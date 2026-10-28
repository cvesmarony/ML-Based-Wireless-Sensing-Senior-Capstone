clear; %test
%% Capture metadata - adjust before every capture
locationName = 'unspecified';   % e.g. 'Halligan 102'
x_lab        = 0;               % RX position in the room grid
y_lab        = 0;
rx_height    = 0.571;           % m
orientation  = 0;               % RX orientation (deg)
captureNotes = '';              % free text: people present, doors open, etc.
data_folder  = 'wifi_data/';
saveFullCapture = false;        % also save the whole raw buffer (~16 MB per 100 ms)

%% Configuration
sampleRate   = 20e6;
CaptureTime  = 0.01;                        % 100 ms burst capture
TotalSamples = round(sampleRate * CaptureTime);
channel      = 6;
centerFreq   = 2407e6 + 5e6 * channel;
rxGain       = 45;
sdrSerial    = '3839223';                  % check with "findsdru"
antennaPort  = 'RX2';

% Packets below minSnrDb are dropped before decoding (even the 6 Mb/s L-SIG
% header is unreliable there). Payloads are only decoded when the SNR meets
% the approximate requirement for the rate given in L-SIG (AWGN, ~10% PER).
minSnrDb      = 5;
legacyRatesMbps = [6 9 12 18 24 36 48 54];
requiredSnrDb   = [5 7 8 11 14 17 21 23];

%% Hardware / Capture Interface
rxObj = comm.SDRuReceiver( ...
    'Platform',           'B210', ...
    'SerialNum',          sdrSerial, ...
    'MasterClockRate',    20e6, ...
    'DecimationFactor',   1, ...
    'SamplesPerFrame',    TotalSamples, ...
    'CenterFrequency',    centerFreq, ...
    'ChannelMapping',     1, ...
    'ReceiveAntennaPort', antennaPort, ...
    'Gain',               rxGain, ...
    'OutputDataType',     'double');

disp('Capturing 20 MHz Wi-Fi frame...');
% Host clock at the start of the capture call; USB latency makes this
% accurate to a few ms, not to the sample.
captureStartTime = datetime('now', 'TimeZone', 'local');
[rxBuffer, ~, overrun] = rxObj();
release(rxObj);

if overrun
    warning('Buffer overrun detected. Some samples may be dropped.');
end

sig = rxBuffer(:, 1);

%% 1. Continuous Packet Detection
% In a 100 ms buffer, there are likely multiple packets.
% Detect all packet candidates above threshold.
threshold = 0.5;
pktStarts = zeros(0, 1);
scanStart = 1;
while scanStart + 159 <= length(sig)
    localPktStart = wlanPacketDetect( ...
        sig(scanStart:end), 'CBW20', 0, threshold);
    if isempty(localPktStart)
        break;
    end

    % wlanPacketDetect returns a 0-based offset into the searched segment.
    pktStart = scanStart + double(localPktStart);
    pktStarts(end+1, 1) = pktStart;

    % Move past the repeated short-training field so this packet is not
    % detected again while leaving later packets in the burst searchable.
    scanStart = max(scanStart + 1, pktStart + 160);
end

if isempty(pktStarts)
    error('No Wi-Fi packets detected in capture buffer.');
end

fprintf('Detected %d packet candidates.\n', length(pktStarts));

% Data containers for training features
csiDataset   = []; % [N_packets x 52] complex
powerDataset = []; % [N_packets x 1] dBm or linear power
snrDataset   = []; % [N_packets x 1] estimated SNR
packetPlotInfo = zeros(length(pktStarts), 6); % boundaries and detection candidate number
plotPacketCount = 0;
packetMacs = repmat({'MAC unavailable'}, length(pktStarts), 1);
packetMacStatus = repmat({'L-SIG was not decoded'}, length(pktStarts), 1);
featureCandidateIds = zeros(length(pktStarts), 1);
featureRouterMacs = cell(length(pktStarts), 1);
featureCfoHz = NaN(length(pktStarts), 1);
packetPowerDb = NaN(length(pktStarts), 1);
packetSnrDb = NaN(length(pktStarts), 1);
packetCfoHz = NaN(length(pktStarts), 1);
packetNoiseVar = NaN(length(pktStarts), 1);
packetRateMbps = NaN(length(pktStarts), 1);
packetPsduLength = NaN(length(pktStarts), 1);
packetFormat = repmat({''}, length(pktStarts), 1);
featureCount = 0;

%% 2. Process Detected Packets
previousPacketEnd = 0;
for p = 1:length(pktStarts)
    pktStart = pktStarts(p);

    % The detector can re-trigger on OFDM data symbols; skip candidates that
    % fall inside a packet whose length was already decoded from its L-SIG.
    if pktStart < previousPacketEnd
        packetMacStatus{p} = 'false detection inside a previous packet';
        continue;
    end

    % Leave room for fine timing to align the packet and for the legacy header.
    if (pktStart + 320) > length(sig)
        packetMacStatus{p} = 'capture ended before the L-SIG field';
        continue;
    end
    
    % Extract preamble segment
    preamble = sig(pktStart : min(pktStart + 400, length(sig)));
    
    % --- Coarse Frequency Offset Correction using L-STF ---
    cfoCoarse = wlanCoarseCFOEstimate(preamble(1:160), 'CBW20');
    preambleCorr = helperFrequencyOffset( ...
        preamble, -cfoCoarse, sampleRate, pktStart);
    
    % --- Fine Symbol Timing / Synchronization ---
    % The returned offset aligns the start of the packet (L-STF).
    fineOffset = wlanSymbolTimingEstimate(preambleCorr, 'CBW20');
    syncStart = pktStart + fineOffset;
    
    if syncStart < 1 || syncStart + 399 > length(sig)
        packetMacStatus{p} = 'capture ended before the synchronized packet header';
        continue;
    end
    
    % L-LTF follows the 160-sample L-STF.
    lltfStart = syncStart + 160;

    % Extract corrected L-LTF (160 samples: 32 CP + 64 T1 + 64 T2)
    lltf_raw = sig(lltfStart : lltfStart + 160 - 1);
    
    % --- Fine Frequency Offset Correction using L-LTF ---
    lltf_coarseCorrected = helperFrequencyOffset( ...
        lltf_raw, -cfoCoarse, sampleRate, lltfStart);
    cfoFine = wlanFineCFOEstimate(lltf_coarseCorrected, 'CBW20');
    cfoTotal = cfoCoarse + cfoFine;
    lltf_corrected = helperFrequencyOffset( ...
        lltf_raw, -cfoTotal, sampleRate, lltfStart);
    
    % --- Demodulate L-LTF and Estimate CSI ---
    % Demodulation performs FFT and discards the cyclic prefix
    lltf_demod = wlanLLTFDemodulate(lltf_corrected, 'CBW20');
    chanEst = wlanLLTFChannelEstimate(lltf_demod, 'CBW20');
    noiseVarEst = wlanLLTFNoiseEstimate(lltf_demod);
    snrEst_dB = 10 * log10((mean(abs(chanEst).^2)) / (noiseVarEst + eps));
    packetSnrDb(p) = snrEst_dB;
    if snrEst_dB < minSnrDb
        packetMacStatus{p} = sprintf('below minimum SNR (%g dB)', minSnrDb);
        continue;
    end

    % Decode L-SIG to estimate the end of this legacy Wi-Fi packet.
    lsigStart = syncStart + 320;
    if lsigStart + 79 <= length(sig)
        lsig = sig(lsigStart : lsigStart + 79);
        lsigCorrected = helperFrequencyOffset( ...
            lsig, -cfoTotal, sampleRate, lsigStart);
        [lsigBits, lsigFail] = wlanLSIGRecover( ...
            lsigCorrected, chanEst, noiseVarEst, 'CBW20');
        if ~lsigFail
            rateCode = sum(double(lsigBits(1:4)) .* (2.^(3:-1:0))');
            psduLength = sum(double(lsigBits(6:17)) .* (2.^(0:11))');
            rateCodes = [13 15 5 7 9 11 1 3];
            dataBitsPerSymbol = [24 36 48 72 96 144 192 216];
            rateIdx = find(rateCodes == rateCode, 1);
            if ~isempty(rateIdx)
                dataSymbols = ceil((16 + 8 * psduLength + 6) / ...
                    dataBitsPerSymbol(rateIdx));

                phyStart = syncStart;
                dataStart = syncStart + 400;
                packetEnd = phyStart + 400 + 80 * dataSymbols;

                % 802.11n/ac/ax packets also start with a legacy L-SIG
                % (6 Mb/s, spoofed length) but their data cannot be decoded
                % as non-HT, so identify them instead of letting them fail
                % the FCS check.
                pktFormat = 'unknown';
                if lsigStart + 239 <= length(sig)
                    formatSamples = helperFrequencyOffset( ...
                        sig(lsigStart:lsigStart + 239), -cfoTotal, ...
                        sampleRate, lsigStart);
                    pktFormat = wlanFormatDetect( ...
                        formatSamples, chanEst, noiseVarEst, 'CBW20');
                end

                % Recover the PSDU and read the BSSID (router radio MAC)
                % from beacons and from data frames to/from an AP.
                if ~strcmp(pktFormat, 'Non-HT') && ~strcmp(pktFormat, 'unknown')
                    packetMacStatus{p} = sprintf( ...
                        '%s packet (only non-HT frames are decoded)', pktFormat);
                elseif snrEst_dB < requiredSnrDb(rateIdx)
                    packetMacStatus{p} = sprintf( ...
                        'SNR too low for %d Mb/s payload (needs ~%d dB)', ...
                        legacyRatesMbps(rateIdx), requiredSnrDb(rateIdx));
                elseif psduLength < 10
                    packetMacStatus{p} = 'L-SIG PSDU length is too short for a MAC header';
                elseif packetEnd > length(sig)
                    packetMacStatus{p} = sprintf( ...
                        'capture ended before packet data (needs sample %d; capture has %d)', ...
                        packetEnd - 1, length(sig));
                else
                    dataSamples = sig(dataStart:packetEnd-1);
                    dataCorrected = helperFrequencyOffset( ...
                        dataSamples, -cfoTotal, sampleRate, dataStart);
                    nonHTConfig = wlanNonHTConfig( ...
                        'ChannelBandwidth', 'CBW20', ...
                        'MCS', rateIdx - 1, ...
                        'PSDULength', psduLength);
                    psduBits = wlanNonHTDataRecover( ...
                        dataCorrected, chanEst, noiseVarEst, nonHTConfig);
                    [packetMacs{p}, packetMacStatus{p}] = ...
                        helperExtractBssid(psduBits);
                end

                packetCfoHz(p) = cfoTotal;
                packetNoiseVar(p) = noiseVarEst;
                packetRateMbps(p) = legacyRatesMbps(rateIdx);
                packetPsduLength(p) = psduLength;
                packetFormat{p} = pktFormat;
                previousPacketEnd = packetEnd;
                plotPacketCount = plotPacketCount + 1;
                packetPlotInfo(plotPacketCount, :) = ...
                    [phyStart, lltfStart, lsigStart, dataStart, packetEnd, p];
            else
                packetMacStatus{p} = 'L-SIG contains an unsupported legacy data rate';
            end
        else
            packetMacStatus{p} = 'L-SIG decoding failed';
        end
    else
        packetMacStatus{p} = 'capture ended before the L-SIG field';
    end
    
    % chanEst is a 52x1 complex vector: H[k]
    % --- Extract Power / RSSI Metrics ---
    % Total packet preamble power (instantaneous power profile)
    pktPwrLinear = mean(abs(lltf_corrected).^2);
    pktPwr_dB = 10 * log10(pktPwrLinear + eps);
    packetPowerDb(p) = pktPwr_dB;

    % Subcarrier-averaged channel power
    csiPower_dB = 10 * log10(mean(abs(chanEst).^2) + eps);

    % Store feature vector
    csiDataset(end+1, :)   = chanEst.';      % 52 complex subcarrier responses
    powerDataset(end+1, :) = [pktPwr_dB, csiPower_dB];
    snrDataset(end+1, :)   = snrEst_dB;
    featureCount = featureCount + 1;
    featureCandidateIds(featureCount) = p;
    featureRouterMacs{featureCount} = packetMacs{p};
    featureCfoHz(featureCount) = cfoTotal;
end
packetPlotInfo = packetPlotInfo(1:plotPacketCount, :);
featureCandidateIds = featureCandidateIds(1:featureCount);
featureRouterMacs = featureRouterMacs(1:featureCount);
featureCfoHz = featureCfoHz(1:featureCount);

%% Report BSSIDs and compare them with the MACs printed on the routers
% The MAC on a router's sticker is usually its Ethernet/base MAC. The BSSID
% it transmits over the air is derived from it (base + a small offset that
% differs per band and per SSID, sometimes with the locally administered
% bit set), so an exact match is not expected.
routerLabelMacs = {'6C:D6:E3:58:34:4C', '6C:D6:E3:58:43:C4', ...
    '6C:D6:E3:58:56:68', '6C:D6:E3:58:35:E4', '6C:D6:E3:58:4A:94'}; % from the router stickers
maxLabelOffset = 64;  % max distance in the low 3 bytes to call a BSSID "related"

[statusNames, ~, statusGroups] = unique(packetMacStatus);
fprintf('MAC decode summary for %d packet candidates (median L-LTF SNR):\n', ...
    length(pktStarts));
for k = 1:length(statusNames)
    fprintf('  %4d x %-55s %6.1f dB\n', sum(statusGroups == k), ...
        statusNames{k}, median(packetSnrDb(statusGroups == k), 'omitnan'));
end

bssidMask = ~strcmp(packetMacs, 'MAC unavailable');
if any(bssidMask)
    observedBssids = packetMacs(bssidMask);
    observedPowers = packetPowerDb(bssidMask);
    [distinctBssids, ~, bssidGroups] = unique(observedBssids);
    bssidPeakPowers = accumarray(bssidGroups, observedPowers, [], @max);
    [rankedPeakPowers, rankedGroups] = sort(bssidPeakPowers, 'descend');
    fprintf('Found %d distinct BSSID(s) in FCS-validated frames:\n', ...
        length(rankedGroups));
    for rank = 1:length(rankedGroups)
        bssid = distinctBssids{rankedGroups(rank)};
        candidateIds = find(strcmp(packetMacs, bssid));
        fprintf('%d. %s (peak power %.2f dB; %d packet(s)) %s\n', ...
            rank, bssid, rankedPeakPowers(rank), length(candidateIds), ...
            helperMatchLabelMac(bssid, routerLabelMacs, maxLabelOffset));
        for candidateIdx = 1:length(candidateIds)
            candidateId = candidateIds(candidateIdx);
            fprintf('   Candidate %d: %.2f dB (%s).\n', ...
                candidateId, packetPowerDb(candidateId), ...
                packetMacStatus{candidateId});
        end
    end
else
    fprintf(['No FCS-validated frames carried a BSSID. Beacons sent at ' ...
        '1 Mb/s (802.11b DSSS) cannot be decoded by this OFDM receiver.\n']);
end

%% Save usable packets for data collection
% Only packets whose BSSID was recovered (FCS-validated, so the whole packet
% is inside the capture) are saved. Times and text are stored as char/double
% (not datetime/string objects) so the file also loads in Python
% (scipy.io.loadmat).
savedRows = find(~strcmp(packetMacs(packetPlotInfo(:, 6)), 'MAC unavailable'));
packetTemplate = struct('iq', [], 'start_sample', [], 'end_sample', [], ...
    'time_iso', '', 'time_posix', [], 'bssid', '', ...
    'mac_status', '', 'format', '', 'rate_mbps', [], 'psdu_length', [], ...
    'snr_db', [], 'power_dbfs', [], 'cfo_hz', [], 'noise_var', [], 'csi', []);
packets = repmat(packetTemplate, length(savedRows), 1);
for k = 1:length(savedRows)
    row = savedRows(k);
    p = packetPlotInfo(row, 6);
    startSample = packetPlotInfo(row, 1);
    endSample = packetPlotInfo(row, 5) - 1;
    packetTime = captureStartTime + seconds((startSample - 1) / sampleRate);

    packets(k).iq = sig(startSample:endSample);   % raw, no CFO correction
    packets(k).start_sample = startSample;         % index into the capture buffer
    packets(k).end_sample = endSample;
    packets(k).time_iso = char(packetTime, 'yyyy-MM-dd''T''HH:mm:ss.SSSSSSXXX');
    packets(k).time_posix = posixtime(packetTime);
    packets(k).bssid = packetMacs{p};
    packets(k).mac_status = packetMacStatus{p};
    packets(k).format = packetFormat{p};
    packets(k).rate_mbps = packetRateMbps(p);      % from L-SIG (6 for HT/VHT/HE)
    packets(k).psdu_length = packetPsduLength(p);  % bytes, from L-SIG
    packets(k).snr_db = packetSnrDb(p);
    packets(k).power_dbfs = packetPowerDb(p);
    packets(k).cfo_hz = packetCfoHz(p);
    packets(k).noise_var = packetNoiseVar(p);
    packets(k).csi = csiDataset(featureCandidateIds == p, :).'; % 52x1 L-LTF estimate
end

wifi_data.packets           = packets;
wifi_data.capture_time_iso  = char(captureStartTime, ...
    'yyyy-MM-dd''T''HH:mm:ss.SSSSSSXXX');
wifi_data.capture_time_posix = posixtime(captureStartTime);
wifi_data.location          = locationName;
wifi_data.rx_pos            = [x_lab, y_lab, rx_height];
wifi_data.rx_orient         = orientation;
wifi_data.notes             = captureNotes;
wifi_data.channel           = channel;
wifi_data.freq              = centerFreq;
wifi_data.Fs                = sampleRate;
wifi_data.CaptureTime       = CaptureTime;
wifi_data.gain_db           = rxGain;
wifi_data.sdr_serial        = sdrSerial;
wifi_data.antenna_port      = antennaPort;
wifi_data.overrun           = logical(overrun);
wifi_data.min_snr_db        = minSnrDb;
wifi_data.num_candidates    = length(pktStarts);
wifi_data.router_label_macs = routerLabelMacs;
wifi_data.matlab_version    = version;
if saveFullCapture
    wifi_data.full_iq = sig;
end

if isempty(packets)
    fprintf('No packets with a recovered BSSID; nothing saved.\n');
else
    if ~isfolder(data_folder)
        mkdir(data_folder);
    end
    filename = sprintf('%swifi_ch%d_%g_%g_%g_%s.mat', data_folder, channel, ...
        x_lab, y_lab, orientation, char(captureStartTime, 'yyyyMMdd_HHmmss'));
    if isfile(filename)
        error('duplicate file name: %s', filename);
    end
    save(filename, 'wifi_data');
    fprintf('Saved %d packet(s) with a BSSID to %s\n', length(packets), filename);
end

%% Raw I/Q waveform with packet-field markers
t = (0:length(sig)-1) / sampleRate;
if ~isempty(packetPlotInfo)
    firstPacketSamples = packetPlotInfo(1, 5) - packetPlotInfo(1, 1);
    firstPacketDurationMs = firstPacketSamples / sampleRate * 1e3;
    fprintf(['First decoded packet duration (legacy OFDM estimate): ' ...
        '%.3f ms (%d samples), derived from L-SIG.\n'], ...
        firstPacketDurationMs, firstPacketSamples);
end

figure('Color', 'w');
plot(t * 1e3, abs(sig));
title('Wi-Fi IQ Magnitude - Full Capture');
xlabel('Time (ms)');
ylabel('|IQ|');
grid on;
hold on;
for p = 1:size(packetPlotInfo, 1)
    startTimeMs = (packetPlotInfo(p, 1) - 1) / sampleRate * 1e3;
    endTimeMs = (packetPlotInfo(p, 5) - 1) / sampleRate * 1e3;
    packetNumber = packetPlotInfo(p, 6);
    xline(startTimeMs, '-', ...
        sprintf('Packet %d | %s', packetNumber, packetMacs{packetNumber}), ...
        'Color', [0 0.45 0.74]);
    xline(endTimeMs, '-', sprintf('Packet %d end', packetNumber), ...
        'Color', [0.85 0.1 0.1], ...
        'LineWidth', 1.5);
end

if ~isempty(packetPlotInfo)
    packetStart = max(1, packetPlotInfo(1, 1));
    packetEndBoundary = packetPlotInfo(1, 5);
    packetEnd = min(length(sig), packetEndBoundary - 1);
    timeRange = packetStart:packetEnd;
    packetStartTimeMs = (packetStart - 1) / sampleRate * 1e3;
    packetEndTimeMs = (packetEndBoundary - 1) / sampleRate * 1e3;
    captureEndTimeMs = length(sig) / sampleRate * 1e3;

    figure('Color', 'w');
    plot(t(timeRange) * 1e3, abs(sig(timeRange)));
    firstCandidateId = packetPlotInfo(1, 6);
    title(sprintf('Packet %d from %s (%.3f ms)', firstCandidateId, ...
        packetMacs{firstCandidateId}, firstPacketDurationMs));
    xlabel('Time (ms)');
    ylabel('|IQ|');
    grid on;
    hold on;

    fieldSamples = packetPlotInfo(1, 1:4);
    fieldTimes = [(fieldSamples - 1) / sampleRate * 1e3, packetEndTimeMs];
    labels = {'Packet start / L-STF', 'L-LTF start', 'L-SIG start', ...
        'DATA start', 'Packet end'};
    for k = 1:length(fieldTimes)
        if k == 1 || k == length(fieldTimes)
            xline(fieldTimes(k), '-', labels{k}, ...
                'Color', [0.85 0.1 0.1], 'LineWidth', 1.5, ...
                'LabelVerticalAlignment', 'middle');
        else
            xline(fieldTimes(k), '--', labels{k}, ...
                'LabelVerticalAlignment', 'middle');
        end
    end
    if packetEndTimeMs <= captureEndTimeMs
        xlim([packetStartTimeMs packetEndTimeMs]);
    else
        xlim([packetStartTimeMs captureEndTimeMs]);
        text(0.99, 0.95, 'Capture ends before estimated packet end', ...
            'Units', 'normalized', 'HorizontalAlignment', 'right', ...
            'Color', [0.85 0.1 0.1]);
    end

    headerEnd = min(length(sig), packetPlotInfo(1, 4) + 159);
    headerRange = packetStart:headerEnd;
    figure('Color', 'w');
    plot(t(headerRange) * 1e3, abs(sig(headerRange)));
    title(sprintf('Preamble and Header - Packet %d from %s', ...
        firstCandidateId, packetMacs{firstCandidateId}));
    xlabel('Time (ms)');
    ylabel('|IQ|');
    grid on;
    detailTimes = fieldTimes(1:4);
    detailLabels = labels(1:4);
    hold on;
    for k = 1:length(detailTimes)
        xline(detailTimes(k), '--', detailLabels{k}, ...
            'LabelVerticalAlignment', 'middle');
    end
else
    warning('Could not decode an L-SIG field; packet end markers are unavailable.');
end

%% 3. Visualization (one figure per plot, sized for slides)
if ~isempty(csiDataset)
    subcarrierIdx = [-26:-1, 1:26];
    nFeatures = size(csiDataset, 1);
    [routerLabels, ~, routerGroup] = unique(featureRouterMacs);
    % Fixed color order; routers past the sixth share gray.
    palette = [0x2a 0x78 0xd6; 0xeb 0x68 0x34; 0x1b 0xaf 0x7a; ...
               0xed 0xa1 0x00; 0xe8 0x7b 0xa4; 0x00 0x83 0x00] / 255;
    routerColors = repmat([0.6 0.6 0.6], length(routerLabels), 1);
    nColored = min(6, length(routerLabels));
    routerColors(1:nColored, :) = palette(1:nColored, :);
    slideSize = [100 100 1100 620];

    % Phase with the linear term (timing offset) and constant term
    % (residual CFO) removed, leaving the channel's own phase shape.
    csiPhase = unwrap(angle(csiDataset), [], 2);
    for packetIdx = 1:nFeatures
        coeffs = polyfit(subcarrierIdx, csiPhase(packetIdx, :), 1);
        csiPhase(packetIdx, :) = csiPhase(packetIdx, :) - ...
            polyval(coeffs, subcarrierIdx);
    end

    % Power delay profile: IFFT of the CSI zero-padded to 512 points
    % (3.2 us span, 6.25 ns bins), aligned to the strongest path since a
    % single unsynchronized receiver has no absolute delay reference.
    nfft = 512;
    dtPdp = 3.2e-6 / nfft;
    leadSamples = round(200e-9 / dtPdp);
    pdpDb = zeros(nFeatures, nfft);
    for packetIdx = 1:nFeatures
        H = zeros(nfft, 1);
        H(mod(subcarrierIdx, nfft) + 1) = csiDataset(packetIdx, :);
        h = ifft(H);
        [~, peakIdx] = max(abs(h));
        h = circshift(h, leadSamples - (peakIdx - 1));
        pdpDb(packetIdx, :) = 20 * log10(abs(h) / max(abs(h)) + eps);
    end
    tauNs = ((0:nfft-1) - leadSamples) * dtPdp * 1e9;
    featureTimeMs = (pktStarts(featureCandidateIds) - 1) / sampleRate * 1e3;

    % CSI magnitude
    figure('Color', 'w', 'Position', slideSize);
    ax = gca; hold(ax, 'on');
    helperPlotByRouter(ax, subcarrierIdx, abs(csiDataset), routerGroup, ...
        routerLabels, routerColors);
    title(ax, 'CSI Magnitude |H[k]| (L-LTF)');
    xlabel(ax, 'Subcarrier Index'); ylabel(ax, '|H|');
    helperStyleSlideAxes(ax);

    % CSI phase
    figure('Color', 'w', 'Position', slideSize);
    ax = gca; hold(ax, 'on');
    helperPlotByRouter(ax, subcarrierIdx, csiPhase, routerGroup, ...
        routerLabels, routerColors);
    title(ax, 'CSI Phase \angleH[k], Linear Trend Removed');
    xlabel(ax, 'Subcarrier Index'); ylabel(ax, 'Phase (rad)');
    helperStyleSlideAxes(ax);

    % Packet power
    figure('Color', 'w', 'Position', slideSize);
    ax = gca; hold(ax, 'on');
    for groupIdx = 1:length(routerLabels)
        groupPackets = routerGroup == groupIdx;
        histogram(ax, powerDataset(groupPackets, 1), 20, ...
            'FaceColor', routerColors(groupIdx, :), 'EdgeColor', 'w', ...
            'DisplayName', sprintf('%s (%d packets)', routerLabels{groupIdx}, ...
            sum(groupPackets)));
    end
    title(ax, 'Received Packet Power Distribution');
    xlabel(ax, 'Power (dBFS)'); ylabel(ax, 'Frame Count');
    helperStyleSlideAxes(ax);

    % Power delay profile
    figure('Color', 'w', 'Position', slideSize);
    ax = gca; hold(ax, 'on');
    helperPlotByRouter(ax, tauNs, pdpDb, routerGroup, routerLabels, routerColors);
    xlim(ax, [-200 1000]); ylim(ax, [-30 1]);
    title(ax, 'Power Delay Profile (IFFT of CSI)');
    xlabel(ax, 'Delay Relative to Strongest Path (ns)');
    ylabel(ax, 'Power (dB, normalized)');
    helperStyleSlideAxes(ax);

    % CFO over time: a transmitter's oscillator offset stays near-constant
    figure('Color', 'w', 'Position', slideSize);
    ax = gca; hold(ax, 'on');
    for groupIdx = 1:length(routerLabels)
        groupPackets = routerGroup == groupIdx;
        scatter(ax, featureTimeMs(groupPackets), featureCfoHz(groupPackets) / 1e3, ...
            80, routerColors(groupIdx, :), 'filled', 'MarkerEdgeColor', 'w', ...
            'DisplayName', routerLabels{groupIdx});
    end
    title(ax, 'Carrier Frequency Offset Over Time');
    xlabel(ax, 'Time in Capture (ms)'); ylabel(ax, 'CFO (kHz)');
    helperStyleSlideAxes(ax);

    % CFO vs power: separate clusters ~ separate transmitters
    figure('Color', 'w', 'Position', slideSize);
    ax = gca; hold(ax, 'on');
    for groupIdx = 1:length(routerLabels)
        groupPackets = routerGroup == groupIdx;
        scatter(ax, powerDataset(groupPackets, 1), featureCfoHz(groupPackets) / 1e3, ...
            80, routerColors(groupIdx, :), 'filled', 'MarkerEdgeColor', 'w', ...
            'DisplayName', routerLabels{groupIdx});
    end
    title(ax, 'CFO vs Packet Power (Transmitter Fingerprint)');
    xlabel(ax, 'Power (dBFS)'); ylabel(ax, 'CFO (kHz)');
    helperStyleSlideAxes(ax);

    % CSI magnitude heatmap, packets in capture order
    figure('Color', 'w', 'Position', slideSize);
    ax = gca;
    imagesc(ax, subcarrierIdx, 1:nFeatures, 20 * log10(abs(csiDataset) + eps));
    blueRamp = [0xcd 0xe2 0xfb; 0x9e 0xc5 0xf4; 0x6d 0xa7 0xec; 0x39 0x87 0xe5; ...
                0x25 0x6a 0xbf; 0x18 0x4f 0x95; 0x0d 0x36 0x6b] / 255;
    colormap(ax, interp1(linspace(0, 1, 7), blueRamp, linspace(0, 1, 256)));
    cb = colorbar(ax); cb.Label.String = '|H| (dB)';
    title(ax, 'CSI Magnitude, All Packets');
    xlabel(ax, 'Subcarrier Index'); ylabel(ax, 'Packet (capture order)');
    ax.FontSize = 16; ax.TitleFontSizeMultiplier = 1.2;
end

%% Helper plotting functions
% One line per packet (rows of y), colored by router; only the first packet
% of each router gets a legend entry.
function helperPlotByRouter(ax, x, y, routerGroup, routerLabels, routerColors)
    for groupIdx = 1:length(routerLabels)
        rows = find(routerGroup == groupIdx);
        for k = 1:length(rows)
            h = plot(ax, x, y(rows(k), :), 'LineWidth', 1.5, ...
                'Color', [routerColors(groupIdx, :) 0.6]);
            if k == 1
                h.DisplayName = sprintf('%s (%d packets)', ...
                    routerLabels{groupIdx}, length(rows));
            else
                h.HandleVisibility = 'off';
            end
        end
    end
end

function helperStyleSlideAxes(ax)
    ax.FontSize = 16;
    ax.TitleFontSizeMultiplier = 1.2;
    grid(ax, 'on'); ax.GridAlpha = 0.15; box(ax, 'off');
    legend(ax, 'Location', 'bestoutside', 'Interpreter', 'none');
end

%% Helper CFO correction function
function out = helperFrequencyOffset(in, fOffset, fs, startSample)
    t = (startSample - 1 + (0:length(in)-1).') / fs;
    out = in .* exp(1j * 2 * pi * fOffset * t);
end

function [macAddress, statusMessage] = helperExtractBssid(psduBits)
    macAddress = 'MAC unavailable';
    [macConfig, ~, decodeStatus] = wlanMPDUDecode( ...
        psduBits, 'Non-HT', 'SuppressWarnings', true);
    if strcmp(string(decodeStatus), 'UnsupportedFrameSubtype')
        % MATLAB does not parse e.g. probe responses or authentication
        % frames, but all management frames carry the BSSID in Address 3.
        [macAddress, statusMessage] = helperParseManagementHeader(psduBits);
        return;
    end
    if decodeStatus ~= 0 || ~macConfig.Decoded
        statusMessage = sprintf('MPDU failed validation (%s)', ...
            char(string(decodeStatus)));
        return;
    end

    % Which address field holds the BSSID depends on the frame type and,
    % for data frames, on the To DS / From DS bits (IEEE 802.11 Table 9-30).
    frameType = char(macConfig.FrameType);
    switch frameType
        case {'Beacon', 'Action', 'Action No Ack'}
            address = macConfig.Address3;
            source = frameType;
        case {'Data', 'Null', 'QoS Data', 'QoS Null'}
            if macConfig.ToDS && ~macConfig.FromDS
                address = macConfig.Address1;
                source = [frameType ' client->AP'];
            elseif ~macConfig.ToDS && macConfig.FromDS
                address = macConfig.Address2;
                source = [frameType ' AP->client'];
            elseif ~macConfig.ToDS && ~macConfig.FromDS
                address = macConfig.Address3;
                source = [frameType ' no-DS'];
            else
                statusMessage = sprintf( ...
                    'decoded 4-address %s frame; BSSID not reported', frameType);
                return;
            end
        otherwise
            statusMessage = sprintf( ...
                'decoded %s frame carries no BSSID', frameType);
            return;
    end

    address = char(address);
    if length(address) ~= 12 || bitand(hex2dec(address(1:2)), 1)
        statusMessage = sprintf( ...
            'decoded %s frame has a group address in the BSSID field', frameType);
        return;
    end

    macAddress = sprintf('%s:%s:%s:%s:%s:%s', ...
        address(1:2), address(3:4), address(5:6), ...
        address(7:8), address(9:10), address(11:12));
    statusMessage = sprintf('BSSID from %s frame', source);
end

function [macAddress, statusMessage] = helperParseManagementHeader(psduBits)
    macAddress = 'MAC unavailable';
    % PSDU bits are sent least significant bit first within each octet.
    octets = (2.^(0:7)) * reshape(double(psduBits), 8, []);
    if length(octets) < 28 || ~helperFcsValid(octets)
        statusMessage = 'unsupported frame subtype failed FCS check';
        return;
    end

    frameType = bitand(bitshift(octets(1), -2), 3);
    subtype = bitshift(octets(1), -4);
    if frameType ~= 0
        statusMessage = sprintf( ...
            'unsupported type %d subtype %d frame carries no BSSID', ...
            frameType, subtype);
        return;
    end

    subtypeNames = {'Assoc Request', 'Assoc Response', 'Reassoc Request', ...
        'Reassoc Response', 'Probe Request', 'Probe Response', 'Timing Advert', ...
        'Reserved', 'Beacon', 'ATIM', 'Disassociation', 'Authentication', ...
        'Deauthentication', 'Action', 'Action No Ack', 'Reserved'};
    frameName = subtypeNames{subtype + 1};
    bssidOctets = octets(17:22);
    if bitand(bssidOctets(1), 1)
        statusMessage = sprintf( ...
            'decoded %s frame has a group address in the BSSID field', frameName);
        return;
    end

    macAddress = strjoin(cellstr(dec2hex(bssidOctets, 2)), ':');
    statusMessage = sprintf('BSSID from %s frame', frameName);
end

function isValid = helperFcsValid(octets)
    % IEEE 802.11 FCS: reflected CRC-32 over the frame, sent LSB first.
    crc = uint32(4294967295);
    for k = 1:length(octets) - 4
        crc = bitxor(crc, uint32(octets(k)));
        for bit = 1:8
            if bitand(crc, 1)
                crc = bitxor(bitshift(crc, -1), uint32(3988292384));
            else
                crc = bitshift(crc, -1);
            end
        end
    end
    receivedFcs = octets(end-3:end) * (256.^(0:3)).';
    isValid = double(bitcmp(crc)) == receivedFcs;
end

function matchText = helperMatchLabelMac(bssid, labelMacs, maxOffset)
    matchText = '';
    bssidBytes = helperMacBytes(bssid);
    for k = 1:length(labelMacs)
        labelBytes = helperMacBytes(labelMacs{k});
        % Ignore the locally administered bit (0x02), which many APs set on
        % their extra virtual BSSIDs.
        sameOui = bitand(bssidBytes(1), 253) == bitand(labelBytes(1), 253) && ...
            all(bssidBytes(2:3) == labelBytes(2:3));
        offset = (bssidBytes(4:6) - labelBytes(4:6)) * [65536; 256; 1];
        if sameOui && offset == 0 && bssidBytes(1) == labelBytes(1)
            matchText = sprintf('<-- exact match to label %s', labelMacs{k});
            return;
        elseif sameOui && abs(offset) <= maxOffset
            matchText = sprintf('<-- likely router with label %s (offset %+d)', ...
                labelMacs{k}, offset);
            return;
        end
    end
end

function bytes = helperMacBytes(mac)
    hexDigits = regexprep(mac, '[^0-9A-Fa-f]', '');
    bytes = hex2dec(reshape(hexDigits, 2, []).').';
end