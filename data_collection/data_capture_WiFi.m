clear; %test
%% Configuration
sampleRate   = 20e6;
CaptureTime  = 0.1;                        % 100 ms burst capture
TotalSamples = round(sampleRate * CaptureTime);
channel      = 6;
centerFreq   = 2407e6 + 5e6 * channel;

%% Hardware / Capture Interface
rxObj = comm.SDRuReceiver( ...
    'Platform',           'B210', ...
    'SerialNum',          '34C78FD', ... #will have to check this with "findsdru"
    'MasterClockRate',    20e6, ...
    'DecimationFactor',   1, ...
    'SamplesPerFrame',    TotalSamples, ...
    'CenterFrequency',    centerFreq, ...
    'ChannelMapping',     1, ...
    'ReceiveAntennaPort', 'RX2', ...
    'Gain',               45, ...
    'OutputDataType',     'double');

disp('Capturing 20 MHz Wi-Fi frame...');
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

    pktStart = scanStart + double(localPktStart) - 1;
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
packetPowerDb = NaN(length(pktStarts), 1);
featureCount = 0;

%% 2. Process Detected Packets
for p = 1:length(pktStarts)
    pktStart = pktStarts(p);
    
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

                % Recover the PSDU to identify beacon frames and their BSSID.
                if psduLength < 10
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
                        helperExtractBeaconBssid(psduBits);
                end

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
    snrEst_dB   = 10 * log10((mean(abs(chanEst).^2)) / (noiseVarEst + eps));
    
    % Store feature vector
    csiDataset(end+1, :)   = chanEst.';      % 52 complex subcarrier responses
    powerDataset(end+1, :) = [pktPwr_dB, csiPower_dB];
    snrDataset(end+1, :)   = snrEst_dB;
    featureCount = featureCount + 1;
    featureCandidateIds(featureCount) = p;
    featureRouterMacs{featureCount} = packetMacs{p};
end
packetPlotInfo = packetPlotInfo(1:plotPacketCount, :);
featureCandidateIds = featureCandidateIds(1:featureCount);
featureRouterMacs = featureRouterMacs(1:featureCount);

macPrefix = '6C:D6:E3:58:';
validatedBssidMask = ~strcmp(packetMacs, 'MAC unavailable') & ...
    isfinite(packetPowerDb) & ...
    strcmp(packetMacStatus, 'validated beacon BSSID');
validBssidMask = validatedBssidMask & ...
    strncmp(packetMacs, macPrefix, length(macPrefix));
if any(validBssidMask)
    observedBssids = packetMacs(validBssidMask);
    observedPowers = packetPowerDb(validBssidMask);
    [distinctBssids, ~, bssidGroups] = unique(observedBssids);
    bssidPeakPowers = accumarray(bssidGroups, observedPowers, [], @max);
    [rankedPeakPowers, rankedGroups] = sort(bssidPeakPowers, 'descend');
    fprintf('Found %d distinct validated beacon BSSID(s):\n', ...
        length(rankedGroups));
    for rank = 1:length(rankedGroups)
        bssid = distinctBssids{rankedGroups(rank)};
        candidateIds = find(strcmp(packetMacs, bssid));
        fprintf('%d. %s (peak power %.2f dB; %d packet(s))\n', ...
            rank, bssid, rankedPeakPowers(rank), length(candidateIds));
        for candidateIdx = 1:length(candidateIds)
            candidateId = candidateIds(candidateIdx);
            fprintf('   Candidate %d: %.2f dB (%s).\n', ...
                candidateId, packetPowerDb(candidateId), ...
                packetMacStatus{candidateId});
        end
    end
else
    if any(validatedBssidMask)
        fprintf(['No validated beacon BSSIDs matched prefix %s ' ...
            '(%d other validated BSSID(s) found).\n'], ...
            macPrefix, sum(validatedBssidMask));
    else
        fprintf(['No beacon BSSIDs passed MPDU validation; ' ...
            'non-beacon frames are not reported as router BSSIDs.\n']);
    end
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
for p = 1:length(pktStarts)
    if any(packetPlotInfo(:, 6) == p)
        continue;
    end
    candidateTime = (pktStarts(p) - 1) / sampleRate * 1e3;
    xline(candidateTime, ':', sprintf('Candidate %d (not decoded)', p));
end
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

%% 3. Visualization
if ~isempty(csiDataset)
    figure('Color', 'w');

    subcarrierIdx = [-26:-1, 1:26];
    [routerLabels, ~, routerGroup] = unique(featureRouterMacs);
    routerColors = lines(length(routerLabels));
    
    subplot(3,1,1);
    hold on;
    for packetIdx = 1:size(csiDataset, 1)
        groupIdx = routerGroup(packetIdx);
        plot(subcarrierIdx, abs(csiDataset(packetIdx, :)), '-o', ...
            'Color', routerColors(groupIdx, :), ...
            'DisplayName', sprintf('Packet %d | %s', ...
            featureCandidateIds(packetIdx), featureRouterMacs{packetIdx}));
    end
    title('Wi-Fi Subcarrier Channel Magnitude (|H[k]|)');
    xlabel('Subcarrier Index'); ylabel('|H|'); grid on; legend('show');
    
    subplot(3,1,2);
    hold on;
    for packetIdx = 1:size(csiDataset, 1)
        groupIdx = routerGroup(packetIdx);
        plot(subcarrierIdx, unwrap(angle(csiDataset(packetIdx, :))), '-o', ...
            'Color', routerColors(groupIdx, :), ...
            'DisplayName', sprintf('Packet %d | %s', ...
            featureCandidateIds(packetIdx), featureRouterMacs{packetIdx}));
    end
    title('Wi-Fi Subcarrier Phase (\angle H[k])');
    xlabel('Subcarrier Index'); ylabel('Phase (rad)'); grid on; legend('show');
    
    subplot(3,1,3);
    hold on;
    for groupIdx = 1:length(routerLabels)
        groupPackets = routerGroup == groupIdx;
        histogram(powerDataset(groupPackets, 1), 20, ...
            'DisplayName', sprintf('%s (%d packets)', routerLabels{groupIdx}, ...
            sum(groupPackets)));
    end
    title('Received Packet Power Distribution (dB)');
    xlabel('Power (dBFS)'); ylabel('Frame Count'); grid on; legend('show');
end

%% Helper CFO correction function
function out = helperFrequencyOffset(in, fOffset, fs, startSample)
    t = (startSample - 1 + (0:length(in)-1).') / fs;
    out = in .* exp(1j * 2 * pi * fOffset * t);
end

function [macAddress, statusMessage] = helperExtractBeaconBssid(psduBits)
    macAddress = 'MAC unavailable';
    [macConfig, ~, decodeStatus] = wlanMPDUDecode( ...
        psduBits, 'Non-HT', 'SuppressWarnings', true);
    if decodeStatus ~= 0 || ~macConfig.Decoded
        statusMessage = sprintf( ...
            'MPDU failed validation (wlanMPDUDecode status %d)', decodeStatus);
        return;
    end

    frameType = char(macConfig.FrameType);
    if ~strcmpi(frameType, 'Beacon')
        statusMessage = sprintf( ...
            'decoded %s frame is not a beacon; BSSID not reported', frameType);
        return;
    end

    address = char(macConfig.Address3);
    if length(address) ~= 12
        statusMessage = 'decoded beacon has no valid BSSID';
        return;
    end

    macAddress = sprintf('%s:%s:%s:%s:%s:%s', ...
        address(1:2), address(3:4), address(5:6), ...
        address(7:8), address(9:10), address(11:12));
    statusMessage = 'validated beacon BSSID';
end