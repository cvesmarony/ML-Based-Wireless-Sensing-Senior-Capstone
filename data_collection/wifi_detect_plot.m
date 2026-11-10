clear;
% Configures the B210, listens for captureMs milliseconds on both antennas,
% runs EVERY detected packet candidate through the Wi-Fi decoder, prints
% where each candidate ended up and why (with the evidence), and plots
% |IQ| vs time with a line where each real packet begins and ends.
% Nothing is saved. Decoding follows data_capture_WiFi.m (MRC over both
% antennas).

%% Configuration
captureMs   = 100;              % how long to listen (ms)
sampleRate  = 20e6;
channel     = 6;
centerFreq  = 2407e6 + 5e6 * channel;
rxGain      = 45;               % dB, same on both channels
sdrSerial   = '3839223';        % check with "findsdru"
antennaPort = 'RX2';
numRx       = 1;                % 2 antennas at 20 MS/s overflows this laptop's USB link

% Below minSnrDb even the most robust part of a packet (the 6 Mb/s L-SIG
% header) cannot be decoded. requiredSnrDb is the approximate SNR each
% payload rate needs (AWGN, ~10% packet error rate); it is only used as
% evidence in the report - every payload is decoded regardless.
minSnrDb        = 5;
legacyRatesMbps = [6 9 12 18 24 36 48 54];
requiredSnrDb   = [5 7 8 11 14 17 21 23];

%% Capture
fprintf('Capturing %g ms on %d antennas, channel %d (%.0f MHz)...\n', ...
    captureMs, numRx, channel, centerFreq / 1e6);
[sig, overrun] = captureB210(sdrSerial, centerFreq, sampleRate, rxGain, ...
    antennaPort, numRx, captureMs / 1e3);     % [nSamples x numRx]
if overrun
    warning('Buffer overrun detected. Some samples may be dropped.');
end
nSamples = size(sig, 1);
% The ADC saturates at |I| or |Q| = 1; clipped samples corrupt packets.
isClipped = any(abs(real(sig)) >= 0.999 | abs(imag(sig)) >= 0.999, 2);
if any(isClipped)
    warning(['%d samples are clipped at ADC full scale: lower rxGain ' ...
        '(currently %g dB).'], sum(isClipped), rxGain);
end

%% Packet detection
% wlanPacketDetect flags anything whose autocorrelation looks like the
% repeating Wi-Fi short preamble. These are CANDIDATES, not packets.
pktStarts = zeros(0, 1);
scanStart = 1;
while scanStart + 159 <= nSamples
    localPktStart = wlanPacketDetect(sig(scanStart:end, :), 'CBW20', 0, 0.5);
    if isempty(localPktStart)
        break;
    end
    pktStart = scanStart + double(localPktStart);   % 0-based offset
    pktStarts(end+1, 1) = pktStart; %#ok<SAGROW>
    scanStart = max(scanStart + 1, pktStart + 160);  % skip past this L-STF
end
nCand = length(pktStarts);
if nCand == 0
    fprintf('No Wi-Fi packet candidates detected in %g ms.\n', captureMs);
end

%% Decode every candidate
% Outcome categories, in the order they are reported.
DECODED_MAC   = 1;  DECODED_NOMAC = 2;  CRC_FAIL = 3;  HT_PAYLOAD = 4;
CUT_OFF       = 5;  RETRIGGER     = 6;  LOW_SNR  = 7;  BAD_HEADER = 8;
FAKE_HEADER   = 9;

category   = zeros(nCand, 1);
reason     = repmat({''}, nCand, 1);
packetMac  = repmat({'MAC unavailable'}, nCand, 1);
snrDb      = NaN(nCand, 1);          % post-MRC L-LTF SNR
pktBounds  = NaN(nCand, 2);          % [first sample, one past last sample]
candMs     = (pktStarts - 1) / sampleRate * 1e3;
hostPacket = zeros(nCand, 1);        % for re-triggers: packet they fall in

previousPacketEnd = 0;
previousPacketId  = 0;
for p = 1:nCand
    pktStart = pktStarts(p);

    % The detector re-fires on OFDM data symbols inside a packet whose end
    % we already know from its header's length field.
    if pktStart < previousPacketEnd
        category(p) = RETRIGGER;
        hostPacket(p) = previousPacketId;
        continue;
    end
    if pktStart + 400 > nSamples
        category(p) = CUT_OFF;
        reason{p} = 'capture ended during the preamble';
        continue;
    end

    % Coarse CFO (L-STF), fine timing
    preamble = sig(pktStart:pktStart + 400, :);
    cfoCoarse = wlanCoarseCFOEstimate(preamble(1:160, :), 'CBW20');
    preambleCorr = helperFrequencyOffset(preamble, -cfoCoarse, sampleRate, pktStart);
    syncStart = pktStart + wlanSymbolTimingEstimate(preambleCorr, 'CBW20');
    if syncStart < 1 || syncStart + 399 > nSamples
        category(p) = CUT_OFF;
        reason{p} = 'capture ended during the preamble';
        continue;
    end

    % Fine CFO, channel and SNR from the L-LTF
    lltfStart = syncStart + 160;
    lltfRaw = sig(lltfStart:lltfStart + 159, :);
    cfoFine = wlanFineCFOEstimate( ...
        helperFrequencyOffset(lltfRaw, -cfoCoarse, sampleRate, lltfStart), 'CBW20');
    cfoTotal = cfoCoarse + cfoFine;
    lltfDemod = wlanLLTFDemodulate( ...
        helperFrequencyOffset(lltfRaw, -cfoTotal, sampleRate, lltfStart), 'CBW20');
    chanEst = wlanLLTFChannelEstimate(lltfDemod, 'CBW20');
    noiseVar = wlanLLTFNoiseEstimate(lltfDemod);
    snrPerRx = squeeze(mean(abs(chanEst).^2, 1)).' / (noiseVar + eps);
    snrDb(p) = 10 * log10(sum(snrPerRx));
    if snrDb(p) < minSnrDb
        category(p) = LOW_SNR;
        continue;
    end

    % L-SIG header (always 6 Mb/s): rate and length of the payload
    lsigStart = syncStart + 320;
    lsig = helperFrequencyOffset(sig(lsigStart:lsigStart + 79, :), ...
        -cfoTotal, sampleRate, lsigStart);
    [lsigBits, lsigFail] = wlanLSIGRecover(lsig, chanEst, noiseVar, 'CBW20');
    if lsigFail
        category(p) = BAD_HEADER;
        reason{p} = sprintf('SNR %.1f dB; header parity/rate check failed', snrDb(p));
        continue;
    end
    rateCode = sum(double(lsigBits(1:4)) .* (2.^(3:-1:0))');
    psduLength = sum(double(lsigBits(6:17)) .* (2.^(0:11))');
    rateIdx = find([13 15 5 7 9 11 1 3] == rateCode, 1);
    if isempty(rateIdx)
        category(p) = BAD_HEADER;
        reason{p} = sprintf('SNR %.1f dB; header has invalid rate code %d', ...
            snrDb(p), rateCode);
        continue;
    end
    if psduLength < 14
        category(p) = BAD_HEADER;
        reason{p} = sprintf(['SNR %.1f dB; header says %d bytes, shorter than ' ...
            'the smallest 802.11 frame (14 B)'], snrDb(p), psduLength);
        continue;
    end
    rateMbps = legacyRatesMbps(rateIdx);
    dataBitsPerSymbol = [24 36 48 72 96 144 192 216];
    dataSymbols = ceil((16 + 8 * psduLength + 6) / dataBitsPerSymbol(rateIdx));
    packetEnd = syncStart + 400 + 80 * dataSymbols;
    pktBounds(p, :) = [syncStart, packetEnd];

    % 802.11n/ac/ax packets start with the same legacy header
    pktFormat = 'unknown';
    if lsigStart + 239 <= nSamples
        pktFormat = wlanFormatDetect(helperFrequencyOffset( ...
            sig(lsigStart:lsigStart + 239, :), -cfoTotal, sampleRate, lsigStart), ...
            chanEst, noiseVar, 'CBW20');
    end
    if strcmp(pktFormat, 'HT-GF')
        % Greenfield packets have no L-SIG at all, so a "valid" L-SIG plus
        % a greenfield verdict is self-contradictory: the header was noise
        % that happened to pass the 1-bit parity check. Its length field is
        % not trusted for re-trigger rejection either.
        category(p) = FAKE_HEADER;
        reason{p} = sprintf(['SNR %.1f dB; format check says HT-greenfield, ' ...
            'which has no L-SIG, so the "header" was noise'], snrDb(p));
        pktBounds(p, :) = NaN;
        continue;
    end
    previousPacketEnd = packetEnd;
    previousPacketId = p;
    headerText = sprintf('%s, %d Mb/s, %d B, SNR %.1f dB', ...
        pktFormat, rateMbps, psduLength, snrDb(p));

    if ~any(strcmp(pktFormat, {'Non-HT', 'unknown'}))
        category(p) = HT_PAYLOAD;
        reason{p} = sprintf('%s, header length %d B, SNR %.1f dB', ...
            pktFormat, psduLength, snrDb(p));
    elseif packetEnd > nSamples
        category(p) = CUT_OFF;
        reason{p} = sprintf('%s; ends at sample %d, capture has %d', ...
            headerText, packetEnd - 1, nSamples);
    else
        % Payload: always attempted, whatever the SNR. The CRC-32 checksum
        % (FCS) at the end of the frame decides whether the bits are right.
        dataStart = syncStart + 400;
        dataCorr = helperFrequencyOffset(sig(dataStart:packetEnd - 1, :), ...
            -cfoTotal, sampleRate, dataStart);
        cfg = wlanNonHTConfig('ChannelBandwidth', 'CBW20', ...
            'MCS', rateIdx - 1, 'PSDULength', psduLength);
        psduBits = wlanNonHTDataRecover(dataCorr, chanEst, noiseVar, cfg);
        [packetMac{p}, frameText, crcOk] = helperExtractBssid(psduBits);
        if ~crcOk
            category(p) = CRC_FAIL;
            reason{p} = sprintf('%s; %s; %d Mb/s needs ~%d dB; %d clipped samples', ...
                headerText, frameText, rateMbps, requiredSnrDb(rateIdx), ...
                sum(isClipped(syncStart:packetEnd - 1)));
        elseif ~strcmp(packetMac{p}, 'MAC unavailable')
            category(p) = DECODED_MAC;
            reason{p} = sprintf('%s; %s -> %s', headerText, frameText, packetMac{p});
        else
            category(p) = DECODED_NOMAC;
            reason{p} = sprintf('%s; %s', headerText, frameText);
        end
    end
end

%% Report
catTitle = {
    'DECODED: checksum passed, router MAC (BSSID) read'
    'DECODED: checksum passed, but this frame type has no router-MAC field'
    'REAL PACKET, NOT DECODED: payload decoded but checksum failed'
    'REAL PACKET, NOT DECODED: 802.11n/ac/ax payload'
    'CUT OFF: packet runs past the edge of the capture'
    'NOT A PACKET: detector re-fired inside an earlier packet'
    'NOT A PACKET: preamble SNR too low (noise)'
    'NOT A PACKET: header failed its validity check'
    'NOT A PACKET: header passed its parity check by chance'};
catWhy = {
    {'Header, payload and 32-bit CRC checksum all passed: every bit of the frame', ...
     'is verified correct. The router MAC is read from the MAC header.'}
    {'Checksum passed, so these are fully decoded real frames. ACK/CTS-type frames', ...
     'only carry a receiver address, so there is no router MAC to report.'}
    {'The 6 Mb/s header decoded, so a real packet was there, and its payload was', ...
     'demodulated. The checksum failed, so the payload bits are wrong. Evidence on', ...
     'each line: SNR below what its data rate needs -> too weak/far; SNR above it ->', ...
     'corrupted mid-packet (clipped samples > 0 -> gain too high; otherwise most', ...
     'likely another transmitter overlapping it).'}
    {'Header decoded and the format check identified 802.11n/ac/ax. The payload uses', ...
     'HT/VHT/HE modulation, which this script does not demodulate (802.11a/g only).'}
    {'The packet started or ended outside the capture window, so it is incomplete.'}
    {'The detector looks for the repeating preamble pattern and also fires on OFDM', ...
     'data symbols. Evidence: each of these starts between the start and end of a', ...
     'packet whose length we already read from its header.'}
    {sprintf(['Preamble SNR below %g dB. Below this not even the most robust part of ' ...
     'a packet'], minSnrDb), ...
     '(the 6 Mb/s header) can be decoded: these are noise bumps that crossed the detector.'}
    {'SNR was high enough, but the header failed its parity/rate/length check:', ...
     'noise or two overlapping packets, not a decodable Wi-Fi header.'}
    {'L-SIG has a single parity bit, so about 1 in 4 noise bursts produce a "valid"', ...
     'header. The format check then contradicts it (greenfield has no L-SIG).'}};

counts = accumarray(category, 1, [numel(catTitle) 1]);
nDecoded = sum(counts([DECODED_MAC DECODED_NOMAC]));
nRealNot = sum(counts([CRC_FAIL HT_PAYLOAD CUT_OFF]));
nFalse   = sum(counts([RETRIGGER LOW_SNR BAD_HEADER FAKE_HEADER]));
rule = repmat('=', 1, 78);

fprintf('\n%s\n RESULTS: %g ms capture, channel %d, %d antennas (MRC)\n%s\n', ...
    rule, captureMs, channel, numRx, rule);
fprintf(['The detector flagged %d candidates: anything that looks like a Wi-Fi ' ...
    'preamble.\nALL %d were run through the decoder. Each stops at the first ' ...
    'stage it fails:\n  preamble SNR -> 6 Mb/s header (L-SIG) -> format check ' ...
    '-> payload -> CRC-32 checksum\n\n'], nCand, nCand);
fprintf('  %-58s %5d\n', 'SUCCESSFULLY DECODED (checksum passed)', nDecoded);
fprintf('  %-58s %5d\n', 'REAL PACKETS NOT DECODED (or cut off)', nRealNot);
fprintf('  %-58s %5d\n', 'NOT PACKETS (detector false alarms)', nFalse);
fprintf('  %s\n  %-58s %5d\n', repmat('-', 1, 64), 'Total candidates', nCand);

for c = 1:numel(catTitle)
    rows = find(category == c);
    if isempty(rows)
        continue;
    end
    fprintf('\n[%d] %s: %d\n', c, catTitle{c}, numel(rows));
    fprintf('    %s\n', catWhy{c}{:});
    switch c
        case RETRIGGER
            hosts = unique(hostPacket(rows));
            for h = hosts'
                n = sum(hostPacket(rows) == h);
                fprintf('    %3d re-trigger(s) inside #%d (%.3f-%.3f ms, %s)\n', ...
                    n, h, (pktBounds(h, :) - 1) / sampleRate * 1e3, ...
                    catTitle{category(h)}(1:find(catTitle{category(h)} == ':', 1) - 1));
            end
        case LOW_SNR
            fprintf('    SNR of these candidates: min %.1f, median %.1f, max %.1f dB\n', ...
                min(snrDb(rows)), median(snrDb(rows)), max(snrDb(rows)));
        otherwise
            for p = rows'
                fprintf('    #%-4d @ %8.3f ms  %s\n', p, candMs(p), reason{p});
            end
    end
end

macRows = find(category == DECODED_MAC);
if ~isempty(macRows)
    [macs, ~, g] = unique(packetMac(macRows));
    fprintf('\nRouter MACs (BSSIDs) heard: %d\n', numel(macs));
    for k = 1:numel(macs)
        fprintf('    %s  %d frame(s), best SNR %.1f dB\n', macs{k}, ...
            sum(g == k), max(snrDb(macRows(g == k))));
    end
end
fprintf('%s\n', rule);

%% Plot: |IQ| vs time with packet start/end lines
% Only real packets (decoded or not) get lines; false alarms do not. Start
% lines carry the candidate number (and MAC if decoded) so they can be
% matched to the report; zoom in to separate packets close together.
t = (0:nSamples - 1) / sampleRate * 1e3;
green = [0 0.6 0.2];  orange = [0.93 0.55 0];  red = [0.85 0.1 0.1];
figure('Color', 'w', 'Position', [100 100 1400 600]);
hold on; grid on;
hTrace = plot(t, abs(sig), 'LineWidth', 0.5);
traceColors = [0.55 0.55 0.55; 0.75 0.80 0.95];
for r = 1:numel(hTrace)
    hTrace(r).Color = traceColors(min(r, end), :);
end
plot(NaN, NaN, '-', 'Color', green, 'LineWidth', 1.5);
plot(NaN, NaN, '-', 'Color', orange, 'LineWidth', 1.5);
plot(NaN, NaN, ':', 'Color', red, 'LineWidth', 1.5);
legend([arrayfun(@(r) sprintf('RX %d', r), 1:numRx, 'UniformOutput', false), ...
    {'Start: decoded', 'Start: real, not decoded', 'Packet end'}], ...
    'Location', 'northeastoutside');
for p = find(category <= CUT_OFF & ~isnan(pktBounds(:, 1)))'
    startMs = (pktBounds(p, 1) - 1) / sampleRate * 1e3;
    if category(p) == DECODED_MAC
        xline(startMs, '-', sprintf('#%d %s', p, packetMac{p}), ...
            'Color', green, 'LineWidth', 1.5, 'HandleVisibility', 'off');
    elseif category(p) == DECODED_NOMAC
        xline(startMs, '-', sprintf('#%d decoded', p), ...
            'Color', green, 'LineWidth', 1.5, 'HandleVisibility', 'off');
    else
        xline(startMs, '-', sprintf('#%d', p), ...
            'Color', orange, 'LineWidth', 1.2, 'HandleVisibility', 'off');
    end
    if pktBounds(p, 2) <= nSamples
        xline((pktBounds(p, 2) - 1) / sampleRate * 1e3, ':', ...
            'Color', red, 'LineWidth', 1.2, 'HandleVisibility', 'off');
    end
end
title(sprintf('Wi-Fi IQ Magnitude - %g ms, channel %d (%d decoded, %d real not decoded)', ...
    captureMs, channel, nDecoded, nRealNot));
xlabel('Time (ms)');
ylabel('|IQ|');
xlim([0 t(end)]);

%% Helpers
function out = helperFrequencyOffset(in, fOffset, fs, startSample)
    t = (startSample - 1 + (0:size(in, 1) - 1).') / fs;
    out = in .* exp(1j * 2 * pi * fOffset * t);
end

% crcOk is true when the frame's CRC-32 (FCS) passed, whether or not the
% frame type carries a BSSID.
function [macAddress, statusMessage, crcOk] = helperExtractBssid(psduBits)
    macAddress = 'MAC unavailable';
    crcOk = false;
    [macConfig, ~, decodeStatus] = wlanMPDUDecode( ...
        psduBits, 'Non-HT', 'SuppressWarnings', true);
    if strcmp(string(decodeStatus), 'UnsupportedFrameSubtype')
        % MATLAB does not parse e.g. probe responses, but all management
        % frames carry the BSSID in Address 3.
        [macAddress, statusMessage, crcOk] = helperParseManagementHeader(psduBits);
        return;
    end
    if decodeStatus ~= 0 || ~macConfig.Decoded
        if strcmp(string(decodeStatus), 'FCSFailed')
            statusMessage = 'CRC-32 checksum FAILED';
        else
            statusMessage = sprintf('frame validation failed (%s)', ...
                char(string(decodeStatus)));
        end
        return;
    end
    crcOk = true;

    % Which address holds the BSSID depends on frame type and To/From DS
    % bits (IEEE 802.11 Table 9-30).
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
                statusMessage = sprintf('checksum OK, 4-address %s frame', frameType);
                return;
            end
        otherwise
            statusMessage = sprintf('checksum OK, %s frame (no BSSID field)', frameType);
            return;
    end

    address = char(address);
    if length(address) ~= 12 || bitand(hex2dec(address(1:2)), 1)
        statusMessage = sprintf('checksum OK, %s frame with group address', frameType);
        return;
    end
    macAddress = sprintf('%s:%s:%s:%s:%s:%s', address(1:2), address(3:4), ...
        address(5:6), address(7:8), address(9:10), address(11:12));
    statusMessage = sprintf('checksum OK, %s frame', source);
end

function [macAddress, statusMessage, crcOk] = helperParseManagementHeader(psduBits)
    macAddress = 'MAC unavailable';
    % PSDU bits are sent least significant bit first within each octet.
    octets = (2.^(0:7)) * reshape(double(psduBits), 8, []);
    crcOk = length(octets) >= 28 && helperFcsValid(octets);
    if ~crcOk
        statusMessage = 'CRC-32 checksum FAILED';
        return;
    end
    frameType = bitand(bitshift(octets(1), -2), 3);
    subtype = bitshift(octets(1), -4);
    if frameType ~= 0
        statusMessage = sprintf('checksum OK, type %d subtype %d frame (no BSSID field)', ...
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
        statusMessage = sprintf('checksum OK, %s frame with group address', frameName);
        return;
    end
    macAddress = strjoin(cellstr(dec2hex(bssidOctets, 2)), ':');
    statusMessage = sprintf('checksum OK, %s frame', frameName);
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
