function stats = wifiCaptureCheck(iq, sampleRate, titleText)
% Fast health check of one raw Wi-Fi capture, run after every capture.
%
% The metric is the per-packet preamble SNR: for every real 802.11 packet
% (its L-SIG header passed the parity, rate, length and format checks) the
% SNR is measured on the L-LTF, the same field the CSI comes from. A good
% capture has many packets at a usable SNR, no dropped samples and no
% clipping. The figure shows
%   top:    SNR of each packet vs time, with dropped samples shaded red
%   bottom: |CSI| (L-LTF channel estimate) of each packet per antenna
% No payloads are decoded, which keeps this quick; use wifi_detect_plot.m
% for the full per-packet diagnosis.
%
% iq is [nSamples x numRx], int16 (as captureB210 returns) or double.
% stats is a small struct that wifi_capture_raw.m stores with the capture.
if nargin < 3
    titleText = '';
end
if sampleRate ~= 20e6
    error('wifiCaptureCheck expects a 20 MS/s capture, got %g MS/s.', sampleRate / 1e6);
end
if isinteger(iq)
    x = double(iq) / 32768;
else
    x = iq;
end
[nSamples, numRx] = size(x);

% Capture integrity. Exact-zero rows are samples the B210 dropped (USB
% overflow); |I| or |Q| at full scale is ADC clipping (gain too high).
dropped = all(x == 0, 2);
nClipped = sum(any(abs(real(x)) >= 0.999 | abs(imag(x)) >= 0.999, 2));

minSnrDb = 5;
dataBitsPerSymbol = [24 36 48 72 96 144 192 216];
rateCodes = [13 15 5 7 9 11 1 3];

pktTimeMs = zeros(0, 1);
snrDb     = zeros(0, 1);
snrRxDb   = zeros(0, numRx);
csiMag    = zeros(52, numRx, 0);

% Energy gate: wlanPacketDetect is the slow step, so only run it where the
% power in a 4 us block is at least 6 dB above the noise floor (taken as
% the 10th-percentile block power). That matches the minSnrDb cut below.
blockLen = 80;
nBlocks = floor(nSamples / blockLen);
blockPow = sum(reshape(abs(x(1:nBlocks * blockLen, :)).^2, blockLen, nBlocks, numRx), [1 3]);
sortedPow = sort(blockPow(blockPow > 0));
if isempty(sortedPow)
    active = false(1, nBlocks);
else
    active = blockPow > 4 * sortedPow(max(1, round(0.1 * end)));
end

% Search in short windows (wlanPacketDetect correlates its whole input)
% that overlap by more than one preamble, so no packet is missed.
searchLen = 8000;
scanStart = 1;
while scanStart + 399 <= nSamples
    % Jump to one block before the next block with signal in it.
    blk = ceil(scanStart / blockLen);
    nextActive = find(active(min(blk, nBlocks):end), 1);
    if isempty(nextActive)
        break;
    end
    scanStart = max(scanStart, (blk + nextActive - 3) * blockLen + 1);
    if scanStart + 399 > nSamples
        break;
    end
    seg = x(scanStart:min(nSamples, scanStart + searchLen - 1), :);
    offset = wlanPacketDetect(seg, 'CBW20', 0, 0.5);
    if isempty(offset)
        scanStart = scanStart + searchLen - 400;
        continue;
    end
    pktStart = scanStart + double(offset);
    scanStart = pktStart + 160;       % default: step past this L-STF
    if pktStart + 400 > nSamples
        break;
    end

    % Coarse CFO and timing on the L-STF, fine CFO and channel on the L-LTF
    preamble = x(pktStart:pktStart + 400, :);
    cfoCoarse = wlanCoarseCFOEstimate(preamble(1:160, :), 'CBW20');
    syncStart = pktStart + wlanSymbolTimingEstimate( ...
        freqShift(preamble, -cfoCoarse, sampleRate, pktStart), 'CBW20');
    if syncStart < 1 || syncStart + 639 > nSamples
        continue;
    end
    lltfStart = syncStart + 160;
    lltf = x(lltfStart:lltfStart + 159, :);
    cfo = cfoCoarse + wlanFineCFOEstimate( ...
        freqShift(lltf, -cfoCoarse, sampleRate, lltfStart), 'CBW20');
    lltfDemod = wlanLLTFDemodulate(freqShift(lltf, -cfo, sampleRate, lltfStart), 'CBW20');
    chanEst = wlanLLTFChannelEstimate(lltfDemod, 'CBW20');
    noiseVar = wlanLLTFNoiseEstimate(lltfDemod);
    snrPerRx = reshape(mean(abs(chanEst).^2, 1), 1, numRx) / (noiseVar + eps);
    if 10 * log10(sum(snrPerRx)) < minSnrDb
        continue;
    end

    % L-SIG header: real packets pass parity and carry a valid rate and a
    % plausible length; a "greenfield" format verdict means the header was
    % noise (greenfield packets have no L-SIG).
    lsigStart = syncStart + 320;
    [lsigBits, lsigFail] = wlanLSIGRecover(freqShift( ...
        x(lsigStart:lsigStart + 79, :), -cfo, sampleRate, lsigStart), ...
        chanEst, noiseVar, 'CBW20');
    if lsigFail
        continue;
    end
    rateIdx = find(rateCodes == sum(double(lsigBits(1:4)) .* (2.^(3:-1:0))'), 1);
    psduLength = sum(double(lsigBits(6:17)) .* (2.^(0:11))');
    if isempty(rateIdx) || psduLength < 14
        continue;
    end
    pktFormat = wlanFormatDetect(freqShift(x(lsigStart:lsigStart + 239, :), ...
        -cfo, sampleRate, lsigStart), chanEst, noiseVar, 'CBW20');
    if strcmp(pktFormat, 'HT-GF')
        continue;
    end

    pktTimeMs(end + 1, 1) = (syncStart - 1) / sampleRate * 1e3; %#ok<AGROW>
    snrDb(end + 1, 1) = 10 * log10(sum(snrPerRx)); %#ok<AGROW>
    snrRxDb(end + 1, :) = 10 * log10(snrPerRx); %#ok<AGROW>
    csiMag(:, :, end + 1) = reshape(abs(chanEst), 52, numRx); %#ok<AGROW>
    % L-SIG length/rate give the packet's duration for every format (HT/VHT
    % use it as legacy spoofing), so skip straight past the packet.
    nSym = ceil((16 + 8 * psduLength + 6) / dataBitsPerSymbol(rateIdx));
    scanStart = syncStart + 400 + 80 * nSym;
end

%% Verdict
nPackets = numel(snrDb);
problems = {};
if mean(dropped) > 0.01
    problems{end + 1} = sprintf('%.0f%% of samples dropped (USB overflow)', 100 * mean(dropped));
end
if nClipped > 0
    problems{end + 1} = sprintf('%d samples clipped: lower rxGain', nClipped);
end
if nPackets == 0
    problems{end + 1} = 'no Wi-Fi packets found';
elseif median(snrDb) < 10
    problems{end + 1} = sprintf('weak signal: median SNR %.1f dB', median(snrDb));
end
if isempty(problems)
    verdict = 'OK';
else
    verdict = ['CHECK: ' strjoin(problems, '; ')];
end

stats.num_packets      = nPackets;
stats.median_snr_db    = median(snrDb);           % NaN when no packets
stats.median_snr_rx_db = median(snrRxDb, 1);      % per antenna
stats.packet_rate_hz   = nPackets / (nSamples / sampleRate);
stats.dropped_fraction = mean(dropped);
stats.num_clipped      = nClipped;
stats.ok               = isempty(problems);
stats.verdict          = verdict;

fprintf('Check: %d packets (%.0f /s), median SNR %.1f dB, %.1f%% dropped, %d clipped -> %s\n', ...
    nPackets, stats.packet_rate_hz, stats.median_snr_db, ...
    100 * stats.dropped_fraction, nClipped, verdict);

%% Plot (one window, reused on every run)
fig = figure(9701);
clf(fig);
set(fig, 'Color', 'w', 'Name', 'Wi-Fi capture check');
if stats.ok
    verdictColor = [0 0.5 0.15];
else
    verdictColor = [0.8 0.1 0.1];
end

ax1 = subplot(2, 1, 1, 'Parent', fig);
hold(ax1, 'on'); grid(ax1, 'on');
captureMs = nSamples / sampleRate * 1e3;
yTop = max([30; snrDb + 5]);
if any(dropped)
    % One rectangle per run of dropped samples (plotting all samples is slow).
    edges = diff([0; dropped; 0]);
    runMs = ([find(edges == 1), find(edges == -1)] - 1) / sampleRate * 1e3;
    patch(ax1, runMs(:, [1 1 2 2])', repmat([0; yTop; yTop; 0], 1, size(runMs, 1)), ...
        [1 0.8 0.8], 'EdgeColor', 'none', 'DisplayName', 'Dropped samples');
end
rxMarkers = {'o', 's'};
for r = 1:numRx
    plot(ax1, pktTimeMs, snrRxDb(:, r), rxMarkers{min(r, end)}, 'MarkerSize', 4, ...
        'DisplayName', sprintf('RX %d', r));
end
if numRx > 1
    plot(ax1, pktTimeMs, snrDb, 'k.', 'MarkerSize', 10, 'DisplayName', 'Combined');
end
yline(ax1, 10, ':', 'HandleVisibility', 'off');
xlim(ax1, [0 captureMs]);
ylim(ax1, [0 yTop]);
xlabel(ax1, 'Time (ms)');
ylabel(ax1, 'Preamble SNR (dB)');
legend(ax1, 'Location', 'eastoutside');
title(ax1, {titleText, verdict}, 'Color', verdictColor, 'Interpreter', 'none');

ax2 = subplot(2, 1, 2, 'Parent', fig);
hold(ax2, 'on'); grid(ax2, 'on');
subcarriers = [-26:-1, 1:26]';
rxColors = [0 0.35 0.8; 0.85 0.35 0];
for r = 1:numRx
    if nPackets == 0
        break;
    end
    magDb = 20 * log10(reshape(csiMag(:, r, :), 52, nPackets) + eps);
    plot(ax2, subcarriers, magDb, 'Color', [rxColors(r, :) 0.15], ...
        'HandleVisibility', 'off');
    plot(ax2, subcarriers, median(magDb, 2), 'Color', rxColors(r, :), ...
        'LineWidth', 2, 'DisplayName', sprintf('RX %d median', r));
end
xlim(ax2, [-27 27]);
xlabel(ax2, 'Subcarrier index');
ylabel(ax2, '|CSI| (dB)');
title(ax2, sprintf('L-LTF channel estimate, %d packets', nPackets));
if nPackets > 0
    legend(ax2, 'Location', 'eastoutside');
end
drawnow;
end

function out = freqShift(in, fOffset, fs, startSample)
t = (startSample - 1 + (0:size(in, 1) - 1).') / fs;
out = in .* exp(1j * 2 * pi * fOffset * t);
end
