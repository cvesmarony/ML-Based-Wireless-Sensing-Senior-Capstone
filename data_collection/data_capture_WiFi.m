clear;

%% Configuration
sampleRate   = 20e6;
CaptureTime  = 1;                       % 100 ms = 2,000,000 samples (lightweight, won't lag)
TotalSamples = round(sampleRate * CaptureTime);

channel    = 6;
centerFreq = 2407e6 + 5e6 * channel;

rxObj = comm.SDRuReceiver( ...
    'Platform',           'B210', ...
    'SerialNum',          '34C78FD', ...
    'MasterClockRate',    20e6, ...
    'DecimationFactor',   1, ...
    'SamplesPerFrame',    TotalSamples, ... % Pulls entire burst via step call
    'CenterFrequency',    centerFreq, ...
    'ChannelMapping',     1, ...
    'ReceiveAntennaPort', 'RX2', ...
    'Gain',               45, ...
    'OutputDataType',     'double');

disp('Capturing 20 MHz Wi-Fi frame...');
[rxBuffer, ~, overrun] = rxObj();

release(rxObj);

% Plot envelope
t = (0:length(rxBuffer)-1) / sampleRate * 1e3;
plot(t, abs(rxBuffer));
grid on;
title('20 MHz Wi-Fi Packet Bursts');
xlabel('Time (ms)');
ylabel('Envelope Amplitude');

sig = rxBuffer(:, 1);

% 1. Detect start of Wi-Fi Packet
pktStart = wlanPacketDetect(sig, 'CBW20');

if ~empty(pktStart)
    fprintf('Packet detected at sample index: %d\n', pktStart);

    % Extract the Non-HT Long Training Field (L-LTF)
    % L-LTF starts 160 samples after the start of packet (after L-STF)
    lltf_start = pktStart + 160;
    lltf_samples = sig(lltf_start : lltf_start + 160 - 1);

    % Demodulate LLTF and compute Channel Frequency Response (CSI)
    lltf_demod = wlanLLTFDemodulate(lltf_samples, 'CBW20');
    [chanEst, noiseVarEst] = wlanLLTFChannelEstimate(lltf_demod, 'CBW20');

    % chanEst is a 52x1 complex vector representing H[k] for each subcarrier!
    figure;
    subplot(2,1,1);
    stem(-26:26, [abs(chanEst(1:26)); 0; abs(chanEst(27:52))]); % Subcarrier magnitudes
    title('Wi-Fi Channel Magnitude (|H[k]|) per Subcarrier');
    xlabel('Subcarrier Index'); ylabel('|H|'); grid on;

    subplot(2,1,2);
    stem(-26:26, [angle(chanEst(1:26)); 0; angle(chanEst(27:52))]);
    title('Wi-Fi Channel Phase (\angle H[k])');
    xlabel('Subcarrier Index'); ylabel('Phase (rad)'); grid on;
else
    disp('No packet detected in this capture window.');
end