clear;

%% Configuration
sampleRate   = 20e6;
CaptureTime  = 0.1;                       % 100 ms = 2,000,000 samples (lightweight, won't lag)
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