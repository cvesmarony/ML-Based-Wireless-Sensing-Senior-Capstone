clear;
% Raw 2-antenna Wi-Fi I/Q capture from the B210. No processing: the buffer
% is saved as-is. Decode and plot it afterwards with data_capture_WiFi.m.

%% Capture metadata - adjust before every capture
locationName = 'unspecified';   % e.g. 'Halligan 102'
x_lab        = 0;               % RX position in the room grid
y_lab        = 0;
rx_height    = 0.571;           % m
orientation  = 0;               % RX orientation (deg)
captureNotes = '';              % free text: people present, doors open, etc.
data_folder  = 'wifi_raw/';

%% Configuration
sampleRate   = 20e6;
CaptureTime  = 0.1;            % s (~3.2 MB per 10 ms per antenna)
channel      = 6;
centerFreq   = 2407e6 + 5e6 * channel;
rxGain       = 45;              % dB, same on both channels
sdrSerial    = '3839223';       % check with "findsdru"
antennaPort  = 'RX2';
numRx        = 2;               % RF A and RF B

%% Capture
rxObj = comm.SDRuReceiver( ...
    'Platform',           'B210', ...
    'SerialNum',          sdrSerial, ...
    'MasterClockRate',    sampleRate, ...
    'DecimationFactor',   1, ...
    'SamplesPerFrame',    round(sampleRate * CaptureTime), ...
    'CenterFrequency',    centerFreq, ...
    'ChannelMapping',     1:numRx, ...
    'ReceiveAntennaPort', antennaPort, ...
    'Gain',               repmat(rxGain, 1, numRx), ...
    'OutputDataType',     'double');

fprintf('Capturing %g ms on %d antennas, channel %d (%.0f MHz)...\n', ...
    CaptureTime * 1e3, numRx, channel, centerFreq / 1e6);
% Host clock at the start of the capture call; USB latency makes this
% accurate to a few ms, not to the sample.
captureStartTime = datetime('now', 'TimeZone', 'local');
[rxBuffer, ~, overrun] = rxObj();
release(rxObj);

if overrun
    warning('Buffer overrun detected. Some samples may be dropped.');
end

%% Save
% Times and text are stored as char/double (not datetime/string objects) so
% the file also loads in Python (scipy.io.loadmat).
raw_wifi.iq                 = rxBuffer;   % [nSamples x numRx] complex double
raw_wifi.Fs                 = sampleRate;
raw_wifi.CaptureTime        = CaptureTime;
raw_wifi.channel            = channel;
raw_wifi.freq               = centerFreq;
raw_wifi.gain_db            = rxGain;
raw_wifi.sdr_serial         = sdrSerial;
raw_wifi.antenna_port       = antennaPort;
raw_wifi.num_rx             = numRx;
raw_wifi.overrun            = logical(overrun);
raw_wifi.capture_time_iso   = char(captureStartTime, ...
    'yyyy-MM-dd''T''HH:mm:ss.SSSSSSXXX');
raw_wifi.capture_time_posix = posixtime(captureStartTime);
raw_wifi.location           = locationName;
raw_wifi.rx_pos             = [x_lab, y_lab, rx_height];
raw_wifi.rx_orient          = orientation;
raw_wifi.notes              = captureNotes;
raw_wifi.matlab_version     = version;

if ~isfolder(data_folder)
    mkdir(data_folder);
end
filename = sprintf('%sraw_ch%d_%g_%g_%g_%s.mat', data_folder, channel, ...
    x_lab, y_lab, orientation, char(captureStartTime, 'yyyyMMdd_HHmmss'));
if isfile(filename)
    error('duplicate file name: %s', filename);
end
save(filename, 'raw_wifi');
fprintf('Saved %d x %d samples to %s\n', size(rxBuffer, 1), size(rxBuffer, 2), filename);
