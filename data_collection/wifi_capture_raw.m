% Raw Wi-Fi I/Q capture from the B210, then a quick health check
% (wifiCaptureCheck: packets found, their SNR and CSI, dropped/clipped
% samples) so every capture is verified before moving on. The buffer is
% saved unprocessed as int16; decode it later with data_capture_WiFi.m.
% Edit only the metadata block between captures.

%% Capture metadata - adjust before every capture
x_lab        = 0;               % RX position in the room grid
y_lab        = 0;
orientation  = 0;               % RX orientation (deg)
captureNotes = '';              % free text: people present, doors open, etc.

%% Fixed setup
locationName = 'unspecified';   % e.g. 'Halligan 102'
rx_height    = 0.571;           % m
data_folder  = 'wifi_raw/';
sampleRate   = 20e6;
CaptureTime  = 0.1;             % s (8 MB per 100 ms per antenna as int16)
channel      = 6;
centerFreq   = 2407e6 + 5e6 * channel;
rxGain       = 76;              % dB; the check flags clipping if this is too high
sdrSerial    = '3839223';       % check with "findsdru"
antennaPort  = 'RX2';
numRx        = 1;               % 2 antennas at 20 MS/s overflows this laptop's USB link

%% Capture
captureStartTime = datetime('now', 'TimeZone', 'local');   % host clock, +/- a few ms
[iq, overrun] = captureB210(sdrSerial, centerFreq, sampleRate, rxGain, ...
    antennaPort, numRx, CaptureTime);
fileTag = sprintf('raw_ch%d_%g_%g_%g_%s', channel, x_lab, y_lab, orientation, ...
    char(captureStartTime, 'yyyyMMdd_HHmmss'));
check = wifiCaptureCheck(iq, sampleRate, fileTag);

%% Save
% Plain numeric/char fields only, saved as uncompressed v7.3 (HDF5): fast
% to write, and readable from Python with h5py.
raw_wifi = struct( ...
    'iq',                 iq, ...          % [nSamples x numRx] complex int16
    'Fs',                 sampleRate, ...
    'CaptureTime',        CaptureTime, ...
    'channel',            channel, ...
    'freq',               centerFreq, ...
    'gain_db',            rxGain, ...
    'sdr_serial',         sdrSerial, ...
    'antenna_port',       antennaPort, ...
    'num_rx',             numRx, ...
    'overrun',            logical(overrun), ...
    'capture_time_iso',   char(captureStartTime, 'yyyy-MM-dd''T''HH:mm:ss.SSSSSSXXX'), ...
    'capture_time_posix', posixtime(captureStartTime), ...
    'location',           locationName, ...
    'rx_pos',             [x_lab, y_lab, rx_height], ...
    'rx_orient',          orientation, ...
    'notes',              captureNotes, ...
    'check',              check);

if ~isfolder(data_folder)
    mkdir(data_folder);
end
filename = fullfile(data_folder, [fileTag '.mat']);
if isfile(filename)
    error('duplicate file name: %s', filename);
end
save(filename, 'raw_wifi', '-v7.3', '-nocompression');
fprintf('Saved %s\n', filename);
