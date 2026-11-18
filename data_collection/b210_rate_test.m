% Finds which sample rates / antenna counts / read modes this laptop can
% stream from the B210 without dropping samples. Each row is a 100 ms
% capture. Row 1 is exactly data_capture.m's setup (known good) as the
% control. "zeros" is the fraction of the buffer that never arrived; the
% overrun flag alone can be raised by a single hiccup.
clear;
sdrSerial  = '3839223';         % check with "findsdru"
centerFreq = 2437e6;
rxGain     = 45;
captureSec = 0.1;

%          MasterClockRate  Decimation  numRx  mode
tests = {   30e6,             60,         2,    'single'   % = data_capture.m
            20e6,             4,          1,    'single'   %  5 MS/s
            20e6,             2,          1,    'single'   % 10 MS/s
            20e6,             1,          1,    'single'   % 20 MS/s, one big frame
            20e6,             1,          1,    'burst'    % 20 MS/s, captureB210
            20e6,             1,          2,    'single'
            20e6,             1,          2,    'burst'};

fprintf('%-8s %-6s %-7s %-8s %-8s %-8s %s\n', 'MS/s', 'numRx', 'mode', 'MB/s', ...
    'overrun', 'zeros', 'error');
for k = 1:size(tests, 1)
    [mcr, dec, numRx, mode] = tests{k, :};
    fs = mcr / dec;
    overrun = NaN; zeroFrac = NaN; errText = '';
    try
        if strcmp(mode, 'burst')
            [iq, overrun] = captureB210(sdrSerial, centerFreq, fs, rxGain, ...
                'RX2', numRx, captureSec);
        else
            rx = comm.SDRuReceiver('Platform', 'B210', 'SerialNum', sdrSerial, ...
                'MasterClockRate', mcr, 'DecimationFactor', dec, ...
                'SamplesPerFrame', round(fs * captureSec), ...
                'CenterFrequency', centerFreq, 'ChannelMapping', 1:numRx, ...
                'Gain', repmat(rxGain, 1, numRx), 'OutputDataType', 'int16');
            len = 0;
            while len == 0      % first call returns empty while the radio starts
                [iq, len, overrun] = rx();
            end
            release(rx);
        end
        zeroFrac = mean(all(iq == 0, 2));
    catch err
        errText = strtok(err.message, newline);
    end
    fprintf('%-8g %-6d %-7s %-8g %-8d %-7.1f%% %s\n', fs / 1e6, numRx, mode, ...
        4 * numRx * fs / 1e6, overrun, 100 * zeroFrac, errText);
    pause(1);
end
