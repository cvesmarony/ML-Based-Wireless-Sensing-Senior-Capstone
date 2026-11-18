function [iq, overrun] = captureB210(sdrSerial, centerFreq, sampleRate, ...
        rxGain, antennaPort, numRx, captureSec)
% Contiguous multi-antenna I/Q capture from a B210. Returns iq as the
% radio's native complex int16, [nSamples x numRx] (full scale = 32768;
% divide by 32768 for doubles in +/-1). int16 is a quarter the memory of
% double and spares the host the conversion while the stream is running.
%
% Requesting the whole capture as one huge frame overruns: the host cannot
% drain the USB stream fast enough, so only the first few ms arrive and the
% rest of the buffer is zeros. Burst mode instead has the driver buffer
% numFrames contiguous 1 ms frames, which are then read out one at a time.
frameLen  = round(sampleRate * captureSec);
numFrames = ceil(captureSec * sampleRate / frameLen);
channel = 6;

rxObj = comm.SDRuReceiver( ...
    'Platform',           'B210', ...
    'SerialNum',          '3839244', ...
    'MasterClockRate',    20e6, ...
    'DecimationFactor',   1, ...
    'SamplesPerFrame',    frameLen, ...
    'CenterFrequency',    2407e6 + 5e6 * channel, ...
    'ChannelMapping',     1:numRx, ...
    'ReceiveAntennaPort', 'RX2', ...
    'Gain',               repmat(rxGain, 1, numRx), ...
    'OutputDataType',     'double', ...
    'EnableBurstMode',    true, ...
    'NumFramesInBurst',   numFrames);
cleanup = onCleanup(@() release(rxObj));

iq = [];        % allocated on the first frame with the radio's output type
overrun = false;
for k = 1:numFrames
    len = 0;
    while len == 0          % the first call returns empty while the radio starts
        try
            [frame, len, frameOverrun] = rxObj();
        catch err
            if contains(err.message, 'Overflow')
                error(['B210 overflowed: the USB link cannot sustain %d antenna(s) ' ...
                    'at %g MS/s (%g MB/s). On this laptop 1 antenna at 20 MS/s ' ...
                    'works and 2 do not; set numRx = 1.'], numRx, sampleRate / 1e6, ...
                    4 * numRx * sampleRate / 1e6);
            end
            rethrow(err);
        end
    end
    if isempty(iq)
        % Same type as the frames (int16, single or double): preallocating
        % int16 and writing +/-1-scaled single/double frames into it rounds
        % every sample to 0, which then looks exactly like dropped samples.
        iq = zeros(numFrames * frameLen, numRx, 'like', frame);
    end
    iq((k - 1) * frameLen + (1:frameLen), :) = frame;
    overrun = overrun || frameOverrun;
end
iq = iq(1:round(captureSec * sampleRate), :);

% Safety net: a stalled stream shows up as long runs of exact zeros.
zeroRows = all(iq == 0, 2);
if mean(zeroRows) > 0.01
    warning(['%.0f%% of the capture is zeros: samples were dropped. Check the ' ...
        'B210 is on a USB 3 port, or lower the sample rate.'], 100 * mean(zeroRows));
end
end
