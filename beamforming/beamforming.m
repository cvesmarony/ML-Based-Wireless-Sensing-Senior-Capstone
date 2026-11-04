clear all;
%% ==============================================================================
% NeRF2 TRANSMITTER - BEAMFORMING
% ==============================================================================
%% BEAMFORMING CONFIGURATION
beam_config = struct();
beam_config.fc = 5.4e9; % Center frequency (Hz)
sps = 5;

%% 1. HARDWARE CONFIG
tx = comm.SDRuTransmitter(...
    'Platform', 'B210', ...
    'SerialNum', '3839244', ...
    'ChannelMapping', [1 2], ...
    'MasterClockRate', 30e6, ...
    'InterpolationFactor', 30, ...
    'CenterFrequency', beam_config.fc, ...
    'Gain', [80 80], ...
    'TransportDataType', 'int16');

Fs = 30e6 / 30; % 3 Msps

%% 2. SIGNAL GENERATION - STREAM S

pn_S0 = comm.PNSequence(...
    'Polynomial', [10 3 0], ...
    'InitialConditions', ones(1,10), ...
    'SamplesPerFrame', 1032);

sounding_bits_S0 = pn_S0();
sounding_bits_S90 = sounding_bits_S0;
for ii = 1:1032/4
    sounding_bits_S90(4*ii-3) = - (sounding_bits_S0(4*ii-1) - 1);
    sounding_bits_S90(4*ii-2) = sounding_bits_S0(4*ii);
    sounding_bits_S90(4*ii-1) = sounding_bits_S0(4*ii-3);
    sounding_bits_S90(4*ii) = - (sounding_bits_S0(4*ii-2) - 1);
end

% QPSK Modulate
qpskMod = comm.QPSKModulator('BitInput',true, 'PhaseOffset', pi/4);

sounding_syms_S0 = qpskMod(sounding_bits_S0);
sounding_syms_S90 = qpskMod(sounding_bits_S90);

% RRC Pulse Shaping
rrcTx = comm.RaisedCosineTransmitFilter(...
    'RolloffFactor', 0.25, ...
    'FilterSpanInSymbols', 6, ...
    'OutputSamplesPerSymbol', sps);

waveform_S0 = rrcTx(sounding_syms_S0);

rrcTx = comm.RaisedCosineTransmitFilter(...
    'RolloffFactor', 0.25, ...
    'FilterSpanInSymbols', 6, ...
    'OutputSamplesPerSymbol', sps);

waveform_S90 = rrcTx(sounding_syms_S90);

waveform_S0 = waveform_S0(3*sps + 1 : 3*sps + 2560);
waveform_S90 = waveform_S90(3*sps + 1 : 3*sps + 2560);


%% 3. NORMALIZE STREAMS
waveform_S0 = waveform_S0 / max(abs(waveform_S0)) * 0.8;
waveform_S90 = waveform_S90 / max(abs(waveform_S90)) * 0.8;

bits_S0 = sounding_bits_S0(1:1024);
bits_S90 = sounding_bits_S90(1:1024);

syms_S0 = sounding_syms_S0(1:512);
syms_S90 = sounding_syms_S90(1:512);

% save('waveform_STTD.mat', 'waveform_S1', 'waveform_S2', 'bits_S1', 'bits_S2', 'syms_S1', 'syms_S2');

%% 4. TWO-ANTENNA BEAMFORMING
% Creates a 0 deg beam and 90 deg beam using the same initial signal S

% Apply equal-power complex weights to form the initial coherent streams
fc = 5.4e9;
c = 299792458;
lambda = c/fc;

% Antenna array parameters
d = lambda/2;
theta1 = deg2rad(0);
theta2 = deg2rad(90);

% Weight each signal
beamPhase0 = exp(1j * ((2*pi * d) / lambda) * cos(theta1));
S0_weight1 = 1 / sqrt(2);
S0_weight2 = beamPhase0 / sqrt(2);
waveform_S0w1 = S0_weight1 * waveform_S0;
waveform_S0w2 = S0_weight2 * waveform_S0;

beamPhase90 = exp(1j * ((2*pi * d) / lambda) * cos(theta2));
S90_weight1 = 1 / sqrt(2);
S90_weight2 = beamPhase90 / sqrt(2);
waveform_S90w1 = S90_weight1 * waveform_S90;
waveform_S90w2 = S90_weight2 * waveform_S90;

% framing
silence = complex(zeros(5000,1));
silence2 = complex(zeros(7560,1));
tx1_0_frame = [silence; waveform_S0w1(:); silence2];
tx1_90_frame = [silence; waveform_S90w1(:); silence2];
tx2_0_frame = [silence; waveform_S0w2(:); silence2];
tx2_90_frame = [silence; waveform_S90w2(:); silence2];

% frames to be transmitted: 10 frames of 0 deg beam, 10 frames of 90 deg
% beam, repeat
tx1_frame = [repmat(tx1_0_frame, 20, 1); repmat(tx1_90_frame, 20, 1)];
tx2_frame = [repmat(tx2_0_frame, 20, 1); repmat(tx2_90_frame, 20, 1)];

% Generate simultaneous antenna signals
tx1 = tx1_frame;
tx2 = tx2_frame;

% Transmit continuously through both channels
disp("Beamforming transmission active. Press Ctrl+C to stop.");

while true
    tx([tx1, tx2]);
end
