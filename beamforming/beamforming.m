clear all;
%% ==============================================================================
% NeRF2 TRANSMITTER - MIMO MODE (turn-taking antennas)
% ==============================================================================
%% MIMO CONFIGURATION
mimo_config = struct();
mimo_config.fc = 5.4e9; % Center frequency (Hz)
sps = 5;

%% 1. HARDWARE CONFIG
tx = comm.SDRuTransmitter(...
    'Platform', 'B210', ...
    'SerialNum', '34C78EF', ...
    'ChannelMapping', [1 2], ...
    'MasterClockRate', 30e6, ...
    'InterpolationFactor', 60, ...
    'CenterFrequency', mimo_config.fc, ...
    'Gain', [80 80], ...
    'TransportDataType', 'int16');

Fs = 30e6 / 60; % 3 Msps

%% 2. SIGNAL GENERATION - STREAM S

pn_S1 = comm.PNSequence(...
    'Polynomial', [10 3 0], ...
    'InitialConditions', ones(1,10), ...
    'SamplesPerFrame', 1032);

sounding_bits_S1 = pn_S1();

sounding_bits_S2 = sounding_bits_S1;
for ii = 1:1032/4
    sounding_bits_S2(4*ii-3) = - (sounding_bits_S1(4*ii-1) - 1);
    sounding_bits_S2(4*ii-2) = sounding_bits_S1(4*ii);
    sounding_bits_S2(4*ii-1) = sounding_bits_S1(4*ii-3);
    sounding_bits_S2(4*ii) = - (sounding_bits_S1(4*ii-2) - 1);
end

% QPSK Modulate
qpskMod = comm.QPSKModulator('BitInput',true, 'PhaseOffset', pi/4);

sounding_syms_S1 = qpskMod(sounding_bits_S1);
sounding_syms_S2 = qpskMod(sounding_bits_S2);

% RRC Pulse Shaping
rrcTx = comm.RaisedCosineTransmitFilter(...
    'RolloffFactor', 0.25, ...
    'FilterSpanInSymbols', 6, ...
    'OutputSamplesPerSymbol', sps);

waveform_S1 = rrcTx(sounding_syms_S1);

rrcTx = comm.RaisedCosineTransmitFilter(...
    'RolloffFactor', 0.25, ...
    'FilterSpanInSymbols', 6, ...
    'OutputSamplesPerSymbol', sps);

waveform_S2 = rrcTx(sounding_syms_S2);

waveform_S1 = waveform_S1(3*sps + 1 : 3*sps + 2560);
waveform_S2 = waveform_S2(3*sps + 1 : 3*sps + 2560);

%% 3. NORMALIZE BOTH STREAMS
waveform_S1 = waveform_S1 / max(abs(waveform_S1)) * 0.8;
waveform_S2 = waveform_S2 / max(abs(waveform_S2)) * 0.8;

bits_S1 = sounding_bits_S1(1:1024);
bits_S2 = sounding_bits_S2(1:1024);

syms_S1 = sounding_syms_S1(1:512);
syms_S2 = sounding_syms_S2(1:512);

save('waveform_STTD.mat', 'waveform_S1', 'waveform_S2', 'bits_S1', 'bits_S2', 'syms_S1', 'syms_S2');


% disp("Alternating between MIMO Streams Continuously... Press Ctrl+C to stop.");

silence = complex(zeros(5000, 1));
silence2 = complex(zeros(7560, 1));

%% 4. TWO-ANTENNA PHASE-STEERED BEAMFORMING

% Carrier frequency
fc = 5.4e9;
c = 299792458;
lambda = c/fc;

% Antenna array parameters
d = 0.03;           % Half-wavelength antenna spacing
theta_deg = 0;         % Desired beam angle from broadside
theta = deg2rad(theta_deg);

% Steering phase for a two-element uniform linear array
% Sign depends on antenna ordering and angle convention.
phi = -2*pi*d*sin(theta)/lambda;

% Use ONE common waveform for coherent beamforming
s = waveform_S1(:);

% Optional framing: simultaneous transmission, not turn-taking
silence = complex(zeros(5000,1));
silence2 = complex(zeros(7560,1));
s_frame = [silence; s; silence2];

% Equal-amplitude weights, normalized to unit total power
w1 = 1/sqrt(2);
w2 = exp(1j*phi)/sqrt(2);

% Generate simultaneous antenna signals
tx1 = w1*s_frame;
tx2 = w2*s_frame;

% Transmit continuously through both channels
disp("Beamforming transmission active. Press Ctrl+C to stop.");

while true
    tx([tx1, tx2]);
end
