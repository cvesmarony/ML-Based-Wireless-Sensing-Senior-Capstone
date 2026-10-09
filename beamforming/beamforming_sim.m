clear all;
%% ==============================================================================
% NeRF2 TRANSMITTER - MIMO MODE (turn-taking antennas)
% ==============================================================================
%% MIMO CONFIGURATION
mimo_config = struct();
mimo_config.fc = 915e6; % Center frequency (Hz)
sps = 5;

%% 1. HARDWARE CONFIG
% tx = comm.SDRuTransmitter(...
%     'Platform', 'B210', ...
%     'SerialNum', '34C78EF', ...
%     'ChannelMapping', [1 2], ...
%     'MasterClockRate', 30e6, ...
%     'InterpolationFactor', 60, ...
%     'CenterFrequency', mimo_config.fc, ...
%     'Gain', [60 60], ...
%     'TransportDataType', 'int16');

Fs = 30e6 / 60; % 3 Msps

%% 2. SIGNAL GENERATION - STREAM S

% pn_S1 = comm.PNSequence(...
%     'Polynomial', [10 3 0], ...
%     'InitialConditions', ones(1,10), ...
%     'SamplesPerFrame', 1032);
% 
% sounding_bits_S1 = pn_S1();
% 
% sounding_bits_S2 = sounding_bits_S1;
% for ii = 1:1032/4
%     sounding_bits_S2(4*ii-3) = - (sounding_bits_S1(4*ii-1) - 1);
%     sounding_bits_S2(4*ii-2) = sounding_bits_S1(4*ii);
%     sounding_bits_S2(4*ii-1) = sounding_bits_S1(4*ii-3);
%     sounding_bits_S2(4*ii) = - (sounding_bits_S1(4*ii-2) - 1);
% end

% QPSK Modulate
qpskMod = comm.QPSKModulator('BitInput',true, 'PhaseOffset', pi/4);
% likely need 2 separate qpsk modulating functions for different phases

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

% save('waveform_STTD.mat', 'waveform_S1', 'waveform_S2', 'bits_S1', 'bits_S2', 'syms_S1', 'syms_S2');
% 
% 
% disp("Alternating between MIMO Streams Continuously... Press Ctrl+C to stop.");

silence = complex(zeros(5000, 1));
silence2 = complex(zeros(7560, 1));

% change this to send waveforms simultaneously
% probably don't need frames since constantly sending data
tx_frame_S1 = [silence; waveform_S1; silence2];
tx_frame_S2 = [silence2; silence; waveform_S2];



% %% 5. CONTINUOUS TRANSMISSION LOOP
% while true
%     tx([tx_frame_S1, tx_frame_S2]);
% end


%% 6. PLOT TRANSMIT FRAMES 1 AND 2
% Time vector in milliseconds
t_ms = (0:length(tx_frame_S1)-1) / Fs * 1000;

% Figure 1: Full Frame Turn-Taking Overview
figure('Name', 'MIMO Transmit Frames Overview', 'NumberTitle', 'off');

subplot(2,1,1);
plot(t_ms, abs(tx_frame_S1), 'LineWidth', 1.2, 'DisplayName', 'Frame 1 (Stream S1)');
hold on;
plot(t_ms, abs(tx_frame_S2), '--', 'LineWidth', 1.2, 'DisplayName', 'Frame 2 (Stream S2)');
hold off;
title('NeRF2 MIMO Transmit Frames Amplitude Envelope');
xlabel('Time (ms)');
ylabel('Magnitude |s(t)|');
ylim([-0.05, 1.0]);
grid on;
legend('Location', 'northeast');

subplot(2,1,2);
plot(t_ms, real(tx_frame_S1), 'LineWidth', 1.0, 'DisplayName', 'Frame 1 In-Phase (I)');
hold on;
plot(t_ms, real(tx_frame_S2), 'LineWidth', 1.0, 'DisplayName', 'Frame 2 In-Phase (I)');
hold off;
title('In-Phase (Real) Waveforms over Time');
xlabel('Time (ms)');
ylabel('Amplitude');
grid on;
legend('Location', 'northeast');

% Figure 2: Detailed 3-Panel Breakdown with Zoomed Bursts
figure('Name', 'MIMO Transmit Frames Detailed Breakdown', 'NumberTitle', 'off');

% Panel 1: Magnitude Overview
subplot(3,1,1);
plot(t_ms, abs(tx_frame_S1), 'Color', [0 0.447 0.741], 'LineWidth', 1.2, 'DisplayName', 'Frame 1 (Stream S1)');
hold on;
plot(t_ms, abs(tx_frame_S2), '--', 'Color', [0.85 0.325 0.098], 'LineWidth', 1.2, 'DisplayName', 'Frame 2 (Stream S2)');
hold off;
title('MIMO Turn-Taking Transmission Frames Overview');
xlabel('Time (ms)');
ylabel('Magnitude |s(t)|');
ylim([-0.05, 1.0]);
grid on;
legend('Location', 'northeast');

% Panel 2: Zoomed Burst 1 (10.0 ms to 15.12 ms)
subplot(3,1,2);
plot(t_ms, real(tx_frame_S1), 'DisplayName', 'Frame 1 In-Phase (I)');
hold on;
plot(t_ms, imag(tx_frame_S1), 'DisplayName', 'Frame 1 Quadrature (Q)');
hold off;
xlim([9.5, 15.6]);
title('Zoomed View: Frame 1 Active Sounding Burst (10.0 ms to 15.12 ms)');
xlabel('Time (ms)');
ylabel('Amplitude');
grid on;
legend('Location', 'northeast');

% Panel 3: Zoomed Burst 2 (25.12 ms to 30.24 ms)
subplot(3,1,3);
plot(t_ms, real(tx_frame_S2), 'DisplayName', 'Frame 2 In-Phase (I)');
hold on;
plot(t_ms, imag(tx_frame_S2), 'DisplayName', 'Frame 2 Quadrature (Q)');
hold off;
xlim([24.5, 30.6]);
title('Zoomed View: Frame 2 Active Sounding Burst (25.12 ms to 30.24 ms)');
xlabel('Time (ms)');
ylabel('Amplitude');
grid on;
legend('Location', 'northeast');
