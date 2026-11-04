files = dir(fullfile('wifi_raw', 'raw_ch6_0_0_0_20261104_173049.mat'));
[~, idx] = max([files.datenum]);
S = load(fullfile(files(idx).folder, files(idx).name), 'raw_wifi');
rxBuffer = S.raw_wifi.iq;

figure;
plot(abs(rxBuffer(:, 1)));
title('Antenna 1 IQ Magnitude');
xlabel('Sample');
ylabel('|IQ|');
grid on;