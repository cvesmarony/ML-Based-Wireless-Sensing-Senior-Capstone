%% DATA CAPTURE
% Requires initialize_rx.m to have already been run.

if ~exist('rxObj', 'var')
    error('Run initialize_rx.m first.');
end

% Adjust to match RX position
x_lab = 7; % 0 to 15
y_lab = 11; % 0 to 11
orientation = 0; % 0, 45, 90, 135

disp(string(x_lab) + " " + string(y_lab) + " " + string(orientation));

% Adjust to match TX position, ALSO UPDATE THE FILENAME TO TX#

tx_num = 1;
tx_pos = [11, 15, 0];

% tx_num = 2;
% tx_pos = [1, 8, 0];

% tx_num = 3;
% tx_pos = [8, 5, 0];

if x_lab == round(x_lab) && y_lab == round(y_lab)
    data_folder = "data_f2/";
elseif x_lab == round(x_lab) || y_lab == round(y_lab)
    disp("Something is mismatched check x y")
else
    data_folder = "test_data_f2/";
end

%% -------------------------------------------------------------------------
% CONFIGURATION
% -------------------------------------------------------------------------
Fs          = 1e6;
CaptureTime = 2;           % seconds
TotalSamples = round(Fs * CaptureTime);

% Change decimation factor -> change Fs
% masterclockrate/decimationfactor = sampleperframe

%% -------------------------------------------------------------------------
% CAPTURE
% -------------------------------------------------------------------------
disp(['Capturing ', num2str(CaptureTime), ' seconds of data...']);

% [rxBuffer, ~, overrun] = rxObj();

[rxBuffer, len, overrun] = rxObj();

disp(['Samples returned: ', num2str(len)]);
disp(['Buffer size: ', num2str(size(rxBuffer,1))]);
disp(['Overrun: ', num2str(overrun)]);
if overrun
    warning('Overrun detected during capture.');
end
if size(rxBuffer,1) < TotalSamples
    warning('Incomplete capture.');
end

disp('Capture complete. Saving raw IQ buffer...');

%% -------------------------------------------------------------------------
% SAVE
% -------------------------------------------------------------------------
raw_data.iq          = rxBuffer;        % [TotalSamples x 2] complex double
raw_data.Fs          = Fs;
raw_data.CaptureTime = CaptureTime;
raw_data.freq        = 5.4e9;
raw_data.rx_pos      = [x_lab, y_lab, 0.571];
raw_data.tx_pos      = tx_pos;
raw_data.mimo        = false;
raw_data.rx_orient   = orientation;
    
if x_lab >= 0
    filename = data_folder + 'beam_tx' + string(tx_num) + '_' + string(x_lab) + '_' + string(y_lab) + '_' + string(orientation) + '.mat';
else
    filename = data_folder + 'beam_tx' + string(tx_num) + '_n1_' + string(y_lab) + '_' + string(orientation) + '.mat';
end

% if you f-up a file then just uncomment this and ignore the error
% generated
%save(filename, 'raw_data');
%disp('Fixer File saved!');

if not(isfile(filename))
    save(filename, 'raw_data'); % dev filename
    disp('File saved!');
else

    error("duplicate file name")
end

%% MUST RELEASE RX IN CLI AFTER
% release(rxObj);
% clear rxObj;