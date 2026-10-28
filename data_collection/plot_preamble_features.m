%% Plot preamble-derived localization features from saved captures
% Loads every wifi_data/*.mat written by data_capture_WiFi.m and plots what
% the L-STF / L-LTF give us: CSI amplitude and phase, power delay profile,
% and CFO / power per transmitter. Packets are colored by BSSID.
clear;
data_folder = 'wifi_data/';
files = dir(fullfile(data_folder, '*.mat'));
if isempty(files)
    error('No .mat files in %s', data_folder);
end

%% Gather packets from every capture
csi = zeros(52, 0); cfoHz = []; powerDb = []; tPosix = []; bssid = {};
for f = 1:numel(files)
    s = load(fullfile(files(f).folder, files(f).name), 'wifi_data');
    pk = s.wifi_data.packets;
    for k = 1:numel(pk)
        if numel(pk(k).csi) ~= 52
            continue;
        end
        csi(:, end+1) = pk(k).csi(:);
        cfoHz(end+1) = pk(k).cfo_hz;
        powerDb(end+1) = pk(k).power_dbfs;
        tPosix(end+1) = pk(k).time_posix;
        bssid{end+1} = pk(k).bssid;
    end
end
[tPosix, order] = sort(tPosix);
csi = csi(:, order); cfoHz = cfoHz(order); powerDb = powerDb(order);
bssid = bssid(order);
nPkt = numel(tPosix);
fprintf('%d packets from %d file(s)\n', nPkt, numel(files));

%% Group by BSSID: six most common get a color, the rest fold into "Other"
% Colors are assigned in alphabetical BSSID order so a BSSID keeps its color.
palette = [0x2a 0x78 0xd6; 0xeb 0x68 0x34; 0x1b 0xaf 0x7a; ...
           0xed 0xa1 0x00; 0xe8 0x7b 0xa4; 0x00 0x83 0x00] / 255;
otherColor = [0.6 0.6 0.6];
[ids, ~, g] = unique(bssid);
counts = accumarray(g(:), 1);
[~, byCount] = sort(counts, 'descend');
kept = sort(byCount(1:min(6, numel(ids))));
groupNames = [ids(kept), {'Other'}];
groupColors = [palette(1:numel(kept), :); otherColor];
grp = (numel(kept) + 1) * ones(1, nPkt);
for i = 1:numel(kept)
    grp(g == kept(i)) = i;
end
usedGroups = unique(grp);

%% Derived features
sc = [-26:-1, 1:26].';                    % L-LTF data/pilot subcarrier indices
ampDb = 20 * log10(abs(csi) + eps);

% Phase: unwrap, then remove the linear term (timing offset) and constant
% term (residual CFO / phase offset) so what is left is the channel shape.
phase = unwrap(angle(csi));
for p = 1:nPkt
    coeffs = polyfit(sc, phase(:, p), 1);
    phase(:, p) = phase(:, p) - polyval(coeffs, sc);
end

% Power delay profile: IFFT of the CSI, zero-padded to 512 points (same
% 312.5 kHz spacing -> 3.2 us span, 6.25 ns bins). Absolute delay is
% unknown with one unsynchronized receiver, so each profile is aligned to
% its strongest path at 0 ns.
nfft = 512;
dt = 3.2e-6 / nfft;
leadSamples = round(200e-9 / dt);
pdpDb = zeros(nfft, nPkt);
for p = 1:nPkt
    H = zeros(nfft, 1);
    H(mod(sc, nfft) + 1) = csi(:, p);
    h = ifft(H);
    [~, peak] = max(abs(h));
    h = circshift(h, leadSamples - (peak - 1));
    pdpDb(:, p) = 20 * log10(abs(h) / max(abs(h)) + eps);
end
tauNs = ((0:nfft-1) - leadSamples) * dt * 1e9;

%% Plot
fig = figure('Name', 'Preamble features', 'Color', 'w', ...
    'Position', [100 100 1400 750]);
tl = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('Preamble features - %d packets', nPkt));
lineAlpha = 0.35;

% 1. CSI amplitude
ax1 = nexttile(tl); hold(ax1, 'on');
for p = 1:nPkt
    plot(ax1, sc, ampDb(:, p), 'Color', [groupColors(grp(p), :) lineAlpha], ...
        'LineWidth', 1);
end
xlabel(ax1, 'Subcarrier index'); ylabel(ax1, '|H| (dB)');
title(ax1, 'CSI amplitude (L-LTF)');

% 2. Sanitized CSI phase
ax2 = nexttile(tl); hold(ax2, 'on');
for p = 1:nPkt
    plot(ax2, sc, phase(:, p), 'Color', [groupColors(grp(p), :) lineAlpha], ...
        'LineWidth', 1);
end
xlabel(ax2, 'Subcarrier index'); ylabel(ax2, 'Phase (rad)');
title(ax2, 'CSI phase, linear trend removed');

% 3. Power delay profile
ax3 = nexttile(tl); hold(ax3, 'on');
for p = 1:nPkt
    plot(ax3, tauNs, pdpDb(:, p), 'Color', [groupColors(grp(p), :) lineAlpha], ...
        'LineWidth', 1);
end
xlim(ax3, [-200 1000]); ylim(ax3, [-30 1]);
xlabel(ax3, 'Delay relative to strongest path (ns)');
ylabel(ax3, 'Power (dB, normalized)');
title(ax3, 'Power delay profile');

% 4. CFO over time: one transmitter's oscillator holds a near-constant offset
ax4 = nexttile(tl); hold(ax4, 'on');
tRel = tPosix - tPosix(1);
for gi = usedGroups
    m = grp == gi;
    scatter(ax4, tRel(m), cfoHz(m) / 1e3, 36, groupColors(gi, :), 'filled', ...
        'MarkerEdgeColor', 'w', 'DisplayName', groupNames{gi});
end
xlabel(ax4, 'Time since first packet (s)'); ylabel(ax4, 'CFO (kHz)');
title(ax4, 'CFO over time');

% 5. CFO vs power: clusters ~ individual transmitters
ax5 = nexttile(tl); hold(ax5, 'on');
for gi = usedGroups
    m = grp == gi;
    scatter(ax5, powerDb(m), cfoHz(m) / 1e3, 36, groupColors(gi, :), ...
        'filled', 'MarkerEdgeColor', 'w');
end
xlabel(ax5, 'Packet power (dBFS)'); ylabel(ax5, 'CFO (kHz)');
title(ax5, 'CFO vs power (transmitter fingerprint)');

% 6. CSI amplitude heatmap, packets in time order
ax6 = nexttile(tl);
imagesc(ax6, sc, 1:nPkt, ampDb.');
blueRamp = [0xcd 0xe2 0xfb; 0x9e 0xc5 0xf4; 0x6d 0xa7 0xec; 0x39 0x87 0xe5; ...
            0x25 0x6a 0xbf; 0x18 0x4f 0x95; 0x0d 0x36 0x6b] / 255;
colormap(ax6, interp1(linspace(0, 1, 7), blueRamp, linspace(0, 1, 256)));
cb = colorbar(ax6); cb.Label.String = '|H| (dB)';
xlabel(ax6, 'Subcarrier index'); ylabel(ax6, 'Packet (time order)');
title(ax6, 'CSI amplitude, all packets');

for ax = [ax1 ax2 ax3 ax4 ax5]
    grid(ax, 'on'); ax.GridAlpha = 0.15; box(ax, 'off');
end

lgd = legend(ax4, 'Orientation', 'horizontal', 'Interpreter', 'none');
lgd.Layout.Tile = 'south';
title(lgd, 'BSSID');
