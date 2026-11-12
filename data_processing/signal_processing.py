"""
Receiver — Offline Signal Processing

Outputs: processed_<tag>.mat  in final_data struct layout.

Usage:
    python signal_processing.py --input raw_capture.mat --input noise_capture.mat --output processed_capture.mat
"""

import argparse
import sys
import numpy as np
import scipy.io as sio
import scipy.signal as sps_sig
import matplotlib.pyplot as plt
from scipy.cluster.vq import kmeans2

parser = argparse.ArgumentParser(description='Offline receiver signal processing')
parser.add_argument('--input',  default='../beamforming/data_f2/beam_tx1_0_0_0.mat', help='Input .mat file')
parser.add_argument('--noise',  default='../beamforming/data_f2/beam_tx2_0_0_0.mat', help='Input .mat file')
parser.add_argument('--output', default='processed.mat',        help='Output .mat file')
args = parser.parse_args()

# ==============================================================================
# 1.  CONFIGURATION
# ==============================================================================
Fs           = 30e6 / 30
CaptureTime  = 2.0
TotalSamples = round(Fs * CaptureTime)

SPS     = 5
ROLLOFF = 0.25
SPAN    = 6

NUM_PACKETS    = None
FRAME_SAMPLES  = 15120  # 5000 silence + 2560 waveform + 7
MIN_PEAK_DIST  = int(0.8 * FRAME_SAMPLES)
MIN_PROM_FRAC  = 0.3
pn_bits_S0 = sio.loadmat('waveform_STTD.mat')['bits_S0'].ravel().astype(np.uint8)[6:]
pn_bits_S90 = sio.loadmat('waveform_STTD.mat')['bits_S90'].ravel().astype(np.uint8)[6:]


def qpsk_modulate(bits, phase_offset=np.pi / 4):
    assert len(bits) % 2 == 0
    b = bits.reshape(-1, 2)
    GRAY_MAP = {(0,0): 0, (0,1): 1, (1,1): 2, (1,0): 3}
    gray_idx = np.array([GRAY_MAP[tuple(bb)] for bb in b], dtype=int)
    return np.exp(1j * (phase_offset + np.pi / 2 * gray_idx))


def qpsk_demodulate(syms, phase_offset=np.pi / 4):
    candidates = np.exp(1j * (phase_offset + np.pi / 2 * np.arange(4)))
    dists      = np.abs(syms[:, None] - candidates[None, :])
    idx        = np.argmin(dists, axis=1)
    gray_to_bits = np.array([[0, 0], [0, 1], [1, 1], [1, 0]], dtype=np.uint8)
    return gray_to_bits[idx].ravel()


def rrc_rx_filter(y, h, sps, ref_syms=None):
    """
    Matched filtering with optional fractional timing synchronization.
    
    y: Received signal
    h: RRC filter coefficients
    sps: Samples per symbol
    ref_syms: (Optional) Preamble symbols used to find the best sampling phase.
    """
    # 1. Perform the matched filtering at the full sampling rate
    yf = np.convolve(y, h, mode='same')
    
    if ref_syms is None:
        # Fallback to default behavior if no reference is provided
        return yf[::sps]

    best_phase = 0
    max_corr_sq = -1
    
    # 2. Search for the optimal fractional phase (k)
    # We iterate through all possible starting samples within one symbol period
    for k in range(sps):
        candidate = yf[k::sps]
        
        # Cross-correlate candidate symbols with the known preamble
        # We use 'valid' to find the peak within the sequence
        correlation = np.convolve(candidate, np.conj(ref_syms[::-1]), mode='valid')
        peak_val = np.max(np.abs(correlation)**2)
        
        if peak_val > max_corr_sq:
            max_corr_sq = peak_val
            best_phase = k
            
    # 3. Return the downsampled signal using the best phase found
    return yf[best_phase::sps]


def coarse_cfo_estimate(x, fs):
    x4   = x ** 4
    n    = len(x4)
    nfft = 1 << (n - 1).bit_length()
    X4   = np.fft.fft(x4, n=nfft)
    freqs = np.fft.fftfreq(nfft, d=1/fs)
    mag  = np.abs(X4)
    mag[0] = 0.0
    return float(freqs[int(np.argmax(mag))]) / 4.0


def xcorr(x, y):
    x = np.asarray(x, dtype=complex)
    y = np.asarray(y, dtype=complex)
    corr = sps_sig.correlate(x, y, mode='full', method='auto')
    lags = sps_sig.correlation_lags(len(x), len(y), mode='full')
    return corr, lags


def fit_to(arr, n):
    if len(arr) >= n:
        return arr[:n]
    return np.pad(arr, (0, n - len(arr)))


def kmeans_cluster_means(syms):
    pts = np.column_stack([syms.real, syms.imag])
    centroids, _ = kmeans2(pts, 4, iter=20, minit='points', missing='warn')
    return centroids[:, 0] + 1j * centroids[:, 1]


def compute_evm(syms):
    try:
        pts = np.column_stack([syms.real, syms.imag])
        centroids, labels = kmeans2(pts, 4, iter=20, minit='points', missing='warn')
        errors = pts - centroids[labels]
        return float(np.mean(np.sum(errors**2, axis=1)))
    except Exception:
        return float('NaN')


def collect_valid_packets(lags, pks_locs, pks_vals, n_pkt, buf_len):
    valid = []
    for k in range(len(pks_locs)):
        lag  = int(lags[pks_locs[k]])
        s_py = lag
        e_py = s_py + n_pkt
        if s_py >= 0 and e_py <= buf_len and pks_vals[k] > 0:
            valid.append({'startIdx': s_py, 'peakVal': float(pks_vals[k]), 'peakTime': s_py / Fs})
        if len(valid) >= NUM_PACKETS:
            break
    return valid


def detect_packets(rx, ref_sig):
    rx = np.asarray(rx)
    if rx.ndim == 1:
        corr_mag = np.abs(sps_sig.correlate(
            rx, ref_sig, mode='full', method='auto'))
        lags = sps_sig.correlation_lags(len(rx), len(ref_sig), mode='full')
    else:
        corr1 = sps_sig.correlate(rx[:, 0], ref_sig, mode='full', method='auto')
        corr2 = sps_sig.correlate(rx[:, 1], ref_sig, mode='full', method='auto')
        corr_mag = np.sqrt(np.abs(corr1)**2 + np.abs(corr2)**2)
        lags = sps_sig.correlation_lags(rx.shape[0], len(ref_sig), mode='full')

    valid_lags = (lags >= 0) & (lags + len(ref_sig) <= rx.shape[0])
    corr_mag = corr_mag * valid_lags
    max_corr = float(np.max(corr_mag)) if corr_mag.size else 0.0
    min_prom = MIN_PROM_FRAC * max_corr
    pks_locs, _ = sps_sig.find_peaks(
        corr_mag, distance=MIN_PEAK_DIST, prominence=min_prom)
    pks_vals = corr_mag[pks_locs]
    sort_idx = np.argsort(pks_vals)[::-1]
    return corr_mag, corr_mag, lags, pks_locs[sort_idx], pks_vals[sort_idx]


def collect_valid_packets(lags, pks_locs, pks_vals, n_pkt, buf_len):
    valid = []
    for k in range(len(pks_locs)):
        lag  = int(lags[pks_locs[k]])
        s_py = lag
        e_py = s_py + n_pkt
        if s_py >= 0 and e_py <= buf_len and pks_vals[k] > 0:
            valid.append({'startIdx': s_py, 'peakVal': float(pks_vals[k]), 'peakTime': s_py / Fs})

        if NUM_PACKETS is not None and len(valid) >= NUM_PACKETS:
            break

    # Detection peaks are initially ranked by strength. Restore capture-time
    # order so the angle sequence can be tracked through the recording.
    valid.sort(key=lambda pkt: pkt["startIdx"])
    return valid


def assign_template_detections(lags, corr_s0, corr_s90, packet_len, buffer_len):
    # Combine both template correlation magnitudes.
    combined_corr = np.maximum(corr_s0, corr_s90)

    # Only consider locations where the complete packet fits.
    valid = (lags >= 0) & (lags + packet_len <= buffer_len)

    combined_corr = combined_corr.copy()
    combined_corr[~valid] = 0

    # Detect packet locations using the combined correlation.
    max_corr = np.max(combined_corr)

    peaks, properties = sps_sig.find_peaks(
        combined_corr,
        distance=MIN_PEAK_DIST,
        prominence=MIN_PROM_FRAC * max_corr
    )

    # Classify each detected packet.
    packets_s0 = []
    packets_s90 = []

    SEARCH_RADIUS = 50

    for peak in peaks:
        lag = int(lags[peak])

        # Compare both templates around the same packet location.
        lo = max(0, peak - SEARCH_RADIUS)
        hi = min(len(lags), peak + SEARCH_RADIUS + 1)

        score_s0 = np.max(corr_s0[lo:hi])
        score_s90 = np.max(corr_s90[lo:hi])

        print(
            f"Packet at {lag/Fs:.5f}s | "
            f"S0 score: {score_s0:.2f} | "
            f"S90 score: {score_s90:.2f}"
        )

        if score_s0 >= score_s90:

            packets_s0.append({
                'startIdx': lag,
                'peakVal': float(score_s0),
                'peakTime': lag / Fs,
                'template': 'S0'
            })

        else:

            packets_s90.append({
                'startIdx': lag,
                'peakVal': float(score_s90),
                'peakTime': lag / Fs,
                'template': 'S90'
            })

    packets_s0.sort(key=lambda p: p['startIdx'])
    packets_s90.sort(key=lambda p: p['startIdx'])

    return packets_s0, packets_s90


def process_packets(valid_packets, rxBuffer, ref_syms, pn_bits, noise_power):
    NUM_REF_SYMS = len(ref_syms)
    num_p       = len(valid_packets)
    all_syms_ant1 = np.zeros((NUM_REF_SYMS, num_p), dtype=complex)
    all_syms_ant2 = np.zeros((NUM_REF_SYMS, num_p), dtype=complex)
    results = []

    for p, pkt in enumerate(valid_packets):
        s = max(0, pkt['startIdx'] - 1000)  # extraction begins before detected packet
        extraction_len = N_PKT - (SPAN // 2 * SPS) + 2000  # include trailing margin


        pkt_ant1 = rxBuffer[s:s+extraction_len, 0]
        pkt_ant2 = rxBuffer[s:s+extraction_len, 1]

        est_cfo = coarse_cfo_estimate(pkt_ant1, Fs)
        t_vec   = np.arange(len(pkt_ant2)) / Fs
        cfo_vec = np.exp(-1j * 2 * np.pi * est_cfo * t_vec)

        # Print the normalized frequency offset
        normalized_cfo = est_cfo / Fs
        print(f"CFO: {est_cfo:.2f} Hz")
        print(f"Phase rotation per symbol: {normalized_cfo * 360:.2f} degrees")

        clean_ant1 = pkt_ant1 * cfo_vec
        clean_ant2 = pkt_ant2 * cfo_vec

        if (len(pkt_ant1) == extraction_len and len(pkt_ant2) == extraction_len):
            
            mf_ant1 = rrc_rx_filter(clean_ant1, H_rrc, SPS, ref_syms=ref_syms)
            mf_ant2 = rrc_rx_filter(clean_ant2, H_rrc, SPS, ref_syms=ref_syms)

            fine_corr1, fine_lags1 = xcorr(mf_ant1, ref_syms)
            fine_idx1              = int(np.argmax(np.abs(fine_corr1)))
            los_phase_ref         = np.angle(fine_corr1[fine_idx1])
            phase_anchor          = np.exp(-1j * los_phase_ref)

            # In your packet detection block:
            peak_val = np.max(np.abs(fine_corr1))
            mean_val = np.mean(np.abs(fine_corr1))
            print(f"Correlation Peak-to-Average Ratio: {peak_val / mean_val:.2f}")

            fine_corr2, fine_lags2 = xcorr(mf_ant2, ref_syms)
            fine_idx2              = int(np.argmax(np.abs(fine_corr2)))

            mf_ant1_anchored = mf_ant1 * phase_anchor
            mf_ant2_anchored = mf_ant2 * phase_anchor

            fine_timing_offset1 = int(fine_lags1[fine_idx1])
            start_idx1 = max(0, min(fine_timing_offset1, len(mf_ant1) - NUM_REF_SYMS))
            end_idx1   = start_idx1 + NUM_REF_SYMS
            fine_timing_offset2 = int(fine_lags2[fine_idx2])
            start_idx2 = max(0, min(fine_timing_offset2, len(mf_ant2) - NUM_REF_SYMS))
            end_idx2   = start_idx2 + NUM_REF_SYMS

            syms_ant1 = mf_ant1_anchored[start_idx1:end_idx1]
            syms_ant2 = mf_ant2_anchored[start_idx2:end_idx2]

            symbol_indices_ant1 = range(len(syms_ant1))#np.where(np.abs(syms_ant1) > 0.5 * np.mean(np.abs(syms_ant1)))[0]
            symbol_indices_ant2 = range(len(syms_ant2))#np.where(np.abs(syms_ant2) > 0.5 * np.mean(np.abs(syms_ant2)))[0]

            syms_ant1 = syms_ant1[symbol_indices_ant1]
            syms_ant2 = syms_ant2[symbol_indices_ant2]

            # if p < 5:  # First 5 packets
            #     # Plot symbol magnitudes vs symbol index
            #     plt.figure(figsize=(10, 4))
            #     plt.plot(np.arange(len(syms_ant1)), np.abs(syms_ant1), marker='o', linestyle='-', label='Antenna 1')
            #     plt.plot(np.arange(len(syms_ant2)), np.abs(syms_ant2), marker='x', linestyle='-', label='Antenna 2')
            #     plt.xlabel('Symbol Index within Packet')
            #     plt.ylabel('Magnitude')
            #     plt.title(f'Packet {p} Symbol Magnitudes')
            #     plt.grid(True)
            #     plt.legend()
            #     plt.tight_layout()
            #     plt.show()

            all_syms_ant1[:, p] = fit_to(syms_ant1, NUM_REF_SYMS)
            all_syms_ant2[:, p] = fit_to(syms_ant2, NUM_REF_SYMS)
        

            # 1. Extract the symbols
            syms_to_demod1 = mf_ant1[start_idx1:end_idx1]
            syms_to_demod2 = mf_ant2[start_idx2:end_idx2]

            # Check for phase drift across the packet
            first_sym_phase = np.angle(syms_to_demod1[0] * np.conj(ref_syms[0]))
            last_sym_phase = np.angle(syms_to_demod1[-1] * np.conj(ref_syms[len(syms_to_demod1)-1]))
            print(f"Phase Drift across preamble: {np.degrees(last_sym_phase - first_sym_phase):.2f} degrees")

            # Estimate packet SNR from time-domain packet power and the
            # independently captured noise-only power for each receive channel.
            # packet_power1 = float(np.mean(np.abs(pkt_ant1)**2))
            # packet_power2 = float(np.mean(np.abs(pkt_ant2)**2))
            # signal_power1 = packet_power1 - noise_power[0]
            # signal_power2 = packet_power2 - noise_power[1]

            packet_start = pkt['startIdx']
            packet_end = packet_start + N_PKT

            actual_packet1 = rxBuffer[packet_start:packet_end, 0]
            actual_packet2 = rxBuffer[packet_start:packet_end, 1]

            packet_power1 = np.mean(np.abs(actual_packet1)**2)
            packet_power2 = np.mean(np.abs(actual_packet2)**2)
            signal_power1 = packet_power1 - noise_power[0]
            signal_power2 = packet_power2 - noise_power[1]
            snr1_db = (10.0 * np.log10(signal_power1 / noise_power[0])
                       if signal_power1 > 0 else np.nan)
            snr2_db = (10.0 * np.log10(signal_power2 / noise_power[1])
                       if signal_power2 > 0 else np.nan)
            print(f"Packet SNR: Ant1={snr1_db:.2f} dB | Ant2={snr2_db:.2f} dB")

            # We compare the first N symbols of our extracted window to the known reference symbols
            num_ref = min(len(syms_to_demod1), len(ref_syms))
            phase_corr1 = np.angle(np.sum(syms_to_demod1[:num_ref] * np.conj(ref_syms[:num_ref])))
            phase_corr2 = np.angle(np.sum(syms_to_demod2[:num_ref] * np.conj(ref_syms[:num_ref])))

            # 3. Demodulate
            bits_out1 = qpsk_demodulate(syms_to_demod1 * np.exp(-1j * phase_corr1))
            bits_out2 = qpsk_demodulate(syms_to_demod2 * np.exp(-1j * phase_corr2))

            print(f"Truth Bits: {pn_bits[:20]}")
            print(f"RX Bits:    {bits_out1[:20]}")

            # 4. The BER calculation
            L = min(len(pn_bits), len(bits_out1))
        
            n_err1 = int(np.sum(pn_bits[:L] != bits_out1[:L]))
            n_err2 = int(np.sum(pn_bits[:L] != bits_out2[:L]))

            ber1 = n_err1 / L
            ber2 = n_err2 / L

            csi1 = np.mean(syms_ant1 / ref_syms[symbol_indices_ant1])
            csi2 = np.mean(syms_ant2 / ref_syms[symbol_indices_ant2])
            if p == 0:
                print(f"CSI1 {csi1}  CSI2 {csi2}")
            results.append({
                'ber1': ber1,
                'ber2': ber2,
                'snr1_db': snr1_db,
                'snr2_db': snr2_db,
                'packet_power1': packet_power1,
                'packet_power2': packet_power2,
                'cfo': est_cfo,
                'csi1': csi1,
                'csi2': csi2,
                'rel_phase': np.angle(csi2),
                'peakVal': pkt['peakVal'],
                'peakTime': pkt['peakTime']
            })

    return all_syms_ant1, all_syms_ant2, results


def cluster_packets(all_syms_ant1, all_syms_ant2):
    num_p = all_syms_ant1.shape[1]
    cluster_means_ant1 = np.full((4, num_p), np.nan, dtype=complex)
    cluster_means_ant2 = np.full((4, num_p), np.nan, dtype=complex)
    for p in range(num_p):
        try:
            cluster_means_ant1[:, p] = kmeans_cluster_means(all_syms_ant1[:, p])
        except Exception:
            pass
        try:
            cluster_means_ant2[:, p] = kmeans_cluster_means(all_syms_ant2[:, p])
        except Exception:
            pass
    return cluster_means_ant1, cluster_means_ant2


def plot_constellations(all_syms_ant1, all_syms_ant2, valid_packets, label):
    sample_idxs = range(min(len(valid_packets), 10))
    n_cols = max(len(sample_idxs), 2)
    fig, axes = plt.subplots(2, n_cols, figsize=(4 * n_cols, 8))
    if n_cols == 1:
        axes = axes.reshape(2, 1)
    for si, p in enumerate(sample_idxs):
        syms0 = all_syms_ant1[:, p]
        syms90 = all_syms_ant2[:, p]
        ax = axes[0, si]
        ax.scatter(syms0.real, syms0.imag, s=8, c='b', alpha=0.3)
        ax.grid(True); ax.set_aspect('equal')
        ax.set_title(f'{label} Pkt {p+1} | Ant1 | EVM={compute_evm(syms0):.3f}', fontsize=8)
        ax.set_xlabel('I'); ax.set_ylabel('Q')
        ax = axes[1, si]
        ax.scatter(syms90.real, syms90.imag, s=8, c='r', alpha=0.3)
        ax.grid(True); ax.set_aspect('equal')
        ax.set_title(f'{label} Pkt {p+1} | Ant2 | EVM={compute_evm(syms90):.3f}', fontsize=8)
        ax.set_xlabel('I'); ax.set_ylabel('Q')
    for si in range(len(sample_idxs), n_cols):
        axes[0, si].set_visible(False)
        axes[1, si].set_visible(False)
    fig.tight_layout()


def plot_packet_snr(results_s1, results_s2, output_path):
    """Plot per-packet SNR (dB) for both receive antennas and packet streams."""
    fig, axes = plt.subplots(2, 1, figsize=(11, 7), sharex=False)
    for ax, results, label in (
        (axes[0], results_s1, "S0"),
        (axes[1], results_s2, "S90"),
    ):
        if not results:
            ax.text(0.5, 0.5, "No packets processed", ha="center", va="center",
                    transform=ax.transAxes)
        else:
            packet_idx = np.arange(1, len(results) + 1)
            snr1 = np.array([r["snr1_db"] for r in results])
            snr2 = np.array([r["snr2_db"] for r in results])
            ax.plot(packet_idx, snr1, "o-", markersize=4, label="Antenna 1")
            ax.plot(packet_idx, snr2, "s-", markersize=4, label="Antenna 2")
            ax.legend()
        ax.set_title(f"{label}: Per-packet SNR")
        ax.set_xlabel("Packet number")
        ax.set_ylabel("SNR (dB)")
        ax.grid(True, alpha=0.3)
    fig.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=160, bbox_inches="tight")
        print(f"SNR plot saved -> {output_path}")


def _smooth_binary_labels(labels, window=3):
    """Remove isolated one-packet label flips while preserving longer blocks."""
    labels = np.asarray(labels, dtype=int)
    if len(labels) < window:
        return labels
    if window % 2 == 0:
        window += 1
    pad = window // 2
    padded = np.pad(labels, (pad, pad), mode="edge")
    return np.array([
        int(np.median(padded[i:i + window]) >= 0.5)
        for i in range(len(labels))
    ])


def classify_packet_angles(results_s1, results_s2, noise_power):
    tagged = [(r, "S0") for r in results_s1] + [(r, "S90") for r in results_s2]
    if len(tagged) < 2:
        raise RuntimeError("Need at least two processed packets to infer two angles.")

    strengths = []
    for result, _stream in tagged:
        # Normalize each channel by its independently measured noise floor,
        # then combine the two receive channels.
        snr1_linear = max(result["packet_power1"] - noise_power[0], 0.0) / noise_power[0]
        snr2_linear = max(result["packet_power2"] - noise_power[1], 0.0) / noise_power[1]
        strength_linear = max(0.5 * (snr1_linear + snr2_linear), 1e-12)
        strengths.append(10.0 * np.log10(strength_linear))

    strengths = np.asarray(strengths)

    # Deterministic 1-D two-means clustering. Initialize at the min and max
    # to avoid random cluster initialization changing the angle assignments.
    centers = np.array([np.min(strengths), np.max(strengths)], dtype=float)
    if np.isclose(centers[0], centers[1]):
        print("WARNING: Packet powers do not show two distinguishable levels; "
              "all packets will be labeled theta=0 for now.")
        labels = np.zeros(len(strengths), dtype=int)
    else:
        for _ in range(100):
            distances = np.abs(strengths[:, None] - centers[None, :])
            new_labels = np.argmin(distances, axis=1)
            new_centers = centers.copy()
            for cluster in (0, 1):
                members = strengths[new_labels == cluster]
                if len(members):
                    new_centers[cluster] = np.mean(members)
            if np.allclose(new_centers, centers):
                break
            centers = new_centers
        labels = np.argmin(np.abs(strengths[:, None] - centers[None, :]), axis=1)

    high_cluster = int(np.argmax(centers))
    # First assign the global two-level result to each stream.
    cursor = 0
    for result, _stream in tagged:
        result["rx_strength_db"] = float(strengths[cursor])
        result["theta_cluster"] = int(labels[cursor])
        cursor += 1

    # Smooth each chronological stream independently. This suppresses a
    # single bad power estimate without assuming a fixed dwell time.
    for results in (results_s1, results_s2):
        results.sort(key=lambda r: r["peakTime"])
        raw_labels = np.array([r["theta_cluster"] for r in results], dtype=int)
        smooth_labels = _smooth_binary_labels(raw_labels, window=3)
        for result, cluster in zip(results, smooth_labels):
            result["theta_deg"] = 0 if int(cluster) == high_cluster else 90

    high_values = strengths[labels == high_cluster]
    low_values = strengths[labels != high_cluster]
    threshold = float(0.5 * (np.mean(high_values) + np.mean(low_values))) \
        if len(high_values) and len(low_values) else float("nan")

    print("\nAutomatic angle classification:")
    print(f"  Higher-power cluster -> theta=0 degrees (mean {np.mean(high_values):.2f} dB)")
    if len(low_values):
        print(f"  Lower-power cluster  -> theta=90 degrees (mean {np.mean(low_values):.2f} dB)")
    print(f"  Approximate power threshold: {threshold:.2f} dB")
    for label, results in (("S0", results_s1), ("S90", results_s2)):
        n0 = sum(r["theta_deg"] == 0 for r in results)
        n90 = sum(r["theta_deg"] == 90 for r in results)
        print(f"  {label}: theta=0: {n0} packets | theta=90: {n90} packets")

    return threshold


def plot_snr_for_angle(results_s1, results_s2, theta_deg, output_path):
    """Create a separate per-packet SNR figure for one inferred angle."""
    fig, axes = plt.subplots(2, 1, figsize=(11, 7), sharex=False)
    for ax, results, stream in (
        (axes[0], results_s1, "S0"),
        (axes[1], results_s2, "S90"),
    ):
        selected = [r for r in results if r.get("theta_deg") == theta_deg]
        if not selected:
            ax.text(0.5, 0.5, f"No theta={theta_deg} packets",
                    ha="center", va="center", transform=ax.transAxes)
        else:
            times = np.array([r["peakTime"] for r in selected])
            ax.plot(times, [r["snr1_db"] for r in selected],
                    "o-", markersize=4, label="Antenna 1")
            ax.plot(times, [r["snr2_db"] for r in selected],
                    "s-", markersize=4, label="Antenna 2")
            ax.legend()
        ax.set_title(f"{stream}: Per-packet SNR at theta={theta_deg} degrees")
        ax.set_xlabel("Capture time (s)")
        ax.set_ylabel("SNR (dB)")
        ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(output_path, dpi=160, bbox_inches="tight")
    print(f"Angle-specific SNR plot saved -> {output_path}")


def plot_angle_strength(results_s1, results_s2, output_path):
    """Show the inferred angle sequence and packet-strength separation."""
    fig, axes = plt.subplots(2, 1, figsize=(12, 7), sharex=True)
    for ax, results, stream in (
        (axes[0], results_s1, "S0"),
        (axes[1], results_s2, "S90"),
    ):
        if results:
            times = np.array([r["peakTime"] for r in results])
            strength = np.array([r["rx_strength_db"] for r in results])
            theta0 = np.array([r["theta_deg"] == 0 for r in results])
            ax.scatter(times[theta0], strength[theta0], label="theta=0 (higher power)")
            ax.scatter(times[~theta0], strength[~theta0], label="theta=90 (lower power)")
            ax.legend()
        ax.set_ylabel("Relative received power (dB)")
        ax.set_title(f"{stream}: Inferred angle by packet")
        ax.grid(True, alpha=0.3)
    axes[-1].set_xlabel("Capture time (s)")
    fig.tight_layout()
    fig.savefig(output_path, dpi=160, bbox_inches="tight")
    print(f"Angle classification plot saved -> {output_path}")


def plot_constellations_by_angle(results, all_syms_ant1, all_syms_ant2,
                                 stream_label, theta_deg, output_path):
    """Plot up to the first 10 constellations belonging to one angle."""
    selected = [r for r in results if r.get("theta_deg") == theta_deg]
    if not selected:
        print(f"No {stream_label} packets for theta={theta_deg}; skipping constellation plot.")
        return

    # packet_index maps each result back to its original symbol-array column.
    indices = [r["packet_index"] for r in selected[:10]]
    syms1 = all_syms_ant1[:, indices]
    syms2 = all_syms_ant2[:, indices]
    pseudo_packets = list(range(len(indices)))
    plot_constellations(syms1, syms2, pseudo_packets,
                        f"{stream_label}, theta={theta_deg}")
    plt.gcf().savefig(output_path, dpi=160, bbox_inches="tight")
    print(f"Angle-specific constellation plot saved -> {output_path}")


def print_summary(label, results):
    ber1_arr = np.array([r['ber1'] for r in results])
    ber2_arr = np.array([r['ber2'] for r in results])
    cfo_arr  = np.array([r['cfo']  for r in results])
    print(f"\n================ {label} SUMMARY ================")
    print(f"Total packets processed: {len(results)}")
    print(f"\nAntenna 1:")
    print(f"  Mean BER:    {np.mean(ber2_arr):.6f}")
    print(f"  Median BER:  {np.median(ber2_arr):.6f}")
    print(f"  Min BER:     {np.min(ber2_arr):.6f}")
    print(f"  Max BER:     {np.max(ber2_arr):.6f}")
    print(f"\nCFO Statistics:")
    print(f"  Mean CFO:    {np.mean(cfo_arr):.1f} Hz")
    print(f"  Std CFO:     {np.std(cfo_arr):.1f} Hz")
    print("=" * (len(label) + 27) + "\n")


# ==============================================================================
# 2.  BUILD REFERENCE SIGNALS
# ==============================================================================
H_rrc = sio.loadmat("rrcTx_coeffs.mat")["coeffs"].ravel()
# H_rrc = sio.loadmat("rrcTx_coeffs.mat")["coeffs"]['Numerator'].ravel()[0].ravel().astype(np.float64)
print(f"\nRRC coeffs: sum={np.sum(H_rrc):.6f} (should be ≈5.0 for RX)")
print(f"RRC energy: {np.sum(np.abs(H_rrc)**2):.6f}")

ref_syms_S0  = qpsk_modulate(pn_bits_S0)
ref_syms_S90  = qpsk_modulate(pn_bits_S90)
NUM_REF_SYMS = len(ref_syms_S0)

ref_sig_S0 = sio.loadmat('waveform_STTD.mat')['waveform_S0'].ravel()
ref_sig_S90 = sio.loadmat('waveform_STTD.mat')['waveform_S90'].ravel()

# Check correlation between the two reference waveforms

corr_00 = np.abs(np.vdot(ref_sig_S0, ref_sig_S0))
corr_90 = np.abs(np.vdot(ref_sig_S90, ref_sig_S90))

corr_0_90 = np.abs(np.vdot(ref_sig_S0, ref_sig_S90))

print("\nReference waveform comparison:")
print(f"S0 autocorrelation:    {corr_00:.2f}")
print(f"S90 autocorrelation:   {corr_90:.2f}")
print(f"S0-S90 cross-corr:     {corr_0_90:.2f}")
print(f"Normalized cross-corr: {corr_0_90 / np.sqrt(corr_00*corr_90):.4f}")

# And verify what qpsk_modulate produces matches ref_syms
ref_syms_check = qpsk_modulate(pn_bits_S0)
print("ref_syms_S0[:5]:", ref_syms_S0[:5])
print("modulated check[:5]:", ref_syms_check[:5])
print("match:", np.allclose(ref_syms_S0, ref_syms_check))

N_PKT = len(ref_sig_S0)
if len(ref_sig_S90) != N_PKT:
    sys.exit(f"S0/S90 waveform lengths differ: {N_PKT} vs {len(ref_sig_S90)}")
if N_PKT != NUM_REF_SYMS * SPS:
    print(f"WARNING: reference waveform has {N_PKT} samples but {NUM_REF_SYMS} symbols at {SPS} SPS")
print(f"Reference signal: {N_PKT} samples | {NUM_REF_SYMS} symbols")

print(len(ref_sig_S90))
# ==============================================================================
# 3.  LOAD IQ DATA
# ==============================================================================
print(f"\\nLoading raw capture {args.input} ...")
mat = sio.loadmat(args.input, squeeze_me=False)
raw_data = mat['raw_data']
iq_data = raw_data['iq'][0, 0]
rx_pos = raw_data['rx_pos'][0, 0].flatten().astype(float)
tx_pos = raw_data['tx_pos'][0, 0].flatten().astype(float)

rxBuffer = np.asarray(iq_data, dtype=np.complex128)
if rxBuffer.ndim != 2 or rxBuffer.shape[1] < 2:
    sys.exit(f"Expected raw IQ data with shape (samples, >=2 channels); got {rxBuffer.shape}")
rxBuffer = rxBuffer[:, :2] - np.mean(rxBuffer[:, :2], axis=0, keepdims=True)
CaptureTime = rxBuffer.shape[0] / Fs
TotalSamples = rxBuffer.shape[0]
print(f"RX sample rate: {Fs/1e6:.3f} MSps | capture: {CaptureTime:.3f} s | samples: {len(rxBuffer)}")

print(f"Loading noise-only capture {args.noise} ...")
noise_mat = sio.loadmat(args.noise, squeeze_me=False)
noise_raw = noise_mat['raw_data']
noise_iq = np.asarray(noise_raw['iq'][0, 0], dtype=np.complex128)
if noise_iq.ndim != 2 or noise_iq.shape[1] < 2:
    sys.exit(f"Expected noise IQ data with shape (samples, >=2 channels); got {noise_iq.shape}")
noise_iq = noise_iq[:, :2] - np.mean(noise_iq[:, :2], axis=0, keepdims=True)
noise_power = np.mean(np.abs(noise_iq)**2, axis=0)
if np.any(noise_power <= 0):
    sys.exit(f"Noise power must be positive; measured {noise_power}")

print(f"Noise power: Ant1={noise_power[0]:.6g}, Ant2={noise_power[1]:.6g}")

# ==============================================================================
# 4.  PACKET DETECTION
# ==============================================================================
print("Cross-correlation for packet detection...")

rx1 = rxBuffer[:, 0]
rx2 = rxBuffer[:, 1]

corrVals_S0, corrMag_S0, lags_S0, pks_locs_S0, pks_vals_S0 = detect_packets(rxBuffer, ref_sig_S0)

corrVals_S90, corrMag_S90, lags_S90, pks_locs_S90, pks_vals_S90 = detect_packets(rxBuffer, ref_sig_S90)

valid_packets_S0, valid_packets_S90 = assign_template_detections(
    lags_S0,
    corrMag_S0,
    corrMag_S90,
    N_PKT,
    rxBuffer.shape[0]
)

if not valid_packets_S0:
    sys.exit("No valid S0 packets found above threshold.")
if not valid_packets_S90:
    sys.exit("No valid S90 packets found above threshold.")

print(f"Found {len(valid_packets_S0)} valid S0 packets")
print(f"Found {len(valid_packets_S90)} valid S90 packets")
print(f"S0 start indices: {[p['startIdx'] for p in valid_packets_S0]}")
print(f"S90 start indices: {[p['startIdx'] for p in valid_packets_S90]}\n")

all_packets = (
    [(p['startIdx'], 'S0') for p in valid_packets_S0] +
    [(p['startIdx'], 'S90') for p in valid_packets_S90]
)

all_packets.sort(key=lambda x: x[0])

print("\nChronological packet sequence:")

for start, label in all_packets:
    print(f"Time: {start/Fs:.5f}s | Template: {label}")
# ==============================================================================
# 5.  CORRELATION PLOTS
# ==============================================================================
for label, lags, corrMag, valid_packets in [
    ('S0', lags_S0, corrMag_S0, valid_packets_S0),
    ('S90', lags_S90, corrMag_S90, valid_packets_S90),
]:
    causal_mask = (lags >= 0) & (lags <= (TotalSamples - N_PKT))
    causal_lags = lags[causal_mask]
    causal_corr = corrMag[causal_mask]
    causal_time = causal_lags / Fs
    causal_norm = causal_corr / np.max(causal_corr)

    fig, ax = plt.subplots(figsize=(14, 4))
    ax.plot(causal_time, causal_norm, color=[0.2, 0.45, 0.8], linewidth=0.8, label='|Xcorr|')
    for k, pkt in enumerate(valid_packets):
        t_pk    = pkt['startIdx'] / Fs
        closest = int(np.argmin(np.abs(causal_time - t_pk)))
        ph      = float(causal_norm[closest])
        ax.plot([t_pk, t_pk], [0, ph], 'r-', linewidth=1.4)
        ax.plot(t_pk, ph, 'ro', markersize=6, markerfacecolor='r')
        # ax.text(t_pk, min(ph + 0.05, 1.12), f'P{k+1}', fontsize=7, ha='center', color='r')
    ax.axhline(0.75, color='k', linestyle='--', linewidth=1, label='75% threshold')
    ax.set_xlim([0, CaptureTime]); ax.set_ylim([0, 1.15])
    ax.grid(True); ax.set_xlabel('Time (s)'); ax.set_ylabel('Normalised |Xcorr|')
    ax.set_title(f'{label} Cross-Correlation vs Time — {len(valid_packets)} Packets Detected')
    ax.legend(loc='upper right')
    fig.tight_layout()

# ==============================================================================
# 6.  PER-PACKET PROCESSING
# ==============================================================================
print("Processing S0 packets...")
all_syms_ant1_S0, all_syms_ant2_S0, results_S0 = process_packets(valid_packets_S0, rxBuffer, ref_syms_S0, pn_bits_S0, noise_power)

print("Processing S90 packets...")
all_syms_ant1_S90, all_syms_ant2_S90, results_S90 = process_packets(valid_packets_S90, rxBuffer, ref_syms_S90, pn_bits_S90, noise_power)

# ==============================================================================
# 7.  CLUSTERING
# ==============================================================================
print("\nComputing cluster means...")
cluster_means_ant1_S0, cluster_means_ant2_S0 = cluster_packets(all_syms_ant1_S0, all_syms_ant2_S0)
cluster_means_ant1_S90, cluster_means_ant2_S90 = cluster_packets(all_syms_ant1_S90, all_syms_ant2_S90)

# ==============================================================================
# 8.  CONSTELLATION PLOTS
# ==============================================================================
plot_constellations(all_syms_ant1_S0, all_syms_ant2_S0, valid_packets_S0, 'S0')
plot_constellations(all_syms_ant1_S90, all_syms_ant2_S90, valid_packets_S90, 'S90')

# ==============================================================================
# 9.  SUMMARY STATISTICS
# ==============================================================================
print_summary('S0', results_S0)
print_summary('S90', results_S90)

snr_plot_path = args.output.rsplit('.', 1)[0] + '_snr.png'
plot_packet_snr(results_S0, results_S90, snr_plot_path)

# ==============================================================================
# 11.  SAVE RESULTS
# ==============================================================================
BER_THRESHOLD = 0.001

def filter_by_ber(results, syms_ant1, syms_ant2, threshold=BER_THRESHOLD):
    """Return results and symbol arrays keeping only packets with BER1 <= threshold."""
    keep = [i for i, r in enumerate(results) if r['ber1'] <= threshold]
    n_orig = len(results)
    print(f"  BER filter: keeping {len(keep)}/{n_orig} packets (BER1 ≤ {threshold})")
    filtered_results = [results[i] for i in keep]
    filtered_ant1    = syms_ant1[:, keep] #if keep else np.zeros((syms_ant1.shape[0], 0), dtype=complex)
    filtered_ant2    = syms_ant2[:, keep] #if keep else np.zeros((syms_ant2.shape[0], 0), dtype=complex)
    return filtered_results, filtered_ant1, filtered_ant2, len(keep)

print("\nApplying BER filter before saving...")
results_S0_f, all_syms_ant1_S0_f, all_syms_ant2_S0_f, num_pkts = filter_by_ber(
    results_S0, all_syms_ant1_S0, all_syms_ant2_S0)
results_S90_f, all_syms_ant1_S90_f, all_syms_ant2_S90_f, num_pkts = filter_by_ber(
    results_S90, all_syms_ant1_S90, all_syms_ant2_S90)

# plot_constellations(all_syms_ant1_S0, all_syms_ant2_S0, range(num_pkts), 'S0')
# plot_constellations(all_syms_ant1_S90, all_syms_ant2_S90, range(num_pkts), 'S90')

# Re-cluster on the filtered symbol sets
cluster_means_ant1_S0_f, cluster_means_ant2_S0_f = cluster_packets(all_syms_ant1_S0_f, all_syms_ant2_S0_f)
cluster_means_ant1_S90_f, cluster_means_ant2_S90_f = cluster_packets(all_syms_ant1_S90_f, all_syms_ant2_S90_f)

def make_results_struct(results):
    return {
        'ber1':     np.array([r['ber1']     for r in results]).reshape(1, -1),
        'cfo':      np.array([r['cfo']      for r in results]).reshape(1, -1),
        'peakVal':  np.array([r['peakVal']  for r in results]).reshape(1, -1),
        'peakTime': np.array([r['peakTime'] for r in results]).reshape(1, -1),
    }

final_data = {
    'results_S0':             np.array([make_results_struct(results_S0_f)]),
    'results_S90':             np.array([make_results_struct(results_S90_f)]),
    'num_packets_S0':         np.array([[len(results_S0_f)]], dtype=float),
    'num_packets_S90':         np.array([[len(results_S90_f)]], dtype=float),
    'all_syms_ant1_S0':       all_syms_ant1_S0_f,
    'all_syms_ant2_S0':       all_syms_ant2_S0_f,
    'all_syms_ant1_S90':       all_syms_ant1_S90_f,
    'all_syms_ant2_S90':       all_syms_ant2_S90_f,
    'cluster_means_ant1_S0':  cluster_means_ant1_S0_f,
    'cluster_means_ant2_S0':  cluster_means_ant2_S0_f,
    'cluster_means_ant1_S90':  cluster_means_ant1_S90_f,
    'cluster_means_ant2_S90':  cluster_means_ant2_S90_f,
    'rx_pos':                 rx_pos.reshape(1, -1),
    'tx_pos':                 tx_pos.reshape(1, -1),
    'freq':                   np.array([[915e6]]),
    'tx_bf_angle':            np.array([[np.nan]]),
    'rx_orientation':         np.array([[np.nan]]),
    'tx_mimo':                np.array([[np.nan]]),
}

sio.savemat(args.output, {'final_data': final_data})
print(f"Results saved → {args.output}")

plt.show()