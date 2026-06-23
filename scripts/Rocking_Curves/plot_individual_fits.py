# -*- coding: utf-8 -*-
"""
Individual Scan Fit Plotter
===========================
Author: Tomasz Stawski, tomasz.stawski@bam.de
Date: 2026-06-19
Version: 0.1

Description:
Plots experimental data and baseline-subtracted fit envelopes for individual rocking curve runs.
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Directories
base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
processed_dir = os.path.join(base_dir, "data/processed/Rocking_Curves")
metrics_path = os.path.join(base_dir, "data/processed/Rocking_Curves/08062026/all_samples_rocking_peaks_vs_phi.csv")
out_dir = os.path.join(base_dir, "data/processed/Rocking_Curves/08062026/Fits_Individual")

# Load metrics
df_metrics = pd.read_csv(metrics_path)

# Samples and their phi values
samples_phi = {
    "SH-124-B3": [0, 30, 60, 90, 120, 150, 180],
    "SH-125-A": [0, 30, 60, 90, 120, 150]
}

# Color palette for peaks
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2']

def gaussian(x, h, x0, w):
    return h * np.exp(-(x - x0)**2 / (2 * w**2))

for sample, phi_list in samples_phi.items():
    sample_out_dir = os.path.join(out_dir, sample)
    os.makedirs(sample_out_dir, exist_ok=True)
    
    for phi in phi_list:
        csv_path = os.path.join(processed_dir, sample, f"{sample}_corrected_rocking_{phi}.csv")
        if not os.path.exists(csv_path):
            print(f"Skipping {sample} Phi={phi}: CSV not found at {csv_path}")
            continue
            
        # Load raw and baseline data
        df_data = pd.read_csv(csv_path)
        theta = df_data['Theta (degrees)'].values
        raw_int = df_data['Raw Intensity'].values
        baseline = df_data['Model Baseline'].values
        
        # Filter metrics for this scan
        df_scan = df_metrics[(df_metrics['Sample'] == sample) & (df_metrics['Phi (degrees)'] == phi)]
        
        # Reconstruct fit
        full_fit = baseline.copy()
        peaks_to_plot = []
        
        for idx, row in df_scan.iterrows():
            name = row['Peak Name']
            center = row['Peak Center (Theta)']
            fwhm = row['FWHM (degrees)']
            height = row['Net Height']
            tilt = row['Tilt Angle (Chi)']
            
            if height > 0 and fwhm > 0:
                w = fwhm / 2.35482
                peak_y = gaussian(theta, height, center, w)
                full_fit += peak_y
                peaks_to_plot.append({
                    "name": name,
                    "center": center,
                    "height": height,
                    "w": w,
                    "tilt": tilt,
                    "curve": peak_y
                })
                
        # Create figure
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 7))
        fig.suptitle(f"{sample} — Phi = {phi}° Rocking Curve Fit", fontsize=16, fontweight='bold')
        
        for ax, is_log in [(ax1, False), (ax2, True)]:
            # Plot raw data
            ax.plot(theta, raw_int, 'o', color='gray', markersize=3, alpha=0.4, label='Raw data')
            # Plot baseline
            ax.plot(theta, baseline, 'r-', linewidth=1.5, label='Background Baseline')
            # Plot full fit envelope
            ax.plot(theta, full_fit, 'k-', linewidth=2, label='Total Fit Envelope')
            
            # Plot individual peaks
            for i, p in enumerate(peaks_to_plot):
                color = colors[i % len(colors)]
                # Peak profile added to baseline
                peak_profile = baseline + p["curve"]
                ax.plot(theta, peak_profile, '--', color=color, linewidth=1.2, label=f"{p['name']} ($\\chi \\approx {p['tilt']:.1f}^\\circ$)")
                
                # Fill between baseline and peak+baseline
                ax.fill_between(theta, baseline, peak_profile, color=color, alpha=0.15)
                
                # Vertical line at peak center
                ax.axvline(p["center"], color=color, linestyle=':', alpha=0.5)
                
            ax.set_xlabel("Theta (degrees)", fontsize=12)
            ax.set_ylabel("Intensity (counts)", fontsize=12)
            if is_log:
                ax.set_yscale('log')
                ax.set_title("Logarithmic Scale", fontsize=13)
                # Set reasonable limits for log scale to focus on the data
                ax.set_ylim(max(10.0, np.min(raw_int)*0.5), np.max(raw_int)*2.0)
            else:
                ax.set_title("Linear Scale", fontsize=13)
                ax.set_ylim(0, np.max(raw_int)*1.1)
                
            ax.grid(True, which='both', linestyle=':', alpha=0.5)
            ax.legend(loc='upper right', fontsize=9, framealpha=0.9)
            
        plt.tight_layout(rect=[0, 0, 1, 0.95])
        out_path = os.path.join(sample_out_dir, f"{sample}_fit_phi_{phi}.png")
        plt.savefig(out_path, dpi=150)
        plt.close()
        print(f"Generated individual plot: {out_path}")
