# -*- coding: utf-8 -*-
"""
Zoomed Peak Fit Plotter
=======================
Author: Tomasz Stawski, tomasz.stawski@bam.de
Date: 2026-06-19
Version: 0.1

Description:
Generates high-resolution plots focused on the peak deconvolution region of rocking curve sweeps.
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
                
        # Define zoom window for Peak 2a & 2b (10.0° to 14.5°)
        zoom_min, zoom_max = 10.0, 14.5
        mask = (theta >= zoom_min) & (theta <= zoom_max)
        
        if not mask.any():
            print(f"Skipping {sample} Phi={phi}: No data in zoom range {zoom_min}° - {zoom_max}°")
            continue
            
        # Create figure for zoomed plot
        plt.figure(figsize=(10, 6))
        
        # Plot raw data (show individual points clearly)
        plt.plot(theta, raw_int, 'o', color='#7f7f7f', markersize=4, alpha=0.6, label='Raw Data Points')
        # Plot baseline
        plt.plot(theta, baseline, 'r-', linewidth=2, label='Background Baseline')
        # Plot full fit envelope
        plt.plot(theta, full_fit, 'k-', linewidth=2.5, label='Total Fit Envelope')
        
        # Plot individual peaks in the zoom region
        peak_idx = 0
        for p in peaks_to_plot:
            # Check if peak is within or close to the zoom window
            if zoom_min - 1.0 <= p["center"] <= zoom_max + 1.0:
                color = colors[peak_idx % len(colors)]
                peak_profile = baseline + p["curve"]
                plt.plot(theta, peak_profile, '--', color=color, linewidth=1.5, label=f"{p['name']} ($\\chi \\approx {p['tilt']:.2f}^\\circ$)")
                plt.fill_between(theta, baseline, peak_profile, color=color, alpha=0.15)
                plt.axvline(p["center"], color=color, linestyle=':', alpha=0.6)
                peak_idx += 1
                
        plt.xlim(zoom_min, zoom_max)
        
        # Adjust y-limit to fit the data in the zoom range
        y_data = raw_int[mask]
        y_fit = full_fit[mask]
        y_max = max(np.max(y_data), np.max(y_fit))
        y_min = min(np.min(y_data), np.min(baseline[mask]))
        plt.ylim(max(0, y_min - (y_max - y_min)*0.1), y_max + (y_max - y_min)*0.1)
        
        plt.xlabel("Theta (degrees)", fontsize=13)
        plt.ylabel("Intensity (counts)", fontsize=13)
        plt.title(f"{sample} Fit Zoom — Phi = {phi}° (Peak 2a/2b Region)", fontsize=14, fontweight='bold')
        plt.grid(True, which='both', linestyle=':', alpha=0.6)
        plt.legend(loc='upper right', fontsize=10, framealpha=0.9)
        
        plt.tight_layout()
        out_path = os.path.join(sample_out_dir, f"{sample}_fit_phi_{phi}_zoom.png")
        plt.savefig(out_path, dpi=150)
        plt.close()
        print(f"Generated zoomed individual plot: {out_path}")
