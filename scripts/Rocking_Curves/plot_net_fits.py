# -*- coding: utf-8 -*-
"""
Net Intensity Fit Plotter
=========================
Author: Tomasz Stawski, tomasz.stawski@bam.de
Date: 2026-06-19
Version: 0.1

Description:
Reconstructs and plots baseline-subtracted rocking curves alongside individual Gaussian components.
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
        
        # Calculate net raw intensity
        net_raw = raw_int - baseline
        
        # Filter metrics for this scan
        df_scan = df_metrics[(df_metrics['Sample'] == sample) & (df_metrics['Phi (degrees)'] == phi)]
        
        # Reconstruct net fit (sum of Gaussians only, no baseline)
        net_fit = np.zeros_like(theta)
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
                net_fit += peak_y
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
            
        # Create figure for net intensity zoomed plot
        plt.figure(figsize=(10, 6))
        
        # Plot raw net data points
        plt.plot(theta, net_raw, 'o', color='#7f7f7f', markersize=5, alpha=0.7, label='Experimental Net Data (Raw - Baseline)')
        
        # Plot flat zero baseline
        plt.axhline(0, color='red', linestyle='-', linewidth=1.5, label='Subtracted Background (y=0)')
        
        # Plot cumulative net fit envelope
        plt.plot(theta, net_fit, 'k-', linewidth=2.5, label='Total Net Fit Envelope')
        
        # Plot individual peaks in the zoom region starting from y=0
        peak_idx = 0
        for p in peaks_to_plot:
            if zoom_min - 1.0 <= p["center"] <= zoom_max + 1.0:
                color = colors[peak_idx % len(colors)]
                plt.plot(theta, p["curve"], '--', color=color, linewidth=1.8, label=f"{p['name']} ($\\chi \\approx {p['tilt']:.2f}^\\circ$)")
                plt.fill_between(theta, 0, p["curve"], color=color, alpha=0.15)
                plt.axvline(p["center"], color=color, linestyle=':', alpha=0.6)
                peak_idx += 1
                
        plt.xlim(zoom_min, zoom_max)
        
        # Adjust y-limit based on net intensity in the range
        y_data = net_raw[mask]
        y_fit = net_fit[mask]
        y_max = max(np.max(y_data), np.max(y_fit))
        y_min = min(np.min(y_data), 0)
        plt.ylim(y_min - (y_max - y_min)*0.1, y_max + (y_max - y_min)*0.1)
        
        plt.xlabel("Theta (degrees)", fontsize=13)
        plt.ylabel("Net Intensity (counts)", fontsize=13)
        plt.title(f"{sample} Net Fit Zoom — Phi = {phi}° (Baseline Subtracted)", fontsize=14, fontweight='bold')
        plt.grid(True, which='both', linestyle=':', alpha=0.6)
        plt.legend(loc='upper right', fontsize=10, framealpha=0.9)
        
        plt.tight_layout()
        out_path = os.path.join(sample_out_dir, f"{sample}_fit_phi_{phi}_net_zoom.png")
        plt.savefig(out_path, dpi=150)
        plt.close()
        print(f"Generated net zoomed individual plot: {out_path}")
