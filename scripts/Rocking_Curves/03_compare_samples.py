#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Sample Comparison Utility
=========================
Author: Tomasz Stawski, tomasz.stawski@bam.de
Date: 2026-06-23
Version: 1.0

Description:
Compares baseline-corrected rocking curves, peak positions, and intensities
across SH-124-B3 and SH-125-A at representative azimuthal orientations.
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Paths
base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
rc_b3_path = os.path.join(base_dir, "data/processed/Rocking_Curves/SH-124-B3/SH-124-B3_corrected_rocking_60.csv")
rc_a_path = os.path.join(base_dir, "data/processed/Rocking_Curves/SH-125-A/SH-125-A_corrected_rocking_60.csv")
output_dir = os.path.join(base_dir, "data/processed/Rocking_Curves/Comparison")
os.makedirs(output_dir, exist_ok=True)

# Load data helper
def load_csv(path):
    """Loads baseline-corrected rocking curve profiles from a CSV file, returning scattering angle and intensity arrays."""
    if not os.path.exists(path):
        raise FileNotFoundError(f"Processed file not found: {path}")
    df = pd.read_csv(path)
    return df["Theta (degrees)"].values, df["Raw Intensity"].values, df["Model Baseline"].values, df["Corrected Net Intensity"].values

try:
    theta_b3, int_b3, baseline_b3, net_b3 = load_csv(rc_b3_path)
    theta_a, int_a, baseline_a, net_a = load_csv(rc_a_path)
except Exception as e:
    print(f"Error loading processed datasets: {e}")
    sys.exit(1)

# Fitting models
def gaussian(t, h, t0, w):
    """Evaluates a Gaussian peak profile."""
    return h * np.exp(-(t-t0)**2 / (2*w**2))

# Nominal 2Theta values for alignment and Chi calculation
two_theta_b3 = 29.44
two_theta_a = 29.3425

# 1. Fit SH-124-B3 peaks
peaks_b3_init = [9.36, 10.60, 12.10, 12.88, 15.72, 16.30, 17.10, 17.50, 22.34]
ranges_b3 = [(9.0, 9.7), (10.2, 11.0), (11.7, 12.3), (12.4, 13.4), (15.2, 16.1), (16.1, 16.7), (16.8, 17.3), (17.3, 17.9), (21.8, 22.8)]
bounds_b3 = [
    ([0, 9.1, 0.05], [10000, 9.5, 0.5]),
    ([0, 10.3, 0.05], [10000, 10.9, 0.5]),
    ([0, 11.8, 0.05], [10000, 12.2, 0.5]),
    ([0, 12.6, 0.05], [10000, 13.1, 0.4]),
    ([0, 15.5, 0.05], [10000, 16.0, 0.6]),
    ([0, 16.0, 0.05], [10000, 16.5, 0.5]),
    ([0, 16.9, 0.05], [10000, 17.2, 0.5]),
    ([0, 17.4, 0.05], [10000, 17.8, 0.5]),
    ([0, 22.1, 0.05], [10000, 22.6, 0.3])
]
names_b3 = ["B3_Peak_9.22", "B3_Peak_10.60", "B3_Peak_12.10", "B3_Peak_12.92", "B3_Peak_15.70", "B3_Peak_16.19", "B3_Peak_17.10", "B3_Peak_17.50", "B3_Peak_22.34"]
results_b3 = []

for name, center, (t_min, t_max), bnd in zip(names_b3, peaks_b3_init, ranges_b3, bounds_b3):
    mask = (theta_b3 >= t_min) & (theta_b3 <= t_max)
    try:
        popt_g, _ = curve_fit(gaussian, theta_b3[mask], net_b3[mask], p0=[100, center, 0.15], bounds=bnd)
        h, t0, w = popt_g
        fwhm = 2.355 * w
        area = h * w * np.sqrt(2 * np.pi)
        results_b3.append({
            'Sample': 'SH-124-B3',
            'Peak ID': name,
            'Peak Center (Theta)': t0,
            'Tilt Angle (Chi)': t0 - two_theta_b3/2,
            'FWHM (deg)': fwhm,
            'Net Height': h,
            'Net Area (cts deg)': area
        })
    except Exception as e:
        print(f"Error fitting SH-124-B3 peak {name}: {e}")

# 2. Fit SH-125-A peaks
peaks_a_init = [9.97, 12.77, 14.35, 22.52]
ranges_a = [(9.5, 10.5), (12.3, 13.3), (13.7, 15.0), (22.0, 23.0)]
bounds_a = [
    ([0, 9.7, 0.05], [200000, 10.2, 0.4]),
    ([0, 12.5, 0.05], [200000, 13.1, 0.4]),
    ([0, 14.0, 0.05], [200000, 14.6, 0.4]),
    ([0, 22.2, 0.05], [200000, 22.8, 0.3])
]
names_a = ["A_Peak_9.99", "A_Peak_12.81", "A_Peak_14.35", "A_Peak_22.53"]
results_a = []

for name, center, (t_min, t_max), bnd in zip(names_a, peaks_a_init, ranges_a, bounds_a):
    mask = (theta_a >= t_min) & (theta_a <= t_max)
    try:
        popt_g, _ = curve_fit(gaussian, theta_a[mask], net_a[mask], p0=[1000, center, 0.15], bounds=bnd)
        h, t0, w = popt_g
        fwhm = 2.355 * w
        area = h * w * np.sqrt(2 * np.pi)
        results_a.append({
            'Sample': 'SH-125-A',
            'Peak ID': name,
            'Peak Center (Theta)': t0,
            'Tilt Angle (Chi)': t0 - two_theta_a/2,
            'FWHM (deg)': fwhm,
            'Net Height': h,
            'Net Area (cts deg)': area
        })
    except Exception as e:
        print(f"Error fitting SH-125-A peak {name}: {e}")

# Combine results and save comparison table
df_compare = pd.DataFrame(results_b3 + results_a)
df_compare.to_csv(os.path.join(output_dir, "samples_peaks_comparison.csv"), index=False)
print("Saved peaks comparison table.")

# Plot comparison
fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

# Plot SH-124-B3 Residuals vs Chi
chi_b3 = theta_b3 - two_theta_b3/2
axes[0].plot(chi_b3, net_b3, 'g-', label='SH-124-B3 Residuals (Net Intensity)', alpha=0.7)
for r in results_b3:
    chi0 = r['Tilt Angle (Chi)']
    h = r['Net Height']
    w = r['FWHM (deg)'] / 2.355
    axes[0].plot(chi_b3, gaussian(theta_b3, h, r['Peak Center (Theta)'], w), 'r--', alpha=0.8)
    axes[0].text(chi0 + 0.1, h * 0.9, f"$\\chi$={chi0:.2f}°", fontsize=9, color='red')
axes[0].axhline(0, color='gray', linestyle='--')
axes[0].set_ylabel('Net Residual Intensity (counts)')
axes[0].legend(loc='upper right')
axes[0].grid(True)
axes[0].set_title('Residual Orientation Peaks Comparison: SH-124-B3 vs. SH-125-A')

# Plot SH-125-A Residuals vs Chi
chi_a = theta_a - two_theta_a/2
axes[1].plot(chi_a, net_a, 'b-', label='SH-125-A Residuals (Net Intensity)', alpha=0.7)
for r in results_a:
    chi0 = r['Tilt Angle (Chi)']
    h = r['Net Height']
    w = r['FWHM (deg)'] / 2.355
    axes[1].plot(chi_a, gaussian(theta_a, h, r['Peak Center (Theta)'], w), 'r--', alpha=0.8)
    axes[1].text(chi0 + 0.1, h * 0.9, f"$\\chi$={chi0:.2f}°", fontsize=9, color='red')
axes[1].axhline(0, color='gray', linestyle='--')
axes[1].set_xlabel('Tilt Angle Chi (degrees)')
axes[1].set_ylabel('Net Residual Intensity (counts)')
axes[1].legend(loc='upper right')
axes[1].grid(True)

plt.tight_layout()
plot_path = os.path.join(output_dir, "samples_comparison.png")
plt.savefig(plot_path, dpi=150)
plt.close()
print(f"Saved comparison plot to {plot_path}")
