# -*- coding: utf-8 -*-
"""
Sample Comparison Utility
=========================
Author: Tomasz Stawski, tomasz.stawski@bam.de
Date: 2026-06-19
Version: 0.1

Description:
Compares rocking curve shapes, peak positions, and intensities across different thin film growth conditions.
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Paths
base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
rc_20250526_xy = os.path.join(base_dir, "data/processed/Rocking_Curves/SH_124_B3/SH_124_B3_rocking-curve_01_exported.xy")
rc_20260601_xy = os.path.join(base_dir, "data/processed/Rocking_Curves/SH-125-A/SH-125-A_rocking-curve_01_exported.xy")
output_dir = os.path.join(base_dir, "data/processed/Rocking_Curves/Comparison")
os.makedirs(output_dir, exist_ok=True)

# Load data helper
def load_xy(path):
    data = []
    with open(path, 'r') as f:
        f.readline() # skip header
        for line in f:
            parts = line.strip().split()
            if len(parts) == 2:
                data.append([float(parts[0]), float(parts[1])])
    arr = np.array(data)
    return arr[:, 0], arr[:, 1]

theta_b3, int_b3 = load_xy(rc_20250526_xy)
theta_a, int_a = load_xy(rc_20260601_xy)

# Fitting models
def bg_model(t, I0, c0, c1, c2, c3):
    return I0 / np.sin(np.radians(t)) + c0 + c1*t + c2*t**2 + c3*t**3

def gaussian(t, h, t0, w):
    return h * np.exp(-(t-t0)**2 / (2*w**2))

# 1. Fit SH_124_B3 Background
# Exclude peaks: 9.0-9.7, 12.2-13.6, 15.0-16.8, and low angles < 7.0
bg_mask_b3 = (theta_b3 >= 7.0) & \
             ((theta_b3 < 9.0) | (theta_b3 > 9.7)) & \
             ((theta_b3 < 12.2) | (theta_b3 > 13.6)) & \
             ((theta_b3 < 15.0) | (theta_b3 > 16.8))

popt_b3, _ = curve_fit(bg_model, theta_b3[bg_mask_b3], int_b3[bg_mask_b3], p0=[1000, 10000, -500, 10, -0.1])
baseline_b3 = bg_model(theta_b3, *popt_b3)
net_b3 = int_b3 - baseline_b3

# 2. Fit SH-125-A Background
bg_mask_a = (theta_a < 9.5) | \
            ((theta_a > 10.5) & (theta_a < 12.3)) | \
            ((theta_a > 13.3) & (theta_a < 13.7)) | \
            ((theta_a > 15.0) & (theta_a < 21.8)) | \
            (theta_a > 23.2)

popt_a, _ = curve_fit(bg_model, theta_a[bg_mask_a], int_a[bg_mask_a], p0=[100000, 5e6, -3e5, 8000, -100])
baseline_a = bg_model(theta_a, *popt_a)
net_a = int_a - baseline_a

# Peak fitting for SH_124_B3 (Nominal 2Theta is 29.44)
two_theta_b3 = 29.44
peaks_b3_init = [9.36, 12.88, 15.72, 16.30, 22.34]
ranges_b3 = [(9.0, 9.7), (12.4, 13.4), (15.2, 16.1), (16.1, 16.7), (21.8, 22.8)]
bounds_b3 = [
    ([0, 9.1, 0.05], [100, 9.5, 0.5]),
    ([0, 12.6, 0.05], [500, 13.1, 0.4]),
    ([0, 15.5, 0.05], [500, 16.0, 0.6]),
    ([0, 16.0, 0.05], [500, 16.5, 0.5]),
    ([0, 22.1, 0.05], [200, 22.6, 0.3])
]
names_b3 = ["B3_Peak_9.22", "B3_Peak_12.92", "B3_Peak_15.70", "B3_Peak_16.19", "B3_Peak_22.34"]
results_b3 = []

for name, center, (t_min, t_max), bnd in zip(names_b3, peaks_b3_init, ranges_b3, bounds_b3):
    mask = (theta_b3 >= t_min) & (theta_b3 <= t_max)
    try:
        popt_g, _ = curve_fit(gaussian, theta_b3[mask], net_b3[mask], p0=[100, center, 0.15], bounds=bnd)
        h, t0, w = popt_g
        fwhm = 2.355 * w
        area = h * w * np.sqrt(2 * np.pi)
        iso_val = popt_b3[0] / np.sin(np.radians(t0))
        
        results_b3.append({
            'Sample': 'SH_124_B3 (2026-05-26)',
            'Peak ID': name,
            'Peak Center (Theta)': t0,
            'Tilt Angle (Chi)': t0 - two_theta_b3/2,
            'FWHM (deg)': fwhm,
            'Net Height': h,
            'Net Area (cts deg)': area,
            'Area/Base Ratio': area / iso_val
        })
    except Exception as e:
        print(f"Error fitting SH_124_B3 {name}: {e}")

# Fit minor peaks for B3
try:
    def double_gaussian_bg(t, h1, t01, w1, h2, t02, w2, c0, c1):
        return (h1 * np.exp(-(t-t01)**2 / (2*w1**2)) + 
                h2 * np.exp(-(t-t02)**2 / (2*w2**2)) + 
                c0 + c1 * (t - 17.4))

    mask_minor = (theta_b3 >= 16.7) & (theta_b3 <= 18.2)
    popt_minor, _ = curve_fit(
        double_gaussian_bg, theta_b3[mask_minor], net_b3[mask_minor],
        p0=[30, 17.10, 0.1, 30, 17.50, 0.1, 40, -10],
        bounds=(
            [1, 16.9, 0.03, 1, 17.3, 0.03, -100, -100],
            [100, 17.3, 0.3, 100, 17.7, 0.3, 200, 100]
        )
    )
    
    # Peak A
    h1, t01, w1 = popt_minor[0], popt_minor[1], popt_minor[2]
    area1 = h1 * w1 * np.sqrt(2 * np.pi)
    iso_val1 = popt_b3[0] / np.sin(np.radians(t01))
    results_b3.append({
        'Sample': 'SH_124_B3 (2026-05-26)',
        'Peak ID': 'B3_Minor_A',
        'Peak Center (Theta)': t01,
        'Tilt Angle (Chi)': t01 - two_theta_b3/2,
        'FWHM (deg)': 2.355 * w1,
        'Net Height': h1,
        'Net Area (cts deg)': area1,
        'Area/Base Ratio': area1 / iso_val1
    })
    
    # Peak B
    h2, t02, w2 = popt_minor[3], popt_minor[4], popt_minor[5]
    area2 = h2 * w2 * np.sqrt(2 * np.pi)
    iso_val2 = popt_b3[0] / np.sin(np.radians(t02))
    results_b3.append({
        'Sample': 'SH_124_B3 (2026-05-26)',
        'Peak ID': 'B3_Minor_B',
        'Peak Center (Theta)': t02,
        'Tilt Angle (Chi)': t02 - two_theta_b3/2,
        'FWHM (deg)': 2.355 * w2,
        'Net Height': h2,
        'Net Area (cts deg)': area2,
        'Area/Base Ratio': area2 / iso_val2
    })
except Exception as e:
    print(f"Error fitting SH_124_B3 minor peaks: {e}")

# Peak fitting for SH-125-A (Nominal 2Theta is 29.3425)
two_theta_a = 29.3425
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
        popt_g, _ = curve_fit(gaussian, theta_a[mask], net_a[mask], p0=[50000, center, 0.15], bounds=bnd)
        h, t0, w = popt_g
        fwhm = 2.355 * w
        area = h * w * np.sqrt(2 * np.pi)
        iso_val = popt_a[0] / np.sin(np.radians(t0))
        
        results_a.append({
            'Sample': 'SH-125-A (2026-06-01)',
            'Peak ID': name,
            'Peak Center (Theta)': t0,
            'Tilt Angle (Chi)': t0 - two_theta_a/2,
            'FWHM (deg)': fwhm,
            'Net Height': h,
            'Net Area (cts deg)': area,
            'Area/Base Ratio': area / iso_val
        })
    except Exception as e:
        print(f"Error fitting SH-125-A {name}: {e}")

# Fit minor peaks for A
try:
    def double_gaussian_bg(t, h1, t01, w1, h2, t02, w2, c0, c1):
        return (h1 * np.exp(-(t-t01)**2 / (2*w1**2)) + 
                h2 * np.exp(-(t-t02)**2 / (2*w2**2)) + 
                c0 + c1 * (t - 17.4))

    mask_minor = (theta_a >= 16.7) & (theta_a <= 18.2)
    popt_minor, _ = curve_fit(
        double_gaussian_bg, theta_a[mask_minor], net_a[mask_minor],
        p0=[35000, 17.15, 0.1, 15000, 17.55, 0.1, 5000, -2000],
        bounds=(
            [1000, 17.0, 0.03, 1000, 17.4, 0.03, -50000, -50000],
            [60000, 17.3, 0.3, 30000, 17.7, 0.3, 50000, 50000]
        )
    )
    
    # Peak A
    h1, t01, w1 = popt_minor[0], popt_minor[1], popt_minor[2]
    area1 = h1 * w1 * np.sqrt(2 * np.pi)
    iso_val1 = popt_a[0] / np.sin(np.radians(t01))
    results_a.append({
        'Sample': 'SH-125-A (2026-06-01)',
        'Peak ID': 'A_Minor_A',
        'Peak Center (Theta)': t01,
        'Tilt Angle (Chi)': t01 - two_theta_a/2,
        'FWHM (deg)': 2.355 * w1,
        'Net Height': h1,
        'Net Area (cts deg)': area1,
        'Area/Base Ratio': area1 / iso_val1
    })
    
    # Peak B
    h2, t02, w2 = popt_minor[3], popt_minor[4], popt_minor[5]
    area2 = h2 * w2 * np.sqrt(2 * np.pi)
    iso_val2 = popt_a[0] / np.sin(np.radians(t02))
    results_a.append({
        'Sample': 'SH-125-A (2026-06-01)',
        'Peak ID': 'A_Minor_B',
        'Peak Center (Theta)': t02,
        'Tilt Angle (Chi)': t02 - two_theta_a/2,
        'FWHM (deg)': 2.355 * w2,
        'Net Height': h2,
        'Net Area (cts deg)': area2,
        'Area/Base Ratio': area2 / iso_val2
    })
except Exception as e:
    print(f"Error fitting SH-125-A minor peaks: {e}")

# Combine results and save
df_compare = pd.DataFrame(results_b3 + results_a)
df_compare.to_csv(os.path.join(output_dir, "samples_peaks_comparison.csv"), index=False)
print("Saved peaks comparison table.")

# Plot comparison
fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

# Plot SH_124_B3 Residuals vs Chi
chi_b3 = theta_b3 - two_theta_b3/2
axes[0].plot(chi_b3, net_b3, 'g-', label='SH_124_B3 Residuals (Net Intensity)', alpha=0.7)
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
axes[0].set_title('Residual Orientation Peaks Comparison: SH_124_B3 vs. SH-125-A')

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
