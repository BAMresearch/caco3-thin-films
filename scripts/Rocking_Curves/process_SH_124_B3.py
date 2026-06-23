# -*- coding: utf-8 -*-
"""
Rocking Curve Processor for SH-124-B3
=====================================
Author: Tomasz Stawski, tomasz.stawski@bam.de
Date: 2026-06-19
Version: 0.1

Description:
Extracts and fits multi-peak Gaussian distributions to baseline-corrected rocking curves of sample SH-124-B3.
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Define paths
base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
processed_dir = os.path.join(base_dir, "data/processed/Rocking_Curves/SH_124_B3")
analysis_dir = os.path.join(base_dir, "data/processed/Rocking_Curves/SH_124_B3")
os.makedirs(processed_dir, exist_ok=True)
os.makedirs(analysis_dir, exist_ok=True)

xy_path = os.path.join(processed_dir, "SH_124_B3_rocking-curve_01_exported.xy")

# Load data
data = []
with open(xy_path, 'r') as f:
    f.readline()  # skip header
    for line in f:
        parts = line.strip().split()
        if len(parts) == 2:
            data.append([float(parts[0]), float(parts[1])])

data = np.array(data)
theta = data[:, 0]
intensity = data[:, 1]

# 1. Background subtraction and volume correction
# Exclude peaks: 9.0-9.7, 12.2-13.6, 15.0-16.8, and low angles < 7.0
bg_mask = (theta >= 7.0) & \
          ((theta < 9.0) | (theta > 9.7)) & \
          ((theta < 12.2) | (theta > 13.6)) & \
          ((theta < 15.0) | (theta > 16.8))

def bg_model(t, I0, c0, c1, c2, c3):
    return I0 / np.sin(np.radians(t)) + c0 + c1*t + c2*t**2 + c3*t**3

popt_bg, _ = curve_fit(bg_model, theta[bg_mask], intensity[bg_mask], p0=[1000, 10000, -500, 10, -0.1])
I0_fit = popt_bg[0]
baseline = bg_model(theta, *popt_bg)
net_intensity = intensity - baseline

# 2. Fit 5 peaks in residuals
peaks_theta = [9.36, 12.88, 15.72, 16.30, 22.34]
ranges = [(9.0, 9.7), (12.4, 13.4), (15.2, 16.1), (16.1, 16.7), (21.8, 22.8)]
bounds = [
    ([0, 9.1, 0.05], [100, 9.5, 0.5]),
    ([0, 12.6, 0.05], [500, 13.1, 0.4]),
    ([0, 15.5, 0.05], [500, 16.0, 0.6]),
    ([0, 16.0, 0.05], [500, 16.5, 0.5]),
    ([0, 22.1, 0.05], [200, 22.6, 0.3])
]
peak_names = ["Peak 1 (Tilt)", "Peak 2 (Tilt)", "Peak 3 (Near-specular)", "Peak 4 (Near-specular)", "Peak 5 (Tilt)"]
two_theta = 29.44
peak_results = []

def gaussian(t, h, t0, w):
    return h * np.exp(-(t-t0)**2 / (2*w**2))

for name, center, (t_min, t_max), bnd in zip(peak_names, peaks_theta, ranges, bounds):
    mask = (theta >= t_min) & (theta <= t_max)
    try:
        popt_g, _ = curve_fit(gaussian, theta[mask], net_intensity[mask], p0=[100, center, 0.15], bounds=bnd)
        h, t0, w = popt_g
        fwhm = 2.355 * w
        area = h * w * np.sqrt(2 * np.pi)
        iso_val = I0_fit / np.sin(np.radians(t0))
        ratio = area / iso_val
        
        peak_results.append({
            'Peak Name': name,
            'Peak Center (Theta)': t0,
            'Tilt Angle (Chi)': t0 - two_theta/2,
            'FWHM (deg)': fwhm,
            'Net Height': h,
            'Net Area (cts deg)': area,
            'Isotropic Base': iso_val,
            'Area/Base Ratio': ratio
        })
    except Exception as e:
        print(f"Error fitting Peak {center}: {e}")

# Fit double Gaussian + linear background for the minor peaks in 16.7 - 18.2 region
try:
    def double_gaussian_bg(t, h1, t01, w1, h2, t02, w2, c0, c1):
        return (h1 * np.exp(-(t-t01)**2 / (2*w1**2)) + 
                h2 * np.exp(-(t-t02)**2 / (2*w2**2)) + 
                c0 + c1 * (t - 17.4))

    mask_minor = (theta >= 16.7) & (theta <= 18.2)
    t_minor = theta[mask_minor]
    y_minor = net_intensity[mask_minor]

    popt_minor, _ = curve_fit(
        double_gaussian_bg, t_minor, y_minor,
        p0=[30, 17.10, 0.1, 30, 17.50, 0.1, 40, -10],
        bounds=(
            [1, 16.9, 0.03, 1, 17.3, 0.03, -100, -100],
            [100, 17.3, 0.3, 100, 17.7, 0.3, 200, 100]
        )
    )

    # Resolve Peak A
    h1, t01, w1 = popt_minor[0], popt_minor[1], popt_minor[2]
    fwhm1 = 2.355 * w1
    area1 = h1 * w1 * np.sqrt(2 * np.pi)
    iso_val1 = I0_fit / np.sin(np.radians(t01))
    peak_results.append({
        'Peak Name': 'Minor Peak A (Tilt)',
        'Peak Center (Theta)': t01,
        'Tilt Angle (Chi)': t01 - two_theta/2,
        'FWHM (deg)': fwhm1,
        'Net Height': h1,
        'Net Area (cts deg)': area1,
        'Isotropic Base': iso_val1,
        'Area/Base Ratio': area1 / iso_val1
    })

    # Resolve Peak B
    h2, t02, w2 = popt_minor[3], popt_minor[4], popt_minor[5]
    fwhm2 = 2.355 * w2
    area2 = h2 * w2 * np.sqrt(2 * np.pi)
    iso_val2 = I0_fit / np.sin(np.radians(t02))
    peak_results.append({
        'Peak Name': 'Minor Peak B (Tilt)',
        'Peak Center (Theta)': t02,
        'Tilt Angle (Chi)': t02 - two_theta/2,
        'FWHM (deg)': fwhm2,
        'Net Height': h2,
        'Net Area (cts deg)': area2,
        'Isotropic Base': iso_val2,
        'Area/Base Ratio': area2 / iso_val2
    })
except Exception as e:
    print(f"Error fitting minor peaks in SH_124_B3: {e}")

df_peaks = pd.DataFrame(peak_results)
df_peaks.to_csv(os.path.join(analysis_dir, "SH_124_B3_rocking_peaks_metrics.csv"), index=False)
print("Saved peak metrics.")

df_corrected = pd.DataFrame({
    'Theta (degrees)': theta,
    'Raw Intensity': intensity,
    'Model Baseline': baseline,
    'Corrected Net Intensity': net_intensity
})
df_corrected.to_csv(os.path.join(processed_dir, "SH_124_B3_corrected_rocking_curve.csv"), index=False)
print("Saved corrected curve.")

# Plot rocking curve analysis
fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
axes[0].plot(theta, intensity, 'g.', label='Experimental data', alpha=0.6)
axes[0].plot(theta, baseline, 'r-', linewidth=2, label='Model Baseline (Isotropic + Background)')
axes[0].set_ylabel('Intensity (counts)')
axes[0].legend()
axes[0].grid(True)
axes[0].set_title('SH_124_B3 Rocking Curve Analysis (2Theta = 29.44°)')

axes[1].plot(theta, net_intensity, 'purple', label='Net Residual Intensity (Corrected)')
for res in peak_results:
    t0 = res['Peak Center (Theta)']
    h = res['Net Height']
    w = res['FWHM (deg)'] / 2.355
    axes[1].plot(theta, gaussian(theta, h, t0, w), 'r--', alpha=0.8)
    axes[1].text(t0 + 0.1, h * 0.9, f"{t0:.2f}°\n$\\chi$={res['Tilt Angle (Chi)']:.2f}°", fontsize=9, color='red')

axes[1].axhline(0, color='gray', linestyle='--')
axes[1].set_xlabel('Theta (degrees)')
axes[1].set_ylabel('Net Residual Intensity (counts)')
axes[1].legend()
axes[1].grid(True)
plt.tight_layout()
plot_rc_path = os.path.join(analysis_dir, "SH_124_B3_rocking_curve_analysis.png")
plt.savefig(plot_rc_path, dpi=150)
plt.close()
print(f"Saved rocking curve analysis plot to {plot_rc_path}")
