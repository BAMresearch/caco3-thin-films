# -*- coding: utf-8 -*-
"""
Rocking Curve Processor for 09 June 2026 Run
============================================
Author: Tomasz Stawski, tomasz.stawski@bam.de
Date: 2026-06-19
Version: 0.1

Description:
Parses Bruker brml files and processes rocking curves from the 09 June 2026 measurement batch.
"""
import os
import zipfile
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Define paths
base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
raw_dir = os.path.join(base_dir, "data/raw/Rocking_Curves/Reference")
processed_dir = os.path.join(base_dir, "data/processed/Rocking_Curves/Reference")
analysis_dir = os.path.join(base_dir, "data/processed/Rocking_Curves/Reference")

os.makedirs(processed_dir, exist_ok=True)
os.makedirs(analysis_dir, exist_ok=True)

calcite_brml = os.path.join(raw_dir, "Calcite single crystal rocking curve.brml")
holder_brml = os.path.join(raw_dir, "sample holder.brml")

# Helper: Extract data from BRML XML
def extract_brml_data(brml_path):
    with zipfile.ZipFile(brml_path, 'r') as z:
        xml_content = z.read("Experiment0/RawData0.xml")
        root = ET.fromstring(xml_content)
        ns = root.tag.split('}')[0] + '}' if '}' in root.tag else ""
        datums = root.findall(f".//{ns}Datum")
        data_list = []
        for d in datums:
            parts = [float(x) for x in d.text.strip().split(',')]
            data_list.append(parts)
        return np.array(data_list)

# Fitting helpers
def bg_model(t, I0, c0, c1, c2, c3):
    return I0 / np.sin(np.radians(t)) + c0 + c1*t + c2*t**2 + c3*t**3

def gaussian(t, h, t0, w):
    return h * np.exp(-(t-t0)**2 / (2*w**2))

# ==============================================================================
# 1. PROCESS CALCITE SINGLE CRYSTAL ROCKING CURVE
# ==============================================================================
print("Processing Calcite Single Crystal Rocking Curve...")
calcite_arr = extract_brml_data(calcite_brml)

# Export raw data to XY
calcite_xy_path = os.path.join(processed_dir, "calcite_single_crystal_rocking_exported.xy")
with open(calcite_xy_path, 'w') as f:
    f.write('Id: "Calcite Single Crystal" Comment: "" Operator: "Lab Manager" Anode: "Cu" Wavelength: "1.5406" Scantype: "Theta (Rocking)" TimePerStep: "2000"\n')
    for row in calcite_arr:
        f.write(f"{row[2]:.5f} {row[3]:.3f}\n")
print(f"Exported raw rocking curve to {calcite_xy_path}")

theta_c = calcite_arr[:, 2]
int_c = calcite_arr[:, 3]

# Fit background: exclude peak region (10.0 to 22.0)
bg_mask_c = (theta_c < 10.0) | (theta_c > 22.0)
popt_bg_c, _ = curve_fit(bg_model, theta_c[bg_mask_c], int_c[bg_mask_c], p0=[100000, 5e5, -1e4, 500, -10])
baseline_c = bg_model(theta_c, *popt_bg_c)
net_int_c = int_c - baseline_c

# Fit Gaussian peak
max_val_c = net_int_c.max()
max_theta_c = theta_c[np.argmax(net_int_c)]

# Find FWHM estimate
half_max_c = max_val_c / 2.0
idx_above_c = np.where(net_int_c >= half_max_c)[0]
fwhm_est_c = theta_c[idx_above_c[-1]] - theta_c[idx_above_c[0]]

popt_g_c, _ = curve_fit(gaussian, theta_c[(theta_c >= 10.0) & (theta_c <= 22.0)], net_int_c[(theta_c >= 10.0) & (theta_c <= 22.0)], p0=[max_val_c, max_theta_c, fwhm_est_c / 2.355])
h_c, t0_c, w_c = popt_g_c
fwhm_c = 2.355 * w_c
area_c = h_c * w_c * np.sqrt(2 * np.pi)

# Isotropic contribution at the peak position
I0_fit_c = popt_bg_c[0]
iso_val_c = I0_fit_c / np.sin(np.radians(t0_c))
ratio_c = area_c / iso_val_c

# Save metrics to CSV
df_metrics_c = pd.DataFrame([{
    'Peak Name': 'Calcite (104) Single Crystal Peak',
    'Peak Center (Theta)': t0_c,
    'Tilt Angle (Chi)': t0_c - 29.3725/2, # Nominal 2Theta = 29.3725 from raw metadata
    'FWHM (deg)': fwhm_c,
    'Net Height': h_c,
    'Net Area (cts deg)': area_c,
    'Isotropic Base': iso_val_c,
    'Area/Base Ratio': ratio_c
}])
metrics_c_path = os.path.join(analysis_dir, "calcite_single_crystal_rocking_peaks_metrics.csv")
df_metrics_c.to_csv(metrics_c_path, index=False)
print(f"Saved rocking peak metrics to {metrics_c_path}")

# Save corrected data to CSV
df_corr_c = pd.DataFrame({
    'Theta (degrees)': theta_c,
    'Raw Intensity': int_c,
    'Model Baseline': baseline_c,
    'Corrected Net Intensity': net_int_c
})
corr_c_path = os.path.join(processed_dir, "calcite_single_crystal_corrected_rocking.csv")
df_corr_c.to_csv(corr_c_path, index=False)
print(f"Saved corrected rocking curve profiles to {corr_c_path}")

# ==============================================================================
# 2. PROCESS SAMPLE HOLDER 2THETA-THETA SCAN
# ==============================================================================
print("Processing Sample Holder 2Theta-Theta Scan...")
holder_arr = extract_brml_data(holder_brml)

# Export raw data to XY
holder_xy_path = os.path.join(processed_dir, "sample_holder_2theta_exported.xy")
with open(holder_xy_path, 'w') as f:
    f.write('Id: "Sample Holder" Comment: "" Operator: "Lab Manager" Anode: "Cu" Wavelength: "1.5406" Scantype: "2Theta-Theta" TimePerStep: "20000"\n')
    for row in holder_arr:
        f.write(f"{row[2]:.5f} {row[4]:.3f}\n")
print(f"Exported raw 2Theta-Theta scan to {holder_xy_path}")

twotheta_h = holder_arr[:, 2]
int_h = holder_arr[:, 4]

# Fit 3rd-order polynomial baseline
poly_coeff_h = np.polyfit(twotheta_h, int_h, 3)
baseline_h = np.polyval(poly_coeff_h, twotheta_h)
net_int_h = int_h - baseline_h

# Save corrected data to CSV
df_corr_h = pd.DataFrame({
    '2Theta (degrees)': twotheta_h,
    'Raw Intensity': int_h,
    'Polynomial Baseline': baseline_h,
    'Corrected Net Intensity': net_int_h
})
corr_h_path = os.path.join(processed_dir, "sample_holder_corrected_2theta.csv")
df_corr_h.to_csv(corr_h_path, index=False)
print(f"Saved corrected sample holder profiles to {corr_h_path}")

# ==============================================================================
# 3. GENERATE PLOTS
# ==============================================================================
print("Generating Plots...")

# Plot 1: Calcite Single Crystal Rocking Curve Analysis
fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
axes[0].plot(theta_c, int_c, '.', color='#1f77b4', label='Experimental data', alpha=0.6)
axes[0].plot(theta_c, baseline_c, '-', color='#d62728', linewidth=2, label='Model Baseline (Isotropic + Background)')
axes[0].set_ylabel('Intensity (counts)')
axes[0].legend()
axes[0].grid(True, linestyle='--', alpha=0.5)
axes[0].set_title('Calcite Single Crystal Rocking Curve Analysis (2Theta = 29.3725°)')

axes[1].plot(theta_c, net_int_c, color='#9467bd', label='Net Residual Intensity (Corrected)')
axes[1].plot(theta_c, gaussian(theta_c, h_c, t0_c, w_c), '--', color='#d62728', linewidth=1.5, label='Gaussian Fit')
axes[1].text(t0_c + 0.2, h_c * 0.9, f"Center = {t0_c:.3f}°\nFWHM = {fwhm_c:.3f}°\n$\\chi$ = {t0_c - 29.3725/2:.3f}°", fontsize=10, color='red', bbox=dict(facecolor='white', alpha=0.8, boxstyle='round,pad=0.3'))
axes[1].axhline(0, color='gray', linestyle='--')
axes[1].set_xlabel('Theta (degrees)')
axes[1].set_ylabel('Net Residual Intensity (counts)')
axes[1].legend()
axes[1].grid(True, linestyle='--', alpha=0.5)

plt.tight_layout()
plot_rc_path = os.path.join(analysis_dir, "calcite_single_crystal_rocking_curve_analysis.png")
plt.savefig(plot_rc_path, dpi=150)
plt.close()
print(f"Saved Calcite RC plot to {plot_rc_path}")

# Plot 2: Sample Holder 2Theta-Theta Scan Analysis
plt.figure(figsize=(10, 6))
plt.plot(twotheta_h, int_h, '.', color='#2ca02c', label='Raw data', alpha=0.6)
plt.plot(twotheta_h, baseline_h, '-', color='#d62728', linewidth=2, label='Polynomial Baseline')
plt.plot(twotheta_h, net_int_h, color='#1f77b4', label='Net Peak Intensity (Residuals)')
plt.axhline(0, color='gray', linestyle='--')
plt.xlabel('2Theta (degrees)')
plt.ylabel('Intensity (counts)')
plt.title('Sample Holder 2Theta-Theta Background Scan')
plt.legend()
plt.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()
plot_sh_path = os.path.join(analysis_dir, "sample_holder_2theta_analysis.png")
plt.savefig(plot_sh_path, dpi=150)
plt.close()
print(f"Saved Sample Holder plot to {plot_sh_path}")

print("All processing and plotting tasks completed successfully!")
