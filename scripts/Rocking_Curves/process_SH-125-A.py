# -*- coding: utf-8 -*-
"""
Rocking Curve Processor for SH-125-A
====================================
Author: Tomasz Stawski, tomasz.stawski@bam.de
Date: 2026-06-19
Version: 0.1

Description:
Performs background correction and peak fitting for rocking curve sweeps of sample SH-125-A.
"""
import os
import shutil
import zipfile
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import find_peaks

# Define paths
base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
raw_dir = os.path.join(base_dir, "data/raw/Rocking_Curves/SH-125-A")
processed_dir = os.path.join(base_dir, "data/processed/Rocking_Curves/SH-125-A")
analysis_dir = os.path.join(base_dir, "data/processed/Rocking_Curves/SH-125-A")
os.makedirs(raw_dir, exist_ok=True)
os.makedirs(processed_dir, exist_ok=True)
os.makedirs(analysis_dir, exist_ok=True)

# 1. Move BRML files to target directory
brml_files = ["SH-125-A_rocking.brml", "SH-125-A_2Theta.brml"]
for f in brml_files:
    src = os.path.join(base_dir, f)
    dst = os.path.join(raw_dir, f)
    if os.path.exists(src):
        shutil.move(src, dst)
        print(f"Moved {f} to {raw_dir}")
    elif os.path.exists(dst):
        print(f"{f} already in {raw_dir}")
    else:
        print(f"Error: {f} not found!")

rocking_brml = os.path.join(raw_dir, "SH-125-A_rocking.brml")
twotheta_brml = os.path.join(raw_dir, "SH-125-A_2Theta.brml")

# 2. Extract data helper
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

# Load data arrays
rocking_arr = extract_brml_data(rocking_brml)
twotheta_arr = extract_brml_data(twotheta_brml)

# 3. Export to .xy files
rocking_xy_path = os.path.join(processed_dir, "SH-125-A_rocking-curve_01_exported.xy")
with open(rocking_xy_path, 'w') as f:
    f.write('Id: "SH-125-A" Comment: "" Operator: "Lab Manager" Anode: "Cu" Wavelength: "1.5406" Scantype: "Theta (Rocking)" TimePerStep: "350" X: "0" Y: "0" Z: "0"\n')
    for row in rocking_arr:
        f.write(f"{row[2]:.5f} {row[3]:.3f}\n")
print(f"Exported rocking curve to {rocking_xy_path}")

twotheta_xy_path = os.path.join(processed_dir, "SH-125-A_2Theta_exported.xy")
with open(twotheta_xy_path, 'w') as f:
    f.write('Id: "SH-125-A" Comment: "" Operator: "Lab Manager" Anode: "Cu" Wavelength: "1.5406" Scantype: "2Theta-Theta" TimePerStep: "1" X: "0" Y: "0" Z: "0"\n')
    for row in twotheta_arr:
        f.write(f"{row[2]:.5f} {row[4]:.3f}\n")
print(f"Exported 2Theta scan to {twotheta_xy_path}")

# 4. Perform analysis on Rocking Curve
theta = rocking_arr[:, 2]
intensity = rocking_arr[:, 3]

# Model fit to background (exclude regions with peaks to fit smooth baseline)
bg_mask = (theta < 9.5) | ((theta > 10.5) & (theta < 12.3)) | ((theta > 13.3) & (theta < 13.7)) | ((theta > 15.0) & (theta < 21.8)) | (theta > 23.2)

def bg_model(t, I0, c0, c1, c2, c3):
    return I0 / np.sin(np.radians(t)) + c0 + c1*t + c2*t**2 + c3*t**3

popt_bg, _ = curve_fit(bg_model, theta[bg_mask], intensity[bg_mask], p0=[100000, 5e6, -3e5, 8000, -100])
I0_fit = popt_bg[0]
baseline = bg_model(theta, *popt_bg)
net_intensity = intensity - baseline

# Fit Gaussian to the four main peaks in the residuals
peaks_theta = [9.97, 12.77, 14.35, 22.52]
peak_ranges = [(9.5, 10.5), (12.3, 13.3), (13.7, 15.0), (22.0, 23.0)]
peak_names = ["Peak 1 (Tilt)", "Peak 2 (Tilt)", "Peak 3 (Specular)", "Peak 4 (Tilt)"]
peak_results = []

def gaussian(t, h, t0, w):
    return h * np.exp(-(t-t0)**2 / (2*w**2))

for name, center, (t_min, t_max) in zip(peak_names, peaks_theta, peak_ranges):
    mask = (theta >= t_min) & (theta <= t_max)
    t_peak = theta[mask]
    r_peak = net_intensity[mask]
    
    try:
        popt_g, _ = curve_fit(gaussian, t_peak, r_peak, p0=[50000, center, 0.15])
        h, t0, w = popt_g
        fwhm = 2.355 * w
        area = h * w * np.sqrt(2 * np.pi)
        
        # Calculate ratio to isotropic contribution at the peak position
        iso_val = I0_fit / np.sin(np.radians(t0))
        ratio = area / iso_val
        
        peak_results.append({
            'Peak Name': name,
            'Peak Center (Theta)': t0,
            'Tilt Angle (Chi)': t0 - 29.3425/2,
            'FWHM (deg)': fwhm,
            'Net Height': h,
            'Net Area (cts deg)': area,
            'Isotropic Base': iso_val,
            'Area/Base Ratio': ratio
        })
    except Exception as e:
        print(f"Error fitting peak at {center}: {e}")

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
        p0=[35000, 17.15, 0.1, 15000, 17.55, 0.1, 5000, -2000],
        bounds=(
            [1000, 17.0, 0.03, 1000, 17.4, 0.03, -50000, -50000],
            [60000, 17.3, 0.3, 30000, 17.7, 0.3, 50000, 50000]
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
        'Tilt Angle (Chi)': t01 - 29.3425/2,
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
        'Tilt Angle (Chi)': t02 - 29.3425/2,
        'FWHM (deg)': fwhm2,
        'Net Height': h2,
        'Net Area (cts deg)': area2,
        'Isotropic Base': iso_val2,
        'Area/Base Ratio': area2 / iso_val2
    })
except Exception as e:
    print(f"Error fitting minor peaks in SH-125-A: {e}")

df_peaks = pd.DataFrame(peak_results)
df_peaks.to_csv(os.path.join(analysis_dir, "SH-125-A_rocking_peaks_metrics.csv"), index=False)
print("Saved rocking peak metrics table.")

# 5. Fit 2Theta scan peaks (Calcite 104 and Vaterite 110)
twotheta_2t = twotheta_arr[:, 2]
intensity_2t = twotheta_arr[:, 4]

# Polynomial background for 2Theta scan
poly_coeff = np.polyfit(twotheta_2t, intensity_2t, 3)
baseline_2t = np.polyval(poly_coeff, twotheta_2t)
net_intensity_2t = intensity_2t - baseline_2t

# Fit Calcite peak (around 29.34)
calcite_mask = (twotheta_2t >= 28.5) & (twotheta_2t <= 30.2)
popt_calcite, _ = curve_fit(gaussian, twotheta_2t[calcite_mask], net_intensity_2t[calcite_mask], p0=[3500, 29.34, 0.15])
h_c, t0_c, w_c = popt_calcite
fwhm_c = 2.355 * w_c
area_c = h_c * w_c * np.sqrt(2 * np.pi)

# Fit Vaterite peak (around 32.8)
vaterite_mask = (twotheta_2t >= 32.0) & (twotheta_2t <= 33.6)
popt_vaterite, _ = curve_fit(gaussian, twotheta_2t[vaterite_mask], net_intensity_2t[vaterite_mask], p0=[2000, 32.8, 0.15])
h_v, t0_v, w_v = popt_vaterite
fwhm_v = 2.355 * w_v
area_v = h_v * w_v * np.sqrt(2 * np.pi)

df_phases = pd.DataFrame([
    {'Phase': 'Calcite (104)', 'Peak Center (2Theta)': t0_c, 'FWHM (deg)': fwhm_c, 'Height (cts)': h_c, 'Area (cts deg)': area_c},
    {'Phase': 'Vaterite (110)', 'Peak Center (2Theta)': t0_v, 'FWHM (deg)': fwhm_v, 'Height (cts)': h_v, 'Area (cts deg)': area_v}
])
df_phases.to_csv(os.path.join(analysis_dir, "SH-125-A_2theta_phase_metrics.csv"), index=False)
print("Saved 2Theta phase metrics table.")

# 6. Save corrected profiles to CSV
df_corrected = pd.DataFrame({
    'Theta (degrees)': theta,
    'Raw Intensity': intensity,
    'Model Baseline': baseline,
    'Corrected Net Intensity': net_intensity
})
df_corrected.to_csv(os.path.join(processed_dir, "SH-125-A_corrected_rocking_curve.csv"), index=False)
print("Saved corrected rocking curve profiles CSV.")

# 7. Generate plots
fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
axes[0].plot(theta, intensity, 'b.', label='Experimental data', alpha=0.6)
axes[0].plot(theta, baseline, 'r-', linewidth=2, label='Model Baseline (Isotropic + Background)')
axes[0].set_ylabel('Intensity (counts)')
axes[0].legend()
axes[0].grid(True)
axes[0].set_title('SH-125-A Rocking Curve Analysis (2Theta = 29.3425°)')

axes[1].plot(theta, net_intensity, 'purple', label='Net Residual Intensity (Corrected)')
# Plot Gaussian fits
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
plot_rc_path = os.path.join(analysis_dir, "SH-125-A_rocking_curve_analysis.png")
plt.savefig(plot_rc_path, dpi=150)
plt.close()
print(f"Saved rocking curve analysis plot to {plot_rc_path}")

# Plot 2: 2Theta phase scan analysis
plt.figure(figsize=(10, 6))
plt.plot(twotheta_2t, intensity_2t, 'g-', label='Raw data', alpha=0.6)
plt.plot(twotheta_2t, baseline_2t, 'r--', label='Polynomial Baseline')
plt.plot(twotheta_2t, net_intensity_2t, 'blue', label='Net Peak Intensity')
plt.plot(twotheta_2t, gaussian(twotheta_2t, h_c, t0_c, w_c), 'm-', linewidth=2, label='Calcite Fit')
plt.plot(twotheta_2t, gaussian(twotheta_2t, h_v, t0_v, w_v), 'c-', linewidth=2, label='Vaterite Fit')

plt.text(t0_c - 0.5, h_c + 2000, f"Calcite (104)\n{t0_c:.3f}°\nFWHM={fwhm_c:.3f}°", color='purple')
plt.text(t0_v - 0.5, h_v + 2000, f"Vaterite (110)\n{t0_v:.3f}°\nFWHM={fwhm_v:.3f}°", color='teal')

plt.xlabel('2Theta (degrees)')
plt.ylabel('Intensity (counts)')
plt.title('SH-125-A 2Theta Scan Phase results')
plt.legend()
plt.grid(True)
plt.tight_layout()
plot_2t_path = os.path.join(analysis_dir, "SH-125-A_2Theta_analysis.png")
plt.savefig(plot_2t_path, dpi=150)
plt.close()
print(f"Saved 2Theta scan analysis plot to {plot_2t_path}")
