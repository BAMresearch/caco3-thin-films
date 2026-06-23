# -*- coding: utf-8 -*-
"""
Rocking Curve Processor for 15 June 2026 Run
============================================
Author: Tomasz Stawski, tomasz.stawski@bam.de
Date: 2026-06-19
Version: 0.1

Description:
Processes and extracts peak metrics from rocking curves measured on 15 June 2026.
"""
import os
import zipfile
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d

# Define directories
base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
raw_data_dir = os.path.join(base_dir, "data/raw/Rocking_Curves/SH-104-1")
processed_dir = os.path.join(base_dir, "data/processed/Rocking_Curves/SH-104-1")
analysis_dir = os.path.join(base_dir, "data/processed/Rocking_Curves/15062026")

os.makedirs(processed_dir, exist_ok=True)
os.makedirs(analysis_dir, exist_ok=True)

# Define sample
sample = "SH-104-1"
sample_processed_dir = os.path.join(processed_dir, sample)
sample_analysis_dir = os.path.join(analysis_dir, sample)
os.makedirs(sample_processed_dir, exist_ok=True)
os.makedirs(sample_analysis_dir, exist_ok=True)

# Helper: Extract data from BRML XML file with RAW fallback
def extract_brml_data(brml_path):
    if not os.path.exists(brml_path) or os.path.getsize(brml_path) == 0:
        raw_path = brml_path.replace(".brml", ".raw")
        if os.path.exists(raw_path):
            print(f"  Warning: {os.path.basename(brml_path)} is missing/empty. Recovering data from {os.path.basename(raw_path)}...")
            is_rocking = "Rocking" in brml_path
            with open(raw_path, 'rb') as f:
                raw_bytes = f.read()
            if is_rocking:
                offset = 1488
                num_points = 1001
                intensities = np.frombuffer(raw_bytes[offset:], dtype=np.float32)
                theta = np.linspace(4.705, 24.705, num_points)
                out_arr = np.zeros((num_points, 5))
                out_arr[:, 2] = theta
                out_arr[:, 3] = intensities
                return out_arr
            else:
                offset = 1384
                num_points = 401
                intensities = np.frombuffer(raw_bytes[offset:], dtype=np.float32)
                twotheta = np.linspace(27.0001, 35.0001, num_points)
                out_arr = np.zeros((num_points, 5))
                out_arr[:, 2] = twotheta
                out_arr[:, 4] = intensities
                return out_arr
        else:
            raise FileNotFoundError(f"Neither {brml_path} nor {raw_path} could be read.")

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

# Fitting helper: Gaussian profile
def gaussian(t, h, t0, w):
    return h * np.exp(-(t-t0)**2 / (2*w**2))

# Fitting helper: Isotropic thin-film volume correction + 3rd-order polynomial background
def bg_model(t, I0, c0, c1, c2, c3):
    return I0 / np.sin(np.radians(t)) + c0 + c1*t + c2*t**2 + c3*t**3

# Peak configuration for SH-104-1 (similar to SH-125-A, as it has a similar Calcite/Vaterite signature)
peaks_config = [
    {"name": "Peak 1 (Tilt)", "init_center": 9.97, "bounds": ([0, 9.0, 0.0425], [200000, 10.8, 0.276])},
    {"name": "Peak 2a (Tilt)", "init_center": 11.70, "bounds": ([0, 11.0, 0.0425], [200000, 12.2, 0.276])},
    {"name": "Peak 2b (Tilt)", "init_center": 12.77, "bounds": ([0, 12.3, 0.0425], [200000, 13.6, 0.276])},
    {"name": "Peak 3 (Specular)", "init_center": 14.35, "bounds": ([0, 13.7, 0.0425], [200000, 15.0, 0.276])},
    {"name": "Minor Peak A (Tilt)", "init_center": 17.10, "bounds": ([0, 16.8, 0.0425], [100000, 17.3, 0.276])},
    {"name": "Minor Peak B (Tilt)", "init_center": 17.50, "bounds": ([0, 17.3, 0.0425], [100000, 17.9, 0.276])},
    {"name": "Peak 4 (Tilt)", "init_center": 22.52, "bounds": ([0, 21.5, 0.0425], [200000, 23.2, 0.276])}
]

# Detect Phi values
files = os.listdir(raw_data_dir)
phi_values = sorted(list(set([int(f.split("_")[1].split(".")[0]) for f in files if (f.startswith("2Theta_") or f.startswith("Rocking_")) and f.endswith((".brml", ".raw"))])))
print(f"Detected Phi values for {sample}: {phi_values}")

all_metrics = []
phi_rocking_data = {}
phi_twotheta_data = {}

for phi in phi_values:
    print(f"\n--- Phi = {phi}° ---")
    
    # 1. Process 2Theta scan (Phase composition & Center check)
    twotheta_file = f"2Theta_{phi}.brml"
    twotheta_path = os.path.join(raw_data_dir, twotheta_file)
    
    arr_2t = extract_brml_data(twotheta_path)
    
    # Export to processed XY
    xy_2t_path = os.path.join(sample_processed_dir, f"{sample}_2Theta_{phi}_exported.xy")
    with open(xy_2t_path, 'w') as f:
        f.write(f'Id: "{sample}" Phi: "{phi}" Scantype: "2Theta-Theta" Anode: "Cu" Wavelength: "1.5406"\n')
        for row in arr_2t:
            f.write(f"{row[2]:.5f} {row[4]:.3f}\n")
            
    twotheta_2t = arr_2t[:, 2]
    intensity_2t = arr_2t[:, 4]
    
    # Fit polynomial baseline
    poly_coeff = np.polyfit(twotheta_2t, intensity_2t, 3)
    baseline_2t = np.polyval(poly_coeff, twotheta_2t)
    net_intensity_2t = intensity_2t - baseline_2t
    
    # Fit Calcite (104) peak
    c_mask = (twotheta_2t >= 28.5) & (twotheta_2t <= 30.5)
    h_c, t0_c, fwhm_c, area_c = 0.0, 29.4, 0.0, 0.0
    try:
        popt_c, _ = curve_fit(gaussian, twotheta_2t[c_mask], net_intensity_2t[c_mask], p0=[intensity_2t.max() - baseline_2t.max(), 29.4, 0.15])
        h_c, t0_c, w_c = popt_c
        fwhm_c = 2.355 * w_c
        area_c = h_c * w_c * np.sqrt(2 * np.pi)
    except Exception as e:
        print(f"  Warning: Calcite 2Theta fit failed for Phi={phi}: {e}")
        
    # Fit Vaterite (110) peak
    h_v, t0_v, fwhm_v, area_v = 0.0, 32.8, 0.0, 0.0
    vaterite_detected = False
    v_mask = (twotheta_2t >= 31.8) & (twotheta_2t <= 33.8)
    try:
        popt_v, _ = curve_fit(gaussian, twotheta_2t[v_mask], net_intensity_2t[v_mask], p0=[10000, 32.8, 0.15])
        h_v_fit, t0_v_fit, w_v_fit = popt_v
        fwhm_v_fit = 2.355 * w_v_fit
        area_v_fit = h_v_fit * w_v_fit * np.sqrt(2 * np.pi)
        
        # Check if the fit is valid (not fitting background noise or showing weird values)
        if h_v_fit > 0.02 * h_c and 32.4 < t0_v_fit < 33.2 and fwhm_v_fit < 1.0:
            h_v, t0_v, fwhm_v, area_v = h_v_fit, t0_v_fit, fwhm_v_fit, area_v_fit
            vaterite_detected = True
    except Exception as e:
        pass
        
    c_v_ratio = area_c / area_v if area_v > 0 else np.nan
    
    print(f"  2Theta Phase Analysis:")
    print(f"    Calcite (104) Peak: center={t0_c:.4f}°, height={h_c:.1f}, FWHM={fwhm_c:.4f}°, area={area_c:.1f}")
    if vaterite_detected:
        print(f"    Vaterite (110) Peak: center={t0_v:.4f}°, height={h_v:.1f}, FWHM={fwhm_v:.4f}°, area={area_v:.1f}")
        print(f"    Calcite/Vaterite Area Ratio: {c_v_ratio:.3f}")
    else:
        print(f"    Vaterite (110) Peak: Not detected / negligible")
        
    phi_twotheta_data[phi] = {
        "2theta": twotheta_2t,
        "intensity": intensity_2t,
        "baseline": baseline_2t,
        "net": net_intensity_2t,
        "calcite": (h_c, t0_c, fwhm_c, area_c),
        "vaterite": (h_v, t0_v, fwhm_v, area_v),
        "vaterite_detected": vaterite_detected
    }
    
    # 2. Process Rocking Curve
    rocking_file = f"Rocking_{phi}.brml"
    rocking_path = os.path.join(raw_data_dir, rocking_file)
    
    arr_rc = extract_brml_data(rocking_path)
    
    # Export to processed XY
    xy_rc_path = os.path.join(sample_processed_dir, f"{sample}_Rocking_{phi}_exported.xy")
    with open(xy_rc_path, 'w') as f:
        f.write(f'Id: "{sample}" Phi: "{phi}" Scantype: "Theta (Rocking)" Anode: "Cu" Wavelength: "1.5406"\n')
        for row in arr_rc:
            f.write(f"{row[2]:.5f} {row[3]:.3f}\n")
            
    theta = arr_rc[:, 2]
    intensity = arr_rc[:, 3]
    
    # Background fitting mask
    bg_mask = (theta >= 4.0) & ((theta < 9.0) | ((theta > 18.5) & (theta < 21.5)) | (theta > 23.5))
    
    # Iteratively fit background to reject outliers
    popt_bg = [100000, 5e6, -3e5, 8000, -100]
    for iteration in range(3):
        try:
            popt_bg, _ = curve_fit(bg_model, theta[bg_mask], intensity[bg_mask], p0=popt_bg)
            baseline = bg_model(theta, *popt_bg)
            residuals = intensity - baseline
            noise = np.std(residuals[bg_mask])
            bg_mask = ((theta >= 4.0) & ((theta < 9.0) | ((theta > 18.5) & (theta < 21.5)) | (theta > 23.5))) & (residuals < 2.5 * noise)
        except Exception as e:
            # Fallback to polynomial
            poly_coeff = np.polyfit(theta[bg_mask], intensity[bg_mask], 3)
            popt_bg = [0.0, poly_coeff[3], poly_coeff[2], poly_coeff[1], poly_coeff[0]]
            break
            
    I0_init, c0_init, c1_init, c2_init, c3_init = popt_bg
    
    # Calculate noise level and peak rejection threshold
    residuals_bg = intensity - bg_model(theta, *popt_bg)
    diffs = np.diff(residuals_bg[bg_mask])
    noise_std = np.std(diffs) / np.sqrt(2) if len(diffs) > 1 else 350.0
    threshold = 1.5 * noise_std
    
    # Build global model for peaks
    net_int_init = intensity - bg_model(theta, *popt_bg)
    flat_guesses = []
    bounds_min = [I0_init*0.5 if I0_init > 0 else -1e5, c0_init - 1e6, c1_init - 1e5, c2_init - 1e4, c3_init - 1e3]
    bounds_max = [I0_init*2.0 if I0_init > 0 else 1e5, c0_init + 1e6, c1_init + 1e5, c2_init + 1e4, c3_init + 1e3]
    
    for p in peaks_config:
        c_min, c_max = p["bounds"][0][1], p["bounds"][1][1]
        p_mask = (theta >= c_min) & (theta <= c_max)
        if p_mask.any():
            center_guess = theta[p_mask][np.argmax(net_int_init[p_mask])]
            h_guess = max(net_int_init[p_mask])
        else:
            center_guess = p["init_center"]
            h_guess = 100.0
            
        h_guess = max(h_guess, 100.0)
        flat_guesses.extend([h_guess, center_guess, 0.15])
        
        bounds_min.extend(p["bounds"][0])
        bounds_max.extend(p["bounds"][1])
        
    p0 = [I0_init, c0_init, c1_init, c2_init, c3_init] + flat_guesses
    bounds = (bounds_min, bounds_max)
    
    # Fit in log10 space
    log_intensity = np.log10(np.clip(intensity, 1.0, None))
    
    def log_fit_func(t, I0, c0, c1, c2, c3, *peak_params):
        baseline_val = I0 / np.sin(np.radians(t)) + c0 + c1*t + c2*t**2 + c3*t**3
        val = baseline_val.copy()
        num_peaks = len(peak_params) // 3
        for i in range(num_peaks):
            h = peak_params[3*i]
            t0 = peak_params[3*i+1]
            w = peak_params[3*i+2]
            val += h * np.exp(-(t-t0)**2 / (2*w**2))
        return np.log10(np.clip(val, 1.0, None))
        
    try:
        popt, _ = curve_fit(log_fit_func, theta, log_intensity, p0=p0, bounds=bounds)
        I0_fit, c0_fit, c1_fit, c2_fit, c3_fit = popt[:5]
        baseline = I0_fit / np.sin(np.radians(theta)) + c0_fit + c1_fit*theta + c2_fit*theta**2 + c3_fit*theta**3
        
        peak_fits = []
        peak_params = popt[5:]
        for i, p in enumerate(peaks_config):
            h = peak_params[3*i]
            t0 = peak_params[3*i+1]
            w = peak_params[3*i+2]
            fwhm = 2.355 * w
            area = h * w * np.sqrt(2 * np.pi)
            
            if h < threshold:
                h, fwhm, area = 0.0, 0.0, 0.0
                
            if I0_fit > 0:
                iso_val = I0_fit / np.sin(np.radians(t0))
                ratio = area / iso_val
            else:
                ratio = 0.0
                
            peak_fits.append({
                "name": p["name"],
                "center": t0,
                "tilt": t0 - t0_c/2 if t0_c > 0 else t0 - 14.7,
                "fwhm": fwhm,
                "height": h,
                "area": area,
                "ratio": ratio,
                "w": w
            })
    except Exception as e:
        print(f"  Warning: Log-scale global fit failed for Phi={phi}: {e}")
        baseline = bg_model(theta, I0_init, c0_init, c1_init, c2_init, c3_init)
        I0_fit = I0_init
        peak_fits = []
        for p in peaks_config:
            peak_fits.append({
                "name": p["name"],
                "center": p["init_center"],
                "tilt": p["init_center"] - t0_c/2 if t0_c > 0 else p["init_center"] - 14.7,
                "fwhm": 0.0,
                "height": 0.0,
                "area": 0.0,
                "ratio": 0.0,
                "w": 0.15
            })
            
    net_intensity = intensity - baseline
    
    # Export corrected curve to CSV
    df_corr = pd.DataFrame({
        'Theta (degrees)': theta,
        'Raw Intensity': intensity,
        'Model Baseline': baseline,
        'Corrected Net Intensity': net_intensity
    })
    df_corr.to_csv(os.path.join(sample_processed_dir, f"{sample}_corrected_rocking_{phi}.csv"), index=False)
    
    print(f"  Rocking Curve Peak Analysis:")
    for pfit in peak_fits:
        print(f"    {pfit['name']}: center={pfit['center']:.4f}° (tilt={pfit['tilt']:.4f}°), height={pfit['height']:.1f}, area={pfit['area']:.1f}, area/base={pfit['ratio']:.5%}")
        
        all_metrics.append({
            "Sample": sample,
            "Phi (degrees)": phi,
            "Peak Name": pfit["name"],
            "Peak Center (Theta)": pfit["center"],
            "Tilt Angle (Chi)": pfit["tilt"],
            "FWHM (degrees)": pfit["fwhm"],
            "Net Height": pfit["height"],
            "Net Area (cts deg)": pfit["area"],
            "Area/Base Ratio": pfit["ratio"]
        })
        
    phi_rocking_data[phi] = {
        "theta": theta,
        "intensity": intensity,
        "baseline": baseline,
        "net": net_intensity,
        "peaks": peak_fits
    }

# Save CSV metrics
df_sample_metrics = pd.DataFrame(all_metrics)
metrics_csv_path = os.path.join(sample_analysis_dir, f"{sample}_rocking_peaks_vs_phi.csv")
df_sample_metrics.to_csv(metrics_csv_path, index=False)
print(f"\nSaved metrics CSV to {metrics_csv_path}")

# ==============================================================================
# GENERATE PLOTS
# ==============================================================================
print("\nGenerating Plots...")

# 1. 2Theta Waterfall Plot
plt.figure(figsize=(10, 8))
offset = 0.0
for phi in phi_values:
    data = phi_twotheta_data.get(phi)
    if data is not None:
        plt.plot(data["2theta"], data["intensity"] / 1e3 + offset, label=f"Phi = {phi}°")
        offset += 30.0
plt.xlabel("2Theta (degrees)")
plt.ylabel("Intensity (kcounts, offset)")
plt.title(f"2Theta Scans Stacked: {sample}")
plt.legend(loc='upper right')
plt.grid(True)
plt.tight_layout()
plot1_path = os.path.join(sample_analysis_dir, f"{sample}_2theta_scans_stacked.png")
plt.savefig(plot1_path, dpi=150)
plt.close()
print(f"Generated: {plot1_path}")

# 2. Rocking Curve Waterfall (Raw & Baseline) - Log Scale
plt.figure(figsize=(11, 9))
factor = 1.0
for phi in phi_values:
    data = phi_rocking_data.get(phi)
    if data is not None:
        plt.plot(data["theta"], data["intensity"] * factor, '.', alpha=0.5, label=f"Phi = {phi}° Raw ($\times${factor:.1e})")
        plt.plot(data["theta"], data["baseline"] * factor, 'r-', linewidth=1.5)
        factor *= 5.0
plt.yscale('log')
plt.xlabel("Theta (degrees)")
plt.ylabel("Intensity (counts, stacked log scale)")
plt.title(f"Rocking Curves and Baselines Stacked (Log Scale): {sample}")
plt.legend(loc='upper right', bbox_to_anchor=(1.15, 1.0))
plt.grid(True, which='both', linestyle=':', alpha=0.5)
plt.tight_layout()
plot2_path = os.path.join(sample_analysis_dir, f"{sample}_rocking_curves_stacked.png")
plt.savefig(plot2_path, dpi=150)
plt.close()
print(f"Generated: {plot2_path}")

# 3. Rocking Curve Peak Fits Stacked (Log Scale)
plt.figure(figsize=(11, 9))
factor = 1.0
for phi in phi_values:
    data = phi_rocking_data.get(phi)
    if data is not None:
        baseline_curve = data["baseline"]
        full_model = baseline_curve.copy()
        for p in data["peaks"]:
            if p["height"] > 0:
                full_model += gaussian(data["theta"], p["height"], p["center"], p["w"])
        
        plt.plot(data["theta"], data["intensity"] * factor, '.', alpha=0.3, label=f"Phi = {phi}° Raw ($\times${factor:.1e})")
        plt.plot(data["theta"], full_model * factor, 'k-', linewidth=1.2, label=f"Phi = {phi}° Fit" if phi == phi_values[0] else "")
        
        for p in data["peaks"]:
            if p["height"] > 0:
                peak_curve = baseline_curve + gaussian(data["theta"], p["height"], p["center"], p["w"])
                plt.plot(data["theta"], peak_curve * factor, '--', linewidth=1.0)
        factor *= 5.0
plt.yscale('log')
plt.xlabel("Theta (degrees)")
plt.ylabel("Intensity (counts, stacked log scale)")
plt.title(f"Rocking Curve Global Fits Stacked (Log Scale): {sample}")
plt.legend(loc='upper right', bbox_to_anchor=(1.15, 1.0))
plt.grid(True, which='both', linestyle=':', alpha=0.5)
plt.tight_layout()
plot3_path = os.path.join(sample_analysis_dir, f"{sample}_rocking_residuals_stacked.png")
plt.savefig(plot3_path, dpi=150)
plt.close()
print(f"Generated: {plot3_path}")

# 4. Peak Intensities (Area) vs. Phi
plt.figure(figsize=(10, 6))
unique_peaks = df_sample_metrics["Peak Name"].unique()
for peak_name in unique_peaks:
    peak_data = df_sample_metrics[df_sample_metrics["Peak Name"] == peak_name].sort_values("Phi (degrees)")
    plt.plot(peak_data["Phi (degrees)"], peak_data["Net Area (cts deg)"] / 1e3, 'o-', linewidth=2, label=peak_name)
plt.xlabel("Phi (degrees)")
plt.ylabel("Net Area (kcounts·deg)")
plt.title(f"Rocking Curve Residual Peak Areas vs. Phi: {sample}")
plt.legend()
plt.grid(True)
plt.tight_layout()
plot4_path = os.path.join(sample_analysis_dir, f"{sample}_peak_areas_vs_phi.png")
plt.savefig(plot4_path, dpi=150)
plt.close()
print(f"Generated: {plot4_path}")

# 5. Peak Tilt Angle (Chi) vs. Phi
plt.figure(figsize=(10, 6))
for peak_name in unique_peaks:
    peak_data = df_sample_metrics[df_sample_metrics["Peak Name"] == peak_name].sort_values("Phi (degrees)")
    detected_mask = peak_data["Net Area (cts deg)"] > 0
    if detected_mask.sum() > 0:
        plt.plot(peak_data.loc[detected_mask, "Phi (degrees)"], peak_data.loc[detected_mask, "Tilt Angle (Chi)"], 'o-', linewidth=2, label=peak_name)
plt.xlabel("Phi (degrees)")
plt.ylabel("Tilt Angle Chi (degrees)")
plt.title(f"Peak Tilt Angle (Chi) vs. Phi: {sample}")
plt.legend()
plt.grid(True)
plt.tight_layout()
plot5_path = os.path.join(sample_analysis_dir, f"{sample}_peak_tilt_vs_phi.png")
plt.savefig(plot5_path, dpi=150)
plt.close()
print(f"Generated: {plot5_path}")

# 6. Phase Composition and Calcite 2Theta shift vs. Phi
fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
phi_list = sorted(phi_twotheta_data.keys())
calcite_centers = [phi_twotheta_data[p]["calcite"][1] for p in phi_list]
axes[0].plot(phi_list, calcite_centers, 'go-', linewidth=2, label='Calcite (104) Center')
axes[0].set_ylabel("2Theta Position (degrees)")
axes[0].set_title(f"Calcite Peak Shift and Phase Area vs. Phi: {sample}")
axes[0].grid(True)
axes[0].legend()

calcite_areas = [phi_twotheta_data[p]["calcite"][3] / 1e3 for p in phi_list]
axes[1].plot(phi_list, calcite_areas, 'bo-', linewidth=2, label='Calcite (104) Area')
has_vaterite = any([phi_twotheta_data[p]["vaterite_detected"] for p in phi_list])
if has_vaterite:
    vaterite_areas = [phi_twotheta_data[p]["vaterite"][3] / 1e3 for p in phi_list]
    axes[1].plot(phi_list, vaterite_areas, 'mo-', linewidth=2, label='Vaterite (110) Area')
axes[1].set_xlabel("Phi (degrees)")
axes[1].set_ylabel("Peak Area (kcounts·deg)")
axes[1].grid(True)
axes[1].legend()
plt.tight_layout()
plot6_path = os.path.join(sample_analysis_dir, f"{sample}_phase_metrics_vs_phi.png")
plt.savefig(plot6_path, dpi=150)
plt.close()
print(f"Generated: {plot6_path}")

# 7 & 8. 2D Polar Texture Plots (Raw & Net)
print("Generating 2D Polar Texture Plots...")
theta_0_dict = {phi: phi_twotheta_data[phi]["calcite"][1] / 2.0 if phi_twotheta_data[phi]["calcite"][3] > 0 else 14.7 for phi in phi_values}
rocking_scans = {}
for phi in phi_values:
    csv_path = os.path.join(sample_processed_dir, f"{sample}_corrected_rocking_{phi}.csv")
    df_c = pd.read_csv(csv_path)
    theta_arr = df_c['Theta (degrees)'].values
    raw_arr = df_c['Raw Intensity'].values
    net_arr = df_c['Corrected Net Intensity'].values
    rocking_scans[phi] = {
        "raw_interp": interp1d(theta_arr, raw_arr, bounds_error=False, fill_value=0.0),
        "net_interp": interp1d(theta_arr, net_arr, bounds_error=False, fill_value=0.0)
    }

phi_polar_deg = np.arange(0, 360 + 30, 30)
phi_polar_rad = np.radians(phi_polar_deg)
r_tilt = np.linspace(0, 10.0, 500)

Z_raw = np.zeros((len(r_tilt), len(phi_polar_deg)))
Z_net = np.zeros((len(r_tilt), len(phi_polar_deg)))

for idx_phi, phi_p in enumerate(phi_polar_deg):
    phi_mapped = phi_p % 360
    raw_profile = np.zeros_like(r_tilt)
    net_profile = np.zeros_like(r_tilt)
    
    if phi_mapped <= 150:
        scan_phi = phi_mapped
        if scan_phi in rocking_scans:
            theta_0 = theta_0_dict[scan_phi]
            raw_profile = rocking_scans[scan_phi]["raw_interp"](theta_0 + r_tilt)
            net_profile = rocking_scans[scan_phi]["net_interp"](theta_0 + r_tilt)
    elif phi_mapped == 180:
        profiles_raw = []
        profiles_net = []
        if 180 in rocking_scans:
            theta_0 = theta_0_dict[180]
            profiles_raw.append(rocking_scans[180]["raw_interp"](theta_0 + r_tilt))
            profiles_net.append(rocking_scans[180]["net_interp"](theta_0 + r_tilt))
        if 0 in rocking_scans:
            theta_0 = theta_0_dict[0]
            profiles_raw.append(rocking_scans[0]["raw_interp"](theta_0 - r_tilt))
            profiles_net.append(rocking_scans[0]["net_interp"](theta_0 - r_tilt))
        if profiles_raw:
            raw_profile = np.mean(profiles_raw, axis=0)
            net_profile = np.mean(profiles_net, axis=0)
    else:
        scan_phi = phi_mapped - 180
        if scan_phi in rocking_scans:
            theta_0 = theta_0_dict[scan_phi]
            raw_profile = rocking_scans[scan_phi]["raw_interp"](theta_0 - r_tilt)
            net_profile = rocking_scans[scan_phi]["net_interp"](theta_0 - r_tilt)
            
    Z_raw[:, idx_phi] = raw_profile
    Z_net[:, idx_phi] = net_profile

Z_raw[0, :] = np.mean(Z_raw[0, :])
Z_net[0, :] = np.mean(Z_net[0, :])

for mode, Z_data in [("Raw", Z_raw), ("Net", Z_net)]:
    fig, ax = plt.subplots(subplot_kw=dict(projection='polar'), figsize=(8, 8))
    Phi_mesh, R_mesh = np.meshgrid(phi_polar_rad, r_tilt)
    Z_plot = np.clip(Z_data, 0, None) if mode == "Net" else Z_data
    contour = ax.contourf(Phi_mesh, R_mesh, Z_plot, levels=60, cmap='plasma')
    cbar = fig.colorbar(contour, ax=ax, pad=0.1, shrink=0.7)
    cbar.set_label("Intensity (counts)")
    
    ax.set_theta_zero_location("E")
    ax.set_theta_direction(1)
    
    ax.set_ylim(0, 10.0)
    ax.set_yticks([2, 4, 6, 8, 10])
    ax.set_yticklabels(["2°", "4°", "6°", "8°", "10°"], color="grey", size=10)
    ax.set_rlabel_position(45)
    
    ax.set_xticks(np.radians(np.arange(0, 360, 30)))
    ax.set_xticklabels([f"{d}°" for d in np.arange(0, 360, 30)], size=10)
    
    ax.grid(True, color="grey", linestyle=":", alpha=0.5)
    ax.set_title(f"Calcite (104) 2D Texture Plot ({mode} Intensity)\nSample: {sample}", y=1.08, fontsize=12, fontweight='bold')
    
    plot_path = os.path.join(sample_analysis_dir, f"{sample}_texture_{mode.lower()}.png")
    plt.tight_layout()
    plt.savefig(plot_path, dpi=150)
    plt.close()
    print(f"Generated: {plot_path}")

# 9. Side-by-Side Plot (Raw Stacked vs Baseline-Corrected Stacked)
print("Generating Side-by-Side Plot...")
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 7))

# Left Subplot: Raw + Baseline (Log)
factor = 1.0
for phi in phi_values:
    csv_path = os.path.join(sample_processed_dir, f"{sample}_corrected_rocking_{phi}.csv")
    df_c = pd.read_csv(csv_path)
    theta_arr = df_c["Theta (degrees)"]
    raw_intensity = df_c["Raw Intensity"]
    baseline_arr = df_c["Model Baseline"]
    
    ax1.plot(theta_arr, raw_intensity * factor, '.', alpha=0.4, label=f"{phi}° ($\times${factor:.1e})")
    ax1.plot(theta_arr, baseline_arr * factor, 'r-', linewidth=1.2)
    factor *= 5.0
    
ax1.set_yscale('log')
ax1.set_xlabel("Theta (degrees)")
ax1.set_ylabel("Intensity (counts, stacked log scale)")
ax1.set_title("Raw Rocking Curves & Baselines (Log Scale)")
ax1.grid(True, which='both', linestyle=':', alpha=0.5)
ax1.legend(loc='upper right', fontsize=8, ncol=2)

# Right Subplot: Net Corrected (Linear)
step_offset = 4000
for idx, phi in enumerate(phi_values):
    csv_path = os.path.join(sample_processed_dir, f"{sample}_corrected_rocking_{phi}.csv")
    df_c = pd.read_csv(csv_path)
    theta_arr = df_c["Theta (degrees)"]
    net_intensity = df_c["Corrected Net Intensity"]
    y_val = net_intensity + idx * step_offset
    
    ax2.plot(theta_arr, y_val, '-', linewidth=1.5, label=f"Phi = {phi}°")
    ax2.axhline(y=idx * step_offset, color='grey', linestyle='--', linewidth=0.8, alpha=0.7)
    ax2.text(theta_arr.max() + 0.2, idx * step_offset, f"{phi}°", 
             verticalalignment='center', fontsize=9, fontweight='bold')
             
ax2.set_xlabel("Theta (degrees)")
ax2.set_ylabel("Net Intensity (counts, stacked linear scale)")
ax2.set_title("Baseline-Corrected Net Curves (Linear Scale)")
ax2.grid(True, which='both', linestyle=':', alpha=0.5)
ax2.set_xlim(theta_arr.min() - 0.5, theta_arr.max() + 1.5)

fig.suptitle(f"{sample} Rocking Curve Analysis: Raw vs. Baseline-Corrected", fontsize=14, fontweight='bold')
plt.tight_layout()
side_by_side_path = os.path.join(sample_analysis_dir, f"{sample}_side_by_side.png")
plt.savefig(side_by_side_path, dpi=150)
plt.close()
print(f"Generated side-by-side plot: {side_by_side_path}")

# Copy side-by-side and texture net plots to conversation artifacts
conversation_artifacts_dir = "/home/tomek/.gemini/antigravity/brain/2974caf8-c2d8-4b77-9163-9cf57c4c82cc/artifacts"
os.makedirs(conversation_artifacts_dir, exist_ok=True)
import shutil
shutil.copy2(side_by_side_path, os.path.join(conversation_artifacts_dir, f"{sample}_side_by_side.png"))
shutil.copy2(os.path.join(sample_analysis_dir, f"{sample}_texture_net.png"), os.path.join(conversation_artifacts_dir, f"{sample}_texture_net.png"))
print("Copied primary figures to conversation artifacts.")

print("\nData processing for SH-104-1 (15062026) completed successfully!")
