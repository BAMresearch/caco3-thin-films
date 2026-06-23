# -*- coding: utf-8 -*-
"""
Azimuthal Variation Processor
=============================
Author: Tomasz Stawski, tomasz.stawski@bam.de
Date: 2026-06-19
Version: 0.1

Description:
Fits rocking curves across all azimuthal orientations to trace the orientation dependence of diffraction intensities.
"""
import os
import shutil
import zipfile
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Define directories
base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
raw_data_dir = os.path.join(base_dir, "data/raw/Rocking_Curves")
processed_dir = os.path.join(base_dir, "data/processed/Rocking_Curves")
analysis_dir = os.path.join(base_dir, "data/processed/Rocking_Curves/08062026")

os.makedirs(processed_dir, exist_ok=True)
os.makedirs(analysis_dir, exist_ok=True)

# Helper: Extract data from BRML XML file
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

# Fitting helper: Gaussian profile
def gaussian(t, h, t0, w):
    return h * np.exp(-(t-t0)**2 / (2*w**2))

# Fitting helper: Double Gaussian + linear background for minor peaks
def double_gaussian_bg(t, h1, t01, w1, h2, t02, w2, c0, c1):
    return (h1 * np.exp(-(t-t01)**2 / (2*w1**2)) + 
            h2 * np.exp(-(t-t02)**2 / (2*w2**2)) + 
            c0 + c1 * (t - 17.4))

# Fitting helper: Isotropic thin-film volume correction + 3rd-order polynomial background
def bg_model(t, I0, c0, c1, c2, c3):
    return I0 / np.sin(np.radians(t)) + c0 + c1*t + c2*t**2 + c3*t**3

# Define samples and their Phi values
samples = ["SH-124-B3", "SH-125-A"]

# Detailed peak parameters from previous analyses
sample_configs = {
    "SH-124-B3": {
        "bg_mask_fn": lambda theta: (theta >= 7.0) & \
                                   ((theta < 7.5) | ((theta > 18.5) & (theta < 21.5)) | (theta > 23.5)),
        "peaks": [
            {"name": "Peak 1 (Tilt)", "init_center": 9.36, "bounds": ([0, 7.5, 0.0425], [200000, 10.3, 0.276])},
            {"name": "Peak 1b (Tilt)", "init_center": 10.60, "bounds": ([0, 10.4, 0.0425], [200000, 11.4, 0.276])},
            {"name": "Peak 2a (Tilt)", "init_center": 12.10, "bounds": ([0, 11.7, 0.0425], [200000, 12.3, 0.276])},
            {"name": "Peak 2b (Tilt)", "init_center": 12.88, "bounds": ([0, 12.4, 0.0425], [200000, 13.5, 0.276])},
            {"name": "Peak 3 (Near-specular)", "init_center": 15.72, "bounds": ([0, 15.0, 0.0425], [200000, 16.0, 0.276])},
            {"name": "Peak 4 (Near-specular)", "init_center": 16.30, "bounds": ([0, 16.1, 0.0425], [200000, 16.7, 0.276])},
            {"name": "Minor Peak A (Tilt)", "init_center": 17.10, "bounds": ([0, 16.8, 0.0425], [100000, 17.3, 0.276])},
            {"name": "Minor Peak B (Tilt)", "init_center": 17.50, "bounds": ([0, 17.3, 0.0425], [100000, 17.9, 0.276])},
            {"name": "Peak 5 (Tilt)", "init_center": 22.34, "bounds": ([0, 21.5, 0.0425], [200000, 23.2, 0.276])}
        ]
    },
    "SH-125-A": {
        "bg_mask_fn": lambda theta: (theta >= 4.0) & \
                                   ((theta < 9.0) | ((theta > 18.5) & (theta < 21.5)) | (theta > 23.5)),
        "peaks": [
            {"name": "Peak 1 (Tilt)", "init_center": 9.97, "bounds": ([0, 9.0, 0.0425], [200000, 10.8, 0.276])},
            {"name": "Peak 2a (Tilt)", "init_center": 11.70, "bounds": ([0, 11.7, 0.0425], [200000, 12.3, 0.276])},
            {"name": "Peak 2b (Tilt)", "init_center": 12.77, "bounds": ([0, 12.4, 0.0425], [200000, 13.6, 0.276])},
            {"name": "Peak 3 (Specular)", "init_center": 14.35, "bounds": ([0, 13.7, 0.0425], [200000, 15.0, 0.276])},
            {"name": "Minor Peak A (Tilt)", "init_center": 17.10, "bounds": ([0, 16.8, 0.0425], [100000, 17.3, 0.276])},
            {"name": "Minor Peak B (Tilt)", "init_center": 17.50, "bounds": ([0, 17.3, 0.0425], [100000, 17.9, 0.276])},
            {"name": "Peak 4 (Tilt)", "init_center": 22.52, "bounds": ([0, 21.5, 0.0425], [200000, 23.2, 0.276])}
        ]
    }
}

# Dictionary to hold all processed data and fitting metrics
all_metrics = []

for sample in samples:
    print(f"\n========================================\nProcessing Sample: {sample}\n========================================")
    sample_raw_dir = os.path.join(raw_data_dir, sample)
    sample_processed_dir = os.path.join(processed_dir, sample)
    sample_analysis_dir = os.path.join(analysis_dir, sample)
    os.makedirs(sample_processed_dir, exist_ok=True)
    os.makedirs(sample_analysis_dir, exist_ok=True)
    
    # List files to determine Phi values
    files = os.listdir(sample_raw_dir)
    phi_values = sorted(list(set([int(f.split("_")[1].split(".")[0]) for f in files if f.endswith(".brml")])))
    print(f"Detected Phi values: {phi_values}")
    
    config = sample_configs[sample]
    
    # We will accumulate data for plotting
    phi_rocking_data = {}
    phi_twotheta_data = {}
    
    for phi in phi_values:
        print(f"\n--- Phi = {phi}° ---")
        
        # 1. Process 2Theta scan (Phase composition & Center check)
        twotheta_file = f"2Theta_{phi}.brml"
        twotheta_path = os.path.join(sample_raw_dir, twotheta_file)
        
        if os.path.exists(twotheta_path):
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
            
            # Fit Vaterite (110) peak (only for SH-125-A, or try on B3 but check height)
            h_v, t0_v, fwhm_v, area_v = 0.0, 32.8, 0.0, 0.0
            vaterite_detected = False
            
            # For B3, we only check if there is an actual peak
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
        else:
            print(f"  Warning: 2Theta file {twotheta_file} not found!")
            t0_c = 29.4 # fallback default
            phi_twotheta_data[phi] = None
            
        # 2. Process Rocking Curve
        rocking_file = f"Rocking_{phi}.brml"
        rocking_path = os.path.join(sample_raw_dir, rocking_file)
        
        if os.path.exists(rocking_path):
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
            mask_fn = config["bg_mask_fn"]
            bg_mask = mask_fn(theta)
            
            # Iteratively fit background to reject outliers
            popt_bg = [100000, 5e6, -3e5, 8000, -100]
            for iteration in range(3):
                try:
                    popt_bg, _ = curve_fit(bg_model, theta[bg_mask], intensity[bg_mask], p0=popt_bg)
                    baseline = bg_model(theta, *popt_bg)
                    residuals = intensity - baseline
                    noise = np.std(residuals[bg_mask])
                    bg_mask = mask_fn(theta) & (residuals < 2.5 * noise)
                except Exception as e:
                    # Fallback to simple polynomial if it fails
                    poly_coeff = np.polyfit(theta[bg_mask], intensity[bg_mask], 3)
                    popt_bg = [0.0, poly_coeff[3], poly_coeff[2], poly_coeff[1], poly_coeff[0]]
                    break
                    
            I0_init, c0_init, c1_init, c2_init, c3_init = popt_bg
            
            # Calculate local noise level and rejection threshold
            residuals_bg = intensity - bg_model(theta, *popt_bg)
            diffs = np.diff(residuals_bg[mask_fn(theta)])
            noise_std = np.std(diffs) / np.sqrt(2) if len(diffs) > 1 else 350.0
            threshold = 1.5 * noise_std
            
            # Build global model for peaks
            peaks_list = list(config["peaks"])
            
            # Form initial guesses and bounds using local maxima of net intensity
            net_int_init = intensity - bg_model(theta, *popt_bg)
            flat_guesses = []
            bounds_min = [I0_init*0.5 if I0_init > 0 else -1e5, c0_init - 1e6, c1_init - 1e5, c2_init - 1e4, c3_init - 1e3]
            bounds_max = [I0_init*2.0 if I0_init > 0 else 1e5, c0_init + 1e6, c1_init + 1e5, c2_init + 1e4, c3_init + 1e3]
            
            for p in peaks_list:
                # Find peak position dynamically in center bounds
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
                
                # Bounds
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
                for i, p in enumerate(peaks_list):
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
                        "tilt": t0 - t0_c/2,
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
                for p in peaks_list:
                    peak_fits.append({
                        "name": p["name"],
                        "center": p["init_center"],
                        "tilt": p["init_center"] - t0_c/2,
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
                
                # Store in master metrics list
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
        else:
            print(f"  Warning: Rocking file {rocking_file} not found!")
            phi_rocking_data[phi] = None
            
    # Save results to a dataframe
    df_sample_metrics = pd.DataFrame([m for m in all_metrics if m["Sample"] == sample])
    df_sample_metrics.to_csv(os.path.join(sample_analysis_dir, f"{sample}_rocking_peaks_vs_phi.csv"), index=False)
    
    # ----------------------------------------------------
    # Generate Plots for current Sample
    # ----------------------------------------------------
    
    # Plot 1: 2Theta Waterfall Plot
    plt.figure(figsize=(10, 8))
    offset = 0.0
    for phi in phi_values:
        data = phi_twotheta_data.get(phi)
        if data is not None:
            plt.plot(data["2theta"], data["intensity"] / 1e3 + offset, label=f"Phi = {phi}°")
            offset += 30.0 # offset in thousands of counts
    plt.xlabel("2Theta (degrees)")
    plt.ylabel("Intensity (kcounts, offset)")
    plt.title(f"2Theta Scans Stacked: {sample}")
    plt.legend(loc='upper right')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(sample_analysis_dir, f"{sample}_2theta_scans_stacked.png"), dpi=150)
    plt.close()
    
    # Plot 2: Rocking Curve Waterfall (Raw & Baseline) - Log Scale
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
    plt.savefig(os.path.join(sample_analysis_dir, f"{sample}_rocking_curves_stacked.png"), dpi=150)
    plt.close()

    # Plot 3: Rocking Curve Peak Fits Stacked (Log Scale)
    plt.figure(figsize=(11, 9))
    factor = 1.0
    for phi in phi_values:
        data = phi_rocking_data.get(phi)
        if data is not None:
            # Calculate full model
            baseline_curve = data["baseline"]
            full_model = baseline_curve.copy()
            for p in data["peaks"]:
                if p["height"] > 0:
                    full_model += gaussian(data["theta"], p["height"], p["center"], p["w"])
            
            plt.plot(data["theta"], data["intensity"] * factor, '.', alpha=0.3, label=f"Phi = {phi}° Raw ($\times${factor:.1e})")
            plt.plot(data["theta"], full_model * factor, 'k-', linewidth=1.2, label=f"Phi = {phi}° Fit" if phi == phi_values[0] else "")
            
            # Plot individual peaks on top of baseline
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
    plt.savefig(os.path.join(sample_analysis_dir, f"{sample}_rocking_residuals_stacked.png"), dpi=150)
    plt.close()

    # Plot 4: Peak Intensities (Area) vs. Phi
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
    plt.savefig(os.path.join(sample_analysis_dir, f"{sample}_peak_areas_vs_phi.png"), dpi=150)
    plt.close()

    # Plot 5: Peak Tilt Angle (Chi) vs. Phi (Checking peak shifts)
    plt.figure(figsize=(10, 6))
    for peak_name in unique_peaks:
        peak_data = df_sample_metrics[df_sample_metrics["Peak Name"] == peak_name].sort_values("Phi (degrees)")
        # Only plot where the peak is actually detected
        detected_mask = peak_data["Net Area (cts deg)"] > 0
        if detected_mask.sum() > 0:
            plt.plot(peak_data.loc[detected_mask, "Phi (degrees)"], peak_data.loc[detected_mask, "Tilt Angle (Chi)"], 'o-', linewidth=2, label=peak_name)
    plt.xlabel("Phi (degrees)")
    plt.ylabel("Tilt Angle Chi (degrees)")
    plt.title(f"Peak Tilt Angle (Chi) vs. Phi: {sample}")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(sample_analysis_dir, f"{sample}_peak_tilt_vs_phi.png"), dpi=150)
    plt.close()

    # Plot 6: Phase Composition and Calcite 2Theta shift vs. Phi
    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    
    # Subplot 0: 2Theta center position
    phi_list = sorted(phi_twotheta_data.keys())
    calcite_centers = [phi_twotheta_data[p]["calcite"][1] for p in phi_list if phi_twotheta_data[p] is not None]
    axes[0].plot(phi_list, calcite_centers, 'go-', linewidth=2, label='Calcite (104) Center')
    axes[0].set_ylabel("2Theta Position (degrees)")
    axes[0].set_title(f"Calcite Peak Shift and Phase Area vs. Phi: {sample}")
    axes[0].grid(True)
    axes[0].legend()
    
    # Subplot 1: Calcite & Vaterite Peak Areas
    calcite_areas = [phi_twotheta_data[p]["calcite"][3] / 1e3 for p in phi_list if phi_twotheta_data[p] is not None]
    axes[1].plot(phi_list, calcite_areas, 'bo-', linewidth=2, label='Calcite (104) Area')
    
    has_vaterite = any([phi_twotheta_data[p]["vaterite_detected"] for p in phi_list if phi_twotheta_data[p] is not None])
    if has_vaterite:
        vaterite_areas = [phi_twotheta_data[p]["vaterite"][3] / 1e3 for p in phi_list if phi_twotheta_data[p] is not None]
        axes[1].plot(phi_list, vaterite_areas, 'mo-', linewidth=2, label='Vaterite (110) Area')
    axes[1].set_xlabel("Phi (degrees)")
    axes[1].set_ylabel("Peak Area (kcounts·deg)")
    axes[1].grid(True)
    axes[1].legend()
    
    plt.tight_layout()
    plt.savefig(os.path.join(sample_analysis_dir, f"{sample}_phase_metrics_vs_phi.png"), dpi=150)
    plt.close()

# Save master CSV
df_master_metrics = pd.DataFrame(all_metrics)
df_master_metrics.to_csv(os.path.join(analysis_dir, "all_samples_rocking_peaks_vs_phi.csv"), index=False)

# Create a comparison plot of peak area vs phi for major tilt peaks in both samples
plt.figure(figsize=(12, 7))

# B3 Peak 2a is tilt ~ -2.8, Peak 2b is tilt ~ -1.8, Peak 5 is tilt ~ +7.7
# A Peak 2a is tilt ~ -2.8, Peak 2b is tilt ~ -1.8, Peak 4 is tilt ~ +7.7
for sample in samples:
    df_s = df_master_metrics[df_master_metrics["Sample"] == sample]
    
    tilt_2_8 = df_s[df_s["Peak Name"].str.contains("Peak 2a")]
    tilt_1_8 = df_s[df_s["Peak Name"].str.contains("Peak 2b")]
    
    # +7.7 tilt is Peak 5 in B3, Peak 4 in A
    if sample == "SH-124-B3":
        tilt_7_7 = df_s[df_s["Peak Name"] == "Peak 5 (Tilt)"]
    else:
        tilt_7_7 = df_s[df_s["Peak Name"] == "Peak 4 (Tilt)"]
        
    if not tilt_2_8.empty and (tilt_2_8["Net Area (cts deg)"] > 0).any():
        # filter out zero entries for cleaner plotting
        df_plot = tilt_2_8[tilt_2_8["Net Area (cts deg)"] > 0]
        plt.plot(df_plot["Phi (degrees)"], df_plot["Net Area (cts deg)"] / 1e3, 'v--', linewidth=2, label=f"{sample} (tilt $\\approx -2.8^\\circ$)")
        
    if not tilt_1_8.empty and (tilt_1_8["Net Area (cts deg)"] > 0).any():
        df_plot = tilt_1_8[tilt_1_8["Net Area (cts deg)"] > 0]
        plt.plot(df_plot["Phi (degrees)"], df_plot["Net Area (cts deg)"] / 1e3, 'o--', linewidth=2, label=f"{sample} (tilt $\\approx -1.8^\\circ$)")
        
    if not tilt_7_7.empty and (tilt_7_7["Net Area (cts deg)"] > 0).any():
        df_plot = tilt_7_7[tilt_7_7["Net Area (cts deg)"] > 0]
        plt.plot(df_plot["Phi (degrees)"], df_plot["Net Area (cts deg)"] / 1e3, 's-', linewidth=2, label=f"{sample} (tilt $\\approx +7.7^\\circ$)")

plt.xlabel("Phi (degrees)")
plt.ylabel("Net Area (kcounts·deg)")
plt.title("Comparison of Azimuthal (Phi) Dependence of Key Tilt Components")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig(os.path.join(analysis_dir, "samples_key_peaks_comparison_vs_phi.png"), dpi=150)
plt.close()

print("\nProcessing completed successfully. All figures and CSV files generated.")
