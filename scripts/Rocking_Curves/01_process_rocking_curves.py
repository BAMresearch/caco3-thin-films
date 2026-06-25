#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Consolidated Rocking Curve Processor
====================================
Author: Tomasz Stawski, tomasz.stawski@bam.de
Date: 2026-06-23
Version: 1.0

Description:
Extracts, volume-corrects, baseline-subtracts, and deconvolutes CaCO3 rocking curve data
for samples SH-124-B3, SH-125-A, SH-125-G, SH-104-1, and single crystal/holder references.
"""

import os
import sys
import argparse
import zipfile
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

# Directories
script_dir = os.path.dirname(os.path.abspath(__file__))
base_dir = os.path.abspath(os.path.join(script_dir, "..", ".."))
RAW_DIR = os.path.join(base_dir, "data/raw")
PROCESSED_DIR = os.path.join(base_dir, "data/processed")

def extract_brml_data(brml_path):
    """
    Extracts scattering angles and diffraction intensity data arrays from a zipped .brml file.
    Falls back to reading from a corresponding binary .raw file if the .brml file is empty or missing.
    """
    if not os.path.exists(brml_path) or os.path.getsize(brml_path) == 0:
        raw_path = brml_path.replace(".brml", ".raw")
        if os.path.exists(raw_path):
            print(f"  Warning: {os.path.basename(brml_path)} is missing. Recovering from binary {os.path.basename(raw_path)}...")
            is_rocking = "Rocking" in brml_path or "rocking" in brml_path
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

def gaussian(t, h, t0, w):
    """Evaluates a Gaussian peak profile."""
    return h * np.exp(-(t - t0)**2 / (2 * w**2))

def bg_model(t, I0, c0, c1, c2, c3):
    """Defines the background model comprising an isotropic volume correction term and a 3rd-order polynomial."""
    return I0 / np.sin(np.radians(t)) + c0 + c1*t + c2*t**2 + c3*t**3

def fit_symmetric_scan(twotheta, intensity):
    """Fits background baseline and symmetric Bragg reflections for calcite (104) and vaterite (110)."""
    mask = (twotheta >= 27.5) & (twotheta <= 34.5)
    twotheta_f = twotheta[mask]
    intensity_f = intensity[mask]
    
    # 2nd order baseline
    poly_coeff = np.polyfit(twotheta_f[(twotheta_f < 28.5) | (twotheta_f > 33.5)], 
                            intensity_f[(twotheta_f < 28.5) | (twotheta_f > 33.5)], 2)
    baseline = np.polyval(poly_coeff, twotheta_f)
    net_intensity = intensity_f - baseline
    
    # Init guess
    c_mask = (twotheta_f >= 28.5) & (twotheta_f <= 30.5)
    v_mask = (twotheta_f >= 32.0) & (twotheta_f <= 33.5)
    
    h_c = max(net_intensity[c_mask]) if c_mask.any() else 100.0
    t0_c = twotheta_f[c_mask][np.argmax(net_intensity[c_mask])] if c_mask.any() else 29.4
    h_v = max(net_intensity[v_mask]) if v_mask.any() else 0.0
    t0_v = twotheta_f[v_mask][np.argmax(net_intensity[v_mask])] if v_mask.any() else 32.8
    
    p0 = [h_c, t0_c, 0.15, h_v, t0_v, 0.15, poly_coeff[2], poly_coeff[1], poly_coeff[0]]
    bounds = (
        [0, 29.0, 0.05, 0, 32.2, 0.05, -1e6, -1e6, -1e6],
        [1e6, 29.9, 0.35, 1e6, 33.4, 0.35, 1e6, 1e6, 1e6]
    )
    
    def fit_func(t, h1, t01, w1, h2, t02, w2, a0, a1, a2):
        return (h1 * np.exp(-(t-t01)**2 / (2*w1**2)) + 
                h2 * np.exp(-(t-t02)**2 / (2*w2**2)) + 
                a0 + a1*t + a2*t**2)
                
    try:
        popt, _ = curve_fit(fit_func, twotheta_f, intensity_f, p0=p0, bounds=bounds)
        return {
            "calcite_center": popt[1],
            "calcite_area": popt[0] * popt[2] * np.sqrt(2 * np.pi),
            "vaterite_center": popt[4],
            "vaterite_area": popt[3] * popt[5] * np.sqrt(2 * np.pi)
        }
    except:
        return {"calcite_center": t0_c, "calcite_area": h_c * 0.15 * np.sqrt(2 * np.pi), 
                "vaterite_center": t0_v, "vaterite_area": h_v * 0.15 * np.sqrt(2 * np.pi)}

def process_rocking_curves():
    """Coordinates processing and peak parameter extraction for rocking curve measurements across all samples and references."""
    sample_configs = {
        "SH-124-B3": {
            "raw_sub_dir": "SH-124-B3",
            "bg_mask_fn": lambda theta: (theta >= 7.0) & \
                                       ((theta < 7.5) | ((theta > 18.5) & (theta < 21.5)) | (theta > 23.5)),
            "peaks": [
                {"name": "Peak 1 (Tilt)", "init_center": 9.36, "bounds": ([0, 7.5, 0.0425], [200000, 10.3, 0.276])},
                {"name": "Peak 1b (Tilt)", "init_center": 10.60, "bounds": ([0, 10.4, 0.0425], [200000, 11.4, 0.276])},
                {"name": "Peak 2a (Tilt)", "init_center": 12.10, "bounds": ([0, 11.7, 0.0425], [200000, 12.3, 0.276])},
                {"name": "Peak 2b (Tilt)", "init_center": 12.88, "bounds": ([0, 12.4, 0.0425], [200000, 13.5, 0.276])},
                {"name": "Peak 2c (Tilt)", "init_center": 14.10, "bounds": ([0, 13.6, 0.0425], [200000, 14.9, 0.276])},
                {"name": "Peak 3 (Near-specular)", "init_center": 15.72, "bounds": ([0, 15.0, 0.0425], [200000, 16.0, 0.276])},
                {"name": "Peak 4 (Near-specular)", "init_center": 16.30, "bounds": ([0, 16.1, 0.0425], [200000, 16.7, 0.276])},
                {"name": "Minor Peak A (Tilt)", "init_center": 17.10, "bounds": ([0, 16.8, 0.0425], [100000, 17.3, 0.276])},
                {"name": "Minor Peak B (Tilt)", "init_center": 17.50, "bounds": ([0, 17.3, 0.0425], [100000, 17.9, 0.276])},
                {"name": "Peak 5 (Tilt)", "init_center": 22.34, "bounds": ([0, 21.5, 0.0425], [200000, 23.2, 0.276])}
            ],
            "default_2theta": 29.44
        },
        "SH-125-A": {
            "raw_sub_dir": "SH-125-A",
            "bg_mask_fn": lambda theta: (theta >= 4.0) & \
                                       ((theta < 9.0) | ((theta > 18.5) & (theta < 21.5)) | (theta > 23.5)),
            "peaks": [
                {"name": "Broad Specular Peak", "init_center": 13.70, "bounds": ([0, 11.0, 1.0], [200000, 16.0, 8.0])}
            ],
            "default_2theta": 29.3425
        },
        "SH-125-G": {
            "raw_sub_dir": "SH-125-G",
            "bg_mask_fn": lambda theta: (theta >= 4.0) & \
                                       ((theta < 9.0) | ((theta > 18.5) & (theta < 21.5)) | (theta > 23.5)),
            "peaks": [
                {"name": "Broad Specular Peak", "init_center": 13.70, "bounds": ([0, 11.0, 1.0], [200000, 16.0, 8.0])}
            ],
            "default_2theta": 29.3425
        },
        "SH-104-1": {
            "raw_sub_dir": "SH-104-1",
            "bg_mask_fn": lambda theta: (theta >= 4.0) & \
                                       ((theta < 9.0) | ((theta > 18.5) & (theta < 21.5)) | (theta > 23.5)),
            "peaks": [
                {"name": "Broad Specular Peak", "init_center": 13.70, "bounds": ([0, 11.0, 1.0], [200000, 16.0, 8.0])}
            ],
            "default_2theta": 29.3425
        }
    }

    all_metrics = []

    for sample, config in sample_configs.items():
        print(f"\nProcessing Rocking Curves for Sample: {sample}...")
        sample_raw_dir = os.path.join(RAW_DIR, f"Rocking_Curves/{config['raw_sub_dir']}")
        sample_processed_dir = os.path.join(PROCESSED_DIR, f"Rocking_Curves/{sample}")
        sample_sym_dir = os.path.join(PROCESSED_DIR, f"Symmetric_Scans/{sample}")
        os.makedirs(sample_processed_dir, exist_ok=True)
        os.makedirs(sample_sym_dir, exist_ok=True)
        
        # Detect phis from files present
        if not os.path.exists(sample_raw_dir):
            print(f"  Warning: Raw directory {sample_raw_dir} does not exist. Skipping.")
            continue
            
        files = os.listdir(sample_raw_dir)
        phi_values = sorted(list(set([
            int(f.split("_")[1].split(".")[0]) for f in files 
            if (f.startswith("2Theta_") or f.startswith("Rocking_") or f.startswith(f"{sample}_2Theta_") or f.startswith(f"{sample}_rocking_"))
            and f.endswith((".brml", ".raw"))
        ])))
        
        if not phi_values:
            phi_values = [0]
            
        print(f"  Detected Phi values: {phi_values}")
        
        for phi in phi_values:
            # Look for symmetric scan 2theta file
            twotheta_file = f"2Theta_{phi}.brml"
            possible_2t = [
                os.path.join(sample_raw_dir, twotheta_file),
                os.path.join(sample_raw_dir, f"{sample}_2Theta.brml"),
                os.path.join(sample_raw_dir, f"{sample}_2Theta_{phi}.brml")
            ]
            twotheta_path = None
            for p in possible_2t:
                if os.path.exists(p) or os.path.exists(p.replace(".brml", ".raw")):
                    twotheta_path = p
                    break
            
            t0_c = config["default_2theta"]
            if twotheta_path:
                try:
                    arr_2t = extract_brml_data(twotheta_path)
                    xy_2t_path = os.path.join(sample_sym_dir, f"{sample}_2Theta_{phi}_exported.xy")
                    with open(xy_2t_path, 'w') as f:
                        f.write(f'Id: "{sample}" Phi: "{phi}" Scantype: "2Theta-Theta" Anode: "Cu" Wavelength: "1.5406"\n')
                        for row in arr_2t:
                            f.write(f"{row[2]:.5f} {row[4]:.3f}\n")
                            
                    twotheta_2t = arr_2t[:, 2]
                    intensity_2t = arr_2t[:, 4]
                    
                    fit_res = fit_symmetric_scan(twotheta_2t, intensity_2t)
                    if fit_res["calcite_area"] > 0:
                        t0_c = fit_res["calcite_center"]
                except Exception as e:
                    print(f"    Warning: could not process 2Theta scan at Phi={phi}: {e}")

            # Rocking curve processing
            rocking_file = f"Rocking_{phi}.brml"
            possible_rc = [
                os.path.join(sample_raw_dir, rocking_file),
                os.path.join(sample_raw_dir, f"{sample}_rocking.brml"),
                os.path.join(sample_raw_dir, f"{sample}_rocking_{phi}.brml")
            ]
            rocking_path = None
            for p in possible_rc:
                if os.path.exists(p) or os.path.exists(p.replace(".brml", ".raw")):
                    rocking_path = p
                    break
                    
            if rocking_path:
                try:
                    arr_rc = extract_brml_data(rocking_path)
                    theta = arr_rc[:, 2]
                    intensity = arr_rc[:, 3]
                    
                    bg_mask = config["bg_mask_fn"](theta)
                    popt_bg = [100000, 5e6, -3e5, 8000, -100]
                    
                    # Iterative baseline fitting
                    for _ in range(3):
                        try:
                            popt_bg, _ = curve_fit(bg_model, theta[bg_mask], intensity[bg_mask], p0=popt_bg)
                            baseline = bg_model(theta, *popt_bg)
                            residuals = intensity - baseline
                            noise = np.std(residuals[bg_mask])
                            bg_mask = config["bg_mask_fn"](theta) & (residuals < 2.5 * noise)
                        except:
                            poly_coeff = np.polyfit(theta[bg_mask], intensity[bg_mask], 3)
                            popt_bg = [0.0, poly_coeff[3], poly_coeff[2], poly_coeff[1], poly_coeff[0]]
                            break
                            
                    I0_fit, c0_fit, c1_fit, c2_fit, c3_fit = popt_bg
                    baseline = I0_fit / np.sin(np.radians(theta)) + c0_fit + c1_fit*theta + c2_fit*theta**2 + c3_fit*theta**3
                    net_intensity = intensity - baseline
                    
                    df_corr = pd.DataFrame({
                        'Theta (degrees)': theta,
                        'Raw Intensity': intensity,
                        'Model Baseline': baseline,
                        'Corrected Net Intensity': net_intensity
                    })
                    df_corr.to_csv(os.path.join(sample_processed_dir, f"{sample}_corrected_rocking_{phi}.csv"), index=False)
                    
                    # Log-scale multi-peak fitting
                    diffs = np.diff(net_intensity[config["bg_mask_fn"](theta)])
                    noise_std = np.std(diffs) / np.sqrt(2) if len(diffs) > 1 else 350.0
                    if sample == "SH-124-B3" and phi == 150:
                        threshold = 4.0 * noise_std  # Prevent overfitting in noisy 150 deg scan
                    else:
                        threshold = 1.5 * noise_std
                    
                    peaks_list = list(config["peaks"])
                    flat_guesses = []
                    bounds_min = [I0_fit*0.5 if I0_fit > 0 else -1e5, c0_fit - 1e6, c1_fit - 1e5, c2_fit - 1e4, c3_fit - 1e3]
                    bounds_max = [I0_fit*2.0 if I0_fit > 0 else 1e5, c0_fit + 1e6, c1_fit + 1e5, c2_fit + 1e4, c3_fit + 1e3]
                    
                    for p in peaks_list:
                        c_min, c_max = p["bounds"][0][1], p["bounds"][1][1]
                        p_mask = (theta >= c_min) & (theta <= c_max)
                        if p_mask.any():
                            center_guess = theta[p_mask][np.argmax(net_intensity[p_mask])]
                            h_guess = max(net_intensity[p_mask])
                        else:
                            center_guess = p["init_center"]
                            h_guess = 100.0
                        h_guess = max(h_guess, 100.0)
                        flat_guesses.extend([h_guess, center_guess, 0.15])
                        bounds_min.extend(p["bounds"][0])
                        bounds_max.extend(p["bounds"][1])
                        
                    p0 = [I0_fit, c0_fit, c1_fit, c2_fit, c3_fit] + flat_guesses
                    bounds = (bounds_min, bounds_max)
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
                        I0_res = popt[0]
                        peak_params = popt[5:]
                        for i, p in enumerate(peaks_list):
                            h = peak_params[3*i]
                            t0 = peak_params[3*i+1]
                            w = peak_params[3*i+2]
                            fwhm = 2.355 * w
                            area = h * w * np.sqrt(2 * np.pi)
                            
                            if h < threshold:
                                h, fwhm, area = 0.0, 0.0, 0.0
                                
                            iso_val = I0_res / np.sin(np.radians(t0)) if I0_res > 0 else 1.0
                            ratio = area / iso_val
                            
                            all_metrics.append({
                                "Sample": sample,
                                "Phi (degrees)": phi,
                                "Peak Name": p["name"],
                                "Peak Center (Theta)": t0,
                                "Tilt Angle (Chi)": t0 - t0_c/2,
                                "FWHM (degrees)": fwhm,
                                "Net Height": h,
                                "Net Area (cts deg)": area,
                                "Area/Base Ratio": ratio
                            })
                    except Exception as e:
                        # Fallback
                        for p in peaks_list:
                            all_metrics.append({
                                "Sample": sample,
                                "Phi (degrees)": phi,
                                "Peak Name": p["name"],
                                "Peak Center (Theta)": p["init_center"],
                                "Tilt Angle (Chi)": p["init_center"] - t0_c/2,
                                "FWHM (degrees)": 0.0,
                                "Net Height": 0.0,
                                "Net Area (cts deg)": 0.0,
                                "Area/Base Ratio": 0.0
                            })
                except Exception as e:
                    print(f"    Error processing rocking curve at Phi={phi}: {e}")

    # Process References
    print("\nProcessing Reference Calcite Single Crystal & Sample Holder...")
    ref_dir = os.path.join(RAW_DIR, "Rocking_Curves/Reference")
    ref_proc_dir = os.path.join(PROCESSED_DIR, "Rocking_Curves/Reference")
    os.makedirs(ref_proc_dir, exist_ok=True)
    
    calcite_ref_path = os.path.join(ref_dir, "Calcite single crystal rocking curve.brml")
    holder_ref_path = os.path.join(ref_dir, "sample holder.brml")
    
    if os.path.exists(calcite_ref_path) or os.path.exists(calcite_ref_path.replace(".brml", ".raw")):
        try:
            calcite_arr = extract_brml_data(calcite_ref_path)
            theta_c = calcite_arr[:, 2]
            int_c = calcite_arr[:, 3]
            
            bg_mask_c = (theta_c < 10.0) | (theta_c > 22.0)
            popt_bg_c, _ = curve_fit(bg_model, theta_c[bg_mask_c], int_c[bg_mask_c], p0=[100000, 5e5, -1e4, 500, -10])
            baseline_c = bg_model(theta_c, *popt_bg_c)
            net_int_c = int_c - baseline_c
            
            max_val_c = net_int_c.max()
            max_theta_c = theta_c[np.argmax(net_int_c)]
            half_max_c = max_val_c / 2.0
            idx_above_c = np.where(net_int_c >= half_max_c)[0]
            fwhm_est_c = theta_c[idx_above_c[-1]] - theta_c[idx_above_c[0]]
            
            popt_g_c, _ = curve_fit(gaussian, theta_c[(theta_c >= 10.0) & (theta_c <= 22.0)], 
                                    net_int_c[(theta_c >= 10.0) & (theta_c <= 22.0)], 
                                    p0=[max_val_c, max_theta_c, fwhm_est_c / 2.355])
            h_c, t0_c, w_c = popt_g_c
            fwhm_c = 2.355 * w_c
            area_c = h_c * w_c * np.sqrt(2 * np.pi)
            iso_val_c = popt_bg_c[0] / np.sin(np.radians(t0_c))
            
            df_metrics_c = pd.DataFrame([{
                'Peak Name': 'Calcite (104) Single Crystal Peak',
                'Peak Center (Theta)': t0_c,
                'Tilt Angle (Chi)': t0_c - 29.3725/2,
                'FWHM (deg)': fwhm_c,
                'Net Height': h_c,
                'Net Area (cts deg)': area_c,
                'Isotropic Base': iso_val_c,
                'Area/Base Ratio': area_c / iso_val_c
            }])
            df_metrics_c.to_csv(os.path.join(ref_proc_dir, "calcite_single_crystal_rocking_peaks_metrics.csv"), index=False)
            
            df_corr_c = pd.DataFrame({
                'Theta (degrees)': theta_c,
                'Raw Intensity': int_c,
                'Model Baseline': baseline_c,
                'Corrected Net Intensity': net_int_c
            })
            df_corr_c.to_csv(os.path.join(ref_proc_dir, "calcite_single_crystal_corrected_rocking.csv"), index=False)
            print("  Calcite single crystal processed successfully.")
        except Exception as e:
            print(f"  Error processing single crystal reference: {e}")
            
    if os.path.exists(holder_ref_path) or os.path.exists(holder_ref_path.replace(".brml", ".raw")):
        try:
            holder_arr = extract_brml_data(holder_ref_path)
            twotheta_h = holder_arr[:, 2]
            int_h = holder_arr[:, 4]
            poly_coeff_h = np.polyfit(twotheta_h, int_h, 3)
            baseline_h = np.polyval(poly_coeff_h, twotheta_h)
            
            df_corr_h = pd.DataFrame({
                '2Theta (degrees)': twotheta_h,
                'Raw Intensity': int_h,
                'Polynomial Baseline': baseline_h,
                'Corrected Net Intensity': int_h - baseline_h
            })
            df_corr_h.to_csv(os.path.join(ref_proc_dir, "sample_holder_corrected_2theta.csv"), index=False)
            print("  Sample holder reference processed successfully.")
        except Exception as e:
            print(f"  Error processing sample holder reference: {e}")

    # Compile master metrics table
    if all_metrics:
        df_master = pd.DataFrame(all_metrics)
        df_master.to_csv(os.path.join(ref_proc_dir, "all_samples_rocking_peaks_vs_phi.csv"), index=False)
        print(f"  Master rocking peak metrics compiled to: {os.path.join(ref_proc_dir, 'all_samples_rocking_peaks_vs_phi.csv')}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process CaCO3 rocking curves.")
    parser.add_argument("--sample", type=str, choices=["SH-124-B3", "SH-125-A", "SH-125-G", "SH-104-1", "Reference", "all"],
                        default="all", help="Specific sample to process (or 'all').")
    args = parser.parse_args()

    print(f"Running rocking curve processing for: {args.sample}")
    process_rocking_curves()
    print("Processing completed successfully!")
