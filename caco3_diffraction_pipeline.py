# -*- coding: utf-8 -*-
"""
Comprehensive Diffraction Processing Pipeline
=============================================
Author: Tomasz Stawski, tomasz.stawski@bam.de
Date: 2026-06-19
Version: 0.1

Description:
Implements the core processing pipeline including detector calibration, baseline correction, and peak deconvolution.
"""
"""
Calcium Carbonate (CaCO3) Thin Film Diffraction Processing Pipeline
===================================================================

This script implements the complete data analysis and visualization pipeline for
characterising thin films of calcium carbonate. It processes stationary 2D-XRD data
and azimuthal rocking curves to evaluate phase stability, polymorphic distribution,
and crystallographic orientation (texture).

This script is fully self-contained and reproducible. It processes raw data files
(.gfrm, .brml, and .raw binary blocks) to generate intermediate processed tables
and publication-quality figures.

Processing Steps Included:
--------------------------
1. 2D-XRD Cake Plotting and 1D Profile Integration
2. Symmetric 2Theta-Theta Stacked Profiles
3. Rocking Curve Stacked Net Intensity Sweeps
4. 2D Polar Projection Pole Figures
5. Phase Metrics (calcite vs vaterite peak areas) vs. Azimuthal Angle
6. Peak Deconvolution (Gaussian fitting) and Background Subtraction (Volume Correction)

Usage:
------
Run this script from the package directory:
    python caco3_diffraction_pipeline.py
"""

import os
import glob
import zipfile
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
import h5py
import scipy.ndimage as ndimage
from scipy.interpolate import griddata, interp1d
from scipy.signal import find_peaks, savgol_filter
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# Check for silx convert
try:
    import silx
    from silx.io.convert import convert
    HAS_SILX = True
except ImportError:
    HAS_SILX = False

# Get script folder to resolve relative data paths
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
RAW_DIR = os.path.join(SCRIPT_DIR, "data", "raw")
PROCESSED_DIR = os.path.join(SCRIPT_DIR, "data", "processed")
PLOT_DIR = os.path.join(SCRIPT_DIR, "results", "figures")

os.makedirs(PROCESSED_DIR, exist_ok=True)
os.makedirs(PLOT_DIR, exist_ok=True)

# Set plot style parameters for publication quality
plt.rcParams.update({
    'font.size': 11,
    'font.family': 'sans-serif',
    'font.sans-serif': ['DejaVu Sans', 'Arial', 'Helvetica'],
    'axes.labelsize': 12,
    'axes.titlesize': 13,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'figure.titlesize': 14,
    'svg.fonttype': 'none'
})

# ==============================================================================
# DATA EXTRACTION HELPERS
# ==============================================================================
def extract_brml_data(brml_path):
    """
    Extracts scattering angles and diffraction intensity data arrays from a zipped .brml file.
    Falls back to reading from a corresponding binary .raw file if the .brml file is empty or missing.
    """
    if not os.path.exists(brml_path) or os.path.getsize(brml_path) == 0:
        raw_path = brml_path.replace(".brml", ".raw")
        if os.path.exists(raw_path):
            print(f"  Warning: {os.path.basename(brml_path)} is missing/empty. Recovering data from binary {os.path.basename(raw_path)}...")
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

def load_reference_peaks(filepath, min_intensity=2.0):
    """Parses crystallographic reference peak positions and intensities from a text file."""
    ref_peaks = []
    if os.path.exists(filepath):
        with open(filepath, 'r') as f:
            for line in f:
                if line.strip().startswith('h') or not line.strip():
                    continue
                parts = line.split()
                if len(parts) >= 9:
                    h, k, l = parts[0], parts[1], parts[2]
                    d_spacing = float(parts[3])
                    two_theta = float(parts[7])
                    intensity = float(parts[8])
                    if intensity >= min_intensity:
                        ref_peaks.append((two_theta, intensity, h, k, l, d_spacing))
    return ref_peaks

# ==============================================================================
# FITTING HELPERS
# ==============================================================================
def gaussian(t, h, t0, w):
    """Evaluates a Gaussian peak profile."""
    return h * np.exp(-(t - t0)**2 / (2 * w**2))

def bg_model(t, I0, c0, c1, c2, c3):
    """Defines the background model comprising an isotropic volume correction term and a 3rd-order polynomial."""
    return I0 / np.sin(np.radians(t)) + c0 + c1*t + c2*t**2 + c3*t**3

def fit_symmetric_scan(twotheta, intensity):
    """Fits background baseline and symmetric Bragg reflections for calcite (104) and vaterite (110)."""
    bg_mask = (twotheta < 28.0) | ((twotheta > 31.0) & (twotheta < 31.8)) | (twotheta > 34.0)
    if not bg_mask.any():
        bg_mask = np.ones_like(twotheta, dtype=bool)
    poly_coeff = np.polyfit(twotheta[bg_mask], intensity[bg_mask], 3)
    baseline = np.polyval(poly_coeff, twotheta)
    net_intensity = intensity - baseline

    c_mask = (twotheta >= 28.2) & (twotheta <= 30.8)
    h_c, t0_c, w_c, area_c = 0.0, 29.4, 0.15, 0.0
    try:
        popt_c, _ = curve_fit(gaussian, twotheta[c_mask], net_intensity[c_mask], p0=[intensity.max() - baseline.max(), 29.4, 0.15])
        h_c, t0_c, w_c = popt_c
        area_c = h_c * w_c * np.sqrt(2 * np.pi)
    except:
        pass
        
    # Vaterite (110) at ~32.8
    v_mask = (twotheta >= 31.8) & (twotheta <= 33.8)
    h_v, t0_v, w_v, area_v = 0.0, 32.8, 0.15, 0.0
    try:
        popt_v, _ = curve_fit(gaussian, twotheta[v_mask], net_intensity[v_mask], p0=[1000.0, 32.8, 0.15])
        h_v_fit, t0_v_fit, w_v_fit = popt_v
        fwhm_v_fit = 2.355 * w_v_fit
        if h_v_fit > 0.01 * h_c and 32.4 < t0_v_fit < 33.4 and fwhm_v_fit < 1.0:
            h_v, t0_v, w_v = popt_v
            area_v = h_v * w_v * np.sqrt(2 * np.pi)
    except:
        pass
        
    return {
        "calcite_area": area_c,
        "calcite_center": t0_c,
        "vaterite_area": area_v,
        "vaterite_center": t0_v,
        "baseline": baseline
    }

# ==============================================================================
# PIPELINE STAGE 1: RAW DATA PROCESSING
# ==============================================================================
def run_data_processing():
    """Coordinates the extraction, volume correction, and peak fitting for all 2D-XRD and rocking curve raw data."""
    
    # 1. 2D-XRD GFRM processing
    print("\nProcessing 2D-XRD Detector Frames...")
    gfrm_files = {
        "SH-125-G": "Universtitaet Erlangen Nuernberg Institut_1345_250124_091118-000.gfrm",
        "SH-124-B3": "Universtitaet Erlangen Nuernberg Institut_1343_250123_162242-000.gfrm",
        "SH-125-A": "Universtitaet Erlangen Nuernberg Institut_1339_250115_125507-000.gfrm",
        "SH-104-1": "Universtitaet Erlangen Nuernberg Institut_1347_250128_084659-000.gfrm"
    }
    
    for sample_name, gfrm_filename in gfrm_files.items():
        print(f"\nProcessing 2D-XRD for {sample_name}...")
        gfrm_file = os.path.join(RAW_DIR, f"2D_XRD/{gfrm_filename}")
        h5_file = gfrm_file + ".h5"
        
        if not os.path.exists(h5_file):
            if HAS_SILX:
                print(f"  Converting {os.path.basename(gfrm_file)} to HDF5 using silx...")
                convert(gfrm_file, h5_file)
            else:
                raise FileNotFoundError(f"silx library not installed, and pre-converted HDF5 file {h5_file} not found!")
        else:
            print(f"  Using pre-converted HDF5 frame file for {sample_name}.")
            
        # Read detector image
        with h5py.File(h5_file, 'r') as file:
            start_angle = float(file['scan_0/instrument/detector_0/others/ANGLES'][()][0].decode().strip().split()[0])
            increment = float(file['scan_0/instrument/detector_0/others/INCREME'][()][0])
            ncols = int(file['scan_0/instrument/detector_0/others/NCOLS'][()][0].decode().strip().split()[0])
            ending_angle = float(file['scan_0/instrument/detector_0/others/ENDING'][()][0].decode().strip().split()[0])
            detector_image = file['scan_0/instrument/detector_0/data'][()]
    
        calculated_end = start_angle + (ncols - 1) * increment
        pixel_size = 0.075  # mm
        D = 305.809         # mm
        num_y, num_x = detector_image.shape
        y_mid = (num_y - 1) / 2
    
        two_theta_x = np.linspace(start_angle, calculated_end, ncols)
        y_pixel_pos = (np.arange(num_y) - y_mid) * pixel_size
        two_theta_y = np.degrees(np.arctan2(y_pixel_pos, D))
        two_theta_X_grid, two_theta_Y_grid = np.meshgrid(two_theta_x, two_theta_y)
        two_theta_total = np.sqrt(two_theta_X_grid**2 + two_theta_Y_grid**2)
        phi = np.degrees(np.arctan2(two_theta_Y_grid, two_theta_X_grid))
    
        max_phi = max(abs(np.degrees(np.arctan2((num_y - y_mid - 1) * pixel_size, D))),
                      abs(np.degrees(np.arctan2((0 - y_mid) * pixel_size, D))))
        phi_lin = np.linspace(-max_phi, max_phi, 500)
        theta_lin = two_theta_x
        theta_grid, phi_grid = np.meshgrid(theta_lin, phi_lin)
        
        points = np.vstack((two_theta_total.ravel(), phi.ravel())).T
        values = detector_image.ravel()
        cake_plot = griddata(points, values, (theta_grid, phi_grid), method='linear', fill_value=np.nan)
    
        # 1D profile integration
        integrated_profile = np.nansum(cake_plot, axis=0)
        smoothed_profile = savgol_filter(integrated_profile, window_length=31, polyorder=3)
        
        # morphological baseline
        window_size = 151
        bg_min = ndimage.minimum_filter1d(smoothed_profile, size=window_size)
        bg_max = ndimage.maximum_filter1d(bg_min, size=window_size)
        baseline = savgol_filter(bg_max, window_length=101, polyorder=2)
        corrected_profile = smoothed_profile - baseline
        corrected_profile[corrected_profile < 0] = 0
    
        # Save Excel sheets
        excel_path = os.path.join(PROCESSED_DIR, f"2D_XRD/{sample_name}_output_data.xlsx")
        os.makedirs(os.path.dirname(excel_path), exist_ok=True)
        with pd.ExcelWriter(excel_path, engine='openpyxl') as writer:
            pd.DataFrame({
                'Parameter': ['Start_Angle', 'Increment', 'NCOLS', 'Ending_Angle', 'Calculated_Ending'],
                'Value': [start_angle, increment, ncols, ending_angle, calculated_end]
            }).to_excel(writer, sheet_name='Scan_Metadata', index=False)
            pd.DataFrame(cake_plot, index=phi_lin, columns=theta_lin).to_excel(writer, sheet_name='Cake_Plot_Data')
            pd.DataFrame(detector_image, index=two_theta_y, columns=two_theta_x).to_excel(writer, sheet_name='2D_Original_2ThetaXY')
            pd.DataFrame({
                '2Theta': theta_lin,
                'Smoothed Profile': smoothed_profile,
                'Baseline': baseline,
                'Corrected Profile': corrected_profile
            }).to_excel(writer, sheet_name='Corrected_Profiles')
        print(f"  Processed Excel data written to: {excel_path}")

    # 2. Rocking curve peak fitting configuration
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
                {"name": "Tilted Peak", "init_center": 12.00, "bounds": ([0, 9.0, 0.2], [200000, 21.0, 4.0])},
                {"name": "Specular Peak", "init_center": 15.00, "bounds": ([0, 9.0, 0.2], [200000, 21.0, 4.0])}
            ],
            "default_2theta": 29.3425
        },
        "SH-125-G": {
            "raw_sub_dir": "SH-125-G",
            "bg_mask_fn": lambda theta: (theta >= 4.0) & \
                                       ((theta < 9.0) | ((theta > 18.5) & (theta < 21.5)) | (theta > 23.5)),
            "peaks": [
                {"name": "Tilted Peak", "init_center": 12.00, "bounds": ([0, 9.0, 0.2], [200000, 21.0, 4.0])},
                {"name": "Specular Peak", "init_center": 15.00, "bounds": ([0, 9.0, 0.2], [200000, 21.0, 4.0])}
            ],
            "default_2theta": 29.3425
        },
        "SH-104-1": {
            "raw_sub_dir": "SH-104-1",
            "bg_mask_fn": lambda theta: (theta >= 4.0) & \
                                       ((theta < 9.0) | ((theta > 18.5) & (theta < 21.5)) | (theta > 23.5)),
            "peaks": [
                {"name": "Tilted Peak", "init_center": 12.00, "bounds": ([0, 9.0, 0.2], [200000, 21.0, 4.0])},
                {"name": "Specular Peak", "init_center": 15.00, "bounds": ([0, 9.0, 0.2], [200000, 21.0, 4.0])}
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
        files = os.listdir(sample_raw_dir)
        phi_values = sorted(list(set([
            int(f.split("_")[1].split(".")[0]) for f in files 
            if (f.startswith("2Theta_") or f.startswith("Rocking_") or f.startswith(f"{sample}_2Theta_") or f.startswith(f"{sample}_rocking_"))
            and f.endswith((".brml", ".raw"))
        ])))
        
        # Fallback if names are different
        if not phi_values:
            # Maybe single scan files? (e.g. SH-125-A)
            phi_values = [0]
            
        print(f"  Detected Phi values: {phi_values}")
        
        for phi in phi_values:
            # Look for symmetric scan 2theta file
            twotheta_file = f"2Theta_{phi}.brml"
            # Fallback file names
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
                    
                    # Fit to find center of calcite
                    fit_res = fit_symmetric_scan(twotheta_2t, intensity_2t)
                    if fit_res["calcite_area"] > 0:
                        t0_c = fit_res["calcite_center"]
                except Exception as e:
                    print(f"    Warning: could not process 2Theta scan at Phi={phi}: {e}")

            # Look for rocking curve file
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
                    
                    # Fit background baseline
                    bg_mask = config["bg_mask_fn"](theta)
                    popt_bg = [100000, 5e6, -3e5, 8000, -100]
                    
                    # Iteratively fit baseline to exclude peaks
                    for _ in range(3):
                        try:
                            popt_bg, _ = curve_fit(bg_model, theta[bg_mask], intensity[bg_mask], p0=popt_bg)
                            baseline = bg_model(theta, *popt_bg)
                            residuals = intensity - baseline
                            noise = np.std(residuals[bg_mask])
                            bg_mask = config["bg_mask_fn"](theta) & (residuals < 2.5 * noise)
                        except:
                            # Fallback polynomial
                            poly_coeff = np.polyfit(theta[bg_mask], intensity[bg_mask], 3)
                            popt_bg = [0.0, poly_coeff[3], poly_coeff[2], poly_coeff[1], poly_coeff[0]]
                            break
                            
                    I0_fit, c0_fit, c1_fit, c2_fit, c3_fit = popt_bg
                    baseline = I0_fit / np.sin(np.radians(theta)) + c0_fit + c1_fit*theta + c2_fit*theta**2 + c3_fit*theta**3
                    net_intensity = intensity - baseline
                    
                    # Save corrected curve
                    df_corr = pd.DataFrame({
                        'Theta (degrees)': theta,
                        'Raw Intensity': intensity,
                        'Model Baseline': baseline,
                        'Corrected Net Intensity': net_intensity
                    })
                    df_corr.to_csv(os.path.join(sample_processed_dir, f"{sample}_corrected_rocking_{phi}.csv"), index=False)
                    
                    # Log-scale fitting for peaks
                    residuals_bg = intensity - baseline
                    diffs = np.diff(residuals_bg[config["bg_mask_fn"](theta)])
                    noise_std = np.std(diffs) / np.sqrt(2) if len(diffs) > 1 else 350.0
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
                        w_guess = 1.5 if sample in ["SH-125-A", "SH-125-G", "SH-104-1"] else 0.15
                        flat_guesses.extend([h_guess, center_guess, w_guess])
                        bounds_min.extend(p["bounds"][0])
                        bounds_max.extend(p["bounds"][1])
                        
                    p0 = [I0_fit, c0_fit, c1_fit, c2_fit, c3_fit] + flat_guesses
                    bounds = (bounds_min, bounds_max)
                    log_intensity = np.log10(np.clip(intensity, 1.0, None))
                    
                    try:
                        if sample in ["SH-125-A", "SH-125-G", "SH-104-1"]:
                            # Limit to [7.5, 22.5] and fix baseline parameters for isotropic samples
                            fit_mask = (theta >= 7.5) & (theta <= 22.5)
                            
                            def log_fit_func_fixed(t, *peak_params):
                                baseline_val = I0_fit / np.sin(np.radians(t)) + c0_fit + c1_fit*t + c2_fit*t**2 + c3_fit*t**3
                                val = baseline_val.copy()
                                num_peaks = len(peak_params) // 3
                                for i in range(num_peaks):
                                    h = peak_params[3*i]
                                    t0 = peak_params[3*i+1]
                                    w = peak_params[3*i+2]
                                    val += h * np.exp(-(t-t0)**2 / (2*w**2))
                                return np.log10(np.clip(val, 1.0, None))
                                
                            p0_peaks = flat_guesses
                            bounds_peaks = (bounds_min[5:], bounds_max[5:])
                            
                            popt_peaks, _ = curve_fit(log_fit_func_fixed, theta[fit_mask], log_intensity[fit_mask], p0=p0_peaks, bounds=bounds_peaks)
                            popt = np.concatenate([[I0_fit, c0_fit, c1_fit, c2_fit, c3_fit], popt_peaks])
                        else:
                            # Full range, floating baseline joint fit for the textured sample (SH-124-B3)
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
                                
                            popt, _ = curve_fit(log_fit_func, theta, log_intensity, p0=p0, bounds=bounds)
                            # Update refined baseline in the CSV file for SH-124-B3
                            I0_res, c0_res, c1_res, c2_res, c3_res = popt[0:5]
                            baseline_refined = I0_res / np.sin(np.radians(theta)) + c0_res + c1_res*theta + c2_res*theta**2 + c3_res*theta**3
                            net_intensity_refined = intensity - baseline_refined
                            df_corr["Model Baseline"] = baseline_refined
                            df_corr["Corrected Net Intensity"] = net_intensity_refined
                            df_corr.to_csv(os.path.join(sample_processed_dir, f"{sample}_corrected_rocking_{phi}.csv"), index=False)
                            
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
                        # Fallback simple metrics - try fitting Gaussian directly to net_intensity
                        for p in peaks_list:
                            h_fit, t0_fit, fwhm_fit, area_fit = 0.0, p["init_center"], 0.0, 0.0
                            try:
                                c_min, c_max = p["bounds"][0][1], p["bounds"][1][1]
                                p_mask = (theta >= c_min) & (theta <= c_max)
                                popt_net, _ = curve_fit(
                                    gaussian, 
                                    theta[p_mask], 
                                    net_intensity[p_mask], 
                                    p0=[max(net_intensity[p_mask]), p["init_center"], 2.0]
                                )
                                h_fit, t0_fit, w_fit = popt_net
                                fwhm_fit = 2.355 * w_fit
                                area_fit = h_fit * w_fit * np.sqrt(2 * np.pi)
                            except:
                                pass
                            
                            iso_val = I0_fit / np.sin(np.radians(t0_fit)) if I0_fit > 0 else 1.0
                            ratio = area_fit / iso_val if area_fit > 0 else 0.0
                            
                            all_metrics.append({
                                "Sample": sample,
                                "Phi (degrees)": phi,
                                "Peak Name": p["name"],
                                "Peak Center (Theta)": t0_fit,
                                "Tilt Angle (Chi)": t0_fit - t0_c/2,
                                "FWHM (degrees)": fwhm_fit if h_fit > 10.0 else 0.0,
                                "Net Height": h_fit if h_fit > 10.0 else 0.0,
                                "Net Area (cts deg)": area_fit if h_fit > 10.0 else 0.0,
                                "Area/Base Ratio": ratio if h_fit > 10.0 else 0.0
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
    
    # Single crystal processing
    if os.path.exists(calcite_ref_path) or os.path.exists(calcite_ref_path.replace(".brml", ".raw")):
        try:
            calcite_arr = extract_brml_data(calcite_ref_path)
            theta_c = calcite_arr[:, 2]
            int_c = calcite_arr[:, 3]
            
            # background
            bg_mask_c = (theta_c < 10.0) | (theta_c > 22.0)
            popt_bg_c, _ = curve_fit(bg_model, theta_c[bg_mask_c], int_c[bg_mask_c], p0=[100000, 5e5, -1e4, 500, -10])
            baseline_c = bg_model(theta_c, *popt_bg_c)
            net_int_c = int_c - baseline_c
            
            # Gaussian
            max_val_c = net_int_c.max()
            max_theta_c = theta_c[np.argmax(net_int_c)]
            half_max_c = max_val_c / 2.0
            idx_above_c = np.where(net_int_c >= half_max_c)[0]
            fwhm_est_c = theta_c[idx_above_c[-1]] - theta_c[idx_above_c[0]]
            
            popt_g_c, _ = curve_fit(gaussian, theta_c[(theta_c >= 10.0) & (theta_c <= 22.0)], net_int_c[(theta_c >= 10.0) & (theta_c <= 22.0)], p0=[max_val_c, max_theta_c, fwhm_est_c / 2.355])
            h_c, t0_c, w_c = popt_g_c
            fwhm_c = 2.355 * w_c
            area_c = h_c * w_c * np.sqrt(2 * np.pi)
            iso_val_c = popt_bg_c[0] / np.sin(np.radians(t0_c))
            
            # Save single crystal metrics
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
            
            # Save corrected profiles
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
            
    # Sample holder processing
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

    # Save master metrics table
    df_master = pd.DataFrame(all_metrics)
    df_master.to_csv(os.path.join(ref_proc_dir, "all_samples_rocking_peaks_vs_phi.csv"), index=False)
    print(f"  Master rocking peak metrics compiled to: {os.path.join(ref_proc_dir, 'all_samples_rocking_peaks_vs_phi.csv')}")
    print("\nSTAGE 1 COMPLETE: All raw measurements processed and saved.")

# ==============================================================================
# PIPELINE STAGE 2: FIGURE GENERATION
# ==============================================================================
def generate_all_plots():
    """Generates all publication-grade figures from the processed datasets and saves them as PNG and SVG files under results/figures/."""
    print("\n======================================================================")
    print("STAGE 2: GENERATING ALL PUBLICATION FIGURES")
    print("======================================================================")
    
    # --------------------------------------------------------------------------
    # FIGURE 1: 2D-XRD DETECTOR FRAMES, CAKE PLOTS, AND INTEGRATED PROFILES
    # --------------------------------------------------------------------------
    print("Generating Figure 1 (2D-XRD analysis for all samples)...")
    samples_2d = ["SH-125-G", "SH-124-B3", "SH-125-A", "SH-104-1"]
    for sname in samples_2d:
        excel_path = os.path.join(PROCESSED_DIR, f"2D_XRD/{sname}_output_data.xlsx")
        if os.path.exists(excel_path):
            try:
                print(f"  Plotting 2D-XRD analysis for {sname}...")
                df_orig = pd.read_excel(excel_path, sheet_name='2D_Original_2ThetaXY', index_col=0)
                df_cake = pd.read_excel(excel_path, sheet_name='Cake_Plot_Data', index_col=0)
                df_profile = pd.read_excel(excel_path, sheet_name='Corrected_Profiles')
                
                orig_x = df_orig.columns.values.astype(float)
                orig_y = df_orig.index.values.astype(float)
                orig_image = df_orig.values
                
                theta_lin = df_cake.columns.values.astype(float)
                phi_lin = df_cake.index.values.astype(float)
                cake_plot = df_cake.values
                
                fig = plt.figure(figsize=(10, 12))
                gs = gridspec.GridSpec(3, 1, height_ratios=[1.2, 1.2, 1.0], hspace=0.35)
                
                # Panel (a): Original detector frame
                ax1 = fig.add_subplot(gs[0])
                vmin = np.percentile(orig_image[~np.isnan(orig_image)], 1)
                vmax = np.percentile(orig_image[~np.isnan(orig_image)], 99)
                im1 = ax1.imshow(orig_image, extent=[orig_x.min(), orig_x.max(), orig_y.min(), orig_y.max()],
                                 aspect='auto', origin='lower', cmap='inferno', vmin=vmin, vmax=vmax)
                ax1.set_xlabel('2$\\theta_X$ (°)')
                ax1.set_ylabel('2$\\theta_Y$ (°)')
                ax1.set_title(f'Original 2D detector frame ({sname})', fontweight='bold', fontsize=11)
                fig.colorbar(im1, ax=ax1, label='Intensity (counts)', pad=0.02)
                ax1.text(-0.08, 1.05, "(a)", transform=ax1.transAxes, fontsize=14, fontweight='bold', va='top')
                
                # Panel (b): Resampled polar cake plot
                ax2 = fig.add_subplot(gs[1])
                vmin_c = np.percentile(cake_plot[~np.isnan(cake_plot)], 1)
                vmax_c = np.percentile(cake_plot[~np.isnan(cake_plot)], 99)
                im2 = ax2.imshow(cake_plot, extent=[theta_lin.min(), theta_lin.max(), phi_lin.min(), phi_lin.max()],
                                 aspect='auto', origin='lower', cmap='inferno', vmin=vmin_c, vmax=vmax_c)
                ax2.set_xlabel('2$\\theta$ (°)')
                ax2.set_ylabel('Azimuthal Angle $\phi$ (°)')
                ax2.set_title('Resampled 2D polar cake plot (2$\\theta$ vs. $\phi$)', fontweight='bold', fontsize=11)
                fig.colorbar(im2, ax=ax2, label='Intensity (counts)', pad=0.02)
                ax2.text(-0.08, 1.05, "(b)", transform=ax2.transAxes, fontsize=14, fontweight='bold', va='top')
                
                # Panel (c): 1D integrated profile
                ax3 = fig.add_subplot(gs[2])
                ax3.plot(df_profile['2Theta'], df_profile['Smoothed Profile'], color='#3182bd', alpha=0.5, label='Integrated profile')
                ax3.plot(df_profile['2Theta'], df_profile['Corrected Profile'], color='#e6550d', linewidth=1.5, label='Baseline-corrected')
                ax3.plot(df_profile['2Theta'], df_profile['Baseline'], color='grey', linestyle=':', label='Morphological baseline')
                ax3.set_xlabel('2$\\theta$ (°)')
                ax3.set_ylabel('Intensity (counts)')
                ax3.set_title('Azimuthally integrated 1D profile and phase identification', fontweight='bold', fontsize=11)
                ax3.grid(True, linestyle=':', alpha=0.5)
                
                y_max = df_profile['Smoothed Profile'].max()
                ax3.axvline(29.4, color='#2ca02c', linestyle='--', alpha=0.7, label='calcite (104) ref')
                ax3.text(29.4, y_max * 0.9, ' calcite (104)\n 3.03 Å', color='#2ca02c', ha='left', va='top', fontsize=9)
                
                ax3.axvline(32.8, color='#9467bd', linestyle='--', alpha=0.7, label='vaterite (110) ref')
                ax3.text(32.8, y_max * 0.9, ' vaterite (110)\n 2.73 Å', color='#9467bd', ha='left', va='top', fontsize=9)
                
                ax3.legend(loc='upper right', framealpha=0.9)
                ax3.set_xlim(20, 55)
                ax3.text(-0.08, 1.05, "(c)", transform=ax3.transAxes, fontsize=14, fontweight='bold', va='top')
                
                fig_name = f"fig1_2d_xrd_analysis_{sname.lower().replace('-', '_')}"
                plt.savefig(os.path.join(PLOT_DIR, f"{fig_name}.png"), dpi=300, bbox_inches='tight')
                plt.savefig(os.path.join(PLOT_DIR, f"{fig_name}.svg"), dpi=300, bbox_inches='tight')
                
                # Keep compatibility with original file name for G
                if sname == "SH-125-G":
                    plt.savefig(os.path.join(PLOT_DIR, "fig1_2d_xrd_analysis.png"), dpi=300, bbox_inches='tight')
                    plt.savefig(os.path.join(PLOT_DIR, "fig1_2d_xrd_analysis.svg"), dpi=300, bbox_inches='tight')
                    
                plt.close()
            except Exception as e:
                print(f"  Error loading/plotting 2D-XRD for {sname}: {e}")
            
    # --------------------------------------------------------------------------
    # FIGURE 2: STACKED 2THETA SCANS FOR ALL MEASURED SAMPLES (2x2 GRID)
    # --------------------------------------------------------------------------
    print("Generating Figure 2: Stacked 2Theta scans for all samples...")
    samples_2theta = {
        "SH-125-G": {"dir": os.path.join(PROCESSED_DIR, "Symmetric_Scans/SH-125-G"), "title": "(a) SH-125-G (mixed calcite-vaterite)"},
        "SH-124-B3": {"dir": os.path.join(PROCESSED_DIR, "Symmetric_Scans/SH-124-B3"), "title": "(b) SH-124-B3 (pure calcite)"},
        "SH-125-A": {"dir": os.path.join(PROCESSED_DIR, "Symmetric_Scans/SH-125-A"), "title": "(c) SH-125-A (mixed calcite-vaterite)"},
        "SH-104-1": {"dir": os.path.join(PROCESSED_DIR, "Symmetric_Scans/SH-104-1"), "title": "(d) SH-104-1 (reference, mainly calcite)"}
    }
    
    try:
        fig, axes = plt.subplots(2, 2, figsize=(14, 13), sharex=True)
        axes_flat = axes.flatten()
        
        for idx_s, (sname, config) in enumerate(samples_2theta.items()):
            ax = axes_flat[idx_s]
            s_dir = config["dir"]
            if not os.path.exists(s_dir):
                continue
            xy_files = sorted(glob.glob(os.path.join(s_dir, "*_2Theta_*_exported.xy")))
            if not xy_files:
                continue
                
            offset = 0.0
            phi_to_file = {}
            for f in xy_files:
                phi_val = int(os.path.basename(f).split("2Theta_")[1].split("_")[0])
                phi_to_file[phi_val] = f
            sorted_phis = sorted(phi_to_file.keys())
            colors = plt.cm.viridis(np.linspace(0, 0.85, len(sorted_phis)))
            
            for idx, phi in enumerate(sorted_phis):
                arr = np.loadtxt(phi_to_file[phi], skiprows=1)
                twotheta = arr[:, 0]
                intensity_k = arr[:, 1] / 1e3
                ax.plot(twotheta, intensity_k + offset, color=colors[idx], linewidth=1.5, label=f"$\phi$ = {phi}°")
                ax.axhline(y=offset, color='grey', linestyle='--', linewidth=0.5, alpha=0.5)
                ax.text(twotheta.max() + 0.1, offset + 10.0, f"{phi}°", fontweight='bold', va='center', color=colors[idx])
                offset += 60.0
                
            ax.set_xlabel('2$\\theta$ (°)')
            ax.set_ylabel('Intensity (kcounts, stacked)')
            ax.set_title(config["title"], fontweight='bold', fontsize=12)
            ax.set_xlim(27.0, 35.0)
            ax.set_ylim(-10.0, offset + 150.0)
            ax.axvspan(32.4, 33.4, color='#9467bd', alpha=0.08)
            ax.axvline(29.4, color='#2ca02c', linestyle=':', alpha=0.5)
            ax.axvline(32.8, color='#9467bd', linestyle=':', alpha=0.5)
            ax.text(29.4, offset + 30.0, 'calcite (104)', color='#2ca02c', ha='center', va='bottom', fontsize=9, fontweight='bold', rotation=90)
            ax.text(32.8, offset + 30.0, 'vaterite (110)', color='#9467bd', ha='center', va='bottom', fontsize=9, fontweight='bold', rotation=90)
            ax.grid(True, which='both', axis='x', linestyle=':', alpha=0.5)
            ax.legend(loc='lower left', framealpha=0.5, fontsize=9)
            
        plt.suptitle("Azimuthal dependence of symmetric 2$\\theta-\\theta$ scans for all samples", fontsize=15, fontweight='bold', y=0.98)
        plt.tight_layout()
        plt.savefig(os.path.join(PLOT_DIR, "fig2_stacked_2theta_all_samples.png"), dpi=300, bbox_inches='tight')
        plt.savefig(os.path.join(PLOT_DIR, "fig2_stacked_2theta_all_samples.svg"), dpi=300, bbox_inches='tight')
        plt.close()
        print("  Saved Figure 2 to results/figures/")
    except Exception as e:
        print(f"  Error plotting Figure 2: {e}")

    # --------------------------------------------------------------------------
    # FIGURE 3: STACKED ROCKING CURVES FOR ALL MEASURED CaCO3 SAMPLES (2x2 GRID)
    # --------------------------------------------------------------------------
    print("Generating Figure 3: Stacked Rocking Curves (2x2 Grid, all samples)...")
    samples_config = {
        "SH-124-B3": {
            "processed_dir": os.path.join(PROCESSED_DIR, "Rocking_Curves/SH-124-B3"),
            "phi_values": [0, 30, 60, 90, 120, 150, 180],
            "net_offset": 5000,
            "title": "Sample SH-124-B3 (pure calcite)"
        },
        "SH-125-A": {
            "processed_dir": os.path.join(PROCESSED_DIR, "Rocking_Curves/SH-125-A"),
            "phi_values": [0, 30, 60, 90, 120, 150],
            "net_offset": 5000,
            "title": "Sample SH-125-A (pure calcite)"
        },
        "SH-104-1": {
            "processed_dir": os.path.join(PROCESSED_DIR, "Rocking_Curves/SH-104-1"),
            "phi_values": [0, 30, 60, 90, 120, 150],
            "net_offset": 3000,
            "title": "Sample SH-104-1 (mixed calcite-vaterite)"
        },
        "SH-125-G": {
            "processed_dir": os.path.join(PROCESSED_DIR, "Rocking_Curves/SH-125-G"),
            "phi_values": [0, 30, 60, 120, 150, 180],
            "net_offset": 5000,
            "title": "Sample SH-125-G (mixed calcite-vaterite)"
        }
    }
    try:
        fig, axes = plt.subplots(2, 2, figsize=(14, 12), sharex=True)
        axes_flat = axes.flatten()
        for idx_s, (sample_name, config) in enumerate(samples_config.items()):
            ax = axes_flat[idx_s]
            phi_values = config["phi_values"]
            step_offset = config["net_offset"]
            colors = plt.cm.plasma(np.linspace(0.1, 0.9, len(phi_values)))
            
            for idx_phi, phi in enumerate(phi_values):
                csv_path = os.path.join(config["processed_dir"], f"{sample_name}_corrected_rocking_{phi}.csv")
                if not os.path.exists(csv_path):
                    continue
                df = pd.read_csv(csv_path)
                theta = df["Theta (degrees)"]
                net_intensity = df["Corrected Net Intensity"]
                y_val = net_intensity + idx_phi * step_offset
                
                ax.plot(theta, y_val, '-', linewidth=1.5, color=colors[idx_phi], label=f"$\phi$ = {phi}°")
                ax.axhline(y=idx_phi * step_offset, color='grey', linestyle='--', linewidth=0.8, alpha=0.5)
                ax.text(theta.max() + 0.3, idx_phi * step_offset, f"{phi}°", 
                         verticalalignment='center', fontsize=9, fontweight='bold', color=colors[idx_phi])
                         
            ax.set_xlabel("Theta $\\theta$ (°)")
            ax.set_ylabel("Net Intensity (counts, stacked)")
            ax.set_title(config["title"], fontweight='bold', fontsize=12)
            ax.grid(True, which='both', linestyle=':', alpha=0.5)
            ax.set_xlim(4.0, 26.0)
            ax.text(-0.08, 1.05, f"({chr(97 + idx_s)})", transform=ax.transAxes, fontsize=14, fontweight='bold', va='top')
            
        plt.suptitle("Baseline-corrected net rocking curves vs. azimuthal angle $\phi$", fontsize=14, fontweight='bold', y=0.98)
        plt.tight_layout()
        plt.savefig(os.path.join(PLOT_DIR, "fig3_stacked_net_rocking_curves.png"), dpi=300, bbox_inches='tight')
        plt.savefig(os.path.join(PLOT_DIR, "fig3_stacked_net_rocking_curves.svg"), dpi=300, bbox_inches='tight')
        plt.close()
        print("  Saved Figure 3 to results/figures/")
    except Exception as e:
        print(f"  Error plotting Figure 3: {e}")

    # --------------------------------------------------------------------------
    # FIGURE 4: 2D TEXTURE POLE FIGURES
    # --------------------------------------------------------------------------
    print("Generating Figure 4: 2D Texture Pole Figures (2x2 Grid)...")
    samples_pole = [
        ("SH-124-B3", os.path.join(PROCESSED_DIR, "Rocking_Curves/SH-124-B3"), [0, 30, 60, 90, 120, 150, 180], os.path.join(PROCESSED_DIR, "Symmetric_Scans/SH-125-G")),
        ("SH-125-A", os.path.join(PROCESSED_DIR, "Rocking_Curves/SH-125-A"), [0, 30, 60, 90, 120, 150], os.path.join(PROCESSED_DIR, "Symmetric_Scans/SH-125-A")),
        ("SH-104-1", os.path.join(PROCESSED_DIR, "Rocking_Curves/SH-104-1"), [0, 30, 60, 90, 120, 150], os.path.join(PROCESSED_DIR, "Symmetric_Scans/SH-104-1")),
        ("SH-125-G", os.path.join(PROCESSED_DIR, "Rocking_Curves/SH-125-G"), [0, 30, 60, 120, 150, 180], os.path.join(PROCESSED_DIR, "Symmetric_Scans/SH-125-G"))
    ]
    try:
        fig = plt.figure(figsize=(12, 11))
        gs = gridspec.GridSpec(2, 2, wspace=0.3, hspace=0.3)
        
        for idx_s, (sample, pdir, phis, sym_dir) in enumerate(samples_pole):
            rocking_scans = {}
            theta_0_dict = {}
            
            for phi in phis:
                xy_2t_path = glob.glob(os.path.join(sym_dir, f"*_2Theta_{phi}_exported.xy"))
                theta_0 = 14.7
                if xy_2t_path:
                    try:
                        arr = np.loadtxt(xy_2t_path[0], skiprows=1)
                        twotheta_arr = arr[:, 0]
                        intensity_arr = arr[:, 1]
                        poly_coeff = np.polyfit(twotheta_arr, intensity_arr, 3)
                        net_intensity = intensity_arr - np.polyval(poly_coeff, twotheta_arr)
                        c_mask = (twotheta_arr >= 28.5) & (twotheta_arr <= 30.5)
                        popt_c, _ = curve_fit(gaussian, twotheta_arr[c_mask], net_intensity[c_mask], p0=[intensity_arr.max(), 29.4, 0.15])
                        theta_0 = popt_c[1] / 2.0
                    except:
                        pass
                theta_0_dict[phi] = theta_0
                
                rc_path = os.path.join(pdir, f"{sample}_corrected_rocking_{phi}.csv")
                if os.path.exists(rc_path):
                    df_rc = pd.read_csv(rc_path)
                    rocking_scans[phi] = interp1d(df_rc['Theta (degrees)'].values, df_rc['Corrected Net Intensity'].values, bounds_error=False, fill_value=0.0)
                    
            phi_polar_deg = np.arange(0, 360 + 30, 30)
            phi_polar_rad = np.radians(phi_polar_deg)
            r_tilt = np.linspace(0, 10.0, 300)
            Z_net = np.zeros((len(r_tilt), len(phi_polar_deg)))
            measured_polar_mask = np.ones((len(r_tilt), len(phi_polar_deg)), dtype=bool)
            
            for idx_phi, phi_p in enumerate(phi_polar_deg):
                phi_mapped = phi_p % 360
                net_profile = np.zeros_like(r_tilt)
                measured = True
                
                if phi_mapped <= 180:
                    scan_phi = phi_mapped
                    if scan_phi in rocking_scans:
                        theta_0 = theta_0_dict.get(scan_phi, 14.7)
                        net_profile = rocking_scans[scan_phi](theta_0 + r_tilt)
                    else:
                        measured = False
                else:
                    scan_phi = phi_mapped - 180
                    if scan_phi in rocking_scans:
                        theta_0 = theta_0_dict.get(scan_phi, 14.7)
                        net_profile = rocking_scans[scan_phi](theta_0 - r_tilt)
                    else:
                        measured = False
                        
                Z_net[:, idx_phi] = net_profile
                if not measured:
                    measured_polar_mask[:, idx_phi] = False
                    
            Z_net[0, :] = np.nanmean(Z_net[0, :])
            Z_net = np.clip(Z_net, 0, None)
            
            Z_net_masked = Z_net.copy()
            Z_net_masked[~measured_polar_mask] = np.nan
            
            ax = fig.add_subplot(gs[idx_s], projection='polar')
            Phi_mesh, R_mesh = np.meshgrid(phi_polar_rad, r_tilt)
            
            z_max = np.nanmax(Z_net_masked) if not np.all(np.isnan(Z_net_masked)) else 1000.0
            levels = np.linspace(0, max(z_max, 10.0), 60)
            contour = ax.contourf(Phi_mesh, R_mesh, Z_net_masked, levels=levels, cmap='plasma')
            
            ax.set_theta_zero_location("E")
            ax.set_theta_direction(1)
            ax.set_ylim(0, 10.0)
            ax.set_yticks([2, 4, 6, 8, 10])
            ax.set_yticklabels(["2°", "4°", "6°", "8°", "10°"], color="#555555", size=9)
            ax.set_rlabel_position(45)
            
            ax.set_xticks(np.radians(np.arange(0, 360, 30)))
            ax.set_xticklabels([f"{d}°" for d in np.arange(0, 360, 30)], size=9)
            ax.grid(True, color="grey", linestyle=":", alpha=0.5)
            ax.set_title(f"({chr(97 + idx_s)}) {sample}", y=1.08, fontweight='bold', fontsize=12)
            
            cbar = fig.colorbar(contour, ax=ax, pad=0.08, shrink=0.7)
            cbar.ax.tick_params(labelsize=9)
            cbar.set_label("Net Intensity (counts)", fontsize=9)
            
        plt.suptitle("Compiled calcite (104) 2D polar texture pole figures", fontsize=15, fontweight='bold', y=0.98)
        plt.savefig(os.path.join(PLOT_DIR, "fig4_texture_pole_figures.png"), dpi=300, bbox_inches='tight')
        plt.savefig(os.path.join(PLOT_DIR, "fig4_texture_pole_figures.svg"), dpi=300, bbox_inches='tight')
        plt.close()
        print("  Saved Figure 4 to results/figures/")
    except Exception as e:
        print(f"  Error plotting Figure 4: {e}")

    # --------------------------------------------------------------------------
    # FIGURE 5: PHASE METRICS VS PHI FOR ALL MEASURED SAMPLES
    # --------------------------------------------------------------------------
    print("Generating Figure 5: Phase Areas vs. Phi...")
    samples_phase = {
        "SH-124-B3": {
            "dir": os.path.join(PROCESSED_DIR, "Symmetric_Scans/SH-124-B3"),
            "phis": [0, 30, 60, 90, 120, 150, 180],
            "color_c": "#7f7f7f", "color_v": "#bcbd22", "marker": "d", "ls": ":"
        },
        "SH-125-A": {
            "dir": os.path.join(PROCESSED_DIR, "Symmetric_Scans/SH-125-A"),
            "phis": [0, 30, 60, 90, 120, 150],
            "color_c": "#2ca02c", "color_v": "#9467bd", "marker": "o", "ls": "-"
        },
        "SH-104-1": {
            "dir": os.path.join(PROCESSED_DIR, "Symmetric_Scans/SH-104-1"),
            "phis": [0, 30, 60, 90, 120, 150],
            "color_c": "#1f77b4", "color_v": "#d62728", "marker": "s", "ls": "--"
        },
        "SH-125-G": {
            "dir": os.path.join(PROCESSED_DIR, "Symmetric_Scans/SH-125-G"),
            "phis": [0, 30, 60, 90, 120, 150, 180],
            "color_c": "#ff7f0e", "color_v": "#8c564b", "marker": "^", "ls": "-."
        }
    }
    try:
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 8), sharex=True)
        
        for sample, config in samples_phase.items():
            pdir = config["dir"]
            phis = config["phis"]
            calcite_areas = []
            vaterite_areas = []
            valid_phis = []
            
            for phi in phis:
                xy_files = glob.glob(os.path.join(pdir, f"*_2Theta_{phi}_exported.xy"))
                if not xy_files:
                    continue
                try:
                    arr = np.loadtxt(xy_files[0], skiprows=1)
                    twotheta = arr[:, 0]
                    intensity = arr[:, 1]
                    
                    fit_res = fit_symmetric_scan(twotheta, intensity)
                    calcite_areas.append(fit_res["calcite_area"] / 1e3)
                    vaterite_areas.append(fit_res["vaterite_area"] / 1e3)
                    valid_phis.append(phi)
                except Exception as e:
                    print(f"    Error fitting {sample} Phi {phi}: {e}")
                    
            ax1.plot(valid_phis, calcite_areas, marker=config["marker"], linestyle=config["ls"],
                     color=config["color_c"], linewidth=2, label=f"{sample} calcite (104)")
            ax2.plot(valid_phis, vaterite_areas, marker=config["marker"], linestyle=config["ls"],
                     color=config["color_v"], linewidth=2, label=f"{sample} vaterite (110)")
                     
        ax1.set_ylabel("calcite (104) peak area\n(kcounts·°)")
        ax1.set_title("Peak areas vs. azimuthal angle $\phi$ (symmetric XRD scans)", fontweight='bold')
        ax1.grid(True, linestyle=':', alpha=0.5)
        ax1.legend(loc='upper right')
        ax1.text(-0.08, 1.05, "(a)", transform=ax1.transAxes, fontsize=14, fontweight='bold', va='top')
        
        ax2.set_xlabel("Azimuthal Angle $\phi$ (°)")
        ax2.set_ylabel("vaterite (110) peak area\n(kcounts·°)")
        ax2.grid(True, linestyle=':', alpha=0.5)
        ax2.legend(loc='upper right')
        ax2.set_xticks(np.arange(0, 210, 30))
        ax2.set_xlim(-5, 185)
        ax2.text(-0.08, 1.05, "(b)", transform=ax2.transAxes, fontsize=14, fontweight='bold', va='top')
        ax2.axvspan(25, 65, color='#9467bd', alpha=0.1)
        ax2.text(45, ax2.get_ylim()[1]*0.8, "Epitaxial vaterite\norientation zone", color='#5c3d75', ha='center', fontweight='bold', fontsize=9)
        
        plt.tight_layout()
        plt.savefig(os.path.join(PLOT_DIR, "fig5_phase_metrics_vs_phi.png"), dpi=300, bbox_inches='tight')
        plt.savefig(os.path.join(PLOT_DIR, "fig5_phase_metrics_vs_phi.svg"), dpi=300, bbox_inches='tight')
        plt.close()
        print("  Saved Figure 5 to results/figures/")
    except Exception as e:
        print(f"  Error plotting Figure 5: {e}")

    # --------------------------------------------------------------------------
    # FIGURE 6: SH-125-A ROCKING CURVE ANALYSIS
    # --------------------------------------------------------------------------
    print("Generating Figure 6: Raw vs. volume-corrected rocking curve for SH-125-A...")
    csv_path = os.path.join(PROCESSED_DIR, "Rocking_Curves/SH-125-A/SH-125-A_corrected_rocking_0.csv")
    metrics_path = os.path.join(PROCESSED_DIR, "Rocking_Curves/Reference/all_samples_rocking_peaks_vs_phi.csv")
    if os.path.exists(csv_path) and os.path.exists(metrics_path):
        try:
            df_rc = pd.read_csv(csv_path)
            theta = df_rc['Theta (degrees)'].values
            intensity = df_rc['Raw Intensity'].values
            baseline = df_rc['Model Baseline'].values
            net_intensity = df_rc['Corrected Net Intensity'].values
            
            df_m = pd.read_csv(metrics_path)
            df_scan = df_m[(df_m['Sample'] == 'SH-125-A') & (df_m['Phi (degrees)'] == 0)]
            
            fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
            axes[0].plot(theta, intensity, 'b.', label='Experimental Data', alpha=0.6)
            axes[0].plot(theta, baseline, 'r-', linewidth=2, label='Model Baseline (Isotropic + Background)')
            axes[0].set_ylabel('Intensity (counts)')
            axes[0].legend()
            axes[0].grid(True)
            axes[0].set_title('SH-125-A Rocking Curve Analysis (2Theta = 29.3425°)')
            
            axes[1].plot(theta, net_intensity, '.', color='#7f7f7f', markersize=4, alpha=0.5, label='Experimental Net Data')
            fit_total = np.zeros_like(theta)
            for idx, row in df_scan.iterrows():
                t0 = row['Peak Center (Theta)']
                h = row['Net Height']
                fwhm = row['FWHM (degrees)']
                w = fwhm / 2.355
                p_name = row['Peak Name']
                if h > 0:
                    y_peak = gaussian(theta, h, t0, w)
                    fit_total += y_peak
                    is_tilt = "Tilt" in p_name or "Tilted" in p_name
                    p_color = '#2ca02c' if is_tilt else '#1f77b4'
                    axes[1].plot(theta, y_peak, color=p_color, linestyle=':', linewidth=1.5, label=f'Fitted {p_name}')
                    axes[1].fill_between(theta, 0, y_peak, color=p_color, alpha=0.08)
                    # Stagger annotations to prevent overlapping
                    text_x = t0 - 1.8 if is_tilt else t0 + 0.3
                    text_y = h * 0.6 if is_tilt else h * 0.8
                    axes[1].text(text_x, text_y, f"{p_name}\nCenter={t0:.2f}°\nFWHM={fwhm:.2f}°\n$\\chi$={row['Tilt Angle (Chi)']:.2f}°", 
                                 fontsize=9, color=p_color, fontweight='bold', 
                                 bbox=dict(facecolor='white', alpha=0.85, boxstyle='round,pad=0.2'))
            axes[1].plot(theta, fit_total, 'k-', linewidth=2, label='Total Fit Envelope')
            axes[1].axhline(0, color='gray', linestyle='--')
            axes[1].set_xlabel('Theta (degrees)')
            axes[1].set_ylabel('Net Residual Intensity (counts)')
            axes[1].legend()
            axes[1].grid(True)
            plt.tight_layout()
            plt.savefig(os.path.join(PLOT_DIR, "SH-125-A_rocking_curve_analysis.png"), dpi=150)
            plt.close()
            print("  Saved Figure 6 to results/figures/")
        except Exception as e:
            print(f"  Error plotting Figure 6: {e}")

    # --------------------------------------------------------------------------
    # FIGURE 7: SH-125-G SIDE-BY-SIDE ROCKING CURVES
    # --------------------------------------------------------------------------
    print("Generating Figure 7: Stacked side-by-side curves for SH-125-G...")
    g_rc_dir = os.path.join(PROCESSED_DIR, "Rocking_Curves/SH-125-G")
    if os.path.exists(g_rc_dir):
        try:
            phi_values = [0, 30, 60, 120, 150, 180]
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 7))
            
            factor = 1.0
            for phi in phi_values:
                csv_path = os.path.join(g_rc_dir, f"SH-125-G_corrected_rocking_{phi}.csv")
                if not os.path.exists(csv_path): continue
                df = pd.read_csv(csv_path)
                theta = df["Theta (degrees)"]
                raw_intensity = df["Raw Intensity"]
                baseline = df["Model Baseline"]
                ax1.plot(theta, raw_intensity * factor, '.', alpha=0.4, label=f"{phi}° ($\times${factor:.1e})")
                ax1.plot(theta, baseline * factor, 'r-', linewidth=1.2)
                factor *= 5.0
            ax1.set_yscale('log')
            ax1.set_xlabel("Theta (degrees)")
            ax1.set_ylabel("Intensity (counts, stacked log scale)")
            ax1.set_title("Raw Rocking Curves & Baselines (Log Scale)")
            ax1.grid(True, which='both', linestyle=':', alpha=0.5)
            ax1.legend(loc='upper right', fontsize=8, ncol=2)
            
            step_offset = 5000
            for idx, phi in enumerate(phi_values):
                csv_path = os.path.join(g_rc_dir, f"SH-125-G_corrected_rocking_{phi}.csv")
                if not os.path.exists(csv_path): continue
                df = pd.read_csv(csv_path)
                theta = df["Theta (degrees)"]
                net_intensity = df["Corrected Net Intensity"]
                y_val = net_intensity + idx * step_offset
                ax2.plot(theta, y_val, '-', linewidth=1.5, label=f"Phi = {phi}°")
                ax2.axhline(y=idx * step_offset, color='grey', linestyle='--', linewidth=0.8, alpha=0.7)
                ax2.text(theta.max() + 0.2, idx * step_offset, f"{phi}°", verticalalignment='center', fontsize=9, fontweight='bold')
            ax2.set_xlabel("Theta (degrees)")
            ax2.set_ylabel("Net Intensity (counts, stacked linear scale)")
            ax2.set_title("Baseline-Corrected Net Curves (Linear Scale)")
            ax2.grid(True, which='both', linestyle=':', alpha=0.5)
            ax2.set_xlim(theta.min() - 0.5, theta.max() + 1.5)
            fig.suptitle("SH-125-G Rocking Curve Analysis: Raw vs. Baseline-Corrected", fontsize=14, fontweight='bold')
            plt.tight_layout()
            plt.savefig(os.path.join(PLOT_DIR, "SH-125-G_side_by_side.png"), dpi=150)
            plt.close()
            print("  Saved Figure 7 to results/figures/")
        except Exception as e:
            print(f"  Error plotting Figure 7: {e}")

    # --------------------------------------------------------------------------
    # FIGURE 8: SH-104-1 SIDE-BY-SIDE ROCKING CURVES
    # --------------------------------------------------------------------------
    print("Generating Figure 8: Stacked side-by-side curves for SH-104-1...")
    ref_rc_dir = os.path.join(PROCESSED_DIR, "Rocking_Curves/SH-104-1")
    if os.path.exists(ref_rc_dir):
        try:
            phi_values = [0, 30, 60, 90, 120, 150]
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 7))
            
            factor = 1.0
            for phi in phi_values:
                csv_path = os.path.join(ref_rc_dir, f"SH-104-1_corrected_rocking_{phi}.csv")
                if not os.path.exists(csv_path): continue
                df = pd.read_csv(csv_path)
                theta = df["Theta (degrees)"]
                raw_intensity = df["Raw Intensity"]
                baseline = df["Model Baseline"]
                ax1.plot(theta, raw_intensity * factor, '.', alpha=0.4, label=f"{phi}° ($\times${factor:.1e})")
                ax1.plot(theta, baseline * factor, 'r-', linewidth=1.2)
                factor *= 5.0
            ax1.set_yscale('log')
            ax1.set_xlabel("Theta (degrees)")
            ax1.set_ylabel("Intensity (counts, stacked log scale)")
            ax1.set_title("Raw Rocking Curves & Baselines (Log Scale)")
            ax1.grid(True, which='both', linestyle=':', alpha=0.5)
            ax1.legend(loc='upper right', fontsize=8, ncol=2)
            
            step_offset = 3000
            for idx, phi in enumerate(phi_values):
                csv_path = os.path.join(ref_rc_dir, f"SH-104-1_corrected_rocking_{phi}.csv")
                if not os.path.exists(csv_path): continue
                df = pd.read_csv(csv_path)
                theta = df["Theta (degrees)"]
                net_intensity = df["Corrected Net Intensity"]
                y_val = net_intensity + idx * step_offset
                ax2.plot(theta, y_val, '-', linewidth=1.5, label=f"Phi = {phi}°")
                ax2.axhline(y=idx * step_offset, color='grey', linestyle='--', linewidth=0.8, alpha=0.7)
                ax2.text(theta.max() + 0.2, idx * step_offset, f"{phi}°", verticalalignment='center', fontsize=9, fontweight='bold')
            ax2.set_xlabel("Theta (degrees)")
            ax2.set_ylabel("Net Intensity (counts, stacked linear scale)")
            ax2.set_title("Baseline-Corrected Net Curves (Linear Scale)")
            ax2.grid(True, which='both', linestyle=':', alpha=0.5)
            ax2.set_xlim(theta.min() - 0.5, theta.max() + 1.5)
            fig.suptitle("SH-104-1 Rocking Curve Analysis: Raw vs. Baseline-Corrected", fontsize=14, fontweight='bold')
            plt.tight_layout()
            plt.savefig(os.path.join(PLOT_DIR, "SH-104-1_side_by_side.png"), dpi=150)
            plt.close()
            print("  Saved Figure 8 to results/figures/")
        except Exception as e:
            print(f"  Error plotting Figure 8: {e}")

    # --------------------------------------------------------------------------
    # FIGURE A1: BACKGROUND SUBTRACTION sequence (SH-124-B3, phi=60)
    # --------------------------------------------------------------------------
    print("Generating Appendix Figure A1: Rocking curve background subtraction sequence...")
    rc_60_csv = os.path.join(PROCESSED_DIR, "Rocking_Curves/SH-124-B3/SH-124-B3_corrected_rocking_60.csv")
    metrics_path = os.path.join(PROCESSED_DIR, "Rocking_Curves/Reference/all_samples_rocking_peaks_vs_phi.csv")
    if os.path.exists(rc_60_csv) and os.path.exists(metrics_path):
        try:
            df_rc = pd.read_csv(rc_60_csv)
            theta = df_rc['Theta (degrees)'].values
            raw_int = df_rc['Raw Intensity'].values
            baseline = df_rc['Model Baseline'].values
            net_raw = raw_int - baseline
            
            df_m = pd.read_csv(metrics_path)
            df_scan = df_m[(df_m['Sample'] == 'SH-124-B3') & (df_m['Phi (degrees)'] == 60)]
            
            full_fit = baseline.copy()
            net_fit = np.zeros_like(theta)
            peaks_to_plot = []
            colors_peaks = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
            
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
                    net_fit += peak_y
                    peaks_to_plot.append({
                        "name": name, "center": center, "height": height, "w": w, "tilt": tilt, "curve": peak_y
                    })
            
            fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 5.5))
            fig.suptitle("Rocking curve background subtraction sequence (SH-124-B3, $\phi = 60^\circ$)", fontsize=14, fontweight='bold', y=0.98)
            
            ax1.plot(theta, raw_int, 'o', color='gray', markersize=3, alpha=0.5, label='Raw data')
            ax1.set_xlabel("Theta $\\theta$ (°)")
            ax1.set_ylabel("Intensity (counts)")
            ax1.set_title("Before: Raw Intensity Profile")
            ax1.grid(True, linestyle=':', alpha=0.5)
            ax1.legend(loc='upper right')
            ax1.set_ylim(0, np.max(raw_int)*1.1)
            ax1.text(0.05, 0.95, "Before", transform=ax1.transAxes, fontsize=12, fontweight='bold', va='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
            
            ax2.plot(theta, raw_int, 'o', color='gray', markersize=3, alpha=0.4, label='Raw data')
            ax2.plot(theta, baseline, 'r-', linewidth=1.5, label='Background baseline')
            ax2.plot(theta, full_fit, 'k-', linewidth=2, label='Total fit envelope')
            ax2.set_xlabel("Theta $\\theta$ (°)")
            ax2.set_ylabel("Intensity (counts)")
            ax2.set_title("During: Background & Envelope Fit")
            ax2.grid(True, linestyle=':', alpha=0.5)
            ax2.legend(loc='upper right')
            ax2.set_ylim(0, np.max(raw_int)*1.1)
            ax2.text(0.05, 0.95, "During", transform=ax2.transAxes, fontsize=12, fontweight='bold', va='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
            
            ax3.plot(theta, net_raw, 'o', color='#7f7f7f', markersize=3, alpha=0.5, label='Experimental net data')
            ax3.plot(theta, net_fit, 'k-', linewidth=2, label='Total net fit envelope')
            ax3.axhline(0, color='red', linestyle='-', linewidth=1.5, label='Subtracted baseline (y=0)')
            for i, p in enumerate(peaks_to_plot):
                color = colors_peaks[i % len(colors_peaks)]
                ax3.plot(theta, p["curve"], '--', color=color, linewidth=1.2, label=f"{p['name']} ($\\chi \\approx {p['tilt']:.1f}^\\circ$)")
                ax3.fill_between(theta, 0, p["curve"], color=color, alpha=0.12)
                ax3.axvline(p["center"], color=color, linestyle=':', alpha=0.4)
            ax3.set_xlabel("Theta $\\theta$ (°)")
            ax3.set_ylabel("Net Intensity (counts)")
            ax3.set_title("After: Deconvoluted Net Peaks")
            ax3.grid(True, linestyle=':', alpha=0.5)
            ax3.legend(loc='upper right', fontsize=8, framealpha=0.9)
            ax3.set_ylim(min(-200, np.min(net_raw)*1.1), np.max(net_raw)*1.1)
            ax3.text(0.05, 0.95, "After", transform=ax3.transAxes, fontsize=12, fontweight='bold', va='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
            
            plt.tight_layout(rect=[0, 0, 1, 0.95])
            plt.savefig(os.path.join(PLOT_DIR, "fig_a1_background_subtraction.png"), dpi=300, bbox_inches='tight')
            plt.close()
            print("  Saved Figure A1 to results/figures/")
        except Exception as e:
            print(f"  Error plotting Figure A1: {e}")

    # --------------------------------------------------------------------------
    # FIGURE A2: PEAK DECONVOLUTION zoom (SH-124-B3, phi=60)
    # --------------------------------------------------------------------------
    print("Generating Appendix Figure A2: Zoomed net peak deconvolution...")
    if os.path.exists(rc_60_csv) and os.path.exists(metrics_path):
        try:
            df_rc = pd.read_csv(rc_60_csv)
            theta = df_rc['Theta (degrees)'].values
            raw_int = df_rc['Raw Intensity'].values
            baseline = df_rc['Model Baseline'].values
            net_raw = raw_int - baseline
            
            df_m = pd.read_csv(metrics_path)
            df_scan = df_m[(df_m['Sample'] == 'SH-124-B3') & (df_m['Phi (degrees)'] == 60)]
            
            net_fit = np.zeros_like(theta)
            peaks_to_plot = []
            colors_peaks = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
            
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
                        "name": name, "center": center, "height": height, "w": w, "tilt": tilt, "curve": peak_y
                    })
                    
            zoom_min, zoom_max = 10.0, 14.5
            mask = (theta >= zoom_min) & (theta <= zoom_max)
            
            plt.figure(figsize=(9, 6))
            plt.plot(theta, net_raw, 'o', color='#7f7f7f', markersize=5, alpha=0.7, label='Experimental net data (raw - baseline)')
            plt.axhline(0, color='red', linestyle='-', linewidth=1.5, label='Subtracted background (y=0)')
            plt.plot(theta, net_fit, 'k-', linewidth=2.5, label='Total net fit envelope')
            
            peak_idx = 0
            for p in peaks_to_plot:
                if zoom_min - 1.0 <= p["center"] <= zoom_max + 1.0:
                    color = colors_peaks[peak_idx % len(colors_peaks)]
                    plt.plot(theta, p["curve"], '--', color=color, linewidth=1.8, label=f"{p['name']} ($\\chi \\approx {p['tilt']:.2f}^\\circ$)")
                    plt.fill_between(theta, 0, p["curve"], color=color, alpha=0.15)
                    plt.axvline(p["center"], color=color, linestyle=':', alpha=0.6)
                    peak_idx += 1
            plt.xlim(zoom_min, zoom_max)
            y_data = net_raw[mask]
            y_fit = net_fit[mask]
            y_max = max(np.max(y_data), np.max(y_fit))
            y_min = min(np.min(y_data), 0)
            plt.ylim(y_min - (y_max - y_min)*0.1, y_max + (y_max - y_min)*0.1)
            plt.xlabel("Theta $\\theta$ (°)", fontsize=12)
            plt.ylabel("Net Intensity (counts)", fontsize=12)
            plt.title("Deconvoluted net rocking curve and tilt components (SH-124-B3, $\phi = 60^\circ$)", fontsize=13, fontweight='bold')
            plt.grid(True, which='both', linestyle=':', alpha=0.5)
            plt.legend(loc='upper right', fontsize=10, framealpha=0.9)
            plt.tight_layout()
            plt.savefig(os.path.join(PLOT_DIR, "fig_a2_peak_deconvolution.png"), dpi=300, bbox_inches='tight')
            plt.close()
            print("  Saved Figure A2 to results/figures/")
        except Exception as e:
            print(f"  Error plotting Figure A2: {e}")

    # --------------------------------------------------------------------------
    # FIGURE A3: ZOOMED FITTING FOR SH-124-B3 AT PHI = 30 (LOG SCALE)
    # --------------------------------------------------------------------------
    print("Generating Appendix Figure A3: Zoomed rocking curve fit at phi=30...")
    rc_30_csv = os.path.join(PROCESSED_DIR, "Rocking_Curves/SH-124-B3/SH-124-B3_corrected_rocking_30.csv")
    if os.path.exists(rc_30_csv) and os.path.exists(metrics_path):
        try:
            df_rc = pd.read_csv(rc_30_csv)
            theta = df_rc['Theta (degrees)'].values
            raw_int = df_rc['Raw Intensity'].values
            baseline = df_rc['Model Baseline'].values
            
            df_m = pd.read_csv(metrics_path)
            df_scan = df_m[(df_m['Sample'] == 'SH-124-B3') & (df_m['Phi (degrees)'] == 30)]
            
            full_fit = baseline.copy()
            peaks_to_plot = []
            colors_peaks = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2']
            
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
                        "name": name, "center": center, "height": height, "w": w, "tilt": tilt, "curve": peak_y
                    })
            
            zoom_min, zoom_max = 10.0, 14.5
            mask = (theta >= zoom_min) & (theta <= zoom_max)
            
            plt.figure(figsize=(10, 6))
            plt.plot(theta, raw_int, 'o', color='#7f7f7f', markersize=4, alpha=0.6, label='Raw Data Points')
            plt.plot(theta, baseline, 'r-', linewidth=2, label='Background Baseline')
            plt.plot(theta, full_fit, 'k-', linewidth=2.5, label='Total Fit Envelope')
            
            peak_idx = 0
            for p in peaks_to_plot:
                if zoom_min - 1.0 <= p["center"] <= zoom_max + 1.0:
                    color = colors_peaks[peak_idx % len(colors_peaks)]
                    peak_profile = baseline + p["curve"]
                    plt.plot(theta, peak_profile, '--', color=color, linewidth=1.5, label=f"{p['name']} ($\\chi \\approx {p['tilt']:.2f}^\\circ$)")
                    plt.fill_between(theta, baseline, peak_profile, color=color, alpha=0.15)
                    plt.axvline(p["center"], color=color, linestyle=':', alpha=0.6)
                    peak_idx += 1
                    
            plt.xlim(zoom_min, zoom_max)
            y_data = raw_int[mask]
            y_fit = full_fit[mask]
            y_max = max(np.max(y_data), np.max(y_fit))
            y_min = min(np.min(y_data), np.min(baseline[mask]))
            plt.ylim(max(0, y_min - (y_max - y_min)*0.1), y_max + (y_max - y_min)*0.1)
            plt.xlabel("Theta (degrees)", fontsize=13)
            plt.ylabel("Intensity (counts)", fontsize=13)
            plt.title("SH-124-B3 Fit Zoom — Phi = 30° (Peak 2a/2b Region)", fontsize=14, fontweight='bold')
            plt.grid(True, which='both', linestyle=':', alpha=0.6)
            plt.legend(loc='upper right', fontsize=10, framealpha=0.9)
            plt.tight_layout()
            plt.savefig(os.path.join(PLOT_DIR, "SH-124-B3_fit_phi_30_zoom.png"), dpi=150)
            plt.close()
            print("  Saved Figure A3 to results/figures/")
        except Exception as e:
            print(f"  Error plotting Figure A3: {e}")

    # --------------------------------------------------------------------------
    # FIGURE A4: ZOOMED NET DECONVOLUTION FOR SH-124-B3 AT PHI = 30
    # --------------------------------------------------------------------------
    print("Generating Appendix Figure A4: Zoomed net peak deconvolution at phi=30...")
    if os.path.exists(rc_30_csv) and os.path.exists(metrics_path):
        try:
            df_rc = pd.read_csv(rc_30_csv)
            theta = df_rc['Theta (degrees)'].values
            raw_int = df_rc['Raw Intensity'].values
            baseline = df_rc['Model Baseline'].values
            net_raw = raw_int - baseline
            
            df_m = pd.read_csv(metrics_path)
            df_scan = df_m[(df_m['Sample'] == 'SH-124-B3') & (df_m['Phi (degrees)'] == 30)]
            
            net_fit = np.zeros_like(theta)
            peaks_to_plot = []
            colors_peaks = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2']
            
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
                        "name": name, "center": center, "height": height, "w": w, "tilt": tilt, "curve": peak_y
                    })
                    
            zoom_min, zoom_max = 10.0, 14.5
            mask = (theta >= zoom_min) & (theta <= zoom_max)
            
            plt.figure(figsize=(10, 6))
            plt.plot(theta, net_raw, 'o', color='#7f7f7f', markersize=5, alpha=0.7, label='Experimental Net Data (Raw - Baseline)')
            plt.axhline(0, color='red', linestyle='-', linewidth=1.5, label='Subtracted Background (y=0)')
            plt.plot(theta, net_fit, 'k-', linewidth=2.5, label='Total Net Fit Envelope')
            
            peak_idx = 0
            for p in peaks_to_plot:
                if zoom_min - 1.0 <= p["center"] <= zoom_max + 1.0:
                    color = colors_peaks[peak_idx % len(colors_peaks)]
                    plt.plot(theta, p["curve"], '--', color=color, linewidth=1.8, label=f"{p['name']} ($\\chi \\approx {p['tilt']:.2f}^\\circ$)")
                    plt.fill_between(theta, 0, p["curve"], color=color, alpha=0.15)
                    plt.axvline(p["center"], color=color, linestyle=':', alpha=0.6)
                    peak_idx += 1
                    
            plt.xlim(zoom_min, zoom_max)
            y_data = net_raw[mask]
            y_fit = net_fit[mask]
            y_max = max(np.max(y_data), np.max(y_fit))
            y_min = min(np.min(y_data), 0)
            plt.ylim(y_min - (y_max - y_min)*0.1, y_max + (y_max - y_min)*0.1)
            plt.xlabel("Theta (degrees)", fontsize=13)
            plt.ylabel("Net Intensity (counts)", fontsize=13)
            plt.title("SH-124-B3 Net Fit Zoom — Phi = 30° (Baseline Subtracted)", fontsize=14, fontweight='bold')
            plt.grid(True, which='both', linestyle=':', alpha=0.6)
            plt.legend(loc='upper right', fontsize=10, framealpha=0.9)
            plt.tight_layout()
            plt.savefig(os.path.join(PLOT_DIR, "SH-124-B3_fit_phi_30_net_zoom.png"), dpi=150)
            plt.close()
            print("  Saved Figure A4 to results/figures/")
        except Exception as e:
            print(f"  Error plotting Figure A4: {e}")

    # --------------------------------------------------------------------------
    # FIGURE A5: CALCITE SINGLE CRYSTAL ROCKING CURVE ANALYSIS
    # --------------------------------------------------------------------------
    print("Generating Appendix Figure A5: Calcite single crystal rocking curve...")
    ref_c_csv = os.path.join(PROCESSED_DIR, "Rocking_Curves/Reference/calcite_single_crystal_corrected_rocking.csv")
    ref_c_metrics = os.path.join(PROCESSED_DIR, "Rocking_Curves/Reference/calcite_single_crystal_rocking_peaks_metrics.csv")
    if os.path.exists(ref_c_csv) and os.path.exists(ref_c_metrics):
        try:
            df_c = pd.read_csv(ref_c_csv)
            theta_c = df_c['Theta (degrees)'].values
            int_c = df_c['Raw Intensity'].values
            baseline_c = df_c['Model Baseline'].values
            net_int_c = df_c['Corrected Net Intensity'].values
            
            df_m = pd.read_csv(ref_c_metrics)
            row = df_m.iloc[0]
            t0_c = row['Peak Center (Theta)']
            h_c = row['Net Height']
            fwhm_c = row['FWHM (deg)']
            w_c = fwhm_c / 2.355
            
            fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
            axes[0].plot(theta_c, int_c, '.', color='#1f77b4', label='Experimental Data', alpha=0.6)
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
            plt.savefig(os.path.join(PLOT_DIR, "calcite_single_crystal_rocking_curve_analysis.png"), dpi=150)
            plt.close()
            print("  Saved Figure A5 to results/figures/")
        except Exception as e:
            print(f"  Error plotting Figure A5: {e}")

    # --------------------------------------------------------------------------
    # FIGURE A6: DEMONSTRATION OF SYMMETRIC DIFFRACTION PEAK FITTING SEQUENCE
    # --------------------------------------------------------------------------
    print("Generating Appendix Figure A6: Symmetric peak fitting demonstration...")
    g_30_files = glob.glob(os.path.join(PROCESSED_DIR, "Symmetric_Scans/SH-125-G/*_2Theta_30_exported.xy"))
    if g_30_files:
        try:
            arr = np.loadtxt(g_30_files[0], skiprows=1)
            twotheta = arr[:, 0]
            intensity = arr[:, 1]
            
            # Subset data strictly to user-requested range [27.5, 34.5] degrees
            mask = (twotheta >= 27.5) & (twotheta <= 34.5)
            twotheta = twotheta[mask]
            intensity = intensity[mask]
            
            # Joint fitting function with a linear background baseline
            def joint_fit_linear(t, c0, c1, *peak_params):
                baseline = c0 + c1*t
                y = baseline.copy()
                n = len(peak_params) // 3
                for i in range(n):
                    h = peak_params[3*i]
                    t0 = peak_params[3*i+1]
                    w = peak_params[3*i+2]
                    y += h * np.exp(-(t - t0)**2 / (2 * w**2))
                return y
                
            # Initial guess for background by fitting to the extremes of the subset
            bg_mask = (twotheta < 28.0) | (twotheta > 34.2)
            poly_coeff = np.polyfit(twotheta[bg_mask], intensity[bg_mask], 1)
            c1_guess, c0_guess = poly_coeff
            
            p0 = [
                c0_guess, c1_guess,   # Background baseline
                2000.0, 28.5, 0.15,   # Peak 1 (28.5)
                18000.0, 29.4, 0.15,  # Peak 2 (29.4, Calcite 104)
                1000.0, 31.5, 0.15,   # Peak 3 (31.5)
                1000.0, 32.8, 0.15,   # Peak 4 (32.8, Vaterite 110)
                1000.0, 33.8, 0.15    # Peak 5 (33.8)
            ]
            
            bounds = (
                [-np.inf, -np.inf,
                 0.0, 28.2, 0.05,
                 0.0, 29.2, 0.05,
                 0.0, 31.2, 0.05,
                 0.0, 32.4, 0.05,
                 0.0, 33.5, 0.05],
                [np.inf, np.inf,
                 20000.0, 28.8, 0.5,
                 50000.0, 29.6, 0.5,
                 20000.0, 31.8, 0.5,
                 20000.0, 33.2, 0.5,
                 20000.0, 33.95, 0.25]
            )
            
            popt, _ = curve_fit(joint_fit_linear, twotheta, intensity, p0=p0, bounds=bounds)
            
            c0_fit, c1_fit = popt[0], popt[1]
            baseline = c0_fit + c1_fit * twotheta
            net_intensity = intensity - baseline
            
            # Extract fitted peak parameters
            h1, t1, w1 = popt[2], popt[3], popt[4]
            h2, t2, w2 = popt[5], popt[6], popt[7]
            h3, t3, w3 = popt[8], popt[9], popt[10]
            h4, t4, w4 = popt[11], popt[12], popt[13]
            h5, t5, w5 = popt[14], popt[15], popt[16]
            
            # Individual peaks
            fit_p1 = gaussian(twotheta, h1, t1, w1)
            fit_calcite = gaussian(twotheta, h2, t2, w2)
            fit_p2 = gaussian(twotheta, h3, t3, w3)
            fit_vaterite = gaussian(twotheta, h4, t4, w4)
            fit_p3 = gaussian(twotheta, h5, t5, w5)
            
            total_fit = baseline + fit_p1 + fit_calcite + fit_p2 + fit_vaterite + fit_p3
            sum_peaks = fit_p1 + fit_calcite + fit_p2 + fit_vaterite + fit_p3
            
            # Areas
            area_p1 = h1 * w1 * np.sqrt(2 * np.pi)
            area_calcite = h2 * w2 * np.sqrt(2 * np.pi)
            area_p2 = h3 * w3 * np.sqrt(2 * np.pi)
            area_vaterite = h4 * w4 * np.sqrt(2 * np.pi)
            area_p3 = h5 * w5 * np.sqrt(2 * np.pi)
            
            fig, axes = plt.subplots(2, 1, figsize=(10, 9), sharex=True)
            
            # Panel (a): Raw data and baseline
            axes[0].plot(twotheta, intensity, 'o', color='#7f7f7f', markersize=4, alpha=0.6, label='Raw Symmetric Scan Data')
            axes[0].plot(twotheta, baseline, 'r-', linewidth=2, label='Fitted Linear Background Baseline')
            axes[0].plot(twotheta, total_fit, 'k--', linewidth=1.5, label='Total Fit (Background + 5 Peaks)')
            axes[0].set_ylabel('Intensity (counts)', fontsize=12)
            axes[0].legend(loc='upper right', framealpha=0.9)
            axes[0].grid(True, linestyle=':', alpha=0.5)
            axes[0].set_title('Symmetric Scan Peak Fitting Sequence (SH-125-G, $\phi = 30^\circ$)', fontsize=14, fontweight='bold')
            axes[0].text(-0.08, 1.05, "(a)", transform=axes[0].transAxes, fontsize=14, fontweight='bold', va='top')
            
            # Panel (b): Net intensity and Gaussian deconvolution
            axes[1].plot(twotheta, net_intensity, 'o', color='#7f7f7f', markersize=4, alpha=0.6, label='Experimental Net Data (Raw - Background)')
            axes[1].plot(twotheta, sum_peaks, 'k-', linewidth=2, label='Sum of Fitted Peak Components')
            
            # Plot individual peaks
            axes[1].plot(twotheta, fit_p1, color='#1f77b4', linewidth=1.5, linestyle='--', label='Peak 28.5°')
            axes[1].fill_between(twotheta, 0, fit_p1, color='#1f77b4', alpha=0.12)
            
            axes[1].plot(twotheta, fit_calcite, color='#2ca02c', linewidth=1.8, linestyle='--', label='Fitted Calcite (104) ~29.4°')
            axes[1].fill_between(twotheta, 0, fit_calcite, color='#2ca02c', alpha=0.15)
            
            axes[1].plot(twotheta, fit_p2, color='#ff7f0e', linewidth=1.5, linestyle='--', label='Peak 31.5°')
            axes[1].fill_between(twotheta, 0, fit_p2, color='#ff7f0e', alpha=0.12)
            
            axes[1].plot(twotheta, fit_vaterite, color='#9467bd', linewidth=1.8, linestyle='--', label='Fitted Vaterite (110) ~32.8°')
            axes[1].fill_between(twotheta, 0, fit_vaterite, color='#9467bd', alpha=0.15)
            
            axes[1].plot(twotheta, fit_p3, color='#e377c2', linewidth=1.5, linestyle='--', label='Peak 33.8°')
            axes[1].fill_between(twotheta, 0, fit_p3, color='#e377c2', alpha=0.12)
            
            axes[1].axhline(0, color='gray', linestyle='--', linewidth=0.8)
            axes[1].set_xlabel('2$\\theta$ (°)', fontsize=12)
            axes[1].set_ylabel('Net Residual Intensity (counts)', fontsize=12)
            axes[1].set_xlim(27.5, 34.5)
            axes[1].legend(loc='upper right', framealpha=0.9, fontsize=9.5)
            axes[1].grid(True, linestyle=':', alpha=0.5)
            axes[1].text(-0.08, 1.05, "(b)", transform=axes[1].transAxes, fontsize=14, fontweight='bold', va='top')
            
            # Add text box with fit results
            fit_text = (
                f"Fit Results (Joint Fit):\n"
                f"  Peak 28.5°: Center = {t1:.3f}°, FWHM = {w1*2.355:.3f}°, Area = {area_p1/1e3:.2f} kcts·°\n"
                f"  Calcite (104): Center = {t2:.3f}°, FWHM = {w2*2.355:.3f}°, Area = {area_calcite/1e3:.2f} kcts·°\n"
                f"  Peak 31.5°: Center = {t3:.3f}°, FWHM = {w3*2.355:.3f}°, Area = {area_p2/1e3:.2f} kcts·°\n"
                f"  Vaterite (110): Center = {t4:.3f}°, FWHM = {w4*2.355:.3f}°, Area = {area_vaterite/1e3:.2f} kcts·°\n"
                f"  Peak 33.8°: Center = {t5:.3f}°, FWHM = {w5*2.355:.3f}°, Area = {area_p3/1e3:.2f} kcts·°"
            )
            # Find max net intensity to place text safely
            max_net_val = np.max(net_intensity)
            axes[1].text(27.6, max_net_val*0.42, fit_text, fontsize=8.5, 
                         bbox=dict(facecolor='white', alpha=0.85, boxstyle='round,pad=0.4'))
            
            plt.tight_layout()
            plt.savefig(os.path.join(PLOT_DIR, "fig_a6_symmetric_peak_fits.png"), dpi=300, bbox_inches='tight')
            plt.savefig(os.path.join(PLOT_DIR, "fig_a6_symmetric_peak_fits.svg"), dpi=300, bbox_inches='tight')
            plt.close()
            print("  Saved Figure A6 to results/figures/")
        except Exception as e:
            print(f"  Error plotting Figure A6: {e}")

    # --------------------------------------------------------------------------
    # FIGURES A7-A10: COMPREHENSIVE PEAK FITS FOR ALL ROCKING CURVES
    # --------------------------------------------------------------------------
    print("Generating Appendix Figures A7-A10: Rocking curve fits for all samples...")
    metrics_path = os.path.join(PROCESSED_DIR, "Rocking_Curves/Reference/all_samples_rocking_peaks_vs_phi.csv")
    if os.path.exists(metrics_path):
        try:
            df_metrics = pd.read_csv(metrics_path)
            
            fit_configs = {
                "SH-124-B3": {
                    "phi_values": [0, 30, 60, 90, 120, 150, 180],
                    "grid": (4, 2),
                    "fig_size": (12, 16),
                    "fig_name": "fig_a7_sh124b3_fits",
                    "title": "Appendix Figure A7: Rocking Curve Fits for Sample SH-124-B3"
                },
                "SH-125-A": {
                    "phi_values": [0, 30, 60, 90, 120, 150],
                    "grid": (3, 2),
                    "fig_size": (12, 12),
                    "fig_name": "fig_a8_sh125a_fits",
                    "title": "Appendix Figure A8: Rocking Curve Fits for Sample SH-125-A"
                },
                "SH-125-G": {
                    "phi_values": [0, 30, 60, 120, 150, 180],
                    "grid": (3, 2),
                    "fig_size": (12, 12),
                    "fig_name": "fig_a9_sh125g_fits",
                    "title": "Appendix Figure A9: Rocking Curve Fits for Sample SH-125-G"
                },
                "SH-104-1": {
                    "phi_values": [0, 30, 60, 90, 120, 150],
                    "grid": (3, 2),
                    "fig_size": (12, 12),
                    "fig_name": "fig_a10_sh1041_fits",
                    "title": "Appendix Figure A10: Rocking Curve Fits for Sample SH-104-1"
                }
            }
            
            for sample, info in fit_configs.items():
                # Pre-calculate global min/max of net intensity for this sample across all phis
                global_min = 0.0
                global_max = 0.0
                for phi in info["phi_values"]:
                    csv_path = os.path.join(PROCESSED_DIR, f"Rocking_Curves/{sample}/{sample}_corrected_rocking_{phi}.csv")
                    if os.path.exists(csv_path):
                        df_rc = pd.read_csv(csv_path)
                        intensity = df_rc["Raw Intensity"].values
                        baseline = df_rc["Model Baseline"].values
                        net_intensity = intensity - baseline
                        global_min = min(global_min, net_intensity.min())
                        global_max = max(global_max, net_intensity.max())
                
                # Dynamic range with padding
                padding = 0.1 * (global_max - global_min)
                ymin = global_min - padding
                ymax = global_max + padding
                # Ensure minimum scale range to look good even for flat curves
                if ymax - ymin < 500:
                    ymin, ymax = -250, 250

                rows, cols = info["grid"]
                fig, axes = plt.subplots(rows, cols, figsize=info["fig_size"], sharex=True)
                axes_flat = axes.flatten()
                
                for idx_phi, phi in enumerate(info["phi_values"]):
                    ax = axes_flat[idx_phi]
                    csv_path = os.path.join(PROCESSED_DIR, f"Rocking_Curves/{sample}/{sample}_corrected_rocking_{phi}.csv")
                    if not os.path.exists(csv_path):
                        ax.text(0.5, 0.5, f"No data for $\\phi$ = {phi}°", ha='center', va='center')
                        continue
                    

                    df_rc = pd.read_csv(csv_path)
                    theta = df_rc["Theta (degrees)"].values
                    intensity = df_rc["Raw Intensity"].values
                    baseline = df_rc["Model Baseline"].values
                    net_intensity = intensity - baseline
                    
                    # Reconstruct fit (in net intensity space)
                    fit_net_intensity = np.zeros_like(theta)
                    
                    # Get peak metrics for this sample and phi
                    df_p = df_metrics[(df_metrics["Sample"] == sample) & (df_metrics["Phi (degrees)"] == phi)]
                    
                    ax.plot(theta, net_intensity, '.', color='#7f7f7f', markersize=3, alpha=0.5, label='Experimental Net Data' if idx_phi == 0 else '')
                    ax.plot(theta, np.zeros_like(theta), 'r--', linewidth=1.2, label='Subtracted Baseline (y=0)' if idx_phi == 0 else '')
                    
                    added_labels = set()
                    for _, row_p in df_p.iterrows():
                        h = row_p["Net Height"]
                        t0 = row_p["Peak Center (Theta)"]
                        fwhm = row_p["FWHM (degrees)"]
                        if h > 0 and fwhm > 0:
                            w = fwhm / 2.355
                            y_peak = h * np.exp(-(theta - t0)**2 / (2 * w**2))
                            fit_net_intensity += y_peak
                            
                            # Plot individual peak
                            is_tilt = "Tilt" in row_p["Peak Name"] or "Tilted" in row_p["Peak Name"]
                            p_color = '#2ca02c' if is_tilt else '#1f77b4'
                            p_name = 'Fitted Tilt Component' if is_tilt else 'Fitted Specular Component'
                            if idx_phi == 0 and p_name not in added_labels:
                                p_label = p_name
                                added_labels.add(p_name)
                            else:
                                p_label = ''
                            ax.plot(theta, y_peak, color=p_color, linewidth=0.8, linestyle=':', label=p_label)
                            ax.fill_between(theta, 0, y_peak, color=p_color, alpha=0.08)
                    ax.plot(theta, fit_net_intensity, 'k-', linewidth=1.5, label='Total Net Fit' if idx_phi == 0 else '')
                    
                    # Set y-limits to the pre-calculated uniform range for this sample
                    ax.set_ylim(ymin, ymax)
                    ax.set_title(f"$\\phi$ = {phi}°", fontweight='bold', fontsize=11)
                    ax.grid(True, linestyle=':', alpha=0.4)
                    
                    if idx_phi % cols == 0:
                        ax.set_ylabel('Net Intensity (counts)')
                    if idx_phi >= (rows - 1) * cols or idx_phi == len(info["phi_values"]) - 1:
                        ax.set_xlabel('Theta $\\theta$ (°)')
                        
                # Hide unused axes
                for idx_ax in range(len(info["phi_values"]), len(axes_flat)):
                    fig.delaxes(axes_flat[idx_ax])
                    
                # Add legend to the first subplot
                axes_flat[0].legend(loc='lower left', fontsize=9, framealpha=0.8)
                
                plt.suptitle(info["title"], fontsize=14, fontweight='bold', y=0.98)
                plt.tight_layout()
                plt.savefig(os.path.join(PLOT_DIR, f"{info['fig_name']}.png"), dpi=150, bbox_inches='tight')
                plt.savefig(os.path.join(PLOT_DIR, f"{info['fig_name']}.svg"), dpi=150, bbox_inches='tight')
                plt.close()
                print(f"  Saved {info['fig_name']} to results/figures/")
        except Exception as e:
            print(f"  Error plotting comprehensive rocking curve fits: {e}")

    print("\nSTAGE 2 COMPLETE: All publication figures generated successfully in results/figures/ directory.")

if __name__ == "__main__":
    print("======================================================================")
    print("    CaCO3 Thin Film Diffraction Processing & Plotting Pipeline")
    print("======================================================================")
    
    # Run Stage 1 (Raw Data Processing)
    run_data_processing()
    
    # Run Stage 2 (Figure Generation)
    generate_all_plots()
    
    print("\n======================================================================")
    print("Process successfully finished!")
    print("Output data:      data/processed/")
    print("Output figures:   results/figures/")
    print("======================================================================")
