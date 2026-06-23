# -*- coding: utf-8 -*-
"""
2D Polar Texture Plotter
========================
Author: Tomasz Stawski, tomasz.stawski@bam.de
Date: 2026-06-19
Version: 0.1

Description:
Maps 1D rocking curves collected at various azimuths to construct 2D polar projection pole figures.
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import zipfile
import xml.etree.ElementTree as ET
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d

# Directories
base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
raw_data_dir = os.path.join(base_dir, "data/raw/Rocking_Curves")
processed_dir = os.path.join(base_dir, "data/processed/Rocking_Curves")
analysis_dir = os.path.join(base_dir, "data/processed/Rocking_Curves/08062026")

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

# Fitting helpers
def gaussian(t, h, t0, w):
    return h * np.exp(-(t-t0)**2 / (2*w**2))

def bg_model(t, I0, c0, c1, c2, c3):
    return I0 / np.sin(np.radians(t)) + c0 + c1*t + c2*t**2 + c3*t**3

# Define samples
samples = ["SH-124-B3", "SH-125-A"]

# Background masks from config
bg_mask_configs = {
    "SH-124-B3": lambda theta: (theta >= 7.0) & \
                               ((theta < 9.0) | (theta > 9.7)) & \
                               ((theta < 12.2) | (theta > 13.6)) & \
                               ((theta < 15.0) | (theta > 16.8)),
    "SH-125-A": lambda theta: (theta < 9.5) | \
                              ((theta > 10.5) & (theta < 12.3)) | \
                              ((theta > 13.3) & (theta < 13.7)) | \
                              ((theta > 15.0) & (theta < 21.8)) | \
                              (theta > 23.2),
}

# We will generate a texture plot for each sample
for sample in samples:
    print(f"\nProcessing texture plots for {sample}...")
    sample_raw_dir = os.path.join(raw_data_dir, sample)
    sample_analysis_dir = os.path.join(analysis_dir, sample)
    
    # List files to determine Phi values
    files = os.listdir(sample_raw_dir)
    phi_values = sorted(list(set([int(f.split("_")[1].split(".")[0]) for f in files if f.endswith(".brml")])))
    
    # 1. Fit Calcite 2Theta scans to get theta_0 for each phi
    theta_0_dict = {}
    for phi in phi_values:
        twotheta_path = os.path.join(sample_raw_dir, f"2Theta_{phi}.brml")
        if os.path.exists(twotheta_path):
            arr_2t = extract_brml_data(twotheta_path)
            twotheta_2t = arr_2t[:, 2]
            intensity_2t = arr_2t[:, 4]
            
            # Fit baseline and peak
            poly_coeff = np.polyfit(twotheta_2t, intensity_2t, 3)
            baseline_2t = np.polyval(poly_coeff, twotheta_2t)
            net_intensity_2t = intensity_2t - baseline_2t
            
            c_mask = (twotheta_2t >= 28.5) & (twotheta_2t <= 30.5)
            try:
                popt_c, _ = curve_fit(gaussian, twotheta_2t[c_mask], net_intensity_2t[c_mask], 
                                      p0=[intensity_2t.max() - baseline_2t.max(), 29.4, 0.15])
                theta_0_dict[phi] = popt_c[1] / 2.0
            except:
                theta_0_dict[phi] = 14.7
        else:
            theta_0_dict[phi] = 14.7
            
    # 2. Load and baseline-subtract rocking curves
    rocking_scans = {}
    for phi in phi_values:
        rocking_path = os.path.join(sample_raw_dir, f"Rocking_{phi}.brml")
        if os.path.exists(rocking_path):
            arr_rc = extract_brml_data(rocking_path)
            theta = arr_rc[:, 2]
            intensity = arr_rc[:, 3]
            
            # Fit baseline
            bg_mask = bg_mask_configs[sample](theta)
            try:
                popt_bg, _ = curve_fit(bg_model, theta[bg_mask], intensity[bg_mask], 
                                      p0=[100000, 5e6, -3e5, 8000, -100])
                baseline = bg_model(theta, *popt_bg)
            except:
                poly_coeff = np.polyfit(theta[bg_mask], intensity[bg_mask], 3)
                baseline = np.polyval(poly_coeff, theta)
            
            net_intensity = intensity - baseline
            
            # Store interpolators for raw and net intensity
            rocking_scans[phi] = {
                "raw_interp": interp1d(theta, intensity, bounds_error=False, fill_value=0.0),
                "net_interp": interp1d(theta, net_intensity, bounds_error=False, fill_value=0.0)
            }
            
    # 3. Construct the Polar Grid
    # Azimuths (Phi) in 30 degree steps: 0, 30, 60, ..., 360
    phi_polar_deg = np.arange(0, 360 + 30, 30)
    phi_polar_rad = np.radians(phi_polar_deg)
    
    # Radii (Tilt Chi) from 0 to 10 degrees
    r_tilt = np.linspace(0, 10.0, 500)
    
    # Grid arrays for Raw and Net intensity
    Z_raw = np.zeros((len(r_tilt), len(phi_polar_deg)))
    Z_net = np.zeros((len(r_tilt), len(phi_polar_deg)))
    
    for idx_phi, phi_p in enumerate(phi_polar_deg):
        # We need to map this polar angle to the measured scan and tilt direction
        # If phi_p is 360, it wraps to 0
        phi_mapped = phi_p % 360
        
        # Let's collect data at this polar angle
        raw_profile = np.zeros_like(r_tilt)
        net_profile = np.zeros_like(r_tilt)
        
        if phi_mapped <= 150:
            # Matches the measured scan at phi_mapped directly (positive tilt side, chi >= 0)
            scan_phi = phi_mapped
            if scan_phi in rocking_scans:
                theta_0 = theta_0_dict[scan_phi]
                raw_profile = rocking_scans[scan_phi]["raw_interp"](theta_0 + r_tilt)
                net_profile = rocking_scans[scan_phi]["net_interp"](theta_0 + r_tilt)
        elif phi_mapped == 180:
            # Can be positive side of scan 180 (if exists) or negative side of scan 0 (chi <= 0)
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
        else: # phi_mapped >= 210
            # Maps to negative tilt side of scan (phi_mapped - 180)
            scan_phi = phi_mapped - 180
            if scan_phi in rocking_scans:
                theta_0 = theta_0_dict[scan_phi]
                raw_profile = rocking_scans[scan_phi]["raw_interp"](theta_0 - r_tilt)
                net_profile = rocking_scans[scan_phi]["net_interp"](theta_0 - r_tilt)
                
        Z_raw[:, idx_phi] = raw_profile
        Z_net[:, idx_phi] = net_profile
        
    # Ensure origin is smooth (average across all phi)
    Z_raw[0, :] = np.mean(Z_raw[0, :])
    Z_net[0, :] = np.mean(Z_net[0, :])
    
    # 4. Plotting raw and net texture plots
    for mode, Z_data in [("Raw", Z_raw), ("Net Residual", Z_net)]:
        fig, ax = plt.subplots(subplot_kw=dict(projection='polar'), figsize=(8, 8))
        
        # Grid mesh
        Phi_mesh, R_mesh = np.meshgrid(phi_polar_rad, r_tilt)
        
        # Use contourf to generate smooth contour plot
        # Set negative net intensities to zero for visualization
        Z_plot = np.clip(Z_data, 0, None) if "Net" in mode else Z_data
        
        contour = ax.contourf(Phi_mesh, R_mesh, Z_plot, levels=60, cmap='plasma')
        cbar = fig.colorbar(contour, ax=ax, pad=0.1, shrink=0.7)
        cbar.set_label("Intensity (counts)")
        
        # Label coordinates
        ax.set_theta_zero_location("E") # Zero at the Right (standard math orientation)
        ax.set_theta_direction(1) # Counter-clockwise
        
        # Radial limits and labels
        ax.set_ylim(0, 10.0)
        ax.set_yticks([2, 4, 6, 8, 10])
        ax.set_yticklabels(["2°", "4°", "6°", "8°", "10°"], color="grey", size=10)
        ax.set_rlabel_position(45) # put radial labels at 45 degrees
        
        # Azimuthal labels (Phi)
        ax.set_xticks(np.radians(np.arange(0, 360, 30)))
        ax.set_xticklabels([f"{d}°" for d in np.arange(0, 360, 30)], size=10)
        
        ax.grid(True, color="grey", linestyle=":", alpha=0.5)
        ax.set_title(f"Calcite (104) 2D Texture Plot ({mode} Intensity)\nSample: {sample}", y=1.08, fontsize=12, fontweight='bold')
        
        # Save figure
        mode_suffix = "raw" if "Raw" in mode else "net"
        plot_path = os.path.join(sample_analysis_dir, f"{sample}_texture_{mode_suffix}.png")
        plt.tight_layout()
        plt.savefig(plot_path, dpi=150)
        print(f"  Generated {plot_path}")
        plt.close()

print("\nTexture plots generation completed successfully.")
