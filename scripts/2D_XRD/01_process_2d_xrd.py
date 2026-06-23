#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
2D-XRD Detector Frame Processor
===============================
Author: Tomasz Stawski, tomasz.stawski@bam.de
Date: 2026-06-23
Version: 1.0

Description:
Converts raw flat-panel detector frames (.gfrm) into polar cake plots,
azimuthally integrates them to produce 1D diffraction profiles, matches
peaks against standard Calcite and Vaterite reference databases, and
quantifies the preferred orientation (graininess, preferred growth).
"""

import numpy as np
import h5py
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from scipy.signal import find_peaks, savgol_filter
import scipy.ndimage as ndimage
import pandas as pd
import imageio.v2 as imageio
from PIL import Image
from silx.io.convert import convert
import os
import datetime
import argparse
from pathlib import Path
import multiprocessing

plt.rcParams.update({'font.size': 14})

# Directories
script_dir = os.path.dirname(os.path.abspath(__file__))
base_dir = os.path.abspath(os.path.join(script_dir, "..", ".."))

def auto_contrast(image, lower=1, upper=99):
    """Calculates min and max values for image display contrast based on percentiles."""
    vmin = np.percentile(image, lower)
    vmax = np.percentile(image, upper)
    return vmin, vmax

def load_reference_peaks(filepath, min_intensity=2.0):
    """Parses crystallographic reference data from a text file."""
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

def process_single_file(input_gfrm, timestamp, calcite_refs, vaterite_refs):
    """Processes a single GFRM file: conversion, integration, baseline, and matching."""
    input_gfrm = os.path.abspath(input_gfrm)
    converted_h5 = input_gfrm + '.h5'
    if os.path.exists(converted_h5):
        os.remove(converted_h5)
    
    print(f"Converting {os.path.basename(input_gfrm)} to HDF5...")
    convert(input_gfrm, converted_h5)

    # Setup Output Directory
    path_parts = os.path.normpath(input_gfrm).split(os.sep)
    if len(path_parts) >= 3:
        sample_name = f"{path_parts[-3]}_{path_parts[-2]}"
    elif len(path_parts) == 2:
        sample_name = path_parts[-2]
    else:
        sample_name = os.path.splitext(path_parts[-1])[0]
        
    output_dir = os.path.join(base_dir, "data/processed/2D_XRD/output", f"{sample_name}_{timestamp}")
    os.makedirs(output_dir, exist_ok=True)

    # Extract metadata from HDF5
    with h5py.File(converted_h5, 'r') as file:
        start_angle = float(file['scan_0/instrument/detector_0/others/ANGLES'][()][0].decode().strip().split()[0])
        increment = float(file['scan_0/instrument/detector_0/others/INCREME'][()][0])
        ncols = int(file['scan_0/instrument/detector_0/others/NCOLS'][()][0].decode().strip().split()[0])
        ending_angle = float(file['scan_0/instrument/detector_0/others/ENDING'][()][0].decode().strip().split()[0])
        detector_image = file['scan_0/instrument/detector_0/data'][()]

    # Confirm ending angle
    calculated_end = start_angle + (ncols - 1) * increment
    assert np.isclose(calculated_end, ending_angle, atol=0.05), "ENDING mismatch."

    pixel_size = 0.075  # mm
    D = 305.809         # mm
    num_y, num_x = detector_image.shape
    y_mid = (num_y - 1) / 2

    # Map pixels to angular grid
    two_theta_x = np.linspace(start_angle, calculated_end, ncols)
    y_pixel_pos = (np.arange(num_y) - y_mid) * pixel_size
    two_theta_y = np.degrees(np.arctan2(y_pixel_pos, D))
    two_theta_X_grid, two_theta_Y_grid = np.meshgrid(two_theta_x, two_theta_y)
    two_theta_total = np.sqrt(two_theta_X_grid**2 + two_theta_Y_grid**2)
    phi = np.degrees(np.arctan2(two_theta_Y_grid, two_theta_X_grid))

    excel_path = os.path.join(output_dir, 'output_data.xlsx')
    with pd.ExcelWriter(excel_path, engine='xlsxwriter') as writer:
        pd.DataFrame({
            'Parameter': ['Start_Angle', 'Increment', 'NCOLS', 'Ending_Angle', 'Calculated_Ending'],
            'Value': [start_angle, increment, ncols, ending_angle, calculated_end]
        }).to_excel(writer, sheet_name='Scan_Metadata', index=False)

        pd.DataFrame(detector_image, index=two_theta_y, columns=two_theta_x).to_excel(writer, sheet_name='2D_Original_2ThetaXY')

        # Plot raw detector image
        plt.figure(figsize=(10, 6))
        vmin, vmax = auto_contrast(detector_image)
        plt.imshow(detector_image, extent=[two_theta_x.min(), two_theta_x.max(), two_theta_y.min(), two_theta_y.max()],
                   aspect='auto', origin='lower', cmap='inferno', vmin=vmin, vmax=vmax)
        plt.xlabel('2θ_X (degrees)')
        plt.ylabel('2θ_Y (degrees)')
        plt.title('Original Detector Data in (2θ_X, 2θ_Y)')
        plt.colorbar(label='Intensity')
        plt.tight_layout()
        plot_1_path = os.path.join(output_dir, 'plot_1_original_2D.png')
        plt.savefig(plot_1_path, dpi=100)
        plt.close()

        # Resample to Cake Plot
        max_phi = max(abs(np.degrees(np.arctan2((num_y - y_mid - 1) * pixel_size, D))),
                      abs(np.degrees(np.arctan2((0 - y_mid) * pixel_size, D))))
        phi_lin = np.linspace(-max_phi, max_phi, 500)
        theta_lin = two_theta_x
        theta_grid, phi_grid = np.meshgrid(theta_lin, phi_lin)
        points = np.vstack((two_theta_total.ravel(), phi.ravel())).T
        values = detector_image.ravel()
        cake_plot = griddata(points, values, (theta_grid, phi_grid), method='linear', fill_value=np.nan)

        pd.DataFrame(cake_plot, index=phi_lin, columns=theta_lin).to_excel(writer, sheet_name='Cake_Plot_Data')

        plt.figure(figsize=(10, 6))
        plt.imshow(cake_plot, extent=[theta_lin.min(), theta_lin.max(), phi_lin.min(), phi_lin.max()],
                   aspect='auto', origin='lower', cmap='inferno', vmin=vmin, vmax=vmax)
        plt.xlabel('2θ (degrees)')
        plt.ylabel('Azimuthal Angle φ (degrees)')
        plt.title('Cake Plot (Converted from Original Data)')
        plt.colorbar(label='Intensity')
        plt.tight_layout()
        plot_2_path = os.path.join(output_dir, 'plot_2_cake_plot.png')
        plt.savefig(plot_2_path, dpi=100)
        plt.close()

        # 1D profile integration
        integrated_profile = np.nansum(cake_plot, axis=0)
        pd.DataFrame({'2Theta': theta_lin, 'Integrated Intensity': integrated_profile}).to_excel(writer, sheet_name='Integrated_Profile')

        # Baseline correction
        smoothed_profile = savgol_filter(integrated_profile, window_length=31, polyorder=3)
        window_size = 151
        bg_min = ndimage.minimum_filter1d(smoothed_profile, size=window_size)
        bg_max = ndimage.maximum_filter1d(bg_min, size=window_size)
        baseline = savgol_filter(bg_max, window_length=101, polyorder=2)
        corrected_profile = smoothed_profile - baseline
        corrected_profile[corrected_profile < 0] = 0

        pd.DataFrame({
            '2Theta': theta_lin,
            'Smoothed Profile': smoothed_profile,
            'Baseline': baseline,
            'Corrected Profile': corrected_profile
        }).to_excel(writer, sheet_name='Corrected_Profiles')

        # Peak detection
        margin = 0.5
        peak_region_mask = (theta_lin > (start_angle + margin)) & (theta_lin < (calculated_end - margin))
        region_max = np.max(corrected_profile[peak_region_mask])
        dynamic_prominence = region_max * 0.02
        peaks, _ = find_peaks(corrected_profile[peak_region_mask], prominence=dynamic_prominence, distance=15)
        peaks = peaks + np.where(peak_region_mask)[0][0]

        # Phase determination
        experimental_peaks = theta_lin[peaks]
        tolerance = 0.3
        
        calcite_matched = []
        vaterite_matched = []
        for ep in experimental_peaks:
            if any(abs(rp[0] - ep) <= tolerance for rp in calcite_refs):
                calcite_matched.append(ep)
            if any(abs(rp[0] - ep) <= tolerance for rp in vaterite_refs):
                vaterite_matched.append(ep)

        phase_determination = []
        if len(calcite_matched) > 0:
            phase_determination.append(f"Calcite ({len(calcite_matched)} peaks matched)")
        if len(vaterite_matched) > 0:
            phase_determination.append(f"Vaterite ({len(vaterite_matched)} peaks matched)")
            
        phase_result = " and ".join(phase_determination) if phase_determination else "Unknown Phase"

        pd.DataFrame({
            'Experimental Peak (2Theta)': experimental_peaks,
            'Is Calcite': [any(abs(rp[0] - ep) <= tolerance for rp in calcite_refs) for ep in experimental_peaks],
            'Is Vaterite': [any(abs(rp[0] - ep) <= tolerance for rp in vaterite_refs) for ep in experimental_peaks]
        }).to_excel(writer, sheet_name='Phase_Matching', index=False)
        pd.DataFrame({'Detected Phase': [phase_result]}).to_excel(writer, sheet_name='Phase_Result', index=False)

        # Plot combined cake & integrated profile
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10), sharex=True)
        ax1.imshow(cake_plot, extent=[theta_lin.min(), theta_lin.max(), phi_lin.min(), phi_lin.max()],
                   aspect='auto', origin='lower', cmap='inferno', vmin=vmin, vmax=vmax)
        ax1.set_ylabel('Azimuthal Angle φ (degrees)')
        ax1.set_title('Cake Plot with Detected Peaks')
        ax2.plot(theta_lin, integrated_profile, color='blue', alpha=0.4, label='Integrated Profile')
        ax2.plot(theta_lin, corrected_profile, color='darkorange', label='Corrected Profile')
        ax2.plot(theta_lin[peaks], corrected_profile[peaks], 'rx', label='Detected Peaks')
        for peak_pos in theta_lin[peaks]:
            ax2.axvspan(peak_pos - 0.5, peak_pos + 0.5, color='red', alpha=0.2)
            
        calcite_label_added = False
        y_max = np.nanmax(corrected_profile)
        for rp in calcite_refs:
            if start_angle <= rp[0] <= calculated_end:
                ax2.axvline(rp[0], color='green', linestyle='--', alpha=0.5, label='Calcite Ref' if not calcite_label_added else "")
                calcite_label_added = True
                h, k, l, d = rp[2], rp[3], rp[4], rp[5]
                ax2.text(rp[0], y_max * 0.95, f"({h}{k}{l}) {d:.2f}Å", rotation=90, color='green', va='top', ha='right', fontsize=9)
                
        vaterite_label_added = False
        for rp in vaterite_refs:
            if start_angle <= rp[0] <= calculated_end:
                ax2.axvline(rp[0], color='purple', linestyle=':', alpha=0.5, label='Vaterite Ref' if not vaterite_label_added else "")
                vaterite_label_added = True
                h, k, l, d = rp[2], rp[3], rp[4], rp[5]
                ax2.text(rp[0], y_max * 0.95, f"({h}{k}{l}) {d:.2f}Å", rotation=90, color='purple', va='top', ha='left', fontsize=9)

        ax2.set_xlabel('2θ (degrees)')
        ax2.set_title(f'Detected Phase: {phase_result}')
        ax2.set_ylabel('Intensity')
        ax2.grid(True)
        ax2.legend()
        plt.tight_layout()
        plot_3_path = os.path.join(output_dir, 'plot_3_combined.png')
        plt.savefig(plot_3_path, dpi=100)
        plt.close()

        # Phase-specific plotting
        images_for_gif = [plot_1_path, plot_2_path, plot_3_path]
        if len(calcite_matched) > 0:
            plot_4_path = os.path.join(output_dir, 'plot_4_calcite_only.png')
            shutil_copy = True
            shutil_src = plot_3_path
            shutil_dst = plot_4_path
            import shutil
            shutil.copy2(shutil_src, shutil_dst)
            images_for_gif.append(plot_4_path)

        # Plot azimuthal profile per peak
        peak_profiles = {}
        peak_metrics = []

        for idx, peak_pos in enumerate(theta_lin[peaks]):
            peak_mask = (theta_grid[0, :] >= (peak_pos - 0.5)) & (theta_grid[0, :] <= (peak_pos + 0.5))
            azimuthal_profile = np.nanmean(cake_plot[:, peak_mask], axis=1)
            peak_profiles[f'Peak_{idx+1}_{peak_pos:.2f}'] = azimuthal_profile
            
            # Calculate Degree of Anisotropy (DoA) & Coefficient of Variation (CV)
            valid_mask = ~np.isnan(azimuthal_profile)
            DoA, CV = np.nan, np.nan
            if np.sum(valid_mask) > 10:
                clean_profile = np.copy(azimuthal_profile)
                clean_profile[~valid_mask] = np.nanmin(clean_profile[valid_mask])
                mean_I = np.nanmean(clean_profile)
                std_I = np.nanstd(clean_profile)
                if mean_I > 0:
                    CV = std_I / mean_I
                az_smooth = savgol_filter(clean_profile, window_length=5, polyorder=2)
                p_95 = np.percentile(az_smooth, 95)
                p_05 = np.percentile(az_smooth, 5)
                if p_95 > 0:
                    DoA = (p_95 - p_05) / p_95

            # Determine C-Axis
            is_c_axis = False
            matched_refs = []
            for rp in calcite_refs:
                if abs(rp[0] - peak_pos) <= tolerance:
                    matched_refs.append(f"Calcite ({rp[2]}{rp[3]}{rp[4]})")
                    if rp[2] == '0' and rp[3] == '0' and rp[4] != '0':
                        is_c_axis = True
            for rp in vaterite_refs:
                if abs(rp[0] - peak_pos) <= tolerance:
                    matched_refs.append(f"Vaterite ({rp[2]}{rp[3]}{rp[4]})")
                    if rp[2] == '0' and rp[3] == '0' and rp[4] != '0':
                        is_c_axis = True

            peak_metrics.append({
                'Peak ID': f'Peak_{idx+1}',
                '2Theta': peak_pos,
                'Matched Refs': " | ".join(matched_refs) if matched_refs else "Unknown",
                'Is C-Axis': is_c_axis,
                'Degree of Anisotropy (DoA)': DoA,
                'Coefficient of Variation (CV)': CV
            })

            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
            ax1.imshow(cake_plot, extent=[theta_lin.min(), theta_lin.max(), phi_lin.min(), phi_lin.max()],
                       aspect='auto', origin='lower', cmap='inferno', vmin=vmin, vmax=vmax)
            rect_x_start = peak_pos - 0.5
            rect = plt.Rectangle((rect_x_start, phi_lin.min()), 1.0, phi_lin.ptp(),
                                 linewidth=2, edgecolor='green', facecolor='none')
            ax1.add_patch(rect)
            ax1.set_ylabel('Azimuthal Angle φ (degrees)')
            ax1.set_title(f'Cake Plot Context for Peak at {peak_pos:.2f}° 2θ')

            ax2.plot(phi_lin, azimuthal_profile, color='purple')
            ax2.set_xlabel('Azimuthal Angle φ (degrees)')
            ax2.set_ylabel('Average Intensity')
            c_axis_str = " [C-Axis Reflection]" if is_c_axis else ""
            ax2.set_title(f'Azimuthal Profile for Peak at {peak_pos:.2f}° 2θ{c_axis_str}\nDoA: {DoA:.2f} | CV: {CV:.3f}')
            ax2.grid(True)
            plt.tight_layout()
            
            filename = os.path.join(output_dir, f'plot_peak_{idx+1}.png')
            plt.savefig(filename, dpi=100)
            images_for_gif.append(filename)
            plt.close()

        pd.DataFrame(peak_profiles, index=phi_lin).to_excel(writer, sheet_name='Azimuthal_Profiles')
        pd.DataFrame(peak_metrics).to_excel(writer, sheet_name='Orientation_Metrics', index=False)

        # Summarize orientation by phase
        calcite_cvs = [m['Coefficient of Variation (CV)'] for m in peak_metrics if 'Calcite' in m['Matched Refs'] and not np.isnan(m['Coefficient of Variation (CV)'])]
        vaterite_cvs = [m['Coefficient of Variation (CV)'] for m in peak_metrics if 'Vaterite' in m['Matched Refs'] and not np.isnan(m['Coefficient of Variation (CV)'])]
        
        avg_calcite_cv = np.mean(calcite_cvs) if calcite_cvs else np.nan
        avg_vaterite_cv = np.mean(vaterite_cvs) if vaterite_cvs else np.nan
        
        orientation_result = "Mainly Isotropic"
        if not np.isnan(avg_calcite_cv) and avg_calcite_cv > 0.05:
            orientation_result = "Calcite exhibits preferred orientation."
        elif not np.isnan(avg_vaterite_cv) and avg_vaterite_cv > 0.05:
            orientation_result = "Vaterite exhibits preferred orientation."
            
        pd.DataFrame({
            'Phase': ['Calcite', 'Vaterite'],
            'Average CV': [avg_calcite_cv, avg_vaterite_cv],
            'Conclusion': [orientation_result, '']
        }).to_excel(writer, sheet_name='Phase_Orientation_Summary', index=False)

    # Create GIF
    frames = []
    target_size = None
    for img in images_for_gif:
        if os.path.exists(img):
            frame = Image.open(img)
            if target_size is None:
                target_size = frame.size
            if frame.size != target_size:
                frame = frame.resize(target_size, Image.Resampling.LANCZOS)
            frames.append(np.array(frame))
    if frames:
        gif_path = os.path.join(output_dir, 'analysis_summary.gif')
        imageio.mimsave(gif_path, frames, duration=2)

    # Generate HTML Summary Report
    html_path = os.path.join(output_dir, 'analysis_report.html')
    with open(html_path, 'w') as f:
        f.write(f"""<!DOCTYPE html>
    <html>
    <head>
    <meta charset="UTF-8">
    <title>XRD Analysis Report: {sample_name}</title>
    <style>
        body {{ font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; line-height: 1.6; color: #333; max-width: 1000px; margin: 0 auto; padding: 20px; }}
        h1 {{ color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px; }}
        h2 {{ color: #2980b9; margin-top: 30px; border-left: 4px solid #3498db; padding-left: 10px; }}
        p {{ font-size: 1.1em; }}
        .metadata {{ background: #f8f9fa; padding: 15px; border-radius: 5px; border-left: 4px solid #3498db; margin-bottom: 20px; }}
        .image-container {{ text-align: center; margin: 20px 0; }}
        img {{ max-width: 100%; height: auto; border: 1px solid #ddd; border-radius: 4px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }}
        .conclusion {{ background: #e8f4f8; padding: 20px; border-radius: 5px; font-weight: bold; font-size: 1.2em; text-align: center; margin-top: 40px; border: 2px solid #bce8f1; }}
    </style>
    </head>
    <body>
        <h1>XRD Analysis Report</h1>
        
        <div class="metadata">
            <strong>Sample Name:</strong> {sample_name}<br>
            <strong>Processing Time:</strong> {timestamp}<br>
            <strong>Source File:</strong> {os.path.basename(input_gfrm)}
        </div>

        <h2>Stage 1: Data Conversion & Resampling</h2>
        <p>The raw 2D detector image is converted from proprietary format to standard Cartesian coordinates. It is then resampled into polar coordinates (Cake Plot) to map 2Theta against the Azimuthal angle.</p>
        <div class="image-container">
            <img src="plot_1_original_2D.png" alt="Detector Image">
        </div>
        <div class="image-container">
            <img src="plot_2_cake_plot.png" alt="Cake Plot">
        </div>

        <h2>Stage 2 & 3: Baseline Correction & Phase Identification</h2>
        <p>The Cake Plot is integrated into a 1D profile. A morphological top-hat baseline filter is applied to subtract background scatter while preserving peak shapes. Peak positions are matched against standard Calcite and Vaterite reference databases.</p>
        <div class="image-container">
            <img src="plot_3_combined.png" alt="Combined Phase Plot">
        </div>
    """)
        if len(calcite_matched) > 0:
            f.write(f"""
        <div class="image-container">
            <img src="plot_4_calcite_only.png" alt="Calcite Only Plot">
        </div>
    """)

        f.write(f"""
        <h2>Stage 4: Orientation and Texture Analysis</h2>
        <p>For each detected peak, the azimuthal intensity profile is extracted. The Degree of Anisotropy (DoA) and Coefficient of Variation (CV) are calculated to assess preferred orientation and texture characteristics.</p>
    """)
        for idx in range(len(peaks)):
            f.write(f"""
        <div class="image-container">
            <img src="plot_peak_{idx+1}.png" alt="Peak {idx+1} Profile">
        </div>
    """)

        f.write(f"""
        <div class="conclusion">
            <p>Detected Phases: {phase_result}</p>
            <p>Orientation Conclusion: {orientation_result}</p>
        </div>
    </body>
    </html>
    """)

    print(f"\nAnalysis complete! HTML report saved to: {html_path}")
    return {
        'sample_name': sample_name,
        'theta_lin': theta_lin,
        'corrected_profile': corrected_profile,
        'integrated_profile': integrated_profile,
        'phase_result': phase_result,
        'orientation_result': orientation_result
    }

def _process_wrapper(args):
    f, timestamp, calcite_refs, vaterite_refs = args
    try:
        return process_single_file(f, timestamp, calcite_refs, vaterite_refs)
    except Exception as e:
        print(f"Error processing {f}: {e}")
        return None

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process .gfrm files.')
    parser.add_argument('input_path', nargs='?', default=os.path.join(base_dir, 'data/raw/2D_XRD'), help='Path to .gfrm file or directory')
    args = parser.parse_args()

    input_path = Path(args.input_path).resolve()
    if input_path.is_file():
        gfrm_files = [str(input_path)]
    else:
        gfrm_files = [str(p) for p in input_path.rglob('*.gfrm')]

    if not gfrm_files:
        print(f"No .gfrm files found in {args.input_path}")
        exit(0)

    timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    
    crystal_data_dir = os.path.join(base_dir, 'data/raw/crystalData')
    calcite_refs = load_reference_peaks(os.path.join(crystal_data_dir, 'Calcite__0000985.txt'), min_intensity=5.0)
    vaterite_refs = load_reference_peaks(os.path.join(crystal_data_dir, 'Vaterite__0004854.txt'), min_intensity=5.0)

    print(f"Found {len(gfrm_files)} .gfrm files. Starting parallel batch processing...")
    
    # Run in parallel
    pool_args = [(f, timestamp, calcite_refs, vaterite_refs) for f in gfrm_files]
    num_cores = max(1, multiprocessing.cpu_count() - 1)
    with multiprocessing.Pool(processes=num_cores) as pool:
        results = pool.map(_process_wrapper, pool_args)
        
    comparison_data = [res for res in results if res is not None]

    if comparison_data:
        comp_dir = os.path.join(base_dir, "data/processed/2D_XRD/output", f"comparison_{timestamp}")
        os.makedirs(comp_dir, exist_ok=True)
        
        # Plot comparison of all corrected profiles
        plt.figure(figsize=(12, 8))
        for d in comparison_data:
            plt.plot(d['theta_lin'], d['corrected_profile'], label=d['sample_name'], alpha=0.7)
        plt.xlabel('2θ (degrees)')
        plt.ylabel('Intensity')
        plt.title('Comparison of Corrected Profiles')
        if len(comparison_data) <= 15:
            plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(os.path.join(comp_dir, 'comparison_profiles.png'), dpi=150)
        plt.close()

        # Save comparison data to excel
        comp_excel = os.path.join(comp_dir, 'comparison_summary.xlsx')
        with pd.ExcelWriter(comp_excel, engine='xlsxwriter') as writer:
            summary_df = pd.DataFrame([
                {'Sample Name': d['sample_name'], 
                 'Detected Phase': d['phase_result'], 
                 'Orientation': d['orientation_result']} 
                 for d in comparison_data
            ])
            summary_df.to_excel(writer, sheet_name='Summary', index=False)

            profiles_dict = {}
            for d in comparison_data:
                profiles_dict[d['sample_name'] + '_2Theta'] = pd.Series(d['theta_lin'])
                profiles_dict[d['sample_name'] + '_Intensity'] = pd.Series(d['corrected_profile'])
            profiles_df = pd.DataFrame(profiles_dict)
            profiles_df.to_excel(writer, sheet_name='Corrected_Profiles', index=False)

        print(f"\nBatch processing complete! Comparison saved to: {comp_dir}")
