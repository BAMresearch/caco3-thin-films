
"""
Created on Fri Mar 21 10:52:49 2025

@author: tomasz.stawski@bam.de

# Manual: Dependencies and Installation Instructions

## Required Dependencies
- numpy
- h5py
- matplotlib
- scipy
- pandas
- imageio
- glob (standard library)
- python-docx
- silx
- xlsxwriter (automatically used by pandas ExcelWriter)

## Installation with Anaconda

You can create a new environment and install all dependencies by running the following commands in your Anaconda prompt:

```
conda create -n cakeplot_env python=3.9
conda activate cakeplot_env
conda install numpy h5py matplotlib scipy pandas imageio
conda install -c conda-forge silx
conda install conda-forge::xlsxwriter
```

## Running the Script
1. Make sure `SH-125-A.gfrm` is in the working directory.
2. Execute the provided script in the activated `cakeplot_env` environment.
3. The script will automatically convert the `.gfrm` file to `.h5`, process the data, generate plots, save them as PNG files, create an Excel file with results, and produce both an animated GIF and a DOCX report with the figures embedded.

## Output Files
- `output_data.xlsx`: contains original data, cake plot data, integrated and corrected profiles, and azimuthal profiles.
- PNG figures (saved as `plot_1_...`, `plot_2_...`, `plot_3_...`, `plot_peak_...`)
- `analysis_summary.gif`: animated GIF showcasing results.
"""




import numpy as np
import h5py
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from scipy.signal import find_peaks, savgol_filter, peak_widths
import pandas as pd
import imageio.v2 as imageio
from PIL import Image
from silx.io.convert import convert
import os

def auto_contrast(image, lower=1, upper=99):
    vmin = np.percentile(image, lower)
    vmax = np.percentile(image, upper)
    return vmin, vmax

plt.rcParams.update({'font.size': 14})

# Convert .gfrm to HDF5 using silx
input_gfrm ='BayBionic Data/XRD_Daten_Udo/SH-125-G/Set1_20-55TT_0.02SS_20SPS/Universtitaet Erlangen Nuernberg Institut_1345_250124_091118-000.gfrm'
converted_h5 = input_gfrm + '.h5'
if os.path.exists(converted_h5):
    os.remove(converted_h5)
convert(input_gfrm, converted_h5)

# Extract metadata from HDF5
with h5py.File(converted_h5, 'r') as file:
    start_angle = float(file['scan_0/instrument/detector_0/others/ANGLES'][()][0].decode().strip().split()[0])
    increment = float(file['scan_0/instrument/detector_0/others/INCREME'][()][0])
    ncols = int(file['scan_0/instrument/detector_0/others/NCOLS'][()][0].decode().strip().split()[0])
    ending_angle = float(file['scan_0/instrument/detector_0/others/ENDING'][()][0].decode().strip().split()[0])

# Confirm ending angle
calculated_end = start_angle + (ncols - 1) * increment
assert np.isclose(calculated_end, ending_angle, atol=0.05), "ENDING mismatch."

# Load detector image
with h5py.File(converted_h5, 'r') as file:
    detector_image = file['scan_0/instrument/detector_0/data'][()]

pixel_size = 0.075  # mm
D = 305.809         # mm
num_y, num_x = detector_image.shape
y_mid = (num_y - 1) / 2

# Dynamically use metadata
two_theta_x = np.linspace(start_angle, calculated_end, ncols)
y_pixel_pos = (np.arange(num_y) - y_mid) * pixel_size
two_theta_y = np.degrees(np.arctan2(y_pixel_pos, D))
two_theta_X_grid, two_theta_Y_grid = np.meshgrid(two_theta_x, two_theta_y)
two_theta_total = np.sqrt(two_theta_X_grid**2 + two_theta_Y_grid**2)
phi = np.degrees(np.arctan2(two_theta_Y_grid, two_theta_X_grid))

with pd.ExcelWriter('output_data.xlsx', engine='xlsxwriter') as writer:
    pd.DataFrame({
        'Parameter': ['Start_Angle', 'Increment', 'NCOLS', 'Ending_Angle', 'Calculated_Ending'],
        'Value': [start_angle, increment, ncols, ending_angle, calculated_end]
    }).to_excel(writer, sheet_name='Scan_Metadata', index=False)

    pd.DataFrame(detector_image, index=two_theta_y, columns=two_theta_x).to_excel(writer, sheet_name='2D_Original_2ThetaXY')

    plt.figure(figsize=(10, 6))
    vmin, vmax = auto_contrast(detector_image)
    plt.imshow(detector_image, extent=[two_theta_x.min(), two_theta_x.max(), two_theta_y.min(), two_theta_y.max()],
               aspect='auto', origin='lower', cmap='inferno', vmin=vmin, vmax=vmax)
    plt.xlabel('2θ_X (degrees)')
    plt.ylabel('2θ_Y (degrees)')
    plt.title('Original Detector Data in (2θ_X, 2θ_Y)')
    plt.colorbar(label='Intensity')
    plt.tight_layout()
    plt.savefig('plot_1_original_2D.png', dpi=100)
    plt.show()

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
    plt.savefig('plot_2_cake_plot.png', dpi=100)
    plt.show()

    integrated_profile = np.nansum(cake_plot, axis=0)
    pd.DataFrame({'2Theta': theta_lin, 'Integrated Intensity': integrated_profile}).to_excel(writer, sheet_name='Integrated_Profile')

    smoothed_profile = savgol_filter(integrated_profile, window_length=31, polyorder=3)
    baseline = savgol_filter(smoothed_profile, window_length=201, polyorder=2)
    corrected_profile = smoothed_profile - baseline

    pd.DataFrame({
        '2Theta': theta_lin,
        'Smoothed Profile': smoothed_profile,
        'Baseline': baseline,
        'Corrected Profile': corrected_profile
    }).to_excel(writer, sheet_name='Corrected_Profiles')

    # Exclude edges for peak detection only
    margin = 0.5
    peak_region_mask = (theta_lin > (start_angle + margin)) & (theta_lin < (calculated_end - margin))
    peaks, _ = find_peaks(corrected_profile[peak_region_mask], prominence=1500000)
    peaks = peaks + np.where(peak_region_mask)[0][0]

    # Load reference data
    def load_reference_peaks(filepath, min_intensity=2.0):
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

    calcite_refs = load_reference_peaks('crystalData/Calcite__0000985.txt', min_intensity=5.0)
    vaterite_refs = load_reference_peaks('crystalData/Vaterite__0004854.txt', min_intensity=5.0)

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
        
    if not phase_determination:
        phase_result = "Unknown Phase"
    else:
        phase_result = " and ".join(phase_determination)

    pd.DataFrame({
        'Experimental Peak (2Theta)': experimental_peaks,
        'Is Calcite': [any(abs(rp[0] - ep) <= tolerance for rp in calcite_refs) for ep in experimental_peaks],
        'Is Vaterite': [any(abs(rp[0] - ep) <= tolerance for rp in vaterite_refs) for ep in experimental_peaks]
    }).to_excel(writer, sheet_name='Phase_Matching', index=False)
    
    pd.DataFrame({'Detected Phase': [phase_result]}).to_excel(writer, sheet_name='Phase_Result', index=False)

    # Plot combined cake and integrated profile with peaks
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
    plt.savefig('plot_3_combined.png', dpi=100)
    plt.show()

    # Save azimuthal profiles per peak
    images_for_gif = ['plot_1_original_2D.png', 'plot_2_cake_plot.png', 'plot_3_combined.png']
    peak_profiles = {}
    peak_metrics = []

    for idx, peak_pos in enumerate(theta_lin[peaks]):
        peak_mask = (theta_grid[0, :] >= (peak_pos - 0.5)) & (theta_grid[0, :] <= (peak_pos + 0.5))
        azimuthal_profile = np.nanmean(cake_plot[:, peak_mask], axis=1)
        peak_profiles[f'Peak_{idx+1}_{peak_pos:.2f}'] = azimuthal_profile
        
        # Calculate Orientation Metrics
        # Because the diffraction rings are mostly isotropic (powder-like) but may contain
        # localized single-crystal or columnar grains, we calculate specific metrics to 
        # quantify these intensity fluctuations:
        # 
        # 1. Coefficient of Variation (CV): std(I) / mean(I).
        #    A perfectly smooth, isotropic powder ring will have a CV near 0. If the ring
        #    contains sharp localized single-crystal reflections, the variance increases relative
        #    to the background mean. This provides a continuous metric for localized intensity variations.
        # 2. Degree of Anisotropy (DoA): (I_95th - I_5th) / I_95th.
        #    A robust metric indicating how much the intensity varies across the ring overall.
        
        valid_mask = ~np.isnan(azimuthal_profile)
        DoA = np.nan
        CV = np.nan
        if np.sum(valid_mask) > 10:
            clean_profile = np.copy(azimuthal_profile)
            clean_profile[~valid_mask] = np.nanmin(clean_profile[valid_mask])
            
            mean_I = np.nanmean(clean_profile)
            std_I = np.nanstd(clean_profile)
            
            if mean_I > 0:
                CV = std_I / mean_I
                
            # Light smoothing to preserve sharp single-crystal discrete reflections but remove basic noise
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
        filename = f'plot_peak_{idx+1}.png'
        plt.savefig(filename, dpi=100)
        images_for_gif.append(filename)
        plt.show()

    pd.DataFrame(peak_profiles, index=phi_lin).to_excel(writer, sheet_name='Azimuthal_Profiles')
    pd.DataFrame(peak_metrics).to_excel(writer, sheet_name='Orientation_Metrics', index=False)

    # Summarize orientation by phase
    calcite_cvs = [m['Coefficient of Variation (CV)'] for m in peak_metrics if 'Calcite' in m['Matched Refs'] and not np.isnan(m['Coefficient of Variation (CV)'])]
    vaterite_cvs = [m['Coefficient of Variation (CV)'] for m in peak_metrics if 'Vaterite' in m['Matched Refs'] and not np.isnan(m['Coefficient of Variation (CV)'])]
    
    avg_calcite_cv = np.mean(calcite_cvs) if calcite_cvs else np.nan
    avg_vaterite_cv = np.mean(vaterite_cvs) if vaterite_cvs else np.nan
    
    orientation_result = "Mainly Isotropic"
    
    if not np.isnan(avg_calcite_cv) and not np.isnan(avg_vaterite_cv):
        if avg_calcite_cv > 0.05 and avg_calcite_cv > avg_vaterite_cv * 1.5:
            orientation_result = "Calcite exhibits strong preferred orientation."
        elif avg_vaterite_cv > 0.05 and avg_vaterite_cv > avg_calcite_cv * 1.5:
            orientation_result = "Vaterite exhibits strong preferred orientation."
        elif avg_calcite_cv > 0.05 or avg_vaterite_cv > 0.05:
            orientation_result = "Both phases exhibit preferred orientation."
    elif not np.isnan(avg_calcite_cv) and avg_calcite_cv > 0.05:
        orientation_result = "Calcite exhibits preferred orientation."
    elif not np.isnan(avg_vaterite_cv) and avg_vaterite_cv > 0.05:
        orientation_result = "Vaterite exhibits preferred orientation."
        
    pd.DataFrame({
        'Phase': ['Calcite', 'Vaterite'],
        'Average CV': [avg_calcite_cv, avg_vaterite_cv],
        'Conclusion': [orientation_result, '']
    }).to_excel(writer, sheet_name='Phase_Orientation_Summary', index=False)
    
    print("\n--- Orientation Analysis Summary ---")
    print(f"Average Calcite CV: {avg_calcite_cv:.4f}")
    print(f"Average Vaterite CV: {avg_vaterite_cv:.4f}")
    print(f"Conclusion: {orientation_result}")
    print("------------------------------------\n")

# Create animated GIF
frames = []
target_size = None
for img in images_for_gif:
    frame = Image.open(img)
    if target_size is None:
        target_size = frame.size
    if frame.size != target_size:
        frame = frame.resize(target_size, Image.Resampling.LANCZOS)
    frames.append(np.array(frame))
imageio.mimsave('analysis_summary.gif', frames, duration=2)
