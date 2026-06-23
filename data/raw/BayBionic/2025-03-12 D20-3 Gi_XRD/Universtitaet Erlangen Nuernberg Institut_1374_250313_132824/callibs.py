
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
from scipy.signal import find_peaks, savgol_filter
import pandas as pd
import imageio
from silx.io.convert import convert
import os

def auto_contrast(image, lower=1, upper=99):
    vmin = np.percentile(image, lower)
    vmax = np.percentile(image, upper)
    return vmin, vmax

plt.rcParams.update({'font.size': 14})

# Convert .gfrm to HDF5 using silx
input_gfrm = 'Universtitaet Erlangen Nuernberg Institut_1374_250313_132824-000.gfrm'
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
    ax2.set_xlabel('2θ (degrees)')
    ax2.set_ylabel('Intensity')
    ax2.grid(True)
    ax2.legend()
    plt.tight_layout()
    plt.savefig('plot_3_combined.png', dpi=100)
    plt.show()

    # Save azimuthal profiles per peak
    images_for_gif = ['plot_1_original_2D.png', 'plot_2_cake_plot.png', 'plot_3_combined.png']
    peak_profiles = {}

    for idx, peak_pos in enumerate(theta_lin[peaks]):
        peak_mask = (theta_grid[0, :] >= (peak_pos - 0.5)) & (theta_grid[0, :] <= (peak_pos + 0.5))
        azimuthal_profile = np.nanmean(cake_plot[:, peak_mask], axis=1)
        peak_profiles[f'Peak_{idx+1}_{peak_pos:.2f}'] = azimuthal_profile

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
        ax2.set_title(f'Azimuthal Profile for Peak at {peak_pos:.2f}° 2θ')
        ax2.grid(True)

        plt.tight_layout()
        filename = f'plot_peak_{idx+1}.png'
        plt.savefig(filename, dpi=100)
        images_for_gif.append(filename)
        plt.show()

    pd.DataFrame(peak_profiles, index=phi_lin).to_excel(writer, sheet_name='Azimuthal_Profiles')

# Create animated GIF
frames = [imageio.imread(img) for img in images_for_gif]
imageio.mimsave('analysis_summary.gif', frames, duration=2)
