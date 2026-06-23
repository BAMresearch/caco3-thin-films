# -*- coding: utf-8 -*-
"""
Summary Plot Generator
======================
Author: Tomasz Stawski, tomasz.stawski@bam.de
Date: 2026-06-19
Version: 0.1

Description:
Creates overall summary plots comparing peak areas, widths, and positions across all thin film samples.
"""
import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit

# Set plot style parameters for publication quality (sans-serif, clean borders, legible fonts)
plt.rcParams.update({
    'font.size': 11,
    'font.family': 'sans-serif',
    'font.sans-serif': ['DejaVu Sans', 'Arial', 'Helvetica'],
    'axes.labelsize': 12,
    'axes.titlesize': 13,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'figure.titlesize': 14,
    'svg.fonttype': 'none'  # Keeps text as text in SVG for editing
})

base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
summary_dir = os.path.join(base_dir, "results/reports")
os.makedirs(summary_dir, exist_ok=True)

# Helper for gaussian fit
def gaussian(t, h, t0, w):
    return h * np.exp(-(t-t0)**2 / (2*w**2))

def fit_symmetric_scan(twotheta, intensity):
    # Fit polynomial baseline
    poly_coeff = np.polyfit(twotheta, intensity, 3)
    baseline = np.polyval(poly_coeff, twotheta)
    net_intensity = intensity - baseline
    
    # Fit calcite (104) peak around 29.4
    c_mask = (twotheta >= 28.5) & (twotheta <= 30.5)
    h_c, t0_c, w_c, area_c = 0.0, 29.4, 0.15, 0.0
    try:
        popt_c, _ = curve_fit(gaussian, twotheta[c_mask], net_intensity[c_mask], p0=[intensity.max() - baseline.max(), 29.4, 0.15])
        h_c, t0_c, w_c = popt_c
        area_c = h_c * w_c * np.sqrt(2 * np.pi)
    except:
        pass
        
    # Fit vaterite (110) peak around 32.8
    v_mask = (twotheta >= 31.8) & (twotheta <= 33.8)
    h_v, t0_v, w_v, area_v = 0.0, 32.8, 0.15, 0.0
    try:
        popt_v, _ = curve_fit(gaussian, twotheta[v_mask], net_intensity[v_mask], p0=[1000, 32.8, 0.15])
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
# FIGURE 1: 2D-XRD CAKE PLOT AND INTEGRATED PROFILE (THE STIMULUS)
# ==============================================================================
print("Generating Figure 1: 2D-XRD Analysis (Cake plot + Integrated profile)...")
# Find SH-125-G Excel output file
excel_files = glob.glob(os.path.join(base_dir, "data/processed/2D_XRD/output/SH-125-G_OK_Set1_*/output_data.xlsx"))
if excel_files:
    excel_path = excel_files[0]
    try:
        # Load sheets
        df_cake = pd.read_excel(excel_path, sheet_name='Cake_Plot_data', index_col=0)
        df_profile = pd.read_excel(excel_path, sheet_name='Corrected_Profiles')
        
        theta_lin = df_cake.columns.values.astype(float)
        phi_lin = df_cake.index.values.astype(float)
        cake_plot = df_cake.values
        
        fig = plt.figure(figsize=(10, 8))
        gs = gridspec.GridSpec(2, 1, height_ratios=[1.2, 1.0], hspace=0.25)
        
        # Subplot A: Cake Plot
        ax1 = fig.add_subplot(gs[0])
        vmin = np.percentile(cake_plot[~np.isnan(cake_plot)], 1)
        vmax = np.percentile(cake_plot[~np.isnan(cake_plot)], 99)
        im = ax1.imshow(cake_plot, extent=[theta_lin.min(), theta_lin.max(), phi_lin.min(), phi_lin.max()],
                        aspect='auto', origin='lower', cmap='inferno', vmin=vmin, vmax=vmax)
        ax1.set_ylabel('Azimuthal Angle $\phi$ (°)')
        ax1.set_title('Resampled 2D-XRD detector cake plot (2$\\theta$ vs. $\phi$)', fontweight='bold', fontsize=11)
        fig.colorbar(im, ax=ax1, label='Intensity (counts)', pad=0.02)
        ax1.text(-0.08, 1.05, "(a)", transform=ax1.transAxes, fontsize=14, fontweight='bold', va='top')
        
        # Subplot B: Integrated Profile
        ax2 = fig.add_subplot(gs[1])
        ax2.plot(df_profile['2Theta'], df_profile['Smoothed Profile'], color='#3182bd', alpha=0.5, label='Integrated profile')
        ax2.plot(df_profile['2Theta'], df_profile['Corrected Profile'], color='#e6550d', linewidth=1.5, label='Baseline-corrected')
        ax2.plot(df_profile['2Theta'], df_profile['Baseline'], color='grey', linestyle=':', label='Morphological baseline')
        
        # Add labels and vertical lines
        ax2.set_xlabel('2$\\theta$ (°)')
        ax2.set_ylabel('Intensity (counts)')
        ax2.set_title('Azimuthally integrated 1D profile and phase identification', fontweight='bold', fontsize=11)
        ax2.grid(True, linestyle=':', alpha=0.5)
        
        # References for calcite and vaterite
        y_max = df_profile['Smoothed Profile'].max()
        ax2.axvline(29.4, color='#2ca02c', linestyle='--', alpha=0.7, label='calcite (104) ref')
        ax2.text(29.4, y_max * 0.9, ' calcite (104)\n 3.03 Å', color='#2ca02c', ha='left', va='top', fontsize=9)
        
        ax2.axvline(32.8, color='#9467bd', linestyle='--', alpha=0.7, label='vaterite (110) ref')
        ax2.text(32.8, y_max * 0.9, ' vaterite (110)\n 2.73 Å', color='#9467bd', ha='left', va='top', fontsize=9)
        
        ax2.legend(loc='upper right', framealpha=0.9)
        ax2.set_xlim(20, 55)
        ax2.text(-0.08, 1.05, "(b)", transform=ax2.transAxes, fontsize=14, fontweight='bold', va='top')
        
        # Save Figure 1
        plt.savefig(os.path.join(summary_dir, "fig1_2d_xrd_analysis.png"), dpi=300, bbox_inches='tight')
        plt.savefig(os.path.join(summary_dir, "fig1_2d_xrd_analysis.svg"), bbox_inches='tight')
        plt.close()
        print("  Successfully generated fig1_2d_xrd_analysis.png and .svg")
    except Exception as e:
        print(f"  Error loading/plotting Figure 1: {e}")
else:
    print("  Could not find SH-125-G Excel file for Figure 1")

# ==============================================================================
# FIGURE 2: STACKED 2THETA SCANS FOR SH-125-G SHOWING VATERITE CONFINEMENT
# ==============================================================================
print("Generating Figure 2: Stacked 2Theta scans for SH-125-G...")
processed_125g_dir = os.path.join(base_dir, "data/processed/Rocking_Curves/SH-125-G")
if os.path.exists(processed_125g_dir):
    xy_files = sorted(glob.glob(os.path.join(processed_125g_dir, "*_2Theta_*_exported.xy")))
    if xy_files:
        fig, ax = plt.subplots(figsize=(8, 9))
        offset = 0.0
        # Sort files by phi angle
        phi_to_file = {}
        for f in xy_files:
            phi_val = int(os.path.basename(f).split("2Theta_")[1].split("_")[0])
            phi_to_file[phi_val] = f
            
        sorted_phis = sorted(phi_to_file.keys())
        colors = plt.cm.viridis(np.linspace(0, 0.85, len(sorted_phis)))
        
        for idx, phi in enumerate(sorted_phis):
            fpath = phi_to_file[phi]
            arr = np.loadtxt(fpath, skiprows=1)
            twotheta = arr[:, 0]
            intensity_k = arr[:, 1] / 1e3  # Convert to kcounts to correct the vertical scale
            
            ax.plot(twotheta, intensity_k + offset, color=colors[idx], linewidth=1.5, label=f"$\phi$ = {phi}°")
            ax.axhline(y=offset, color='grey', linestyle='--', linewidth=0.5, alpha=0.5)
            # Label each scan
            ax.text(twotheta.max() + 0.1, offset + 10.0, f"{phi}°", fontweight='bold', va='center', color=colors[idx])
            offset += 60.0  # Stack offset (kcounts)
            
        ax.set_xlabel('2$\\theta$ (°)')
        ax.set_ylabel('Intensity (kcounts, stacked)')
        ax.set_title('Azimuthal dependence of symmetric 2$\\theta-\\theta$ diffraction scans (SH-125-G)', fontweight='bold')
        ax.set_xlim(27.0, 35.0)
        ax.set_ylim(-10.0, offset + 150.0)  # Correct limit based on stacked intensity range
        
        # Highlight vaterite (110) region
        ax.axvspan(32.4, 33.4, color='#9467bd', alpha=0.08, label='vaterite (110) range')
        ax.axvline(29.4, color='#2ca02c', linestyle=':', alpha=0.5)
        ax.axvline(32.8, color='#9467bd', linestyle=':', alpha=0.5)
        
        # Add peak index labels
        ax.text(29.4, offset + 50.0, 'calcite (104)', color='#2ca02c', ha='center', va='bottom', fontsize=9, fontweight='bold', rotation=90)
        ax.text(32.8, offset + 50.0, 'vaterite (110)', color='#9467bd', ha='center', va='bottom', fontsize=9, fontweight='bold', rotation=90)
        
        ax.grid(True, which='both', axis='x', linestyle=':', alpha=0.5)
        ax.legend(loc='lower left', framealpha=0.5)
        plt.tight_layout()
        plt.savefig(os.path.join(summary_dir, "fig2_stacked_2theta_sh125g.png"), dpi=300, bbox_inches='tight')
        plt.savefig(os.path.join(summary_dir, "fig2_stacked_2theta_sh125g.svg"), bbox_inches='tight')
        plt.close()
        print("  Successfully generated fig2_stacked_2theta_sh125g.png and .svg")
else:
    print("  Could not find SH-125-G processed directory for Figure 2")

# ==============================================================================
# FIGURE 3: STACKED ROCKING CURVES (LINEAR SCALE) FOR B3 AND G
# ==============================================================================
print("Generating Figure 3: Stacked Rocking Curves (B3 and G side-by-side)...")
samples_config = {
    "SH-124-B3": {
        "processed_dir": os.path.join(base_dir, "data/processed/Rocking_Curves/SH-124-B3"),
        "phi_values": [0, 30, 60, 90, 120, 150, 180],
        "net_offset": 5000,
        "title": "Sample SH-124-B3 (pure calcite)"
    },
    "SH-125-G": {
        "processed_dir": os.path.join(base_dir, "data/processed/Rocking_Curves/SH-125-G"),
        "phi_values": [0, 30, 60, 120, 150, 180], # Skipped missing 90
        "net_offset": 5000,
        "title": "Sample SH-125-G (mixed calcite-vaterite)"
    }
}

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 7), sharex=True)
axes = [ax1, ax2]

for idx_s, (sample_name, config) in enumerate(samples_config.items()):
    ax = axes[idx_s]
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
    ax.set_ylabel("Net Intensity (counts, stacked linear scale)")
    ax.set_title(config["title"], fontweight='bold', fontsize=12)
    ax.grid(True, which='both', linestyle=':', alpha=0.5)
    ax.set_xlim(4.0, 26.0)
    ax.text(-0.08, 1.05, f"({chr(97 + idx_s)})", transform=ax.transAxes, fontsize=14, fontweight='bold', va='top')

plt.suptitle("Baseline-corrected net rocking curves vs. azimuthal angle $\phi$", fontsize=14, fontweight='bold', y=0.98)
plt.tight_layout()
plt.savefig(os.path.join(summary_dir, "fig3_stacked_net_rocking_curves.png"), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(summary_dir, "fig3_stacked_net_rocking_curves.svg"), bbox_inches='tight')
plt.close()
print("  Successfully generated fig3_stacked_net_rocking_curves.png and .svg")

# ==============================================================================
# FIGURE 4: 2D TEXTURE POLE FIGURES (2x2 GRID)
# ==============================================================================
print("Generating Figure 4: 2D Texture Pole Figures (2x2 Grid)...")
samples_pole = [
    ("SH-124-B3", os.path.join(base_dir, "data/processed/Rocking_Curves/SH-124-B3"), [0, 30, 60, 90, 120, 150, 180]),
    ("SH-125-A", os.path.join(base_dir, "data/processed/Rocking_Curves/SH-125-A"), [0, 30, 60, 90, 120, 150]),
    ("SH-104-1", os.path.join(base_dir, "data/processed/Rocking_Curves/SH-104-1"), [0, 30, 60, 90, 120, 150]),
    ("SH-125-G", os.path.join(base_dir, "data/processed/Rocking_Curves/SH-125-G"), [0, 30, 60, 120, 150, 180]) # Omit 90
]

fig = plt.figure(figsize=(12, 11))
gs = gridspec.GridSpec(2, 2, wspace=0.3, hspace=0.3)

for idx_s, (sample, pdir, phis) in enumerate(samples_pole):
    # 1. Load rocking curves and interpolate
    rocking_scans = {}
    theta_0_dict = {}
    
    for phi in phis:
        # Check if 2Theta scan exists in pdir
        xy_2t_path = glob.glob(os.path.join(pdir, f"*_2Theta_{phi}_exported.xy"))
        theta_0 = 14.7
        if xy_2t_path:
            arr = np.loadtxt(xy_2t_path[0], skiprows=1)
            twotheta_arr = arr[:, 0]
            intensity_arr = arr[:, 1]
            poly_coeff = np.polyfit(twotheta_arr, intensity_arr, 3)
            net_intensity = intensity_arr - np.polyval(poly_coeff, twotheta_arr)
            c_mask = (twotheta_arr >= 28.5) & (twotheta_arr <= 30.5)
            try:
                popt_c, _ = curve_fit(gaussian, twotheta_arr[c_mask], net_intensity[c_mask], p0=[intensity_arr.max(), 29.4, 0.15])
                theta_0 = popt_c[1] / 2.0
            except:
                pass
        theta_0_dict[phi] = theta_0
        
        # Load rocking curve
        rc_path = os.path.join(pdir, f"{sample}_corrected_rocking_{phi}.csv")
        if os.path.exists(rc_path):
            df_rc = pd.read_csv(rc_path)
            rocking_scans[phi] = interp1d(df_rc['Theta (degrees)'].values, df_rc['Corrected Net Intensity'].values, bounds_error=False, fill_value=0.0)
            
    # 2. Build polar grid
    phi_polar_deg = np.arange(0, 360 + 30, 30)
    phi_polar_rad = np.radians(phi_polar_deg)
    r_tilt = np.linspace(0, 10.0, 300)
    Z_net = np.zeros((len(r_tilt), len(phi_polar_deg)))
    
    # Track which angles are measured (for missing scans like G 90)
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
    
    # For missing data, set to NaN to display as white gap in contourf
    Z_net_masked = Z_net.copy()
    Z_net_masked[~measured_polar_mask] = np.nan
    
    # Plot in sub-panel
    ax = fig.add_subplot(gs[idx_s], projection='polar')
    Phi_mesh, R_mesh = np.meshgrid(phi_polar_rad, r_tilt)
    
    # Handle levels based on data range
    z_max = np.nanmax(Z_net_masked) if not np.all(np.isnan(Z_net_masked)) else 1000.0
    levels = np.linspace(0, max(z_max, 10.0), 60)
    
    contour = ax.contourf(Phi_mesh, R_mesh, Z_net_masked, levels=levels, cmap='plasma')
    
    # Labels
    ax.set_theta_zero_location("E")  # 0° on the right
    ax.set_theta_direction(1)       # Counter-clockwise
    ax.set_ylim(0, 10.0)
    ax.set_yticks([2, 4, 6, 8, 10])
    ax.set_yticklabels(["2°", "4°", "6°", "8°", "10°"], color="#555555", size=9)
    ax.set_rlabel_position(45)
    
    ax.set_xticks(np.radians(np.arange(0, 360, 30)))
    ax.set_xticklabels([f"{d}°" for d in np.arange(0, 360, 30)], size=9)
    ax.grid(True, color="grey", linestyle=":", alpha=0.5)
    
    ax.set_title(f"({chr(97 + idx_s)}) {sample}", y=1.08, fontweight='bold', fontsize=12)
    
    # Add colorbar for each panel
    cbar = fig.colorbar(contour, ax=ax, pad=0.08, shrink=0.7)
    cbar.ax.tick_params(labelsize=9)
    cbar.set_label("Net Intensity (counts)", fontsize=9)

plt.suptitle("Compiled calcite (104) 2D polar texture pole figures", fontsize=15, fontweight='bold', y=0.98)
plt.savefig(os.path.join(summary_dir, "fig4_texture_pole_figures.png"), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(summary_dir, "fig4_texture_pole_figures.svg"), bbox_inches='tight')
plt.close()
print("  Successfully generated fig4_texture_pole_figures.png and .svg")

# ==============================================================================
# FIGURE 5: PHASE METRICS (CALCITE VS VATERITE PEAK AREAS) VS PHI
# ==============================================================================
print("Generating Figure 5: Phase Areas vs. Phi (vaterite confinement)...")
# We will read the .xy files for SH-125-A, SH-104-1, and SH-125-G, fit the calcite and vaterite peaks,
# and plot them as a function of phi to show the confinement of vaterite in-plane.

samples_phase = {
    "SH-125-A": {
        "dir": os.path.join(base_dir, "data/processed/Rocking_Curves/SH-125-A"),
        "phis": [0, 30, 60, 90, 120, 150],
        "color_c": "#2ca02c", "color_v": "#9467bd", "marker": "o", "ls": "-"
    },
    "SH-104-1": {
        "dir": os.path.join(base_dir, "data/processed/Rocking_Curves/SH-104-1"),
        "phis": [0, 30, 60, 90, 120, 150],
        "color_c": "#1f77b4", "color_v": "#d62728", "marker": "s", "ls": "--"
    },
    "SH-125-G": {
        "dir": os.path.join(base_dir, "data/processed/Rocking_Curves/SH-125-G"),
        "phis": [0, 30, 60, 90, 120, 150, 180],
        "color_c": "#ff7f0e", "color_v": "#8c564b", "marker": "^", "ls": "-."
    }
}

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
            
    # Plot calcite in top panel
    ax1.plot(valid_phis, calcite_areas, marker=config["marker"], linestyle=config["ls"],
             color=config["color_c"], linewidth=2, label=f"{sample} calcite (104)")
             
    # Plot vaterite in bottom panel
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

# Highlight the selective vaterite orientation zone
ax2.axvspan(25, 65, color='#9467bd', alpha=0.1)
ax2.text(45, ax2.get_ylim()[1]*0.8, "Epitaxial vaterite\norientation zone", color='#5c3d75', ha='center', fontweight='bold', fontsize=9)

plt.tight_layout()
plt.savefig(os.path.join(summary_dir, "fig5_phase_metrics_vs_phi.png"), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(summary_dir, "fig5_phase_metrics_vs_phi.svg"), bbox_inches='tight')
plt.close()
print("  Successfully generated fig5_phase_metrics_vs_phi.png and .svg")

# ==============================================================================
# APPENDIX FIGURES: BACKGROUND SUBTRACTION AND PEAK DECONVOLUTION
# ==============================================================================
print("Generating Appendix Figure A1: Rocking curve background subtraction...")
rc_60_csv = os.path.join(base_dir, "data/processed/Rocking_Curves/SH-124-B3/SH-124-B3_corrected_rocking_60.csv")
metrics_csv = os.path.join(base_dir, "data/processed/Rocking_Curves/Reference/all_samples_rocking_peaks_vs_phi.csv")

if os.path.exists(rc_60_csv) and os.path.exists(metrics_csv):
    try:
        # Load raw data and baseline
        df_rc = pd.read_csv(rc_60_csv)
        theta = df_rc['Theta (degrees)'].values
        raw_int = df_rc['Raw Intensity'].values
        baseline = df_rc['Model Baseline'].values
        net_raw = raw_int - baseline
        
        # Load metrics to reconstruct fits
        df_m = pd.read_csv(metrics_csv)
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
                    "name": name,
                    "center": center,
                    "height": height,
                    "w": w,
                    "tilt": tilt,
                    "curve": peak_y
                })
        
        # Plot Figure A1: Before, During, and After background subtraction
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 5.5))
        fig.suptitle("Rocking curve background subtraction sequence (SH-124-B3, $\phi = 60^\circ$)", fontsize=14, fontweight='bold', y=0.98)
        
        # Panel 1: BEFORE (Raw Experimental Data)
        ax1.plot(theta, raw_int, 'o', color='gray', markersize=3, alpha=0.5, label='Raw data')
        ax1.set_xlabel("Theta $\\theta$ (°)")
        ax1.set_ylabel("Intensity (counts)")
        ax1.set_title("Before: Raw Intensity Profile")
        ax1.grid(True, linestyle=':', alpha=0.5)
        ax1.legend(loc='upper right')
        ax1.set_ylim(0, np.max(raw_int)*1.1)
        ax1.text(0.05, 0.95, "Before", transform=ax1.transAxes, fontsize=12, fontweight='bold', va='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        # Panel 2: DURING (Fitting Background & Peak Envelope)
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
        
        # Panel 3: AFTER (Subtracted Background showing Net Gaussian Components)
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
        plt.savefig(os.path.join(summary_dir, "fig_a1_background_subtraction.png"), dpi=300, bbox_inches='tight')
        plt.savefig(os.path.join(summary_dir, "fig_a1_background_subtraction.svg"), bbox_inches='tight')
        plt.close()
        print("  Successfully generated fig_a1_background_subtraction.png and .svg")
        
        # Plot Figure A2: Zoomed Net Peak Deconvolution
        print("Generating Appendix Figure A2: Zoomed net peak deconvolution...")
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
        plt.savefig(os.path.join(summary_dir, "fig_a2_peak_deconvolution.png"), dpi=300, bbox_inches='tight')
        plt.savefig(os.path.join(summary_dir, "fig_a2_peak_deconvolution.svg"), bbox_inches='tight')
        plt.close()
        print("  Successfully generated fig_a2_peak_deconvolution.png and .svg")
        
    except Exception as e:
        print(f"  Error generating Appendix figures: {e}")
else:
    print("  Required files for Appendix figures not found")

print("\nAll publication figures generated successfully in: ", summary_dir)
