#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Publication Figure Generator for CaCO3 Thin Films
=================================================
Author: Tomasz Stawski, tomasz.stawski@bam.de
Date: 2026-06-23
Version: 1.0

Description:
Generates all publication-grade figures (Figures 1 to 8 and Appendix Figures A1 to A5)
from processed 2D-XRD and rocking curve datasets.
"""

import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit

# Directories
script_dir = os.path.dirname(os.path.abspath(__file__))
base_dir = os.path.abspath(os.path.join(script_dir, "..", ".."))
PROCESSED_DIR = os.path.join(base_dir, "data/processed")
PLOT_DIR = os.path.join(base_dir, "results/figures")
os.makedirs(PLOT_DIR, exist_ok=True)

def gaussian(t, h, t0, w):
    """Evaluates a Gaussian peak profile."""
    return h * np.exp(-(t - t0)**2 / (2 * w**2))

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

def generate_all_plots():
    """Generates all publication-grade plots and saves them as PNG/SVG files under results/figures/."""
    print("\n======================================================================")
    print("STAGE 2: GENERATING ALL PUBLICATION FIGURES")
    print("======================================================================")
    
    # FIGURE 1: 2D-XRD detector cake plot + 1D profile
    print("Generating Figure 1: 2D-XRD detector cake plot + 1D profile...")
    excel_path = os.path.join(PROCESSED_DIR, "2D_XRD/SH-125-G_output_data.xlsx")
    if os.path.exists(excel_path):
        try:
            df_cake = pd.read_excel(excel_path, sheet_name='Cake_Plot_data', index_col=0)
            df_profile = pd.read_excel(excel_path, sheet_name='Corrected_Profiles')
            theta_lin = df_cake.columns.values.astype(float)
            phi_lin = df_cake.index.values.astype(float)
            cake_plot = df_cake.values
            
            fig = plt.figure(figsize=(10, 8))
            gs = gridspec.GridSpec(2, 1, height_ratios=[1.2, 1.0], hspace=0.25)
            
            ax1 = fig.add_subplot(gs[0])
            vmin = np.percentile(cake_plot[~np.isnan(cake_plot)], 1)
            vmax = np.percentile(cake_plot[~np.isnan(cake_plot)], 99)
            im = ax1.imshow(cake_plot, extent=[theta_lin.min(), theta_lin.max(), phi_lin.min(), phi_lin.max()],
                            aspect='auto', origin='lower', cmap='inferno', vmin=vmin, vmax=vmax)
            ax1.set_ylabel('Azimuthal Angle $\phi$ (°)')
            ax1.set_title('Resampled 2D-XRD detector cake plot (2$\\theta$ vs. $\phi$)', fontweight='bold', fontsize=11)
            fig.colorbar(im, ax=ax1, label='Intensity (counts)', pad=0.02)
            ax1.text(-0.08, 1.05, "(a)", transform=ax1.transAxes, fontsize=14, fontweight='bold', va='top')
            
            ax2 = fig.add_subplot(gs[1])
            ax2.plot(df_profile['2Theta'], df_profile['Smoothed Profile'], color='#3182bd', alpha=0.5, label='Integrated profile')
            ax2.plot(df_profile['2Theta'], df_profile['Corrected Profile'], color='#e6550d', linewidth=1.5, label='Baseline-corrected')
            ax2.plot(df_profile['2Theta'], df_profile['Baseline'], color='grey', linestyle=':', label='Morphological baseline')
            ax2.set_xlabel('2$\\theta$ (°)')
            ax2.set_ylabel('Intensity (counts)')
            ax2.set_title('Azimuthally integrated 1D profile and phase identification', fontweight='bold', fontsize=11)
            ax2.grid(True, linestyle=':', alpha=0.5)
            
            y_max = df_profile['Smoothed Profile'].max()
            ax2.axvline(29.4, color='#2ca02c', linestyle='--', alpha=0.7, label='calcite (104) ref')
            ax2.text(29.4, y_max * 0.9, ' calcite (104)\n 3.03 Å', color='#2ca02c', ha='left', va='top', fontsize=9)
            ax2.axvline(32.8, color='#9467bd', linestyle='--', alpha=0.7, label='vaterite (110) ref')
            ax2.text(32.8, y_max * 0.9, ' vaterite (110)\n 2.73 Å', color='#9467bd', ha='left', va='top', fontsize=9)
            ax2.legend(loc='upper right', framealpha=0.9)
            ax2.set_xlim(20, 55)
            ax2.text(-0.08, 1.05, "(b)", transform=ax2.transAxes, fontsize=14, fontweight='bold', va='top')
            
            plt.savefig(os.path.join(PLOT_DIR, "fig1_2d_xrd_analysis.png"), dpi=300, bbox_inches='tight')
            plt.close()
            print("  Saved Figure 1 to results/figures/")
        except Exception as e:
            print(f"  Error loading/plotting Figure 1: {e}")

    # FIGURE 2: Stacked 2Theta scans for SH-125-G
    print("Generating Figure 2: Stacked 2Theta scans for SH-125-G...")
    g_sym_dir = os.path.join(PROCESSED_DIR, "Symmetric_Scans/SH-125-G")
    if os.path.exists(g_sym_dir):
        try:
            xy_files = sorted(glob.glob(os.path.join(g_sym_dir, "*_2Theta_*_exported.xy")))
            if xy_files:
                fig, ax = plt.subplots(figsize=(8, 9))
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
                ax.set_title('Azimuthal dependence of symmetric 2$\\theta-\\theta$ diffraction scans (SH-125-G)', fontweight='bold')
                ax.set_xlim(27.0, 35.0)
                ax.set_ylim(-10.0, offset + 150.0)
                ax.axvspan(32.4, 33.4, color='#9467bd', alpha=0.08, label='vaterite (110) range')
                ax.axvline(29.4, color='#2ca02c', linestyle=':', alpha=0.5)
                ax.axvline(32.8, color='#9467bd', linestyle=':', alpha=0.5)
                ax.text(29.4, offset + 50.0, 'calcite (104)', color='#2ca02c', ha='center', va='bottom', fontsize=9, fontweight='bold', rotation=90)
                ax.text(32.8, offset + 50.0, 'vaterite (110)', color='#9467bd', ha='center', va='bottom', fontsize=9, fontweight='bold', rotation=90)
                ax.grid(True, which='both', axis='x', linestyle=':', alpha=0.5)
                ax.legend(loc='lower left', framealpha=0.5)
                plt.tight_layout()
                plt.savefig(os.path.join(PLOT_DIR, "fig2_stacked_2theta_sh125g.png"), dpi=300, bbox_inches='tight')
                plt.close()
                print("  Saved Figure 2 to results/figures/")
        except Exception as e:
            print(f"  Error plotting Figure 2: {e}")

    # FIGURE 3: Stacked Rocking Curves (B3 and G)
    print("Generating Figure 3: Stacked Rocking Curves (B3 and G)...")
    samples_config = {
        "SH-124-B3": {
            "processed_dir": os.path.join(PROCESSED_DIR, "Rocking_Curves/SH-124-B3"),
            "phi_values": [0, 30, 60, 90, 120, 150, 180],
            "net_offset": 5000,
            "title": "Sample SH-124-B3 (pure calcite)"
        },
        "SH-125-G": {
            "processed_dir": os.path.join(PROCESSED_DIR, "Rocking_Curves/SH-125-G"),
            "phi_values": [0, 30, 60, 120, 150, 180],
            "net_offset": 5000,
            "title": "Sample SH-125-G (mixed calcite-vaterite)"
        }
    }
    try:
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
        plt.savefig(os.path.join(PLOT_DIR, "fig3_stacked_net_rocking_curves.png"), dpi=300, bbox_inches='tight')
        plt.close()
        print("  Saved Figure 3 to results/figures/")
    except Exception as e:
        print(f"  Error plotting Figure 3: {e}")

    # FIGURE 4: 2D Texture Pole Figures
    print("Generating Figure 4: 2D Texture Pole Figures (2x2 Grid)...")
    samples_pole = [
        ("SH-124-B3", os.path.join(PROCESSED_DIR, "Rocking_Curves/SH-124-B3"), [0, 30, 60, 90, 120, 150, 180], os.path.join(PROCESSED_DIR, "Symmetric_Scans/SH-124-B3")),
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
        plt.close()
        print("  Saved Figure 4 to results/figures/")
    except Exception as e:
        print(f"  Error plotting Figure 4: {e}")

    # FIGURE 5: Phase Areas vs. Phi
    print("Generating Figure 5: Phase Areas vs. Phi...")
    samples_phase = {
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
        plt.close()
        print("  Saved Figure 5 to results/figures/")
    except Exception as e:
        print(f"  Error plotting Figure 5: {e}")

    # FIGURE 6: Raw vs. volume-corrected rocking curve for SH-125-A
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
            
            axes[1].plot(theta, net_intensity, 'purple', label='Net Residual Intensity (Corrected)')
            for idx, row in df_scan.iterrows():
                t0 = row['Peak Center (Theta)']
                h = row['Net Height']
                w = row['FWHM (degrees)'] / 2.355
                if h > 0:
                    axes[1].plot(theta, gaussian(theta, h, t0, w), 'r--', alpha=0.8)
                    axes[1].text(t0 + 0.1, h * 0.9, f"{t0:.2f}°\n$\\chi$={row['Tilt Angle (Chi)']:.2f}°", fontsize=9, color='red')
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

    # FIGURE 7: SH-125-G Side-By-Side Rocking Curves
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

    # FIGURE 8: SH-104-1 Side-By-Side Rocking Curves
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

    # APPENDIX FIGURE A1: Background subtraction sequence
    print("Generating Appendix Figure A1: Rocking curve background subtraction sequence...")
    rc_60_csv = os.path.join(PROCESSED_DIR, "Rocking_Curves/SH-124-B3/SH-124-B3_corrected_rocking_60.csv")
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

    # APPENDIX FIGURE A2: Zoomed net peak deconvolution
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

    # APPENDIX FIGURE A3: Zoomed rocking curve fit at phi=30
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

    # APPENDIX FIGURE A4: Zoomed net peak deconvolution at phi=30
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

    # APPENDIX FIGURE A5: Calcite single crystal rocking curve
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

    print("\nSTAGE 2 COMPLETE: All publication figures generated successfully in results/figures/ directory.")

if __name__ == "__main__":
    generate_all_plots()
