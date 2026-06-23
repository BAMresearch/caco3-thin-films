# -*- coding: utf-8 -*-
"""
Side-by-Side Rocking Curve Plotter
==================================
Author: Tomasz Stawski, tomasz.stawski@bam.de
Date: 2026-06-19
Version: 0.1

Description:
Generates side-by-side comparison plots of rocking curve intensities and peak fitting results.
"""
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
processed_dir = os.path.join(base_dir, "data/processed/Rocking_Curves")
analysis_dir = os.path.join(base_dir, "data/processed/Rocking_Curves/08062026")
artifacts_dir = "/home/tomek/.gemini/antigravity/brain/5bd7f6a4-d046-484c-b15a-5104b17f4fb0/artifacts"

samples_config = {
    "SH-124-B3": {
        "phi_values": [0, 30, 60, 90, 120, 150, 180],
        "net_offset": 6000,
        "raw_factor": 5.0
    },
    "SH-125-A": {
        "phi_values": [0, 30, 60, 90, 120, 150],
        "net_offset": 4000,
        "raw_factor": 5.0
    }
}

for sample, config in samples_config.items():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 7))
    phi_values = config["phi_values"]
    
    # ------------------ Left Subplot: Raw + Baseline (Log) ------------------
    factor = 1.0
    for phi in phi_values:
        csv_path = os.path.join(processed_dir, sample, f"{sample}_corrected_rocking_{phi}.csv")
        if not os.path.exists(csv_path):
            continue
        df = pd.read_csv(csv_path)
        theta = df["Theta (degrees)"]
        raw_intensity = df["Raw Intensity"]
        baseline = df["Model Baseline"]
        
        ax1.plot(theta, raw_intensity * factor, '.', alpha=0.4, label=f"{phi}° ($\times${factor:.1e})")
        ax1.plot(theta, baseline * factor, 'r-', linewidth=1.2)
        factor *= config["raw_factor"]
        
    ax1.set_yscale('log')
    ax1.set_xlabel("Theta (degrees)")
    ax1.set_ylabel("Intensity (counts, stacked log scale)")
    ax1.set_title("Raw Rocking Curves & Baselines (Log Scale)")
    ax1.grid(True, which='both', linestyle=':', alpha=0.5)
    ax1.legend(loc='upper right', fontsize=8, ncol=2)
    
    # ------------------ Right Subplot: Net Corrected (Linear) ------------------
    step_offset = config["net_offset"]
    for idx, phi in enumerate(phi_values):
        csv_path = os.path.join(processed_dir, sample, f"{sample}_corrected_rocking_{phi}.csv")
        if not os.path.exists(csv_path):
            continue
        df = pd.read_csv(csv_path)
        theta = df["Theta (degrees)"]
        net_intensity = df["Corrected Net Intensity"]
        y_val = net_intensity + idx * step_offset
        
        ax2.plot(theta, y_val, '-', linewidth=1.5, label=f"Phi = {phi}°")
        ax2.axhline(y=idx * step_offset, color='grey', linestyle='--', linewidth=0.8, alpha=0.7)
        ax2.text(theta.max() + 0.2, idx * step_offset, f"{phi}°", 
                 verticalalignment='center', fontsize=9, fontweight='bold')
                 
    ax2.set_xlabel("Theta (degrees)")
    ax2.set_ylabel("Net Intensity (counts, stacked linear scale)")
    ax2.set_title("Baseline-Corrected Net Curves (Linear Scale)")
    ax2.grid(True, which='both', linestyle=':', alpha=0.5)
    ax2.set_xlim(theta.min() - 0.5, theta.max() + 1.5)
    
    fig.suptitle(f"{sample} Rocking Curve Analysis: Raw vs. Baseline-Corrected", fontsize=14, fontweight='bold')
    plt.tight_layout()
    
    # Save locally
    out_dir = os.path.join(analysis_dir, sample)
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, f"{sample}_side_by_side.png")
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"Generated side-by-side plot: {out_path}")
    
    # Copy to artifacts
    art_path = os.path.join(artifacts_dir, f"{sample}_side_by_side.png")
    import shutil
    shutil.copy2(out_path, art_path)
    print(f"Copied to artifacts: {art_path}")
