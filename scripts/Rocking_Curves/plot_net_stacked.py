# -*- coding: utf-8 -*-
"""
Stacked Net Intensity Plotter
=============================
Author: Tomasz Stawski, tomasz.stawski@bam.de
Date: 2026-06-19
Version: 0.1

Description:
Plots baseline-corrected rocking curves stacked vertically to visualise azimuthal variations.
"""
import os
import pandas as pd
import matplotlib.pyplot as plt

base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
processed_dir = os.path.join(base_dir, "data/processed/Rocking_Curves")
analysis_dir = os.path.join(base_dir, "data/processed/Rocking_Curves/08062026")
artifacts_dir = "/home/tomek/.gemini/antigravity/brain/5bd7f6a4-d046-484c-b15a-5104b17f4fb0/artifacts"

samples_config = {
    "SH-124-B3": {
        "phi_values": [0, 30, 60, 90, 120, 150, 180],
        "offset": 6000
    },
    "SH-125-A": {
        "phi_values": [0, 30, 60, 90, 120, 150],
        "offset": 4000
    }
}

for sample, config in samples_config.items():
    plt.figure(figsize=(10, 8))
    phi_values = config["phi_values"]
    step_offset = config["offset"]
    
    for idx, phi in enumerate(phi_values):
        csv_path = os.path.join(processed_dir, sample, f"{sample}_corrected_rocking_{phi}.csv")
        if not os.path.exists(csv_path):
            print(f"Warning: {csv_path} not found.")
            continue
            
        df = pd.read_csv(csv_path)
        theta = df["Theta (degrees)"]
        net_intensity = df["Corrected Net Intensity"]
        y_val = net_intensity + idx * step_offset
        
        # Plot net curve
        plt.plot(theta, y_val, '-', linewidth=1.5, label=f"Phi = {phi}°")
        
        # Plot horizontal line at baseline level (y = idx * step_offset)
        plt.axhline(y=idx * step_offset, color='grey', linestyle='--', linewidth=0.8, alpha=0.7)
        
        # Label each curve on the right
        plt.text(theta.max() + 0.2, idx * step_offset, f"Phi = {phi}°", 
                 verticalalignment='center', fontsize=9, fontweight='bold')
                 
    plt.xlabel("Theta (degrees)")
    plt.ylabel("Baseline-Corrected Net Intensity (counts, stacked linear scale)")
    plt.title(f"Baseline-Corrected Net Rocking Curves: {sample}")
    plt.grid(True, which='both', linestyle=':', alpha=0.5)
    plt.xlim(theta.min() - 0.5, theta.max() + 2.0)
    plt.tight_layout()
    
    # Save locally
    out_dir = os.path.join(analysis_dir, sample)
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, f"{sample}_net_rocking_curves_stacked.png")
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"Generated stacked net rocking curves plot: {out_path}")
    
    # Copy to artifacts
    art_path = os.path.join(artifacts_dir, f"{sample}_net_rocking_curves_stacked.png")
    import shutil
    shutil.copy2(out_path, art_path)
    print(f"Copied to artifacts: {art_path}")
