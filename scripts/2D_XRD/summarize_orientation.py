# -*- coding: utf-8 -*-
"""
2D-XRD Orientation Summary
==========================
Author: Tomasz Stawski, tomasz.stawski@bam.de
Date: 2026-06-19
Version: 0.1

Description:
Compiles and summarises orientation parameters and phase classification metrics from processed 2D-XRD data.
"""
import pandas as pd
import glob
import os

# Find the latest comparison directory
comp_dirs = sorted(glob.glob('data/processed/2D_XRD/output/comparison_*'))
if not comp_dirs:
    print("No comparison directories found.")
    exit()
latest_comp = comp_dirs[-1]
comp_excel = os.path.join(latest_comp, 'comparison_summary.xlsx')
print(f"Reading summary from: {comp_excel}")
summary_df = pd.read_excel(comp_excel, sheet_name='Summary')
print("\n--- Overall Summary ---")
print(summary_df.to_string())

print("\n--- Detailed Orientation Metrics per Sample ---")
# Find the latest individual sample directories
sample_dirs = sorted(glob.glob('data/processed/2D_XRD/output/*_20260428_*'))
sample_dirs = [d for d in sample_dirs if 'comparison' not in d]

for d in sample_dirs:
    data_excel = os.path.join(d, 'output_data.xlsx')
    if os.path.exists(data_excel):
        metrics_df = pd.read_excel(data_excel, sheet_name='Orientation_Metrics')
        # We look for peaks with high DoA or CV (e.g. DoA > 0.1 or CV > 0.05)
        # However, let's just print the peaks that have the highest CV/DoA for each phase
        print(f"\nSample: {os.path.basename(d)}")
        metrics_df = metrics_df.dropna(subset=['Coefficient of Variation (CV)'])
        if not metrics_df.empty:
            metrics_df = metrics_df.sort_values(by='Coefficient of Variation (CV)', ascending=False)
            print(metrics_df[['Peak ID', 'Matched Refs', 'Is C-Axis', 'Degree of Anisotropy (DoA)', 'Coefficient of Variation (CV)']].head(3).to_string(index=False))
        else:
            print("No valid orientation metrics.")

