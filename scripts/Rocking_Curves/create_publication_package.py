# -*- coding: utf-8 -*-
"""
Publication Package Exporter
============================
Author: Tomasz Stawski, tomasz.stawski@bam.de
Date: 2026-06-19
Version: 0.1

Description:
Packages key vector figures and data tables into a zip archive for submission.
"""
import os
import shutil
import glob
import zipfile

base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
package_dir = os.path.join(base_dir, "results")
data_dir = os.path.join(package_dir, "data")

os.makedirs(data_dir, exist_ok=True)

# 1. Create subfolders for data
os.makedirs(os.path.join(data_dir, "2D_XRD"), exist_ok=True)
os.makedirs(os.path.join(data_dir, "Symmetric_Scans/SH-125-G"), exist_ok=True)
os.makedirs(os.path.join(data_dir, "Symmetric_Scans/SH-125-A"), exist_ok=True)
os.makedirs(os.path.join(data_dir, "Symmetric_Scans/SH-104-1"), exist_ok=True)

os.makedirs(os.path.join(data_dir, "Rocking_Curves/SH-124-B3"), exist_ok=True)
os.makedirs(os.path.join(data_dir, "Rocking_Curves/SH-125-G"), exist_ok=True)
os.makedirs(os.path.join(data_dir, "Rocking_Curves/SH-125-A"), exist_ok=True)
os.makedirs(os.path.join(data_dir, "Rocking_Curves/SH-104-1"), exist_ok=True)
os.makedirs(os.path.join(data_dir, "Rocking_Curves/Reference"), exist_ok=True)

# 2. Copy the files
# Figure 1: 2D-XRD Excel sheet
excel_files = glob.glob(os.path.join(base_dir, "data/processed/2D_XRD/output/SH-125-G_OK_Set1_*/output_data.xlsx"))
if excel_files:
    shutil.copy(excel_files[0], os.path.join(data_dir, "2D_XRD/SH-125-G_output_data.xlsx"))

# Figure 2: Symmetric scans for SH-125-G
sh_125g_2t_files = glob.glob(os.path.join(base_dir, "data/processed/Rocking_Curves/SH-125-G/*_2Theta_*_exported.xy"))
for f in sh_125g_2t_files:
    shutil.copy(f, os.path.join(data_dir, "Symmetric_Scans/SH-125-G/"))

# Figure 5: Symmetric scans for SH-125-A and SH-104-1
sh_125a_2t_files = glob.glob(os.path.join(base_dir, "data/processed/Rocking_Curves/SH-125-A/*_2Theta_*_exported.xy"))
for f in sh_125a_2t_files:
    shutil.copy(f, os.path.join(data_dir, "Symmetric_Scans/SH-125-A/"))

sh_1041_2t_files = glob.glob(os.path.join(base_dir, "data/processed/Rocking_Curves/SH-104-1/*_2Theta_*_exported.xy"))
for f in sh_1041_2t_files:
    shutil.copy(f, os.path.join(data_dir, "Symmetric_Scans/SH-104-1/"))

# Rocking curves CSVs
for s, sdir in [("SH-124-B3", "08062026"), ("SH-125-A", "08062026"), ("SH-104-1", "15062026"), ("SH-125-G", "16062026")]:
    rc_files = glob.glob(os.path.join(base_dir, f"data/{sdir}/Processed/{s}/*_corrected_rocking_*.csv"))
    for f in rc_files:
        shutil.copy(f, os.path.join(data_dir, f"Rocking_Curves/{s}/"))

# Single crystal metrics reference CSV and general rocking curve metrics
sc_csv_path = os.path.join(base_dir, "data/processed/Rocking_Curves/Reference/calcite_single_crystal_rocking_peaks_metrics.csv")
if os.path.exists(sc_csv_path):
    shutil.copy(sc_csv_path, os.path.join(data_dir, "Rocking_Curves/Reference/calcite_single_crystal_rocking_peaks_metrics.csv"))

all_peaks_csv = os.path.join(base_dir, "data/processed/Rocking_Curves/Reference/all_samples_rocking_peaks_vs_phi.csv")
if os.path.exists(all_peaks_csv):
    shutil.copy(all_peaks_csv, os.path.join(data_dir, "Rocking_Curves/Reference/all_samples_rocking_peaks_vs_phi.csv"))

print("Copied all data files successfully.")

# 3. Create the ZIP package
zip_path = os.path.join(base_dir, "results.zip")
with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
    for root, dirs, files in os.walk(package_dir):
        for file in files:
            file_path = os.path.join(root, file)
            # Add file with relative path inside the zip
            arcname = os.path.relpath(file_path, package_dir)
            zipf.write(file_path, arcname)

print(f"Created self-contained ZIP package at {zip_path}")
