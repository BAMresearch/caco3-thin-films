#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Publication Package Exporter
============================
Author: Tomasz Stawski, tomasz.stawski@bam.de
Date: 2026-06-23
Version: 1.0

Description:
Packages key processed CSV files, figures, and technical reports into a zip archive for submission.
"""
import os
import shutil
import zipfile

base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
results_dir = os.path.join(base_dir, "results")
zip_path = os.path.join(base_dir, "results_publication_package.zip")

print("Creating publication package...")

# Create a zip archive of the entire results/ folder
if os.path.exists(results_dir):
    with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
        for root, dirs, files in os.walk(results_dir):
            for file in files:
                file_path = os.path.join(root, file)
                arcname = os.path.relpath(file_path, results_dir)
                zipf.write(file_path, arcname)
    print(f"Created publication package zip at: {zip_path}")
else:
    print(f"Error: results directory not found at {results_dir}")
