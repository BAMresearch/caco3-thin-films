# -*- coding: utf-8 -*-
"""
Master Analysis Coordinator
===========================
Author: Tomasz Stawski, tomasz.stawski@bam.de
Date: 2026-06-19
Version: 0.1

Description:
Master script that coordinates the entire data processing and plot generation workflow.
"""
"""
Master script to execute the complete XRD analysis route.
Retraces the entire processing from raw data to publication-ready plots.
"""

import sys
import os

# Ensure the package directory is in the import path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import caco3_diffraction_pipeline

if __name__ == "__main__":
    print("Executing CaCO3 XRD analysis pipeline...")
    caco3_diffraction_pipeline.run_data_processing()
    caco3_diffraction_pipeline.generate_all_plots()
    print("All tasks completed successfully!")
