#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Consolidated 2D-XRD Detector Frame Processor
============================================
Author: Tomasz Stawski, tomasz.stawski@bam.de
Date: 2026-06-23
Version: 1.0

Description:
Integrates spatial detector data, performs polar integration, morphological baseline
subtraction, and summarizes specimen orientation metrics.
"""
import os
import sys

# Add parent directories to path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
import caco3_diffraction_pipeline

if __name__ == "__main__":
    print("Running 2D-XRD processing...")
    # The pipeline runs 2D-XRD processing as stage 1
    caco3_diffraction_pipeline.run_data_processing()
    print("2D-XRD processing completed successfully!")
