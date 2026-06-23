#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Consolidated Rocking Curve Processor
====================================
Author: Tomasz Stawski, tomasz.stawski@bam.de
Date: 2026-06-23
Version: 1.0

Description:
Consolidates raw-data extraction, volume correction, and baseline fitting
for CaCO3 rocking curves.
"""
import os
import sys
import argparse

# Add parent directories to path so we can import caco3_diffraction_pipeline
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
import caco3_diffraction_pipeline

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process CaCO3 rocking curves.")
    parser.add_argument("--sample", type=str, choices=["SH-124-B3", "SH-125-A", "SH-125-G", "SH-104-1", "Reference", "all"],
                        default="all", help="Specific sample to process (or 'all').")
    args = parser.parse_args()

    print(f"Running rocking curve processing for: {args.sample}")
    if args.sample == "all":
        caco3_diffraction_pipeline.run_data_processing()
    else:
        # We run the main data processing which handles all samples
        caco3_diffraction_pipeline.run_data_processing()
    print("Processing completed successfully!")
