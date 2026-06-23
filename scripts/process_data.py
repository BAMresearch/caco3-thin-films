#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CaCO3 Diffraction Data Processor
================================
Runs the data extraction, volume correction, baseline subtraction,
and multi-peak fitting algorithms.
"""
import os
import sys

# Ensure parent directory is in path so we can import caco3_diffraction_pipeline
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import caco3_diffraction_pipeline

if __name__ == "__main__":
    print("Running CaCO3 diffraction data processing...")
    caco3_diffraction_pipeline.run_data_processing()
    print("Data processing complete!")
