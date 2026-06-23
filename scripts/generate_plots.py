#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CaCO3 Publication Figure Generator
==================================
Generates publication-quality 2D texture pole figures, symmetric scan offsets,
and rocking curves fits.
"""
import os
import sys

# Ensure parent directory is in path so we can import caco3_diffraction_pipeline
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import caco3_diffraction_pipeline

if __name__ == "__main__":
    print("Generating publication figures...")
    caco3_diffraction_pipeline.generate_all_plots()
    print("Figures generated successfully!")
