#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Consolidated Rocking Curve Plotter
==================================
Author: Tomasz Stawski, tomasz.stawski@bam.de
Date: 2026-06-23
Version: 1.0

Description:
Consolidates plotting functionalities for net intensity curves, stacked curves,
zoom deconvolution fits, and polar 2D texture pole figures.
"""
import os
import sys
import argparse

# Add parent directories to path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
import caco3_diffraction_pipeline

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate plots for CaCO3 diffraction analysis.")
    parser.add_argument("--plot-type", type=str, choices=["all", "stacked", "texture", "fits", "side-by-side"],
                        default="all", help="Type of plots to generate.")
    args = parser.parse_args()

    print(f"Generating plots: {args.plot_type}")
    caco3_diffraction_pipeline.generate_all_plots()
    print("Plot generation completed successfully!")
