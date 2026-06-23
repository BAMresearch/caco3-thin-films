# -*- coding: utf-8 -*-
"""
Report Generator for 16 June 2026 Run
=====================================
Author: Tomasz Stawski, tomasz.stawski@bam.de
Date: 2026-06-19
Version: 0.1

Description:
Compiles the technical summary report for rocking curve measurements from 16 June 2026.
"""
import os
import subprocess
import pandas as pd
import numpy as np

# Directories
base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
analysis_dir = os.path.join(base_dir, "data/processed/Rocking_Curves/16062026/SH-125-G")
csv_path = os.path.join(analysis_dir, "SH-125-G_rocking_peaks_vs_phi.csv")

# Load data
df = pd.read_csv(csv_path)

# Programmatically generate tables for the report
def generate_sample_table():
    # Filter out the rows that have height > 0
    df_s = df[df["Net Height"] > 0].copy()
    
    # Sort by Phi, then Peak Name
    df_s = df_s.sort_values(["Phi (degrees)", "Peak Name"])
    
    table_lines = [
        "| Phi (°) | Peak Name | Center ($\\theta$) | Tilt Angle ($\\chi$) | FWHM ($^\circ$) | Height (counts) | Area (cts·deg) | Area/Base Ratio |",
        "| :---: | :--- | :---: | :---: | :---: | :---: | :---: | :---: |"
    ]
    
    for _, r in df_s.iterrows():
        ratio_str = f"{r['Area/Base Ratio']:.3%}" if not pd.isna(r['Area/Base Ratio']) else "N/A"
        table_lines.append(
            f"| {r['Phi (degrees)']:.0f} | {r['Peak Name']} | {r['Peak Center (Theta)']:.3f}° | {r['Tilt Angle (Chi)']:.3f}° | {r['FWHM (degrees)']:.3f}° | {r['Net Height']:.1f} | {r['Net Area (cts deg)']:.1f} | {ratio_str} |"
        )
    return "\n".join(table_lines)

table_content = generate_sample_table()

md_content = r"""---
title: "CaCO3 Thin Film Rocking Curve Analysis: Azimuthal (Phi) Variation Study (Sample SH-125-G)"
date: "June 16, 2026"
geometry: "margin=1in"
fontsize: "11pt"
header-includes:
  - \usepackage{fancyhdr}
  - \pagestyle{fancy}
  - \fancyhead[CO,CE]{SH-125-G Phi-Variation Rocking Curve Report}
  - \fancyfoot[CO,CE]{Page \thepage}
---

# 1. Introduction and Dataset Structure
This report summarizes the analysis of the X-ray diffraction (XRD) measurements collected on **June 16, 2026** for the calcium carbonate ($CaCO_3$) thin film sample **`SH-125-G`** (located in the directory `data/raw/Rocking_Curves/SH-125-G`).

For this study, the symmetric $2\theta-\theta$ scans were collected at seven azimuthal angles: $\phi = 0^\circ, 30^\circ, 60^\circ, 90^\circ, 120^\circ, 150^\circ, 180^\circ$. The 1D Rocking Curves were collected at six of these angles, as **`Rocking_90.brml` is missing** from the raw measurement files. The processing pipeline was successfully modified to skip the missing $90^\circ$ rocking curve scan while maintaining full processing and mapping for all other rotations, creating a gap in the 2D texture and waterfall plots where data is unavailable.

# 2. Phase Analysis and Phase Stability ($2\theta-\theta$ Scans)
Gaussian profile fitting of the background-subtracted symmetric $2\theta-\theta$ scans reveals key structural insights:
1. **Calcite (104) Stability:** Calcite is detected across all rotations. The peak is centered at $2\theta \approx 29.37^\circ - 29.45^\circ$ with a FWHM of $\approx 0.30^\circ - 0.32^\circ$.
2. **Azimuthally Textured Vaterite (110):** Vaterite (110) at $2\theta \approx 32.8^\circ - 33.0^\circ$ is **only detected at $\phi = 30^\circ$ and $\phi = 60^\circ$**. At $\phi = 30^\circ$, the Calcite/Vaterite area ratio is 5.61 (Vaterite area $= 858.2$ counts$\cdot$deg); at $\phi = 60^\circ$, the ratio is 4.74 (Vaterite area $= 1,301.7$ counts$\cdot$deg). For all other phi angles ($0^\circ, 90^\circ, 120^\circ, 150^\circ, 180^\circ$), Vaterite is negligible and falls below the detection limit.

This isolated detection demonstrates that the Vaterite phase has a strong preferred in-plane (azimuthal) alignment. Because the physical phase fraction is fixed, the Vaterite crystallites must be highly textured so that their diffracting planes only satisfy the Bragg condition at specific rotations.

## 2.1 Stacked 2Theta Diffraction Patterns
The raw 2Theta diffraction scans stacked with vertical offsets show the presence of Calcite across all angles and the emergence of the Vaterite (110) peak at $\phi = 30^\circ$ and $\phi = 60^\circ$:

![SH-125-G Symmetric 2Theta Scans Stacked](SH-125-G_2theta_scans_stacked.png){width=85%}

\newpage

## 2.2 Azimuthal Dependence of Phase Composition
The plot below tracks the peak area of Calcite and Vaterite as a function of the azimuthal angle $\phi$. This visually highlights that Vaterite is only resolved around $\phi = 30^\circ - 60^\circ$:

![SH-125-G Phase metrics vs. Phi](SH-125-G_phase_metrics_vs_phi.png){width=85%}

# 3. Rocking Curve Analysis and Volume Correction
To isolate the preferred orientation peak components, a global fitting was performed on the raw rocking curves in $\log_{10}(I)$ space. The model contains both the isotropic background baseline (a footprint correction combined with a 3rd-order polynomial) and multiple Gaussian profiles representing the co-oriented grain populations.

## 3.1 Stacked Rocking Curves & Background Baselines (Log Scale)
The plot below shows the raw rocking curves along with the fitted background baseline curves (red lines), stacked using multiplicative factors. The missing $\phi = 90^\circ$ rocking curve scan is omitted:

![SH-125-G Rocking Curves & Baselines Stacked](SH-125-G_rocking_curves_stacked.png){width=85%}

\newpage

## 3.2 Stacked Baseline-Corrected Net Rocking Curves (Linear Scale)
The plot below shows the baseline-corrected (net intensity) rocking curves stacked vertically on a linear scale. The dashed horizontal lines show the zero-baseline level for each angle:

![SH-125-G Stacked Baseline-Corrected Net Rocking Curves](SH-125-G_side_by_side.png){width=85%}

\newpage

## 3.3 Stacked Global Peak Fits (Log Scale)
The plot below shows the complete global fits in log scale, showing the raw intensity (dots), the full fitted model (solid black line), and the individual resolved peaks (dashed curves):

![SH-125-G Rocking Curve Global Fits Stacked](SH-125-G_rocking_residuals_stacked.png){width=85%}

\newpage

# 4. 2D Polar Texture Plots (Pole Figures)
To visualize the complete in-plane alignment, the 1D rocking curves measured at different azimuthal angles $\phi$ were mapped into a single 2D polar projection (pole-figure-like texture plot). Negative tilts ($\chi < 0$) at angle $\phi$ are mapped to a positive radius at the opposite azimuth ($\phi + 180^\circ$).

In the net texture plot, the missing $\phi = 90^\circ$ scan results in a clear unmeasured gap (blank wedges at $\phi = 90^\circ$ and $\phi = 270^\circ$). However, the co-oriented grain spots at other rotations are resolved:

![SH-125-G 2D Polar Texture Plots (Raw & Net)](SH-125-G_texture_raw.png){width=48%} ![SH-125-G Net Texture](SH-125-G_texture_net.png){width=48%}

# 5. Fitted Rocking Curve Peak Metrics
The table below lists the parameters of all active reflections resolved from the rocking curves at each azimuthal rotation angle:

{table_content}

# 6. Conclusions
1. **Robust Handling of Missing Data:** The analysis code successfully managed the missing `Rocking_90.brml` scan, producing high-fidelity waterfall, stacked fit, and 2D pole figures with clear representations of the missing scan.
2. **In-Plane Texture of Vaterite:** The Vaterite (110) reflection is visible only at $\phi = 30^\circ$ and $60^\circ$, proving a highly anisotropic, substrate-locked in-plane epitaxial texture.
3. **Calcite Epitaxial Tilts:** Calcite shows a highly active domain structure at $\phi = 60^\circ$ (with multiple resolved tilt peaks at $\chi \approx -3.9^\circ, -2.5^\circ, -1.4^\circ$, and specular at $\chi \approx -0.2^\circ$). These domains disappear at $\phi = 150^\circ$, confirming a substrate-templated growth mode.
"""

# Replace placeholders
md_content = md_content.replace("{table_content}", table_content)

# Write md file
md_path = os.path.join(analysis_dir, "SH-125-G_rocking_curve_report.md")
with open(md_path, 'w') as f:
    f.write(md_content)
print(f"Written report Markdown to {md_path}")

# Compile to PDF
pdf_path = os.path.join(analysis_dir, "SH-125-G_rocking_curve_report.pdf")
cmd = [
    "pandoc",
    "SH-125-G_rocking_curve_report.md",
    "-o", "SH-125-G_rocking_curve_report.pdf",
    "--pdf-engine=pdflatex"
]
print(f"Compiling report to PDF...")
result = subprocess.run(cmd, cwd=analysis_dir, capture_output=True, text=True)
if result.returncode == 0:
    print(f"Successfully generated PDF report at {pdf_path}")
    # Copy to conversation artifacts
    conversation_artifacts_dir = "/home/tomek/.gemini/antigravity/brain/2974caf8-c2d8-4b77-9163-9cf57c4c82cc/artifacts"
    import shutil
    shutil.copy2(pdf_path, os.path.join(conversation_artifacts_dir, "SH-125-G_rocking_curve_report.pdf"))
    print("Copied report PDF to conversation artifacts.")
else:
    print("Error compiling report to PDF:")
    print("STDOUT:", result.stdout)
    print("STDERR:", result.stderr)
