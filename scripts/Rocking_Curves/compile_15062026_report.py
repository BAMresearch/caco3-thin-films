# -*- coding: utf-8 -*-
"""
Report Generator for 15 June 2026 Run
=====================================
Author: Tomasz Stawski, tomasz.stawski@bam.de
Date: 2026-06-19
Version: 0.1

Description:
Compiles the technical summary report for rocking curve measurements from 15 June 2026.
"""
import os
import subprocess
import pandas as pd
import numpy as np

# Directories
base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
analysis_dir = os.path.join(base_dir, "data/processed/Rocking_Curves/15062026/SH-104-1")
csv_path = os.path.join(analysis_dir, "SH-104-1_rocking_peaks_vs_phi.csv")

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
title: "CaCO3 Thin Film Rocking Curve Analysis: Azimuthal (Phi) Variation Study (Sample SH-104-1)"
date: "June 15, 2026"
geometry: "margin=1in"
fontsize: "11pt"
header-includes:
  - \usepackage{fancyhdr}
  - \pagestyle{fancy}
  - \fancyhead[CO,CE]{SH-104-1 Phi-Variation Rocking Curve Report}
  - \fancyfoot[CO,CE]{Page \thepage}
---

# 1. Introduction and Dataset Structure
This report summarizes the analysis of a new set of X-ray diffraction (XRD) measurements collected on **June 15, 2026** for the calcium carbonate ($CaCO_3$) thin film sample **`SH-104-1`**. 

The measurement was carried out at six azimuthal angles: $\phi = 0^\circ, 30^\circ, 60^\circ, 90^\circ, 120^\circ, 150^\circ$. For each angle, two types of scans were collected:
* **Symmetric $2\theta-\theta$ scans** (covering $2\theta$ from $27.0001^\circ$ to $35.0001^\circ$).
* **1D Rocking Curves** (covering $\theta$ from $4.705^\circ$ to $24.705^\circ$).

During the analysis, the file `Rocking_60.brml` was found to be empty (0 bytes). To prevent loss of information, a custom binary parser was implemented to recover the 1001 intensity points directly from the raw Bruker binary file `Rocking_60.raw` starting at offset 1488 (verified by comparing with non-empty files). This successful recovery allowed a complete 6-angle azimuthal study.

# 2. Phase Analysis and Phase Stability ($2\theta-\theta$ Scans)
Gaussian profile fitting of the background-subtracted symmetric $2\theta-\theta$ scans reveals several key structural insights:
1. **Calcite (104) Stability:** Calcite is highly stable and detected at all phi angles. The peak is centered at $2\theta \approx 29.40^\circ - 29.43^\circ$ with a narrow FWHM of $\approx 0.26^\circ - 0.28^\circ$.
2. **Azimuthally Textured Vaterite (110):** Vaterite (110) at $2\theta \approx 33.05^\circ$ is **only detected at $\phi = 60^\circ$** (with an area of $1,585.8$ counts$\cdot$deg and a Calcite/Vaterite area ratio of 5.79). At all other phi angles, the Vaterite signal is negligible and falls below the detection threshold. 

The fact that Vaterite is visible exclusively at a single azimuthal rotation ($\phi = 60^\circ$) suggests a highly anisotropic, strongly in-plane oriented (epitaxial) Vaterite phase, whose reflecting planes only align with the diffraction condition at this specific rotation.

## 2.1 Stacked 2Theta Diffraction Patterns
The raw 2Theta diffraction scans stacked with vertical offsets show the presence of Calcite across all angles and the isolated Vaterite reflection at $\phi = 60^\circ$:

![SH-104-1 Symmetric 2Theta Scans Stacked](SH-104-1_2theta_scans_stacked.png){width=85%}

\newpage

## 2.2 Azimuthal Dependence of Phase Composition
The plot below tracks the peak area of Calcite and Vaterite as a function of the azimuthal angle $\phi$. This visually highlights the emergence of the Vaterite phase at $\phi = 60^\circ$:

![SH-104-1 Phase metrics vs. Phi](SH-104-1_phase_metrics_vs_phi.png){width=85%}

# 3. Rocking Curve Analysis and Volume Correction
To isolate the preferred orientation peak components, a global fitting was performed on the raw rocking curves in $\log_{10}(I)$ space. The model contains both the isotropic background baseline (a footprint correction combined with a 3rd-order polynomial) and multiple Gaussian profiles representing the co-oriented grain populations.

## 3.1 Stacked Rocking Curves & Background Baselines (Log Scale)
The plot below shows the raw rocking curves along with the fitted background baseline curves (red lines), stacked using multiplicative factors for each azimuthal angle $\phi$:

![SH-104-1 Rocking Curves & Baselines Stacked](SH-104-1_rocking_curves_stacked.png){width=85%}

\newpage

## 3.2 Stacked Baseline-Corrected Net Rocking Curves (Linear Scale)
The plot below shows the baseline-corrected (net intensity) rocking curves stacked vertically on a linear scale. The dashed horizontal lines show the zero-baseline level for each azimuthal angle $\phi$:

![SH-104-1 Stacked Baseline-Corrected Net Rocking Curves](SH-104-1_side_by_side.png){width=85%}

\newpage

## 3.3 Stacked Global Peak Fits (Log Scale)
The plot below shows the complete global fits in log scale, showing the raw intensity (dots), the full fitted model (solid black line), and the individual resolved peaks (dashed curves):

![SH-104-1 Rocking Curve Global Fits Stacked](SH-104-1_rocking_residuals_stacked.png){width=85%}

\newpage

# 4. 2D Polar Texture Plots (Pole Figures)
To visualize the complete in-plane alignment, the 1D rocking curves measured at different azimuthal angles $\phi$ were mapped into a single 2D polar projection (pole-figure-like texture plot). Negative tilts ($\chi < 0$) at angle $\phi$ are mapped to a positive radius at the opposite azimuth ($\phi + 180^\circ$).

The raw texture plot shows a strong specular orientation at the center ($\chi = 0$), whereas the baseline-subtracted (net residual) texture plot clearly resolves the preferred off-axis orientation spots representing the co-oriented grain populations:

![SH-104-1 2D Polar Texture Plots (Raw & Net)](SH-104-1_texture_raw.png){width=48%} ![SH-104-1 Net Texture](SH-104-1_texture_net.png){width=48%}

# 5. Fitted Rocking Curve Peak Metrics
The table below lists the parameters of all active reflections resolved from the rocking curves at each azimuthal rotation angle:

{table_content}

# 6. Conclusions
1. **Successful Binary Recovery:** Implementation of a Bruker RAW binary parser successfully recovered the missing intensity dataset for $\phi = 60^\circ$ (`Rocking_60.brml`), preserving the integrity of this study.
2. **Azimuthally Textured Vaterite Co-growth:** A highly oriented Vaterite phase (diffracting solely at $\phi = 60^\circ$) coexists with Calcite, indicating a strong template-guided epitaxial growth of Vaterite on the substrate.
3. **Epitaxial Calcite Domains:** The modulation of the rocking curve peaks with azimuthal rotation shows a strong preferred orientation. Specifically, the negative tilt reflections (e.g. Peak 2b at $\chi \approx -1.1^\circ$ and Peak 3 specular at $\chi \approx 0.1^\circ$) are dominant and active at specific rotations, confirming that the crystallographic axes are locked in-plane by substrate epitaxy.
"""

# Replace placeholders
md_content = md_content.replace("{table_content}", table_content)

# Write md file
md_path = os.path.join(analysis_dir, "SH-104-1_rocking_curve_report.md")
with open(md_path, 'w') as f:
    f.write(md_content)
print(f"Written report Markdown to {md_path}")

# Compile to PDF
pdf_path = os.path.join(analysis_dir, "SH-104-1_rocking_curve_report.pdf")
cmd = [
    "pandoc",
    "SH-104-1_rocking_curve_report.md",
    "-o", "SH-104-1_rocking_curve_report.pdf",
    "--pdf-engine=pdflatex"
]
print(f"Compiling report to PDF...")
result = subprocess.run(cmd, cwd=analysis_dir, capture_output=True, text=True)
if result.returncode == 0:
    print(f"Successfully generated PDF report at {pdf_path}")
    # Copy to conversation artifacts
    conversation_artifacts_dir = "/home/tomek/.gemini/antigravity/brain/2974caf8-c2d8-4b77-9163-9cf57c4c82cc/artifacts"
    import shutil
    shutil.copy2(pdf_path, os.path.join(conversation_artifacts_dir, "SH-104-1_rocking_curve_report.pdf"))
    print("Copied report PDF to conversation artifacts.")
else:
    print("Error compiling report to PDF:")
    print("STDOUT:", result.stdout)
    print("STDERR:", result.stderr)
