# -*- coding: utf-8 -*-
"""
Azimuthal Report Generator
==========================
Author: Tomasz Stawski, tomasz.stawski@bam.de
Date: 2026-06-19
Version: 0.1

Description:
Compiles a technical report focusing on the phi-dependence of rocking curves and phase fractions.
"""
import os
import subprocess
import pandas as pd
import numpy as np

# Directories
base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
analysis_dir = os.path.join(base_dir, "data/processed/Rocking_Curves/08062026")
csv_path = os.path.join(analysis_dir, "all_samples_rocking_peaks_vs_phi.csv")

# Load data
df = pd.read_csv(csv_path)

# Programmatically generate tables for the report
def generate_sample_table(sample_name):
    df_s = df[df["Sample"] == sample_name].copy()
    
    # We want to display key peaks: Peak 2 (tilt ~ -1.8) and Peak 4/5 (tilt ~ 7.7)
    # Let's filter out the rows that have height > 0
    df_s = df_s[df_s["Net Height"] > 0]
    
    # Sort by Phi, then Peak Name
    df_s = df_s.sort_values(["Phi (degrees)", "Peak Name"])
    
    table_lines = [
        "| Phi (°C) | Peak Name | Center ($\\theta$) | Tilt Angle ($\\chi$) | FWHM ($^\circ$) | Height (counts) | Area (cts·deg) | Area/Base Ratio |",
        "| :---: | :--- | :---: | :---: | :---: | :---: | :---: | :---: |"
    ]
    
    for _, r in df_s.iterrows():
        ratio_str = f"{r['Area/Base Ratio']:.3%}" if not pd.isna(r['Area/Base Ratio']) else "N/A"
        table_lines.append(
            f"| {r['Phi (degrees)']:.0f} | {r['Peak Name']} | {r['Peak Center (Theta)']:.3f}° | {r['Tilt Angle (Chi)']:.3f}° | {r['FWHM (degrees)']:.3f}° | {r['Net Height']:.1f} | {r['Net Area (cts deg)']:.1f} | {ratio_str} |"
        )
    return "\n".join(table_lines)

table_b3 = generate_sample_table("SH-124-B3")
table_a = generate_sample_table("SH-125-A")

md_content = r"""---
title: "CaCO3 Thin Film Rocking Curve Analysis: Azimuthal (Phi) Variation Study"
date: "June 8, 2026"
geometry: "margin=1in"
fontsize: "11pt"
header-includes:
  - \usepackage{fancyhdr}
  - \pagestyle{fancy}
  - \fancyhead[CO,CE]{Phi-Variation Rocking Curve Report}
  - \fancyfoot[CO,CE]{Page \thepage}
---

# 1. Introduction and Dataset Structure
This report summarizes the analysis of a new set of X-ray diffraction (XRD) measurements carried out on **June 8, 2026** for two calcium carbonate ($CaCO_3$) thin film samples:
1. **`SH-124-B3`** (measured with $\phi$ rotations of $0^\circ, 30^\circ, 60^\circ, 90^\circ, 120^\circ, 150^\circ, 180^\circ$)
2. **`SH-125-A`** (measured with $\phi$ rotations of $0^\circ, 30^\circ, 60^\circ, 90^\circ, 120^\circ, 150^\circ$)

For each sample and at each azimuthal angle $\phi$, two types of scans were collected:
* **Symmetric $2\theta-\theta$ scans** (covering $2\theta$ from $27.0001^\circ$ to $35.0001^\circ$ in steps of $0.02^\circ$).
* **1D Rocking Curves** (covering $\theta$ from $4.705^\circ$ to $24.705^\circ$ in steps of $0.02^\circ$).

This azimuthal study addresses a major limitation of prior 1D rocking curve analyses: the inability to distinguish between an isotropic fiber texture (random in-plane orientation) and a preferred in-plane (azimuthal) alignment.

# 2. Phase Analysis and Phase Stability ($2\theta-\theta$ Scans)
Fitting Gaussian profiles to the background-subtracted symmetric scans yields the following structural findings:
* **`SH-124-B3`** is a **pure Calcite (104) film**. The Vaterite (110) peak at $2\theta \approx 32.8^\circ$ is absent or negligible across all $\phi$ angles. The Calcite (104) peak is centered at $2\theta \approx 29.41^\circ - 29.46^\circ$ with a FWHM of $\approx 0.28^\circ$.
* **`SH-125-A`** is a **mixed-phase Calcite-Vaterite film**. Both Calcite (104) (centered at $2\theta \approx 29.36^\circ - 29.43^\circ$, FWHM $\approx 0.22^\circ - 0.35^\circ$) and Vaterite (110) (centered at $2\theta \approx 32.76^\circ - 32.86^\circ$, FWHM $\approx 0.38^\circ - 0.45^\circ$) are resolved.

## 2.1 Stacked 2Theta Diffraction Patterns
The raw 2Theta diffraction scans stacked with offsets for both samples show the phase stability and presence of peaks as a function of the azimuthal angle $\phi$:

![SH-124-B3 Symmetric 2Theta Scans Stacked](SH-124-B3/SH-124-B3_2theta_scans_stacked.png){width=80%}

\newpage

![SH-125-A Symmetric 2Theta Scans Stacked](SH-125-A/SH-125-A_2theta_scans_stacked.png){width=80%}

## 2.2 Azimuthal Anisotropy of Phase Composition (SH-125-A)
In `SH-125-A`, the Calcite (104) peak area varies significantly as the sample is rotated in $\phi$, ranging from a minimum of $2,661.5$ counts$\cdot$deg at $\phi = 120^\circ$ to a maximum of $5,440.6$ counts$\cdot$deg at $\phi = 90^\circ$ (a factor of 2 variation). In contrast, the Vaterite (110) peak area remains relatively stable between $1,430$ and $1,950$ counts$\cdot$deg. 

This results in the **Calcite-to-Vaterite area ratio** shifting from $1.38$ (at $\phi = 120^\circ$) to $3.33$ (at $\phi = 90^\circ$). Because the physical phase fraction of a film cannot change upon simple rotation, this variation indicates a **strong in-plane preferred orientation (texture) in the Calcite crystallites**, while the Vaterite phase appears much more isotropic.

![SH-125-A Phase metrics vs. Phi](SH-125-A/SH-125-A_phase_metrics_vs_phi.png){width=80%}

\newpage

# 3. Rocking Curve Analysis and Volume Correction
To isolate the preferred orientation peak components, a global fitting was performed on the raw rocking curves in $\log_{10}(I)$ space. The model contains both the isotropic background baseline (a foot-print correction combined with a 3rd-order polynomial) and multiple Gaussian profiles representing the co-oriented grain populations.

## 3.1 Stacked Rocking Curves & Background Baselines (Log Scale)
The plots below show the raw rocking curves (dots) along with the fitted background baseline curves (red lines), stacked using multiplicative factors for each azimuthal angle $\phi$. This demonstrates the alignment and shape of the background model across all scans:

![SH-124-B3 Rocking Curves & Baselines Stacked](SH-124-B3/SH-124-B3_rocking_curves_stacked.png){width=80%}

\newpage

![SH-125-A Rocking Curves & Baselines Stacked](SH-125-A/SH-125-A_rocking_curves_stacked.png){width=80%}

\newpage

## 3.2 Stacked Baseline-Corrected Net Rocking Curves (Linear Scale)
The plots below show the baseline-corrected (net intensity) rocking curves stacked vertically on a linear scale. For each curve, the background baseline has been subtracted ($I_{\text{net}} = I_{\text{raw}} - I_{\text{baseline}}$), isolating the pure peak contributions. The dashed horizontal lines show the zero-baseline level for each azimuthal angle $\phi$:

![SH-124-B3 Stacked Baseline-Corrected Net Rocking Curves](SH-124-B3/SH-124-B3_net_rocking_curves_stacked.png){width=80%}

\newpage

![SH-125-A Stacked Baseline-Corrected Net Rocking Curves](SH-125-A/SH-125-A_net_rocking_curves_stacked.png){width=80%}

\newpage

## 3.3 Stacked Global Peak Fits (Log Scale)
The plots below show the complete global fits in log scale. The raw intensity is represented by dots, the full fitted model (baseline + all Gaussian peaks) is represented by the solid black line, and the individual resolved peaks are plotted as colored dashed curves rising above the background baseline. This representation clearly displays how each peak, even a minor one, is resolved relative to the background:

![SH-124-B3 Rocking Curve Global Fits Stacked](SH-124-B3/SH-124-B3_rocking_residuals_stacked.png){width=80%}

\newpage

![SH-125-A Rocking Curve Global Fits Stacked](SH-125-A/SH-125-A_rocking_residuals_stacked.png){width=80%}

\newpage

## 3.4 Azimuthal Dependence of Residual Peaks (In-Plane Texture)
The intensities of the residual tilt peaks exhibit dramatic changes as a function of the azimuthal angle $\phi$. This is visually summarized in the peak areas vs. $\phi$ plot:

![Comparison of Key Tilt Peak Areas vs. Phi](samples_key_peaks_comparison_vs_phi.png){width=80%}

### 1. The $\chi \approx -2.8^\circ$ and $\chi \approx -1.8^\circ$ Components (Peak 2a and Peak 2b)
Our updated fitting uncovers two distinct, co-existing tilt states in the negative tilt region:
* **`SH-124-B3`**: Both components are highly active in the $\phi = 30^\circ - 120^\circ$ range. The $\chi \approx -1.8^\circ$ component (Peak 2b) peaks at $\phi = 60^\circ$ (height $= 2,350.9$), while the newly resolved $\chi \approx -2.8^\circ$ component (Peak 2a) peaks at $\phi = 120^\circ$ (height $= 2,159.9$).
* **`SH-125-A`**: A similar trend is observed, but both components peak synchronously at $\phi = 60^\circ$, with Peak 2a reaching height $= 1,022.1$ counts and Peak 2b reaching height $= 1,046.1$ counts.

### 2. The $\chi \approx +7.7^\circ$ Component (Conical Tilt / Helical Columns)
* **`SH-124-B3`**: The $+7.7^\circ$ component (Peak 5) peaks strongly at $\phi = 120^\circ$ (height $= 957.4$ counts), while remaining moderately strong at $\phi = 0^\circ, 30^\circ, 60^\circ, 90^\circ$.
* **`SH-125-A`**: The $+7.7^\circ$ component (Peak 4) is resolved only at $\phi = 0^\circ$ (height $= 594.6$) and $\phi = 30^\circ$ (height $= 146.1$), disappearing completely at all other $\phi$ angles.

## 3.5 Physical and Crystallographic Interpretation
The strong modulation of peak intensities with $\phi$ confirms that the co-oriented grain populations (crystallite statistics) are **not azimuthally isotropic**. 
1. **Absence of Fiber Texture:** If the tilted crystallites (such as the helical columns observed in SEM) were randomly rotated around the surface normal, the rocking curves would be identical at all $\phi$ angles. Instead, we observe specific azimuthal "spotting".
2. **In-Plane Epitaxy/Alignment:** The substrate template dictates not only the tilt angle $\chi$ relative to the surface normal, but also imposes a preferred azimuthal alignment $\phi$. This suggests a substrate-guided epitaxial relation where specific crystallographic directions of the calcite nuclei align with the substrate's in-plane directions.
3. **Helical Cones and Asymmetric Chirality:** The different $\phi$ dependencies of the negative tilt ($-1.8^\circ$) and positive tilt ($+7.7^\circ$) components suggest a complex spatial arrangement, likely related to the helical columnar growth twisting in a preferred direction or a tilt in the columnar growth axis itself.

\newpage

# 4. 2D Polar Texture Plots (Pole Figures)
To visualize the complete in-plane alignment, the 1D rocking curves measured at different azimuthal angles $\phi$ were mapped into a single 2D polar projection (pole-figure-like texture plot). Negative tilts ($\chi < 0$) at angle $\phi$ are mapped to a positive radius at the opposite azimuth ($\phi + 180^\circ$).

The raw texture plots show a strong specular orientation at the center ($\chi = 0$), whereas the baseline-subtracted (net residual) texture plots clearly resolve the preferred off-axis orientation spots representing the co-oriented grain populations:

![SH-124-B3 2D Polar Texture Plots (Raw & Net)](SH-124-B3/SH-124-B3_texture_raw.png){width=45%} ![SH-124-B3 Net Texture](SH-124-B3/SH-124-B3_texture_net.png){width=45%}

\newpage

![SH-125-A 2D Polar Texture Plots (Raw & Net)](SH-125-A/SH-125-A_texture_raw.png){width=45%} ![SH-125-A Net Texture](SH-125-A/SH-125-A_texture_net.png){width=45%}

\newpage

# 5. Fitted Rocking Curve Peak Metrics

### 5.1 Sample `SH-124-B3`
{table_b3}

### 5.2 Sample `SH-125-A`
{table_a}

# 6. Conclusions
1. **In-Plane Anisotropy Resolved:** The addition of $\phi$ rotation measurements successfully resolves the azimuthal anisotropy of the film. The preferred orientation components are highly localized in-plane, confirming a template-guided nucleation.
2. **Phase Anisotropy in SH-125-A:** In the Calcite-Vaterite mixed film, Calcite shows a high degree of in-plane alignment, leading to a strong modulation of the apparent Calcite/Vaterite ratio with $\phi$ in symmetric scans.
3. **Helical columnar alignment:** The helical columns (associated with the $\approx -1.8^\circ$ and $\approx 7.7^\circ$ tilt states) grow along preferred azimuthal directions relative to the substrate template.
"""

# Replace placeholders
md_content = md_content.replace("{table_b3}", table_b3)
md_content = md_content.replace("{table_a}", table_a)

# Write md file
md_path = os.path.join(analysis_dir, "rocking_curve_phi_variation_report.md")
with open(md_path, 'w') as f:
    f.write(md_content)
print(f"Written report Markdown to {md_path}")

# Compile to PDF
pdf_path = os.path.join(analysis_dir, "rocking_curve_phi_variation_report.pdf")
cmd = [
    "pandoc",
    "rocking_curve_phi_variation_report.md",
    "-o", "rocking_curve_phi_variation_report.pdf",
    "--pdf-engine=pdflatex"
]
print(f"Compiling report to PDF...")
result = subprocess.run(cmd, cwd=analysis_dir, capture_output=True, text=True)
if result.returncode == 0:
    print(f"Successfully generated PDF report at {pdf_path}")
else:
    print("Error compiling report to PDF:")
    print("STDOUT:", result.stdout)
    print("STDERR:", result.stderr)
