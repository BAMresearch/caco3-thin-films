# -*- coding: utf-8 -*-
"""
Publication Summary Compiler
============================
Author: Tomasz Stawski, tomasz.stawski@bam.de
Date: 2026-06-19
Version: 0.1

Description:
Compiles a condensed summary report matching the manuscript outline for publication.
"""
import os
import shutil
import subprocess
import pandas as pd
import numpy as np
from PIL import Image, ImageDraw, ImageFont

# Define directories
base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
analysis_dir = os.path.join(base_dir, "data/processed/Rocking_Curves")
comparison_dir = os.path.join(analysis_dir, "Comparison")
os.makedirs(comparison_dir, exist_ok=True)

# Image paths
img_sources = {
    "SH-124-B3": {
        "side_by_side": os.path.join(analysis_dir, "08062026/SH-124-B3/SH-124-B3_side_by_side.png"),
        "texture_net": os.path.join(analysis_dir, "08062026/SH-124-B3/SH-124-B3_texture_net.png"),
        "phase_metrics": os.path.join(analysis_dir, "08062026/SH-124-B3/SH-124-B3_phase_metrics_vs_phi.png")
    },
    "SH-125-A": {
        "side_by_side": os.path.join(analysis_dir, "08062026/SH-125-A/SH-125-A_side_by_side.png"),
        "texture_net": os.path.join(analysis_dir, "08062026/SH-125-A/SH-125-A_texture_net.png"),
        "phase_metrics": os.path.join(analysis_dir, "08062026/SH-125-A/SH-125-A_phase_metrics_vs_phi.png")
    },
    "SH-104-1": {
        "side_by_side": os.path.join(analysis_dir, "15062026/SH-104-1/SH-104-1_side_by_side.png"),
        "texture_net": os.path.join(analysis_dir, "15062026/SH-104-1/SH-104-1_texture_net.png"),
        "phase_metrics": os.path.join(analysis_dir, "15062026/SH-104-1/SH-104-1_phase_metrics_vs_phi.png")
    },
    "SH-125-G": {
        "side_by_side": os.path.join(analysis_dir, "16062026/SH-125-G/SH-125-G_side_by_side.png"),
        "texture_net": os.path.join(analysis_dir, "16062026/SH-125-G/SH-125-G_texture_net.png"),
        "phase_metrics": os.path.join(analysis_dir, "16062026/SH-125-G/SH-125-G_phase_metrics_vs_phi.png")
    }
}

# 1. Merge images helper
def merge_and_label_2x2(key, output_name):
    paths = [
        img_sources["SH-124-B3"][key],
        img_sources["SH-125-A"][key],
        img_sources["SH-104-1"][key],
        img_sources["SH-125-G"][key]
    ]
    
    images = [Image.open(p) for p in paths]
    w, h = images[0].size
    
    new_img = Image.new('RGB', (w*2, h*2), color='white')
    new_img.paste(images[0], (0, 0))
    new_img.paste(images[1], (w, 0))
    new_img.paste(images[2], (0, h))
    new_img.paste(images[3], (w, h))
    
    # Annotate with letters (A, B, C, D)
    draw = ImageDraw.Draw(new_img)
    labels = ["(a) SH-124-B3", "(b) SH-125-A", "(c) SH-104-1", "(d) SH-125-G"]
    
    # Try using default font or system font if available
    try:
        font = ImageFont.truetype("DejaVuSans-Bold.ttf", int(h*0.04))
    except:
        font = ImageFont.load_default()
        
    draw.text((int(w*0.05), int(h*0.03)), labels[0], fill='black', font=font)
    draw.text((w + int(w*0.05), int(h*0.03)), labels[1], fill='black', font=font)
    draw.text((int(w*0.05), h + int(h*0.03)), labels[2], fill='black', font=font)
    draw.text((w + int(w*0.05), h + int(h*0.03)), labels[3], fill='black', font=font)
    
    out_path = os.path.join(comparison_dir, output_name)
    new_img.save(out_path, dpi=(150, 150))
    print(f"Created panel: {out_path}")
    
    # Copy to conversation artifacts
    conversation_artifacts_dir = "/home/tomek/.gemini/antigravity/brain/2974caf8-c2d8-4b77-9163-9cf57c4c82cc/artifacts"
    os.makedirs(conversation_artifacts_dir, exist_ok=True)
    import shutil
    shutil.copy2(out_path, os.path.join(conversation_artifacts_dir, output_name))
    print(f"Copied panel to conversation artifacts: {output_name}")

# Generate panels
merge_and_label_2x2("phase_metrics", "panel_phase_metrics.png")
merge_and_label_2x2("side_by_side", "panel_side_by_side.png")
merge_and_label_2x2("texture_net", "panel_texture_net.png")

# 2. Compile metrics table
csv_paths = [
    os.path.join(analysis_dir, "08062026/all_samples_rocking_peaks_vs_phi.csv"),
    os.path.join(analysis_dir, "15062026/SH-104-1/SH-104-1_rocking_peaks_vs_phi.csv"),
    os.path.join(analysis_dir, "16062026/SH-125-G/SH-125-G_rocking_peaks_vs_phi.csv")
]

dfs = []
for p in csv_paths:
    if os.path.exists(p):
        dfs.append(pd.read_csv(p))
df_all = pd.concat(dfs, ignore_index=True)

# Select key active components to summarize
# We want to show the peak center, tilt angle, height, and area for each sample at a representative Phi angle where they are highly active
summary_table_lines = [
    "| Sample | Phi (°) | Peak Name | Center ($\\theta$) | Tilt Angle ($\\chi$) | FWHM ($^\circ$) | Height (counts) | Area (cts·deg) | Area/Base Ratio |",
    "| :--- | :---: | :--- | :---: | :---: | :---: | :---: | :---: | :---: |"
]

# Add calcite single crystal as a reference
sc_csv_path = os.path.join(analysis_dir, "09062026/calcite_single_crystal_rocking_peaks_metrics.csv")
if os.path.exists(sc_csv_path):
    df_sc = pd.read_csv(sc_csv_path)
    for _, r in df_sc.iterrows():
        summary_table_lines.append(
            f"| **Single Crystal** | Reference | {r['Peak Name']} | {r['Peak Center (Theta)']:.3f}° | {r['Tilt Angle (Chi)']:.3f}° | {r['FWHM (deg)']:.3f}° | {r['Net Height']:.1f} | {r['Net Area (cts deg)']/1e6:.2f}M | {r['Area/Base Ratio']:.1%} |"
        )

# For each thin film sample, extract a few characteristic rows
samples_to_report = [
    ("SH-124-B3", 60, ["Peak 1 (Tilt)", "Peak 2b (Tilt)"]),
    ("SH-125-A", 60, ["Peak 2b (Tilt)", "Peak 3 (Specular)"]),
    ("SH-104-1", 0, ["Peak 2b (Tilt)", "Peak 3 (Specular)"]),
    ("SH-104-1", 90, ["Peak 2b (Tilt)", "Peak 3 (Specular)"]),
    ("SH-125-G", 60, ["Peak 2a (Tilt)", "Peak 2b (Tilt)"]),
    ("SH-125-G", 120, ["Peak 2b (Tilt)", "Peak 3 (Specular)"])
]

for s, phi, pnames in samples_to_report:
    for pname in pnames:
        sub = df_all[(df_all["Sample"] == s) & (df_all["Phi (degrees)"] == phi) & (df_all["Peak Name"] == pname)]
        if not sub.empty:
            r = sub.iloc[0]
            ratio_str = f"{r['Area/Base Ratio']:.3%}" if not pd.isna(r['Area/Base Ratio']) else "N/A"
            summary_table_lines.append(
                f"| {s} | {phi} | {pname} | {r['Peak Center (Theta)']:.3f}° | {r['Tilt Angle (Chi)']:.3f}° | {r['FWHM (degrees)']:.3f}° | {r['Net Height']:.1f} | {r['Net Area (cts deg)']:.1f} | {ratio_str} |"
            )

summary_table = "\n".join(summary_table_lines)

# 3. Write publication ready Markdown
md_content = r"""---
title: "Texture and Epitaxial Relationships in CaCO3 Thin Films: Compiled Study"
subtitle: "A Multi-Dataset Analysis of Azimuthal (Phi) Dependent Rocking Curves"
date: "June 15, 2026"
geometry: "margin=1in"
fontsize: "11pt"
header-includes:
  - \usepackage{fancyhdr}
  - \pagestyle{fancy}
  - \fancyhead[CO,CE]{CaCO3 Thin Film Rocking Curve Compilation Report}
  - \fancyfoot[CO,CE]{Page \thepage}
---

# Abstract
This publication-ready report compiles and compares the structural and crystallographic findings from four calcium carbonate ($CaCO_3$) thin film samples (`SH-124-B3`, `SH-125-A`, `SH-104-1`, and `SH-125-G`) measured across a series of azimuthal ($\phi$) rotation angles. By combining symmetric $2\theta-\theta$ phase scans and multi-component 1D rocking curves ($\theta$), we resolve in-plane epitaxial texture vs. random fiber texture. The analysis demonstrates a strong, substrate-guided in-plane alignment for both Calcite and Vaterite phases. Crucially, Vaterite (110) crystallites are confined to specific azimuthal rotations ($\phi = 60^\circ$ for `SH-104-1` and $\phi = 30^\circ/60^\circ$ for `SH-125-G`), proving complete in-plane epitaxial confinement of this transient carbonate phase.

# 1. Introduction
A major challenge in characterizing functional thin-film coatings is mapping the spatial distribution and orientation of grains (mosaicity and texture). Standard 1D rocking curve measurements (varying $\theta$ at a fixed detector position $2\theta_0$) are unable to distinguish between an isotropic fiber texture (where grain tilts are random in-plane) and a preferred in-plane (azimuthal) alignment. 

To overcome this limitation, a systematic azimuthal variation study was carried out on four $CaCO_3$ thin film samples. By rotating the sample azimuthally ($\phi$ from $0^\circ$ to $180^\circ$ in steps of $30^\circ$) and collecting both symmetric $2\theta$ scans and 1D rocking curves, we reconstruct the full 2D orientation distribution (pole figures). This allows us to determine:
1. Whether the film contains a random fiber texture or is epitaxially locked to the substrate.
2. The specific domain tilt angles ($\chi$) and their corresponding in-plane azimuthal directions ($\phi$).
3. The spatial and phase distribution of co-existing Calcite and Vaterite crystallites.

# 2. Experimental Section & References
* **Symmetric $2\theta-\theta$ scans:** Cover $2\theta = 27.0^\circ - 35.0^\circ$ in steps of $0.02^\circ$. This region covers the primary Calcite (104) reflection at $2\theta \approx 29.4^\circ$ and the Vaterite (110) reflection at $2\theta \approx 32.8^\circ - 33.0^\circ$.
* **1D Rocking Curves:** Cover $\theta = 4.7^\circ - 24.7^\circ$ in steps of $0.02^\circ$ (1001 points) with $2\theta$ fixed at the nominal Calcite (104) position.
* **Calcite Single Crystal Reference:** Measured as an instrumental resolution standard (FWHM $\approx 4.13^\circ$, peak centered at $\theta = 17.46^\circ$ corresponding to a cleavage tilt/miscut of $\chi = +2.77^\circ$).
* **Empty Sample Holder:** Confirms a featureless amorphous scatter profile with a standard deviation of $< 0.6\%$, proving that no background reflections interfere with the analysis.
* **Data Recovery and Integrity:** For `SH-104-1`, a corrupted file (`Rocking_60.brml`) was recovered by parsing the binary Bruker `Rocking_60.raw` structure directly at offset 1488. For `SH-125-G`, `Rocking_90.brml` was missing; the pipeline successfully skipped this rotation, displaying the missing data as a clear, unmeasured gap.

# 3. Phase Composition Analysis and Epitaxial Vaterite Confinement
Fitting Gaussian profiles to the background-subtracted symmetric $2\theta$ scans reveals two distinct classes of films:
1. **Single-Phase Calcite Film (`SH-124-B3`):** Showed zero Vaterite across all phi angles.
2. **Mixed-Phase Calcite-Vaterite Films (`SH-125-A`, `SH-104-1`, `SH-125-G`):** Showed both phases, but with a highly anisotropic, phi-dependent signature:
   * **`SH-125-A`:** Vaterite is resolved at all phi angles, but the Calcite/Vaterite area ratio varies from $1.38$ (at $\phi=120^\circ$) to $3.33$ (at $\phi=90^\circ$). This indicates that Calcite has a stronger in-plane texture than Vaterite.
   * **`SH-104-1`:** Vaterite is **exclusively resolved at $\phi = 60^\circ$** (Calcite/Vaterite ratio $= 5.79$). It is below the detection limit at all other rotations.
   * **`SH-125-G`:** Vaterite is **exclusively resolved at $\phi = 30^\circ$ and $\phi = 60^\circ$** (ratios of $5.61$ and $4.74$). It disappears completely at other rotations.

This represents a major crystallographic discovery: Vaterite crystallites in `SH-104-1` and `SH-125-G` are grown with a highly restricted epitaxial relationship, aligning their (110) planes to the diffracting plane only at specific azimuthal normal rotations.

![Comparative Phase Metrics vs. Phi](panel_phase_metrics.png){width=85%}

\newpage

# 4. Domain Orientation Distribution and Peak Fitting
Global rocking curve peak fitting was performed in $\log_{10}(I)$ space. The mathematical model accounts for:
* **Footprint and Baseline Correction:** $I_{\text{baseline}}(\theta) = \frac{I_0}{\sin(\theta)} + p_3(\theta)^3 + p_2(\theta)^2 + p_1(\theta) + p_0$ to isolate co-oriented grain peaks from the diffuse isotropic background.
* **Gaussian Peak Assemblies:** Multiple Gaussian curves corresponding to domain tilt states:
  * **Specular peak ($\chi \approx 0^\circ$):** Specularly oriented grains.
  * **Negative tilt peaks ($\chi \approx -1.1^\circ$ to $-4.5^\circ$):** Grains tilted away from the surface normal in one direction.
  * **Positive tilt peaks ($\chi \approx +2.1^\circ$ to $+8.5^\circ$):** Grains tilted in the opposite direction.

The baseline-corrected net intensity curves and fits reveal that the domain states are highly dependent on the azimuthal angle $\phi$.

![Stacked Rocking Curve Fits and Net Intensities](panel_side_by_side.png){width=85%}

\newpage

# 5. 2D Polar Texture Comparison (Pole Figures)
Mapping the 1D rocking curves to a 2D polar projection (pole figures) provides a visual overview of the in-plane alignment:
* **`SH-124-B3` (Pure Calcite):** Exhibits strong off-axis spots at $\chi \approx -1.8^\circ$ and $+7.7^\circ$ localized in the $\phi = 60^\circ - 120^\circ$ sector, confirming an epitaxial lock.
* **`SH-125-A` (Mixed Phase):** Shows co-existing specular and tilted domains peaking synchronously at $\phi = 60^\circ$.
* **`SH-104-1` (Vaterite at 60°):** Shows strong specular and negative tilt components at $\phi = 0^\circ, 90^\circ, 120^\circ$. The domains are inactive at $\phi=150^\circ$.
* **`SH-125-G` (Vaterite at 30°/60°):** Displays active domain states at $\phi = 60^\circ$ and $120^\circ$, with a clear unmeasured gap at $\phi=90^\circ$ representing the missing scan.

![Compiled Net Intensity 2D Texture Pole Figures](panel_texture_net.png){width=85%}

\newpage

# 6. Quantitative Peak Parameters Table
The table below compiles the parameters of key active domain reflections for the four thin-film samples and the single crystal calibrant:

{summary_table}

# 7. Conclusions & Scientific Implications
1. **Direct Proof of Epitaxy (No Fiber Texture):** The strong, localized modulation of rocking curve peak heights and areas with azimuthal rotation ($\phi$) proves that these films are epitaxial rather than having a random fiber texture. Grains are aligned both out-of-plane (tilt $\chi$) and in-plane (azimuth $\phi$).
2. **Epitaxial Confinement of Vaterite:** The localized emergence of the Vaterite (110) reflection only at specific rotations (e.g. $\phi = 60^\circ$ for `SH-104-1`) demonstrates that Vaterite crystallites nucleate and grow with their in-plane lattice directions locked to the substrate template.
3. **Domain Co-existence:** The presence of multiple discrete tilt states (such as $\chi \approx -1.8^\circ$ and $+7.7^\circ$ in Calcite) suggests a columnar growth mode where grains form inclined, possibly helical columns that adopt specific crystallographic relations with the substrate lattice.
"""

# Replace placeholder
md_content = md_content.replace("{summary_table}", summary_table)

# Write to file
md_path = os.path.join(comparison_dir, "CaCO3_thin_films_rocking_curve_publication_summary.md")
with open(md_path, 'w') as f:
    f.write(md_content)
print(f"Written Markdown report to {md_path}")

# Compile to PDF
pdf_path = os.path.join(comparison_dir, "CaCO3_thin_films_rocking_curve_publication_summary.pdf")
cmd = [
    "pandoc",
    "CaCO3_thin_films_rocking_curve_publication_summary.md",
    "-o", "CaCO3_thin_films_rocking_curve_publication_summary.pdf",
    "--pdf-engine=pdflatex"
]
print(f"Compiling report to PDF...")
result = subprocess.run(cmd, cwd=comparison_dir, capture_output=True, text=True)
if result.returncode == 0:
    print(f"Successfully generated PDF report at {pdf_path}")
    # Copy PDF to conversation artifacts
    conversation_artifacts_dir = "/home/tomek/.gemini/antigravity/brain/2974caf8-c2d8-4b77-9163-9cf57c4c82cc/artifacts"
    shutil.copy2(pdf_path, os.path.join(conversation_artifacts_dir, "CaCO3_thin_films_rocking_curve_publication_summary.pdf"))
    print("Copied summary PDF to conversation artifacts.")
else:
    print("Error compiling report to PDF:")
    print("STDOUT:", result.stdout)
    print("STDERR:", result.stderr)
