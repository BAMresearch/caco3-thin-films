# -*- coding: utf-8 -*-
"""
Multi-Batch Data Integrator
===========================
Author: Tomasz Stawski, tomasz.stawski@bam.de
Date: 2026-06-19
Version: 0.1

Description:
Combines processed peak metrics from the 08 June and 09 June 2026 measurement series.
"""
import os
import shutil
import subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Paths
base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
src_analysis = os.path.join(base_dir, "data/processed/Rocking_Curves/08062026")
dst_analysis = os.path.join(base_dir, "data/processed/Rocking_Curves/Reference")
data_0806_dir = os.path.join(base_dir, "data/processed/Rocking_Curves")
data_0906_dir = os.path.join(base_dir, "data/processed/Rocking_Curves/Reference")
artifacts_dir = "/home/tomek/.gemini/antigravity/brain/297360fa-0f2b-4ba6-96f5-6c818b856b29/artifacts"

# 1. Copy 08062026 analysis to 08062026_integrated to preserve original analysis
print("Copying analysis directory to 08062026_integrated...")
if os.path.exists(dst_analysis):
    shutil.rmtree(dst_analysis)
shutil.copytree(src_analysis, dst_analysis)
print(f"Copied to {dst_analysis}")

# 2. Load data for profile shape comparison (Plot A)
print("Loading profile data for shape comparison...")
df_crystal = pd.read_csv(os.path.join(data_0906_dir, "calcite_single_crystal_corrected_rocking.csv"))
df_b3 = pd.read_csv(os.path.join(data_0806_dir, "SH-124-B3/SH-124-B3_corrected_rocking_60.csv"))
df_a = pd.read_csv(os.path.join(data_0806_dir, "SH-125-A/SH-125-A_corrected_rocking_60.csv"))

# Calculate Chi
chi_crystal = df_crystal["Theta (degrees)"] - 14.68625
chi_b3 = df_b3["Theta (degrees)"] - 14.716665
chi_a = df_a["Theta (degrees)"] - 14.685524

# Normalize net intensity
net_crystal_norm = df_crystal["Corrected Net Intensity"] / df_crystal["Corrected Net Intensity"].max()
net_b3_norm = df_b3["Corrected Net Intensity"] / df_b3["Corrected Net Intensity"].max()
net_a_norm = df_a["Corrected Net Intensity"] / df_a["Corrected Net Intensity"].max()

# Plot A: Rocking Curve Profile Shape Comparison
plt.figure(figsize=(10, 6))
plt.plot(chi_crystal, net_crystal_norm, '-', color='#800080', linewidth=3.0, label='Calcite Single Crystal (FWHM = 4.13°)')
plt.plot(chi_b3, net_b3_norm, '--', color='#2ca02c', linewidth=1.5, label='SH-124-B3 at Phi = 60° (Net normalized)')
plt.plot(chi_a, net_a_norm, '--', color='#1f77b4', linewidth=1.5, label='SH-125-A at Phi = 60° (Net normalized)')

plt.xlabel("Tilt Angle Chi (degrees)", fontsize=12)
plt.ylabel("Normalized Net Intensity (a.u.)", fontsize=12)
plt.title("Rocking Curve Shape Comparison: Mosaicity vs. Thin-Film Texture", fontsize=13, fontweight='bold')
plt.xlim(-10, 10)
plt.ylim(-0.1, 1.1)
plt.grid(True, linestyle='--', alpha=0.5)
plt.legend(fontsize=10, loc='upper right')
plt.tight_layout()

plot_shape_path = os.path.join(dst_analysis, "rocking_curve_shape_comparison.png")
plt.savefig(plot_shape_path, dpi=150)
plt.close()
print(f"Saved profile shape comparison plot to {plot_shape_path}")

# 3. Load data for 2Theta symmetric scan comparison (Plot B)
print("Loading 2Theta symmetric scan data...")
def load_xy(path):
    data = []
    with open(path, 'r') as f:
        f.readline()  # skip header
        for line in f:
            parts = line.strip().split()
            if len(parts) == 2:
                data.append([float(parts[0]), float(parts[1])])
    return np.array(data)

arr_holder = load_xy(os.path.join(data_0906_dir, "sample_holder_2theta_exported.xy"))
arr_b3 = load_xy(os.path.join(data_0806_dir, "SH-124-B3/SH-124-B3_2Theta_0_exported.xy"))
arr_a = load_xy(os.path.join(data_0806_dir, "SH-125-A/SH-125-A_2Theta_0_exported.xy"))

# Plot B: Symmetric 2Theta-Theta Scan Comparison (Log Scale)
plt.figure(figsize=(10, 6))
plt.plot(arr_holder[:, 0], arr_holder[:, 1], '-', color='#7f7f7f', linewidth=2.0, label='Sample Holder Background')
plt.plot(arr_b3[:, 0], arr_b3[:, 1] + 3e5, '-', color='#2ca02c', linewidth=1.5, label='SH-124-B3 at Phi = 0° (+3e5 offset)')
plt.plot(arr_a[:, 0], arr_a[:, 1] + 6e5, '-', color='#1f77b4', linewidth=1.5, label='SH-125-A at Phi = 0° (+6e5 offset)')

plt.xlabel("2Theta (degrees)", fontsize=12)
plt.ylabel("Intensity (counts, stacked log scale)", fontsize=12)
plt.title("Symmetric 2Theta-Theta Scans: Background vs. Thin Films", fontsize=13, fontweight='bold')
plt.yscale('log')
plt.grid(True, which='both', linestyle='--', alpha=0.5)
plt.legend(fontsize=10, loc='upper right')
plt.tight_layout()

plot_2t_path = os.path.join(dst_analysis, "symmetric_scan_comparison.png")
plt.savefig(plot_2t_path, dpi=150)
plt.close()
print(f"Saved 2Theta scan comparison plot to {plot_2t_path}")

# Copy plots to artifacts
shutil.copy2(plot_shape_path, os.path.join(artifacts_dir, "rocking_curve_shape_comparison.png"))
shutil.copy2(plot_2t_path, os.path.join(artifacts_dir, "symmetric_scan_comparison.png"))
print("Copied comparison plots to artifacts.")

# 4. Read the original June 8 report text
print("Constructing integrated report Markdown...")
orig_report_md_path = os.path.join(src_analysis, "rocking_curve_phi_variation_report.md")
with open(orig_report_md_path, 'r') as f:
    orig_md = f.read()

# We will search for Section 5 (Fitted Rocking Curve Peak Metrics) and insert our new sections before it
search_str = "# 5. Fitted Rocking Curve Peak Metrics"

integrated_sections = r"""# 5. Calibration and Verification of the Substrate Background
To ensure that the residual peak components resolved in the rocking curves of the thin films (Section 3) are true orientation features of the CaCO3 crystallites rather than instrumental or substrate scatter artifacts, a baseline calibration study was performed on **June 9, 2026**.

## 5.1 Sample Holder symmetric scan
A symmetric $2\theta-\theta$ scan of the empty sample holder was collected in the $2\theta = 27^\circ - 35^\circ$ range. Fitting a 3rd-order polynomial to the raw intensity data shows a flat, slowly varying amorphous scattering profile with no sharp Bragg reflections. The residuals from this fit show a standard deviation of only $3837.7$ counts (less than $0.6\%$ of the total signal), confirming that the sample holder behaves purely as an amorphous background scatterer.

## 5.2 Stacked Comparison of Symmetric Scans
The plot below compares the empty sample holder profile with the symmetric phase scans of the two thin-film samples at $\phi = 0^\circ$. While the thin films exhibit sharp crystalline peaks (Calcite (104) at $2\theta \approx 29.4^\circ$ and Vaterite (110) at $2\theta \approx 32.8^\circ$), the sample holder background is smooth, verifying that the background subtraction methodology is physically sound:

![Symmetric 2Theta-Theta Scans Comparison: Background vs. Thin Films](symmetric_scan_comparison.png){width=80%}

\newpage

# 6. Mosaicity & Tilt State Comparison: Single Crystal vs. Thin Films
To benchmark the orientation features of our thin-film samples, a rocking curve of a **Calcite Single Crystal** was collected on **June 9, 2026** at the Calcite (104) Bragg condition ($2\theta = 29.3725^\circ$).

## 6.1 Fit Analysis of the Calcite Single Crystal Rocking Curve
The single crystal rocking curve was fitted using the same physical volume correction and baseline subtraction methodology as the thin-film scans. The peak metrics extracted from the Gaussian fit to the residual intensity are:
* **Peak Center ($\theta_0$):** $17.4550^\circ$
* **Out-of-Plane Tilt Angle ($\chi = \theta_0 - 2\theta_0/2$):** $+2.7688^\circ$
* **Peak FWHM:** $4.1325^\circ$
* **Net Peak Height:** $2.656 \times 10^6$ counts
* **Net Peak Area:** $1.168 \times 10^7$ counts$\cdot$deg
* **Area/Base Ratio:** $3528.6\%$

## 6.2 Rocking Curve Shape Comparison
The plot below shows the normalized net residual profiles of the Calcite single crystal compared to the thin films `SH-124-B3` and `SH-125-A` (each at the highly active $\phi = 60^\circ$ rotation) as a function of the tilt angle $\chi$:

![Rocking Curve Shape Comparison: Mosaicity vs. Thin-Film Texture](rocking_curve_shape_comparison.png){width=80%}

## 6.3 Scientific Discussion & Insights
1. **Correlation of Tilt States:**
   The calcite single crystal is tilted relative to the specular condition by $\chi = +2.7688^\circ$. This tilt state aligns remarkably well with the minor tilt peaks resolved in the thin films:
   * **`SH-124-B3`** has minor tilt components at $\chi \approx +2.51^\circ$ (Minor Peak A) and $\chi \approx +2.93^\circ$ (Minor Peak B).
   * **`SH-125-A`** has similar components at $\chi \approx +2.49^\circ$ (Minor Peak A) and $\chi \approx +2.89^\circ$ (Minor Peak B).
   This close correlation suggests that the single crystal mounting or cleavage plane possesses a miscut tilt that matches the specific epitaxial template orientation pathways preferred by the thin-film calcite crystallites during nucleation and growth.
   
2. **Bragg Peak Broadening (Mosaicity vs. Film Grain Spread):**
   Surprisingly, the FWHM of the single crystal Bragg peak ($4.1325^\circ$) is much broader than the individual tilt peaks of the thin-film samples ($0.15^\circ - 0.45^\circ$). 
   While thin-film samples represent a polycrystalline ensemble where individual grain populations are co-oriented with a very narrow tilt spread, the single crystal under study exhibits a broad profile. This indicates a high density of defects, low-angle mosaic grain boundaries, or lattice curvature within the crystal, combined with experimental footprint broadening.

\newpage
"""

# Insert the new sections before Section 5
parts = orig_md.split(search_str)
if len(parts) == 2:
    integrated_md = parts[0] + integrated_sections + search_str + parts[1]
    
    # Update title in YAML header of the MD content to indicate Integrated Report
    integrated_md = integrated_md.replace(
        'title: "CaCO3 Thin Film Rocking Curve Analysis: Azimuthal (Phi) Variation Study"',
        'title: "CaCO3 Thin Film Rocking Curve Analysis: Integrated Study with Single Crystal & Background Calibration"'
    )
    
    # Save the new MD report
    dst_md_path = os.path.join(dst_analysis, "rocking_curve_phi_variation_report_integrated.md")
    with open(dst_md_path, 'w') as f:
        f.write(integrated_md)
    print(f"Written integrated MD report to {dst_md_path}")
    
    # Compile to PDF using pandoc
    dst_pdf_path = os.path.join(dst_analysis, "rocking_curve_phi_variation_report_integrated.pdf")
    cmd = [
        "pandoc",
        "rocking_curve_phi_variation_report_integrated.md",
        "-o", "rocking_curve_phi_variation_report_integrated.pdf",
        "--pdf-engine=pdflatex"
    ]
    print("Compiling integrated report to PDF...")
    result = subprocess.run(cmd, cwd=dst_analysis, capture_output=True, text=True)
    if result.returncode == 0:
        print(f"Successfully compiled integrated PDF report at {dst_pdf_path}")
        
        # Copy to artifacts
        art_pdf_path = os.path.join(artifacts_dir, "rocking_curve_phi_variation_report_integrated.pdf")
        shutil.copy2(dst_pdf_path, art_pdf_path)
        print(f"Copied integrated PDF report to artifacts: {art_pdf_path}")
    else:
        print("Error compiling report to PDF:")
        print("STDOUT:", result.stdout)
        print("STDERR:", result.stderr)
else:
    print("Error: Could not locate Section 5 in the original report markdown!")
