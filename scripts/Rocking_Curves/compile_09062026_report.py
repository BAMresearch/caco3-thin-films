# -*- coding: utf-8 -*-
"""
Report Generator for 09 June 2026 Run
=====================================
Author: Tomasz Stawski, tomasz.stawski@bam.de
Date: 2026-06-19
Version: 0.1

Description:
Formulates tables and compiles the technical summary report for rocking curve measurements from 09 June 2026.
"""
import os
import subprocess
import shutil

# Paths
base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
analysis_dir = os.path.join(base_dir, "data/processed/Rocking_Curves/Reference")
artifacts_dir = "/home/tomek/.gemini/antigravity/brain/297360fa-0f2b-4ba6-96f5-6c818b856b29/artifacts"

os.makedirs(analysis_dir, exist_ok=True)
os.makedirs(artifacts_dir, exist_ok=True)

md_content = r"""---
title: "Calcite Single Crystal and Sample Holder Background results"
date: "June 9, 2026"
geometry: "margin=1in"
fontsize: "11pt"
header-includes:
  - \usepackage{fancyhdr}
  - \pagestyle{fancy}
  - \fancyhead[CO,CE]{Calcite Single Crystal and Background Study}
  - \fancyfoot[CO,CE]{Page \thepage}
---

# 1. Executive Summary
This report presents the analysis of the new datasets collected on **June 9, 2026**:
1. **Calcite Single Crystal Rocking Curve** (`Calcite single crystal rocking curve.brml`)
2. **Sample Holder 2Theta-Theta Background Scan** (`sample holder.brml`)

Using the same baseline subtraction and fitting methodology developed for thin-film samples, we analyze the single crystal rocking curve to determine its tilt state and peak broadening. Additionally, we evaluate the sample holder scan to confirm the absence of background reflections in the region of interest.

---

# 2. Dataset Structure
* **Calcite Rocking Curve:** 401 points of raw intensity vs. $\theta$ (ranging from $4.7218^\circ$ to $24.7218^\circ$ in steps of $0.05^\circ$). The nominal 2Theta detector position was fixed at $29.3725^\circ$.
* **Sample Holder Scan:** 401 points of intensity vs. $2\theta$ (ranging from $27.0001^\circ$ to $35.0001^\circ$ in steps of $0.02^\circ$). The $\theta$ axis tracks $2\theta/2$ (from $13.50005^\circ$ to $17.50005^\circ$).

---

# 3. Calcite (104) Single Crystal Rocking Curve Analysis
## 3.1 Background Subtraction & Peak Fitting
The rocking curve was modeled using our thin-film background model (isotropic footprint correction + 3rd-order polynomial). The peak region ($10.0^\circ \le \theta \le 22.0^\circ$) was masked during baseline fitting. The background model is given by:
$$I_{\text{baseline}}(\theta) = \frac{I_0}{\sin(\theta)} + c_0 + c_1\theta + c_2\theta^2 + c_3\theta^3$$

After baseline subtraction, the net residual intensity was fitted with a Gaussian profile to extract the peak metrics:

| Parameter | Value |
| :--- | :---: |
| **Peak Center ($\theta$)** | $17.455^\circ$ |
| **Nominal Specular Center ($2\theta_0 / 2$)** | $14.686^\circ$ |
| **Tilt Angle ($\chi = \theta - 2\theta_0/2$)** | $+2.769^\circ$ |
| **Peak FWHM** | $4.133^\circ$ |
| **Net Peak Height** | $2.66 \times 10^6$ counts |
| **Net Peak Area** | $1.17 \times 10^7$ counts$\cdot$deg |
| **Isotropic Base at Peak** | $3.31 \times 10^5$ counts |
| **Area/Base Ratio** | $3528.6\%$ |

## 3.2 Rocking Curve Analysis Plot
The plot below displays the experimental data, the model baseline, the background-corrected net intensity, and the Gaussian peak fit.

![Calcite Single Crystal Rocking Curve Analysis](calcite_single_crystal_rocking_curve_analysis.png){width=85%}

\newpage

# 4. Sample Holder 2Theta-Theta Background Scan
To verify the contribution of the sample holder to the measurements, we analyzed its symmetric scan in the $2\theta$ range of $27.0^\circ - 35.0^\circ$. A 3rd-order polynomial was fitted to the data to model the background scatter:
$$I_{\text{baseline}}(2\theta) = p_3(2\theta)^3 + p_2(2\theta)^2 + p_1(2\theta) + p_0$$

The raw data and residuals are shown in the plot below. The residuals show a standard deviation of only $3837$ counts (less than $0.6\%$ of the average signal of $6.84 \times 10^5$ counts) and contain no sharp Bragg reflections. This confirms that the sample holder behaves as a featureless amorphous scatterer and does not introduce spurious crystalline peaks.

![Sample Holder 2Theta-Theta Background Scan](sample_holder_2theta_analysis.png){width=85%}

---

# 5. Scientific Discussion and Insights
1. **Tilt State Comparison with Thin Films:**
   The calcite single crystal exhibits a clear tilt angle of $\chi = +2.769^\circ$. This position matches very closely with the minor tilt states observed in our thin film samples:
   * **`SH-124-B3`** showed a minor tilt component at $\chi \approx +2.51^\circ$ (Minor Peak A) and $\chi \approx +2.93^\circ$ (Minor Peak B).
   * **`SH-125-A`** showed similar components at $\chi \approx +2.49^\circ$ (Minor Peak A) and $\chi \approx +2.89^\circ$ (Minor Peak B).
   The alignment of the single crystal Bragg peak near $\chi \approx +2.77^\circ$ indicates that the single crystal mounting or its cleavage plane possesses a miscut or tilt that mimics the specific growth alignment planes preferred by the columnar grains in the thin films.
   
2. **Bragg Peak Width (Mosaicity):**
   The FWHM of the single crystal peak is $4.133^\circ$, which is significantly broader than typical values for high-quality single crystals (which are usually $< 0.1^\circ$). This broad profile suggests a high density of defects, multiple mosaic blocks (mosaicity), or a significant curvature in the crystal lattice. Alternatively, it could reflect experimental broadening from a wide beam footprint or divergent slit configurations.
"""

md_path = os.path.join(analysis_dir, "calcite_crystal_and_holder_report.md")
with open(md_path, 'w') as f:
    f.write(md_content)
print(f"Written report Markdown to {md_path}")

pdf_path = os.path.join(analysis_dir, "calcite_crystal_and_holder_report.pdf")
cmd = [
    "pandoc",
    "calcite_crystal_and_holder_report.md",
    "-o", "calcite_crystal_and_holder_report.pdf",
    "--pdf-engine=pdflatex"
]
print("Compiling report to PDF...")
result = subprocess.run(cmd, cwd=analysis_dir, capture_output=True, text=True)
if result.returncode == 0:
    print(f"Successfully generated PDF report at {pdf_path}")
    
    # Copy PDF to artifacts
    art_pdf_path = os.path.join(artifacts_dir, "calcite_crystal_and_holder_report.pdf")
    shutil.copy2(pdf_path, art_pdf_path)
    print(f"Copied PDF report to artifacts: {art_pdf_path}")
    
    # Also copy png plots to artifacts so they can be embedded in model response/artifacts
    shutil.copy2(os.path.join(analysis_dir, "calcite_single_crystal_rocking_curve_analysis.png"), 
                 os.path.join(artifacts_dir, "calcite_single_crystal_rocking_curve_analysis.png"))
    shutil.copy2(os.path.join(analysis_dir, "sample_holder_2theta_analysis.png"), 
                 os.path.join(artifacts_dir, "sample_holder_2theta_analysis.png"))
    print("Copied plots to artifacts.")
else:
    print("Error compiling report to PDF:")
    print("STDOUT:", result.stdout)
    print("STDERR:", result.stderr)
