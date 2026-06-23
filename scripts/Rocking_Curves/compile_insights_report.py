# -*- coding: utf-8 -*-
"""
Insights Report Compiler
========================
Author: Tomasz Stawski, tomasz.stawski@bam.de
Date: 2026-06-19
Version: 0.1

Description:
Formulates observations and generates summaries on growth dynamics from processed diffraction metrics.
"""
import os
import subprocess

base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
analysis_dir = os.path.join(base_dir, "data/processed/Rocking_Curves/08062026")
artifacts_dir = os.path.join(base_dir, "results/reports")

md_content = """---
title: "Structural and Epitaxial Insights: Azimuthal Rocking Curve Study of $CaCO_3$ Thin Films"
date: "June 9, 2026"
geometry: "margin=1in"
fontsize: "11pt"
header-includes:
  - \\usepackage{fancyhdr}
  - \\pagestyle{fancy}
  - \\fancyhead[CO,CE]{$CaCO_3$ Azimuthal Rocking Curve Insights}
  - \\fancyfoot[CO,CE]{Page \\thepage}
---

# 1. Executive Summary
This report analyzes the domain orientation distribution and template-guided alignment of Calcite (104) grains in samples **SH-124-B3** (pure calcite) and **SH-125-A** (mixed calcite/vaterite). By comparing the raw rocking curves against the baseline-corrected (net intensity) curves across different azimuthal angles ($\\phi$), we distinguish between continuous mosaic textures and discrete single-crystal-like epitaxial domains.

---

# 2. Sample SH-124-B3: Discrete Epitaxial Domains
The rocking curve behavior of sample **SH-124-B3** is characterized by highly localized, sharp reflections and the absence of a broad, continuous background envelope:
1. **Negligible Isotropic Matrix:** The baseline-corrected net rocking curves return to zero across most of the angular range, indicating that the film contains a negligible fraction of randomly oriented (isotropic) calcite crystallites.
2. **Discrete Crystalline Islands:** The diffraction patterns are dominated by sharp, isolated peaks at highly specific tilt and rotation coordinates (e.g., the distinct tilt reflection at $\\theta \\approx 9.15^\\circ$ and $13.0^\\circ$ for $\\phi = 30^\\circ$, and the strong tilt peak shifting down to $\\theta \\approx 7.94^\\circ$ at $\\phi = 180^\\circ$).
3. **Low Mosaicity:** The FWHM of these individual reflections is narrow ($\\approx 0.25^\\circ - 0.30^\\circ$), confirming high crystallographic alignment within these isolated domains. This indicates growth consisting of discrete, single-crystal-like epitaxial columns/islands template-matched to the substrate lattice.

![SH-124-B3 Rocking Curve Analysis: Raw Stacked vs. Baseline-Corrected Net Stacked](SH-124-B3/SH-124-B3_side_by_side.png){width=100%}

\\newpage

# 3. Sample SH-125-A: Mosaic Texture and Orthogonal Symmetry
Sample **SH-125-A** exhibits a combination of a broad orientation distribution (mosaic spread) and highly co-oriented, low-mosaicity tilted domains:
1. **Broad Specular Envelope with Fourfold Symmetry:** A broad specular maximum (FWHM $\\ge 0.6^\\circ$) is visible near $\\theta \\approx 15^\\circ$ exclusively for the $\\phi = 60^\\circ$ and $\\phi = 150^\\circ$ curves. The exact $90^\\circ$ azimuthal separation between these active scans ($\\Delta\\phi = 90^\\circ$) indicates a two-fold or four-fold orthogonal in-plane symmetry of this near-specular grain population, dictated by substrate-guided epitaxy.
2. **Coexisting Low-Mosaicity Tilted Domains:** Superimposed on the broad background are sharp, narrow reflections (such as `Peak 2a` at $\\theta \\approx 11.99^\\circ$, FWHM $\\approx 0.26^\\circ$ at $\\phi = 60^\\circ$). These represent highly oriented, low-mosaicity tilted columnar domains (e.g. helical columns) growing along specific template-guided directions.

![SH-125-A Rocking Curve Analysis: Raw Stacked vs. Baseline-Corrected Net Stacked](SH-125-A/SH-125-A_side_by_side.png){width=100%}

---

# 4. Summary of Key Findings
* **SH-124-B3** displays a growth mode consisting of isolated, highly aligned epitaxial grains (low mosaicity, no isotropic matrix) with discrete tilt/azimuth orientations.
* **SH-125-A** represents a mosaic film with a broad, near-specular orientation distribution exhibiting orthogonal in-plane symmetry, coexisting with highly ordered, narrow off-axis tilted domains.
"""

# Write md file
md_path = os.path.join(analysis_dir, "rocking_curve_insights_report.md")
with open(md_path, 'w') as f:
    f.write(md_content)
print(f"Written insights report Markdown to {md_path}")

# Compile to PDF
pdf_path = os.path.join(analysis_dir, "rocking_curve_insights_report.pdf")
cmd = [
    "pandoc",
    "rocking_curve_insights_report.md",
    "-o", "rocking_curve_insights_report.pdf",
    "--pdf-engine=pdflatex"
]
print("Compiling insights report to PDF...")
result = subprocess.run(cmd, cwd=analysis_dir, capture_output=True, text=True)
if result.returncode == 0:
    print(f"Successfully generated PDF report at {pdf_path}")
    # Copy to artifacts
    import shutil
    art_path = os.path.join(artifacts_dir, "rocking_curve_insights_report.pdf")
    shutil.copy2(pdf_path, art_path)
    print(f"Copied PDF report to artifacts: {art_path}")
else:
    print("Error compiling report to PDF:")
    print("STDOUT:", result.stdout)
    print("STDERR:", result.stderr)
