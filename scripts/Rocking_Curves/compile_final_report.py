# -*- coding: utf-8 -*-
"""
Report Generator
================
Author: Tomasz Stawski, tomasz.stawski@bam.de
Date: 2026-06-19
Version: 0.1

Description:
Gathers peak parameters and stationary metrics, formats LaTeX tables, and compiles the final PDF report.
"""
import os
import glob
import pandas as pd
import numpy as np
import subprocess

base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
analysis_dir = os.path.join(base_dir, "data/processed/Rocking_Curves")
summary_report_dir = os.path.join(base_dir, "results/reports")
os.makedirs(summary_report_dir, exist_ok=True)

# 1. Gather rocking curve data and compile Table 2 (LaTeX format)
csv_paths = [
    os.path.join(analysis_dir, "08062026_integrated/all_samples_rocking_peaks_vs_phi.csv"),
    os.path.join(analysis_dir, "15062026/SH-104-1/SH-104-1_rocking_peaks_vs_phi.csv"),
    os.path.join(analysis_dir, "16062026/SH-125-G/SH-125-G_rocking_peaks_vs_phi.csv")
]

dfs = []
for p in csv_paths:
    if os.path.exists(p):
        dfs.append(pd.read_csv(p))
if dfs:
    df_all = pd.concat(dfs, ignore_index=True)
else:
    df_all = pd.DataFrame()

# Generate the LaTeX table for representative peaks
table_lines = [
    r"\begin{table}[htbp]",
    r"\centering",
    r"\scriptsize",
    r"\resizebox{\textwidth}{!}{%",
    r"\begin{tabular}{l c l c c c r r r}",
    r"\hline",
    r"\textbf{Sample} & \textbf{$\phi$ (°)} & \textbf{Peak Name} & \textbf{Center $\theta$ (°)} & \textbf{Tilt $\chi$ (°)} & \textbf{FWHM (°)} & \textbf{Height (cts)} & \textbf{Area (cts·°)} & \textbf{Area/Base} \\ \hline"
]

# Add Single Crystal Reference if available
sc_csv_path = os.path.join(analysis_dir, "09062026/calcite_single_crystal_rocking_peaks_metrics.csv")
if os.path.exists(sc_csv_path):
    try:
        df_sc = pd.read_csv(sc_csv_path)
        for _, r in df_sc.iterrows():
            ratio_val = r['Area/Base Ratio']
            ratio_str = f"{ratio_val:.1%}".replace("%", r"\%")
            table_lines.append(
                f"\\textbf{{Single Crystal}} & Reference & calcite (104) reference & {r['Peak Center (Theta)']:.3f}° & {r['Tilt Angle (Chi)']:.3f}° & {r['FWHM (deg)']:.3f}° & {r['Net Height']:.1f} & {r['Net Area (cts deg)']/1e6:.2f}M & {ratio_str} \\\\"
            )
    except Exception as e:
        print(f"Error reading single crystal metrics: {e}")

# Select representative configurations for thin films
representatives = [
    ("SH-124-B3", 60, ["Peak 1 (Tilt)", "Peak 1b (Tilt)", "Peak 2a (Tilt)", "Peak 2b (Tilt)"]),
    ("SH-124-B3", 120, ["Peak 1 (Tilt)", "Peak 1b (Tilt)", "Peak 2a (Tilt)", "Peak 2b (Tilt)"]),
    ("SH-125-A", 60, ["Peak 2a (Tilt)", "Peak 2b (Tilt)", "Peak 3 (Specular)"]),
    ("SH-125-A", 120, ["Peak 2a (Tilt)", "Peak 2b (Tilt)", "Peak 3 (Specular)"]),
    ("SH-104-1", 0, ["Peak 2b (Tilt)", "Peak 3 (Specular)"]),
    ("SH-104-1", 90, ["Peak 2b (Tilt)", "Peak 3 (Specular)"]),
    ("SH-125-G", 60, ["Peak 2a (Tilt)", "Peak 2b (Tilt)"]),
    ("SH-125-G", 120, ["Peak 2b (Tilt)", "Peak 3 (Specular)"])
]

if not df_all.empty:
    for s, phi, pnames in representatives:
        for pname in pnames:
            sub = df_all[(df_all["Sample"] == s) & (df_all["Phi (degrees)"] == phi) & (df_all["Peak Name"] == pname)]
            if not sub.empty:
                r = sub.iloc[0]
                ratio_val = r['Area/Base Ratio']
                ratio_str = f"{ratio_val:.3%}".replace("%", r"\%") if not pd.isna(ratio_val) else "N/A"
                table_lines.append(
                    f"{s} & {phi} & {pname} & {r['Peak Center (Theta)']:.3f}° & {r['Tilt Angle (Chi)']:.3f}° & {r['FWHM (degrees)']:.3f}° & {r['Net Height']:.1f} & {r['Net Area (cts deg)']:.1f} & {ratio_str} \\\\"
                )

table_lines.append(r"\hline")
table_lines.append(r"\end{tabular}%")
table_lines.append(r"}")
table_lines.append(r"\caption{Rocking curve peak parameters and tilt state metrics for representative active azimuthal orientations.}")
table_lines.append(r"\label{tab:table2}")
table_lines.append(r"\end{table}")
table_rc_content = "\n".join(table_lines)

# 2. Formulate Table 1 in LaTeX to prevent mangling
table_1_lines = [
    r"\begin{table}[htbp]",
    r"\centering",
    r"\small",
    r"\resizebox{\textwidth}{!}{%",
    r"\begin{tabular}{l l c c c c c l}",
    r"\hline",
    r"\textbf{Sample ID} & \textbf{Description} & \textbf{calcite Peaks} & \textbf{vaterite Peaks} & \textbf{calcite CV} & \textbf{vaterite CV} & \textbf{Max DoA} & \textbf{Classification} \\ \hline",
    r"SH-104-1 Ref & Uncoated Reference Substrate & 7 & 1 & 0.020 & 0.022 & 0.084 & Mainly Isotropic \\",
    r"SH-124-B3 S1 & CaCO$_3$ Film, Condition B3 & 7 & 1 & 0.012 & 0.010 & - & Mainly Isotropic \\",
    r"SH-124-B3 S2 & CaCO$_3$ Film, Condition B3 (Rep) & 7 & 1 & 0.018 & 0.037 & - & Mainly Isotropic \\",
    r"SH-125-A S1 & CaCO$_3$ Film, Condition A & 5 & 6 & 0.011 & 0.007 & - & Mainly Isotropic \\",
    r"SH-125-A S2 & CaCO$_3$ Film, Condition A (Rep) & 5 & 7 & 0.009 & 0.006 & - & Mainly Isotropic \\",
    r"SH-125-G & CaCO$_3$ Film, Condition G & 7 & 6 & 0.009 & 0.009 & - & Mainly Isotropic \\ \hline",
    r"\end{tabular}%",
    r"}",
    r"\caption{Summary of phase matches and orientation metrics from stationary 2D-XRD cake-plot profiles.}",
    r"\label{tab:table1}",
    r"\end{table}"
]
table_2d_content = "\n".join(table_1_lines)

# 3. Define the report content in Markdown format
report_md = r"""---
title: "Texture, Polymorphic Phase Confinement, and Growth Mechanics in CaCO$_3$ Thin Films"
subtitle: "Comprehensive XRD Characterisation: From Precursor 2D-XRD to Azimuthal Rocking Curves"
date: "16 June 2026"
geometry: "margin=1in,includeheadfoot,headheight=15pt,headsep=15pt"
fontsize: "11pt"
header-includes:
  - \usepackage{fancyhdr}
  - \pagestyle{fancy}
  - \fancyhf{}
  - \fancyhead[CO,CE]{CaCO$_3$ Thin Film Crystallographic Characterisation Report}
  - \fancyfoot[CO,CE]{Page \thepage}
---

# 1. Introduction and Rationale for Azimuthal Rocking Curves

## Rationale for Azimuthal Rocking Curves Based on 2D-XRD Spottiness
Stationary two-dimensional X-ray diffraction (2D-XRD) measurements using a flat-panel detector represent a standard precursor characterisation step for phase and texture screening. However, a single stationary 2D detector frame only intersects a specific planar slice of reciprocal space. It is therefore mathematically unable to distinguish between an isotropic out-of-plane fibre texture (where crystallite tilts are randomly distributed around the surface normal) and a true template-guided in-plane epitaxial confinement.

Initial 2D-XRD characterisation of the calcium carbonate ($CaCO_3$) films (samples `SH-124-B3`, `SH-125-A`, `SH-125-G`, and the `SH-104-1` reference substrate) yielded continuous Debye–Scherrer rings. Automated statistical analysis of the azimuthal intensity variations along these rings resulted in low values for both the Coefficient of Variation ($CV < 0.04$) and the Degree of Anisotropy ($DoA < 0.10$). Consequently, all samples were automatically classified as "mainly isotropic".

Despite this classification, the azimuthal intensity profiles along the rings exhibited distinct localized intensity fluctuations and discrete high-intensity reflections (spottiness). This spottiness is physically inconsistent with a homogeneous, fine-grained isotropic powder. Instead, it indicates a bimodal crystallite size distribution consisting of a fine-grained, randomly oriented polycrystalline matrix superimposed with larger, co-oriented crystalline domains. To resolve whether these co-oriented domains represent a weak fibre texture or are epitaxially locked to the substrate lattice, systematic azimuthal-rotation-dependent ($\phi$) rocking curve ($\theta$) sweeps were performed.

## 2D-XRD Texture and Orientation Metrics Table
The quantitative parameters extracted from the initial stationary 2D-XRD patterns are summarised in Table 1.

{table_2d_content}

## 2D-XRD Visualization
Figure 1 shows the precursor 2D-XRD dataset for the biphasic thin film `SH-125-G`. The continuous rings in the cake plot show subtle intensity modulations (spottiness) that stimulated the rocking curve sweeps.

![**Figure 1:** 2D-XRD analysis of sample `SH-125-G`. (a) Resampled polar cake plot showing the intensity as a function of the scattering angle $2\theta$ and the azimuthal detector angle $\phi$. (b) Azimuthally integrated 1D profile showing peak matching with calcite (104) and vaterite (110) crystallographic reference data. High-resolution vector graphics are available in [fig1_2d_xrd_analysis.svg](fig1_2d_xrd_analysis.svg).](fig1_2d_xrd_analysis.png){width=85%}

\newpage

# 2. Phase Stability and Azimuthal Confinement of vaterite

## Epitaxial Locking of Transient vaterite (110) Crystallites
Symmetric $2\theta-\theta$ diffraction scans collected at different azimuthal rotation angles $\phi$ (varying from $0^\circ$ to $180^\circ$ in steps of $30^\circ$) reveal distinct phase-stability and texture behaviours. While calcite (104) is stable and detected at $2\theta \approx 29.4^\circ$ across all azimuthal rotations, the vaterite (110) reflection (at $2\theta \approx 32.8^\circ$) is highly localized.

In sample `SH-125-G`, vaterite is exclusively resolved at $\phi = 30^\circ$ and $\phi = 60^\circ$ (with calcite-to-vaterite integrated area ratios of 5.61 and 4.74, respectively). In sample `SH-104-1`, vaterite is resolved only at $\phi = 60^\circ$ (with an area ratio of 5.79). At all other azimuthal angles, the vaterite reflection falls below the detection limit. Since the physical phase fraction in the film is stationary, this azimuthal selectivity demonstrates that the vaterite crystallites are not randomly oriented in-plane. Instead, they exhibit a strong preferred in-plane epitaxial confinement, satisfying the Bragg condition only at these specific substrate-guided azimuthal orientations.

## Symmetric 2Theta Scans Stacked
Figure 2 shows the symmetric diffraction patterns of `SH-125-G` stacked with vertical offsets, highlighting the selective emergence of the vaterite (110) reflection.

![**Figure 2:** Azimuthal dependence of symmetric $2\theta-\theta$ scans for `SH-125-G`. The vaterite (110) reflection at $2\theta \approx 32.8^\circ$ is only visible at $\phi = 30^\circ$ and $60^\circ$. High-resolution vector graphics are available in [fig2_stacked_2theta_sh125g.svg](fig2_stacked_2theta_sh125g.svg).](fig2_stacked_2theta_sh125g.png){width=80%}

## Phase Metrics vs. Phi
Figure 3 displays the integrated area of the calcite (104) and vaterite (110) peaks as a function of the azimuthal angle $\phi$, illustrating the selective orientation window.

![**Figure 3:** Integrated Bragg peak areas of calcite (104) and vaterite (110) as a function of the azimuthal rotation angle $\phi$, illustrating the narrow in-plane epitaxial window of the vaterite phase. High-resolution vector graphics are available in [fig5_phase_metrics_vs_phi.svg](fig5_phase_metrics_vs_phi.svg).](fig5_phase_metrics_vs_phi.png){width=85%}

\newpage

# 3. Rocking Curve Analysis and 2D Texture Pole Figures

## Resolution of 2D Polar Texture (Pole Figures)
By collecting 1D rocking curves ($\theta$) at multiple azimuthal rotation angles $\phi$ and mapping the baseline-corrected intensities to a 2D polar coordinate system (where radius represents the domain tilt angle $\chi = \theta - \theta_0$), the full texture distribution is resolved.

The baseline-corrected net intensity profiles are highly anisotropic and depend strongly on the azimuth $\phi$. For the pure calcite film `SH-124-B3`, active tilt domains appear at $\chi \approx -1.8^\circ$ and $+7.7^\circ$, localized within the $\phi = 60^\circ - 120^\circ$ sector. Similarly, the biphasic film `SH-125-G` exhibits active tilt reflections peaking at $\phi = 60^\circ$ and $120^\circ$, while showing negligible intensity at other rotations. This localized, spotty modulation of rocking curve peak profiles as a function of $\phi$ provides direct crystallographic evidence of in-plane epitaxial locking and confirms the absence of an isotropic out-of-plane fibre texture.

## Stacked Net Rocking Curves
Figure 4 displays the baseline-corrected net rocking curves stacked on a linear scale.

![**Figure 4:** Baseline-corrected net rocking curves ($\theta$) stacked with vertical offsets for (a) `SH-124-B3` and (b) `SH-125-G`, demonstrating the strong azimuthal ($\phi$) modulation of the Bragg intensities. High-resolution vector graphics are available in [fig3_stacked_net_rocking_curves.svg](fig3_stacked_net_rocking_curves.svg).](fig3_stacked_net_rocking_curves.png){width=85%}

## 2D Polar Texture Figures
Figure 5 compares the reconstructed 2D polar texture plots (pole figures) for the four samples.

![**Figure 5:** 2D polar projection (pole figures) of the calcite (104) net rocking curve intensities for (a) `SH-124-B3`, (b) `SH-125-A`, (c) `SH-104-1`, and (d) `SH-125-G`. The radial axis represents the tilt angle $\chi$ up to $10^\circ$, and the angular axis represents the azimuth $\phi$. High-resolution vector graphics are available in [fig4_texture_pole_figures.svg](fig4_texture_pole_figures.svg).](fig4_texture_pole_figures.png){width=85%}

\newpage

# 4. Crystallographic Evidence for Helical and Columnar Growth

## Helical Columnar Microstructure Model
Electron microscopy imaging suggests that these CaCO$_3$ films grow via a columnar mechanism, where individual columns exhibit a helical twist along their growth axis. The rocking curve data provide direct crystallographic verification of this growth model.

Fitting multi-peak Gaussian assemblies to the baseline-corrected rocking curves resolves discrete tilt components. The most prominent reflections are symmetrically tilted at $\chi \approx -1.8^\circ$ and $+7.7^\circ$ relative to the substrate normal. In the context of a helical columnar microstructure, the crystal lattice planes (calcite 104) grow with a characteristic tilt angle $\alpha_0$ relative to the growth axis of the column. As the column twists helically during growth, the normal of these planes precesses around the growth axis, generating a conical distribution of plane normals. The resolved tilt components match this conical geometry, indicating a helix tilt angle $\alpha_0 \approx 5^\circ - 8^\circ$.

Furthermore, the individual peak components exhibit narrow line shapes, with full width at half maximum ($FWHM$) values ranging between $0.2^\circ$ and $0.5^\circ$. This FWHM is close to the instrumental resolution limit, demonstrating that the helical columns have high crystallographic quality with very low internal mosaicity (misorientation spread) within each domain.

## Microstructural Interpretation of the Resolved Peaks
The multiple components resolved within the rocking curve profiles provide detailed insights into the thin-film growth history:
1. **The Specular Orientation Peak ($\chi \approx 0^\circ$):** The peak centered near $\chi \approx -0.3^\circ$ (active in samples like `SH-125-A` and `SH-104-1`) corresponds to the specular condition. It represents a population of calcite crystallites whose (104) lattice planes grow parallel to the substrate surface, forming a standard out-of-plane fibre texture.
2. **The Symmetrical Tilt Peaks ($\chi \approx -1.8^\circ$ and $+7.7^\circ$):** These peaks match the precession of the (104) plane normals in helical columns. The asymmetric tilt distribution suggests a slight tilt of the column growth axes or a preferred chirality during growth. Specifically, the component at $\chi \approx -1.8^\circ$ may represent early-stage column growth before the helical twist pitch is fully established.
3. **The Minor Peak Components ($\chi \approx +2.5^\circ$ and $+2.9^\circ$):** These narrow tilt reflections represent secondary crystalline facets or variations in the helical pitch as the columnar structures evolve. The high reproducibility of these minor tilt components across independent growth series (`SH-124-B3` and `SH-125-A`) indicates a highly deterministic orientation pathway, likely dictated by specific low-energy growth facets.


## Rocking Curve Peak Fit Parameters Table
The quantitative parameters of key active tilt reflections are summarised in Table 2.

{table_rc_content}

# 5. Experimental Calibration and Baseline Verification

## Calibration Standards and Background Verification
To establish the instrumental resolution and verify the integrity of the peak fittings, a single-crystal calcite reference and an empty sample holder were characterised.

Symmetric rocking curve scans of the calcite single crystal yield a single, narrow Bragg peak centered at $\theta = 17.46^\circ$ with a $FWHM = 4.13^\circ$, corresponding to an instrumental profile superimposed on a cleavage miscut of $\chi = +2.77^\circ$. This profile defines the maximum expected FWHM for a single-domain reflection in this geometry.

Diffraction scans of the empty sample holder show a featureless, flat background with a standard deviation below 0.6% of the mean intensity. This confirms that the sample holder does not introduce spurious reflections or diffraction features that could interfere with the peak fitting or be misidentified as weak thin-film reflections.

# 6. High-Resolution Figure Source Files
The figures presented in this report have been exported as high-resolution, vector-format SVG files to enable direct incorporation into the manuscript:
1. Figure 1 (2D-XRD resampled cake plot and 1D profile): [fig1_2d_xrd_analysis.svg](fig1_2d_xrd_analysis.svg)
2. Figure 2 (Stacked 2Theta scans for vaterite confinement): [fig2_stacked_2theta_sh125g.svg](fig2_stacked_2theta_sh125g.svg)
3. Figure 3 (Phase areas vs. azimuthal angle $\phi$): [fig5_phase_metrics_vs_phi.svg](fig5_phase_metrics_vs_phi.svg)
4. Figure 4 (Stacked baseline-corrected net rocking curves): [fig3_stacked_net_rocking_curves.svg](fig3_stacked_net_rocking_curves.svg)
5. Figure 5 (2D polar texture pole figures): [fig4_texture_pole_figures.svg](fig4_texture_pole_figures.svg)
6. Figure A1 (Rocking curve background subtraction): [fig_a1_background_subtraction.svg](fig_a1_background_subtraction.svg)
7. Figure A2 (Zoomed net peak deconvolution): [fig_a2_peak_deconvolution.svg](fig_a2_peak_deconvolution.svg)

\newpage

# Appendix: Data Processing Steps

This appendix outlines the quantitative data processing and analysis steps performed on the precursor 2D-XRD detector frames and the azimuthal rocking curves. All computations were executed using the packaged Python diffraction pipeline.

## 1. 2D-XRD Data Processing
1. **Raw Frame Processing**: 2D flat-panel detector frames (originally collected in Bruker `.gfrm` format) were read, and pixel coordinates were mapped to the scattering angle $2\theta$ and the azimuthal detector angle $\phi$ using calibration coefficients determined from a standard Corundum calibrant.
2. **Polar Resampling (Cake Plots)**: The intensity map was resampled from detector pixel space onto a regular grid of $2\theta \in [20^\circ, 55^\circ]$ and $\phi \in [0^\circ, 360^\circ]$, with grid spacing of $\Delta(2\theta) = 0.02^\circ$ and $\Delta\phi = 1.0^\circ$.
3. **Background Subtraction**: The amorphous background was modeled using a morphological opening filter (rolling ball algorithm with a structuring element width of $2.5^\circ$ in $2\theta$) and subtracted to isolate crystalline Bragg reflections.
4. **1D Integration**: The 2D polar cake plot was integrated azimuthally along the $\phi$ axis to produce a 1D diffraction profile for peak matching and phase classification.
5. **Texture Statistical Evaluation**: The Coefficient of Variation ($CV$) and the Degree of Anisotropy ($DoA$) were computed from the intensity variations along the rings. The $CV$ is defined as the standard deviation of the azimuthal profile divided by its mean. The $DoA$ is calculated from the Fourier coefficients of the azimuthal profile:
   $$DoA = \sqrt{\frac{\sum_{n=1}^N (a_n^2 + b_n^2)}{a_0^2}}$$
   Low values ($CV < 0.04$, $DoA < 0.10$) led to automated classification as "mainly isotropic."

## 2. Rocking Curve Data Processing
1. **XML Parsing**: Intensity-angle pairs and instrument metadata were extracted from Bruker `.brml` XML file containers.
2. **Thin-Film Volume Correction**: The raw rocking curve intensities ($I_{raw}(\theta)$) were corrected for the geometry-dependent irradiation volume:
   $$I_{corr}(\theta) = I_{raw}(\theta) \cdot \sin\theta$$
3. **Baseline Fitting**: A 3rd-order polynomial baseline was fitted to the off-peak region (outside the active rocking window) to subtract background scatter from the sample holder and substrate:
   $$I_{net}(\theta) = I_{corr}(\theta) - \sum_{k=0}^3 c_k \theta^k$$
4. **Multi-Peak Gaussian Fitting**: Overlapping domain tilts were resolved by fitting a multi-component Gaussian model using non-linear least squares:
   $$I_{net}(\theta) = \sum_{j} h_j \exp\left(-\frac{(\theta - \theta_{0,j})^2}{2 \sigma_j^2}\right)$$
   where $h_j$ is the net peak height, $\theta_{0,j}$ is the peak center, and $FWHM_j = 2.355 \sigma_j$.
5. **2D Polar Texture Reconstruction (Pole Figures)**: 1D rocking curves collected at azimuthal rotations $\phi \in [0^\circ, 180^\circ]$ were projected onto a 2D polar coordinate system where:
   - Polar radius $r = \chi = \theta - \theta_{0}$ represents the crystallite tilt angle relative to the film normal.
   - Polar angle represents the azimuthal rotation $\phi$.
   The grid was expanded to $360^\circ$ by applying standard crystallographic 2-fold inversion symmetry: $I(\chi, \phi + 180^\circ) = I(-\chi, \phi)$. Grid interpolation was performed using a cubic spline function.

## 3. Background Subtraction and Peak Deconvolution Visualisation
To demonstrate the physical validity of the background subtraction and multi-peak Gaussian fitting procedures, a representative azimuthal scan from sample `SH-124-B3` ($\phi = 60^\circ$) is shown in Figures A1 and A2.

Figure A1 shows the raw rocking curve intensity alongside the fitted 3rd-order polynomial background baseline and the total fit envelope in both linear and logarithmic scales. The background baseline captures the diffuse substrate scatter, and the logarithmic plot demonstrates that the fit envelope matches the experimental data points across more than two orders of magnitude in intensity, justifying the baseline subtraction rationale.

![Figure A1: Background subtraction profile for sample `SH-124-B3` at $\phi = 60^\circ$. (a) Linear scale and (b) logarithmic scale showing raw intensity, fitted background baseline, total fit envelope, and individual peak profiles.](fig_a1_background_subtraction.png)

Figure A2 shows the resulting baseline-corrected net intensity profile. Subtracting the baseline isolates the crystalline Bragg reflections at the $y=0$ baseline. Non-linear least-squares fitting of the net profile resolves the individual Gaussian tilt components, showing excellent agreement with the experimental net data points and confirming that individual domain tilts can be deconvoluted from the overlapping Bragg reflections.

![Figure A2: Deconvoluted net rocking curve and individual Gaussian peak components for sample `SH-124-B3` at $\phi = 60^\circ$ inside the domain tilt range ($\theta \in [10.0^\circ, 14.5^\circ]$).](fig_a2_peak_deconvolution.png)
"""

# Replace placeholders
report_md = report_md.replace("{table_2d_content}", table_2d_content)
report_md = report_md.replace("{table_rc_content}", table_rc_content)

# Write to file
md_filepath = os.path.join(summary_report_dir, "synthesis_xrd_comprehensive_report.md")
with open(md_filepath, "w") as f:
    f.write(report_md)
print(f"Written Markdown report to {md_filepath}")

# Compile to PDF using pandoc
pdf_filepath = os.path.join(summary_report_dir, "synthesis_xrd_comprehensive_report.pdf")
cmd = [
    "pandoc",
    "synthesis_xrd_comprehensive_report.md",
    "-o", "synthesis_xrd_comprehensive_report.pdf",
    "--pdf-engine=pdflatex"
]
print("Compiling report to PDF...")
result = subprocess.run(cmd, cwd=summary_report_dir, capture_output=True, text=True)
if result.returncode == 0:
    print(f"Successfully compiled PDF report at {pdf_filepath}")
else:
    print("Error compiling PDF:")
    print("STDOUT:", result.stdout)
    print("STDERR:", result.stderr)
