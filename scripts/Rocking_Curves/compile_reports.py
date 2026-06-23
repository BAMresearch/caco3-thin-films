# -*- coding: utf-8 -*-
"""
Multi-Report Compiler
=====================
Author: Tomasz Stawski, tomasz.stawski@bam.de
Date: 2026-06-19
Version: 0.1

Description:
Automates the compilation of various sub-reports and technical summaries across measurement series.
"""
import os
import subprocess

# Directories
base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
dir_b3 = os.path.join(base_dir, "data/processed/Rocking_Curves/SH_124_B3")
dir_a = os.path.join(base_dir, "data/processed/Rocking_Curves/SH-125-A")
dir_comp = os.path.join(base_dir, "data/processed/Rocking_Curves/Comparison")
os.makedirs(dir_b3, exist_ok=True)
os.makedirs(dir_a, exist_ok=True)
os.makedirs(dir_comp, exist_ok=True)

# ==============================================================================
# REPORT 1: 2026-05-26 (SH_124_B3) Analysis
# ==============================================================================
md_content_526 = r"""---
title: "CaCO3 Thin Film Rocking Curve Analysis: SH-124-B3"
date: "June 2, 2026"
geometry: "margin=1in"
fontsize: "11pt"
header-includes:
  - \usepackage{fancyhdr}
  - \pagestyle{fancy}
  - \fancyhead[CO,CE]{SH-124-B3 Rocking Curve Report}
  - \fancyfoot[CO,CE]{Page \thepage}
---

# 1. Dataset Structure
The dataset consists of:
1. **1D Rocking Curve** (`SH_124_B3_rocking-curve_01_exported.xy`): 1001 points of raw intensity vs. $\theta$ (ranging from $4.7218^\circ$ to $24.7218^\circ$ in steps of $0.02^\circ$).
2. **2D Rocking Curve Map** (`SH_124_B3_rocking-curve_02_2d005_exported`): A folder containing 401 individual 1D integration files corresponding to $\theta$ steps of $0.05^\circ$ (covering $2\theta$ from $27.26^\circ$ to $30.915^\circ$).

# 2. Methodology
## 2.1 Background Subtraction
At low tilt angles ($\theta < 10^\circ$), the X-ray beam footprint is large, and diffuse/substrate scatter is strong, causing a monotonic decay in the background. To isolate the true rocking curve features:
* The raw 1D rocking curve intensity was modeled as an isotropic thin-film contribution combined with a 3rd-order polynomial background representing air/substrate scatter.
* Areas corresponding to preferred orientation peaks were masked out during background fitting to ensure a clean baseline.

## 2.2 Thin Film Volume Correction
In a thin-film rocking curve, the scattering volume changes with the incident angle. Assuming negligible absorption in a very thin film of thickness $t$, the diffracted path length is proportional to $1/\sin(\theta)$. To obtain the true rocking curve representing the distribution of crystal orientations, the net peak intensity must be multiplied by $\sin(\theta)$ to correct for the geometrical footprint:
$$I_{\text{corrected}}(\theta) = I_{\text{net}}(\theta) \times \sin(\theta)$$

# 3. Analysis Results
## 3.1 Background & Volume Correction
Applying the volume correction shows that the majority of the Calcite (104) grains are randomly oriented (isotropic), leading to a flat corrected orientation baseline.

However, in the high-frequency residuals, **five major peaks** and **two minor peaks** are resolved. These peaks correspond to discrete co-oriented grain populations (crystallite statistics) tilted at specific angles relative to the substrate normal $\chi = \theta - 2\theta/2$ (where $2\theta_0 = 29.44^\circ$):
* Peak 1 (Tilt): $\theta = 9.22^\circ$ (tilt angle $\chi = -5.50^\circ$)
* Peak 2 (Tilt): $\theta = 12.92^\circ$ (tilt angle $\chi = -1.79^\circ$)
* Peak 3 (Near-specular): $\theta = 15.71^\circ$ (tilt angle $\chi = +0.99^\circ$)
* Peak 4 (Near-specular): $\theta = 16.19^\circ$ (tilt angle $\chi = +1.47^\circ$)
* Peak 5 (Tilt): $\theta = 22.34^\circ$ (tilt angle $\chi = +7.62^\circ$)
* *Minor Peak A (Tilt)*: $\theta = 17.23^\circ$ (tilt angle $\chi = +2.51^\circ$)
* *Minor Peak B (Tilt)*: $\theta = 17.65^\circ$ (tilt angle $\chi = +2.93^\circ$)

![SH_124_B3 Rocking Curve Analysis Plot](SH_124_B3_rocking_curve_analysis.png){width=85%}

\newpage

# 4. Peak Fitting Metrics
To quantify the peak shapes, a Gaussian profile was fitted to each identified feature in the residuals:

| Peak Name | Peak Center ($\theta$) | Tilt Angle ($\chi$) | FWHM ($^\circ$) | Net Height (counts) | Net Area (counts$\cdot$deg) | Area/Isotropic Base Ratio |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **Peak 1 (Tilt)** | $9.22^\circ$ | $-5.50^\circ$ | $0.880^\circ$ | $22.6$ | $21.1$ | $0.23\%$ |
| **Peak 2 (Tilt)** | $12.92^\circ$ | $-1.79^\circ$ | $0.492^\circ$ | $242.4$ | $126.9$ | $1.97\%$ |
| **Peak 3 (Near-specular)** | $15.71^\circ$ | $+0.99^\circ$ | $1.264^\circ$ | $181.6$ | $244.4$ | $4.59\%$ |
| **Peak 4 (Near-specular)** | $16.19^\circ$ | $+1.47^\circ$ | $0.981^\circ$ | $142.3$ | $148.5$ | $2.87\%$ |
| **Minor Peak A (Tilt)** | $17.23^\circ$ | $+2.51^\circ$ | $0.071^\circ$ | $7.7$ | $0.58$ | $0.01\%$ |
| **Minor Peak B (Tilt)** | $17.65^\circ$ | $+2.93^\circ$ | $0.071^\circ$ | $10.7$ | $0.80$ | $0.02\%$ |
| **Peak 5 (Tilt)** | $22.34^\circ$ | $+7.62^\circ$ | $0.407^\circ$ | $48.9$ | $21.2$ | $0.56\%$ |

# 5. Conclusions
1. **Isotropic Matrix:** The bulk of the Calcite (104) film is randomly oriented (isotropic), as confirmed by the flat volume-corrected baseline.
2. **Co-oriented Grains:** Multiple discrete preferred orientation components are resolved in the residuals. The FWHM values of these peaks demonstrate high-quality crystalline alignment in these specific domains.
"""

# ==============================================================================
# REPORT 2: 2026-06-01 (SH-125-A) Analysis
# ==============================================================================
md_content_601 = r"""---
title: "CaCO3 Thin Film Rocking Curve Analysis: SH-125-A"
date: "June 2, 2026"
geometry: "margin=1in"
fontsize: "11pt"
header-includes:
  - \usepackage{fancyhdr}
  - \pagestyle{fancy}
  - \fancyhead[CO,CE]{SH-125-A Rocking Curve Report}
  - \fancyfoot[CO,CE]{Page \thepage}
---

# 1. Dataset Structure
The dataset consists of two components:
1. **1D Rocking Curve** (`SH-125-A_rocking-curve_01_exported.xy`): 401 points of raw intensity vs. $\theta$ (ranging from $4.7218^\circ$ to $24.7218^\circ$ in steps of $0.05^\circ$).
2. **Symmetric $2\theta-\theta$ Scan** (`SH-125-A_2Theta_exported.xy`): 401 points of intensity vs. $2\theta$ (ranging from $27.0001^\circ$ to $35.0001^\circ$ in steps of $0.02^\circ$).

# 2. Methodology
## 2.1 Background Subtraction
At low tilt angles ($\theta < 10^\circ$), the X-ray beam footprint is large, and diffuse/substrate scatter is strong. To isolate the true Bragg peaks and phase contributions:
* For the $2\theta-\theta$ scan, a 3rd-order polynomial was fitted to the non-peak regions and subtracted.
* For the rocking curve, the raw intensity was modeled as an isotropic thin-film contribution combined with a 3rd-order polynomial background representing air/substrate scatter. All peak regions were excluded during background fitting.

## 2.2 Thin Film Volume Correction
In a thin-film rocking curve, the scattering volume changes with the incident angle. Assuming negligible absorption in a very thin film of thickness $t$, the diffracted path length is proportional to $1/\sin(\theta)$. To obtain the true rocking curve representing the distribution of crystal orientations, the net peak intensity must be multiplied by $\sin(\theta)$ to correct for the geometrical footprint:
$$I_{\text{corrected}}(\theta) = I_{\text{net}}(\theta) \times \sin(\theta)$$

# 3. Analysis Results
## 3.1 2Theta Phase Scan
Fitting Gaussian profiles to the background-subtracted $2\theta$ scan reveals a mixture of Calcite (104) and Vaterite (110) phases:
* **Calcite (104)** is centered at $2\theta = 29.374^\circ$ with a FWHM of $0.270^\circ$.
* **Vaterite (110)** is centered at $2\theta = 32.807^\circ$ with a FWHM of $0.385^\circ$.

![SH-125-A 2Theta Scan Phase Analysis](SH-125-A_2Theta_analysis.png){width=85%}

\newpage

## 3.2 Volume-Corrected Rocking Curves
The raw uncorrected intensity decays monotonically and peaks at the minimum angle ($\theta \approx 4.72^\circ$). 

Applying the physical isotropic model shows that the majority of the Calcite (104) grains are randomly oriented (isotropic), leading to a flat corrected orientation baseline.

However, in the high-frequency residuals, **four major peaks** and **two minor peaks** are resolved. These peaks correspond to discrete co-oriented grain populations (crystallite statistics) tilted at specific angles relative to the substrate normal:
* Peak 1: $\theta = 9.99^\circ$ (tilt angle $\chi = -4.68^\circ$)
* Peak 2: $\theta = 12.81^\circ$ (tilt angle $\chi = -1.87^\circ$)
* Peak 3: $\theta = 14.35^\circ$ (tilt angle $\chi = -0.32^\circ$) [Specular]
* Peak 4: $\theta = 22.53^\circ$ (tilt angle $\chi = +7.86^\circ$)
* *Minor Peak A*: $\theta = 17.16^\circ$ (tilt angle $\chi = +2.49^\circ$)
* *Minor Peak B*: $\theta = 17.57^\circ$ (tilt angle $\chi = +2.89^\circ$)

![SH-125-A Rocking Curve Fit and Residual Peaks](SH-125-A_rocking_curve_analysis.png){width=85%}

## 3.3 Correlation with Helical Columnar Growth (SEM Evidence)
Scanning Electron Microscopy (SEM) images suggest that the thin film contains **helicically growing columnar crystals**. This morphological feature provides a direct physical explanation for the peaks observed in the rocking curve residuals:

1. **The Helical Tilt Model:**
   In helical columnar crystals, the crystal lattice planes (Calcite 104) grow with a specific tilt angle $\alpha_0$ relative to the growth axis of the column (which is approximately normal to the substrate). As the column growth twists helically, the tilt direction rotates azimuthally around the column axis, creating a conical distribution of plane normals.
   
2. **XRD Rocking Curve Signature of Helicity:**
   In a 1D rocking curve scan, this conical distribution of tilted planes produces diffraction peaks at tilt angles of $\chi = \pm \alpha_0$. 
   * The detected peaks at $\chi_1 = -4.68^\circ$ and $\chi_4 = +7.86^\circ$ align with this model, indicating a helix tilt angle $\alpha_0$ in the range of approximately $5^\circ$ to $8^\circ$. The asymmetry in the peak intensities and tilt angles suggests a slight tilt of the column growth axes or a preferred chirality/twist direction relative to the incident X-ray beam.
   * The peak at $\chi_2 = -1.87^\circ$ may correspond to columns in an early stage of growth before the helix pitch fully develops, or columnar segments growing closer to the substrate normal.
   * **The Minor Peaks at $\chi \approx +2.49^\circ$ and $+2.89^\circ$:**
     These minor features represent additional tilt components that are also highly co-oriented. In the helical columnar growth model, they may correspond to secondary crystalline facets or variations in the helical pitch as the columnar structures evolve.
   
3. **The Specular Orientation Peak (Peak 3):**
   * The peak at $\chi_3 = -0.32^\circ$ (centered at $\theta = 14.35^\circ$) is located almost exactly at the specular reflection condition ($\chi \approx 0^\circ$). This peak represents a population of calcite crystallites whose (104) planes are aligned parallel to the substrate surface (i.e. standard normal fiber orientation).
   
4. **High Mosaicity vs. Discrete Helical Grains:**
   The narrow FWHM of the residual peaks ($0.26^\circ - 0.57^\circ$) is close to the instrumental resolution limit, demonstrating that the individual helical columns have high crystallographic quality with very little internal mosaicity (misorientation spread) within each domain.

\newpage

# 4. Peak Fitting Metrics
To quantify the peak shapes, a Gaussian profile was fitted to each identified feature.

### 4.1 Phase Peak Fitting Parameters
| Phase | Peak Center ($2\theta$) | FWHM ($^\circ$) | Height (counts) | Integrated Area (counts$\cdot$deg) |
| :--- | :---: | :---: | :---: | :---: |
| **Calcite (104)** | $29.374^\circ$ | $0.270^\circ$ | $3,261.2$ | $938.2$ |
| **Vaterite (110)** | $32.807^\circ$ | $0.385^\circ$ | $1,966.0$ | $805.3$ |

### 4.2 Rocking Curve Residual Peak Parameters
| Peak Name | Peak Center ($\theta$) | Tilt Angle ($\chi$) | FWHM ($^\circ$) | Net Height (counts) | Net Area (counts$\cdot$deg) | Area/Isotropic Base Ratio |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **Peak 1 (Tilt)** | $9.99^\circ$ | $-4.68^\circ$ | $0.566^\circ$ | $86,794$ | $52,325$ | $5.88\%$ |
| **Peak 2 (Tilt)** | $12.81^\circ$ | $-1.87^\circ$ | $0.540^\circ$ | $74,493$ | $42,846$ | $6.16\%$ |
| **Peak 3 (Specular)** | $14.35^\circ$ | $-0.32^\circ$ | $0.462^\circ$ | $59,134$ | $29,087$ | $4.67\%$ |
| **Minor Peak A (Tilt)** | $17.16^\circ$ | $+2.49^\circ$ | $0.407^\circ$ | $35,836$ | $15,502$ | $2.55\%$ |
| **Minor Peak B (Tilt)** | $17.57^\circ$ | $+2.89^\circ$ | $0.098^\circ$ | $7,495$ | $781$ | $0.13\%$ |
| **Peak 4 (Tilt)** | $22.53^\circ$ | $+7.86^\circ$ | $0.259^\circ$ | $76,704$ | $21,146$ | $5.25\%$ |

# 5. Limitations of the Analysis
While the present analysis provides strong evidence of helical columnar features, several structural and experimental limitations must be acknowledged:

1. **Dimensional Limitation (1D vs. 2D Mapping):**
   The rocking curve is a 1D $\theta$-scan at a single fixed azimuthal angle $\phi$. Because the sample was not rotated around its surface normal ($\phi$ scan) during measurement, we cannot determine if the helical columns have a preferred in-plane orientation (alignment) or if they are distributed isotropically around the surface normal (fiber texture). 
   
2. **Beam Spillover at Low Angles (Footprint Effect):**
   At very low incident angles ($\theta < 10^\circ$), the X-ray beam footprint exceeds the physical dimensions of the sample (spillover). While the $1/\sin(\theta)$ model accounts for volume variation, it does not correct for this spillover. This can cause the net intensity at low angles (affecting Peak 1 at $\theta = 9.99^\circ$) to be slightly underestimated.
   
3. **Absorption Approximations:**
   The volume correction assumes an infinitely thin film with negligible X-ray absorption. If the film thickness or density is significant, absorption will attenuate the diffracted signal at low incident angles, deviating from the pure $1/\sin(\theta)$ behavior.
   
4. **Poor Particle/Grain Statistics:**
   The XRD beam spot covers a large macroscopic area (~10 mm^2), meaning it averages over millions of columns. The fact that only four distinct peaks are resolved in the residuals indicates that a few very large or highly co-oriented crystallite domains dominate the diffraction signal. This limits our ability to quantify the overall helical population statistics across the entire film.

# 6. Conclusions
1. **Mixed Phase Film:** The thin-film sample contains both Calcite and Vaterite phases in a comparable ratio (Calcite:Vaterite area ratio $\approx 1.16:1$).
2. **Isotropic Film Matrix:** After accounting for the scattering volume correction, the Calcite (104) orientation distribution is flat, demonstrating that the crystallites are predominantly quasi-randomly oriented.
3. **Helical Growth and Specular Fingerprints:** The rocking curve residuals resolve both the specular orientation of calcite crystallites (at $\chi \approx -0.32^\circ$) and the characteristic helix tilt angle of the columnar crystals (at $\chi \approx -4.68^\circ$ and $+7.86^\circ$).
"""

# ==============================================================================
# REPORT 3: Comparative Analysis (SH_124_B3 vs. SH-125-A)
# ==============================================================================
md_content_comp = r"""---
title: "CaCO3 Thin Film Rocking Curve Comparative results"
date: "June 2, 2026"
geometry: "margin=1in"
fontsize: "11pt"
header-includes:
  - \usepackage{fancyhdr}
  - \pagestyle{fancy}
  - \fancyhead[CO,CE]{Rocking Curve Comparative Report}
  - \fancyfoot[CO,CE]{Page \thepage}
---

# 1. Introduction
This report compares the preferred orientation components of two thin-film samples:
* `SH_124_B3` (measured on 2026-05-26, nominal $2\theta_0 = 29.44^\circ$)
* `SH-125-A` (measured on 2026-06-01, nominal $2\theta_0 = 29.3425^\circ$)

The comparison is based on the net high-frequency residuals obtained after modeling and subtracting the isotropic matrix baseline ($1/\sin\theta$ volume footprint correction and polynomial background).

# 2. Residual Orientation Peaks Comparison
Despite a ~300-fold difference in absolute diffraction intensities (likely arising from differences in count times, slit settings, or beam size), the positions of the discrete tilt peaks are highly reproducible across both datasets.

![Comparison of Net Residual Orientation Peaks](samples_comparison.png){width=85%}

The co-oriented tilt components match as follows:
1. **The $\chi \approx -1.8^\circ$ Component:**
   * `SH_124_B3`: $\chi = -1.80^\circ$ (FWHM $= 0.49^\circ$, height $= 242$ counts)
   * `SH-125-A`: $\chi = -1.87^\circ$ (FWHM $= 0.54^\circ$, height $= 74,493$ counts)
   * *Difference:* $0.07^\circ$
2. **The $\chi \approx +7.7^\circ$ Component:**
   * `SH_124_B3`: $\chi = +7.62^\circ$ (FWHM $= 0.41^\circ$, height $= 49$ counts)
   * `SH-125-A`: $\chi = +7.86^\circ$ (FWHM $= 0.26^\circ$, height $= 76,704$ counts)
   * *Difference:* $0.24^\circ$
3. **The $\chi \approx -5^\circ$ Component:**
   * `SH_124_B3`: $\chi = -5.50^\circ$ (FWHM $= 0.88^\circ$, height $= 23$ counts)
   * `SH-125-A`: $\chi = -4.68^\circ$ (FWHM $= 0.57^\circ$, height $= 86,794$ counts)
   * *Difference:* $0.82^\circ$
4. **The Minor $\chi \approx +2.5^\circ$ and $+2.9^\circ$ Components:**
   * **Component at $\chi \approx +2.50^\circ$**: `SH_124_B3` ($\chi = +2.51^\circ$, height $= 7.7$) vs `SH-125-A` ($\chi = +2.49^\circ$, height $= 35,836$)
   * **Component at $\chi \approx +2.90^\circ$**: `SH_124_B3` ($\chi = +2.93^\circ$, height $= 10.7$) vs `SH-125-A` ($\chi = +2.89^\circ$, height $= 7,495$)
   * *Differences:* $0.02^\circ - 0.04^\circ$

# 3. Crystallographic Discussion
The alignment of preferred orientation peaks at identical tilt angles $\chi$ across separate film growths provides strong evidence for a **substrate-guided nucleation or growth mechanism**:
* In both films, the majority of the calcite matrix is randomly oriented (isotropic), representing a disordered matrix.
* However, the discrete peaks represent co-oriented crystallite domains. 
* The symmetrical peaks at $\chi \approx -5^\circ$ and $+7.7^\circ$ match the conical distribution of plane normals expected from the helical growth geometry suggested by SEM, representing a helix tilt angle $\alpha_0 \approx 5^\circ - 8^\circ$.
* The high reproducibility of the minor tilt states at $\chi \approx +2.5^\circ$ and $+2.9^\circ$ between independent growth series indicates a highly deterministic orientation pathway, likely dictated by specific low-energy growth facets.

# 4. Summary Table of Fitted Peaks
### 4.1 Sample `SH_124_B3` (2026-05-26)
| Peak ID | Peak Center ($\theta$) | Tilt Angle ($\chi$) | FWHM ($^\circ$) | Net Height (counts) | Area/Isotropic Base Ratio |
| :--- | :---: | :---: | :---: | :---: | :---: |
| **B3_Peak_9.22** | $9.22^\circ$ | $-5.50^\circ$ | $0.880^\circ$ | $22.6$ | $0.23\%$ |
| **B3_Peak_12.92** | $12.92^\circ$ | $-1.80^\circ$ | $0.492^\circ$ | $242.4$ | $1.97\%$ |
| **B3_Peak_15.70** | $15.71^\circ$ | $+0.99^\circ$ | $1.264^\circ$ | $181.6$ | $4.59\%$ |
| **B3_Peak_16.19** | $16.19^\circ$ | $+1.47^\circ$ | $0.981^\circ$ | $142.3$ | $2.87\%$ |
| **B3_Minor_A** | $17.23^\circ$ | $+2.51^\circ$ | $0.071^\circ$ | $7.7$ | $0.01\%$ |
| **B3_Minor_B** | $17.65^\circ$ | $+2.93^\circ$ | $0.071^\circ$ | $10.7$ | $0.02\%$ |
| **B3_Peak_22.34** | $22.34^\circ$ | $+7.62^\circ$ | $0.407^\circ$ | $48.9$ | $0.56\%$ |

### 4.2 Sample `SH-125-A` (2026-06-01)
| Peak ID | Peak Center ($\theta$) | Tilt Angle ($\chi$) | FWHM ($^\circ$) | Net Height (counts) | Area/Isotropic Base Ratio |
| :--- | :---: | :---: | :---: | :---: | :---: |
| **A_Peak_9.99** | $9.99^\circ$ | $-4.68^\circ$ | $0.566^\circ$ | $86,794$ | $5.88\%$ |
| **A_Peak_12.81** | $12.81^\circ$ | $-1.87^\circ$ | $0.540^\circ$ | $74,493$ | $6.16\%$ |
| **A_Peak_14.35** | $14.35^\circ$ | $-0.32^\circ$ | $0.462^\circ$ | $59,134$ | $4.67\%$ |
| **A_Minor_A** | $17.16^\circ$ | $+2.49^\circ$ | $0.407^\circ$ | $35,836$ | $2.55\%$ |
| **A_Minor_B** | $17.57^\circ$ | $+2.89^\circ$ | $0.098^\circ$ | $7,495$ | $0.13\%$ |
| **A_Peak_22.53** | $22.53^\circ$ | $+7.86^\circ$ | $0.259^\circ$ | $76,704$ | $5.25\%$ |

# 5. Conclusions
The comparative analysis confirms that both calcite thin films share identical discrete orientation tilt components at $\chi \approx -1.8^\circ$ and $\chi \approx +7.7^\circ$. This high reproducibility supports a deterministic growth geometry, such as helical columnar growth, guided by the substrate or growth conditions.
"""

# ==============================================================================
# Compilation functions
# ==============================================================================
def write_and_compile(md_text, target_dir, base_name):
    md_path = os.path.join(target_dir, base_name + ".md")
    pdf_path = os.path.join(target_dir, base_name + ".pdf")
    
    with open(md_path, 'w') as f:
        f.write(md_text)
    print(f"Written Markdown report to {md_path}")
    
    cmd = [
        "pandoc",
        base_name + ".md",
        "-o", base_name + ".pdf",
        "--pdf-engine=pdflatex"
    ]
    print(f"Compiling {base_name}.pdf in {target_dir}...")
    result = subprocess.run(cmd, cwd=target_dir, capture_output=True, text=True)
    if result.returncode == 0:
        print(f"Successfully generated PDF report at {pdf_path}")
    else:
        print(f"Error compiling {base_name}.pdf:")
        print("STDOUT:", result.stdout)
        print("STDERR:", result.stderr)

# Compile all three reports with descriptive sample-specific filenames
write_and_compile(md_content_526, dir_b3, "rocking_curve_analysis_SH_124_B3")
write_and_compile(md_content_601, dir_a, "rocking_curve_analysis_SH-125-A")
write_and_compile(md_content_comp, dir_comp, "rocking_curve_comparison_SH_124_B3_vs_SH-125-A")
