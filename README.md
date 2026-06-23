# CaCO3 Thin Film Crystallographic Characterisation and Analysis

## Authors & Credits

- **Author**: Tomasz Stawski (tomasz.stawski@bam.de)
- **Version**: 1.0
- **Date**: 2026-06-23

---

This repository contains the data processing, peak deconvolution, and plotting scripts for analysing 2D-XRD detector frames and azimuthal rocking curves of calcium carbonate ($CaCO_3$) thin films.

The project is structured to transition raw experimental data (from stationary detector frames and rocking curve scans) into baseline-corrected, multi-peak fitted datasets, and ultimately compile them into high-resolution publication figures and technical reports.

---

## Directory Structure

The repository is organised as follows:

```
thin films/
├── .gitignore                      # Standard ignore rules (pycache, envs, IDE configs)
├── LICENSE                         # MIT License
├── README.md                       # Project documentation and guide
├── requirements.txt                # Python dependencies
├── run_entire_analysis.py          # Master analysis coordinator entry point
├── caco3_diffraction_pipeline.py   # Core processing and plotting pipeline script
├── data/
│   ├── raw/                        # Original experimental files (unmodified)
│   │   ├── 2D_XRD/                 # Bruker 2D flat-panel detector frames (.gfrm, .h5, .brml)
│   │   ├── crystalData/            # Vaterite and calcite crystallographic database reference files
│   │   └── Rocking_Curves/         # Azimuthal theta rocking scans (.brml and .raw binary blocks)
│   │       ├── SH-104-1/
│   │       ├── SH-124-B3/
│   │       ├── SH-125-A/
│   │       ├── SH-125-G/
│   │       └── Reference/          # Reference single crystal scans
│   └── processed/                  # Output tables, baseline-corrected profiles, metrics
│       ├── 2D_XRD/                 # Integrated 1D profile Excel sheets
│       ├── Symmetric_Scans/        # Exported symmetric 2theta scan xy tables
│       └── Rocking_Curves/         # Deconvoluted rocking curve peak parameters and CSV files
│           ├── SH-104-1/
│           ├── SH-124-B3/
│           ├── SH-125-A/
│           ├── SH-125-G/
│           └── Reference/
├── results/
│   ├── figures/                    # Publication-quality plots (PNG and vector SVG formats)
│   └── reports/                    # Final PDF and Markdown comprehensive reports
└── scripts/                        # Archives of individual processing and plotting scripts
    ├── 2D_XRD/
    └── Rocking_Curves/
```

---

## Data Structure

The dataset spans multiple file formats representing different stages of data acquisition and reduction:

### 1. Raw Data (`data/raw/`)
* **`.gfrm`**: Bruker 2D flat-panel detector frames containing spatial intensity data from stationary measurements.
* **`.h5`**: Hierarchical Data Format (HDF5) files containing resampled arrays, converted from `.gfrm` files using the `silx` library for open access.
* **`.brml`**: Compressed Bruker XML-wrapped files containing 1D scanning data, including symmetric $2\theta-\theta$ scans and azimuthal rocking curves ($\theta$ scans).
* **`.raw`**: Legacy Bruker binary data blocks containing raw measurement arrays, used as fallbacks if `.brml` files are corrupted or missing.
* **`.txt`**: Tab-delimited crystal reference data (e.g. `Calcite__0000985.txt`, `Vaterite__0004854.txt`) listing reflection indices ($hkl$), $d$-spacings, and relative intensities.

### 2. Processed Data (`data/processed/`)
* **`.xlsx`**: Multi-sheet Excel workbooks containing scan metadata, interpolated grid points, and corrected cake profiles.
* **`.xy`**: Clean, space-delimited text profiles (2theta and intensity) suitable for external plotting or refinement software.
* **`.csv`**: Tabular data containing fitted peak parameters, net residual intensities, and model baselines for each azimuthal angle ($\phi$).

---

## Specimen Details

The analysis covers four distinct specimens representing different growth conditions and mineralogical compositions of calcium carbonate ($CaCO_3$) thin films on substrates:

* **`SH-104-1` (Reference Substrate)**: An uncoated reference substrate. The stationary 2D-XRD profile exhibits a mainly isotropic character, providing a baseline measurement of the substrate lattice contributions and instrumental background.
* **`SH-124-B3` (Pure Calcite Film)**: A $CaCO_3$ thin film grown under Condition B3, consisting of pure calcite phase. Rocking curve deconvolution reveals distinct, narrow domain tilts ($\chi \approx -1.8^\circ$ and $+7.7^\circ$) that suggest a helical columnar microstructure growth mechanism.
* **`SH-125-A` (Mixed Calcite-Vaterite Film)**: A thin film grown under Condition A, exhibiting a biphasic mineralogical composition with coexisting calcite and vaterite phases. This specimen is used to validate background-subtraction stability and sample-absorption volume corrections.
* **`SH-125-G` (Mixed Calcite-Vaterite Film)**: A thin film grown under Condition G, featuring a mixed calcite-vaterite composition. The transient vaterite (110) reflection exhibits a highly localized azimuthal selectiveness, satisfying the Bragg condition only at specific rotations ($\phi = 30^\circ$ and $60^\circ$), which indicates template-guided in-plane epitaxial locking.

---

## Functional Capabilities

The core analysis pipeline (`caco3_diffraction_pipeline.py`) implements the following processing operations:

### 1. 2D-XRD Detector Frame Integration
* **Grid Interpolation**: Converts flat-panel Bruker detector coordinates into polar cake plots ($2\theta$ vs. $\phi$) using bivariate linear grid interpolation.
* **Azimuthal Integration**: Integrates the 2D cake plot over the azimuthal angle range to produce 1D profiles.
* **Morphological Baseline Correction**: Estimates the background signal using rolling minimum-maximum filters and subtracts it to isolate Bragg diffraction peaks.

### 2. Symmetric Scan Extraction and Fitting
* **Profile Extraction**: Automatically parses symmetric $2\theta-\theta$ scans for each specimen across multiple azimuthal angles ($\phi$).
* **Polynomial Baseline Correction**: Models diffuse background scatter with a third-order polynomial.
* **Peak Parameterisation**: Fits Gaussian functions to the calcite (104) and vaterite (110) reflections to extract peak centers, heights, and integrated areas.

### 3. Rocking Curve Background Correction
* **Isotropic Volume Correction**: Accounts for sample-geometry absorption effects using an isotropic $I_0 / \sin(\theta)$ model.
* **Iterative Baseline Fitting**: Combines the volume correction with a polynomial background, iteratively fitting the model to regions free of Bragg reflections.

### 4. Bragg Peak Deconvolution
* **Multi-Peak Log-Scale Fitting**: Fits overlapping rocking curve features in log-intensity space to resolve distinct tilt components.
* **Parameter Extraction**: Computes tilt angles ($\chi$), peak heights, full width at half maximum (FWHM) values, and relative calcite/vaterite phase distributions.

### 5. Automated Publication Plot Generation
* Renders 20 publication-ready PNG and vector SVG diagrams, including waterfall plots, stacked net rocking curves, phase composition maps, and polar projection pole figures (texture maps).

---

## Installation & Requirements

The scripts require Python 3.8+ and standard scientific packages. Install the requirements using:

```bash
pip install -r requirements.txt
```

### Optional Dependency for GFRM Conversion
To convert binary `.gfrm` Bruker frame files directly, `silx` is used. If `silx` is not installed, the pipeline automatically falls back to reading the pre-converted `.gfrm.h5` files located in the raw data directory.

---

## Usage

To re-run the entire data processing pipeline and regenerate all publication-ready plots:

```bash
python run_entire_analysis.py
```

This runs the core routines from `caco3_diffraction_pipeline.py`, which populates `data/processed/` with baseline-corrected tables and exports publication-ready figures to `results/figures/`.

---

## Authors & Credits

- **Author**: Tomasz Stawski (tomasz.stawski@bam.de)
- **Version**: 1.0
- **Date**: 2026-06-23
