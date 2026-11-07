# Pevonedistat Kinetic Modeling Pipeline

A Python-based analytical pipeline designed to process UV-Vis absorbance data from drug release experiments, apply cumulative sampling volume and dilution corrections, perform kinetic model fitting, and generate publication-ready plots.

This script evaluates drug release kinetics across four fundamental models: **Zero-Order**, **First-Order**, **Higuchi**, and **Korsmeyer–Peppas**.

---

## Disclaimer & Data Availability

**Note:** Raw experimental release data files (`release_data.csv`) are **not included** in this public repository to protect unpublished research findings and institutional intellectual property. 

To run the pipeline on your own dataset, provide a CSV file named `release_data.csv` in the root directory formatted according to the schema in the [Input Data Format](#-input-data-format) section below.

---

## Features

* **Dilution & Volume Correction:** Automatically accounts for dilution factors ($\times 10$) and cumulative mass removal during periodic sample replacement ($V_\text{sample} = 0.5\text{ mL}$, $V_\text{total} = 5.0\text{ mL}$).
* **Baseline Noise Handling:** Applies lower-bound clipping (`.clip(lower=0)`) to eliminate non-physical negative concentrations resulting from spectrophotometric baseline noise, preventing $\ln(x)$ and $\sqrt{x}$ execution errors.
* **Methodological Validity Constraints:** Enforces kinetic domain restrictions, filtering data to $0 < M_t/M_\infty \le 0.60$ for Higuchi and Korsmeyer–Peppas models to remain aligned with theoretical diffusion boundaries.
* **Automated Parameter Extraction:** Computes correlation coefficients ($R^2$), rate constants ($k_0, k_1, k_H, k_{KP}$), and the release exponent ($n$).
* **High-Resolution Figures:** Exports $300\text{ DPI}$ graphical outputs for primary cumulative release curves and 4-panel kinetic regression grids.

---

## Mathematical Framework

The pipeline computes fractional release ($M_t/M_\infty$) using the maximum observed experimental cumulative release ($M_\infty$).

| Model | Linearized Equation | Parameters | Model Validity / Domain |
| :--- | :--- | :--- | :--- |
| **Zero-Order** | $\frac{M_t}{M_\infty} = q_0 + k_0 t$ | $k_0$: Rate constant ($\text{h}^{-1}$) | Full timeline ($t > 0$) |
| **First-Order** | $\ln\left(1 - \frac{M_t}{M_\infty}\right) = -k_1 t$ | $k_1$: Rate constant ($\text{h}^{-1}$) | $0 < \frac{M_t}{M_\infty} < 1$ |
| **Higuchi** | $\frac{M_t}{M_\infty} = k_H \sqrt{t}$ | $k_H$: Higuchi constant ($\text{h}^{-0.5}$) | $0 < \frac{M_t}{M_\infty} \le 0.60$ |
| **Korsmeyer–Peppas** | $\ln\left(\frac{M_t}{M_\infty}\right) = \ln(k_{KP}) + n \ln(t)$ | $k_{KP}$: Rate constant, $n$: Exponent | $0 < \frac{M_t}{M_\infty} \le 0.60$ |

### Interpretation of Exponent ($n$)
For spherical/particulate geometries:
* $n \le 0.43$: Fickian diffusion dominant.
* $0.43 < n < 0.85$: Consistent with anomalous (non-Fickian) transport under standard model assumptions.
* $n \ge 0.85$: Case II transport (relaxational/swelling-controlled).

---

## Input Data Format

The pipeline expects an input CSV file formatted as follows:

```csv
time_min,absorbance
0,0.025
60,0.045
120,0.078
```

| Column Name | Type | Description |
| :--- | :--- | :--- |
| `time_min` | Numeric | Sampling time point in minutes. |
| `absorbance` | Numeric | Uncorrected raw UV-Vis absorbance ($\text{A.U.}$). |

---

## Execution & Usage

### 1. Prerequisites
Ensure Python 3.8+ is installed along with the required dependencies:

```bash
pip install numpy pandas matplotlib scipy
```

### 2. Experimental Parameters
If your experimental setup differs from the default script parameters, adjust the global variables at the top of `kinetic_analysis.py`:

```python
TOTAL_VOLUME_ML = 5.00       # Release medium total volume (mL)
SAMPLE_VOLUME_ML = 0.50      # Sample replacement volume (mL)
LOADED_MASS_MG = 0.107296137 # Total initial drug mass loaded (mg)
CALIBRATION_SLOPE = 662.29   # UV-Vis standard calibration slope
CALIBRATION_INTERCEPT = 0.025# UV-Vis calibration intercept
DILUTION_FACTOR = 10.0       # Sample dilution factor
```

### 3. Running the Pipeline
Place your `release_data.csv` file in the root directory and run:

```bash
python kinetic_analysis.py
```

---

## Output Files Generated

Upon execution, the script exports:

1. **`kinetic_parameters_summary.csv`** — Exported table containing $R^2$ values, rate constants ($k$), and the Korsmeyer–Peppas exponent ($n$).
2. **`pevonedistat_release_curve.png`** — High-resolution cumulative release profile.
3. **`kinetic_models_fit.png`** — 4-panel linear regression grid displaying model fits.
