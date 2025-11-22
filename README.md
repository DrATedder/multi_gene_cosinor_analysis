# Multi-gene Cosinor Analysis & Plotting

## Overview

`multi_gene_cosinor_plot.py` performs batch 24-hour single-component cosinor analysis for multiple genes from rhythmic expression data. It fits a cosinor model to each gene, computes rhythm parameters (MESOR, amplitude, acrophase, R²), performs a zero-amplitude (F-test) to assess rhythmic significance, and generates publication-quality grayscale figures.

For each gene, the script:

- Fits a standard 24-h cosinor model  
- Computes MESOR, amplitude, acrophase (in hours), R², amplitude SE  
- Computes a **zero-amplitude test p-value (F-test)**  
- Plots raw data (open circles), mean ± SD (closed circles), and the fitted curve  
- Exports each figure to a **PDF**  
- Produces a combined summary table (`.csv`)

---

## Method Summary

Each gene’s expression is modeled using the classical single-component cosinor:

\[
Y(t) = M + A \cos\left(\frac{2\pi t}{24} - \phi\right)
\]

where:

- \( M \): MESOR (midline estimate)  
- \( A \): amplitude  
- \( \phi \): acrophase (radians), converted to hours for reporting  
- \( t \): time (ZT)  
- \( 24 \): assumed circadian period  

The model is fit using nonlinear least squares (SciPy `curve_fit`).  
Rhythmic significance is evaluated via the **zero-amplitude test**:

\[
F = \frac{A^{2} / \mathrm{SE}(A)^{2}}{2}
\]

with degrees of freedom \(df_1 = 2\), \(df_2 = N - 3\).  
The p-value is computed from the F distribution.

---

## Script Description

**Filename:** `multi_gene_cosinor_plot.py`

This script:

1. Reads all `.tsv` files in a specified directory.  
2. Extracts gene names from filenames (e.g., `Bmal1.tsv` → `Bmal1`).  
3. Fits a 24-hour cosinor model to each dataset.  
4. Computes MESOR, amplitude, SE, acrophase (radians + hours), R².  
5. Performs a zero-amplitude F-test to obtain a p-value.  
6. Generates publication-quality grayscale figures:  
   - Raw data (open circles)  
   - Mean ± SD (closed circles)  
   - Fitted curve (solid black)  
7. Saves each figure as **`<gene>_cosinor.pdf`**.  
8. Produces a summary CSV (`cosinor_summary.csv`).

### Outputs

- **Per-gene PDF figures** (`<gene>_cosinor.pdf`)
- **Cosinor summary table** (`cosinor_summary.csv`)
- Console output containing parameter estimates for each gene

**Note**. Each figure will be produced on the same Y-axis scale to ensure meaningful comparison.

---

## Input Data Format

Each gene must be in a `.tsv` file with the following columns:

| ZT | expression ratio (goi/hk) |
|----|---------------------------|

**Required column names (exact):**

- `ZT` → time points (e.g., 0, 4, 8 … 24)  
- `expression ratio (goi/hk)` → normalized expression  

Example files:

Bmal1.tsv
Per2.tsv
Cry1.tsv


**Missing values** are allowed and automatically handled.

---

## Configuration Parameters

Edit these parameters in the script’s **SETTINGS** section:

```python
folder = "/path/to/tsv/files/"  # directory containing .tsv files
show_points = True              # show raw open-circle points
period = 24                     # assumed circadian period
save_summary = True             # export cosinor_summary.csv
save_figures = True             # export PDF figures
```

---
## Dependencies

Requires Python 3.9+ and:
```python
numpy
pandas
scipy
matplotlib
```

---

## Install

```bash
pip3 install numpy pandas scipy matplotlib
```

if `Numpy` version causes issues:

```bash
pip install "numpy<2"
```

---

## Usage

1. Place all .tsv files in your target directory.

2. Edit the folder variable in the script to point to that directory.

3. Run:

```bash
python3 multi_gene_cosinor_plot.py
```
The script will automatically:

* Fit cosinor models for each gene

* Export publication-quality plots

* Produce a parameter summary table

* Print results to the terminal

---

## Outputs

| File                  | Description                                  |
| --------------------- | -------------------------------------------- |
| `<gene>_cosinor.pdf`  | Publication-quality cosinor plot             |
| `cosinor_summary.csv` | Parameter estimates + p-values for all genes |

## Summary Tables

| column            | meaning                                  |
| ----------------- | ---------------------------------------- |
| MESOR             | estimated midline level                  |
| Amplitude         | oscillation amplitude                    |
| Amplitude_SE      | standard error of the amplitude estimate |
| Acrophase_radians | phase (radians)                          |
| Acrophase_hours   | phase converted to hours                 |
| R2                | variance explained by model              |
| ZeroAmp_pvalue    | F-test p-value for non-zero amplitude    |

---

## Interpretation of Results

**MESOR**: average expression level

**Amplitude**: oscillation strength

**Acrophase (hours)**: predicted peak time

**p-value**: statistical evidence of rhythmicity

**R²**: goodness of fit

Low amplitude or large residual variance → higher p-values.
Clear rhythmic genes show large amplitude and p ≪ 0.05.

---

## Notes

* Assumes a single-component 24-h cosinor.

* Timesteps do **not** need to be equally spaced.

* ZT24 is automatically included even if not in the raw data.

* Figures use strict grayscale for publication compatibility.

* Raw data use open circles; means use closed circles.

---

## References

Nelson W., Tong Y.L., Lee J.K., Halberg F. (1979). Methods for cosinor-rhythmometry. Chronobiologia, 6(4), 305–323.

Refinetti R., Cornélissen G., Halberg F. (2007). Procedures for numerical analysis of circadian rhythms. Biological Rhythm Research, 38(4), 275–325.
