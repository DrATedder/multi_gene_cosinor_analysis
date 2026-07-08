# Multi-gene Cosinor Analysis

Python scripts for batch circadian cosinor analysis of rhythmic gene expression data using classical single-component (24-hour) cosinor models. The repository provides three complementary workflows, ranging from single-gene analysis to multi-tissue comparative analyses with multiple-testing correction and bootstrap-based inference.

| Script | Description |
|--------|-------------|
| [`multi_gene_cosinor_plot.py`](#multi_gene_cosinor_plotpy) | Batch cosinor analysis and publication-quality plots for multiple genes. |
| [`cosinor_batch_by_tissue.py`](#cosinor_batch_by_tissuepy) | Multi-tissue cosinor analysis with FDR correction, phase comparisons and clustering. |
| [`cosinor_bootstrap_FDR.py`](#cosinor_bootstrap_fdrpy) | Bootstrap-based cosinor analysis with empirical significance testing and exploratory visualisations. |

---

# Common Methodology

All scripts fit the classical single-component 24-hour cosinor model:

\[
Y(t) = M + A \cos\left(\frac{2\pi t}{24} - \phi\right)
\]

where

- \(M\): MESOR (midline estimate)
- \(A\): amplitude
- \(\phi\): acrophase
- \(t\): Zeitgeber Time (ZT)

Model fitting is performed using nonlinear least squares (`scipy.optimize.curve_fit`).

All scripts report rhythm parameters including amplitude, acrophase and model goodness-of-fit (R²). Depending on the workflow, rhythmic significance is assessed using either the classical zero-amplitude F-test or a bootstrap permutation test, with optional Benjamini-Hochberg False Discovery Rate (FDR) correction.

---

# Common Input Format

Input datasets are tab-separated (`.tsv`) files.

Accepted time columns:

- `ZT`
- `Time`

Accepted expression columns:

- `expression ratio (goi/hk)`
- `fold change (goi/hk)`

Column matching is case-insensitive where applicable.

Missing values are ignored automatically.

---

# Installation

Requires **Python 3.9+**.

Install dependencies with:

```bash
pip install numpy pandas scipy matplotlib statsmodels scikit-learn networkx
```

If NumPy 2 causes compatibility issues:

```bash
pip install "numpy<2"
```

---

# `multi_gene_cosinor_plot.py`

Performs batch 24-hour cosinor analysis for multiple genes contained within a single directory.

### Features

- Fits independent cosinor models for each gene
- Computes MESOR, amplitude, amplitude SE, acrophase and R²
- Performs a zero-amplitude F-test
- Generates publication-quality grayscale figures
- Exports a summary CSV containing rhythm parameters

### Input

One `.tsv` file per gene.

Example:

```text
Bmal1.tsv
Per2.tsv
Cry1.tsv
```

### Configuration

```python
folder = "/path/to/tsv/files/"
show_points = True
period = 24
save_summary = True
save_figures = True
```

### Usage

```bash
python3 multi_gene_cosinor_plot.py
```

or

```bash
python3 multi_gene_cosinor_plot.py /path/to/data/
```

### Outputs

| File | Description |
|------|-------------|
| `<gene>_cosinor.pdf` | Publication-quality cosinor plot |
| `cosinor_summary.csv` | Rhythm parameters for all genes |

### Summary Table

| Column | Meaning |
|---------|---------|
| MESOR | Midline estimate |
| Amplitude | Oscillation amplitude |
| Amplitude_SE | Standard error of amplitude |
| Acrophase_radians | Phase (radians) |
| Acrophase_hours | Phase (hours) |
| R2 | Goodness of fit |
| ZeroAmp_pvalue | Zero-amplitude F-test p-value |

**Notes**

- Uses a classical analytical zero-amplitude F-test.
- Figures share a common Y-axis scale.
- Raw observations are shown as open circles; mean ± SD as filled circles.

---

# `cosinor_batch_by_tissue.py`

Performs independent cosinor analysis across multiple tissues for multiple genes and compares rhythmic behaviour between tissues.

### Features

- Fits independent cosinor models for every gene–tissue combination
- Performs zero-amplitude testing with Benjamini-Hochberg FDR correction
- Identifies rhythmic profiles using FDR and minimum R² thresholds
- Computes pairwise tissue phase differences
- Performs circular K-means clustering of phases
- Generates comparative multi-tissue plots and amplitude heatmaps

### Directory Structure

```text
Circadian_Gene_Expression_Data/

├── Liver_data/
├── Heart_data/
├── Kidney_data/
└── ...
```

Each tissue directory contains one `.tsv` file per gene.

### Configuration

```python
folder = "/path/to/Circadian_Gene_Expression_Data/"

show_points = True
period = 24

save_summary = True
save_figures = True

min_r2_for_rhythmic = 0.20
```

### Usage

```bash
python3 cosinor_batch_by_tissue.py
```

### Outputs

| File | Description |
|------|-------------|
| `<gene>_cosinor.pdf` | Comparative multi-tissue cosinor plot |
| `amplitude_heatmap.pdf` | Gene × tissue amplitude heatmap |
| `cosinor_summary.csv` | Rhythm statistics with FDR correction |
| `phase_shift_stats.csv` | Pairwise tissue phase differences |

### Summary Table

| Column | Meaning |
|---------|---------|
| Gene | Gene name |
| Tissue | Tissue |
| Amplitude | Oscillation amplitude |
| Acrophase | Peak time (hours) |
| p | Zero-amplitude p-value |
| FDR | Benjamini-Hochberg adjusted p-value |
| R2 | Goodness of fit |
| Rhythmic | FDR < 0.05 and R² exceeds threshold |
| PhaseCluster | Circular K-means cluster |

### Phase Shift Table

| Column | Meaning |
|---------|---------|
| Gene | Gene name |
| Tissue1 | First tissue |
| Tissue2 | Second tissue |
| PhaseDiff_hours | Circular phase difference (hours) |

**Notes**

- Tissues are analysed independently.
- Shared Y-axis limits allow direct comparison between genes.
- Colours distinguish tissues within each figure.

---

# `cosinor_bootstrap_FDR.py`

Performs multi-tissue cosinor analysis using empirical bootstrap permutation testing instead of analytical significance testing.

### Features

- Fits cosinor models independently for every gene–tissue combination
- Estimates empirical p-values using bootstrap permutation
- Computes bootstrap confidence intervals for amplitude
- Applies Benjamini-Hochberg FDR correction
- Produces exploratory visualisations of rhythmic structure
- Generates hierarchical clustering and tissue similarity networks

### Bootstrap Method

For each dataset:

1. Fit the observed cosinor model.
2. Randomly permute expression values while preserving sampling times.
3. Refit the cosinor model to each permutation.
4. Estimate the empirical p-value from the null amplitude distribution.

Default: **300 bootstrap iterations** per dataset.

### Configuration

```python
folder = "/path/to/Circadian_Gene_Expression_Data/"

period = 24
BOOTSTRAPS = 300

save_figures = True
```

### Usage

```bash
python3 cosinor_bootstrap_FDR.py
```

### Outputs

| File | Description |
|------|-------------|
| `panel_A.pdf` | Gene × tissue amplitude heatmap |
| `panel_B.pdf` | Polar phase distributions |
| `panel_C.pdf` | Hierarchical clustering dendrogram |
| `panel_D.pdf` | Tissue similarity network |
| `cosinor_bootstrap_FDR_summary.csv` | Bootstrap statistics with FDR correction |

### Summary Table

| Column | Meaning |
|---------|---------|
| Gene | Gene name |
| Tissue | Tissue |
| Amplitude | Oscillation amplitude |
| Phase | Peak time (hours) |
| p | Empirical bootstrap p-value |
| CI_low | Lower bootstrap confidence limit |
| CI_high | Upper bootstrap confidence limit |
| R2 | Goodness of fit |
| FDR | Benjamini-Hochberg adjusted p-value |
| Rhythmic | FDR < 0.05 |

**Notes**

- Uses bootstrap permutation testing rather than the analytical F-test.
- Bootstrap inference is generally more robust for small sample sizes and non-normal data.
- Hierarchical clustering combines amplitude with circular phase information.
- Tissue similarity networks are exploratory and based on phase correlations.

---

# Interpretation of Results

| Metric | Meaning |
|--------|---------|
| **MESOR** | Average expression level |
| **Amplitude** | Oscillation strength |
| **Acrophase / Phase** | Predicted peak time |
| **R²** | Goodness of fit |
| **p-value** | Statistical evidence for rhythmicity |
| **FDR** | Multiple-testing corrected significance |

Higher amplitudes and better model fits generally provide stronger evidence for rhythmic expression.

---

# General Notes

- All analyses assume a single-component 24-hour cosinor model.
- Time points do not need to be equally spaced.
- Missing observations are ignored automatically.
- Figures are exported as publication-quality PDF files.

---

# References

Nelson W., Tong Y.L., Lee J.K., Halberg F. (1979). *Methods for cosinor-rhythmometry.* Chronobiologia, 6(4), 305–323.

Refinetti R., Cornélissen G., Halberg F. (2007). *Procedures for numerical analysis of circadian rhythms.* Biological Rhythm Research, 38(4), 275–325.

Benjamini Y., Hochberg Y. (1995). *Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing.* Journal of the Royal Statistical Society Series B, 57(1), 289–300.

Efron B., Tibshirani R.J. (1993). *An Introduction to the Bootstrap.* Chapman & Hall.
