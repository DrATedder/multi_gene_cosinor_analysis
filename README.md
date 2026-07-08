# Multi-gene Cosinor Analysis & Plotting

* [`multi_gene_cosinor_plot.py`](https://github.com/DrATedder/multi_gene_cosinor_analysis/blob/main/README.md#multi_gene_cosinor_plotpy)
* [`cosinor_batch_by_tissue.py`](https://github.com/DrATedder/multi_gene_cosinor_analysis/blob/main/README.md#cosinor_batch_by_tissuepy)
* [`cosinor_bootstrap_FDR.py`](https://github.com/DrATedder/multi_gene_cosinor_analysis/blob/main/README.md#cosinor_bootstrap_fdrpy)

# `multi_gene_cosinor_plot.py`

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

* Bmal1.tsv
* Per2.tsv
* Cry1.tsv


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
pip3 install "numpy<2"
```

---

## Usage

1. Place all .tsv files in your target directory.

2. Edit the folder variable in the script to point to that directory.

3. Run:

```bash
python3 multi_gene_cosinor_plot.py
```
Alternatively, you can leave the `folder` variable unchanged in the script (i.e. ignore step 2 above), and simply give the destination folder when calling the script:
```bash
python3 multi_gene_cosinor_plot.py /path/to/target_directory/
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

# `cosinor_batch_by_tissue.py`

`cosinor_batch_by_tissue.py` performs batch 24-hour single-component cosinor analysis across multiple tissues for multiple genes. In addition to fitting independent cosinor models for every gene–tissue combination, the script performs multiple-testing correction, identifies rhythmic profiles, compares acrophase shifts between tissues, clusters rhythmic phases, and generates publication-quality comparative figures.

For each gene–tissue dataset, the script:

- Fits a standard 24-h cosinor model
- Computes MESOR, amplitude, acrophase, and R²
- Performs a zero-amplitude (F-test) for rhythmic significance
- Applies Benjamini-Hochberg False Discovery Rate (FDR) correction across all tests
- Classifies rhythmic profiles using both FDR and minimum R² thresholds
- Compares phase differences between tissues for each gene
- Performs phase clustering using circular K-means
- Generates comparative multi-tissue cosinor plots
- Produces an amplitude heatmap summarising rhythmic strength across tissues
- Exports summary tables for downstream analysis

---

## Method Summary

Each gene within each tissue is modeled using the classical single-component cosinor equation:

\[
Y(t) = M + A \cos\left(\frac{2\pi t}{24} - \phi\right)
\]

where

- \(M\): MESOR (midline estimate)
- \(A\): amplitude
- \(\phi\): acrophase
- \(t\): Zeitgeber Time (ZT)
- 24 h: assumed circadian period

The model is fit using nonlinear least squares (`scipy.optimize.curve_fit`).

Model quality is quantified using the coefficient of determination (R²).

Rhythmic significance is assessed using the classical zero-amplitude F-test:

\[
F = \frac{A^{2}/SE(A)^2}{2}
\]

with

- \(df_1 = 2\)
- \(df_2 = N-3\)

P-values are subsequently corrected for multiple testing using the Benjamini-Hochberg False Discovery Rate (FDR) procedure.

A gene–tissue combination is classified as **rhythmic** when:

- FDR < 0.05
- R² exceeds a user-defined threshold (default = 0.20)

---

## Script Description

**Filename:** `multi_tissue_cosinor_analysis.py`

This script:

1. Searches all tissue subdirectories (`*_data`) within a specified folder.
2. Reads every `.tsv` expression file.
3. Automatically extracts gene names from filenames.
4. Automatically extracts tissue names from directory names.
5. Fits an independent 24-hour cosinor model for every gene–tissue combination.
6. Computes:
   - MESOR
   - Amplitude
   - Acrophase (radians and hours)
   - R²
   - Zero-amplitude p-value
7. Applies Benjamini-Hochberg FDR correction across all fitted models.
8. Identifies rhythmic gene–tissue combinations using FDR and minimum R² thresholds.
9. Computes pairwise phase shifts between tissues for each gene.
10. Performs circular phase clustering using K-means.
11. Generates publication-quality comparative plots showing cosinor fits for all tissues for each gene.
12. Produces an amplitude heatmap summarising oscillation amplitudes across tissues.
13. Exports summary tables for rhythmic statistics and tissue phase differences.

---

## Outputs

The script produces:

- **Per-gene comparative cosinor figures** (`<gene>_cosinor.pdf`)
- **Amplitude heatmap** (`amplitude_heatmap.pdf`)
- **Cosinor summary table** (`cosinor_summary.csv`)
- **Pairwise tissue phase-shift table** (`phase_shift_stats.csv`)
- Console output showing summary statistics

Each gene figure displays:

- Raw observations (open circles)
- Mean ± SD at each time point
- Independent fitted cosinor curves for every tissue
- Shared Y-axis scaling across all genes for direct comparison

---

## Input Directory Structure

The script expects one folder per tissue.

Example:

```text
Circadian_Gene_Expression_Data/

├── Liver_data/
│   ├── Bmal1_expression.tsv
│   ├── Per2_expression.tsv
│   └── Cry1_expression.tsv
│
├── Heart_data/
│   ├── Bmal1_expression.tsv
│   ├── Per2_expression.tsv
│   └── Cry1_expression.tsv
│
├── Kidney_data/
│   ├── Bmal1_expression.tsv
│   ├── Per2_expression.tsv
│   └── Cry1_expression.tsv
```

Gene names are extracted from filenames by removing the final underscore suffix.

For example:

```text
Bmal1_expression.tsv
```

becomes

```text
Bmal1
```

Tissue names are extracted automatically from directory names:

```text
Heart_data
```

becomes

```text
Heart
```

---

## Input Data Format

Each `.tsv` file should contain either:

| ZT | expression ratio (goi/hk) |
|----|---------------------------|

or

| Time | fold change (goi/hk) |
|------|-----------------------|

Recognised time columns:

- `ZT`
- `Time`

Recognised expression columns:

- `expression ratio (goi/hk)`
- `fold change (goi/hk)`

Column names are matched automatically (case insensitive).

Missing values are automatically ignored.

---

## Configuration Parameters

Edit the **SETTINGS** section near the top of the script:

```python
folder = "/path/to/Circadian_Gene_Expression_Data/"

show_points = True
period = 24

save_summary = True
save_figures = True

min_r2_for_rhythmic = 0.20
```

### Parameter Description

| Parameter | Description |
|-----------|-------------|
| `folder` | Parent directory containing tissue folders |
| `show_points` | Plot individual observations |
| `period` | Assumed circadian period (hours) |
| `save_summary` | Export CSV summary tables |
| `save_figures` | Export PDF figures |
| `min_r2_for_rhythmic` | Minimum R² required for rhythmic classification |

---

## Dependencies

Requires Python 3.9+ and:

```python
numpy
pandas
scipy
matplotlib
statsmodels
scikit-learn
```

---

## Install

```bash
pip3 install numpy pandas scipy matplotlib statsmodels scikit-learn
```

If `NumPy` version conflicts occur:

```bash
pip3 install "numpy<2"
```

---

## Usage

1. Arrange your data into one folder per tissue.

2. Edit the `folder` variable in the script.

3. Run:

```bash
python3 multi_tissue_cosinor_analysis.py
```

Alternatively, modify the script to accept a directory argument if desired.

The script will automatically:

- Fit cosinor models for every gene in every tissue
- Correct p-values using FDR
- Identify rhythmic profiles
- Compare tissue-specific phases
- Cluster rhythmic phases
- Generate comparative publication-quality figures
- Produce summary tables

---

## Outputs

| File | Description |
|------|-------------|
| `<gene>_cosinor.pdf` | Multi-tissue cosinor comparison figure |
| `amplitude_heatmap.pdf` | Heatmap of amplitudes across genes and tissues |
| `cosinor_summary.csv` | Complete cosinor statistics with FDR correction |
| `phase_shift_stats.csv` | Pairwise tissue acrophase differences |

---

## Cosinor Summary Columns

| Column | Meaning |
|---------|---------|
| Gene | Gene name |
| Tissue | Tissue |
| Amplitude | Estimated oscillation amplitude |
| Acrophase | Peak time (hours) |
| p | Zero-amplitude F-test p-value |
| FDR | Benjamini-Hochberg adjusted p-value |
| R2 | Model goodness of fit |
| Rhythmic | TRUE if FDR < 0.05 and R² exceeds threshold |
| PhaseCluster | Circular K-means phase cluster |

---

## Phase Shift Table

| Column | Meaning |
|---------|---------|
| Gene | Gene name |
| Tissue1 | First tissue |
| Tissue2 | Second tissue |
| PhaseDiff_hours | Circular acrophase difference (hours) |

---

## Interpretation of Results

**MESOR**

Average expression level across the circadian cycle.

**Amplitude**

Strength of rhythmic oscillation.

**Acrophase**

Predicted peak expression time.

**R²**

Proportion of variance explained by the fitted cosinor model.

**p-value**

Evidence against the null hypothesis of zero amplitude.

**FDR**

Multiple-testing corrected significance.

**Rhythmic**

Indicates statistically significant rhythmicity after FDR correction while also requiring adequate model fit.

**PhaseCluster**

Groups gene–tissue combinations exhibiting similar circadian phases.

**PhaseDiff_hours**

Estimated phase shift between tissues for the same gene, accounting for circular time.

---

## Notes

- Assumes a single-component 24-hour cosinor model.
- Each tissue is fit independently.
- Time points do not need to be equally spaced.
- Missing observations are automatically ignored.
- Shared Y-axis limits allow direct visual comparison between genes.
- Colours distinguish tissues within each figure.
- Phase clustering is performed using circular embedding followed by K-means clustering.
- Amplitude heatmaps provide a global overview of tissue-specific rhythmic strength.

---

## References

Nelson W., Tong Y.L., Lee J.K., Halberg F. (1979). *Methods for cosinor-rhythmometry.* Chronobiologia, 6(4), 305–323.

Refinetti R., Cornélissen G., Halberg F. (2007). *Procedures for numerical analysis of circadian rhythms.* Biological Rhythm Research, 38(4), 275–325.

Benjamini Y., Hochberg Y. (1995). *Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing.* Journal of the Royal Statistical Society Series B, 57(1), 289–300.

# `cosinor_bootstrap_FDR.py`

`cosinor_bootstrap_FDR.py` performs robust 24-hour single-component cosinor analysis across multiple tissues using a non-parametric bootstrap approach. Rather than relying solely on analytical significance tests, the script estimates empirical p-values from bootstrap resampling, applies False Discovery Rate (FDR) correction, and generates a suite of exploratory visualisations summarising rhythmic properties across tissues.

For each gene–tissue dataset, the script:

- Fits a standard 24-hour cosinor model
- Computes MESOR, amplitude, acrophase, and R²
- Estimates empirical significance using bootstrap permutation testing
- Calculates bootstrap confidence intervals for amplitude
- Applies Benjamini-Hochberg False Discovery Rate (FDR) correction
- Identifies rhythmic gene–tissue combinations
- Generates amplitude heatmaps
- Visualises phase distributions using polar histograms
- Clusters genes based on rhythmic similarity
- Constructs tissue similarity networks based on phase correlations
- Exports publication-quality figures and summary tables

---

## Method Summary

Each gene within each tissue is modeled using the classical single-component cosinor model:

\[
Y(t) = M + A \cos\left(\frac{2\pi t}{24} - \phi\right)
\]

where

- \(M\): MESOR (midline estimate)
- \(A\): amplitude
- \(\phi\): acrophase
- \(t\): Zeitgeber Time (ZT)
- 24 h: assumed circadian period

Model fitting is performed using nonlinear least squares (`scipy.optimize.curve_fit`).

Goodness of fit is quantified using the coefficient of determination (R²).

Unlike traditional cosinor analysis, rhythmic significance is estimated empirically using bootstrap permutation testing.

For each fitted dataset:

1. Expression values are randomly permuted while preserving sampling times.
2. The cosinor model is refitted to each permuted dataset.
3. A null distribution of amplitudes is generated.
4. The empirical p-value is calculated as the proportion of bootstrap amplitudes greater than or equal to the observed amplitude.

The default analysis performs **300 bootstrap iterations** per dataset.

Following bootstrap analysis, p-values are corrected using the Benjamini-Hochberg False Discovery Rate (FDR) procedure.

---

## Script Description

**Filename:** `multi_tissue_bootstrap_cosinor_analysis.py`

This script:

1. Searches all tissue subdirectories (`*_data`) within a specified folder.
2. Reads every `.tsv` expression dataset.
3. Automatically extracts gene names from filenames.
4. Automatically extracts tissue names from directory names.
5. Fits a 24-hour cosinor model independently for every gene–tissue combination.
6. Computes:
   - MESOR
   - Amplitude
   - Acrophase (hours)
   - R²
7. Performs bootstrap permutation testing to estimate:
   - empirical p-values
   - 95% bootstrap confidence intervals for amplitude
8. Applies Benjamini-Hochberg FDR correction.
9. Generates four publication-quality summary panels:
   - Panel A: Amplitude heatmap
   - Panel B: Polar phase distributions
   - Panel C: Hierarchical clustering dendrogram
   - Panel D: Tissue similarity network
10. Exports a comprehensive summary table.

---

## Outputs

The script produces:

- **Panel A:** `panel_A.pdf` — Amplitude heatmap
- **Panel B:** `panel_B.pdf` — Polar phase distributions
- **Panel C:** `panel_C.pdf` — Hierarchical clustering dendrogram
- **Panel D:** `panel_D.pdf` — Tissue similarity network
- **Bootstrap summary table:** `cosinor_bootstrap_FDR_summary.csv`

Console output displays the first rows of the summary table.

---

## Input Directory Structure

The script expects one directory per tissue.

Example:

```text
Circadian_Gene_Expression_Data/

├── Liver_data/
│   ├── Bmal1_expression.tsv
│   ├── Per2_expression.tsv
│   └── Cry1_expression.tsv
│
├── Heart_data/
│   ├── Bmal1_expression.tsv
│   ├── Per2_expression.tsv
│   └── Cry1_expression.tsv
│
├── Kidney_data/
│   ├── Bmal1_expression.tsv
│   ├── Per2_expression.tsv
│   └── Cry1_expression.tsv
```

Gene names are extracted from filenames by removing the final underscore suffix.

For example:

```text
Bmal1_expression.tsv
```

becomes

```text
Bmal1
```

Tissue names are extracted automatically from directory names.

For example:

```text
Liver_data
```

becomes

```text
Liver
```

---

## Input Data Format

Each `.tsv` file should contain either:

| ZT | expression ratio (goi/hk) |
|----|---------------------------|

or

| Time | expression ratio (goi/hk) |
|------|---------------------------|

Recognised time columns:

- `ZT`
- `Time`

Expression columns are automatically identified by the presence of:

- `GOI/HK`

Missing observations are automatically ignored.

---

## Configuration Parameters

Edit the **SETTINGS** section near the beginning of the script:

```python
folder = "/path/to/Circadian_Gene_Expression_Data/"

period = 24

BOOTSTRAPS = 300

save_figures = True
```

### Parameter Description

| Parameter | Description |
|-----------|-------------|
| `folder` | Parent directory containing tissue folders |
| `period` | Assumed circadian period (hours) |
| `BOOTSTRAPS` | Number of bootstrap permutations |
| `save_figures` | Save publication-quality PDF figures |

---

## Dependencies

Requires Python 3.9+ and:

```python
numpy
pandas
scipy
matplotlib
statsmodels
networkx
```

---

## Install

```bash
pip3 install numpy pandas scipy matplotlib statsmodels networkx
```

If `NumPy` version conflicts occur:

```bash
pip3 install "numpy<2"
```

---

## Usage

1. Arrange datasets into one folder per tissue.

2. Edit the `folder` variable.

3. Run:

```bash
python3 multi_tissue_bootstrap_cosinor_analysis.py
```

The script will automatically:

- Fit bootstrap cosinor models
- Estimate empirical significance
- Correct p-values using FDR
- Generate four publication-quality figures
- Produce a comprehensive summary table

---

## Outputs

| File | Description |
|------|-------------|
| `panel_A.pdf` | Heatmap of oscillation amplitudes |
| `panel_B.pdf` | Polar histograms of tissue phases |
| `panel_C.pdf` | Hierarchical clustering of genes |
| `panel_D.pdf` | Tissue phase-correlation network |
| `cosinor_bootstrap_FDR_summary.csv` | Bootstrap statistics with FDR correction |

---

## Summary Table Columns

| Column | Meaning |
|---------|---------|
| Gene | Gene name |
| Tissue | Tissue |
| Amplitude | Estimated oscillation amplitude |
| Phase | Estimated acrophase (hours) |
| p | Empirical bootstrap p-value |
| CI_low | Lower 95% bootstrap confidence limit |
| CI_high | Upper 95% bootstrap confidence limit |
| R2 | Model goodness of fit |
| FDR | Benjamini-Hochberg adjusted p-value |
| Rhythmic | TRUE if FDR < 0.05 |

---

## Figure Descriptions

### Panel A — Amplitude Heatmap

Displays oscillation amplitudes for every gene across all tissues, providing a global overview of rhythmic strength.

### Panel B — Polar Phase Distributions

Shows polar histograms of acrophase values for each gene, illustrating phase dispersion across tissues.

### Panel C — Hierarchical Clustering

Clusters genes according to combined amplitude and phase characteristics using Ward's hierarchical clustering method.

Phase information is represented using sine and cosine transformations to preserve circular relationships.

### Panel D — Tissue Similarity Network

Constructs a network of tissues connected by correlations between gene acrophases.

Edge widths are proportional to the magnitude of phase correlation, providing a visual summary of tissue synchrony.

---

## Interpretation of Results

**MESOR**

Average expression level across the circadian cycle.

**Amplitude**

Strength of rhythmic oscillation.

**Phase**

Predicted peak expression time (hours).

**Bootstrap p-value**

Empirical probability of observing an amplitude at least as large under the null hypothesis.

**Confidence Interval**

Bootstrap-derived 95% interval for the null amplitude distribution.

**FDR**

Multiple-testing corrected significance.

**R²**

Goodness of fit of the cosinor model.

**Rhythmic**

Indicates statistically significant rhythmicity after FDR correction.

---

## Notes

- Assumes a single-component 24-hour cosinor model.
- Uses empirical bootstrap permutation testing instead of analytical significance tests.
- Bootstrap p-values are generally more robust for small sample sizes or non-normal data.
- Each tissue is analysed independently.
- Missing observations are automatically ignored.
- Phase-based analyses account for the circular nature of circadian time.
- Hierarchical clustering combines amplitude with circular phase information.
- Tissue similarity networks are exploratory and based on phase correlations between tissues.

---

## References

Nelson W., Tong Y.L., Lee J.K., Halberg F. (1979). *Methods for cosinor-rhythmometry.* Chronobiologia, 6(4), 305–323.

Refinetti R., Cornélissen G., Halberg F. (2007). *Procedures for numerical analysis of circadian rhythms.* Biological Rhythm Research, 38(4), 275–325.

Benjamini Y., Hochberg Y. (1995). *Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing.* Journal of the Royal Statistical Society Series B, 57(1), 289–300.

Efron B., Tibshirani R.J. (1993). *An Introduction to the Bootstrap.* Chapman & Hall.
