import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import f
import glob
import os
import sys

# ===============================================================
# SETTINGS
# ===============================================================

import sys
folder = sys.argv[1] if len(sys.argv) > 1 else ""                  # folder with gene .tsv files
show_points = True                                                 # show raw data
period = 24                                                        # circadian period
save_summary = True                                                # save cosinor summary CSV
save_figures = True                                                # export each plot as PDF

# ===============================================================
# FIGURE STYLE — publication-quality grayscale
# ===============================================================

plt.style.use("default")
plt.rcParams.update({
    "font.size": 14,
    "axes.labelsize": 16,
    "axes.titlesize": 18,
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
    "lines.linewidth": 2,
    "figure.dpi": 150
})

# ===============================================================
# COSINOR FUNCTION
# ===============================================================

def cosinor(t, mesor, amp, acro):
    return mesor + amp * np.cos(2 * np.pi * t / period - acro)

def compute_r2(y_obs, y_pred):
    ss_res = np.sum((y_obs - y_pred)**2)
    ss_tot = np.sum((y_obs - np.mean(y_obs))**2)
    return 1 - ss_res/ss_tot

# ===============================================================
# LOAD ALL GENES + GLOBAL Y LIMITS
# ===============================================================

files = glob.glob(os.path.join(folder, "*.tsv"))
datasets = {}
ymin, ymax = float("inf"), float("-inf")

for file_path in files:
    df = pd.read_csv(file_path, sep="\t")
    df = df.rename(columns={
        "ZT": "Time",
        "expression ratio (goi/hk)": "Value"
    })
    datasets[file_path] = df
    ymin = min(ymin, df["Value"].min())
    ymax = max(ymax, df["Value"].max())

padding = (ymax - ymin) * 0.1
ymin -= padding
ymax += padding

# ===============================================================
# COSINOR SUMMARY TABLE
# ===============================================================

summary_rows = []

# ===============================================================
# PROCESS EACH GENE
# ===============================================================

for file_path, df in datasets.items():
    gene_name = os.path.splitext(os.path.basename(file_path))[0]

    times = df["Time"].values
    values = df["Value"].values
    N = len(values)

    # Initial guesses
    guess_mesor = np.mean(values)
    guess_amp   = (values.max() - values.min()) / 2
    guess_acro  = 0

    # Fit
    popt, pcov = curve_fit(
        cosinor,
        times,
        values,
        p0=[guess_mesor, guess_amp, guess_acro]
    )

    mesor, amp, acro = popt

    # Standard error of amplitude
    amp_se = np.sqrt(pcov[1, 1])

    # R²
    predicted = cosinor(times, *popt)
    r2 = compute_r2(values, predicted)

    # Acrophase hours
    acro_hours = ((acro * period) / (2 * np.pi)) % 24

    # ===============================================================
    # ZERO-AMPLITUDE TEST (p-value)
    # ===============================================================

    F_stat = (amp / amp_se)**2 / 2
    df1 = 2
    df2 = N - 3
    p_value = 1 - f.cdf(F_stat, df1, df2)

    # Save stats
    summary_rows.append({
        "Gene": gene_name,
        "MESOR": mesor,
        "Amplitude": amp,
        "Amplitude_SE": amp_se,
        "Acrophase_radians": acro,
        "Acrophase_hours": acro_hours,
        "R2": r2,
        "ZeroAmp_pvalue": p_value
    })

    # ===============================================================
    # COSINOR CURVE (0–24 h)
    # ===============================================================

    t_fit = np.linspace(0, 24, 500)
    y_fit = cosinor(t_fit, *popt)

    # Summary stats
    sumdf = df.groupby("Time").agg(
        mean=("Value", "mean"),
        sd=("Value", "std")
    ).reset_index()

    # ===============================================================
    # PLOT — Publication-quality grayscale
    # ===============================================================

    fig, ax = plt.subplots(figsize=(9, 5))

    # Raw points (open circles)
    if show_points:
        ax.scatter(
            df["Time"], df["Value"],
            facecolors='none',
            edgecolors='black',
            s=70,
            label="Raw data"
        )

    # Mean ± SD (closed circles)
    ax.errorbar(
        sumdf["Time"],
        sumdf["mean"],
        yerr=sumdf["sd"],
        fmt='o',
        color='black',
        capsize=4,
        markersize=8,
        label="Mean ± SD"
    )

    # Cosinor curve
    ax.plot(t_fit, y_fit, color='black', label="Cosinor fit")

    # X-axis ticks including ZT24
    xticks_sorted = sorted(set(df["Time"].unique()).union({24}))
    ax.set_xticks(xticks_sorted)

    ax.set_title(f"{gene_name} — Cosinor Analysis")
    ax.set_xlabel("ZT (hours)")
    ax.set_ylabel("Expression Ratio (GOI/HK)")
    ax.set_ylim(ymin, ymax)
    ax.legend()
    fig.tight_layout()

    # ===============================================================
    # SAVE FIGURE AS PDF
    # ===============================================================

    if save_figures:
        out_pdf = f"{folder}{gene_name}_cosinor.pdf"
        fig.savefig(out_pdf, bbox_inches="tight")
        print(f"Saved figure: {out_pdf}")

    plt.close(fig)

# ===============================================================
# SAVE COSINOR SUMMARY
# ===============================================================

summary_df = pd.DataFrame(summary_rows)
print("\n=== COSINOR SUMMARY ===")
print(summary_df)

if save_summary:
    out_csv = f"{folder}cosinor_summary.csv"
    summary_df.to_csv(out_csv, index=False)
    print(f"\nSaved cosinor summary to: {out_csv}")

