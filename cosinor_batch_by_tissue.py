import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import f
from statsmodels.stats.multitest import multipletests
from sklearn.cluster import KMeans
import glob
import os

# ===============================================================
# SETTINGS
# ===============================================================

folder = '/PATH/TO/YOUR/Circadian_Gene_Expression_Data/'

show_points = True
period = 24
save_summary = True
save_figures = True

# analysis toggles
min_r2_for_rhythmic = 0.2

# ===============================================================
# COSINOR MODEL
# ===============================================================

def cosinor(t, mesor, amp, acro):
    return mesor + amp * np.cos(2 * np.pi * t / period - acro)

def compute_r2(y_obs, y_pred):
    ss_res = np.sum((y_obs - y_pred) ** 2)
    ss_tot = np.sum((y_obs - np.mean(y_obs)) ** 2)
    return 1 - ss_res / ss_tot

# ===============================================================
# PARSING
# ===============================================================

def extract_gene(file_path):
    return os.path.splitext(os.path.basename(file_path))[0].rsplit("_", 1)[0]

def extract_tissue(file_path):
    return os.path.basename(os.path.dirname(file_path)).replace("_data", "")

# ===============================================================
# LOAD DATA
# ===============================================================

files = glob.glob(os.path.join(folder, "*_data", "*.tsv"))

data = {}
all_fits = []   # for stats aggregation

ymin, ymax = float("inf"), float("-inf")

for file_path in sorted(files):

    gene = extract_gene(file_path)
    tissue = extract_tissue(file_path)

    df = pd.read_csv(file_path, sep="\t")
    df.columns = df.columns.str.strip()

    time_col = next((c for c in df.columns if c.lower() in ["zt", "time"]), None)

    value_col = None
    for key in ["expression ratio (goi/hk)", "fold change (goi/hk)"]:
        for c in df.columns:
            if c.strip().lower() == key:
                value_col = c
                break
        if value_col:
            break

    df = df.rename(columns={time_col: "Time", value_col: "Value"})
    df = df.dropna(subset=["Time", "Value"])

    t = df["Time"].values
    y = df["Value"].values

    p0 = [np.mean(y), (np.max(y)-np.min(y))/2, 0]

    popt, pcov = curve_fit(cosinor, t, y, p0=p0, maxfev=10000)

    mesor, amp, acro = popt
    pred = cosinor(t, *popt)

    r2 = compute_r2(y, pred)

    amp_se = np.sqrt(max(pcov[1,1], 0))

    N = len(y)
    if amp_se > 0 and N > 3:
        F = (amp / amp_se) ** 2 / 2
        p_value = 1 - f.cdf(F, 2, N - 3)
    else:
        p_value = np.nan

    acro_hours = (acro * period / (2*np.pi)) % 24

    if gene not in data:
        data[gene] = {}

    data[gene][tissue] = {
        "df": df,
        "mesor": mesor,
        "amp": amp,
        "acro": acro,
        "acro_hours": acro_hours,
        "r2": r2,
        "p": p_value
    }

    all_fits.append({
        "Gene": gene,
        "Tissue": tissue,
        "Amplitude": amp,
        "Acrophase": acro_hours,
        "p": p_value,
        "r2": r2
    })

    ymin, ymax = min(ymin, y.min()), max(ymax, y.max())

# ===============================================================
# FDR CORRECTION
# ===============================================================

fit_df = pd.DataFrame(all_fits)
fit_df["FDR"] = multipletests(fit_df["p"].fillna(1), method="fdr_bh")[1]
fit_df["Rhythmic"] = (fit_df["FDR"] < 0.05) & (fit_df["r2"] > min_r2_for_rhythmic)

# ===============================================================
# AMPLITUDE HEATMAP (GENE × TISSUE)
# ===============================================================

amp_matrix = fit_df.pivot(index="Gene", columns="Tissue", values="Amplitude")

plt.figure(figsize=(8, 5))
plt.imshow(amp_matrix, aspect="auto")
plt.xticks(range(len(amp_matrix.columns)), amp_matrix.columns, rotation=45)
plt.yticks(range(len(amp_matrix.index)), amp_matrix.index)
plt.title("Amplitude Heatmap (Gene × Tissue)")
plt.colorbar(label="Amplitude")
plt.tight_layout()

if save_figures:
    plt.savefig(os.path.join(folder, "amplitude_heatmap.pdf"))
plt.close()

# ===============================================================
# PHASE CLUSTERING (CIRCULAR EMBEDDING)
# ===============================================================

phase_vecs = np.array([
    [np.cos(np.deg2rad(r["Acrophase"])), np.sin(np.deg2rad(r["Acrophase"]))]
    for r in all_fits
])

k = min(3, len(phase_vecs))
kmeans = KMeans(n_clusters=k, n_init=10, random_state=0).fit(phase_vecs)
fit_df["PhaseCluster"] = kmeans.labels_

# ===============================================================
# PHASE SHIFT COMPARISON FUNCTION
# ===============================================================

def circ_diff(a, b):
    return np.angle(np.exp(1j*(np.deg2rad(a)-np.deg2rad(b)))) * 180/np.pi

phase_stats = []

for gene in data:
    tissues = list(data[gene].keys())
    for i in range(len(tissues)):
        for j in range(i+1, len(tissues)):

            t1, t2 = tissues[i], tissues[j]

            if t1 in data[gene] and t2 in data[gene]:
                a1 = data[gene][t1]["acro_hours"]
                a2 = data[gene][t2]["acro_hours"]

                phase_stats.append({
                    "Gene": gene,
                    "Tissue1": t1,
                    "Tissue2": t2,
                    "PhaseDiff_hours": circ_diff(a1*15, a2*15)/15
                })

phase_df = pd.DataFrame(phase_stats)

# ===============================================================
# PLOT EACH GENE
# ===============================================================

tissues = sorted({t for g in data for t in data[g].keys()})
colors = plt.cm.tab10(np.linspace(0, 1, len(tissues)))

for gene, tissue_dict in data.items():

    fig, ax = plt.subplots(figsize=(9, 5.5))
    xticks = {0, 24}

    for color, tissue in zip(colors, tissues):

        if tissue not in tissue_dict:
            continue

        df = tissue_dict[tissue]["df"]
        t = df["Time"].values
        y = df["Value"].values

        xticks.update(t)

        p0 = [np.mean(y), (np.max(y)-np.min(y))/2, 0]
        popt, _ = curve_fit(cosinor, t, y, p0=p0)

        t_fit = np.linspace(0, 24, 500)
        y_fit = cosinor(t_fit, *popt)

        grouped = df.groupby("Time").agg(mean=("Value","mean"), sd=("Value","std")).reset_index()

        if show_points:
            ax.scatter(df["Time"], df["Value"], facecolors="none", edgecolors=[color], s=35)

        ax.errorbar(grouped["Time"], grouped["mean"], yerr=grouped["sd"],
                    fmt="o", color=color, capsize=3, markersize=5)

        ax.plot(t_fit, y_fit, color=color, label=tissue)

    ax.set_xlim(0, 24)
    ax.set_ylim(ymin, ymax)
    ax.set_xticks(sorted(xticks))

    ax.set_title(f"{gene} — Cosinor Across Tissues")
    ax.set_xlabel("ZT (hours)")
    ax.set_ylabel("Expression Ratio (GOI/HK)")
    ax.legend()

    if save_figures:
        plt.savefig(os.path.join(folder, f"{gene}_cosinor.pdf"), bbox_inches="tight")

    plt.close()

# ===============================================================
# SAVE OUTPUTS
# ===============================================================

summary_df = fit_df

if save_summary:
    summary_df.to_csv(os.path.join(folder, "cosinor_summary.csv"), index=False)
    phase_df.to_csv(os.path.join(folder, "phase_shift_stats.csv"), index=False)

print(summary_df.head())
print(phase_df.head())
