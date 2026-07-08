import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from statsmodels.stats.multitest import multipletests
from scipy.cluster.hierarchy import linkage, dendrogram
import networkx as nx
import glob
import os

# ===============================================================
# SETTINGS
# ===============================================================

folder = '/PATH/TO/YOUR/Circadian_Gene_Expression_Data/'
period = 24
BOOTSTRAPS = 300
save_figures = True

# ===============================================================
# MODEL
# ===============================================================

def cosinor(t, mesor, amp, acro):
    return mesor + amp * np.cos(2*np.pi*t/period - acro)

def r2(y, yhat):
    ss_res = np.sum((y - yhat)**2)
    ss_tot = np.sum((y - np.mean(y))**2)
    return 1 - ss_res/ss_tot

# ===============================================================
# PARSING
# ===============================================================

def gene(fp):
    return os.path.splitext(os.path.basename(fp))[0].rsplit("_",1)[0]

def tissue(fp):
    return os.path.basename(os.path.dirname(fp)).replace("_data","")

# ===============================================================
# BOOTSTRAP COSINOR
# ===============================================================

def bootstrap_fit(t, y):

    t = np.array(t)
    y = np.array(y)

    try:
        popt,_ = curve_fit(cosinor, t, y,
                           p0=[np.mean(y),(y.max()-y.min())/2,0],
                           maxfev=10000)
    except:
        return None

    mesor, amp, acro = popt

    yhat = cosinor(t,*popt)

    # bootstrap null distribution
    boot_amps = []

    for _ in range(BOOTSTRAPS):
        y_perm = np.random.permutation(y)

        try:
            popt_b,_ = curve_fit(cosinor, t, y_perm,
                                 p0=[np.mean(y),(y.max()-y.min())/2,0],
                                 maxfev=3000)
            boot_amps.append(popt_b[1])
        except:
            continue

    boot_amps = np.array(boot_amps)

    if len(boot_amps) < 10:
        return None

    p_val = np.mean(np.abs(boot_amps) >= np.abs(amp))
    ci = np.percentile(boot_amps,[2.5,97.5])

    return mesor, amp, acro, p_val, ci, r2(y,yhat)

# ===============================================================
# LOAD + FIT
# ===============================================================

files = glob.glob(os.path.join(folder,"*_data","*.tsv"))

raw = {}
results = []

for fp in files:

    g = gene(fp)
    tiss = tissue(fp)

    df = pd.read_csv(fp, sep="\t")
    df.columns = df.columns.str.strip()

    tcol = next(c for c in df.columns if c.lower() in ["zt","time"])
    vcol = next(c for c in df.columns if "goi/hk" in c.lower())

    df = df.rename(columns={tcol:"t", vcol:"y"}).dropna()

    out = bootstrap_fit(df["t"].values, df["y"].values)

    if out is None:
        continue

    mesor, amp, acro, p_val, ci, r2v = out

    phase_h = (acro*period/(2*np.pi)) % 24

    results.append({
        "Gene": g,
        "Tissue": tiss,
        "Amplitude": amp,
        "Phase": phase_h,
        "p": p_val,
        "CI_low": ci[0],
        "CI_high": ci[1],
        "R2": r2v
    })

    if g not in raw:
        raw[g] = {}
    raw[g][tiss] = df

# ===============================================================
# FDR (RESTORED)
# ===============================================================

df = pd.DataFrame(results)

df["p"] = df["p"].fillna(1.0)
df["FDR"] = multipletests(df["p"], method="fdr_bh")[1]
df["Rhythmic"] = df["FDR"] < 0.05

# ===============================================================
# MATRICES
# ===============================================================

amp_mat = df.pivot(index="Gene", columns="Tissue", values="Amplitude")
phase_mat = df.pivot(index="Gene", columns="Tissue", values="Phase")

# ===============================================================
# PANEL A — AMPLITUDE HEATMAP
# ===============================================================

plt.figure(figsize=(8,5))
plt.imshow(amp_mat, aspect="auto")
plt.xticks(range(len(amp_mat.columns)), amp_mat.columns, rotation=45)
plt.yticks(range(len(amp_mat.index)), amp_mat.index)
plt.title("Panel A — Amplitude (FDR-corrected)")
plt.colorbar()
plt.tight_layout()

if save_figures:
    plt.savefig(os.path.join(folder,"panel_A.pdf"))
plt.close()

# ===============================================================
# PANEL B — PHASE DISTRIBUTION
# ===============================================================

fig = plt.figure(figsize=(12,4))
genes = list(raw.keys())

for i,g in enumerate(genes):
    ax = plt.subplot(1,len(genes),i+1, projection="polar")

    phases = df[df.Gene==g]["Phase"].values
    ax.hist(np.deg2rad(phases), bins=6, alpha=0.7)
    ax.set_title(g)

plt.tight_layout()

if save_figures:
    plt.savefig(os.path.join(folder,"panel_B.pdf"))
plt.close()

# ===============================================================
# PANEL C — CLUSTERING
# ===============================================================

feat = np.hstack([
    amp_mat.fillna(0).values,
    np.sin(np.deg2rad(phase_mat.fillna(0).values)),
    np.cos(np.deg2rad(phase_mat.fillna(0).values))
])

Z = linkage(feat, method="ward")

plt.figure(figsize=(8,5))
dendrogram(Z, labels=amp_mat.index)
plt.title("Panel C — Gene Clustering")
plt.xticks(rotation=45)
plt.tight_layout()

if save_figures:
    plt.savefig(os.path.join(folder,"panel_C.pdf"))
plt.close()

# ===============================================================
# PANEL D — NETWORK
# ===============================================================

tissues = df.Tissue.unique()

G = nx.Graph()
G.add_nodes_from(tissues)

for i,t1 in enumerate(tissues):
    for j,t2 in enumerate(tissues):
        if i>=j: continue

        v1 = df[df.Tissue==t1]["Phase"].values
        v2 = df[df.Tissue==t2]["Phase"].values

        m = min(len(v1),len(v2))
        if m < 3:
            continue

        corr = np.corrcoef(v1[:m], v2[:m])[0,1]

        if np.isfinite(corr):
            G.add_edge(t1,t2,weight=corr)

plt.figure(figsize=(5,5))
pos = nx.spring_layout(G, seed=1)

weights = [G[u][v]["weight"] for u,v in G.edges]
nx.draw(G, pos, with_labels=True, width=np.abs(weights)*3)

plt.title("Panel D — Tissue Network")

if save_figures:
    plt.savefig(os.path.join(folder,"panel_D.pdf"))
plt.close()

# ===============================================================
# SAVE
# ===============================================================

df.to_csv(os.path.join(folder,"cosinor_bootstrap_FDR_summary.csv"), index=False)

print(df.head())
