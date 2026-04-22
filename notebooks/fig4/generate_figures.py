#!/usr/bin/env python3
"""
Merged figure generation for IDH-astro telomere manuscript.

Input files (relative to project root):
  data/v5/all_df_aug_final_18_03_2026.txt        ��� Master dataset: per-sample, per-arm TL with SV/CNV/ecDNA/HAR counts
  data/v5/telogator_telomere_data_merged_v3.csv  — Telogator2 allele-level telomere length estimates (TL_p75)
  data/v5/clinical_data_merged.txt               — Clinical data (grade, survival)

Output figures (written to figures/v5/):
  Main:
    tl_zscore_distribution.pdf/.png              — Figure 4A: density + pie chart of TL z-score categories
    figure4_tl_association.pdf/.png               — Figure 4 composite (panels A–E)
  Supplementary:
    tl_median_clustered_heatmap.pdf/.png          — Clustered heatmap of median TL per patient x arm
    tl_by_event_all_measures.pdf/.png             — TL distributions by event type (TL Min/Max/Median)

Intermediate output (saved for inspection):
  figures/v5/dichotomized_TL_association_results.tsv — Mixed-model results for all event–TL associations

Usage:
  python scripts/v5/generate_figures.py
"""

import pandas as pd
import numpy as np
from scipy import stats
from scipy.cluster.hierarchy import linkage, fcluster, leaves_list
from scipy.spatial.distance import pdist
from scipy.stats import gaussian_kde, fisher_exact as fisher_exact_test
from statsmodels.stats.multitest import multipletests
import statsmodels.formula.api as smf
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import matplotlib.gridspec as gridspec
import matplotlib.font_manager as fm
import seaborn as sns
from pathlib import Path
import time
import warnings
warnings.filterwarnings("ignore")

# =====================================================================
#  Configuration
# =====================================================================
ROOTDIR = Path(__file__).resolve().parents[2]
DATADIR = ROOTDIR / "data" / "v5"
FIGDIR = ROOTDIR / "figures" / "v5"
FIGDIR.mkdir(parents=True, exist_ok=True)

# Register fonts
FONT_DIR = Path.home() / "Documents" / "Barthel-Custom-Powerpoint-Theme" / "fonts"
if FONT_DIR.exists():
    for f in FONT_DIR.glob("BasicSans-*.otf"):
        fm.fontManager.addfont(str(f))

plt.rcParams["font.family"] = "Arial"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["font.size"] = 12
plt.rcParams["axes.labelsize"] = 14
plt.rcParams["axes.titlesize"] = 15
plt.rcParams["xtick.labelsize"] = 11
plt.rcParams["ytick.labelsize"] = 11

# Colors
CLR_BRK = "#2E8B57"
CLR_BRK_LIGHT = "#A8D5BA"
CLR_AMP = "#C0392B"
CLR_AMP_LIGHT = "#F1948A"
CLR_NEUTRAL = "#B0B0B0"
CLR_GRAY = "#4A4A4A"
CLR_NA = "#D3D3D3"

CLR_CRIT_SHORT = "#1B5E37"
CLR_SHORT = "#6DB88A"
CLR_NORMAL = "#B0B0B0"
CLR_LONG = "#E8A09A"
CLR_VERY_LONG = "#C0392B"

Z_CUTS = [-1.5, -0.5, 0.5, 1.5]
CAT_NAMES = ["Critical Short", "Short", "Normal", "Long", "Very Long"]
CAT_COLORS = [CLR_CRIT_SHORT, CLR_SHORT, CLR_NORMAL, CLR_LONG, CLR_VERY_LONG]

MIN_OBS = 5


# =====================================================================
#  Fast random-intercept mixed model solver
# =====================================================================
def fast_lmer(y, X_fixed, groups):
    """
    REML solution for: y ~ X_fixed + (1|group)
    X_fixed: 2D array of fixed-effect predictors (WITHOUT intercept).
    Returns: betas, se, z, p arrays for each predictor.
    """
    y = np.asarray(y, dtype=np.float64)
    X_fixed = np.asarray(X_fixed, dtype=np.float64)
    if X_fixed.ndim == 1:
        X_fixed = X_fixed.reshape(-1, 1)

    n = len(y)
    p = X_fixed.shape[1]
    X = np.column_stack([np.ones(n), X_fixed])

    unique_groups, group_idx = np.unique(groups, return_inverse=True)
    n_groups = len(unique_groups)
    ni = np.bincount(group_idx)

    try:
        XtX_inv = np.linalg.inv(X.T @ X)
    except np.linalg.LinAlgError:
        nans = np.full(p, np.nan)
        return nans, nans, nans, nans
    beta = XtX_inv @ (X.T @ y)
    resid = y - X @ beta
    sigma2_e = np.var(resid)
    sigma2_u = max(sigma2_e * 0.5, 1e-8)

    for _ in range(200):
        s2e_old, s2u_old = sigma2_e, sigma2_u

        theta = sigma2_u / (sigma2_u + sigma2_e / ni)
        resid = y - X @ beta
        gs = np.bincount(group_idx, weights=resid, minlength=n_groups)
        gm = gs / ni
        u_hat = theta * gm

        sigma2_u = max(np.mean(u_hat**2 + sigma2_u * (1 - theta)), 1e-12)
        sigma2_e = max(np.sum((resid - u_hat[group_idx])**2) / n, 1e-12)

        d_inv = 1.0 / (sigma2_e + ni * sigma2_u)

        k = X.shape[1]
        H = np.zeros((k, k))
        g_vec = np.zeros(k)
        for g in range(n_groups):
            mask = group_idx == g
            Xg, yg = X[mask], y[mask]
            sd = 1.0 / sigma2_e
            so = sigma2_u * d_inv[g] / sigma2_e
            Xs = Xg.sum(axis=0)
            H += sd * (Xg.T @ Xg) - so * np.outer(Xs, Xs)
            g_vec += sd * (Xg.T @ yg) - so * Xs * yg.sum()

        try:
            H_inv = np.linalg.inv(H)
        except np.linalg.LinAlgError:
            nans = np.full(p, np.nan)
            return nans, nans, nans, nans
        beta = H_inv @ g_vec

        if abs(sigma2_e - s2e_old) < 1e-8 * max(abs(sigma2_e), 1) and \
           abs(sigma2_u - s2u_old) < 1e-8 * max(abs(sigma2_u), 1):
            break

    betas = beta[1:]
    ses = np.sqrt(np.diag(H_inv)[1:])
    bad = ses < 1e-15
    z_vals = np.where(bad, np.nan, betas / ses)
    p_vals = np.where(bad, np.nan, 2 * stats.norm.sf(np.abs(z_vals)))

    return betas, ses, z_vals, p_vals


# =====================================================================
#  Plotting helpers
# =====================================================================
def violin_box_sina(ax, data_list, positions, colors, width=0.35):
    parts = ax.violinplot(data_list, positions=positions, showextrema=False)
    for pc, color in zip(parts["bodies"], colors):
        pc.set_facecolor(color)
        pc.set_edgecolor(CLR_GRAY)
        pc.set_alpha(0.7)
    bp = ax.boxplot(data_list, positions=positions, widths=width * 0.5,
                     patch_artist=True, showfliers=False, zorder=4)
    for patch in bp["boxes"]:
        patch.set_facecolor("white")
        patch.set_edgecolor(CLR_GRAY)
    for element in ["whiskers", "caps", "medians"]:
        for line in bp[element]:
            line.set_color(CLR_GRAY)
    for i, (vals, color) in enumerate(zip(data_list, colors)):
        jitter = np.random.normal(0, 0.05, len(vals))
        ax.scatter(np.full(len(vals), positions[i]) + jitter, vals,
                   s=5, color=color, alpha=0.25, zorder=3)


def add_bracket(ax, x1, x2, y, text, fontsize=11):
    ax.plot([x1, x1, x2, x2], [y * 0.98, y, y, y * 0.98], color=CLR_GRAY, linewidth=1)
    ax.text((x1 + x2) / 2, y * 1.01, text, ha="center", va="bottom",
            fontsize=fontsize, color=CLR_GRAY)


def star(p):
    return "***" if p < 0.001 else ("**" if p < 0.01 else ("*" if p < 0.05 else "ns"))


def assign_category(z):
    if np.isnan(z):
        return "NA"
    if z < Z_CUTS[0]:
        return CAT_NAMES[0]
    elif z < Z_CUTS[1]:
        return CAT_NAMES[1]
    elif z < Z_CUTS[2]:
        return CAT_NAMES[2]
    elif z < Z_CUTS[3]:
        return CAT_NAMES[3]
    else:
        return CAT_NAMES[4]


# =====================================================================
#  Step 1: Load all data
# =====================================================================
print("=" * 70)
print("  IDH-Astro Telomere Manuscript — Figure Generation (v5)")
print("=" * 70)

df = pd.read_csv(DATADIR / "all_df_aug_final_18_03_2026.txt", sep="\t")
print(f"Master data: {len(df)} rows, {df['Sample_ID'].nunique()} samples")

telo = pd.read_csv(DATADIR / "telogator_telomere_data_merged_v3.csv", sep="\t")
telo["Sample_ID"] = telo["Sample"].str.replace(r"_run\d+$", "", regex=True)
telo.rename(columns={"chr": "chr_arm"}, inplace=True)
telo["TL_kb"] = telo["TL_p75"] / 1000.0

chr_arms = [
    "chr1p", "chr1q", "chr2p", "chr2q", "chr3p", "chr3q", "chr4p", "chr4q",
    "chr5p", "chr5q", "chr6p", "chr6q", "chr7p", "chr7q", "chr8p", "chr8q",
    "chr9p", "chr9q", "chr10p", "chr10q", "chr11p", "chr11q", "chr12p", "chr12q",
    "chr13q", "chr14q", "chr15q", "chr16p", "chr16q", "chr17p", "chr17q",
    "chr18p", "chr18q", "chr19p", "chr19q", "chr20p", "chr20q", "chr21q", "chr22q",
]
properties = [
    "sv_count_DUP", "sv_count_DEL", "Chromothripsis",
    "Single_Chr_Complex_Rearrangement", "Templated_Insertion_Chains",
    "Multiple_Chr_Complex_Rearrangement", "Chromoplexy",
    "auto_complex_inv_Count", "auto_Templated_ins_Count",
    "auto_dup_inv_segment_Count", "auto_foldback_Count",
    "auto_inv_tra_Count", "auto_Templated_ins_inv_Count",
    "ecdna_simple_circle", "ecdna_multi_region_intra_chr",
]

grade_df = pd.read_csv(DATADIR / "clinical_data_merged.txt", sep="\t")
grade_map = dict(zip(grade_df["ID"].astype(str), grade_df["Grade"].astype(str)))
print(f"Clinical data: grade available for {len(grade_map)} samples")


# =====================================================================
#  Step 2: Compute min/max/median TL from telogator
# =====================================================================
print("\n--- Computing allelic TL measures ---")

# Drop interstitial telomere reads
n_before = len(telo)
telo_filt = telo[~telo["allele_id"].astype(str).str.endswith("i")].copy()
print(f"Dropped {n_before - len(telo_filt)} interstitial reads ({n_before} -> {len(telo_filt)})")

telo_agg = (telo_filt.groupby(["Sample_ID", "chr_arm"])["TL_kb"]
            .agg(tl_min="min", tl_max="max", tl_median="median")
            .reset_index())

# Also compute bp-scale aggregation for dichotomized analysis compatibility
telo["TL_p75_kb"] = telo["TL_p75"] / 1000.0
telo_agg_all = (telo.groupby(["Sample_ID", "chr_arm"])["TL_p75_kb"]
                .agg(min_TL_p75_kb="min", max_TL_p75_kb="max", median_TL_p75_kb="median")
                .reset_index())

# Merge into master data
m = df.merge(telo_agg, on=["Sample_ID", "chr_arm"], how="left")
m = m[m["chr_arm"].isin(chr_arms)].copy()
print(f"Merged data: {len(m)} rows, {m['Sample_ID'].nunique()} samples, {len(chr_arms)} arms")

# Also build df_dichot (for dichotomized analysis, uses unfiltered telogator)
df_dichot = df.merge(telo_agg_all, on=["Sample_ID", "chr_arm"], how="left")
df_dichot = df_dichot[df_dichot["chr_arm"].isin(chr_arms)].copy()

# Set Age=0 to NA
n_age_zero = (df_dichot["Age"] == 0).sum()
df_dichot.loc[df_dichot["Age"] == 0, "Age"] = np.nan


# =====================================================================
#  Step 3: Dichotomized TL association analysis (for Panel B forest plot)
# =====================================================================
print("\n" + "=" * 70)
print("  Dichotomized TL Association Analysis")
print("=" * 70)

# Verify properties exist
missing = [p for p in properties if p not in df_dichot.columns]
if missing:
    print(f"WARNING: Properties not in data: {missing}")
    properties = [p for p in properties if p in df_dichot.columns]

nice_names = {}
for col in properties:
    name = col
    name = name.replace("auto_", "").replace("_Count", "").replace("_Rate_per_Mb", " Rate")
    name = name.replace("sv_count_", "SV ").replace("sv_size_sum_", "SV Size ")
    name = name.replace("aneu_", "Aneu ").replace("has_ecdna", "Has ecDNA")
    name = name.replace("HAR_ecDNA_count", "HAR ecDNA")
    name = name.replace("ecdna_", "ecDNA ")
    name = name.replace("_", " ")
    nice_names[col] = name

tl_measures_dichot = {
    "final_average_TL_p75_kb": "TL Avg",
    "min_TL_p75_kb": "TL Min",
    "max_TL_p75_kb": "TL Max",
    "median_TL_p75_kb": "TL Median",
}

groups = df_dichot["Sample_ID"].values

t0 = time.time()
results = []

for prop in properties:
    bin_vals = (df_dichot[prop] > 0).astype(float).values
    n_present = int(bin_vals.sum())
    n_absent = int((bin_vals == 0).sum())

    if n_present < MIN_OBS or n_absent < MIN_OBS:
        continue

    for tl_col, tl_label in tl_measures_dichot.items():
        y = df_dichot[tl_col].values
        valid = ~np.isnan(y)
        if valid.sum() < 10:
            continue

        y_v = y[valid]
        bin_v = bin_vals[valid]
        groups_v = groups[valid]

        row = {
            "property": prop,
            "property_name": nice_names[prop],
            "TL_measure": tl_label,
            "TL_col": tl_col,
            "n_absent": int((bin_v == 0).sum()),
            "n_present": int((bin_v == 1).sum()),
            "mean_absent": y_v[bin_v == 0].mean(),
            "mean_present": y_v[bin_v == 1].mean(),
        }

        # Mixed model (no age)
        betas, ses, z_vals, p_vals = fast_lmer(y_v, bin_v, groups_v)
        row["beta_lmer"] = betas[0]
        row["se_lmer"] = ses[0]
        row["z_lmer"] = z_vals[0]
        row["p_lmer"] = p_vals[0]

        results.append(row)

dichot_df = pd.DataFrame(results)
print(f"{len(dichot_df)} tests in {time.time() - t0:.1f}s")

# FDR correction
if len(dichot_df) > 0 and "p_lmer" in dichot_df.columns:
    valid_mask = dichot_df["p_lmer"].notna()
    if valid_mask.sum() > 0:
        _, fdr_p, _, _ = multipletests(dichot_df.loc[valid_mask, "p_lmer"], method="fdr_bh")
        dichot_df.loc[valid_mask, "p_fdr_lmer"] = fdr_p

# Save intermediate results
dichot_path = FIGDIR / "dichotomized_TL_association_results.tsv"
dichot_df.to_csv(dichot_path, sep="\t", index=False)
print(f"Saved {dichot_path}")


# =====================================================================
#  Figure 1: TL z-score distribution (density + pie chart)
# =====================================================================
print("\n" + "=" * 70)
print("  Figure: tl_zscore_distribution")
print("=" * 70)

# Build TL matrix for heatmap and z-score computation
tl_matrix = m.pivot_table(index="Sample_ID", columns="chr_arm", values="tl_median")
tl_matrix = tl_matrix[chr_arms]

# Impute NaN with column median
raw_matrix = tl_matrix.copy()
for col in tl_matrix.columns:
    tl_matrix[col] = tl_matrix[col].fillna(tl_matrix[col].median())

global_mean = tl_matrix.values.ravel().mean()
global_sd = tl_matrix.values.ravel().std()
z_matrix = (tl_matrix - global_mean) / global_sd
print(f"Global z-score: mean TL = {global_mean:.2f} kb, SD = {global_sd:.2f} kb")

z_vals_all = z_matrix.values.ravel()
z_vals_valid = z_vals_all[~np.isnan(z_vals_all)]
z_cats = [assign_category(z) for z in z_vals_all]
cat_counts = pd.Series(z_cats).value_counts()
total_n = len(z_vals_all)

for cat in CAT_NAMES:
    n = cat_counts.get(cat, 0)
    print(f"  {cat}: {n} ({n/total_n*100:.1f}%)")

fig_dist, (ax_pie, ax_dens) = plt.subplots(1, 2, figsize=(14, 5),
                                             gridspec_kw={"width_ratios": [1, 1.5]})
fig_dist.suptitle("Telomere length categories across chromosome arms",
                  fontsize=15, fontweight="bold", y=1.02)

# Pie chart
pie_order = CAT_NAMES.copy()
na_count = cat_counts.get("NA", 0)
if na_count > 0:
    pie_order.append("NA")
pie_sizes = [cat_counts.get(c, 0) for c in pie_order]
pie_colors_list = CAT_COLORS.copy()
if na_count > 0:
    pie_colors_list.append(CLR_NA)
pie_labels = [f"{c}\n{cat_counts.get(c,0)/total_n*100:.1f}%" for c in pie_order]

wedges, texts = ax_pie.pie(pie_sizes, labels=pie_labels, colors=pie_colors_list,
                            startangle=90, textprops={"fontsize": 9})
ax_pie.set_title("Z-score distribution", fontsize=13, fontweight="bold")

# Density plot with colored regions
kde = gaussian_kde(z_vals_valid, bw_method=0.3)
x_grid = np.linspace(z_vals_valid.min() - 0.5, z_vals_valid.max() + 0.5, 500)
y_kde = kde(x_grid)

boundaries = [x_grid[0]] + Z_CUTS + [x_grid[-1]]
for i, (lo, hi) in enumerate(zip(boundaries[:-1], boundaries[1:])):
    mask = (x_grid >= lo) & (x_grid <= hi)
    ax_dens.fill_between(x_grid[mask], y_kde[mask], alpha=0.7, color=CAT_COLORS[i])
    ax_dens.plot(x_grid[mask], y_kde[mask], color=CLR_GRAY, linewidth=0.8)

for cut in Z_CUTS:
    ax_dens.axvline(cut, color="#666666", linestyle="--", linewidth=0.8)

cat_centers = [
    (boundaries[0] + boundaries[1]) / 2,
    (boundaries[1] + boundaries[2]) / 2,
    (boundaries[2] + boundaries[3]) / 2,
    (boundaries[3] + boundaries[4]) / 2,
    (boundaries[4] + boundaries[5]) / 2,
]
cat_short_names = ["Critical\nshort", "Short", "Normal", "Long", "Very\nlong"]
y_top = y_kde.max() * 1.05
for cx, label in zip(cat_centers, cat_short_names):
    ax_dens.text(cx, y_top, label, ha="center", va="bottom", fontsize=9, color=CLR_GRAY)

ax_dens.set_xlabel("Z-scored telomere length", fontsize=13)
ax_dens.set_ylabel("Density", fontsize=13)
ax_dens.spines["top"].set_visible(False)
ax_dens.spines["right"].set_visible(False)
ax_dens.set_ylim(0, y_top * 1.25)

legend_patches = [mpatches.Patch(facecolor=c, edgecolor="gray", label=n)
                  for c, n in zip(CAT_COLORS, CAT_NAMES)]
ax_dens.legend(handles=legend_patches, loc="upper right", fontsize=9,
               frameon=True, edgecolor="gray", title="TL Category", title_fontsize=10)

plt.tight_layout()
fig_dist.savefig(FIGDIR / "tl_zscore_distribution.pdf", dpi=300, bbox_inches="tight")
fig_dist.savefig(FIGDIR / "tl_zscore_distribution.png", dpi=300, bbox_inches="tight")
plt.close(fig_dist)
print("Saved tl_zscore_distribution.pdf/.png")


# =====================================================================
#  Figure 2: Clustered heatmap (tl_median_clustered_heatmap)
# =====================================================================
print("\n" + "=" * 70)
print("  Figure: tl_median_clustered_heatmap")
print("=" * 70)

samples = tl_matrix.index.tolist()

# Row annotations
grade_colors_map = {"2": "#60A5FA", "3": "#F07167"}
sample_grade = pd.Series({s: grade_map.get(s, "N/A") for s in samples}, name="Grade")
grade_palette = {**grade_colors_map, "N/A": "#E0E0E0"}
grade_colors = sample_grade.map(grade_palette)

har_per_sample = m.groupby("Sample_ID")["HAR_ecDNA_count"].sum()
har_series = pd.Series({s: har_per_sample.get(s, 0) for s in samples}, name="HAR/ecDNA")
mean_tl = tl_matrix.mean(axis=1)
mean_tl.name = "Mean TL"

har_binary = har_series.map(lambda x: "#C0392B" if x > 0 else "#E0E0E0")

def continuous_colors(series, cmap_name, vmin=None, vmax=None):
    cmap = plt.get_cmap(cmap_name)
    v0 = vmin if vmin is not None else series.min()
    v1 = vmax if vmax is not None else series.max()
    norm = mcolors.Normalize(vmin=v0, vmax=v1)
    return series.map(lambda x: mcolors.to_hex(cmap(norm(x))))

tl_colors = continuous_colors(mean_tl, "YlOrRd")

# Event matrices for markers
m["has_foldback"] = ((m["auto_foldback_Count"] > 0) | (m["Foldback_Translocation"] > 0)).astype(int)
m["has_chromothripsis"] = (m["Chromothripsis"] > 0).astype(int)
m["has_ecdna_multi"] = (m["ecdna_multi_region_intra_chr"] > 0).astype(int)
m["has_dup_inv"] = (m["auto_dup_inv_segment_Count"] > 0).astype(int)
m["has_recip_inv"] = (m["sv_count_INV"] > 0).astype(int)

event_matrices = {}
for event_col, label in [("has_foldback", "Foldback"),
                          ("has_chromothripsis", "Chromothripsis"),
                          ("has_ecdna_multi", "ecDNA multi-region"),
                          ("has_dup_inv", "Dup inv"),
                          ("has_recip_inv", "Reciprocal inv")]:
    ev_mat = m.pivot_table(index="Sample_ID", columns="chr_arm", values=event_col, fill_value=0)
    ev_mat = ev_mat.reindex(index=tl_matrix.index, columns=chr_arms, fill_value=0)
    event_matrices[label] = ev_mat

# Clustering
row_linkage = linkage(pdist(z_matrix.values, metric="euclidean"),
                      method="ward", optimal_ordering=True)
col_linkage = linkage(pdist(z_matrix.values.T, metric="euclidean"),
                      method="ward", optimal_ordering=True)

heatmap_cmap = mcolors.LinearSegmentedColormap.from_list(
    "tl_heatmap", [CLR_BRK, "#FFFFFF", CLR_AMP])

g = sns.clustermap(
    z_matrix,
    row_linkage=row_linkage,
    col_linkage=col_linkage,
    cmap=heatmap_cmap,
    vmin=-3, vmax=3,
    center=0,
    figsize=(20, 10),
    linewidths=0.3,
    linecolor="#E8E8E8",
    dendrogram_ratio=(0.12, 0.08),
    cbar_pos=(0.02, 0.82, 0.03, 0.15),
    cbar_kws={"label": "TL z-score", "orientation": "vertical"},
    xticklabels=True,
    yticklabels=True,
)

g.ax_heatmap.set_xlabel("Chromosome Arm", fontsize=14)
g.ax_heatmap.set_ylabel("")
g.ax_heatmap.tick_params(axis="x", rotation=90, labelsize=9)
g.ax_heatmap.tick_params(axis="y", labelsize=10)

g.fig.suptitle("Median Telomere Length per Patient and Chromosome Arm\n"
               "(global z-score, Ward clustering)",
               fontsize=16, fontweight="bold", y=1.02)

# Overlay event markers
row_order = g.dendrogram_row.reordered_ind
col_order = g.dendrogram_col.reordered_ind

event_markers = [
    ("Foldback", "v", "#000000", 5),
    ("Chromothripsis", "D", "#FBBF24", 5),
    ("ecDNA multi-region", "o", "#40D392", 5),
    ("Dup inv", "s", "#C4A9E8", 5),
    ("Reciprocal inv", "^", "#3B82F6", 5),
]
offsets = [(-0.2, -0.15), (0.2, -0.15), (0, 0.18), (-0.2, 0.15), (0.2, 0.15)]

for (label, marker, color, size), (dx, dy) in zip(event_markers, offsets):
    ev_mat = event_matrices[label]
    ev_reordered = ev_mat.iloc[row_order, col_order]
    for i in range(ev_reordered.shape[0]):
        for j in range(ev_reordered.shape[1]):
            if ev_reordered.iloc[i, j] > 0:
                g.ax_heatmap.scatter(j + 0.5 + dx, i + 0.5 + dy,
                                      marker=marker, s=size * 6,
                                      color=color, edgecolors="white",
                                      linewidths=0.4, zorder=5)

# Legend
legend_elements = []
legend_elements.append(mpatches.Patch(facecolor="#60A5FA", edgecolor="gray", label="Grade 2"))
legend_elements.append(mpatches.Patch(facecolor="#F07167", edgecolor="gray", label="Grade 3"))
legend_elements.append(mpatches.Patch(facecolor="#E0E0E0", edgecolor="gray", label="Grade N/A"))
legend_elements.append(mpatches.Patch(facecolor="none", edgecolor="none", label=""))
legend_elements.append(mlines.Line2D([], [], marker="v", color="none", markerfacecolor="#000000",
                                      markeredgecolor="white", markersize=7, label="Foldback"))
legend_elements.append(mlines.Line2D([], [], marker="D", color="none", markerfacecolor="#FBBF24",
                                      markeredgecolor="white", markersize=7, label="Chromothripsis"))
legend_elements.append(mlines.Line2D([], [], marker="o", color="none", markerfacecolor="#40D392",
                                      markeredgecolor="white", markersize=7, label="ecDNA multi-region"))
legend_elements.append(mlines.Line2D([], [], marker="s", color="none", markerfacecolor="#C4A9E8",
                                      markeredgecolor="white", markersize=7, label="Dup inv"))
legend_elements.append(mlines.Line2D([], [], marker="^", color="none", markerfacecolor="#3B82F6",
                                      markeredgecolor="white", markersize=7, label="Reciprocal inv"))

g.ax_heatmap.legend(handles=legend_elements, loc="upper left",
                     bbox_to_anchor=(1.18, 1.0), fontsize=10,
                     frameon=True, edgecolor="gray", title="Annotations",
                     title_fontsize=11)

# Cluster assignments
clusters = fcluster(row_linkage, t=2, criterion="maxclust")
cluster_means = {}
for c in [1, 2]:
    mask_c = clusters == c
    cluster_means[c] = mean_tl[z_matrix.index[mask_c]].mean() if mask_c.any() else 0

label_map = {1: "A", 2: "B"} if cluster_means[1] >= cluster_means[2] else {1: "B", 2: "A"}
cluster_labels = pd.Series([label_map[c] for c in clusters], index=z_matrix.index, name="Cluster")

print(f"\n{'Sample_ID':<12s}  {'Cluster':>7s}  {'Mean TL':>8s}  {'Grade':>5s}  {'HAR/ecDNA':>10s}")
print("-" * 55)
for s in sorted(samples):
    print(f"{s:<12s}  {cluster_labels[s]:>7s}  {mean_tl[s]:8.2f}  "
          f"{sample_grade[s]:>5s}  {har_series[s]:10.0f}")

for cl in ["A", "B"]:
    members = cluster_labels[cluster_labels == cl].index
    mtl_cl = mean_tl[members].mean()
    grades = [sample_grade[s] for s in members]
    print(f"Cluster {cl}: n={len(members)}, mean TL={mtl_cl:.2f} kb, "
          f"Grade 2={grades.count('2')}, Grade 3={grades.count('3')}")

g.savefig(FIGDIR / "tl_median_clustered_heatmap.pdf", dpi=300, bbox_inches="tight")
g.savefig(FIGDIR / "tl_median_clustered_heatmap.png", dpi=300, bbox_inches="tight")
plt.close()
print("Saved tl_median_clustered_heatmap.pdf/.png")


# =====================================================================
#  Figure 3: TL by event — all measures (tl_by_event_all_measures)
# =====================================================================
print("\n" + "=" * 70)
print("  Figure: tl_by_event_all_measures")
print("=" * 70)

N_PERMS = 10000
rng = np.random.default_rng(42)

event_binary_cols = {
    "Foldback": "has_foldback",
    "Chromothripsis": "has_chromothripsis",
    "ecDNA multi-region": "has_ecdna_multi",
    "Dup inv": "has_dup_inv",
    "Reciprocal inv": "has_recip_inv",
}
event_sig_type = {
    "Foldback": "breakage",
    "Chromothripsis": "breakage",
    "ecDNA multi-region": "amplification",
    "Dup inv": "amplification",
    "Reciprocal inv": "breakage",
}

tl_measures_event = [
    ("tl_min", "TL Min (shortest allele)"),
    ("tl_max", "TL Max (longest allele)"),
    ("tl_median", "TL Median"),
]

# Build min/max matrices
tl_min_matrix = m.pivot_table(index="Sample_ID", columns="chr_arm", values="tl_min")
tl_min_matrix = tl_min_matrix.reindex(columns=chr_arms)
for col in tl_min_matrix.columns:
    tl_min_matrix[col] = tl_min_matrix[col].fillna(tl_min_matrix[col].median())

tl_max_matrix = m.pivot_table(index="Sample_ID", columns="chr_arm", values="tl_max")
tl_max_matrix = tl_max_matrix.reindex(columns=chr_arms)
for col in tl_max_matrix.columns:
    tl_max_matrix[col] = tl_max_matrix[col].fillna(tl_max_matrix[col].median())

measure_matrices = {"tl_min": tl_min_matrix, "tl_max": tl_max_matrix, "tl_median": tl_matrix}

# Global z-score each matrix
measure_matrices_z = {}
for k, mat in measure_matrices.items():
    flat = mat.values.ravel()
    g_mean = np.nanmean(flat)
    g_sd = np.nanstd(flat)
    measure_matrices_z[k] = (mat - g_mean) / g_sd

event_cols_to_merge = list(event_binary_cols.values())

fig2, all_axes = plt.subplots(3, 5, figsize=(20, 14), sharey="row",
                               gridspec_kw={"hspace": 0.45, "wspace": 0.15})
fig2.suptitle("Global z-scored TL Distribution by Event Presence\n"
              "(permutation test, 10k permutations)",
              fontsize=15, fontweight="bold", y=1.01)

for row_idx, (tl_col, tl_label) in enumerate(tl_measures_event):
    axes_row = all_axes[row_idx]
    mat = measure_matrices_z[tl_col]

    tl_long = mat.reset_index().melt(id_vars="Sample_ID", var_name="chr_arm", value_name="tl_val")
    tl_long["log2_tl"] = tl_long["tl_val"]
    tl_long = tl_long.merge(m[["Sample_ID", "chr_arm"] + event_cols_to_merge],
                             on=["Sample_ID", "chr_arm"], how="left")

    all_vals = tl_long["log2_tl"].values
    cohort_med = np.nanmedian(all_vals)
    valid_idx = np.where(~np.isnan(all_vals))[0]

    print(f"\n  --- {tl_label} (median z={cohort_med:.2f}) ---")

    for col_idx, (event_name, event_col) in enumerate(event_binary_cols.items()):
        ax = axes_row[col_idx]
        sig = event_sig_type[event_name]
        color = CLR_BRK if sig == "breakage" else CLR_AMP

        event_mask = tl_long[event_col].values == 1
        no_vals = tl_long.loc[~event_mask, "log2_tl"].dropna().values
        yes_vals = tl_long.loc[event_mask, "log2_tl"].dropna().values
        n_pos = int(event_mask.sum())

        if n_pos < 2:
            ax.text(0.5, 0.5, f"n={n_pos}\ntoo few", ha="center", va="center",
                    transform=ax.transAxes, fontsize=11, color=CLR_GRAY)
            if row_idx == 0:
                ax.set_title(event_name, fontsize=12, fontweight="bold", color=color)
            continue

        observed_median = np.median(yes_vals)

        # Global permutation test
        null_global = np.zeros(N_PERMS)
        for p_i in range(N_PERMS):
            sampled = rng.choice(valid_idx, size=n_pos, replace=False)
            null_global[p_i] = np.median(all_vals[sampled])
        obs_abs = abs(observed_median - cohort_med)
        p_global = (np.sum(np.abs(null_global - cohort_med) >= obs_abs) + 1) / (N_PERMS + 1)

        # Fisher sign test
        no_below = int(((~event_mask) & (all_vals < cohort_med)).sum())
        no_above = int(((~event_mask) & (all_vals >= cohort_med)).sum())
        yes_below = int((event_mask & (all_vals < cohort_med)).sum())
        yes_above = int((event_mask & (all_vals >= cohort_med)).sum())
        _, p_fisher = fisher_exact_test([[yes_below, yes_above], [no_below, no_above]])

        # Violin
        parts = ax.violinplot([no_vals, yes_vals], positions=[0, 1], showextrema=False)
        for pc, c in zip(parts["bodies"], ["#B0B0B0", color]):
            pc.set_facecolor(c)
            pc.set_edgecolor(CLR_GRAY)
            pc.set_alpha(0.7)

        # Boxplot
        bp = ax.boxplot([no_vals, yes_vals], positions=[0, 1], widths=0.18,
                         patch_artist=True, showfliers=False, zorder=4)
        for patch in bp["boxes"]:
            patch.set_facecolor("white")
            patch.set_edgecolor(CLR_GRAY)
        for element in ["whiskers", "caps", "medians"]:
            for line in bp[element]:
                line.set_color(CLR_GRAY)

        # Sina
        for i, (vals, c) in enumerate([(no_vals, "#B0B0B0"), (yes_vals, color)]):
            jitter = np.random.normal(0, 0.04, len(vals))
            ax.scatter(np.full(len(vals), i) + jitter, vals,
                       s=4, color=c, alpha=0.2, zorder=3)

        # p-value annotation
        ymax = max(no_vals.max(), yes_vals.max()) * 1.03
        ax.plot([0, 0, 1, 1], [ymax * 0.96, ymax, ymax, ymax * 0.96], color=CLR_GRAY, linewidth=1)
        ax.text(0.5, ymax * 1.01,
                f"p={p_global:.3f} {star(p_global)}",
                ha="center", va="bottom", fontsize=9, color=CLR_GRAY,
                fontweight="bold" if p_global < 0.05 else "normal")

        ax.set_xticks([0, 1])
        if row_idx == 2:
            ax.set_xticklabels([f"No\n({len(no_vals)})", f"Yes\n({len(yes_vals)})"], fontsize=10)
        else:
            ax.set_xticklabels(["No", "Yes"], fontsize=10)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.axhline(cohort_med, color=CLR_GRAY, linewidth=0.5, linestyle="--", alpha=0.4)

        if row_idx == 0:
            ax.set_title(event_name, fontsize=12, fontweight="bold", color=color)

        print(f"  {event_name:<22s}  obs_med={observed_median:+.3f}  p_perm={p_global:.4f} {star(p_global):>3s}  "
              f"p_fisher={p_fisher:.4f} {star(p_fisher):>3s}")

    short_label = tl_label.split(" (")[0]
    axes_row[0].set_ylabel(f"Global z-score\n{short_label}", fontsize=12)

plt.tight_layout()
fig2.savefig(FIGDIR / "tl_by_event_all_measures.pdf", dpi=300, bbox_inches="tight")
fig2.savefig(FIGDIR / "tl_by_event_all_measures.png", dpi=300, bbox_inches="tight")
plt.close()
print("Saved tl_by_event_all_measures.pdf/.png")


# =====================================================================
#  Figure 4: Main composite figure (figure4_tl_association)
# =====================================================================
print("\n" + "=" * 70)
print("  Figure: figure4_tl_association (main composite)")
print("=" * 70)

# Define signatures
m["brk_core"] = ((m["auto_foldback_Count"] > 0) | (m["Foldback_Translocation"] > 0) |
                  (m["sv_count_INV"] > 0) | (m["Chromothripsis"] > 0)).astype(int)
m["brk_full"] = (m["brk_core"] | (m["cnv_rate_per_mb_Loss"] > 0) | (m["aneu_loss"] > 0)).astype(int)

m["amp_core"] = ((m["ecdna_multi_region_intra_chr"] > 0) |
                  (m["auto_dup_inv_segment_Count"] > 0)).astype(int)
m["amp_full"] = (m["amp_core"] | (m["HAR_ecDNA_count"] > 0) | (m["auto_complex_inv_Count"] > 0) |
                  (m["Single_Chr_Complex_Rearrangement"] > 0) |
                  (m["Multiple_Chr_Complex_Rearrangement"] > 0)).astype(int)

for label, col in [("Breakage core", "brk_core"), ("Breakage full", "brk_full"),
                    ("Amplification core", "amp_core"), ("Amplification full", "amp_full")]:
    print(f"  {label}: {m[col].sum()} arms ({m[col].mean()*100:.1f}%)")

# Sample-level for Panel A
sample_tl = telo_filt.groupby("Sample_ID")["TL_kb"].median().reset_index()
sample_tl.columns = ["Sample_ID", "sample_tl_median"]
sample_events = df.groupby("Sample_ID")["HAR_ecDNA_count"].sum().reset_index()
sample_events.columns = ["Sample_ID", "total_HAR_ecDNA"]
sample = sample_tl.merge(sample_events, on="Sample_ID")
sample["has_HAR"] = (sample["total_HAR_ecDNA"] > 0).astype(int)
tl_no_har = sample.loc[sample["has_HAR"] == 0, "sample_tl_median"]
tl_has_har = sample.loc[sample["has_HAR"] == 1, "sample_tl_median"]
stat_w, p_w = stats.mannwhitneyu(tl_has_har, tl_no_har, alternative="two-sided")

# ---- Figure layout ----
fig = plt.figure(figsize=(22, 15))
gs_main = gridspec.GridSpec(2, 12, hspace=0.45, wspace=0.8, height_ratios=[1, 1])

# Panel A: Sample-level boxplot
ax_a = fig.add_subplot(gs_main[0, 0:4])

colors_a = [CLR_NEUTRAL, CLR_AMP]
violin_box_sina(ax_a, [tl_no_har.values, tl_has_har.values], [0, 1], colors_a)
for i, (vals, color) in enumerate([(tl_no_har, CLR_NEUTRAL), (tl_has_har, CLR_AMP)]):
    jitter = np.random.normal(0, 0.05, len(vals))
    ax_a.scatter(np.full(len(vals), i) + jitter, vals.values,
                 s=30, color=color, edgecolor="white", linewidth=0.5, alpha=0.8, zorder=5)

ax_a.set_xticks([0, 1])
ax_a.set_xticklabels([f"No HAR/ecDNA\n(n={len(tl_no_har)})",
                       f"HAR/ecDNA > 0\n(n={len(tl_has_har)})"], fontsize=12)
ax_a.set_ylabel("Sample Median Telomere Length (kb)", fontsize=14)
ax_a.set_title("Tumors with Amplification Events\nHave Longer Telomeres",
               fontsize=14, fontweight="bold")

ymax_a = max(tl_no_har.max(), tl_has_har.max()) * 1.08
stars_a = "**" if p_w < 0.01 else ("*" if p_w < 0.05 else "ns")
add_bracket(ax_a, 0, 1, ymax_a, f"p = {p_w:.3f} {stars_a}")
ax_a.text(-0.12, 1.05, "A", transform=ax_a.transAxes, fontsize=18, fontweight="bold", va="top")
ax_a.spines["top"].set_visible(False)
ax_a.spines["right"].set_visible(False)
print(f"Panel A: Wilcoxon p={p_w:.4f}")

# Panel B: Forest plot — all signature components
ax_b = fig.add_subplot(gs_main[0, 5:12])

forest_events = [
    ("auto_dup_inv_segment_Count", "TL Median", "Dup inv segment", "amp", "primary"),
    ("ecdna_multi_region_intra_chr", "TL Median", "ecDNA multi-region", "amp", "primary"),
    ("auto_complex_inv_Count", "TL Median", "Complex inv", "amp", "secondary"),
    ("Single_Chr_Complex_Rearrangement", "TL Median", "Single chr complex", "amp", "secondary"),
    ("HAR_ecDNA_count", "TL Median", "HAR ecDNA", "amp", "secondary"),
    ("Multiple_Chr_Complex_Rearrangement", "TL Median", "Multi chr complex", "amp", "secondary"),
    ("Chromothripsis", "TL Median", "Chromothripsis", "brk", "primary"),
    ("auto_foldback_Count", "TL Median", "Foldback inversion", "brk", "primary"),
    ("Foldback_Translocation", "TL Median", "Foldback translocation", "brk", "primary"),
    ("sv_count_INV", "TL Median", "Reciprocal inversion", "brk", "primary"),
    ("aneu_loss", "TL Median", "Arm-level loss", "brk", "secondary"),
    ("cnv_rate_per_mb_Loss", "TL Median", "Segmental CNV loss", "brk", "secondary"),
]

tl_col_map = {"TL Min": "tl_min", "TL Max": "tl_max", "TL Median": "tl_median",
              "TL Avg": "final_average_TL_p75_kb"}

forest_data = []
for prop, tl_measure, nice_name, sig_type, tier in forest_events:
    row = dichot_df[(dichot_df["property"] == prop) & (dichot_df["TL_measure"] == tl_measure)]
    if len(row) > 0:
        r = row.iloc[0]
        forest_data.append({
            "name": nice_name, "tl_measure": tl_measure,
            "beta": r["beta_lmer"], "se": r["se_lmer"],
            "p": r["p_lmer"], "n": int(r["n_present"]),
            "sig_type": sig_type, "tier": tier,
        })
    else:
        tl_col = tl_col_map.get(tl_measure, tl_measure)
        m["_tmp_event"] = (m[prop] > 0).astype(int)
        n_present = m["_tmp_event"].sum()
        try:
            mod = smf.mixedlm(f"{tl_col} ~ _tmp_event", m.dropna(subset=[tl_col]), groups="Sample_ID")
            res = mod.fit(reml=True)
            forest_data.append({
                "name": nice_name, "tl_measure": tl_measure,
                "beta": res.fe_params["_tmp_event"], "se": res.bse["_tmp_event"],
                "p": res.pvalues["_tmp_event"], "n": int(n_present),
                "sig_type": sig_type, "tier": tier,
            })
        except Exception as e:
            print(f"  WARNING: Could not compute {nice_name}: {e}")

amp_events = sorted([d for d in forest_data if d["sig_type"] == "amp"], key=lambda x: x["beta"], reverse=True)
brk_events = sorted([d for d in forest_data if d["sig_type"] == "brk"], key=lambda x: x["beta"], reverse=True)
forest_data = amp_events + brk_events

for i, d in enumerate(forest_data):
    color = CLR_AMP if d["sig_type"] == "amp" else CLR_BRK
    mfc = color if d["tier"] == "primary" else "white"
    ci_lo = d["beta"] - 1.96 * d["se"]
    ci_hi = d["beta"] + 1.96 * d["se"]
    ax_b.plot([ci_lo, ci_hi], [i, i], color=color, linewidth=1.5, zorder=4)
    ax_b.plot(d["beta"], i, marker="o", color=color, markerfacecolor=mfc,
              markersize=7, markeredgewidth=1.5, zorder=5)

ax_b.axvline(0, color=CLR_GRAY, linewidth=0.8, linestyle="--", zorder=1)

sep_y = len(amp_events) - 0.5
ax_b.axhline(sep_y, color=CLR_GRAY, linewidth=0.5, linestyle=":", alpha=0.5)

ax_b.set_yticks(range(len(forest_data)))
ax_b.set_yticklabels([d["name"] for d in forest_data], fontsize=12)
ax_b.set_xlabel("Effect on Telomere Length (kb)\npresent − absent, mixed model", fontsize=13)
ax_b.set_title("Within-Patient TL Difference by Genomic Event\n(univariate mixed models)",
               fontsize=14, fontweight="bold")
ax_b.invert_yaxis()
ax_b.spines["top"].set_visible(False)
ax_b.spines["right"].set_visible(False)

xlim = ax_b.get_xlim()
ax_b.set_xlim(xlim[0] - 0.3, max(xlim[1], 2.0) + 3.0)
for i, d in enumerate(forest_data):
    sig_str = "***" if d["p"] < 0.001 else ("**" if d["p"] < 0.01 else ("*" if d["p"] < 0.05 else ""))
    tl_tag = d["tl_measure"].replace("TL ", "")
    tier_tag = "" if d["tier"] == "primary" else " (2°)"
    ax_b.text(ax_b.get_xlim()[1] - 0.05, i,
              f"p={d['p']:.3f}{sig_str} [{tl_tag}, n={d['n']}]{tier_tag}",
              va="center", ha="right", fontsize=10, color=CLR_GRAY)

ax_b.text(ax_b.get_xlim()[0] + 0.1, 0.0, "Amplification", fontsize=12,
          fontweight="bold", color=CLR_AMP, va="center")
ax_b.text(ax_b.get_xlim()[0] + 0.1, sep_y + 0.5, "Breakage", fontsize=12,
          fontweight="bold", color=CLR_BRK, va="center")
ax_b.text(-0.12, 1.05, "B", transform=ax_b.transAxes, fontsize=18, fontweight="bold", va="top")
print(f"Panel B: Forest plot with {len(forest_data)} events")

# Panel C: Breakage — TL Max
ax_c = fig.add_subplot(gs_main[1, 0:4])

brk_core_no = m.loc[m["brk_core"] == 0, "tl_max"].dropna()
brk_core_yes = m.loc[m["brk_core"] == 1, "tl_max"].dropna()
brk_full_no = m.loc[m["brk_full"] == 0, "tl_max"].dropna()
brk_full_yes = m.loc[m["brk_full"] == 1, "tl_max"].dropna()

mod_brk_core = smf.mixedlm("tl_max ~ brk_core", m.dropna(subset=["tl_max"]), groups="Sample_ID").fit(reml=True)
mod_brk_full = smf.mixedlm("tl_max ~ brk_full", m.dropna(subset=["tl_max"]), groups="Sample_ID").fit(reml=True)
p_bc = mod_brk_core.pvalues["brk_core"]; b_bc = mod_brk_core.fe_params["brk_core"]
p_bf = mod_brk_full.pvalues["brk_full"]; b_bf = mod_brk_full.fe_params["brk_full"]

violin_box_sina(ax_c,
    [brk_core_no.values, brk_core_yes.values, brk_full_no.values, brk_full_yes.values],
    [0, 1, 2.5, 3.5],
    [CLR_NEUTRAL, CLR_BRK, CLR_NEUTRAL, CLR_BRK])

ax_c.set_xticks([0, 1, 2.5, 3.5])
ax_c.set_xticklabels([f"No\n(n={len(brk_core_no)})", f"Yes\n(n={len(brk_core_yes)})",
                       f"No\n(n={len(brk_full_no)})", f"Yes\n(n={len(brk_full_yes)})"], fontsize=11)

ymax_b = max(brk_core_no.max(), brk_core_yes.max()) * 0.93
stars_bc = "**" if p_bc < 0.01 else ("*" if p_bc < 0.05 else "ns")
stars_bf = "**" if p_bf < 0.01 else ("*" if p_bf < 0.05 else "ns")
add_bracket(ax_c, 0, 1, ymax_b, f"p={p_bc:.3f}{stars_bc}\n({b_bc:+.2f} kb)", fontsize=10)
add_bracket(ax_c, 2.5, 3.5, ymax_b, f"p={p_bf:.3f}{stars_bf}\n({b_bf:+.2f} kb)", fontsize=10)

ax_c.axvline(1.75, color=CLR_GRAY, linewidth=0.5, linestyle=":", alpha=0.4)
ax_c.text(0.5, -0.18, "Core\n(foldback+INV+chromo)", ha="center", va="top",
          fontsize=11, color=CLR_BRK, transform=ax_c.get_xaxis_transform())
ax_c.text(3.0, -0.18, "Full\n(core+CNV loss+aneu loss)", ha="center", va="top",
          fontsize=11, color=CLR_BRK, transform=ax_c.get_xaxis_transform())

ax_c.set_ylabel("Maximum Allelic Telomere Length (kb)", fontsize=13)
ax_c.set_title("Breakage Signature\n(TL Max — longest allele)", fontsize=14, fontweight="bold", color=CLR_BRK)
ax_c.text(-0.12, 1.05, "C", transform=ax_c.transAxes, fontsize=18, fontweight="bold", va="top")
ax_c.spines["top"].set_visible(False)
ax_c.spines["right"].set_visible(False)
print(f"Panel C: Breakage core p={p_bc:.4f}, full p={p_bf:.4f}")

# Panel D: Amplification — TL Min
ax_d = fig.add_subplot(gs_main[1, 4:8])

amp_core_no = m.loc[m["amp_core"] == 0, "tl_min"].dropna()
amp_core_yes = m.loc[m["amp_core"] == 1, "tl_min"].dropna()
amp_full_no = m.loc[m["amp_full"] == 0, "tl_min"].dropna()
amp_full_yes = m.loc[m["amp_full"] == 1, "tl_min"].dropna()

mod_amp_core = smf.mixedlm("tl_min ~ amp_core", m.dropna(subset=["tl_min"]), groups="Sample_ID").fit(reml=True)
mod_amp_full = smf.mixedlm("tl_min ~ amp_full", m.dropna(subset=["tl_min"]), groups="Sample_ID").fit(reml=True)
p_ac = mod_amp_core.pvalues["amp_core"]; b_ac = mod_amp_core.fe_params["amp_core"]
p_af = mod_amp_full.pvalues["amp_full"]; b_af = mod_amp_full.fe_params["amp_full"]

violin_box_sina(ax_d,
    [amp_core_no.values, amp_core_yes.values, amp_full_no.values, amp_full_yes.values],
    [0, 1, 2.5, 3.5],
    [CLR_NEUTRAL, CLR_AMP, CLR_NEUTRAL, CLR_AMP])

ax_d.set_xticks([0, 1, 2.5, 3.5])
ax_d.set_xticklabels([f"No\n(n={len(amp_core_no)})", f"Yes\n(n={len(amp_core_yes)})",
                       f"No\n(n={len(amp_full_no)})", f"Yes\n(n={len(amp_full_yes)})"], fontsize=11)

ymax_c = max(amp_core_no.max(), amp_core_yes.max()) * 0.85
stars_ac = "**" if p_ac < 0.01 else ("*" if p_ac < 0.05 else "ns")
stars_af = "**" if p_af < 0.01 else ("*" if p_af < 0.05 else "ns")
add_bracket(ax_d, 0, 1, ymax_c, f"p={p_ac:.4f}{stars_ac}\n({b_ac:+.2f} kb)", fontsize=10)
add_bracket(ax_d, 2.5, 3.5, ymax_c, f"p={p_af:.4f}{stars_af}\n({b_af:+.2f} kb)", fontsize=10)

ax_d.axvline(1.75, color=CLR_GRAY, linewidth=0.5, linestyle=":", alpha=0.4)
ax_d.text(0.5, -0.18, "Core\n(ecDNA multi+dup inv)", ha="center", va="top",
          fontsize=11, color=CLR_AMP, transform=ax_d.get_xaxis_transform())
ax_d.text(3.0, -0.18, "Full\n(core+HAR+cpx inv+cpx rearr)", ha="center", va="top",
          fontsize=11, color=CLR_AMP, transform=ax_d.get_xaxis_transform())

ax_d.set_ylabel("Minimum Allelic Telomere Length (kb)", fontsize=13)
ax_d.set_title("Amplification Signature\n(TL Min — shortest allele)", fontsize=14, fontweight="bold", color=CLR_AMP)
ax_d.text(-0.12, 1.05, "D", transform=ax_d.transAxes, fontsize=18, fontweight="bold", va="top")
ax_d.spines["top"].set_visible(False)
ax_d.spines["right"].set_visible(False)
print(f"Panel D: Amplification core p={p_ac:.4f}, full p={p_af:.4f}")

# Panel E: TL measure specificity heatmap (4x3)
ax_e = fig.add_subplot(gs_main[1, 8:12])

sig_cols = ["brk_core", "brk_full", "amp_core", "amp_full"]
sig_labels = ["Breakage\nCore", "Breakage\nFull", "Amplification\nCore", "Amplification\nFull"]
tl_measures_e = ["tl_min", "tl_max", "tl_median"]
tl_labels_e = ["TL Min\n(shortest)", "TL Max\n(longest)", "TL Median"]

pval_matrix = np.zeros((4, 3))
beta_matrix = np.zeros((4, 3))

for i, sig in enumerate(sig_cols):
    for j, tl in enumerate(tl_measures_e):
        mod = smf.mixedlm(f"{tl} ~ {sig}", m.dropna(subset=[tl]), groups="Sample_ID")
        res = mod.fit(reml=True)
        pval_matrix[i, j] = res.pvalues[sig]
        beta_matrix[i, j] = res.fe_params[sig]

neglog = -np.log10(pval_matrix)
signed = neglog * np.sign(beta_matrix)

from matplotlib.colors import LinearSegmentedColormap, Normalize
panel_e_cmap = LinearSegmentedColormap.from_list("brk_amp", [CLR_BRK, "white", CLR_AMP])
panel_e_norm = Normalize(vmin=-3.5, vmax=3.5)

for i in range(4):
    for j in range(3):
        cell_color = panel_e_cmap(panel_e_norm(signed[i, j]))
        rect = plt.Rectangle((j - 0.5, i - 0.5), 1, 1, facecolor=cell_color,
                               edgecolor="white", linewidth=0.5)
        ax_e.add_patch(rect)

        p = pval_matrix[i, j]
        b = beta_matrix[i, j]
        sig_str = "***" if p < 0.001 else ("**" if p < 0.01 else ("*" if p < 0.05 else ""))
        txt_color = "white" if abs(signed[i, j]) > 1.5 else "black"
        ax_e.text(j, i, f"{b:+.2f} kb\np={p:.3f}{sig_str}", ha="center", va="center",
                  fontsize=11, color=txt_color, fontweight="bold" if sig_str else "normal")

ax_e.set_xlim(-0.5, 2.5)
ax_e.set_ylim(3.5, -0.5)
ax_e.set_xticks(range(3))
ax_e.set_xticklabels(tl_labels_e, fontsize=12)
ax_e.set_yticks(range(4))
ax_e.set_yticklabels(sig_labels, fontsize=12)
ax_e.plot([-0.5, 2.5], [1.5, 1.5], color="white", linewidth=3)

sm_colorbar = plt.cm.ScalarMappable(cmap=panel_e_cmap, norm=panel_e_norm)
sm_colorbar.set_array([])
cbar = plt.colorbar(sm_colorbar, ax=ax_e, shrink=0.7, label="Signed −log₁₀(p)", pad=0.08)
cbar.ax.axhline(-np.log10(0.05), color="black", lw=1, ls="--")
cbar.ax.axhline(np.log10(0.05), color="black", lw=1, ls="--")

ax_e.set_title("TL Measure Specificity\n(mixed model, corrected for patient)",
               fontsize=14, fontweight="bold")
ax_e.text(-0.15, 1.05, "E", transform=ax_e.transAxes, fontsize=18, fontweight="bold", va="top")
print("Panel E: TL measure specificity heatmap (4x3)")

plt.savefig(FIGDIR / "figure4_tl_association.pdf", dpi=300, bbox_inches="tight")
plt.savefig(FIGDIR / "figure4_tl_association.png", dpi=300, bbox_inches="tight")
plt.close()
print("Saved figure4_tl_association.pdf/.png")

# =====================================================================
#  Done
# =====================================================================
print(f"\nTotal time: {time.time() - t0:.1f}s")
print(f"All figures saved to {FIGDIR}")
print("Done!")
