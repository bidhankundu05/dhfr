"""
Module 10: PCA and pathway-level heatmap on the GSE11440 expression matrix.

Two complementary views of the same six-sample dataset used by Module 7/9:
  - An unsupervised PCA on the most variable probes, to see whether resistant
    vs. sensitive HT29 samples separate on their own (no gene list, no
    hypothesis) -- a sanity check on how strong the transcriptional
    difference between groups actually is.
  - A supervised, hypothesis-driven view: a clustered, per-gene z-scored
    heatmap of DHFR plus six folate-pathway genes (TYMS, MTHFR, SLC19A1,
    GGH, FPGS, ATIC), to see whether DHFR's resistant-vs-sensitive shift is
    an isolated effect or part of a coordinated pathway-level response.

Reuses scripts/07_geo_expression.py's own sample-grouping, expression-matrix,
and probe-annotation functions directly (loaded dynamically, since a numeric
filename prefix isn't a valid Python import target) -- so this is the same
matrix and grouping Module 7 (and Module 9) used, not a re-derived copy.
"""

import importlib.util
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import config  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.patches as mpatches  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import seaborn as sns  # noqa: E402
from sklearn.decomposition import PCA  # noqa: E402

OUT_DIR = config.RESULTS_EXPRESSION_DIR
PCA_PLOT_OUT = OUT_DIR / "pca_plot.png"
HEATMAP_OUT = OUT_DIR / "pathway_heatmap.png"
NOTES_OUT = OUT_DIR / "pca_heatmap_notes.txt"

N_TOP_VARIABLE_PROBES = 2000
PATHWAY_GENES = ["DHFR", "TYMS", "MTHFR", "SLC19A1", "GGH", "FPGS", "ATIC"]
GROUP_COLORS = {"resistant": "#d62728", "sensitive": "#4c72b0"}


def load_module_07():
    """Dynamically loads scripts/07_geo_expression.py so its
    parse_sample_groups/build_expression_matrix/annotate_probes/
    download_soft_file functions can be reused verbatim -- guaranteeing
    the exact same matrix and grouping as Module 7/9, not a hand-copied
    (and possibly drifting) reimplementation."""
    path = Path(__file__).resolve().parent / "07_geo_expression.py"
    spec = importlib.util.spec_from_file_location("geo_mod", path)
    geo_mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(geo_mod)
    return geo_mod


def sample_labels(sample_df):
    """Short, readable per-sample labels, e.g. 'R1', 'S2'."""
    tag = sample_df["resistance_status"].str[0].str.upper()
    return (tag + sample_df["biological_replicate"].astype("Int64").astype(str)).tolist()


# ---------------------------------------------------------------------------
# Part 1: PCA on the most variable probes
# ---------------------------------------------------------------------------

def run_pca(log2_matrix, n_top=N_TOP_VARIABLE_PROBES):
    variances = log2_matrix.var(axis=1)
    top_probes = variances.sort_values(ascending=False).head(n_top).index
    X = log2_matrix.loc[top_probes].T  # samples x probes

    pca = PCA(n_components=2)
    pcs = pca.fit_transform(X.values)
    pcs_df = pd.DataFrame(pcs, index=X.index, columns=["PC1", "PC2"])
    return pcs_df, pca.explained_variance_ratio_, len(top_probes)


def render_pca_plot(pcs_df, sample_df, explained_var, n_probes, out_path):
    labels = pd.Series(sample_labels(sample_df), index=sample_df["sample_id"]).loc[pcs_df.index]
    groups = sample_df.set_index("sample_id")["resistance_status"].loc[pcs_df.index]

    fig, ax = plt.subplots(figsize=(7, 6))
    for group, color in GROUP_COLORS.items():
        mask = groups == group
        ax.scatter(pcs_df.loc[mask, "PC1"], pcs_df.loc[mask, "PC2"],
                   s=140, color=color, edgecolor="black", linewidth=0.8,
                   label=group, zorder=3)
    for sample_id, row in pcs_df.iterrows():
        ax.annotate(labels[sample_id], (row["PC1"], row["PC2"]),
                    xytext=(7, 7), textcoords="offset points", fontsize=9)

    ax.axhline(0, color="gray", linewidth=0.6, zorder=1)
    ax.axvline(0, color="gray", linewidth=0.6, zorder=1)
    ax.set_xlabel(f"PC1 ({explained_var[0] * 100:.1f}% variance)")
    ax.set_ylabel(f"PC2 ({explained_var[1] * 100:.1f}% variance)")
    ax.set_title(f"{config.GEO_ACCESSION}: PCA on top {n_probes} most-variable probes\n"
                 "(log2 MAS5 signal, all samples)")
    ax.legend(title="Group", loc="best")
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)


def describe_pca_separation(pcs_df, sample_df):
    groups = sample_df.set_index("sample_id")["resistance_status"].loc[pcs_df.index]
    lines = []
    for pc in ("PC1", "PC2"):
        resistant_vals = pcs_df.loc[groups == "resistant", pc]
        sensitive_vals = pcs_df.loc[groups == "sensitive", pc]
        separates = resistant_vals.min() > sensitive_vals.max() or resistant_vals.max() < sensitive_vals.min()
        lines.append((pc, separates))
    return lines


# ---------------------------------------------------------------------------
# Part 2: pathway-gene clustered heatmap
# ---------------------------------------------------------------------------

def representative_probe(gene, log2_matrix, annotation):
    """Same tie-break rule used throughout this pipeline (Modules 7/9):
    among probes annotated to this exact gene symbol, prefer ones without
    extensive cross-hybridization risk ('_x_at'), then take the
    highest-mean-expressed one."""
    mask = annotation["Gene Symbol"].fillna("").apply(lambda s: gene in s.split(" /// "))
    probes = [p for p in annotation[mask].index if p in log2_matrix.index]
    if not probes:
        return None, []
    means = log2_matrix.loc[probes].mean(axis=1)
    non_cross_hyb = [p for p in probes if not p.endswith("_x_at")]
    pool = non_cross_hyb if non_cross_hyb else probes
    rep = means.loc[pool].idxmax()
    return rep, probes


def build_pathway_matrix(log2_matrix, annotation, genes):
    rep_probes = {}
    all_probes_by_gene = {}
    for gene in genes:
        rep, all_probes = representative_probe(gene, log2_matrix, annotation)
        if rep is None:
            print(f"    WARNING: no probe on this platform maps to gene symbol '{gene}' -- skipping.")
            continue
        rep_probes[gene] = rep
        all_probes_by_gene[gene] = all_probes

    gene_matrix = log2_matrix.loc[list(rep_probes.values())].copy()
    gene_matrix.index = list(rep_probes.keys())
    return gene_matrix, rep_probes, all_probes_by_gene


def render_pathway_heatmap(gene_matrix, sample_df, out_path):
    z = gene_matrix.sub(gene_matrix.mean(axis=1), axis=0).div(gene_matrix.std(axis=1), axis=0)

    labels = pd.Series(sample_labels(sample_df), index=sample_df["sample_id"]).loc[z.columns]
    groups = sample_df.set_index("sample_id")["resistance_status"].loc[z.columns]
    col_colors = groups.map(GROUP_COLORS).rename("Group")
    z_display = z.copy()
    z_display.columns = labels.values
    col_colors.index = labels.values

    g = sns.clustermap(
        z_display, cmap="vlag", center=0, col_colors=col_colors,
        figsize=(7, 6), linewidths=0.5, linecolor="white",
        cbar_kws={"label": "z-score (log2 signal, per gene)"},
        dendrogram_ratio=(0.18, 0.14),
    )
    g.ax_heatmap.set_xlabel("Sample")
    g.ax_heatmap.set_ylabel("Gene")
    handles = [mpatches.Patch(color=c, label=grp) for grp, c in GROUP_COLORS.items()]
    g.ax_heatmap.legend(handles=handles, title="Group", bbox_to_anchor=(1.25, 1.15), loc="upper left")
    g.figure.suptitle(f"{config.GEO_ACCESSION}: DHFR + folate-pathway genes\n"
                       "(z-scored log2 MAS5 signal, hierarchically clustered)", y=1.04)
    g.figure.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(g.figure)
    return z


# ---------------------------------------------------------------------------
# Part 3: notes
# ---------------------------------------------------------------------------

def write_notes(pca_sep, z, sample_df, out_path):
    lines = []

    pc1_sep = dict(pca_sep)["PC1"]
    pc2_sep = dict(pca_sep)["PC2"]
    if pc1_sep and not pc2_sep:
        pca_sentence = ("Resistant and sensitive samples separate cleanly along PC1 with no "
                         "overlap between groups, while PC2 does not track group identity.")
    elif pc1_sep and pc2_sep:
        pca_sentence = "Resistant and sensitive samples separate cleanly along both PC1 and PC2."
    elif pc2_sep and not pc1_sep:
        pca_sentence = ("Resistant and sensitive samples separate cleanly along PC2 (not PC1), "
                         "with no overlap between groups.")
    else:
        pca_sentence = "Resistant and sensitive samples do not cleanly separate along either PC1 or PC2."
    lines.append(pca_sentence)
    lines.append("")

    dhfr_z = z.loc["DHFR"]
    correlations = {}
    for gene in z.index:
        if gene == "DHFR":
            continue
        correlations[gene] = float(np.corrcoef(dhfr_z.values, z.loc[gene].values)[0, 1])

    concordant = sorted([g for g, r in correlations.items() if r >= 0.7], key=lambda g: -correlations[g])
    anticoncordant = sorted([g for g, r in correlations.items() if r <= -0.7], key=lambda g: correlations[g])
    weak = [g for g in correlations if g not in concordant and g not in anticoncordant]

    corr_str = ", ".join(f"{g}: r={correlations[g]:+.2f}" for g in
                          sorted(correlations, key=lambda g: -abs(correlations[g])))

    moved_with = concordant + anticoncordant
    if moved_with:
        direction_bits = []
        if concordant:
            direction_bits.append(f"tracking it directly ({', '.join(concordant)})")
        if anticoncordant:
            direction_bits.append(f"moving in the opposite direction ({', '.join(anticoncordant)})")
        lines.append(
            f"DHFR's resistant-vs-sensitive shift is not an isolated effect: {len(moved_with)} of the "
            f"other {len(correlations)} pathway genes shift in a coordinated way alongside it across "
            f"the same six samples -- {' and '.join(direction_bits)} (Pearson r on per-gene z-scores: "
            f"{corr_str}) -- consistent with a broader pathway-level transcriptional response to "
            "chronic methotrexate selection, not a DHFR-specific change."
        )
    else:
        lines.append(
            f"DHFR's resistant-vs-sensitive shift looks isolated: none of the other {len(correlations)} "
            f"pathway genes track it strongly in either direction (Pearson r on per-gene z-scores: "
            f"{corr_str}), which argues against a broad pathway-level co-response and for a more "
            "DHFR-specific change."
        )

    if weak:
        lines.append(
            f"By contrast, {', '.join(weak)} show only weak or inconsistent association with DHFR's "
            "pattern across these six samples, so the coordinated response is not uniform across "
            "every folate-pathway gene tested."
        )

    out_path.write_text("\n".join(lines) + "\n")
    return pca_sentence, correlations


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=" * 90)
    print("Module 10: PCA + pathway-gene heatmap on GSE11440 expression data")
    print("=" * 90)

    print("\n[10.1] Rebuilding the exact Module 7 expression matrix + grouping "
          "(scripts/07_geo_expression.py's own functions, same cached GEO SOFT file)...")
    geo_mod = load_module_07()
    soft_path = geo_mod.download_soft_file()
    import GEOparse
    gse = GEOparse.get_GEO(filepath=str(soft_path), silent=True)
    sample_df = geo_mod.parse_sample_groups(gse)
    matrix = geo_mod.build_expression_matrix(gse, sample_df)
    annotation = geo_mod.annotate_probes(gse)
    log2_matrix = np.log2(matrix.clip(lower=1))
    print(f"    {log2_matrix.shape[0]} probes x {log2_matrix.shape[1]} samples")

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print(f"\n[10.2] PCA on the top {N_TOP_VARIABLE_PROBES} most-variable probes...")
    pcs_df, explained_var, n_probes = run_pca(log2_matrix)
    render_pca_plot(pcs_df, sample_df, explained_var, n_probes, PCA_PLOT_OUT)
    pca_sep = describe_pca_separation(pcs_df, sample_df)
    print(f"    PC1 explains {explained_var[0] * 100:.1f}% variance, "
          f"PC2 explains {explained_var[1] * 100:.1f}% variance")
    for pc, sep in pca_sep:
        print(f"    {pc}: groups {'separate cleanly' if sep else 'do NOT separate'}")
    print(f"    Wrote {PCA_PLOT_OUT}")

    print(f"\n[10.3] Building pathway-gene expression matrix "
          f"({', '.join(PATHWAY_GENES)}) via platform annotation lookup...")
    gene_matrix, rep_probes, all_probes_by_gene = build_pathway_matrix(log2_matrix, annotation, PATHWAY_GENES)
    for gene, rep in rep_probes.items():
        others = [p for p in all_probes_by_gene[gene] if p != rep]
        note = f" (of {len(all_probes_by_gene[gene])} probes on this platform; others: {', '.join(others)})" if others else ""
        print(f"    {gene}: representative probe {rep}{note}")

    print("\n[10.4] Rendering clustered, per-gene z-scored heatmap...")
    z = render_pathway_heatmap(gene_matrix, sample_df, HEATMAP_OUT)
    print(f"    Wrote {HEATMAP_OUT}")

    print("\n[10.5] Writing notes...")
    pca_sentence, correlations = write_notes(pca_sep, z, sample_df, NOTES_OUT)
    print(f"    Wrote {NOTES_OUT}")
    print("\n" + NOTES_OUT.read_text())

    sys.exit(0)


if __name__ == "__main__":
    main()
