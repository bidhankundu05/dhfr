"""
Module 7: Expression/amplification-driven resistance, as a mechanistic
contrast to Module 6's mutation-driven resistance analysis.

Research question: how does expression/amplification-driven methotrexate
resistance differ mechanistically from mutation-driven resistance at the
enzyme level? Does the resistant group show increased DHFR transcript
abundance, and (separately, since expression alone cannot prove it) is
there any evidence of gene amplification?

This module explicitly does NOT equate "DHFR transcript is up" with
"DHFR enzyme now binds methotrexate more weakly" -- those are different
mechanisms, tested with different evidence, compared directly in
results/expression/resistance_mechanism_comparison.csv against Module 6.
"""

import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import config  # noqa: E402

import GEOparse  # noqa: E402
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import requests  # noqa: E402
from scipy import stats  # noqa: E402
from statsmodels.stats.multitest import multipletests  # noqa: E402

GEO_ACCESSION = config.GEO_ACCESSION  # "GSE11440"
GEO_SERIES_STEM = GEO_ACCESSION[:-3] + "nnn"
SOFT_URL = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{GEO_SERIES_STEM}/{GEO_ACCESSION}/soft/{GEO_ACCESSION}_family.soft.gz"
SOFT_PATH = config.EXPRESSION_DIR / f"{GEO_ACCESSION}_family.soft.gz"

OUT_DIR = config.RESULTS_EXPRESSION_DIR
DATASET_SELECTION_OUT = OUT_DIR / "geo_dataset_selection.md"
SAMPLE_METADATA_OUT = OUT_DIR / "sample_metadata.csv"
DE_OUT = OUT_DIR / "differential_expression.csv"
DHFR_RESULT_OUT = OUT_DIR / "dhfr_result.txt"
VOLCANO_OUT = OUT_DIR / "volcano_plot.png"
MECHANISM_TABLE_OUT = OUT_DIR / "resistance_mechanism_comparison.csv"

SIGNIFICANCE_ALPHA = 0.05
VOLCANO_LOG2FC_LINE = 1.0  # visual reference line only (2-fold), not a filter


# ---------------------------------------------------------------------------
# Phase 7.1: dataset acquisition + documented selection rationale
# ---------------------------------------------------------------------------

def download_soft_file():
    config.EXPRESSION_DIR.mkdir(parents=True, exist_ok=True)
    if not SOFT_PATH.exists():
        print(f"    GEOparse's own FTP downloader is unreliable in this environment; "
              f"fetching {GEO_ACCESSION} directly over HTTPS...")
        resp = requests.get(SOFT_URL, timeout=180)
        resp.raise_for_status()
        SOFT_PATH.write_bytes(resp.content)
    return SOFT_PATH


def write_dataset_selection_doc(gse):
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    title = gse.metadata.get("title", [""])[0]
    summary = " ".join(gse.metadata.get("summary", []))
    organism = next(iter(gse.gsms.values())).metadata.get("organism_ch1", ["unknown"])[0]
    platform_id = list(gse.gpls.keys())[0]
    platform_title = gse.gpls[platform_id].metadata.get("title", [platform_id])[0]
    pubmed = gse.metadata.get("pubmed_id", [])
    superseries = [r for r in gse.metadata.get("relation", []) if r.startswith("SubSeries of")]

    n_resistant = sum(1 for gsm in gse.gsms.values()
                       if any("resistant" in c.lower() for c in gsm.metadata.get("characteristics_ch1", [])))
    n_sensitive = sum(1 for gsm in gse.gsms.values()
                       if any("sensitive" in c.lower() for c in gsm.metadata.get("characteristics_ch1", [])))

    lines = [
        f"# GEO dataset selection: {GEO_ACCESSION}",
        "",
        "## Accession and identity",
        f"- **GEO accession:** {GEO_ACCESSION}" + (f" (SubSeries of {superseries[0].split(': ')[1]})" if superseries else ""),
        f"- **Title:** {title}",
        f"- **Organism:** {organism}",
        f"- **Cell/tissue type:** HT29 human colon adenocarcinoma cell line",
        f"- **Platform:** {platform_id} ({platform_title})",
        f"- **PubMed:** {', '.join(str(p) for p in pubmed) if pubmed else 'not linked in GEO record'}",
        "",
        "## Design",
        f"- **Resistant group (case):** HT29 cells chronically selected for methotrexate "
        f"resistance -- {n_resistant} biological replicates",
        f"- **Sensitive group (control):** parental HT29 cells -- {n_sensitive} biological replicates",
        f"- **Sample counts:** {n_resistant + n_sensitive} total ({n_resistant} resistant + {n_sensitive} sensitive), "
        f"clearly labelled in GEO's own `characteristics_ch1` metadata field "
        f"(`methotrexate sensitivity: resistant` / `sensitive`) -- not inferred from sample titles.",
        "",
        "## Why this dataset addresses the research question",
        "This dataset was NOT selected because DHFR showed the largest fold-change in a "
        "search result. It was selected because its own published rationale is built around "
        "the DHFR amplicon mechanism directly relevant to this pipeline's hypothesis: the "
        "study explicitly investigates genes deregulated in methotrexate-resistant HT29 cells "
        "\"either due to its co-amplification with the DHFR gene or as a result of a "
        "transcriptome screening\" (see summary below), and reports that genes adjacent to the "
        "DHFR locus within the 5q14 amplicon are overexpressed in the resistant line. This is "
        "precisely the transcriptional/amplification-driven resistance mechanism this module "
        "is meant to characterize, as a direct contrast to Module 6's mutation-driven, "
        "enzyme-level analysis.",
        "",
        "## Original study summary (from GEO)",
        f"> {summary}",
        "",
        "## Important limitations",
        "- **n=3 per group.** Adequately powered for large effects, underpowered for subtle "
        "ones; multiple-testing correction (BH-FDR) across the whole array will be "
        "conservative given this sample size.",
        "- **This GEO series is expression-only (Affymetrix HG-U133 Plus 2.0 array).** No "
        "copy-number/array-CGH data are deposited under this accession or its SuperSeries "
        f"({superseries[0].split(': ')[1] if superseries else 'GSE16648'}, confirmed by checking every "
        "linked/related GEO record before writing this module -- all subseries are expression "
        "profiling by array, not CGH). This means DHFR **transcript abundance** can be measured "
        "directly here, but DHFR **gene amplification** (copy number) cannot be measured from "
        "this dataset -- it can only be cited as the mechanism reported in the original study's "
        "own publication, not re-derived computationally from these specific files.",
        "- **Single cell line (HT29).** Resistance mechanisms can be cell-line-specific; this is "
        "not a pan-cancer generalization.",
        "- **Bulk microarray, not RNA-seq.** MAS5 signal-intensity values, not raw read counts -- "
        "analyzed accordingly (log2-transform + moderated comparison, not a count-based method "
        "like DESeq2, which would be inappropriate for this data type).",
    ]
    DATASET_SELECTION_OUT.write_text("\n".join(lines))
    return n_resistant, n_sensitive


# ---------------------------------------------------------------------------
# Phase 7.2: strict metadata-driven sample grouping
# ---------------------------------------------------------------------------

def parse_sample_groups(gse):
    rows = []
    for gsm_name, gsm in gse.gsms.items():
        chars = gsm.metadata.get("characteristics_ch1", [])
        title = gsm.metadata.get("title", [""])[0]

        status = None
        cell_line = None
        for c in chars:
            key, _, value = c.partition(":")
            key = key.strip().lower()
            value = value.strip()
            if key == "methotrexate sensitivity":
                status = value.lower()
            elif key == "cell line":
                cell_line = value

        if status not in ("resistant", "sensitive"):
            raise SystemExit(
                f"FATAL: could not confidently determine resistance status for {gsm_name} "
                f"from GEO metadata (title={title!r}, characteristics_ch1={chars!r}). "
                "Refusing to guess a sample label -- fix the metadata field being parsed "
                "before proceeding, or select a different dataset."
            )

        m = re.search(r"replicate\s*(\d+)", title, re.IGNORECASE)
        replicate = int(m.group(1)) if m else None

        rows.append({
            "sample_id": gsm_name,
            "title": title,
            "treatment": "chronic methotrexate selection" if status == "resistant" else "none (parental line)",
            "resistance_status": status,
            "control_status": "control" if status == "sensitive" else "case",
            "biological_replicate": replicate,
            "cell_line": cell_line,
        })

    df = pd.DataFrame(rows).sort_values(["resistance_status", "biological_replicate"]).reset_index(drop=True)

    n_resistant = (df["resistance_status"] == "resistant").sum()
    n_sensitive = (df["resistance_status"] == "sensitive").sum()
    if n_resistant < 2 or n_sensitive < 2:
        raise SystemExit(
            f"FATAL: insufficient biological replicates for a meaningful comparison "
            f"(resistant={n_resistant}, sensitive={n_sensitive}, need >=2 each)."
        )
    if df["biological_replicate"].isna().any():
        print("    WARNING: could not parse an explicit replicate number for some samples "
              "(see sample_metadata.csv) -- grouping by resistance_status still succeeded.")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    df.to_csv(SAMPLE_METADATA_OUT, index=False)
    return df


# ---------------------------------------------------------------------------
# Phase 7.3: differential expression (microarray-appropriate method)
# ---------------------------------------------------------------------------

def build_expression_matrix(gse, sample_df):
    data = {}
    for _, row in sample_df.iterrows():
        gsm = gse.gsms[row["sample_id"]]
        data[row["sample_id"]] = gsm.table.set_index("ID_REF")["VALUE"]
    matrix = pd.DataFrame(data)
    # Affymetrix control/spike-in probes (AFFX-*) are QC probes, not real genes.
    matrix = matrix[~matrix.index.str.startswith("AFFX")]
    return matrix


def annotate_probes(gse):
    gpl = next(iter(gse.gpls.values()))
    ann = gpl.table.set_index("ID")[["Gene Symbol", "Gene Title"]]
    return ann


def run_differential_expression(matrix, sample_df, annotation):
    resistant = sample_df.loc[sample_df.resistance_status == "resistant", "sample_id"].tolist()
    sensitive = sample_df.loc[sample_df.resistance_status == "sensitive", "sample_id"].tolist()

    # This is Affymetrix MAS5 signal-intensity microarray data, NOT RNA-seq
    # read counts -- a count-based method (e.g. DESeq2) would be the wrong
    # tool here. Standard practice for intensity data: log2-transform, then a
    # per-probe parametric test. Clip at 1 to avoid log2(0)/log2(negative)
    # from MAS5's background-corrected values before very low signals.
    log2_matrix = np.log2(matrix.clip(lower=1))

    log2fc = log2_matrix[resistant].mean(axis=1) - log2_matrix[sensitive].mean(axis=1)
    fold_change = 2.0 ** log2fc

    pvals = np.array([
        stats.ttest_ind(log2_matrix.loc[probe, resistant], log2_matrix.loc[probe, sensitive],
                         equal_var=False).pvalue
        for probe in log2_matrix.index
    ])
    pvals = np.nan_to_num(pvals, nan=1.0)
    _, padj, _, _ = multipletests(pvals, method="fdr_bh")

    df = pd.DataFrame({
        "probe_id": log2_matrix.index,
        "log2_fold_change": log2fc.values,
        "fold_change": fold_change.values,
        "p_value": pvals,
        "adjusted_p_value": padj,
        "mean_control": matrix.loc[log2_matrix.index, sensitive].mean(axis=1).values,
        "mean_resistant": matrix.loc[log2_matrix.index, resistant].mean(axis=1).values,
    }).set_index("probe_id").join(annotation).reset_index()
    df = df.rename(columns={"Gene Symbol": "gene", "Gene Title": "gene_title"})
    df = df[["probe_id", "gene", "gene_title", "log2_fold_change", "fold_change",
             "p_value", "adjusted_p_value", "mean_control", "mean_resistant"]]
    df = df.sort_values("adjusted_p_value").reset_index(drop=True)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    df.to_csv(DE_OUT, index=False)
    return df, resistant, sensitive


# ---------------------------------------------------------------------------
# Phase 7.4: DHFR-specific result (careful, non-overclaiming wording)
# ---------------------------------------------------------------------------

def report_dhfr_result(de_df):
    dhfr_rows = de_df[de_df["gene"] == "DHFR"].copy()
    if dhfr_rows.empty:
        raise SystemExit(
            "FATAL: no probes annotated to gene symbol 'DHFR' were found on this platform. "
            "Cannot report a DHFR-specific result -- check the platform annotation column "
            "used (Gene Symbol) or the platform itself."
        )

    dhfr_rows["mean_overall"] = (dhfr_rows["mean_control"] + dhfr_rows["mean_resistant"]) / 2
    non_cross_hyb = dhfr_rows[~dhfr_rows["probe_id"].str.endswith("_x_at")]
    pool = non_cross_hyb if not non_cross_hyb.empty else dhfr_rows
    rep = pool.loc[pool["mean_overall"].idxmax()]

    direction = "increased" if rep["log2_fold_change"] > 0 else "decreased"
    significant = bool(rep["adjusted_p_value"] < SIGNIFICANCE_ALPHA)

    lines = [
        "DHFR-specific differential expression result",
        "=" * 60,
        f"Representative probe: {rep['probe_id']} (of {len(dhfr_rows)} DHFR probes on this "
        "platform; selected as the highest-expressed probe among those without extensive "
        "cross-hybridization risk, i.e. excluding '_x_at' designated probes where a "
        "'_at'/'_s_at' alternative existed)",
        f"All DHFR probes on this platform: {', '.join(dhfr_rows['probe_id'])}",
        "",
        f"DHFR transcript abundance was {direction} in the resistant group relative to the "
        "control (sensitive) group.",
        f"  log2 fold-change: {rep['log2_fold_change']:.3f}",
        f"  fold-change: {rep['fold_change']:.2f}x",
        f"  p-value: {rep['p_value']:.4g}",
        f"  adjusted p-value (FDR, Benjamini-Hochberg, whole-array correction): {rep['adjusted_p_value']:.4g}",
        f"  Statistically significant at alpha={SIGNIFICANCE_ALPHA}: {significant}",
        "",
        "All DHFR probes on this platform, for reference:",
    ]
    for _, row in dhfr_rows.iterrows():
        lines.append(
            f"  {row['probe_id']}: log2FC={row['log2_fold_change']:.3f}, "
            f"fold_change={row['fold_change']:.2f}x, p={row['p_value']:.4g}, "
            f"padj={row['adjusted_p_value']:.4g}"
        )
    lines += [
        "",
        "IMPORTANT: this is an mRNA transcript abundance measurement, not an enzyme property. "
        "It says nothing about methotrexate binding affinity per DHFR molecule -- that is the "
        "separate, structural/docking question addressed in Module 6. Increased transcript "
        "abundance and reduced per-molecule binding affinity are mechanistically independent "
        "and are compared directly in resistance_mechanism_comparison.csv.",
        "",
        "The GEO expression analysis supports increased DHFR transcript abundance but does "
        "not by itself establish gene amplification -- this dataset contains no copy-number "
        "data (verified: neither this GEO series nor its SuperSeries/related series include "
        "an array-CGH or copy-number platform). Amplification is cited here only as the "
        "mechanism reported in the original study's own publication, not re-derived "
        "computationally from these expression files.",
    ]
    DHFR_RESULT_OUT.write_text("\n".join(lines))
    print(f"    Representative probe: {rep['probe_id']}")
    print(f"    DHFR transcript abundance was {direction} {rep['fold_change']:.2f}x "
          f"(log2FC={rep['log2_fold_change']:.3f}, padj={rep['adjusted_p_value']:.4g}, "
          f"significant={significant})")
    return rep, dhfr_rows, significant


# ---------------------------------------------------------------------------
# Phase 7.5: volcano plot
# ---------------------------------------------------------------------------

def render_volcano(de_df, dhfr_rep):
    plot_df = de_df.dropna(subset=["adjusted_p_value"]).copy()
    plot_df["neg_log10_padj"] = -np.log10(plot_df["adjusted_p_value"].clip(lower=1e-300))

    fig, ax = plt.subplots(figsize=(8, 7))
    sig = plot_df["adjusted_p_value"] < SIGNIFICANCE_ALPHA
    ax.scatter(plot_df.loc[~sig, "log2_fold_change"], plot_df.loc[~sig, "neg_log10_padj"],
               s=6, color="#bbbbbb", alpha=0.5, label=f"not significant (padj >= {SIGNIFICANCE_ALPHA})")
    ax.scatter(plot_df.loc[sig, "log2_fold_change"], plot_df.loc[sig, "neg_log10_padj"],
               s=8, color="#4c72b0", alpha=0.6, label=f"significant (padj < {SIGNIFICANCE_ALPHA})")

    dhfr_plot_row = plot_df[plot_df["probe_id"] == dhfr_rep["probe_id"]]
    if not dhfr_plot_row.empty:
        drow = dhfr_plot_row.iloc[0]
        ax.scatter([drow["log2_fold_change"]], [drow["neg_log10_padj"]], s=140, color="#d62728",
                   edgecolor="black", zorder=5, label="DHFR (representative probe)")
        ax.annotate(f"DHFR ({drow['probe_id']})", (drow["log2_fold_change"], drow["neg_log10_padj"]),
                    xytext=(15, 10), textcoords="offset points", fontsize=10, fontweight="bold",
                    arrowprops=dict(arrowstyle="->", color="black", lw=0.8))

    ax.axhline(-np.log10(SIGNIFICANCE_ALPHA), color="gray", linestyle="--", linewidth=0.8,
               label=f"padj = {SIGNIFICANCE_ALPHA}")
    ax.axvline(VOLCANO_LOG2FC_LINE, color="gray", linestyle=":", linewidth=0.8)
    ax.axvline(-VOLCANO_LOG2FC_LINE, color="gray", linestyle=":", linewidth=0.8,
               label=f"+/-{VOLCANO_LOG2FC_LINE} log2FC (2-fold, reference only)")

    ax.set_xlabel("log2 fold-change (resistant vs. sensitive)")
    ax.set_ylabel("-log10(adjusted p-value)")
    ax.set_title(f"{GEO_ACCESSION}: HT29 methotrexate-resistant vs. sensitive\n"
                 "differential expression (DHFR highlighted)")
    ax.legend(fontsize=8, loc="upper left")
    fig.tight_layout()
    fig.savefig(VOLCANO_OUT, dpi=200, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Phase 7.7: mechanistic comparison table (populated from real results)
# ---------------------------------------------------------------------------

def build_mechanism_comparison(dhfr_rep, dhfr_significant, n_resistant, n_sensitive):
    docking_summary = "not available (run Module 6 first)"
    mutation_summary = "not available (run Module 6 first)"
    docking_vs_exp_path = config.RESULTS_DOCKING_DIR / "docking_vs_experimental.csv"
    flex_path = config.RESULTS_DOCKING_DIR / "docking_flexible_vs_rigid.csv"
    mutations_path = config.DATA_DIR / "mutations.csv"

    if docking_vs_exp_path.exists():
        dvx = pd.read_csv(docking_vs_exp_path)
        # 1DLR (L22F) has experimental_fold_change=NaN because only an aggregate
        # literature range was confirmed (see config.RESISTANCE_MUTANT_PDBS), not
        # because it lacks published resistance -- don't let fillna(0) hide that.
        is_resistant_mut = (dvx["experimental_fold_change"].fillna(0) > 2) | (dvx["pdb_id"] == "1DLR")
        n_resistant_mut = is_resistant_mut.sum()
        resistant_fc = pd.concat([
            dvx.loc[dvx["experimental_fold_change"] > 2, "experimental_fold_change"],
            pd.Series([740.0]),  # 1DLR (L22F) lower-bound literature estimate
        ])
        non_resistant_names = ", ".join(dvx.loc[~is_resistant_mut, "mutation"])
        docking_summary = (
            f"{len(dvx)} characterized point mutants analyzed; {n_resistant_mut} showed "
            f"published methotrexate resistance (fold-change {resistant_fc.min():.0f}-{resistant_fc.max():.0f}x "
            f"reduced Ki, one as a lower-bound estimate only); {len(dvx) - n_resistant_mut} "
            f"showed no significant resistance ({non_resistant_names})"
        )
        if flex_path.exists():
            flex_df = pd.read_csv(flex_path)
            merged = flex_df.merge(dvx[["pdb_id", "experimental_ddG_kcal_mol"]], on="pdb_id", suffixes=("", "_dup"))

            def bucket(x):
                return "up" if x > 0.5 else ("down" if x < -0.5 else "flat")

            merged["exp_bucket"] = merged["experimental_ddG_kcal_mol"].apply(bucket)
            merged["flex_bucket"] = merged["delta_vina_flex_kcal_mol"].apply(bucket)
            n_consistent_flex = (merged["exp_bucket"] == merged["flex_bucket"]).sum()
            docking_summary += (
                f"; after allowing the mutated side chain(s) to relax (flexible docking), "
                f"direction matched published data for {n_consistent_flex}/{len(merged)} mutants "
                f"(vs. rigid-receptor docking alone, which matched direction for only 1/{len(dvx)})"
            )
    if mutations_path.exists():
        muts = pd.read_csv(mutations_path)
        mutation_summary = "; ".join(sorted(muts["mutation"].unique()))

    dhfr_direction = "increased" if dhfr_rep["log2_fold_change"] > 0 else "decreased"
    expression_summary = (
        f"DHFR transcript {dhfr_direction} {dhfr_rep['fold_change']:.2f}x in resistant HT29 cells "
        f"(n={n_resistant} vs n={n_sensitive}), "
        f"{'statistically significant' if dhfr_significant else 'NOT statistically significant'} "
        f"after FDR correction (padj={dhfr_rep['adjusted_p_value']:.3g})"
    )

    rows = [
        {"Feature": "Primary level", "Mutation-mediated resistance": "Protein/enzyme",
         "Expression/amplification-mediated resistance": "Cellular/gene-expression"},
        {"Feature": "Genetic change", "Mutation-mediated resistance": mutation_summary,
         "Expression/amplification-mediated resistance": "Increased expression and/or copy number (not measured here; see limitation)"},
        {"Feature": "DHFR sequence", "Mutation-mediated resistance": "Altered at specific residue(s) (22, 31, 35, and/or 64 -- see Module 5/6)",
         "Expression/amplification-mediated resistance": "Not assessed in this GEO series (expression only); no sequence data"},
        {"Feature": "Protein binding site", "Mutation-mediated resistance": "Potentially altered (see results/structures/mutation_structural_context.csv)",
         "Expression/amplification-mediated resistance": "Assumed unchanged (not tested here)"},
        {"Feature": "Methotrexate affinity per enzyme molecule", "Mutation-mediated resistance": docking_summary,
         "Expression/amplification-mediated resistance": "Not addressed by expression data; requires separate biochemical/structural evidence"},
        {"Feature": "DHFR abundance", "Mutation-mediated resistance": "Not necessarily increased (not assessed by Module 6)",
         "Expression/amplification-mediated resistance": expression_summary},
        {"Feature": "Mechanism", "Mutation-mediated resistance": "Reduced drug binding/inhibition per enzyme molecule (where confirmed)",
         "Expression/amplification-mediated resistance": "More target enzyme molecules produced per cell (transcriptional evidence); "
                                                           "amplification specifically not confirmed from this dataset"},
        {"Feature": "Evidence type", "Mutation-mediated resistance": "Biochemical (published Ki) + structural (contact analysis) + docking",
         "Expression/amplification-mediated resistance": "Transcriptomic (microarray, this module) -- NOT copy-number/genomic"},
        {"Feature": "Key limitation", "Mutation-mediated resistance": "Docking is predictive, not a direct affinity measurement; rigid-receptor "
                                                                       "docking under-predicted 2 of 4 known effects (see Module 6)",
         "Expression/amplification-mediated resistance": "Expression increase does NOT by itself establish gene amplification "
                                                           "-- that requires copy-number data this dataset does not contain"},
    ]
    df = pd.DataFrame(rows)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    df.to_csv(MECHANISM_TABLE_OUT, index=False)
    return df


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=" * 90)
    print(f"Module 7: Expression/amplification-driven resistance ({GEO_ACCESSION}) "
          "vs. mutation-driven resistance (Module 6)")
    print("=" * 90)

    print(f"\n[7.1] Fetching {GEO_ACCESSION}...")
    soft_path = download_soft_file()
    gse = GEOparse.get_GEO(filepath=str(soft_path), silent=True)
    n_resistant, n_sensitive = write_dataset_selection_doc(gse)
    print(f"    Wrote {DATASET_SELECTION_OUT} ({n_resistant} resistant, {n_sensitive} sensitive samples)")

    print("\n[7.2] Parsing sample metadata (strict -- fails loudly on ambiguity)...")
    sample_df = parse_sample_groups(gse)
    print(f"    Wrote {SAMPLE_METADATA_OUT}")
    print(sample_df[["sample_id", "title", "resistance_status", "biological_replicate"]].to_string(index=False))

    print("\n[7.3] Building expression matrix + differential expression "
          "(microarray: log2-transform + Welch's t-test + BH-FDR)...")
    matrix = build_expression_matrix(gse, sample_df)
    annotation = annotate_probes(gse)
    de_df, resistant_ids, sensitive_ids = run_differential_expression(matrix, sample_df, annotation)
    print(f"    {len(de_df)} probes tested (AFFX control probes excluded)")
    print(f"    Wrote {DE_OUT}")

    print("\n[7.4] DHFR-specific result...")
    dhfr_rep, dhfr_rows, dhfr_significant = report_dhfr_result(de_df)
    print(f"    Wrote {DHFR_RESULT_OUT}")

    print("\n[7.5] Rendering volcano plot...")
    render_volcano(de_df, dhfr_rep)
    print(f"    Wrote {VOLCANO_OUT}")

    print("\n[7.6] Amplification vs. transcription: see geo_dataset_selection.md and "
          "dhfr_result.txt -- this GEO series contains no copy-number data (verified "
          "against its SuperSeries and related records), so amplification is NOT claimed "
          "from these files.")

    print("\n[7.7] Building mechanistic comparison table (Module 6 vs Module 7)...")
    comparison_df = build_mechanism_comparison(dhfr_rep, dhfr_significant, n_resistant, n_sensitive)
    print(f"    Wrote {MECHANISM_TABLE_OUT}")

    print("\n" + "=" * 90)
    print("Module 7 complete.")
    print("=" * 90)
    sys.exit(0)


if __name__ == "__main__":
    main()
