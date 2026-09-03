"""
Module 9: limma reanalysis of the Module 7 (GSE11440) differential expression.

Module 7's differential expression used a per-probe Welch's t-test -- a
reasonable baseline, but not the field-standard tool for small-n microarray
studies. This module re-analyzes the exact same expression matrix and
parental (sensitive) vs. resistant grouping with limma: lmFit() to fit a
per-probe linear model, eBayes() to moderate the per-probe variance
estimates by borrowing information across all probes, and topTable() to
extract results. This is NOT a from-scratch reimplementation of limma's
empirical Bayes moderation -- it calls the real Bioconductor limma package
in R via rpy2, exactly as the field does.

Requires R (>=4.0) with Bioconductor's limma package installed, and the
Python rpy2 package. If either is missing this script exits with a clear
explanation of what to install, rather than reimplementing the statistics
in pure Python -- reimplementing eBayes moderation would defeat the point
of using the standard tool.

Reuses scripts/07_geo_expression.py's own sample-grouping and expression-
matrix-building functions directly (loaded dynamically since a numeric
filename prefix isn't a valid Python import target), so this is provably
the same matrix and grouping Module 7 used -- not a re-derived, possibly
divergent copy.
"""

import importlib.util
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import config  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

OUT_DIR = config.RESULTS_EXPRESSION_DIR
WELCH_DE_PATH = OUT_DIR / "differential_expression.csv"
LIMMA_RESULTS_OUT = OUT_DIR / "limma_results.csv"
COMPARISON_OUT = OUT_DIR / "welch_vs_limma_comparison.csv"

SIGNIFICANCE_ALPHA = 0.05
TARGET_GENE = "DHFR"


def load_module_07():
    """Dynamically loads scripts/07_geo_expression.py so its
    parse_sample_groups/build_expression_matrix/annotate_probes/
    download_soft_file functions can be reused verbatim -- guaranteeing
    the exact same matrix and grouping as Module 7, not a hand-copied
    (and possibly drifting) reimplementation."""
    path = Path(__file__).resolve().parent / "07_geo_expression.py"
    spec = importlib.util.spec_from_file_location("geo_mod", path)
    geo_mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(geo_mod)
    return geo_mod


def require_rpy2_and_limma():
    """Imports rpy2 and loads limma, or exits with a clear, actionable
    error -- per the brief, falling back to a from-scratch Python
    reimplementation of eBayes moderation is explicitly not an option."""
    try:
        import rpy2.robjects as ro  # noqa: F401
        from rpy2.robjects.packages import importr
    except ImportError as exc:
        print(
            "ERROR: the 'rpy2' Python package is not installed, so this module cannot "
            "call R's limma directly (and will not reimplement empirical Bayes "
            "moderation from scratch in Python -- that would defeat the point of using "
            "the standard tool).\n"
            f"  Import error: {exc}\n\n"
            "To fix: install R (>=4.0) and Bioconductor's limma package, then:\n"
            "    pip install rpy2\n\n"
            "If R/Bioconductor limma cannot be installed in this environment at all, "
            "the limma analysis will need to be run separately on a machine that has "
            "R + Bioconductor limma installed, using an equivalent lmFit()/eBayes()/"
            "topTable() script."
        )
        sys.exit(1)

    try:
        limma = importr("limma")
    except Exception as exc:  # noqa: BLE001
        print(
            "ERROR: R is available via rpy2, but the Bioconductor 'limma' package is "
            "not installed in that R installation.\n"
            f"  Error: {exc}\n\n"
            "To fix, from a shell with that R on PATH:\n"
            '    Rscript -e \'if (!requireNamespace("BiocManager", quietly=TRUE)) '
            'install.packages("BiocManager"); BiocManager::install("limma")\''
        )
        sys.exit(1)

    return limma


def run_limma(log2_matrix, groups):
    """log2_matrix: probes x samples DataFrame (log2 intensity).
    groups: Series aligned to log2_matrix.columns, values "resistant"/"sensitive".
    Returns a DataFrame indexed by probe_id with limma's topTable() columns.
    """
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri

    with (ro.default_converter + pandas2ri.converter).context():
        r_expr_df = ro.conversion.get_conversion().py2rpy(log2_matrix)
        r_group_raw = ro.conversion.get_conversion().py2rpy(
            pd.Series(groups.values, name="group")
        )
    ro.r.assign("expr_df", r_expr_df)
    ro.r.assign("group_raw", r_group_raw)

    # ~0+group (no intercept) + an explicit contrast is the standard,
    # unambiguous way to set up a two-group limma comparison; the
    # "resistant - sensitive" contrast matches Module 7's own
    # (resistant mean) - (sensitive mean) sign convention exactly, so
    # log2FC signs are directly comparable between the two methods.
    ro.r("""
        expr <- as.matrix(expr_df)
        group <- factor(group_raw, levels = c("resistant", "sensitive"))
        design <- model.matrix(~0 + group)
        colnames(design) <- levels(group)
        fit <- lmFit(expr, design)
        contrast <- makeContrasts(ResistantVsSensitive = resistant - sensitive, levels = design)
        fit2 <- contrasts.fit(fit, contrast)
        fit2 <- eBayes(fit2)
        tt <- topTable(fit2, number = Inf, sort.by = "P")
    """)

    with (ro.default_converter + pandas2ri.converter).context():
        tt = ro.conversion.get_conversion().rpy2py(ro.r("tt"))
    tt.index.name = "probe_id"
    return tt.reset_index()


def select_representative_probe(dhfr_rows, log2fc_col, mean_col_a, mean_col_b):
    """Same tie-break rule as Module 7's report_dhfr_result(): prefer
    probes without extensive cross-hybridization risk (i.e. not '_x_at'
    when a '_at'/'_s_at' alternative exists), then pick the
    highest-expressed one among those."""
    mean_overall = (dhfr_rows[mean_col_a] + dhfr_rows[mean_col_b]) / 2
    non_cross_hyb = dhfr_rows[~dhfr_rows["probe_id"].str.endswith("_x_at")]
    pool = non_cross_hyb if not non_cross_hyb.empty else dhfr_rows
    return pool.loc[mean_overall.loc[pool.index].idxmax(), "probe_id"]


def main():
    print("=" * 90)
    print("Module 9: limma reanalysis of Module 7's GSE11440 differential expression")
    print("=" * 90)

    if not WELCH_DE_PATH.exists():
        print(f"ERROR: {WELCH_DE_PATH} not found. Run scripts/07_geo_expression.py first "
              "-- this module reuses its Welch's t-test results for the comparison table.")
        sys.exit(1)

    print("\n[9.1] Checking for R + Bioconductor limma (via rpy2)...")
    require_rpy2_and_limma()
    print("    OK -- limma available.")

    print("\n[9.2] Rebuilding the exact Module 7 expression matrix + grouping "
          "(scripts/07_geo_expression.py's own functions, same cached GEO SOFT file)...")
    geo_mod = load_module_07()
    soft_path = geo_mod.download_soft_file()
    import GEOparse
    gse = GEOparse.get_GEO(filepath=str(soft_path), silent=True)
    sample_df = geo_mod.parse_sample_groups(gse)
    matrix = geo_mod.build_expression_matrix(gse, sample_df)
    annotation = geo_mod.annotate_probes(gse)
    print(f"    {matrix.shape[0]} probes x {matrix.shape[1]} samples "
          f"({(sample_df.resistance_status == 'resistant').sum()} resistant, "
          f"{(sample_df.resistance_status == 'sensitive').sum()} sensitive) -- "
          "same AFFX-exclusion and grouping logic as Module 7.")

    # Same clip(lower=1) + log2 transform as Module 7's Welch's t-test, so
    # the two methods are compared on identical input values -- the only
    # thing that differs is the statistical test itself.
    log2_matrix = np.log2(matrix.clip(lower=1))
    groups = sample_df.set_index("sample_id").loc[log2_matrix.columns, "resistance_status"]

    print("\n[9.3] Running limma: lmFit() -> eBayes() -> topTable()...")
    print("    (a 'zero sample variances ... offset away from zero' warning from R below, "
          "if any, is limma's own routine handling of a small-n edge case, not an error)")
    tt = run_limma(log2_matrix, groups)
    tt = tt.merge(annotation, left_on="probe_id", right_index=True, how="left")
    tt = tt.rename(columns={"Gene Symbol": "gene", "Gene Title": "gene_title"})
    tt = tt[["probe_id", "gene", "gene_title", "logFC", "t", "P.Value", "adj.P.Val"]]
    tt = tt.sort_values("adj.P.Val").reset_index(drop=True)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    tt.to_csv(LIMMA_RESULTS_OUT, index=False)
    print(f"    {len(tt)} probes tested. Wrote {LIMMA_RESULTS_OUT}")
    n_limma_sig = (tt["adj.P.Val"] < SIGNIFICANCE_ALPHA).sum()
    print(f"    {n_limma_sig} / {len(tt)} probes significant at adj.P.Val < {SIGNIFICANCE_ALPHA}")

    print(f"\n[9.4] Comparing {TARGET_GENE} between Welch's t-test (Module 7) and limma...")
    welch_df = pd.read_csv(WELCH_DE_PATH)
    welch_dhfr = welch_df[welch_df["gene"] == TARGET_GENE].copy()
    limma_dhfr = tt[tt["gene"] == TARGET_GENE].copy()
    if welch_dhfr.empty or limma_dhfr.empty:
        print(f"ERROR: no probes annotated to gene symbol '{TARGET_GENE}' found in one or "
              "both result sets -- cannot build the comparison table.")
        sys.exit(1)

    representative = select_representative_probe(
        welch_dhfr, "log2_fold_change", "mean_control", "mean_resistant"
    )

    merged = welch_dhfr[["probe_id", "gene", "log2_fold_change", "adjusted_p_value"]].merge(
        limma_dhfr[["probe_id", "logFC", "P.Value", "adj.P.Val"]], on="probe_id", how="inner"
    )
    merged = merged.rename(columns={
        "log2_fold_change": "welch_log2fc",
        "adjusted_p_value": "welch_adj_p_value",
        "logFC": "limma_log2fc",
        "P.Value": "limma_p_value",
        "adj.P.Val": "limma_adj_p_value",
    })
    merged["representative_probe"] = merged["probe_id"] == representative
    merged["welch_significant"] = merged["welch_adj_p_value"] < SIGNIFICANCE_ALPHA
    merged["limma_significant"] = merged["limma_adj_p_value"] < SIGNIFICANCE_ALPHA
    merged = merged[["probe_id", "gene", "representative_probe", "welch_log2fc",
                      "welch_adj_p_value", "welch_significant", "limma_log2fc",
                      "limma_p_value", "limma_adj_p_value", "limma_significant"]]
    merged = merged.sort_values("representative_probe", ascending=False).reset_index(drop=True)

    merged.to_csv(COMPARISON_OUT, index=False)
    print(f"    Wrote {COMPARISON_OUT}")
    print("\n" + merged.to_string(index=False))

    rep_row = merged[merged["probe_id"] == representative].iloc[0]
    print("\n" + "=" * 90)
    print(f"{TARGET_GENE} result (representative probe {representative}):")
    print(f"  Welch's t-test : log2FC={rep_row['welch_log2fc']:.3f}, "
          f"adj.P.Val={rep_row['welch_adj_p_value']:.4g}  "
          f"({'SIGNIFICANT' if rep_row['welch_significant'] else 'not significant'} at alpha={SIGNIFICANCE_ALPHA})")
    print(f"  limma          : log2FC={rep_row['limma_log2fc']:.3f}, "
          f"p.value={rep_row['limma_p_value']:.4g}, adj.P.Val={rep_row['limma_adj_p_value']:.4g}  "
          f"({'SIGNIFICANT' if rep_row['limma_significant'] else 'not significant'} at alpha={SIGNIFICANCE_ALPHA})")
    print(f"\n  ==> DHFR {'DOES' if rep_row['limma_significant'] else 'does NOT'} come out "
          f"significant under limma (adj.P.Val < {SIGNIFICANCE_ALPHA}).")
    print("=" * 90)

    sys.exit(0)


if __name__ == "__main__":
    main()
