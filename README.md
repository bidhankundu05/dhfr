# DHFR Sequence-to-Structure-to-Function Pipeline

A computational biology pipeline that tests whether sequence conservation,
structural ligand-contact analysis, and molecular docking — validated
directly against published biochemical affinity data — can predict the
methotrexate-binding consequences of four well-characterized human
dihydrofolate reductase (DHFR) resistance mutations, and how that
enzyme-level mechanism compares to the transcriptional/amplification-driven
resistance captured by independent expression analysis.

It started as a 7-module bioinformatics-course pipeline (one protein, one
story, every syllabus topic covered: sequence retrieval, BLAST, MSA,
phylogenetics, structure, docking, expression) and grew into a 10-module
validation study once the original 6-species conservation analysis and
Welch's-t-test expression analysis were judged under-powered and were
each independently re-run at proper scale.

## Research question

*Do standard computational structural descriptors actually predict
experimentally measured methotrexate affinity changes — not just
rationalize them after the fact?*

## Key results

- **Docking is a qualitative, not quantitative, predictor.** Rigid-receptor
  AutoDock Vina docking of the four resistance mutants recovered the
  correct *direction* of binding disruption in only 1 of 4 cases (Pearson
  r = −0.54 vs. published affinity changes, not significant at n = 4).
  Allowing the mutated side chain(s) to relax improved this to 2 of 4
  (Spearman ρ = 0.80) — useful for *explaining* individual mutants after
  the fact, poorly calibrated as a predictor before the fact.
- **Conservation is graded, not binary.** Across an expanded 200-sequence
  ortholog panel (UniProt Swiss-Prot + TrEMBL, replacing the original
  6-species set), Leu22 sits above the protein-wide conservation average
  (76th percentile) while Phe31, Gln35, and Asn64 sit below it (20th–34th
  percentile) — and the physicochemical size of naturally tolerated
  substitutions at Phe31 is clearly smaller than the engineered F31R
  resistance mutation, but comparable in magnitude at Leu22 and Gln35.
- **DHFR transcript upregulation in resistant cells is real and gets
  stronger under proper statistics.** Reanalyzing GEO dataset GSE11440
  with limma (empirical-Bayes moderated t-test) instead of a per-probe
  Welch's t-test moved DHFR from marginally significant (adjusted
  p = 0.048) to solidly significant (adjusted p = 0.0083), with an
  unchanged effect size — the textbook signature of real signal recovered
  by variance moderation, not a method-dependent artifact.
- **That expression shift isn't DHFR-specific.** A clustered, per-gene
  heatmap of six other folate-pathway transcripts shows DHFR's shift
  co-occurring with coordinated changes in SLC19A1 (the methotrexate
  transporter, moving in the opposite direction — consistent with reduced
  drug uptake as an independent resistance mechanism) and GGH.

Full numbers, figures, and discussion are in [`manuscript/`](manuscript/).

## Structure

```
dhfr-pipeline/
├── config.py                       # constants: gene, species, PDB IDs, GEO accession, Entrez creds
├── data/                           # raw/intermediate data, one subfolder per module
│   ├── sequences/  blast/  alignment/  tree/  structures/  docking/  expression/  conservation/*
├── scripts/
│   ├── 00_check_environment.py     # verifies Python/R packages + CLI tools are installed
│   ├── 01_fetch_sequences.py       # Module 1: six-species DHFR retrieval (NCBI/UniProt) — legacy*
│   ├── 02_run_blast.py             # Module 2: BLASTp human DHFR vs. fetched set — legacy*
│   ├── 03_align_mafft.py           # Module 3: MAFFT MSA + six-species conservation — legacy*
│   ├── 04_build_tree.py            # Module 4: phylogenetic tree from the six-species alignment — legacy*
│   ├── 05_fetch_structure.py       # Module 5: RCSB PDB structure retrieval + rendering
│   ├── 06_dock_vina.py             # Module 6: AutoDock Vina docking, rigid + flexible, vs. published Ki
│   ├── 07_geo_expression.py        # Module 7: GSE11440 differential expression (Welch's t-test)
│   ├── 08_expand_conservation.py   # Module 8: 200-sequence UniProt conservation panel (primary analysis)
│   ├── 09_limma_reanalysis.py      # Module 9: limma re-analysis of Module 7's expression matrix
│   └── 10_expression_pca_heatmap.py# Module 10: PCA + folate-pathway co-expression heatmap
├── results/                        # final figures/tables, same subfolder layout as data/
├── manuscript/                     # write-up (Markdown/DOCX/PDF), references, viva slide source material
└── requirements.txt
```

\* Modules 1–4's six-species conservation and phylogenetic-tree analysis
were superseded by Module 8 for the primary conservation claims once a
larger, statistically powered comparison set was needed. They're kept in
the pipeline and reported as supplementary material (BLASTp screen, tree
topology sanity-check) rather than removed, since their output is still
real, correct, and load-bearing for the tree in particular.

## Setup

1. Create and activate a Python virtual environment, then install
   dependencies:

   ```bash
   python3 -m venv .venv
   source .venv/bin/activate
   pip install -r requirements.txt
   ```

2. Set your NCBI Entrez email (required by NCBI's usage policy) and,
   optionally, an API key for higher rate limits:

   ```bash
   export ENTREZ_EMAIL="you@example.com"
   export ENTREZ_API_KEY="your_key_here"
   ```

   Or edit the defaults directly in `config.py`.

3. Install R and Bioconductor's `limma` (needed only for Module 9):

   ```bash
   brew install r   # or your OS's R package/installer
   Rscript -e 'if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager"); BiocManager::install("limma")'
   pip install rpy2
   ```

   `rpy2`'s prebuilt wheel is linked against the official CRAN
   `R.framework` layout; if you installed R via Homebrew (different
   layout) and see an `ImportError`/fallback-to-ABI-mode warning on
   import, rebuild it against your actual R with
   `pip install --no-binary rpy2-rinterface --force-reinstall rpy2` —
   it auto-detects R via `R`/`Rscript` on `PATH`, no manual `R_HOME`
   needed either way.

4. Verify the environment (Python + R packages, CLI tools):

   ```bash
   python scripts/00_check_environment.py
   ```

   This checks Biopython, GEOparse, pandas, scikit-learn, pubchempy, and
   looks for `mafft`, `vina` (AutoDock Vina), `obabel` (Open Babel), and
   `fasttree` on `PATH`. Any tool reported "NOT FOUND" needs to be
   installed separately — see the install commands the script prints, or
   your OS package manager (Homebrew/apt/conda).

## Running the pipeline

Run scripts in order from the project root:

```bash
python scripts/01_fetch_sequences.py
python scripts/02_run_blast.py
python scripts/03_align_mafft.py
python scripts/04_build_tree.py
python scripts/05_fetch_structure.py
python scripts/06_dock_vina.py
python scripts/07_geo_expression.py
python scripts/08_expand_conservation.py
python scripts/09_limma_reanalysis.py
python scripts/10_expression_pca_heatmap.py
```

Each script reads its inputs from `data/<module>/` (or fetches them
fresh) and writes intermediate output there, with final figures/tables
also saved to `results/<module>/`. Modules 9 and 10 reuse Module 7's
cached GEO download and Module 7's own sample-grouping/matrix-building
functions directly (dynamically imported), rather than re-deriving them —
Module 7 must be run at least once first. Module 8 is independent (its
own UniProt query), and doesn't require Modules 1–4 or 7 to have run.

**Network access:** Modules 1, 5, 7, and 8 require live internet access
(NCBI Entrez/UniProt, RCSB PDB, GEO, UniProt REST respectively). Modules
2, 3, 4, 6, 9, and 10 run entirely on locally cached inputs from earlier
steps. If you're running in a sandboxed environment without outbound
network access, the network-dependent steps will need to be run on a
machine with internet access (your laptop or Google Colab) instead.

## Notes on tool substitutions

- **FastTree → VeryFastTree.** FastTree itself has no macOS Homebrew formula.
  `scripts/00_check_environment.py` and the phylogenetics module accept
  [VeryFastTree](https://github.com/citiususc/veryfasttree) as a disclosed,
  CLI-compatible substitute (same flags: `-nt`, `-gtr`, `-quiet`, `-log`,
  `-pseudo`, ...). Installed here via `brew install veryfasttree`.
- **AutoDock Vina.** Not packaged in Homebrew or as a prebuilt `pip` wheel for
  this Python version. Installed here from the official precompiled binary at
  the [ccsb-scripps/AutoDock-Vina releases page](https://github.com/ccsb-scripps/AutoDock-Vina/releases)
  (`vina_1.2.7_mac_aarch64` → `/opt/homebrew/bin/vina`).
- **limma via `rpy2`, not reimplemented.** Module 9 calls R's actual
  Bioconductor `limma` package through `rpy2` rather than reimplementing
  empirical-Bayes variance moderation in Python — the field-standard tool
  should do the field-standard statistics.

## Notes on Module 5 (structure)

SWISS-Model, Phyre2, and MODELLER are GUI/web-submission tools that aren't
practically scriptable end-to-end. As an automatable, equally valid
substitute, this pipeline pulls the real experimental structure from RCSB PDB
(human DHFR + methotrexate, e.g. PDB 1U72) rather than doing homology
modelling from scratch. For species lacking a crystal structure, an
AlphaFold DB predicted model can be used instead.

## Notes on Module 8 (expanded conservation)

A six-sequence alignment supports at most seven possible per-residue
identity values — too few to rank individual positions against a
meaningful background distribution. Module 8 instead retrieves every
reviewed (Swiss-Prot) UniProt entry for `gene:DHFR`, length-filters to
150–260 aa to exclude bifunctional DHFR–thymidylate-synthase fusion
proteins (plants, kinetoplastid parasites, apicomplexans), and — since the
reviewed-only set falls short of a useful minimum — tops up with
taxonomically diverse unreviewed (TrEMBL) entries capped at three
sequences per genus, for a final panel of up to 200 sequences. Per-residue
Shannon entropy, percent identity, and Grantham physicochemical distance
are then computed at every position, so the four resistance positions can
be ranked against the rest of the protein rather than reported in
isolation.

## Manuscript

[`manuscript/DHFR_pipeline_report.docx`](manuscript/DHFR_pipeline_report.md)
is the source of truth; `.docx` and `.pdf` copies (all figures/tables
embedded) and a Zotero-importable `.ris` reference file are generated from
it. Structure, section set, and numbered Vancouver-style citations are
formatted for submission to *Computational and Structural Biotechnology
Journal*.
