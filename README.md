# DHFR Sequence-to-Structure-to-Function Pipeline

A 7-module bioinformatics pipeline that traces Dihydrofolate Reductase (DHFR)
from sequence retrieval through BLAST, multiple sequence alignment,
phylogenetics, 3D structure, molecular docking, and gene expression analysis.
Built for a bioinformatics lab viva — one protein, one story, every syllabus
topic covered.

## Project story

*How does DHFR structure and expression relate to drug resistance across
species and in cancer?*

## Structure

```
dhfr-pipeline/
├── config.py                  # all constants: gene, species, PDB IDs, GEO accession, Entrez creds
├── data/                      # raw/intermediate data, one subfolder per module
│   ├── sequences/  blast/  alignment/  tree/  structures/  docking/  expression/
├── scripts/
│   ├── 00_check_environment.py  # verifies Python packages + CLI tools are installed
│   ├── 01_fetch_sequences.py    # Module 1: fetch DHFR FASTA across species (NCBI/UniProt)
│   ├── 02_run_blast.py          # Module 2: BLASTp human DHFR vs fetched set / nr
│   ├── 03_align_mafft.py        # Module 3: MAFFT multiple sequence alignment
│   ├── 04_build_tree.py         # Module 4: phylogenetic tree from alignment
│   ├── 05_fetch_structure.py    # Module 5: RCSB PDB structure retrieval + rendering
│   ├── 06_dock_vina.py          # Module 6: AutoDock Vina docking of methotrexate + analogs
│   └── 07_geo_expression.py     # Module 7: GEO differential expression + volcano plot
├── results/                   # final figures/tables, same subfolder layout as data/
├── manuscript/                # writeup / viva slide source material
└── requirements.txt
```

## Setup

1. Create and activate a virtual environment, then install dependencies:

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

3. Verify the environment (Python packages + CLI tools):

   ```bash
   python scripts/00_check_environment.py
   ```

   This checks Biopython, GEOparse, pandas, pubchempy, and looks for `mafft`,
   `vina` (AutoDock Vina), `obabel` (Open Babel), and `fasttree` on PATH.
   Any tool reported "NOT FOUND" needs to be installed separately — see the
   install commands the script prints, or your OS package manager
   (Homebrew/apt/conda).

## Running the pipeline

Run scripts in order from the project root, e.g.:

```bash
python scripts/01_fetch_sequences.py
python scripts/02_run_blast.py
python scripts/03_align_mafft.py
python scripts/04_build_tree.py
python scripts/05_fetch_structure.py
python scripts/06_dock_vina.py
python scripts/07_geo_expression.py
```

Each script reads its inputs from `data/<module>/` (or fetches them fresh)
and writes intermediate output to `data/<module>/`, with final figures/tables
also saved to `results/<module>/`.

**Network access:** Modules 1, 2, 5, and 7 require live internet access
(NCBI Entrez, RCSB PDB/AlphaFold, GEO). If you're running in a sandboxed
environment without outbound network access, those steps will need to be run
on a machine with internet access (your laptop or Google Colab) instead.

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

## Notes on Module 5 (structure)

SWISS-Model, Phyre2, and MODELLER are GUI/web-submission tools that aren't
practically scriptable end-to-end. As an automatable, equally valid
substitute, this pipeline pulls the real experimental structure from RCSB PDB
(human DHFR + methotrexate, e.g. PDB 1U72) rather than doing homology
modelling from scratch. For species lacking a crystal structure, an
AlphaFold DB predicted model can be used instead.
