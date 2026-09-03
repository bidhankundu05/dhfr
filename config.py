"""
Central configuration for the DHFR sequence-to-structure-to-function pipeline.

All scripts under scripts/ import constants from here so that species lists,
accessions, and file paths only need to be changed in one place.
"""

import os
from pathlib import Path

# ---------------------------------------------------------------------------
# Project paths
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parent
DATA_DIR = PROJECT_ROOT / "data"
RESULTS_DIR = PROJECT_ROOT / "results"

SEQUENCES_DIR = DATA_DIR / "sequences"
BLAST_DIR = DATA_DIR / "blast"
ALIGNMENT_DIR = DATA_DIR / "alignment"
TREE_DIR = DATA_DIR / "tree"
STRUCTURES_DIR = DATA_DIR / "structures"
DOCKING_DIR = DATA_DIR / "docking"
EXPRESSION_DIR = DATA_DIR / "expression"

RESULTS_SEQUENCES_DIR = RESULTS_DIR / "sequences"
RESULTS_BLAST_DIR = RESULTS_DIR / "blast"
RESULTS_ALIGNMENT_DIR = RESULTS_DIR / "alignment"
RESULTS_TREE_DIR = RESULTS_DIR / "tree"
RESULTS_STRUCTURES_DIR = RESULTS_DIR / "structures"
RESULTS_DOCKING_DIR = RESULTS_DIR / "docking"
RESULTS_EXPRESSION_DIR = RESULTS_DIR / "expression"
RESULTS_CONSERVATION_DIR = RESULTS_DIR / "conservation"

# ---------------------------------------------------------------------------
# Gene / protein of interest
# ---------------------------------------------------------------------------
GENE_SYMBOL = "DHFR"
HUMAN_UNIPROT = "P00374"

# Species set for cross-species sequence retrieval, BLAST, MSA, and phylogenetics.
# Mix of mammals and non-mammals so the tree has real topology to discuss.
SPECIES_LIST = [
    "Homo sapiens",
    "Mus musculus",
    "Rattus norvegicus",
    "Escherichia coli",
    "Saccharomyces cerevisiae",
    "Danio rerio",
]

# "DHFR" is a mammal-centric symbol. Some species carry the same enzyme under
# a different standard gene name (per NCBI Gene / UniProt), so a literal
# "DHFR[Gene Name]" search silently finds nothing (or the wrong gene) for
# them. Listed aliases are tried, in order, after GENE_SYMBOL itself.
GENE_SYMBOL_ALIASES = {
    "Saccharomyces cerevisiae": ["DFR1"],  # SGD standard name; DHFR is an informal synonym only
    "Escherichia coli": ["folA"],  # EcoCyc/UniProt standard name for chromosomal DHFR
}

# ---------------------------------------------------------------------------
# Structure / docking targets: WT + four well-characterized human DHFR
# methotrexate-resistance point mutants with solved crystal structures.
# ---------------------------------------------------------------------------
# Central research question this pipeline is built around: can sequence
# conservation + structural ligand-contact analysis + docking (validated
# against published biochemical affinity data) explain reduced methotrexate
# binding in these mutants, vs. the amplification/overexpression mechanism
# captured separately in Module 7?
#
# Mutation identities were verified computationally (direct sequence diff of
# each PDB entry's SEQRES against WT 1U72), not taken from text descriptions.
# Numbering is mature-protein/PDB numbering (Met1 cleaved) -- the same
# HUMAN_NUMBERING_OFFSET convention used in scripts/03_align_mafft.py.
#
# Ki_MTX = inhibition constant for methotrexate binding to the E:NADPH
# ternary complex. fold_change is computed against each source paper's OWN
# wild-type control (absolute Ki varies somewhat between labs/assay
# batches, but the within-paper mutant:WT ratio is internally consistent
# and is what matters for a thermodynamic dG = RT*ln(Ki) comparison).
REFERENCE_PDB = "1U72"   # WT human DHFR bound to methotrexate (crystal reference)

RESISTANCE_MUTANT_PDBS = {
    "1DLR": {
        "mutations": [(22, "L", "F")],
        "fold_change": 740.0,
        "fold_change_is_lower_bound": True,
        "note": (
            "Exact isolated Ki for L22F not confirmed from an accessible primary "
            "source (paper predates PMC deposit; publisher blocks automated "
            "fetch). Lewis et al. 1995 reports an aggregate 740- to 28,000-fold "
            "decrease across all four Leu22 substitutions tested (Tyr/Phe/Trp/Arg) "
            "vs WT (Ki < 0.031 nM); using the 740x lower bound as a conservative "
            "point estimate for this variant specifically."
        ),
        "source": "Lewis et al. 1995, J Biol Chem 270:5057-5064 (PMID 7890613)",
    },
    "1DLS": {
        "mutations": [(22, "L", "Y")],
        "fold_change": 11.0 / 0.031,
        "fold_change_is_lower_bound": False,
        "note": "Ki_MTX = 11 nM vs WT Ki_MTX < 0.031 nM.",
        "source": (
            "Lewis et al. 1995, J Biol Chem 270:5057-5064 (PMID 7890613); value as "
            "cited in Volpato et al. 2009, J Biol Chem 284:20079-20089"
        ),
    },
    "3EIG": {
        "mutations": [(31, "F", "R"), (35, "Q", "E")],
        "fold_change": 21.0 / 0.031,
        "fold_change_is_lower_bound": False,
        "note": (
            "Ki_MTX = 21 +/- 11 nM vs WT Ki_MTX < 0.031 nM. The double mutant's "
            "effect is synergistic: F31R alone is only ~35x, Q35E alone only "
            "~1.5x -- combined, >650x. Strongest resistance case of the four."
        ),
        "source": "Volpato et al. 2009, J Biol Chem 284:20079-20089 (PMID 19478082), Table 2",
    },
    "3F8Z": {
        "mutations": [(35, "Q", "S"), (64, "N", "S")],
        "fold_change": 0.047 / 0.093,
        "fold_change_is_lower_bound": False,
        "note": (
            "Ki_MTX = 0.047 nM vs this paper's own WT control of 0.093 nM -- "
            "NOT a significant resistance mutation (paper explicitly notes the "
            "difference is not statistically significant). This mutant was "
            "originally studied for human-vs-Pneumocystis-jirovecii DHFR active "
            "site selectivity, not isolated as a clinical resistance variant. "
            "Included here as specified, and deliberately kept as a real "
            "negative case rather than dropped."
        ),
        "source": "Cody et al. 2009, Biochemistry 48:1702-1711 (PMID 19196009), Table 3",
    },
}

# Distinct mutated positions across all four mutants (mature-protein numbering).
# Three of the four (31, 35, 64) are already in the literature active-site set
# used by scripts/03_align_mafft.py; 22 (Lewis et al. 1995) is added there too.
RESISTANCE_POSITIONS = sorted({m[0] for v in RESISTANCE_MUTANT_PDBS.values() for m in v["mutations"]})

# ---------------------------------------------------------------------------
# Expression dataset
# ---------------------------------------------------------------------------
GEO_ACCESSION = "GSE11440"

# ---------------------------------------------------------------------------
# NCBI Entrez access
# ---------------------------------------------------------------------------
# NCBI's usage policy requires a real email address on every Entrez request,
# and recommends an API key for higher rate limits. Fill in your own values
# below (or set them as environment variables) before running scripts/01.
ENTREZ_EMAIL = os.environ.get("ENTREZ_EMAIL", "bidhan.kundu05@gmail.com")
ENTREZ_API_KEY = os.environ.get("ENTREZ_API_KEY", "8933dac6288fbc2114318bee84bc41626408")
