"""
Module 2: BLASTp homolog search.

BLASTs the human DHFR sequence (query) against a local protein database
built from the other species fetched in Module 1 (subjects), using NCBI
BLAST+ (blastp/makeblastdb). Saves a hit table with %identity, alignment
length, e-value, bit score, and query coverage per species.

Human is excluded from the subject database since blasting it against
itself is a trivial 100% self-hit.
"""

import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import config  # noqa: E402

import pandas as pd  # noqa: E402
from Bio import SeqIO  # noqa: E402

INPUT_FASTA = config.SEQUENCES_DIR / "dhfr_multispecies.fasta"
QUERY_FASTA = config.BLAST_DIR / "human_query.fasta"
SUBJECT_FASTA = config.BLAST_DIR / "subject_sequences.fasta"
BLASTDB_PREFIX = config.BLAST_DIR / "blastdb" / "dhfr_subjects"
RAW_OUTPUT = config.BLAST_DIR / "blastp_raw_output.tsv"
HIT_TABLE_OUT = config.RESULTS_BLAST_DIR / "blast_hits.csv"

# Custom outfmt 6 columns: qcovs = query coverage per subject.
OUTFMT_COLUMNS = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "qcovs",
]
OUTFMT_SPEC = "6 " + " ".join(OUTFMT_COLUMNS)


def check_blast_tools():
    missing = [tool for tool in ("makeblastdb", "blastp") if shutil.which(tool) is None]
    if missing:
        print(f"ERROR: required BLAST+ tool(s) not found on PATH: {', '.join(missing)}")
        print("Install with: brew install blast   (macOS)   |   sudo apt-get install ncbi-blast+   (Linux)")
        sys.exit(1)


def load_sequences():
    if not INPUT_FASTA.exists():
        print(f"ERROR: {INPUT_FASTA} not found. Run scripts/01_fetch_sequences.py first.")
        sys.exit(1)
    records = list(SeqIO.parse(INPUT_FASTA, "fasta"))
    if not records:
        print(f"ERROR: {INPUT_FASTA} contains no sequences.")
        sys.exit(1)
    return records


def split_query_and_subjects(records):
    human = [r for r in records if r.id.startswith("Homo_sapiens|")]
    subjects = [r for r in records if not r.id.startswith("Homo_sapiens|")]
    if not human:
        print("ERROR: no Homo_sapiens record found in the fetched FASTA to use as BLAST query.")
        sys.exit(1)
    if not subjects:
        print("ERROR: no non-human sequences available to build a subject database.")
        sys.exit(1)
    return human[0], subjects


def run(cmd):
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"ERROR running: {' '.join(cmd)}")
        print(result.stdout)
        print(result.stderr)
        sys.exit(1)
    return result


def build_blast_db():
    BLASTDB_PREFIX.parent.mkdir(parents=True, exist_ok=True)
    run(
        [
            "makeblastdb",
            "-in", str(SUBJECT_FASTA),
            "-dbtype", "prot",
            "-out", str(BLASTDB_PREFIX),
            "-title", "DHFR cross-species subjects",
        ]
    )


def run_blastp():
    result = run(
        [
            "blastp",
            "-query", str(QUERY_FASTA),
            "-db", str(BLASTDB_PREFIX),
            "-outfmt", OUTFMT_SPEC,
            "-evalue", "10",
        ]
    )
    RAW_OUTPUT.write_text(result.stdout)
    return result.stdout


def parse_hits(raw_tsv_text):
    rows = [line.split("\t") for line in raw_tsv_text.strip().splitlines() if line.strip()]
    df = pd.DataFrame(rows, columns=OUTFMT_COLUMNS)

    numeric_cols = ["pident", "length", "mismatch", "gapopen", "qstart", "qend",
                     "sstart", "send", "evalue", "bitscore", "qcovs"]
    for col in numeric_cols:
        df[col] = pd.to_numeric(df[col])

    # sseqid is "Species_Name|Accession" (the header format written by Module 1).
    df["species"] = df["sseqid"].str.split("|").str[0].str.replace("_", " ")
    df["subject_accession"] = df["sseqid"].str.split("|").str[1]

    df = df.rename(
        columns={
            "pident": "percent_identity",
            "length": "alignment_length",
            "mismatch": "mismatches",
            "gapopen": "gap_opens",
            "evalue": "e_value",
            "bitscore": "bit_score",
            "qcovs": "query_coverage_pct",
        }
    )

    ordered_cols = [
        "species", "subject_accession", "percent_identity", "query_coverage_pct",
        "alignment_length", "mismatches", "gap_opens", "e_value", "bit_score",
    ]
    df = df[ordered_cols].sort_values("percent_identity", ascending=False).reset_index(drop=True)
    return df


def print_summary(df, n_subjects):
    print("\n" + "=" * 90)
    print("SUMMARY")
    print("=" * 90)
    if df.empty:
        print("No BLAST hits found against the subject database.")
    else:
        species_width = max([len(s) for s in df["species"]] + [len("Species")]) + 2
        header = (
            f"{'Species':<{species_width}}{'Accession':<16}{'%Identity':<11}{'Coverage%':<11}"
            f"{'AlnLen':<8}{'E-value':<12}{'BitScore'}"
        )
        print(header)
        print("-" * len(header))
        for _, row in df.iterrows():
            print(
                f"{row['species']:<{species_width}}"
                f"{row['subject_accession']:<16}"
                f"{row['percent_identity']:<11.1f}"
                f"{row['query_coverage_pct']:<11.1f}"
                f"{row['alignment_length']:<8}"
                f"{row['e_value']:<12.2e}"
                f"{row['bit_score']:.1f}"
            )
        print("-" * len(header))
    print(f"{len(df)} hit(s) from {n_subjects} subject sequence(s) searched")


def main():
    print("=" * 90)
    print("Module 2: BLASTp of human DHFR against fetched cross-species set")
    print("=" * 90)

    check_blast_tools()

    config.BLAST_DIR.mkdir(parents=True, exist_ok=True)
    config.RESULTS_BLAST_DIR.mkdir(parents=True, exist_ok=True)

    records = load_sequences()
    human, subjects = split_query_and_subjects(records)

    SeqIO.write([human], QUERY_FASTA, "fasta")
    SeqIO.write(subjects, SUBJECT_FASTA, "fasta")
    print(f"Query: {human.id} ({len(human.seq)} aa)")
    print(f"Subject database: {len(subjects)} sequences -> {SUBJECT_FASTA}")

    print("\nBuilding local BLAST protein database...")
    build_blast_db()

    print("Running blastp...")
    raw_output = run_blastp()

    df = parse_hits(raw_output)
    df.to_csv(HIT_TABLE_OUT, index=False)
    print(f"\nWrote hit table: {HIT_TABLE_OUT}")
    print(f"Wrote raw BLAST output: {RAW_OUTPUT}")

    print_summary(df, len(subjects))

    if df.empty or len(df) < len(subjects):
        n_missing = len(subjects) - len(df)
        if n_missing > 0:
            hit_species = set(df["species"]) if not df.empty else set()
            missing_species = [s.id.split("|")[0].replace("_", " ") for s in subjects
                                if s.id.split("|")[0].replace("_", " ") not in hit_species]
            print(f"\nNote: {n_missing} subject species had no BLAST hit at all: {', '.join(missing_species)}")

    sys.exit(0)


if __name__ == "__main__":
    main()
