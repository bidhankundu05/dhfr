"""
Module 3: Multiple sequence alignment + conservation analysis.

Aligns data/sequences/dhfr_multispecies.fasta with MAFFT (subprocess), falling
back to the EBI MUSCLE REST API if MAFFT isn't installed. Then:
  - scores all pairwise sequence similarities with the BLOSUM62 matrix
  - identifies alignment columns that are 100% conserved across all species
  - cross-checks a handful of those against literature-known human DHFR
    active-site residues
  - renders a per-column conservation track for slides
"""

import shutil
import subprocess
import sys
import time
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import config  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.patches as mpatches  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import requests  # noqa: E402
from Bio import SeqIO  # noqa: E402
from Bio.Align import substitution_matrices  # noqa: E402

INPUT_FASTA = config.SEQUENCES_DIR / "dhfr_multispecies.fasta"
RAW_ALIGNMENT = config.ALIGNMENT_DIR / "dhfr_aligned.fasta"
FINAL_ALIGNMENT = config.RESULTS_ALIGNMENT_DIR / "dhfr_aligned.fasta"
SIMILARITY_MATRIX_OUT = config.RESULTS_ALIGNMENT_DIR / "blosum62_similarity_matrix.csv"
CONSERVED_NOTES_OUT = config.RESULTS_ALIGNMENT_DIR / "conserved_residues_notes.txt"
CONSERVATION_IMAGE_OUT = config.RESULTS_ALIGNMENT_DIR / "alignment_conservation.png"
RESISTANCE_SITE_CONSERVATION_OUT = config.RESULTS_ALIGNMENT_DIR / "resistance_site_conservation.csv"

EBI_MUSCLE_BASE = "https://www.ebi.ac.uk/Tools/services/rest/muscle"
EBI_POLL_INTERVAL_S = 5
EBI_POLL_TIMEOUT_S = 300

# Literature-reported human DHFR active-site residues (mature-protein / PDB
# numbering, i.e. after cleavage of the initiator Met1): Glu30, Phe31, Gln35,
# Ile60, Asn64, Arg70 -- see e.g. structural studies of active-site mutants
# (PMC3176622, ScienceDirect S0022283607010546). Our fetched UniProt sequence
# (P00374) retains Met1, so literature position N == our 1-indexed position
# N + 1; this was verified directly against the fetched sequence before use.
HUMAN_NUMBERING_OFFSET = 1
KNOWN_ACTIVE_SITE_RESIDUES = {
    22: ("L", "Leu22 -- hydrophobic pocket wall; site of characterized methotrexate-resistance mutations (Lewis et al. 1995)"),
    30: ("E", "Glu30 -- catalytic acid; protonates N5 of the pteridine ring during hydride transfer"),
    31: ("F", "Phe31 -- hydrophobic wall of the substrate/inhibitor binding pocket"),
    35: ("Q", "Gln35 -- alpha-helix 1 residue contacting the p-ABA moiety of folate/methotrexate"),
    60: ("I", "Ile60 -- substrate-binding pocket"),
    64: ("N", "Asn64 -- substrate-binding pocket"),
    70: ("R", "Arg70 -- substrate-binding pocket"),
}


# ---------------------------------------------------------------------------
# Alignment (MAFFT, falling back to EBI MUSCLE)
# ---------------------------------------------------------------------------

def run_mafft(input_fasta, output_fasta):
    with open(output_fasta, "w") as out_f:
        result = subprocess.run(
            ["mafft", "--auto", str(input_fasta)],
            stdout=out_f,
            stderr=subprocess.PIPE,
            text=True,
        )
    if result.returncode != 0:
        raise RuntimeError(f"mafft exited {result.returncode}: {result.stderr}")


def run_ebi_muscle(input_fasta, output_fasta, email):
    if not email or email == "your_email@example.com":
        raise RuntimeError(
            "config.ENTREZ_EMAIL must be a real address to use the EBI MUSCLE "
            "REST API fallback (EBI's usage policy requires a contact email, "
            "same as NCBI's)."
        )

    sequence_text = input_fasta.read_text()
    resp = requests.post(
        f"{EBI_MUSCLE_BASE}/run",
        data={"email": email, "format": "fasta", "sequence": sequence_text},
        timeout=30,
    )
    resp.raise_for_status()
    job_id = resp.text.strip()
    print(f"    Submitted EBI MUSCLE job: {job_id}")

    elapsed = 0
    status = None
    while elapsed < EBI_POLL_TIMEOUT_S:
        time.sleep(EBI_POLL_INTERVAL_S)
        elapsed += EBI_POLL_INTERVAL_S
        status_resp = requests.get(f"{EBI_MUSCLE_BASE}/status/{job_id}", timeout=30)
        status_resp.raise_for_status()
        status = status_resp.text.strip()
        print(f"    [{elapsed:>3}s] job status: {status}")
        if status == "FINISHED":
            break
        if status in ("FAILURE", "ERROR", "NOT_FOUND"):
            raise RuntimeError(f"EBI MUSCLE job {job_id} ended with status {status}")
    else:
        raise RuntimeError(f"EBI MUSCLE job {job_id} did not finish within {EBI_POLL_TIMEOUT_S}s")

    result_resp = requests.get(f"{EBI_MUSCLE_BASE}/result/{job_id}/aln-fasta", timeout=30)
    result_resp.raise_for_status()
    output_fasta.write_text(result_resp.text)


def get_alignment():
    config.ALIGNMENT_DIR.mkdir(parents=True, exist_ok=True)
    if shutil.which("mafft"):
        print("Using local MAFFT (mafft --auto)...")
        run_mafft(INPUT_FASTA, RAW_ALIGNMENT)
    else:
        print("mafft not found on PATH -- falling back to EBI MUSCLE REST API...")
        run_ebi_muscle(INPUT_FASTA, RAW_ALIGNMENT, config.ENTREZ_EMAIL)

    records = list(SeqIO.parse(RAW_ALIGNMENT, "fasta"))
    if not records:
        raise RuntimeError(f"Alignment produced no sequences ({RAW_ALIGNMENT})")

    config.RESULTS_ALIGNMENT_DIR.mkdir(parents=True, exist_ok=True)
    SeqIO.write(records, FINAL_ALIGNMENT, "fasta")
    return records


# ---------------------------------------------------------------------------
# Conservation analysis
# ---------------------------------------------------------------------------

def species_label(record):
    species, _, accession = record.id.partition("|")
    return species.replace("_", " "), accession


def per_column_stats(records):
    """Returns a list of dicts, one per alignment column: majority char,
    conservation score (fraction of N sequences matching the majority,
    non-gap character), and whether it's fully conserved."""
    n = len(records)
    aln_len = len(records[0].seq)
    columns = []
    for col in range(aln_len):
        chars = [str(r.seq)[col] for r in records]
        non_gap = [c for c in chars if c != "-"]
        if not non_gap:
            columns.append({"majority": "-", "score": 0.0, "fully_conserved": False})
            continue
        majority_char, majority_count = Counter(non_gap).most_common(1)[0]
        score = majority_count / n
        fully_conserved = score == 1.0
        columns.append({"majority": majority_char, "score": score, "fully_conserved": fully_conserved})
    return columns


def blosum62_similarity_matrix(records):
    matrix = substitution_matrices.load("BLOSUM62")
    labels = [f"{sp} ({acc})" for sp, acc in (species_label(r) for r in records)]
    n = len(records)
    scores = np.zeros((n, n))
    for i in range(n):
        seq_i = str(records[i].seq)
        for j in range(n):
            seq_j = str(records[j].seq)
            total = 0
            for a, b in zip(seq_i, seq_j):
                if a == "-" or b == "-":
                    continue
                total += matrix[a, b]
            scores[i, j] = total
    return pd.DataFrame(scores, index=labels, columns=labels)


def human_col_for_position(records, literature_position):
    """Maps a literature (mature-protein) residue number to an alignment
    column, via the aligned human sequence."""
    human = next(r for r in records if r.id.startswith("Homo_sapiens|"))
    target_ungapped_pos = literature_position + HUMAN_NUMBERING_OFFSET  # 1-indexed, full (Met1-inclusive) sequence
    ungapped_count = 0
    for col, ch in enumerate(str(human.seq)):
        if ch != "-":
            ungapped_count += 1
            if ungapped_count == target_ungapped_pos:
                return col
    return None


def write_conserved_residues_notes(records, columns):
    fully_conserved_cols = [i for i, c in enumerate(columns) if c["fully_conserved"]]

    lines = []
    lines.append("DHFR alignment -- conserved residue analysis")
    lines.append("=" * 60)
    lines.append(f"Alignment length: {len(columns)} columns, {len(records)} sequences")
    lines.append(
        f"Columns 100% conserved across all {len(records)} species "
        f"(identical residue, no gaps): {len(fully_conserved_cols)}"
    )
    lines.append("Conserved column indices (1-indexed, alignment coordinates):")
    lines.append(", ".join(str(i + 1) for i in fully_conserved_cols) if fully_conserved_cols else "(none)")
    lines.append("")
    lines.append("Cross-check against literature-known human DHFR active-site residues")
    lines.append("-" * 70)
    lines.append(
        "Positions below use mature-protein/PDB-style numbering as reported in "
        "the structural literature (Met1 cleaved); our fetched UniProt sequence "
        f"retains Met1, so alignment mapping applies a +{HUMAN_NUMBERING_OFFSET} offset "
        "(verified directly against the fetched sequence before use)."
    )
    lines.append("")

    for lit_pos, (expected_res, description) in KNOWN_ACTIVE_SITE_RESIDUES.items():
        col = human_col_for_position(records, lit_pos)
        if col is None:
            lines.append(f"  {description}: could not map (position beyond human sequence length)")
            continue
        residues_here = {}
        for r in records:
            sp, _ = species_label(r)
            residues_here[sp] = str(r.seq)[col]
        human_res = residues_here.get("Homo sapiens", "?")
        match_note = "matches expected residue" if human_res == expected_res else (
            f"WARNING: expected {expected_res}, found {human_res}"
        )
        is_conserved = columns[col]["fully_conserved"]
        residue_summary = ", ".join(f"{sp}={res}" for sp, res in residues_here.items())

        lines.append(f"  {description}")
        lines.append(f"    Alignment column {col + 1}, human residue = {human_res} ({match_note})")
        lines.append(f"    100% conserved across all species: {is_conserved}")
        lines.append(f"    Residues by species: {residue_summary}")
        lines.append("")

    # Highlight the catalytic acid finding explicitly -- this is the most
    # interesting cross-species result and worth spelling out for the viva.
    catalytic_col = human_col_for_position(records, 30)
    if catalytic_col is not None:
        residues_here = {species_label(r)[0]: str(r.seq)[catalytic_col] for r in records}
        ecoli_res = residues_here.get("Escherichia coli", "?")
        lines.append("Discussion")
        lines.append("-" * 70)
        lines.append(
            "The catalytic acid of DHFR (human Glu30, mature-protein numbering) sits "
            f"at alignment column {catalytic_col + 1}. Across the fetched species this "
            f"column reads: {', '.join(f'{sp}={res}' for sp, res in residues_here.items())}. "
        )
        if ecoli_res == "D" and residues_here.get("Homo sapiens") == "E":
            lines.append(
                "This matches published biochemistry: eukaryotic DHFRs use a glutamate "
                "(Glu30 in human) at this position, while E. coli DHFR uses an aspartate "
                "(Asp27 in its own numbering) as the functional equivalent -- both are "
                "the sole ionizable side chain in the active site that protonates N5 of "
                "the pteridine ring during hydride transfer (Asp27 role documented in "
                "PMC4280594 / PNAS 10.1073/pnas.1415940111). Because Glu->Asp is a "
                "conservative substitution (both acidic, both able to protonate the "
                "substrate) rather than an identical residue, this column does NOT pass "
                "our strict 100%-identity conservation filter above -- a good illustration "
                "that catalytic conservation is about chemistry, not always literal "
                "sequence identity."
            )
        else:
            lines.append(
                "(Species set at this column did not reproduce the classic Glu/Asp "
                "catalytic-acid pattern reported in the literature -- worth double-"
                "checking the mapped column against the source sequences directly.)"
            )
        lines.append("")

    CONSERVED_NOTES_OUT.write_text("\n".join(lines))
    return fully_conserved_cols


def write_resistance_site_conservation(records, columns):
    """The sequence-conservation leg of the hypothesis: an explicit per-site
    conservation score, across all fetched species, at each of the four
    positions actually mutated in the resistance structures (Modules 5-6).
    """
    rows = []
    for lit_pos in config.RESISTANCE_POSITIONS:
        expected_res, description = KNOWN_ACTIVE_SITE_RESIDUES[lit_pos]
        col = human_col_for_position(records, lit_pos)
        row = {
            "position_mature_numbering": lit_pos,
            "role": description,
            "alignment_column": (col + 1) if col is not None else None,
            "conservation_score": columns[col]["score"] if col is not None else None,
            "fully_conserved_across_all_species": columns[col]["fully_conserved"] if col is not None else None,
        }
        for r in records:
            sp, _ = species_label(r)
            row[sp] = str(r.seq)[col] if col is not None else None
        rows.append(row)

    df = pd.DataFrame(rows)
    df.to_csv(RESISTANCE_SITE_CONSERVATION_OUT, index=False)
    return df


# ---------------------------------------------------------------------------
# Conservation image
# ---------------------------------------------------------------------------

def render_conservation_image(records, columns):
    scores = np.array([c["score"] for c in columns])
    n_cols = len(columns)

    fig_width = max(10, n_cols * 0.045)
    fig, ax = plt.subplots(figsize=(fig_width, 3.6))

    cmap = plt.get_cmap("RdYlGn")
    ax.imshow(
        scores.reshape(1, -1),
        aspect="auto",
        cmap=cmap,
        vmin=0,
        vmax=1,
        extent=[0.5, n_cols + 0.5, 0, 1],
    )
    ax.set_yticks([])
    ax.set_ylim(-0.9, 1.0)
    ax.set_xlabel("Alignment position")
    ax.set_title(
        f"DHFR multiple sequence alignment conservation ({len(records)} species, "
        f"{n_cols} columns)",
        pad=12,
    )

    # Mark the mapped literature active-site columns, labelled below the strip
    # so the rotated text has clear room and never collides with the title.
    # Labels whose columns fall close together are staggered vertically so
    # the rotated text doesn't overlap (e.g. Glu30/Phe31 are one column apart).
    label_cols = sorted(
        (human_col_for_position(records, lit_pos), description.split(" -- ")[0])
        for lit_pos, (expected_res, description) in KNOWN_ACTIVE_SITE_RESIDUES.items()
        if human_col_for_position(records, lit_pos) is not None
    )
    prev_col, depth = None, 0
    for col, short_label in label_cols:
        depth = depth + 1 if (prev_col is not None and col - prev_col < 3) else 0
        prev_col = col
        label_y = -0.1 - depth * 0.32
        # The line runs all the way to the label so it's unambiguous which
        # label belongs to which column; the label's own white backing (drawn
        # on top, below) masks the short overlap so the line never visibly
        # crosses the letters themselves.
        ax.plot(
            [col + 1, col + 1], [label_y - 0.05, 1],
            color="black", linewidth=0.8, linestyle="--", alpha=0.6, zorder=1,
        )
        ax.text(
            col + 1, label_y, short_label, rotation=90, fontsize=7,
            ha="center", va="top", zorder=2,
            bbox=dict(facecolor="white", edgecolor="none", pad=0.5, alpha=0.95),
        )

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=0, vmax=1))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, pad=0.015, fraction=0.03)
    cbar.set_label("Fraction of species\nsharing majority residue", fontsize=8)

    fig.savefig(CONSERVATION_IMAGE_OUT, dpi=200, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=" * 90)
    print("Module 3: Multiple sequence alignment + conservation analysis")
    print("=" * 90)

    if not INPUT_FASTA.exists():
        print(f"ERROR: {INPUT_FASTA} not found. Run scripts/01_fetch_sequences.py first.")
        sys.exit(1)

    records = get_alignment()
    print(f"\nAlignment: {len(records)} sequences x {len(records[0].seq)} columns")
    print(f"Wrote alignment: {FINAL_ALIGNMENT}")

    columns = per_column_stats(records)
    n_fully_conserved = sum(1 for c in columns if c["fully_conserved"])
    print(f"Fully conserved columns (100% identical, no gaps): {n_fully_conserved} / {len(columns)}")

    print("\nComputing BLOSUM62 pairwise similarity matrix...")
    sim_df = blosum62_similarity_matrix(records)
    sim_df.to_csv(SIMILARITY_MATRIX_OUT)
    print(f"Wrote similarity matrix: {SIMILARITY_MATRIX_OUT}")
    print(sim_df.round(0).to_string())

    print("\nCross-checking conserved columns against literature active-site residues...")
    write_conserved_residues_notes(records, columns)
    print(f"Wrote notes: {CONSERVED_NOTES_OUT}")

    print("\nScoring conservation at the four resistance-mutation positions "
          f"{config.RESISTANCE_POSITIONS} (sequence-conservation leg of the hypothesis)...")
    site_df = write_resistance_site_conservation(records, columns)
    print(f"Wrote resistance-site conservation table: {RESISTANCE_SITE_CONSERVATION_OUT}")
    print(site_df[["position_mature_numbering", "role", "conservation_score",
                    "fully_conserved_across_all_species"]].to_string(index=False))

    print("\nRendering conservation image...")
    render_conservation_image(records, columns)
    print(f"Wrote image: {CONSERVATION_IMAGE_OUT}")

    sys.exit(0)


if __name__ == "__main__":
    main()
