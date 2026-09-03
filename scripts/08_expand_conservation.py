"""
Module 8: Expanded, properly-powered conservation analysis.

The original 6-species alignment (scripts/03_align_mafft.py) is not powered
to say anything statistically meaningful about conservation at individual
residues -- with n=6 there are only 7 possible "fraction identical" values.
This module replaces that analysis by pulling every reviewed (Swiss-Prot)
DHFR sequence from UniProt, topping up with taxonomically diverse TrEMBL
entries if needed to clear n=50, aligning the result with MAFFT, and
computing real per-column statistics (Shannon entropy, %identity to human
WT, the observed substitution spectrum, and mean Grantham physicochemical
distance) for every position -- not just the four resistance-mutation
sites, so those four can be ranked against a real distribution instead of
reported in isolation.
"""

import math
import re
import shutil
import subprocess
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import config  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
import requests  # noqa: E402
from Bio import SeqIO  # noqa: E402

# ---------------------------------------------------------------------------
# UniProt REST access
# ---------------------------------------------------------------------------
UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
UNIPROT_FIELDS = "accession,organism_name,length,lineage,sequence"
PAGE_SIZE = 500
REQUEST_TIMEOUT = 300
MAX_RETRIES = 4
NEXT_LINK_RE = re.compile(r'<([^>]+)>;\s*rel="next"')

# ---------------------------------------------------------------------------
# Dataset construction
# ---------------------------------------------------------------------------
MIN_SEQ_LEN = 150
MAX_SEQ_LEN = 260
MIN_TOTAL_SEQUENCES = 50
MAX_TOTAL_SEQUENCES = 200
GENUS_CAP = 3
FASTA_WRAP = 70

EXPANDED_FASTA_OUT = config.SEQUENCES_DIR / "dhfr_expanded.fasta"
EXCLUSION_LOG_OUT = config.RESULTS_CONSERVATION_DIR / "exclusion_log.txt"
TAXA_METADATA_OUT = config.RESULTS_CONSERVATION_DIR / "taxa_metadata.csv"
ALIGNED_FASTA_OUT = config.RESULTS_CONSERVATION_DIR / "dhfr_expanded_aligned.fasta"
SITE_TABLE_OUT = config.RESULTS_CONSERVATION_DIR / "site_conservation_table.csv"
FULL_PROFILE_OUT = config.RESULTS_CONSERVATION_DIR / "full_conservation_profile.csv"
FIGURE_OUT = config.RESULTS_CONSERVATION_DIR / "conservation_profile.png"
SUMMARY_OUT = config.RESULTS_CONSERVATION_DIR / "conservation_summary.txt"

# ---------------------------------------------------------------------------
# Human numbering: our fetched UniProt sequences retain Met1, but the
# resistance-mutation literature (and config.RESISTANCE_MUTANT_PDBS) uses
# mature-protein/PDB numbering (Met1 cleaved) -- same convention, same
# +1 offset, as scripts/03_align_mafft.py. Verified against the human
# sequence directly: raw position 23 = L (Leu22 mature), 32 = F (Phe31),
# 36 = Q (Gln35), 65 = N (Asn64).
# ---------------------------------------------------------------------------
HUMAN_NUMBERING_OFFSET = 1
FLAGGED_POSITIONS = {22: "L", 31: "F", 35: "Q", 64: "N"}

# Known methotrexate-resistance substitutions at three of the four flagged
# sites (literature set named in the task; no established resistance
# mutation is catalogued at position 64 in this pipeline's reference set --
# see config.RESISTANCE_MUTANT_PDBS's note on 3F8Z/N64S being non-significant).
KNOWN_RESISTANCE_MUTATIONS = [
    (22, "L", "R"), (22, "L", "F"), (22, "L", "Y"), (22, "L", "W"),
    (31, "F", "R"),
    (35, "Q", "E"), (35, "Q", "S"), (35, "Q", "K"),
]

# ---------------------------------------------------------------------------
# Grantham (1974) physicochemical distance matrix -- composition, polarity,
# and molecular volume of amino acid side chains. Sourced from
# patterninstitute/grantham (github.com/patterninstitute/grantham,
# data-raw/grantham_distance_matrix.csv), cross-checked for symmetry and
# zero diagonal before use.
# ---------------------------------------------------------------------------
GRANTHAM_MATRIX = {
    "A": {"A": 0, "C": 195, "D": 126, "E": 107, "F": 113, "G": 60, "H": 86, "I": 94, "K": 106, "L": 96, "M": 84, "N": 111, "P": 27, "Q": 91, "R": 112, "S": 99, "T": 58, "V": 64, "W": 148, "Y": 112},
    "C": {"A": 195, "C": 0, "D": 154, "E": 170, "F": 205, "G": 159, "H": 174, "I": 198, "K": 202, "L": 198, "M": 196, "N": 139, "P": 169, "Q": 154, "R": 180, "S": 112, "T": 149, "V": 192, "W": 215, "Y": 194},
    "D": {"A": 126, "C": 154, "D": 0, "E": 45, "F": 177, "G": 94, "H": 81, "I": 168, "K": 101, "L": 172, "M": 160, "N": 23, "P": 108, "Q": 61, "R": 96, "S": 65, "T": 85, "V": 152, "W": 181, "Y": 160},
    "E": {"A": 107, "C": 170, "D": 45, "E": 0, "F": 140, "G": 98, "H": 40, "I": 134, "K": 56, "L": 138, "M": 126, "N": 42, "P": 93, "Q": 29, "R": 54, "S": 80, "T": 65, "V": 121, "W": 152, "Y": 122},
    "F": {"A": 113, "C": 205, "D": 177, "E": 140, "F": 0, "G": 153, "H": 100, "I": 21, "K": 102, "L": 22, "M": 28, "N": 158, "P": 114, "Q": 116, "R": 97, "S": 155, "T": 103, "V": 50, "W": 40, "Y": 22},
    "G": {"A": 60, "C": 159, "D": 94, "E": 98, "F": 153, "G": 0, "H": 98, "I": 135, "K": 127, "L": 138, "M": 127, "N": 80, "P": 42, "Q": 87, "R": 125, "S": 56, "T": 59, "V": 109, "W": 184, "Y": 147},
    "H": {"A": 86, "C": 174, "D": 81, "E": 40, "F": 100, "G": 98, "H": 0, "I": 94, "K": 32, "L": 99, "M": 87, "N": 68, "P": 77, "Q": 24, "R": 29, "S": 89, "T": 47, "V": 84, "W": 115, "Y": 83},
    "I": {"A": 94, "C": 198, "D": 168, "E": 134, "F": 21, "G": 135, "H": 94, "I": 0, "K": 102, "L": 5, "M": 10, "N": 149, "P": 95, "Q": 109, "R": 97, "S": 142, "T": 89, "V": 29, "W": 61, "Y": 33},
    "K": {"A": 106, "C": 202, "D": 101, "E": 56, "F": 102, "G": 127, "H": 32, "I": 102, "K": 0, "L": 107, "M": 95, "N": 94, "P": 103, "Q": 53, "R": 26, "S": 121, "T": 78, "V": 97, "W": 110, "Y": 85},
    "L": {"A": 96, "C": 198, "D": 172, "E": 138, "F": 22, "G": 138, "H": 99, "I": 5, "K": 107, "L": 0, "M": 15, "N": 153, "P": 98, "Q": 113, "R": 102, "S": 145, "T": 92, "V": 32, "W": 61, "Y": 36},
    "M": {"A": 84, "C": 196, "D": 160, "E": 126, "F": 28, "G": 127, "H": 87, "I": 10, "K": 95, "L": 15, "M": 0, "N": 142, "P": 87, "Q": 101, "R": 91, "S": 135, "T": 81, "V": 21, "W": 67, "Y": 36},
    "N": {"A": 111, "C": 139, "D": 23, "E": 42, "F": 158, "G": 80, "H": 68, "I": 149, "K": 94, "L": 153, "M": 142, "N": 0, "P": 91, "Q": 46, "R": 86, "S": 46, "T": 65, "V": 133, "W": 174, "Y": 143},
    "P": {"A": 27, "C": 169, "D": 108, "E": 93, "F": 114, "G": 42, "H": 77, "I": 95, "K": 103, "L": 98, "M": 87, "N": 91, "P": 0, "Q": 76, "R": 103, "S": 74, "T": 38, "V": 68, "W": 147, "Y": 110},
    "Q": {"A": 91, "C": 154, "D": 61, "E": 29, "F": 116, "G": 87, "H": 24, "I": 109, "K": 53, "L": 113, "M": 101, "N": 46, "P": 76, "Q": 0, "R": 43, "S": 68, "T": 42, "V": 96, "W": 130, "Y": 99},
    "R": {"A": 112, "C": 180, "D": 96, "E": 54, "F": 97, "G": 125, "H": 29, "I": 97, "K": 26, "L": 102, "M": 91, "N": 86, "P": 103, "Q": 43, "R": 0, "S": 110, "T": 71, "V": 96, "W": 101, "Y": 77},
    "S": {"A": 99, "C": 112, "D": 65, "E": 80, "F": 155, "G": 56, "H": 89, "I": 142, "K": 121, "L": 145, "M": 135, "N": 46, "P": 74, "Q": 68, "R": 110, "S": 0, "T": 58, "V": 124, "W": 177, "Y": 144},
    "T": {"A": 58, "C": 149, "D": 85, "E": 65, "F": 103, "G": 59, "H": 47, "I": 89, "K": 78, "L": 92, "M": 81, "N": 65, "P": 38, "Q": 42, "R": 71, "S": 58, "T": 0, "V": 69, "W": 128, "Y": 92},
    "V": {"A": 64, "C": 192, "D": 152, "E": 121, "F": 50, "G": 109, "H": 84, "I": 29, "K": 97, "L": 32, "M": 21, "N": 133, "P": 68, "Q": 96, "R": 96, "S": 124, "T": 69, "V": 0, "W": 88, "Y": 55},
    "W": {"A": 148, "C": 215, "D": 181, "E": 152, "F": 40, "G": 184, "H": 115, "I": 61, "K": 110, "L": 61, "M": 67, "N": 174, "P": 147, "Q": 130, "R": 101, "S": 177, "T": 128, "V": 88, "W": 0, "Y": 37},
    "Y": {"A": 112, "C": 194, "D": 160, "E": 122, "F": 22, "G": 147, "H": 83, "I": 33, "K": 85, "L": 36, "M": 36, "N": 143, "P": 110, "Q": 99, "R": 77, "S": 144, "T": 92, "V": 55, "W": 37, "Y": 0},
}


# ---------------------------------------------------------------------------
# Step 1-3: UniProt retrieval
# ---------------------------------------------------------------------------

def uniprot_paginated_search(query, label, session):
    """Fetches every hit for `query`, following UniProt's cursor-based Link
    header pagination. Returns a list of raw per-entry JSON dicts."""
    url = UNIPROT_SEARCH_URL
    params = {"query": query, "format": "json", "size": PAGE_SIZE, "fields": UNIPROT_FIELDS}
    results = []
    page = 0
    while url:
        page += 1
        attempt = 0
        while True:
            attempt += 1
            try:
                resp = session.get(url, params=params, timeout=REQUEST_TIMEOUT)
                resp.raise_for_status()
                break
            except requests.exceptions.RequestException as exc:
                if attempt > MAX_RETRIES:
                    raise RuntimeError(
                        f"UniProt query ({label}) failed on page {page} after "
                        f"{MAX_RETRIES} retries: {exc}"
                    ) from exc
                wait = min(60, 5 * (2 ** (attempt - 1)))
                print(f"    [{label} p{page}] request failed ({exc}); retrying in {wait}s "
                      f"(attempt {attempt}/{MAX_RETRIES})...")
                time.sleep(wait)

        batch = resp.json().get("results", [])
        results.extend(batch)
        print(f"    [{label}] page {page}: +{len(batch)} entries (running total {len(results)})")

        params = None  # the next link already carries the full query string
        match = NEXT_LINK_RE.search(resp.headers.get("Link", "") or "")
        url = match.group(1) if match else None
        if url:
            time.sleep(1.0)
    return results


def parse_entry(raw, reviewed):
    organism = raw.get("organism", {})
    sequence = raw.get("sequence", {})
    return {
        "accession": raw.get("primaryAccession"),
        "organism": organism.get("scientificName"),
        "lineage": organism.get("lineage", []),
        "length": sequence.get("length"),
        "sequence": sequence.get("value"),
        "reviewed": reviewed,
    }


def filter_by_length(entries, min_len=MIN_SEQ_LEN, max_len=MAX_SEQ_LEN):
    kept, excluded = [], []
    for e in entries:
        length, seq = e.get("length"), e.get("sequence")
        if length is None or not seq:
            excluded.append({**e, "exclusion_reason": "missing length/sequence data"})
        elif length < min_len:
            excluded.append({**e, "exclusion_reason": f"too short (<{min_len} aa)"})
        elif length > max_len:
            excluded.append({**e, "exclusion_reason": f"too long (>{max_len} aa, likely multi-domain/fusion)"})
        else:
            kept.append(e)
    return kept, excluded


def genus_of(organism_name):
    return (organism_name or "?").split()[0]


def select_diverse_unreviewed(candidates, exclude_organisms, target_total, genus_cap=GENUS_CAP):
    """Round-robins across genera (alphabetically, deterministic) so no
    single well-studied genus (e.g. many Mus musculus strain records)
    crowds out taxonomic breadth, capping each genus at `genus_cap`."""
    by_organism = {}
    for e in sorted(candidates, key=lambda e: (e["organism"] or "", e["accession"] or "")):
        org = e["organism"]
        if not org or org in exclude_organisms or org in by_organism:
            continue
        by_organism[org] = e
    unique_candidates = list(by_organism.values())

    by_genus = defaultdict(list)
    for e in unique_candidates:
        by_genus[genus_of(e["organism"])].append(e)
    for genus_entries in by_genus.values():
        genus_entries.sort(key=lambda e: e["accession"] or "")

    genera_sorted = sorted(by_genus.keys())
    idx_per_genus = {g: 0 for g in genera_sorted}
    genus_counts = Counter()
    selected = []
    progressed = True
    while progressed and len(selected) < target_total:
        progressed = False
        for genus in genera_sorted:
            if len(selected) >= target_total:
                break
            if genus_counts[genus] >= genus_cap:
                continue
            i = idx_per_genus[genus]
            if i >= len(by_genus[genus]):
                continue
            selected.append(by_genus[genus][i])
            idx_per_genus[genus] += 1
            genus_counts[genus] += 1
            progressed = True
    return selected


def dedup_by_organism(entries):
    """One sequence per organism. Reviewed entries win ties (shouldn't
    occur in practice since unreviewed candidates already exclude
    reviewed-set organisms, but kept as a defensive guarantee)."""
    by_org = {}
    for e in sorted(entries, key=lambda e: (e["organism"] or "", not e["reviewed"], e["accession"] or "")):
        by_org.setdefault(e["organism"], e)
    return list(by_org.values())


def write_exclusion_log(reviewed_all, reviewed_kept, reviewed_excluded,
                         unreviewed_all, unreviewed_excluded, final_entries):
    lines = []
    lines.append("DHFR expanded conservation dataset -- sequence exclusion log")
    lines.append("=" * 70)
    lines.append(f"Query: gene:{config.GENE_SYMBOL} AND reviewed:true (Swiss-Prot)")
    lines.append(f"Reviewed hits fetched: {len(reviewed_all)}")
    lines.append("")
    lines.append(f"Length filter: keep {MIN_SEQ_LEN}-{MAX_SEQ_LEN} aa (single-domain, monomeric DHFR).")
    lines.append(
        "This deliberately excludes bifunctional DHFR-thymidylate synthase (DHFR-TS)"
    )
    lines.append(
        "fusion proteins found in plants, kinetoplastid parasites (e.g. Leishmania,"
    )
    lines.append(
        "Trypanosoma), and apicomplexan parasites (e.g. Plasmodium) -- these run"
    )
    lines.append(
        "roughly 2-4x longer than a monofunctional DHFR domain and would badly"
    )
    lines.append("distort a global multiple sequence alignment if included alongside single-domain sequences.")
    lines.append("")
    lines.append("Reviewed (Swiss-Prot) set:")
    lines.append(f"  kept ({MIN_SEQ_LEN}-{MAX_SEQ_LEN} aa):  {len(reviewed_kept)}")
    n_short = sum(1 for e in reviewed_excluded if "too short" in e["exclusion_reason"])
    n_long = sum(1 for e in reviewed_excluded if "too long" in e["exclusion_reason"])
    n_bad = len(reviewed_excluded) - n_short - n_long
    lines.append(f"  excluded, too short (<{MIN_SEQ_LEN} aa):  {n_short}")
    lines.append(f"  excluded, too long (>{MAX_SEQ_LEN} aa, likely fusion/multi-domain):  {n_long}")
    if n_bad:
        lines.append(f"  excluded, missing data:  {n_bad}")
    if reviewed_excluded:
        lines.append("  examples of excluded reviewed entries (accession, organism, length):")
        for e in sorted(reviewed_excluded, key=lambda e: -(e.get("length") or 0))[:10]:
            lines.append(f"    {e['accession']}  {e['organism']}  {e.get('length')} aa  ({e['exclusion_reason']})")
    lines.append("")

    if unreviewed_all is not None:
        lines.append(
            f"Reviewed-only kept count ({len(reviewed_kept)}) was below the minimum of "
            f"{MIN_TOTAL_SEQUENCES}, so unreviewed (TrEMBL) entries were also queried"
        )
        lines.append(f"(gene:{config.GENE_SYMBOL}, reviewed:true filter dropped) and length-filtered the same way:")
        lines.append(f"  unreviewed hits fetched: {len(unreviewed_all)}")
        n_kept_unrev = len(unreviewed_all) - len(unreviewed_excluded)
        n_short_u = sum(1 for e in unreviewed_excluded if "too short" in e["exclusion_reason"])
        n_long_u = sum(1 for e in unreviewed_excluded if "too long" in e["exclusion_reason"])
        n_bad_u = len(unreviewed_excluded) - n_short_u - n_long_u
        lines.append(f"  in length range ({MIN_SEQ_LEN}-{MAX_SEQ_LEN} aa): {n_kept_unrev}")
        lines.append(f"  excluded, too short (<{MIN_SEQ_LEN} aa): {n_short_u}")
        lines.append(f"  excluded, too long (>{MAX_SEQ_LEN} aa, likely fusion/multi-domain): {n_long_u}")
        if n_bad_u:
            lines.append(f"  excluded, missing data: {n_bad_u}")
        if unreviewed_excluded:
            lines.append("  examples of excluded unreviewed entries (accession, organism, length):")
            for e in sorted(unreviewed_excluded, key=lambda e: -(e.get("length") or 0))[:10]:
                lines.append(f"    {e['accession']}  {e['organism']}  {e.get('length')} aa  ({e['exclusion_reason']})")
        lines.append(
            f"  Of those in range, entries were added prioritizing taxonomic diversity"
        )
        lines.append(
            f"  (round-robin across genera, capped at {GENUS_CAP} sequences/genus, so no single"
        )
        lines.append(
            "  well-studied genus -- e.g. multiple Mus musculus strain records -- crowds out breadth)."
        )
        lines.append("")
    else:
        lines.append(
            f"Reviewed-only kept count ({len(reviewed_kept)}) already met the minimum of "
            f"{MIN_TOTAL_SEQUENCES}; unreviewed (TrEMBL) entries were not queried."
        )
        lines.append("")

    lines.append(
        f"Final dataset after deduplication by organism: {len(final_entries)} sequences "
        f"(target range {MIN_TOTAL_SEQUENCES}-{MAX_TOTAL_SEQUENCES})."
    )

    EXCLUSION_LOG_OUT.write_text("\n".join(lines) + "\n")


# ---------------------------------------------------------------------------
# Step 4: FASTA + taxa metadata
# ---------------------------------------------------------------------------

def classify_taxon(lineage):
    lineage_set = set(lineage or [])
    if "Mammalia" in lineage_set:
        return "mammal"
    if "Vertebrata" in lineage_set:
        return "other vertebrate"
    if "Metazoa" in lineage_set:
        return "invertebrate"
    if "Fungi" in lineage_set:
        return "fungi"
    if "Bacteria" in lineage_set:
        return "bacteria"
    return "other"


def sort_for_output(entries):
    def key(e):
        return (0 if e["organism"] == "Homo sapiens" else 1, e["organism"] or "")
    return sorted(entries, key=key)


def write_fasta(entries, out_path):
    with open(out_path, "w") as fh:
        for e in entries:
            header = e["organism"].replace(" ", "_")
            fh.write(f">{header}|{e['accession']}\n")
            seq = e["sequence"]
            for i in range(0, len(seq), FASTA_WRAP):
                fh.write(seq[i:i + FASTA_WRAP] + "\n")


def write_taxa_metadata(entries, out_path):
    df = pd.DataFrame([
        {
            "accession": e["accession"],
            "organism": e["organism"],
            "taxonomic_group": classify_taxon(e["lineage"]),
        }
        for e in entries
    ])
    df.to_csv(out_path, index=False)
    return df


# ---------------------------------------------------------------------------
# Step 5: MAFFT alignment
# ---------------------------------------------------------------------------

def run_mafft(input_fasta, output_fasta):
    if not shutil.which("mafft"):
        print("ERROR: 'mafft' not found on PATH. Install it (e.g. `brew install mafft` "
              "on macOS, `sudo apt-get install mafft` on Linux) and re-run.")
        sys.exit(1)
    with open(output_fasta, "w") as out_f:
        result = subprocess.run(
            ["mafft", "--auto", str(input_fasta)],
            stdout=out_f, stderr=subprocess.PIPE, text=True,
        )
    if result.returncode != 0:
        raise RuntimeError(f"mafft exited {result.returncode}: {result.stderr}")
    records = list(SeqIO.parse(output_fasta, "fasta"))
    if not records:
        raise RuntimeError(f"Alignment produced no sequences ({output_fasta})")
    return records


# ---------------------------------------------------------------------------
# Steps 6-7: per-column conservation statistics
# ---------------------------------------------------------------------------

def shannon_entropy(counts):
    total = sum(counts.values())
    if total == 0:
        return float("nan")
    entropy = 0.0
    for c in counts.values():
        p = c / total
        entropy -= p * math.log2(p)
    return entropy


def format_observed(counts, total):
    parts = [f"{aa}:{c} ({100 * c / total:.1f}%)" for aa, c in counts.most_common()]
    return "; ".join(parts)


def mean_grantham_distance(counts, wt_residue):
    if wt_residue not in GRANTHAM_MATRIX:
        return float("nan")
    diffs_total, diffs_n = 0, 0
    for aa, c in counts.items():
        if aa == wt_residue or aa not in GRANTHAM_MATRIX:
            continue
        diffs_total += GRANTHAM_MATRIX[wt_residue][aa] * c
        diffs_n += c
    return (diffs_total / diffs_n) if diffs_n else float("nan")


def column_stats(records, col_idx, wt_residue):
    chars = [str(r.seq)[col_idx] for r in records]
    non_gap = [c for c in chars if c != "-"]
    n = len(non_gap)
    counts = Counter(non_gap)
    return {
        "shannon_entropy": shannon_entropy(counts),
        "pct_identity_to_wt": (100.0 * counts.get(wt_residue, 0) / n) if n else float("nan"),
        "observed_residues": format_observed(counts, n) if n else "",
        "mean_grantham_distance": mean_grantham_distance(counts, wt_residue),
        "n_sequences": n,
    }


def human_column_to_mature_position(human_aligned_seq):
    """Maps every non-gap alignment column of the aligned human sequence to
    its mature-protein (Met1-cleaved) residue number. Raw position 1 (Met1
    itself) has no mature-numbering equivalent and is skipped."""
    mapping = {}
    raw_pos = 0
    for col, ch in enumerate(human_aligned_seq):
        if ch == "-":
            continue
        raw_pos += 1
        mature_pos = raw_pos - HUMAN_NUMBERING_OFFSET
        if mature_pos >= 1:
            mapping[col] = mature_pos
    return mapping


def compute_conservation_tables(records):
    human = next((r for r in records if r.id.startswith("Homo_sapiens|")), None)
    if human is None:
        print("ERROR: no Homo_sapiens record found in the alignment -- cannot map "
              "human-numbered positions.")
        sys.exit(1)
    human_seq = str(human.seq)

    col_to_mature = human_column_to_mature_position(human_seq)
    mature_to_col = {v: k for k, v in col_to_mature.items()}

    site_rows = []
    for lit_pos, expected_wt in sorted(FLAGGED_POSITIONS.items()):
        col = mature_to_col.get(lit_pos)
        if col is None:
            print(f"WARNING: could not map flagged position {lit_pos} to an alignment column "
                  "(beyond human sequence length?) -- skipping.")
            continue
        wt = human_seq[col]
        if wt != expected_wt:
            print(f"WARNING: flagged position {lit_pos} -- expected WT {expected_wt}, "
                  f"found {wt} in the fetched human sequence.")
        stats = column_stats(records, col, wt)
        site_rows.append({"position": lit_pos, "wt_residue": wt, **stats})
    site_df = pd.DataFrame(site_rows)

    full_rows = []
    for col, mature_pos in sorted(col_to_mature.items(), key=lambda kv: kv[1]):
        if mature_pos in FLAGGED_POSITIONS:
            continue
        wt = human_seq[col]
        stats = column_stats(records, col, wt)
        full_rows.append({"position": mature_pos, "wt_residue": wt, **stats})
    full_df = pd.DataFrame(full_rows)

    for df in (site_df, full_df):
        df["shannon_entropy"] = df["shannon_entropy"].round(4)
        df["pct_identity_to_wt"] = df["pct_identity_to_wt"].round(2)
        df["mean_grantham_distance"] = df["mean_grantham_distance"].round(2)

    return site_df, full_df, len(records)


# ---------------------------------------------------------------------------
# Step 8: figure
# ---------------------------------------------------------------------------

def render_figure(site_df, full_df, n_sequences, out_path):
    combined = pd.concat([site_df, full_df], ignore_index=True).sort_values("position")
    mean_entropy = full_df["shannon_entropy"].mean()

    fig, ax = plt.subplots(figsize=(14, 5.5))
    ax.plot(combined["position"], combined["shannon_entropy"], color="#4C72B0",
            linewidth=1.1, zorder=2)
    ax.axhline(mean_entropy, color="gray", linestyle=":", linewidth=1.2, zorder=1,
               label=f"Mean entropy, all other positions = {mean_entropy:.2f} bits")

    # Stagger labels vertically when two flagged positions fall close together
    # on the x-axis (e.g. 31 and 35 are only 4 residues apart) so the bold
    # "F31"/"Q35"-style labels don't overlap.
    y_top = combined["shannon_entropy"].max()
    x_span = combined["position"].max() - combined["position"].min() or 1
    close_threshold = x_span * 0.03
    prev_pos, depth = None, 0
    for _, row in site_df.sort_values("position").iterrows():
        pos, ent = row["position"], row["shannon_entropy"]
        depth = depth + 1 if (prev_pos is not None and pos - prev_pos < close_threshold) else 0
        prev_pos = pos
        label_y = y_top * 1.08 + 0.15 + depth * (y_top * 0.14 + 0.1)
        ax.scatter([pos], [ent], color="crimson", zorder=3, s=45, edgecolor="white", linewidth=0.6)
        ax.annotate(
            f"{row['wt_residue']}{pos}",
            xy=(pos, ent), xytext=(pos, label_y),
            ha="center", fontsize=9, fontweight="bold", color="crimson",
            arrowprops=dict(arrowstyle="-", color="crimson", linewidth=0.8, alpha=0.7),
        )

    ax.set_ylim(bottom=-0.05, top=y_top * 1.22 + 0.2 + depth * (y_top * 0.14 + 0.1))
    ax.set_xlabel("Residue position (human / mature-protein numbering)")
    ax.set_ylabel("Shannon entropy (bits)")
    ax.set_title(
        f"DHFR alignment conservation profile -- {n_sequences} sequences, "
        f"{len(combined)} human-mapped positions"
    )
    ax.legend(loc="upper left", fontsize=8, frameon=True)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Step 9: summary
# ---------------------------------------------------------------------------

def conservation_percentile(entropy_value, reference_entropies):
    """% of reference positions at least as variable (>= entropy) as this
    one -- i.e. how conserved this position is relative to the rest of the
    protein. 100 = the single most conserved position; 0 = the least."""
    n = len(reference_entropies)
    if n == 0:
        return float("nan")
    n_at_least_as_variable = sum(1 for e in reference_entropies if e >= entropy_value)
    return 100.0 * n_at_least_as_variable / n


def describe_relative(value, reference_mean):
    """Lower Shannon entropy = more conserved."""
    if abs(value - reference_mean) < 1e-9:
        return "about as conserved as"
    return "more conserved than" if value < reference_mean else "less conserved than"


def write_summary(site_df, full_df, taxa_df, n_sequences, aln_len, out_path):
    reference_entropies = full_df["shannon_entropy"].tolist()
    protein_mean_entropy = full_df["shannon_entropy"].mean()
    n_groups = taxa_df["taxonomic_group"].nunique()

    site_by_pos = site_df.set_index("position")

    sentences = []

    sentences.append(
        f"The expanded dataset comprises {n_sequences} DHFR sequences spanning "
        f"{n_groups} taxonomic groups ({', '.join(sorted(taxa_df['taxonomic_group'].unique()))}), "
        f"aligned to {aln_len} columns -- versus 6 sequences previously -- giving each "
        f"per-position statistic real statistical support instead of a handful of "
        f"discrete possible values."
    )

    conservation_clauses = []
    for pos in sorted(FLAGGED_POSITIONS):
        if pos not in site_by_pos.index:
            continue
        row = site_by_pos.loc[pos]
        pct = conservation_percentile(row["shannon_entropy"], reference_entropies)
        rel = describe_relative(row["shannon_entropy"], protein_mean_entropy)
        conservation_clauses.append(
            f"position {pos} ({row['wt_residue']}{pos}: entropy {row['shannon_entropy']:.2f} bits, "
            f"{row['pct_identity_to_wt']:.1f}% identity to WT, {rel} the protein average "
            f"of {protein_mean_entropy:.2f} bits -- {pct:.0f}th percentile of conservation "
            f"across the other {len(full_df)} aligned positions)"
        )
    sentences.append(
        "Relative to the rest of the protein, " + "; ".join(conservation_clauses) + "."
    )

    # Natural substitution character vs. the known resistance mutations at
    # the three sites the literature set actually names (22, 31, 35).
    mutations_by_pos = defaultdict(list)
    for pos, wt, mut in KNOWN_RESISTANCE_MUTATIONS:
        mutations_by_pos[pos].append((wt, mut, GRANTHAM_MATRIX[wt][mut]))

    natural_clauses = []
    verdicts = {}  # pos -> "smaller" | "comparable" | "larger" | None (no natural variation)
    for pos in (22, 31, 35):
        if pos not in site_by_pos.index:
            continue
        row = site_by_pos.loc[pos]
        natural_dist = row["mean_grantham_distance"]
        res_dists = mutations_by_pos[pos]
        res_dist_str = ", ".join(f"{wt}{pos}{mut}={d}" for wt, mut, d in res_dists)
        res_mean = sum(d for _, _, d in res_dists) / len(res_dists)
        if pd.isna(natural_dist):
            verdicts[pos] = None
            natural_clauses.append(
                f"position {pos} shows no natural substitutions in this dataset (100% identity), "
                f"so no physicochemical comparison to the resistance mutations ({res_dist_str}) is possible"
            )
        else:
            if natural_dist < res_mean * 0.7:
                verdict, verdicts[pos] = "substantially smaller than", "smaller"
            elif natural_dist > res_mean * 1.3:
                verdict, verdicts[pos] = "as large as or larger than", "larger"
            else:
                verdict, verdicts[pos] = "comparable in magnitude to", "comparable"
            natural_clauses.append(
                f"position {pos} tolerates natural substitutions averaging Grantham distance "
                f"{natural_dist:.0f} ({row['observed_residues']}), {verdict} the mean Grantham "
                f"distance of the known resistance mutations there ({res_dist_str}, mean {res_mean:.0f})"
            )

    sentences.append(
        "Comparing the physicochemical size of tolerated natural variation to the known "
        "resistance mutations (L22R/F/Y/W, F31R, Q35E/S/K): " + "; ".join(natural_clauses) + "."
    )

    pos64 = site_by_pos.loc[64] if 64 in site_by_pos.index else None
    if pos64 is not None:
        if pd.isna(pos64["mean_grantham_distance"]):
            pos64_clause = (
                f"Position 64 ({pos64['wt_residue']}64) is 100% identical to the human WT residue "
                f"across all {int(pos64['n_sequences'])} sequences in this dataset, and no "
                f"established resistance mutation is catalogued there in this pipeline's reference "
                f"set (the N64S variant tested alongside Q35S in PDB 3F8Z was reported as not "
                f"statistically significant)."
            )
        else:
            pos64_clause = (
                f"Position 64 ({pos64['wt_residue']}64) tolerates natural substitutions averaging "
                f"Grantham distance {pos64['mean_grantham_distance']:.0f} ({pos64['observed_residues']}); "
                f"no established resistance mutation is catalogued there in this pipeline's reference "
                f"set (the N64S variant tested alongside Q35S in PDB 3F8Z was reported as not "
                f"statistically significant)."
            )
        sentences.append(pos64_clause)

    # Build the closing verdict from what was actually found at 22/31/35,
    # rather than collapsing three independent comparisons into one
    # any()/else -- a single "substantially smaller" hit next to two
    # "comparable" ones is a mixed result, not uniform support for either
    # extreme, and the summary should say so.
    smaller_pos = [p for p, v in verdicts.items() if v == "smaller"]
    comparable_pos = [p for p, v in verdicts.items() if v == "comparable"]
    larger_pos = [p for p, v in verdicts.items() if v == "larger"]
    scored_pos = smaller_pos + comparable_pos + larger_pos

    if not scored_pos:
        conclusion = (
            "Overall, none of positions 22, 31, or 35 show any natural substitution in this dataset "
            "(all 100% identical to human WT), so this expanded alignment cannot speak to whether "
            "nature tolerates smaller physicochemical changes there than the resistance mutations do."
        )
    elif smaller_pos and not comparable_pos and not larger_pos:
        conclusion = (
            f"Overall, at every site with natural variation (position{'s' if len(smaller_pos) > 1 else ''} "
            f"{', '.join(map(str, smaller_pos))}), evolution tolerates substitutions that are "
            "physicochemically smaller than the resistance mutations seen there -- consistent with these "
            "residues being under functional constraint that the resistance mutations specifically violate."
        )
    elif larger_pos and not smaller_pos:
        conclusion = (
            "Overall, where natural substitutions are observed at these sites their physicochemical "
            "magnitude is not smaller than the known resistance mutations, so conservation alone does "
            "not cleanly separate 'tolerated' from 'resistance-causing' changes at every flagged "
            "position -- structural/functional context (Modules 5-6) remains necessary."
        )
    else:
        def position_list(positions):
            return " and ".join(map(str, positions)) if len(positions) <= 2 else (
                ", ".join(map(str, positions[:-1])) + f", and {positions[-1]}"
            )

        clause_bits = []
        if smaller_pos:
            clause_bits.append(
                f"at position{'s' if len(smaller_pos) > 1 else ''} {position_list(smaller_pos)}, "
                f"natural variation is clearly milder than the resistance mutation"
                f"{'s' if len(smaller_pos) > 1 else ''} there"
            )
        if comparable_pos:
            clause_bits.append(
                f"at position{'s' if len(comparable_pos) > 1 else ''} {position_list(comparable_pos)}, "
                f"natural variation is roughly the same physicochemical size as the resistance mutation"
                f"{'s' if len(comparable_pos) > 1 else ''} there"
            )
        if larger_pos:
            clause_bits.append(
                f"at position{'s' if len(larger_pos) > 1 else ''} {position_list(larger_pos)}, "
                f"natural variation is at least as large as the resistance mutation"
                f"{'s' if len(larger_pos) > 1 else ''} there"
            )
        conclusion = (
            "Overall the picture is mixed rather than uniform: " + "; ".join(clause_bits) +
            " -- so conservation only partly explains why the specific resistance substitutions are "
            "disruptive, and structural/functional context (Modules 5-6) is still needed for the rest."
        )
    sentences.append(conclusion)

    out_path.write_text(" ".join(sentences) + "\n")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=" * 90)
    print("Module 8: Expanded, properly-powered conservation analysis")
    print("=" * 90)

    config.RESULTS_CONSERVATION_DIR.mkdir(parents=True, exist_ok=True)
    config.SEQUENCES_DIR.mkdir(parents=True, exist_ok=True)

    session = requests.Session()
    session.headers.update({"User-Agent": f"dhfr-pipeline/1.0 ({config.ENTREZ_EMAIL})"})

    # --- Step 1: reviewed entries -----------------------------------------
    print(f"\nStep 1: querying UniProt for reviewed gene:{config.GENE_SYMBOL} entries...")
    reviewed_raw = uniprot_paginated_search(
        f"gene:{config.GENE_SYMBOL} AND reviewed:true", "reviewed", session
    )
    reviewed_all = [parse_entry(r, reviewed=True) for r in reviewed_raw]
    print(f"  {len(reviewed_all)} reviewed entries fetched")

    # --- Step 2: length filter ---------------------------------------------
    print(f"\nStep 2: filtering to single-domain length range {MIN_SEQ_LEN}-{MAX_SEQ_LEN} aa...")
    reviewed_kept, reviewed_excluded = filter_by_length(reviewed_all)
    print(f"  kept {len(reviewed_kept)} / {len(reviewed_all)} reviewed entries "
          f"({len(reviewed_excluded)} excluded)")

    # --- Step 3: top up with unreviewed entries if needed ------------------
    unreviewed_all, unreviewed_excluded = None, []
    unreviewed_selected = []
    if len(reviewed_kept) < MIN_TOTAL_SEQUENCES:
        print(f"\nStep 3: reviewed-only kept count ({len(reviewed_kept)}) < {MIN_TOTAL_SEQUENCES}; "
              "adding unreviewed (TrEMBL) entries...")
        unreviewed_raw = uniprot_paginated_search(
            f"gene:{config.GENE_SYMBOL} AND reviewed:false", "unreviewed", session
        )
        unreviewed_all = [parse_entry(r, reviewed=False) for r in unreviewed_raw]
        print(f"  {len(unreviewed_all)} unreviewed entries fetched")
        unreviewed_in_range, unreviewed_excluded = filter_by_length(unreviewed_all)
        print(f"  {len(unreviewed_in_range)} / {len(unreviewed_all)} in length range")

        reviewed_organisms = {e["organism"] for e in reviewed_kept}
        target_total = MAX_TOTAL_SEQUENCES - len(reviewed_kept)
        unreviewed_selected = select_diverse_unreviewed(
            unreviewed_in_range, exclude_organisms=reviewed_organisms, target_total=target_total,
        )
        n_genera = len({genus_of(e["organism"]) for e in unreviewed_selected})
        print(f"  selected {len(unreviewed_selected)} taxonomically diverse unreviewed entries "
              f"({n_genera} genera, cap {GENUS_CAP}/genus)")
    else:
        print(f"\nStep 3: reviewed-only kept count ({len(reviewed_kept)}) already >= "
              f"{MIN_TOTAL_SEQUENCES}; skipping unreviewed entries.")

    write_exclusion_log(reviewed_all, reviewed_kept, reviewed_excluded,
                         unreviewed_all, unreviewed_excluded,
                         final_entries=dedup_by_organism(reviewed_kept + unreviewed_selected))
    print(f"  Wrote exclusion log: {EXCLUSION_LOG_OUT}")

    # --- Step 4: dedup, FASTA, metadata ------------------------------------
    print("\nStep 4: deduplicating by organism, writing FASTA + taxa metadata...")
    final_entries = dedup_by_organism(reviewed_kept + unreviewed_selected)
    if not any(e["organism"] == "Homo sapiens" for e in final_entries):
        print("ERROR: Homo sapiens dropped out of the final dataset -- cannot proceed "
              "without the human WT reference sequence.")
        sys.exit(1)
    if not (MIN_TOTAL_SEQUENCES <= len(final_entries) <= MAX_TOTAL_SEQUENCES):
        print(f"  NOTE: final count {len(final_entries)} is outside the target range "
              f"{MIN_TOTAL_SEQUENCES}-{MAX_TOTAL_SEQUENCES} (proceeding anyway; see exclusion log).")

    final_entries = sort_for_output(final_entries)
    write_fasta(final_entries, EXPANDED_FASTA_OUT)
    print(f"  Wrote {len(final_entries)} sequences: {EXPANDED_FASTA_OUT}")
    taxa_df = write_taxa_metadata(final_entries, TAXA_METADATA_OUT)
    print(f"  Wrote taxa metadata: {TAXA_METADATA_OUT}")
    print(taxa_df["taxonomic_group"].value_counts().to_string())

    # --- Step 5: MAFFT alignment --------------------------------------------
    print("\nStep 5: aligning with MAFFT (--auto)...")
    records = run_mafft(EXPANDED_FASTA_OUT, ALIGNED_FASTA_OUT)
    print(f"  Alignment: {len(records)} sequences x {len(records[0].seq)} columns")
    print(f"  Wrote alignment: {ALIGNED_FASTA_OUT}")

    # --- Steps 6-7: conservation tables -------------------------------------
    print("\nSteps 6-7: computing per-column conservation statistics...")
    site_df, full_df, n_sequences = compute_conservation_tables(records)
    site_df.to_csv(SITE_TABLE_OUT, index=False)
    print(f"  Wrote site conservation table: {SITE_TABLE_OUT}")
    print(site_df.to_string(index=False))
    full_df.to_csv(FULL_PROFILE_OUT, index=False)
    print(f"  Wrote full conservation profile ({len(full_df)} positions): {FULL_PROFILE_OUT}")

    # --- Step 8: figure -------------------------------------------------------
    print("\nStep 8: rendering conservation profile figure...")
    render_figure(site_df, full_df, n_sequences, FIGURE_OUT)
    print(f"  Wrote figure: {FIGURE_OUT}")

    # --- Step 9: summary --------------------------------------------------
    print("\nStep 9: writing conservation summary...")
    write_summary(site_df, full_df, taxa_df, n_sequences, len(records[0].seq), SUMMARY_OUT)
    print(f"  Wrote summary: {SUMMARY_OUT}")
    print("\n" + SUMMARY_OUT.read_text())

    sys.exit(0)


if __name__ == "__main__":
    main()
