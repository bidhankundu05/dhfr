"""
Module 5: Structure retrieval + mutant ligand-contact analysis.

Fetches WT human DHFR (1U72) plus the four resistance-mutant structures
(1DLR, 1DLS, 3EIG, 3F8Z) from RCSB PDB, then measures every one of the four
resistance positions (22, 31, 35, 64 -- see config.RESISTANCE_POSITIONS)
against its structure's own co-crystallized ligand: minimum heavy-atom
distance and a simple contact-type classification.

This is the structural leg of the hypothesis: does a mutation at a given
position actually sit close enough to the ligand to plausibly disrupt
binding, and does mutating it shift the contact geometry at OTHER
resistance positions too (allosteric effects)?

Note: not every structure is crystallized with literal methotrexate --
1DLR has the close analog MXA, and 3F8Z has an unrelated inhibitor (DH1 /
PY957 in the source paper), not MTX. This is reported explicitly rather
than assumed uniform; Module 6 handles this properly by re-docking the
actual methotrexate molecule into every receptor for a fair comparison.
"""

import sys
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import config  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import requests  # noqa: E402
from Bio.PDB import PDBParser  # noqa: E402
from Bio.PDB.PDBExceptions import PDBConstructionWarning  # noqa: E402

warnings.simplefilter("ignore", PDBConstructionWarning)

RCSB_DOWNLOAD_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"
STRUCTURE_IDS = [config.REFERENCE_PDB] + list(config.RESISTANCE_MUTANT_PDBS.keys())

CONTACT_TABLE_OUT = config.RESULTS_STRUCTURES_DIR / "mutant_ligand_contacts.csv"
CONTACT_IMAGE_OUT = config.RESULTS_STRUCTURES_DIR / "mutant_ligand_distances.png"
STRUCTURE_NOTES_OUT = config.RESULTS_STRUCTURES_DIR / "structure_notes.txt"

# Ligand HETATM residue code actually co-crystallized in each structure,
# verified directly against each downloaded PDB file's own HETATM records.
STRUCTURE_LIGANDS = {
    "1U72": ("MTX", "methotrexate"),
    "1DLR": ("MXA", "6-(2,5-dimethoxybenzyl)-5-methylpyrido[2,3-d]pyrimidine-2,4-diamine (methotrexate analog)"),
    "1DLS": ("MTX", "methotrexate"),
    "3EIG": ("MTX", "methotrexate"),
    "3F8Z": ("DH1", "PY957 (non-methotrexate inhibitor)"),
}

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}

# WT residue at each of the four resistance positions (mature-protein
# numbering), used as the reference identity for the WT structure row.
WT_RESIDUE_AT_POSITION = {22: "L", 31: "F", 35: "Q", 64: "N"}


def download_structure(pdb_id):
    out_path = config.STRUCTURES_DIR / f"{pdb_id}.pdb"
    if not out_path.exists():
        resp = requests.get(RCSB_DOWNLOAD_URL.format(pdb_id=pdb_id), timeout=30)
        resp.raise_for_status()
        out_path.write_text(resp.text)
    return out_path


def find_chain_with_residue(model, resnum):
    for chain in model:
        for residue in chain:
            if residue.id[1] == resnum and residue.id[0] == " ":
                return chain
    return None


def get_residue(chain, resnum):
    for residue in chain:
        if residue.id[1] == resnum and residue.id[0] == " ":
            return residue
    return None


def find_ligand_residue(model, ligand_code):
    for chain in model:
        for residue in chain:
            if residue.resname == ligand_code:
                return residue
    return None


def min_distance(residue, ligand_atoms):
    if not ligand_atoms:
        return None
    min_dist = None
    for atom in residue:
        if atom.element == "H":
            continue
        for lig_atom in ligand_atoms:
            d = atom - lig_atom
            if min_dist is None or d < min_dist:
                min_dist = d
    return min_dist


def classify_contact(dist):
    if dist is None:
        return "ligand not found"
    if dist < 4.0:
        return "direct contact (<4.0 A)"
    if dist < 5.0:
        return "van der Waals contact (4.0-5.0 A)"
    return "no direct contact (>5.0 A)"


def analyze_structure(pdb_id):
    pdb_path = download_structure(pdb_id)
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_id, pdb_path)
    model = structure[0]

    ligand_code, ligand_name = STRUCTURE_LIGANDS[pdb_id]
    ligand_residue = find_ligand_residue(model, ligand_code)
    ligand_atoms = [a for a in ligand_residue if a.element != "H"] if ligand_residue else []

    mutated_positions = {
        pos: mut for pos, wt, mut in config.RESISTANCE_MUTANT_PDBS.get(pdb_id, {}).get("mutations", [])
    }

    rows = []
    for pos in config.RESISTANCE_POSITIONS:
        chain = find_chain_with_residue(model, pos)
        residue = get_residue(chain, pos) if chain else None
        observed = THREE_TO_ONE.get(residue.resname, "?") if residue else None
        dist = min_distance(residue, ligand_atoms) if residue else None

        rows.append({
            "pdb_id": pdb_id,
            "position": pos,
            "expected_wt_residue": WT_RESIDUE_AT_POSITION[pos],
            "observed_residue": observed,
            "is_mutated_at_this_position": pos in mutated_positions,
            "ligand_code": ligand_code,
            "ligand_name": ligand_name,
            "min_distance_to_ligand_angstrom": round(dist, 2) if dist is not None else None,
            "contact_type": classify_contact(dist),
        })
    return rows, ligand_residue is not None


def render_distance_chart(df):
    positions = config.RESISTANCE_POSITIONS
    structures = STRUCTURE_IDS
    x = np.arange(len(positions))
    width = 0.15

    fig, ax = plt.subplots(figsize=(10, 5.5))
    for i, pdb_id in enumerate(structures):
        sub = df[df["pdb_id"] == pdb_id].set_index("position").reindex(positions)
        heights = sub["min_distance_to_ligand_angstrom"].fillna(0).values
        mutated_flags = sub["is_mutated_at_this_position"].fillna(False).values
        colors = ["#d62728" if m else "#4c72b0" for m in mutated_flags]
        bars = ax.bar(x + i * width, heights, width, label=pdb_id, color=colors, alpha=0.5 + 0.1 * i)
        for bar, h, m in zip(bars, heights, mutated_flags):
            if h > 0:
                ax.text(bar.get_x() + bar.get_width() / 2, h + 0.1, f"{h:.1f}",
                        ha="center", va="bottom", fontsize=6, rotation=90)

    ax.set_xticks(x + width * (len(structures) - 1) / 2)
    ax.set_xticklabels([f"Position {p}\n(WT={WT_RESIDUE_AT_POSITION[p]})" for p in positions])
    ax.set_ylabel("Min. heavy-atom distance to bound ligand (angstrom)")
    ax.set_title("Resistance-position distance to co-crystallized ligand, WT vs mutants")
    ax.axhline(4.0, color="gray", linestyle="--", linewidth=0.8, alpha=0.6)
    ax.axhline(5.0, color="gray", linestyle=":", linewidth=0.8, alpha=0.6)

    from matplotlib.patches import Patch
    color_legend = [
        Patch(facecolor="#d62728", alpha=0.7, label="mutated at this position in this structure"),
        Patch(facecolor="#4c72b0", alpha=0.7, label="WT residue at this position"),
    ]
    ax.legend(handles=color_legend, loc="upper right", fontsize=8)
    fig.text(0.99, 0.02, "Bars left-to-right per group: " + ", ".join(structures),
              ha="right", fontsize=7, color="gray")

    fig.tight_layout()
    fig.savefig(CONTACT_IMAGE_OUT, dpi=200, bbox_inches="tight")
    plt.close(fig)


def write_structure_notes(df, ligand_found):
    lines = ["Module 5 -- structure retrieval and ligand-contact analysis", "=" * 70, ""]
    lines.append("Structures fetched:")
    for pdb_id in STRUCTURE_IDS:
        ligand_code, ligand_name = STRUCTURE_LIGANDS[pdb_id]
        found = "OK" if ligand_found.get(pdb_id) else "LIGAND NOT FOUND IN STRUCTURE"
        mutations = config.RESISTANCE_MUTANT_PDBS.get(pdb_id, {}).get("mutations", [])
        mut_desc = ", ".join(f"{wt}{pos}{mut}" for pos, wt, mut in mutations) if mutations else "wild-type"
        lines.append(f"  {pdb_id}: {mut_desc} -- ligand {ligand_code} ({ligand_name}) [{found}]")
    lines.append("")
    lines.append(
        "IMPORTANT: not every structure is co-crystallized with literal methotrexate. "
        "1U72, 1DLS, and 3EIG have MTX bound directly; 1DLR has the close pyridopyrimidine "
        "analog MXA; 3F8Z has an unrelated inhibitor (DH1/PY957). Distances below are each "
        "measured against whatever ligand is actually present in that specific crystal -- "
        "Module 6 re-docks the real methotrexate molecule into every receptor uniformly so "
        "the ΔΔG comparison there is apples-to-apples."
    )
    lines.append("")
    lines.append("Per-position contact summary (mutated position highlighted with *):")
    lines.append("-" * 70)
    for pdb_id in STRUCTURE_IDS:
        lines.append(f"\n{pdb_id}:")
        sub = df[df["pdb_id"] == pdb_id]
        for _, row in sub.iterrows():
            marker = "*" if row["is_mutated_at_this_position"] else " "
            lines.append(
                f"  {marker} position {row['position']} "
                f"(WT={row['expected_wt_residue']}, observed={row['observed_residue']}): "
                f"{row['min_distance_to_ligand_angstrom']} A -- {row['contact_type']}"
            )
    STRUCTURE_NOTES_OUT.write_text("\n".join(lines))


def main():
    print("=" * 90)
    print("Module 5: Structure retrieval + mutant ligand-contact analysis")
    print("=" * 90)

    config.STRUCTURES_DIR.mkdir(parents=True, exist_ok=True)
    config.RESULTS_STRUCTURES_DIR.mkdir(parents=True, exist_ok=True)

    all_rows = []
    ligand_found = {}
    for pdb_id in STRUCTURE_IDS:
        print(f"\n[{pdb_id}]")
        rows, found = analyze_structure(pdb_id)
        ligand_found[pdb_id] = found
        ligand_code, ligand_name = STRUCTURE_LIGANDS[pdb_id]
        print(f"    Downloaded, ligand {ligand_code} ({ligand_name}): {'found' if found else 'NOT FOUND'}")
        for row in rows:
            marker = "*" if row["is_mutated_at_this_position"] else " "
            print(f"    {marker} pos {row['position']}: {row['observed_residue']} "
                  f"-> {row['min_distance_to_ligand_angstrom']} A ({row['contact_type']})")
        all_rows.extend(rows)

    df = pd.DataFrame(all_rows)
    df.to_csv(CONTACT_TABLE_OUT, index=False)
    print(f"\nWrote contact table: {CONTACT_TABLE_OUT}")

    render_distance_chart(df)
    print(f"Wrote distance chart: {CONTACT_IMAGE_OUT}")

    write_structure_notes(df, ligand_found)
    print(f"Wrote structure notes: {STRUCTURE_NOTES_OUT}")

    sys.exit(0)


if __name__ == "__main__":
    main()
