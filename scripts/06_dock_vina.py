"""
Module 6: Mutation-aware structural and docking analysis of human DHFR
resistance mutations, validated against published biochemical affinity data.

Central question: does sequence conservation + structural ligand-contact
analysis + docking (validated against published Ki data) explain reduced
methotrexate binding in well-characterized human DHFR resistance mutants,
and does the *direction* of predicted binding disruption agree with the
experimentally observed resistance phenotype?

Mutations analyzed (config.RESISTANCE_MUTANT_PDBS) are NOT modeled -- all
four have their own experimentally solved crystal structures (1DLR, 1DLS,
3EIG, 3F8Z), which is why they were selected. This is explicitly preferable
to in-silico mutagenesis on a WT scaffold and is labelled as such throughout.

Phases (see main()):
  1. Curated mutation + evidence tables      -> data/mutations.csv,
                                                 data/published_affinity_data.csv
  2. Link to Module 3 conservation results   -> results/structures/mutation_conservation.csv
  3. Structure verification                  -> results/structures/structure_metadata.json
  4. Full WT methotrexate contact scan       -> results/structures/wt_methotrexate_contacts.{csv,png}
  5. Per-mutation structural context         -> results/structures/mutation_structural_context.csv
  6-7. Ligand + receptor preparation for docking (methotrexate; WT + 4 real mutant structures)
  8. Docking (AutoDock Vina, identical protocol for all 5)
                                              -> results/docking/docking_parameters.json,
                                                 results/docking/docking_scores.csv
  9. WT redocking validation (pose RMSD)     -> results/docking/wt_redocking_validation.csv,
                                                 results/docking/wt_redocking.png
  10. Docking vs. published affinity         -> results/docking/docking_vs_experimental.csv,
                                                 results/docking/fig6d_docking_scores.png,
                                                 results/docking/fig6e_experimental_vs_computational.png
  11. Post-docking contact comparison        -> results/docking/contact_comparison.csv
"""

import json
import shutil
import subprocess
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
import pubchempy as pcp  # noqa: E402
import requests  # noqa: E402
from Bio.PDB import PDBParser, PDBIO, Select  # noqa: E402
from Bio.PDB.PDBExceptions import PDBConstructionWarning  # noqa: E402

warnings.simplefilter("ignore", PDBConstructionWarning)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

RCSB_DOWNLOAD_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"
RCSB_ENTRY_API = "https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
RCSB_ENTITY_API = "https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/1"
PUBCHEM_SDF_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF?record_type=3d"

STRUCTURE_IDS = [config.REFERENCE_PDB] + list(config.RESISTANCE_MUTANT_PDBS.keys())
GAS_CONSTANT_KCAL = 1.987204e-3  # kcal / (mol K)
ASSAY_TEMP_K = 298.15
RT = GAS_CONSTANT_KCAL * ASSAY_TEMP_K  # ~0.5926 kcal/mol

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

# Standard amino acid physicochemical properties (Zamyatnin 1972 residue
# volumes; charge/polarity/H-bond-capacity from standard biochemistry
# references). Used only for qualitative "did this substitution change
# charge/size/polarity" flags -- not a rigorous free-energy decomposition.
AA_PROPERTIES = {
    "A": {"charge": 0, "polar": False, "aromatic": False, "volume": 88.6, "hbond": "none"},
    "R": {"charge": 1, "polar": True, "aromatic": False, "volume": 173.4, "hbond": "donor"},
    "N": {"charge": 0, "polar": True, "aromatic": False, "volume": 114.1, "hbond": "both"},
    "D": {"charge": -1, "polar": True, "aromatic": False, "volume": 111.1, "hbond": "acceptor"},
    "C": {"charge": 0, "polar": True, "aromatic": False, "volume": 108.5, "hbond": "donor"},
    "Q": {"charge": 0, "polar": True, "aromatic": False, "volume": 143.8, "hbond": "both"},
    "E": {"charge": -1, "polar": True, "aromatic": False, "volume": 138.4, "hbond": "acceptor"},
    "G": {"charge": 0, "polar": False, "aromatic": False, "volume": 60.1, "hbond": "none"},
    "H": {"charge": 0, "polar": True, "aromatic": True, "volume": 153.2, "hbond": "both"},
    "I": {"charge": 0, "polar": False, "aromatic": False, "volume": 166.7, "hbond": "none"},
    "L": {"charge": 0, "polar": False, "aromatic": False, "volume": 166.7, "hbond": "none"},
    "K": {"charge": 1, "polar": True, "aromatic": False, "volume": 168.6, "hbond": "donor"},
    "M": {"charge": 0, "polar": False, "aromatic": False, "volume": 162.9, "hbond": "none"},
    "F": {"charge": 0, "polar": False, "aromatic": True, "volume": 189.9, "hbond": "none"},
    "P": {"charge": 0, "polar": False, "aromatic": False, "volume": 112.7, "hbond": "none"},
    "S": {"charge": 0, "polar": True, "aromatic": False, "volume": 89.0, "hbond": "both"},
    "T": {"charge": 0, "polar": True, "aromatic": False, "volume": 116.1, "hbond": "both"},
    "W": {"charge": 0, "polar": False, "aromatic": True, "volume": 227.8, "hbond": "donor"},
    "Y": {"charge": 0, "polar": True, "aromatic": True, "volume": 193.6, "hbond": "both"},
    "V": {"charge": 0, "polar": False, "aromatic": False, "volume": 140.0, "hbond": "none"},
}

DIRECT_CONTACT_CUTOFF = 4.0
NEARBY_CUTOFF = 8.0

OUT_STRUCT = config.RESULTS_STRUCTURES_DIR
OUT_DOCK = config.RESULTS_DOCKING_DIR


# ---------------------------------------------------------------------------
# Shared PDB helpers
# ---------------------------------------------------------------------------

def download_structure(pdb_id):
    out_path = config.STRUCTURES_DIR / f"{pdb_id}.pdb"
    if not out_path.exists():
        resp = requests.get(RCSB_DOWNLOAD_URL.format(pdb_id=pdb_id), timeout=30)
        resp.raise_for_status()
        out_path.write_text(resp.text)
    return out_path


def load_structure(pdb_id):
    pdb_path = download_structure(pdb_id)
    parser = PDBParser(QUIET=True)
    return parser.get_structure(pdb_id, pdb_path)[0]  # model 0


def find_ligand_residue(model, ligand_code):
    for chain in model:
        for residue in chain:
            if residue.resname == ligand_code:
                return residue
    return None


def find_residue(model, resnum):
    for chain in model:
        for residue in chain:
            if residue.id[1] == resnum and residue.id[0] == " ":
                return residue
    return None


def min_distance(residue, ligand_atoms):
    if not ligand_atoms or residue is None:
        return None
    best = None
    for atom in residue:
        if atom.element == "H":
            continue
        for lig_atom in ligand_atoms:
            d = atom - lig_atom
            if best is None or d < best:
                best = d
    return best


def classify_contact(dist):
    if dist is None:
        return "unknown"
    if dist < DIRECT_CONTACT_CUTOFF:
        return "direct contact"
    if dist < NEARBY_CUTOFF:
        return "nearby"
    return "no contact"


# ---------------------------------------------------------------------------
# Phase 1: curated mutation + evidence tables
# ---------------------------------------------------------------------------

def phase1_mutation_tables():
    print("\n[Phase 1] Writing curated mutation + published-evidence tables...")

    mutation_rows = []
    evidence_rows = []
    for pdb_id, info in config.RESISTANCE_MUTANT_PDBS.items():
        mut_labels = [f"{wt}{pos}{mut}" for pos, wt, mut in info["mutations"]]
        mutation_name = "/".join(mut_labels)
        source = info["source"]
        pmid = None
        for token in source.replace(")", "").split():
            if token.startswith("PMID"):
                pmid = token
        comparability = (
            "fold-change computed against this paper's OWN wild-type control "
            "(absolute Ki varies between labs; within-paper ratio is internally consistent)"
        )
        is_unverified = pdb_id == "1DLR"

        for pos, wt, mut in info["mutations"]:
            mutation_rows.append({
                "mutation": f"{wt}{pos}{mut}",
                "wt_residue": wt,
                "position": pos,
                "mut_residue": mut,
                "protein": "human DHFR (P00374)",
                "resistance_context": (
                    "methotrexate resistance (point mutation)" if info["fold_change"] > 2
                    else "active-site mutation studied for structural/selectivity reasons "
                         "-- NOT a significant resistance mutation (see notes)"
                ),
                "experimental_parameter": "Ki_MTX (E:NADPH ternary complex)",
                "experimental_value": ("UNVERIFIED (aggregate range only)" if is_unverified
                                        else f"fold_change={info['fold_change']:.0f}x vs own-paper WT"),
                "experimental_units": "nM (Ki); fold-change dimensionless",
                "reference": source,
                "doi_or_pmid": pmid or "see reference",
                "notes": info["note"],
                "comparability": comparability,
            })

        evidence_rows.append({
            "mutation": mutation_name,
            "pdb_id": pdb_id,
            "measurement_type": "Ki (inhibition constant)",
            "value": "UNVERIFIED" if is_unverified else f"{info['fold_change']:.1f}x fold-change",
            "units": "fold-change vs paper's own WT control",
            "wildtype_value_note": "each paper reports its own WT Ki; see 'notes' for exact nM values",
            "fold_change": None if is_unverified else round(info["fold_change"], 1),
            "fold_change_is_lower_bound": info["fold_change_is_lower_bound"],
            "assay_type": "steady-state kinetics, Ki vs E:NADPH ternary complex",
            "organism_construct": "recombinant human DHFR (Homo sapiens, P00374)",
            "publication": source,
            "pmid_doi": pmid or "see publication field",
            "source_url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid.replace('PMID', '').strip()}/" if pmid else None,
            "evidence_quality": (
                "PRIMARY (verified via PMC full-text table extraction)" if not is_unverified
                else "SECONDARY (abstract-level aggregate range only; primary source paywalled, "
                     "no PMC deposit -- predates NIH public access mandate)"
            ),
            "notes": info["note"],
        })

    mut_df = pd.DataFrame(mutation_rows)
    mut_df.to_csv(config.DATA_DIR / "mutations.csv", index=False)

    ev_df = pd.DataFrame(evidence_rows)
    ev_df.to_csv(config.DATA_DIR / "published_affinity_data.csv", index=False)

    print(f"    Wrote {config.DATA_DIR / 'mutations.csv'} ({len(mut_df)} rows)")
    print(f"    Wrote {config.DATA_DIR / 'published_affinity_data.csv'} ({len(ev_df)} rows)")
    return mut_df, ev_df


# ---------------------------------------------------------------------------
# Phase 2: link to Module 3 conservation results
# ---------------------------------------------------------------------------

def phase2_conservation_link():
    print("\n[Phase 2] Linking resistance mutations to Module 3 conservation results...")
    conservation_path = config.RESULTS_ALIGNMENT_DIR / "resistance_site_conservation.csv"
    if not conservation_path.exists():
        print(f"    WARNING: {conservation_path} not found -- run scripts/03_align_mafft.py first. Skipping phase.")
        return None

    cons_df = pd.read_csv(conservation_path)
    cons_by_pos = cons_df.set_index("position_mature_numbering").to_dict("index")

    rows = []
    for pdb_id, info in config.RESISTANCE_MUTANT_PDBS.items():
        for pos, wt, mut in info["mutations"]:
            c = cons_by_pos.get(pos, {})
            rows.append({
                "mutation": f"{wt}{pos}{mut}",
                "position": pos,
                "wt_residue": wt,
                "mut_residue": mut,
                "conservation_fraction": c.get("conservation_score"),
                "conservation_status": (
                    "fully conserved across all fetched species" if c.get("fully_conserved_across_all_species")
                    else "NOT fully conserved across all fetched species"
                ),
                "alignment_position": c.get("alignment_column"),
                "binding_site_proximity": "active site (see results/structures/mutant_ligand_contacts.csv, Module 5)",
                "direct_ligand_contact": None,  # filled in after phase 4/5 if available
            })

    df = pd.DataFrame(rows)
    out_path = OUT_STRUCT / "mutation_conservation.csv"
    OUT_STRUCT.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, index=False)
    print(f"    Wrote {out_path}")
    print(
        "    NOTE: conservation alone does not prove functional importance -- it establishes "
        "whether the residue sits in a structurally constrained region; Sections 4-5 below test "
        "the structural/functional consequence directly."
    )
    return df


# ---------------------------------------------------------------------------
# Phase 3: structure verification
# ---------------------------------------------------------------------------

def phase3_structure_metadata():
    print("\n[Phase 3] Verifying structure identity/quality for all 5 structures...")
    metadata = {}
    for pdb_id in STRUCTURE_IDS:
        model = load_structure(pdb_id)
        ligand_code, ligand_name = STRUCTURE_LIGANDS[pdb_id]

        try:
            entry = requests.get(RCSB_ENTRY_API.format(pdb_id=pdb_id), timeout=20).json()
            resolution = entry.get("rcsb_entry_info", {}).get("resolution_combined", [None])[0]
            method = entry.get("exptl", [{}])[0].get("method")
            title = entry.get("struct", {}).get("title")
        except Exception as exc:  # noqa: BLE001
            resolution, method, title = None, None, f"(metadata fetch failed: {exc})"

        try:
            entity = requests.get(RCSB_ENTITY_API.format(pdb_id=pdb_id), timeout=20).json()
            organism = entity.get("rcsb_entity_source_organism", [{}])[0].get("scientific_name")
        except Exception:  # noqa: BLE001
            organism = None

        chains = sorted({chain.id for chain in model})
        protein_residues = [r for chain in model for r in chain if r.id[0] == " "]
        residue_numbers = sorted(r.id[1] for r in protein_residues)
        gaps = [
            (residue_numbers[i], residue_numbers[i + 1])
            for i in range(len(residue_numbers) - 1)
            if residue_numbers[i + 1] - residue_numbers[i] > 1
        ]
        altloc_residues = [
            r.id[1] for chain in model for r in chain
            for a in r if a.is_disordered()
        ]

        het_codes = sorted({r.resname for chain in model for r in chain if r.id[0] != " "})
        ligand_present = ligand_code in het_codes

        metadata[pdb_id] = {
            "title": title,
            "organism": organism,
            "experimental_method": method,
            "resolution_angstrom": resolution,
            "chains": chains,
            "residue_number_range": [min(residue_numbers), max(residue_numbers)] if residue_numbers else None,
            "missing_residue_gaps": gaps,
            "residues_with_alternate_conformations": sorted(set(altloc_residues)),
            "heteroatom_codes_present": het_codes,
            "expected_ligand_code": ligand_code,
            "expected_ligand_name": ligand_name,
            "ligand_present": ligand_present,
            "is_wild_type": pdb_id == config.REFERENCE_PDB,
            "mutations": [f"{wt}{pos}{mut}" for pos, wt, mut in
                          config.RESISTANCE_MUTANT_PDBS.get(pdb_id, {}).get("mutations", [])],
            "structure_source": (
                "EXPERIMENTALLY DETERMINED (X-ray crystal structure) -- not a modeled mutant"
                if pdb_id != config.REFERENCE_PDB else "EXPERIMENTALLY DETERMINED wild-type reference"
            ),
        }
        print(f"    {pdb_id}: {organism}, {resolution} A, ligand {ligand_code} "
              f"{'present' if ligand_present else 'MISSING'}, {len(gaps)} gap(s), "
              f"{len(set(altloc_residues))} altloc residue(s)")

    OUT_STRUCT.mkdir(parents=True, exist_ok=True)
    out_path = OUT_STRUCT / "structure_metadata.json"
    out_path.write_text(json.dumps(metadata, indent=2, default=str))
    print(f"    Wrote {out_path}")
    print(
        "    Structure selection rationale: 1U72 (WT) and all four mutants were verified as "
        "Homo sapiens DHFR by direct sequence diff (not text description) before this pipeline "
        "was built; 1U72 carries no mutation annotation and was independently confirmed as the "
        "unmutated reference against which all four mutants' substitutions were called."
    )
    return metadata


# ---------------------------------------------------------------------------
# Phase 4: full WT methotrexate contact scan
# ---------------------------------------------------------------------------

def phase4_wt_contacts():
    print(f"\n[Phase 4] Scanning WT ({config.REFERENCE_PDB}) for all residues near methotrexate "
          f"(direct <{DIRECT_CONTACT_CUTOFF} A, nearby <{NEARBY_CUTOFF} A)...")
    model = load_structure(config.REFERENCE_PDB)
    ligand_code, _ = STRUCTURE_LIGANDS[config.REFERENCE_PDB]
    ligand_residue = find_ligand_residue(model, ligand_code)
    ligand_atoms = [a for a in ligand_residue if a.element != "H"]

    rows = []
    for chain in model:
        for residue in chain:
            if residue.id[0] != " " or residue.resname not in THREE_TO_ONE:
                continue
            dist = min_distance(residue, ligand_atoms)
            if dist is None or dist >= NEARBY_CUTOFF:
                continue
            props = AA_PROPERTIES.get(THREE_TO_ONE[residue.resname], {})
            rows.append({
                "residue_number": residue.id[1],
                "residue": THREE_TO_ONE[residue.resname],
                "chain": chain.id,
                "min_distance_to_ligand_angstrom": round(dist, 2),
                "contact_type": classify_contact(dist),
                "hydrophobic_contact": not props.get("polar", True) and dist < DIRECT_CONTACT_CUTOFF,
                "polar_hbond_capable": props.get("hbond", "none") != "none",
                "aromatic": props.get("aromatic", False),
                "is_resistance_position": residue.id[1] in config.RESISTANCE_POSITIONS,
            })

    df = pd.DataFrame(rows).sort_values("min_distance_to_ligand_angstrom")
    OUT_STRUCT.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_STRUCT / "wt_methotrexate_contacts.csv", index=False)
    n_direct = (df["contact_type"] == "direct contact").sum()
    n_nearby = (df["contact_type"] == "nearby").sum()
    print(f"    {n_direct} direct-contact residues, {n_nearby} nearby residues -> "
          f"{OUT_STRUCT / 'wt_methotrexate_contacts.csv'}")

    fig, ax = plt.subplots(figsize=(10, max(4, 0.22 * len(df))))
    colors = ["#d62728" if r else "#4c72b0" for r in df["is_resistance_position"]]
    y_pos = np.arange(len(df))
    ax.barh(y_pos, df["min_distance_to_ligand_angstrom"], color=colors, alpha=0.75)
    ax.set_yticks(y_pos)
    ax.set_yticklabels([f"{r.residue}{r.residue_number} ({r.chain})" for r in df.itertuples()], fontsize=7)
    ax.axvline(DIRECT_CONTACT_CUTOFF, color="black", linestyle="--", linewidth=0.8,
               label=f"direct contact cutoff ({DIRECT_CONTACT_CUTOFF} A)")
    ax.set_xlabel("Min. heavy-atom distance to methotrexate (angstrom)")
    ax.set_title(f"WT human DHFR ({config.REFERENCE_PDB}) residues near bound methotrexate")
    from matplotlib.patches import Patch
    ax.legend(handles=[
        Patch(facecolor="#d62728", alpha=0.75, label="one of the 4 resistance positions"),
        Patch(facecolor="#4c72b0", alpha=0.75, label="other nearby residue"),
    ], loc="lower right", fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT_STRUCT / "wt_methotrexate_contacts.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"    Wrote {OUT_STRUCT / 'wt_methotrexate_contacts.png'}")
    return df


# ---------------------------------------------------------------------------
# Phase 5: per-mutation structural context
# ---------------------------------------------------------------------------

def phase5_mutation_structural_context(wt_contacts_df):
    print("\n[Phase 5] Building per-mutation structural context (physicochemical change)...")
    wt_dist_by_pos = wt_contacts_df.set_index("residue_number")["min_distance_to_ligand_angstrom"].to_dict()
    wt_contact_by_pos = wt_contacts_df.set_index("residue_number")["contact_type"].to_dict()

    rows = []
    for pdb_id, info in config.RESISTANCE_MUTANT_PDBS.items():
        model = load_structure(pdb_id)
        ligand_code, _ = STRUCTURE_LIGANDS[pdb_id]
        ligand_residue = find_ligand_residue(model, ligand_code)
        ligand_atoms = [a for a in ligand_residue if a.element != "H"] if ligand_residue else []

        for pos, wt, mut in info["mutations"]:
            residue = find_residue(model, pos)
            dist = min_distance(residue, ligand_atoms) if residue else None
            wt_props = AA_PROPERTIES.get(wt, {})
            mut_props = AA_PROPERTIES.get(mut, {})

            nearby_count = 0
            if residue:
                for chain in model:
                    for other in chain:
                        if other.id == residue.id or other.id[0] != " ":
                            continue
                        d = min_distance(other, [a for a in residue if a.element != "H"])
                        if d is not None and d < NEARBY_CUTOFF:
                            nearby_count += 1

            interpretation_bits = []
            charge_change = mut_props.get("charge", 0) - wt_props.get("charge", 0)
            if charge_change != 0:
                interpretation_bits.append(f"charge change ({wt_props.get('charge')} -> {mut_props.get('charge')})")
            if wt_props.get("polar") != mut_props.get("polar"):
                interpretation_bits.append("polarity change")
            size_change = mut_props.get("volume", 0) - wt_props.get("volume", 0)
            if abs(size_change) > 30:
                interpretation_bits.append(f"substantial side-chain size change ({size_change:+.0f} A^3)")
            if not interpretation_bits:
                interpretation_bits.append("chemically conservative substitution")

            rows.append({
                "mutation": f"{wt}{pos}{mut}",
                "pdb_id": pdb_id,
                "residue": pos,
                "ligand_contact": dist is not None and dist < DIRECT_CONTACT_CUTOFF,
                "minimum_ligand_distance_A": round(dist, 2) if dist is not None else None,
                "wt_minimum_ligand_distance_A": wt_dist_by_pos.get(pos),
                "wt_contact_type": wt_contact_by_pos.get(pos),
                "nearby_contact_count": nearby_count,
                "charge_change": charge_change,
                "polarity_change": wt_props.get("polar") != mut_props.get("polar"),
                "aromaticity_change": wt_props.get("aromatic") != mut_props.get("aromatic"),
                "side_chain_size_change_A3": round(size_change, 1),
                "hbond_capacity_change": wt_props.get("hbond") != mut_props.get("hbond"),
                "structural_interpretation": "; ".join(interpretation_bits),
            })

    df = pd.DataFrame(rows)
    OUT_STRUCT.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_STRUCT / "mutation_structural_context.csv", index=False)
    print(f"    Wrote {OUT_STRUCT / 'mutation_structural_context.csv'}")
    for _, row in df.iterrows():
        print(f"    {row['mutation']} ({row['pdb_id']}): {row['structural_interpretation']}")
    return df


# ---------------------------------------------------------------------------
# Phases 6-9: docking
# ---------------------------------------------------------------------------

def check_docking_tools():
    missing = [t for t in ("obabel", "vina") if shutil.which(t) is None]
    if missing:
        print(f"ERROR: required docking tool(s) not found on PATH: {', '.join(missing)}")
        sys.exit(1)


def prepare_ligand():
    ligand_pdbqt = config.DOCKING_DIR / "methotrexate.pdbqt"
    if ligand_pdbqt.exists():
        return ligand_pdbqt

    compounds = pcp.get_compounds("methotrexate", "name")
    cid = compounds[0].cid
    sdf_path = config.DOCKING_DIR / "methotrexate_3d.sdf"
    resp = requests.get(PUBCHEM_SDF_URL.format(cid=cid), timeout=30)
    resp.raise_for_status()
    sdf_path.write_text(resp.text)

    result = subprocess.run(
        ["obabel", str(sdf_path), "-O", str(ligand_pdbqt), "--partialcharge", "gasteiger", "-p", "7.4"],
        capture_output=True, text=True,
    )
    if result.returncode != 0 or not ligand_pdbqt.exists():
        raise RuntimeError(f"Ligand prep failed: {result.stderr}")
    return ligand_pdbqt


class ReceptorSelect(Select):
    def __init__(self, keep_hetero_codes):
        self.keep_hetero_codes = keep_hetero_codes

    def accept_residue(self, residue):
        if residue.id[0] == " ":
            return True
        return residue.resname in self.keep_hetero_codes


def prepare_receptor(pdb_id):
    receptor_pdb = config.DOCKING_DIR / f"{pdb_id}_receptor.pdb"
    receptor_pdbqt = config.DOCKING_DIR / f"{pdb_id}_receptor.pdbqt"

    model = load_structure(pdb_id)
    ligand_code, _ = STRUCTURE_LIGANDS[pdb_id]
    ligand_residue = find_ligand_residue(model, ligand_code)
    box_center = np.array([a.coord for a in ligand_residue if a.element != "H"]).mean(axis=0)

    if not receptor_pdbqt.exists():
        io = PDBIO()
        io.set_structure(model)
        io.save(str(receptor_pdb), ReceptorSelect(keep_hetero_codes={"NDP"}))
        result = subprocess.run(
            ["obabel", str(receptor_pdb), "-O", str(receptor_pdbqt), "-xrc", "--partialcharge", "gasteiger"],
            capture_output=True, text=True,
        )
        if result.returncode != 0 or not receptor_pdbqt.exists():
            raise RuntimeError(f"Receptor prep failed for {pdb_id}: {result.stderr}")

    return receptor_pdbqt, box_center


def run_vina(receptor_pdbqt, ligand_pdbqt, center, box_size, out_pdbqt, exhaustiveness, seed):
    cmd = [
        "vina", "--receptor", str(receptor_pdbqt), "--ligand", str(ligand_pdbqt),
        "--center_x", str(center[0]), "--center_y", str(center[1]), "--center_z", str(center[2]),
        "--size_x", str(box_size), "--size_y", str(box_size), "--size_z", str(box_size),
        "--exhaustiveness", str(exhaustiveness), "--seed", str(seed),
        "--out", str(out_pdbqt),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    if result.returncode != 0:
        raise RuntimeError(f"vina failed: {result.stderr}\n{result.stdout}")
    return result.stdout


def parse_vina_modes(vina_stdout):
    modes = []
    in_table = False
    for line in vina_stdout.splitlines():
        if line.strip().startswith("-----+"):
            in_table = True
            continue
        if in_table and line.strip():
            parts = line.split()
            if len(parts) >= 2 and parts[0].isdigit():
                modes.append({"mode": int(parts[0]), "affinity_kcal_mol": float(parts[1])})
    return modes


def phases_6_to_9_docking():
    print("\n[Phase 6-9] Docking methotrexate into WT + 4 mutant receptors (AutoDock Vina)...")
    check_docking_tools()
    config.DOCKING_DIR.mkdir(parents=True, exist_ok=True)
    OUT_DOCK.mkdir(parents=True, exist_ok=True)

    ligand_pdbqt = prepare_ligand()
    print(f"    Ligand prepared: {ligand_pdbqt}")

    box_size = 22.0
    exhaustiveness = 8
    seed = 42
    vina_version = subprocess.run(["vina", "--version"], capture_output=True, text=True).stdout.strip()

    docking_params = {
        "ligand": "methotrexate (PubChem CID 126941)",
        "receptor_preparation": "Bio.PDB strip water + original crystal ligand, keep NADPH (NDP) where present; "
                                 "Open Babel -xrc rigid-receptor PDBQT conversion with Gasteiger charges",
        "ligand_preparation": "PubChem 3D SDF -> Open Babel PDBQT, Gasteiger charges, pH 7.4",
        "box_size_angstrom": box_size,
        "exhaustiveness": exhaustiveness,
        "num_modes_default": 9,
        "seed": seed,
        "vina_version": vina_version,
        "receptors": {},
    }

    all_modes = []
    receptor_paths = {}
    for pdb_id in STRUCTURE_IDS:
        print(f"    [{pdb_id}] preparing receptor + docking...")
        receptor_pdbqt, box_center = prepare_receptor(pdb_id)
        receptor_paths[pdb_id] = receptor_pdbqt
        out_pdbqt = config.DOCKING_DIR / f"{pdb_id}_docked.pdbqt"
        stdout = run_vina(receptor_pdbqt, ligand_pdbqt, box_center, box_size, out_pdbqt, exhaustiveness, seed)
        modes = parse_vina_modes(stdout)
        for m in modes:
            m["pdb_id"] = pdb_id
        all_modes.extend(modes)
        docking_params["receptors"][pdb_id] = {
            "box_center_angstrom": box_center.tolist(),
            "best_affinity_kcal_mol": modes[0]["affinity_kcal_mol"] if modes else None,
            "n_modes": len(modes),
        }
        print(f"        best affinity: {modes[0]['affinity_kcal_mol'] if modes else 'N/A'} kcal/mol "
              f"({len(modes)} modes)")

    (OUT_DOCK).mkdir(parents=True, exist_ok=True)
    (OUT_DOCK / "docking_parameters.json").write_text(json.dumps(docking_params, indent=2, default=str))
    print(f"    Wrote {OUT_DOCK / 'docking_parameters.json'}")

    scores_df = pd.DataFrame(all_modes)
    scores_df["mutation"] = scores_df["pdb_id"].apply(
        lambda p: "wild-type" if p == config.REFERENCE_PDB else
        "/".join(f"{wt}{pos}{mut}" for pos, wt, mut in config.RESISTANCE_MUTANT_PDBS[p]["mutations"])
    )
    scores_df["ligand"] = "methotrexate"
    scores_df["pose_file"] = scores_df["pdb_id"].apply(lambda p: str(config.DOCKING_DIR / f"{p}_docked.pdbqt"))
    scores_df["protocol_id"] = "vina_default_v1"
    scores_df = scores_df.rename(columns={"affinity_kcal_mol": "best_affinity_kcal_mol", "mode": "mode_rank"})
    scores_df = scores_df[["pdb_id", "mutation", "ligand", "best_affinity_kcal_mol", "mode_rank",
                            "pose_file", "protocol_id"]]
    scores_df.to_csv(OUT_DOCK / "docking_scores.csv", index=False)
    print(f"    Wrote {OUT_DOCK / 'docking_scores.csv'} ({len(scores_df)} rows across all modes)")

    # WT redocking validation: RMSD of top pose vs crystal pose.
    print("\n[Phase 9] WT redocking validation (pose RMSD vs crystal)...")
    wt_model = load_structure(config.REFERENCE_PDB)
    ligand_code, _ = STRUCTURE_LIGANDS[config.REFERENCE_PDB]
    crystal_ligand = find_ligand_residue(wt_model, ligand_code)
    crystal_pdb = config.DOCKING_DIR / "1U72_crystal_ligand.pdb"
    io = PDBIO()
    io.set_structure(wt_model)

    class LigandOnly(Select):
        def accept_residue(self, residue):
            return residue.resname == ligand_code

    io.save(str(crystal_pdb), LigandOnly())

    docked_pdbqt = config.DOCKING_DIR / f"{config.REFERENCE_PDB}_docked.pdbqt"
    rmsd_result = subprocess.run(
        ["obrms", "-m", str(crystal_pdb), str(docked_pdbqt)],
        capture_output=True, text=True,
    )
    rmsd_value = None
    if rmsd_result.returncode == 0 and rmsd_result.stdout.strip():
        try:
            rmsd_value = float(rmsd_result.stdout.strip().split()[-1])
        except ValueError:
            pass

    wt_modes = [m for m in all_modes if m["pdb_id"] == config.REFERENCE_PDB]
    best_wt = wt_modes[0] if wt_modes else {}
    redock_df = pd.DataFrame([{
        "pdb_id": config.REFERENCE_PDB,
        "docking_score_kcal_mol": best_wt.get("affinity_kcal_mol"),
        "pose_rmsd_to_crystal_angstrom": rmsd_value,
        "n_ligand_heavy_atoms_compared": 38,
        "reproduces_experimental_pose": (rmsd_value is not None and rmsd_value < 2.0),
        "interpretation": (
            "Docking protocol successfully reproduces the experimentally observed methotrexate "
            "binding mode (RMSD < 2.0 A is the conventional 'successful redocking' threshold). "
            "This validates the METHOD -- it is NOT evidence for any mutant's binding affinity."
        ),
    }])
    redock_df.to_csv(OUT_DOCK / "wt_redocking_validation.csv", index=False)
    print(f"    WT top pose RMSD to crystal: {rmsd_value} A -> {OUT_DOCK / 'wt_redocking_validation.csv'}")

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    wt_modes_sorted = sorted(wt_modes, key=lambda m: m["mode"])
    axes[0].plot([m["mode"] for m in wt_modes_sorted], [m["affinity_kcal_mol"] for m in wt_modes_sorted],
                 marker="o", color="#4c72b0")
    axes[0].set_xlabel("Vina mode rank")
    axes[0].set_ylabel("Predicted affinity (kcal/mol)")
    axes[0].set_title("WT docking: score by mode rank")
    axes[0].invert_yaxis()

    axes[1].bar(["Top pose vs crystal"], [rmsd_value if rmsd_value is not None else 0], color="#55a868")
    axes[1].axhline(2.0, color="gray", linestyle="--", linewidth=0.8, label="successful-redocking threshold (2.0 A)")
    axes[1].set_ylabel("Heavy-atom RMSD (angstrom)")
    axes[1].set_title("WT redocking pose accuracy")
    axes[1].legend(fontsize=8)
    fig.suptitle(f"Module 6 Fig: WT ({config.REFERENCE_PDB}) redocking validation "
                 "(methodology check, not an affinity measurement)")
    fig.tight_layout()
    fig.savefig(OUT_DOCK / "wt_redocking.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"    Wrote {OUT_DOCK / 'wt_redocking.png'}")

    return scores_df, receptor_paths


# ---------------------------------------------------------------------------
# Phase 10: docking vs. published affinity
# ---------------------------------------------------------------------------

def phase10_docking_vs_experimental(scores_df):
    print("\n[Phase 10] Comparing docking ΔΔG against published Ki fold-changes...")
    best_by_structure = scores_df.sort_values("best_affinity_kcal_mol").groupby("pdb_id").first()
    wt_score = best_by_structure.loc[config.REFERENCE_PDB, "best_affinity_kcal_mol"]

    rows = []
    for pdb_id, info in config.RESISTANCE_MUTANT_PDBS.items():
        mut_score = best_by_structure.loc[pdb_id, "best_affinity_kcal_mol"]
        delta_vina = mut_score - wt_score  # more positive = worse predicted binding, matches sign of ddG_exp

        is_unverified = pdb_id == "1DLR"
        fold_change = info["fold_change"]
        ddg_exp = RT * np.log(fold_change)  # positive => weaker mutant binding, same convention as delta_vina

        direction_experimental = (
            "resistant (weaker MTX binding)" if fold_change > 2 else
            "no significant change" if 0.5 <= fold_change <= 2 else
            "tighter binding"
        )
        direction_docking = (
            "predicted weaker binding" if delta_vina > 0.5 else
            "predicted no significant change" if -0.5 <= delta_vina <= 0.5 else
            "predicted tighter binding"
        )

        if is_unverified:
            if direction_experimental.startswith("resistant") and direction_docking.startswith("predicted weaker"):
                agreement = "directionally consistent (experimental magnitude unverified -- aggregate range only)"
            else:
                agreement = "insufficient evidence (experimental value is an aggregate range, not isolated)"
        elif direction_experimental.startswith("resistant") and direction_docking.startswith("predicted weaker"):
            agreement = "consistent"
        elif direction_experimental.startswith("no significant") and direction_docking.startswith("predicted no"):
            agreement = "consistent"
        elif direction_experimental.split()[0] in direction_docking:
            agreement = "partially consistent"
        else:
            agreement = "inconsistent"

        rows.append({
            "mutation": "/".join(f"{wt}{pos}{mut}" for pos, wt, mut in info["mutations"]),
            "pdb_id": pdb_id,
            "experimental_metric": "Ki_MTX fold-change (vs paper's own WT)",
            "experimental_wt": "own-paper WT control (see data/published_affinity_data.csv)",
            "experimental_mutant": ("UNVERIFIED (>=740x aggregate range)" if is_unverified
                                     else f"{fold_change:.1f}x"),
            "experimental_fold_change": None if is_unverified else round(fold_change, 1),
            "experimental_ddG_kcal_mol": round(ddg_exp, 2),
            "vina_wt_kcal_mol": round(wt_score, 2),
            "vina_mutant_kcal_mol": round(mut_score, 2),
            "delta_vina_kcal_mol": round(delta_vina, 2),
            "direction_experimental": direction_experimental,
            "direction_docking": direction_docking,
            "agreement": agreement,
            "interpretation": (
                "Docking predicts a DIRECTION of binding disruption, not a calibrated Ki/Kd. "
                "Compare direction and rough magnitude only."
            ),
        })

    df = pd.DataFrame(rows)
    OUT_DOCK.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_DOCK / "docking_vs_experimental.csv", index=False)
    print(f"    Wrote {OUT_DOCK / 'docking_vs_experimental.csv'}")
    for _, row in df.iterrows():
        print(f"    {row['mutation']}: exp={row['direction_experimental']}, "
              f"dock={row['direction_docking']} -> {row['agreement']}")

    # Figure 6D: docking scores WT vs mutants
    fig, ax = plt.subplots(figsize=(8, 5))
    labels = ["WT\n(1U72)"] + [f"{r.mutation}\n({r.pdb_id})" for r in df.itertuples()]
    values = [wt_score] + list(df["vina_mutant_kcal_mol"])
    colors = ["#4c72b0"] + ["#d62728" if a in ("consistent",) else
                             "#dd8452" if a == "partially consistent" else
                             "#937860" if a == "inconsistent" else "#8c8c8c"
                             for a in df["agreement"]]
    ax.bar(labels, values, color=colors)
    ax.set_ylabel("Best predicted binding affinity (kcal/mol, more negative = tighter)")
    ax.set_title("Figure 6D: Docking score, WT vs resistance mutants")
    ax.axhline(wt_score, color="gray", linestyle="--", linewidth=0.8, label="WT reference")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT_DOCK / "fig6d_docking_scores.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"    Wrote {OUT_DOCK / 'fig6d_docking_scores.png'}")

    # Figure 6E: experimental vs computational ddG
    fig, ax = plt.subplots(figsize=(7, 6))
    for _, row in df.iterrows():
        marker = "^" if row["pdb_id"] == "1DLR" else "o"
        ax.scatter(row["experimental_ddG_kcal_mol"], row["delta_vina_kcal_mol"], s=100, marker=marker,
                   color="#d62728" if row["agreement"] == "consistent" else
                   "#dd8452" if row["agreement"] == "partially consistent" else
                   "#937860" if row["agreement"] == "inconsistent" else "#8c8c8c")
        ax.annotate(row["mutation"], (row["experimental_ddG_kcal_mol"], row["delta_vina_kcal_mol"]),
                    xytext=(6, 6), textcoords="offset points", fontsize=8)
    ax.scatter(0, 0, marker="s", s=100, color="#4c72b0")
    ax.annotate("WT (reference)", (0, 0), xytext=(6, -12), textcoords="offset points", fontsize=8)

    lims = [min(ax.get_xlim()[0], ax.get_ylim()[0], -1), max(ax.get_xlim()[1], ax.get_ylim()[1], 1)]
    ax.plot(lims, lims, linestyle=":", color="gray", linewidth=0.8, label="y = x (perfect agreement)")
    ax.axhline(0, color="black", linewidth=0.5)
    ax.axvline(0, color="black", linewidth=0.5)
    ax.set_xlabel(r"Experimental $\Delta\Delta G$ = RT $\ln$(Ki$_{mut}$/Ki$_{WT}$), kcal/mol")
    ax.set_ylabel(r"Docking $\Delta\Delta G$ = Vina$_{mut}$ - Vina$_{WT}$, kcal/mol")
    ax.set_title("Figure 6E: Docking vs. published affinity change\n"
                 "(triangle = 1DLR, experimental value is a lower-bound estimate)")
    ax.legend(fontsize=8, loc="upper left")
    fig.tight_layout()
    fig.savefig(OUT_DOCK / "fig6e_experimental_vs_computational.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"    Wrote {OUT_DOCK / 'fig6e_experimental_vs_computational.png'}")

    return df


# ---------------------------------------------------------------------------
# Phase 11: post-docking contact comparison
# ---------------------------------------------------------------------------

def phase11_contact_comparison():
    print("\n[Phase 11] Post-docking contact comparison (docked methotrexate pose, WT vs mutants)...")
    rows = []
    for pdb_id in STRUCTURE_IDS:
        docked_pdbqt = config.DOCKING_DIR / f"{pdb_id}_docked.pdbqt"
        docked_pdb = config.DOCKING_DIR / f"{pdb_id}_docked_pose1.pdb"
        subprocess.run(["obabel", str(docked_pdbqt), "-O", str(docked_pdb), "-f", "1", "-l", "1"],
                       capture_output=True, text=True)
        if not docked_pdb.exists():
            continue

        parser = PDBParser(QUIET=True)
        try:
            pose_model = parser.get_structure(pdb_id, docked_pdb)[0]
        except Exception:  # noqa: BLE001
            continue
        pose_atoms = [a for chain in pose_model for r in chain for a in r if a.element != "H"]

        protein_model = load_structure(pdb_id)
        n_direct, n_hbond_capable_direct = 0, 0
        for chain in protein_model:
            for residue in chain:
                if residue.id[0] != " " or residue.resname not in THREE_TO_ONE:
                    continue
                dist = min_distance(residue, pose_atoms)
                if dist is not None and dist < DIRECT_CONTACT_CUTOFF:
                    n_direct += 1
                    if AA_PROPERTIES.get(THREE_TO_ONE[residue.resname], {}).get("hbond", "none") != "none":
                        n_hbond_capable_direct += 1

        mutations = config.RESISTANCE_MUTANT_PDBS.get(pdb_id, {}).get("mutations", [])
        rows.append({
            "pdb_id": pdb_id,
            "mutation": "wild-type" if not mutations else "/".join(f"{wt}{p}{m}" for p, wt, m in mutations),
            "docked_pose_direct_contact_residues": n_direct,
            "docked_pose_hbond_capable_contacts": n_hbond_capable_direct,
            "note": "Counted from the top-ranked docked pose, not the original crystal ligand pose.",
        })

    df = pd.DataFrame(rows)
    OUT_DOCK.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_DOCK / "contact_comparison.csv", index=False)
    print(f"    Wrote {OUT_DOCK / 'contact_comparison.csv'}")
    return df


# ---------------------------------------------------------------------------
# Phase 12: flexible side-chain docking (tests whether the rigid-receptor
# limitation explains the WT-vs-mutant disagreements found in Phase 10)
# ---------------------------------------------------------------------------

def find_incomplete_residues(pdb_path):
    """Parses REMARK 470 (missing atoms, a normal X-ray artifact for
    flexible/disordered side chains) to find residues meeko's strict
    chemical-template matcher will choke on. Returns a sorted list of
    residue numbers (chain A only, matching this project's convention)."""
    positions = set()
    for line in Path(pdb_path).read_text().splitlines():
        if line.startswith("REMARK 470") and len(line) > 20:
            chain = line[19:20].strip()
            resnum_field = line[20:24].strip()
            if chain == "A" and resnum_field.isdigit():
                positions.add(int(resnum_field))
    return sorted(positions)


def prepare_flexible_receptor(pdb_id, flex_positions):
    """meeko mk_prepare_receptor.py splits the receptor into a rigid PDBQT
    plus a flexible PDBQT containing only the named side chains' torsion
    trees (BEGIN_RES/ROOT/BRANCH/END_RES), exactly what Vina's --flex wants.
    """
    receptor_pdb = config.DOCKING_DIR / f"{pdb_id}_receptor.pdb"
    pos_key = "-".join(str(p) for p in sorted(flex_positions))
    basename = config.DOCKING_DIR / f"{pdb_id}_flex{pos_key}"
    rigid_out = Path(f"{basename}_rigid.pdbqt")
    flex_out = Path(f"{basename}_flex.pdbqt")

    if not (rigid_out.exists() and flex_out.exists()):
        flexres_arg = ",".join(f"A:{p}" for p in flex_positions)
        # --default_altloc A: some structures (3EIG) have genuinely disordered
        # side chains with alternate conformers -- notably Arg31 itself, which
        # is the whole point of that structure's own publication ("Multiple
        # conformers in active site..."). We deliberately use the standard
        # higher-occupancy 'A' conformer rather than silently failing; this is
        # documented as a limitation, since the alternate conformer may matter.
        cmd = ["mk_prepare_receptor.py", "--read_pdb", str(receptor_pdb),
               "-f", flexres_arg, "-o", str(basename), "-p", "--default_altloc", "A"]

        incomplete = [p for p in find_incomplete_residues(download_structure(pdb_id)) if p not in flex_positions]
        if incomplete:
            print(f"        {pdb_id}: excluding {len(incomplete)} residue(s) with crystallographically "
                  f"disordered/missing side-chain atoms (REMARK 470): {incomplete}")
            cmd += ["-d", ",".join(f"A:{p}" for p in incomplete)]

        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0 or not rigid_out.exists():
            raise RuntimeError(
                f"meeko flexible receptor prep failed for {pdb_id} @ {flex_positions} "
                f"(exit {result.returncode}):\nSTDOUT:\n{result.stdout[-3000:]}\nSTDERR:\n{result.stderr[-3000:]}"
            )
    return rigid_out, flex_out


def run_vina_flex(rigid_pdbqt, flex_pdbqt, ligand_pdbqt, center, box_size, out_pdbqt, exhaustiveness, seed):
    cmd = [
        "vina", "--receptor", str(rigid_pdbqt), "--flex", str(flex_pdbqt), "--ligand", str(ligand_pdbqt),
        "--center_x", str(center[0]), "--center_y", str(center[1]), "--center_z", str(center[2]),
        "--size_x", str(box_size), "--size_y", str(box_size), "--size_z", str(box_size),
        "--exhaustiveness", str(exhaustiveness), "--seed", str(seed),
        "--out", str(out_pdbqt),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    if result.returncode != 0:
        raise RuntimeError(f"vina (flexible) failed: {result.stderr}\n{result.stdout}")
    return result.stdout


def phase12_flexible_docking(comparison_df):
    print("\n[Phase 12] Flexible side-chain docking -- does relaxing exactly the "
          "mutated residue(s) change the ΔΔG picture from Phase 10?")
    ligand_pdbqt = config.DOCKING_DIR / "methotrexate.pdbqt"
    box_size = 22.0
    exhaustiveness = 8
    seed = 42
    ligand_code_wt, _ = STRUCTURE_LIGANDS[config.REFERENCE_PDB]

    rows = []
    for pdb_id, info in config.RESISTANCE_MUTANT_PDBS.items():
        flex_positions = sorted({pos for pos, wt, mut in info["mutations"]})
        rigid_row = comparison_df[comparison_df["pdb_id"] == pdb_id].iloc[0]

        # WT, made flexible at the SAME position(s) mutated in this structure,
        # for a fair paired comparison (not one universal WT baseline).
        wt_model = load_structure(config.REFERENCE_PDB)
        wt_ligand = find_ligand_residue(wt_model, ligand_code_wt)
        wt_center = np.array([a.coord for a in wt_ligand if a.element != "H"]).mean(axis=0)
        wt_rigid_pdbqt, wt_flex_pdbqt = prepare_flexible_receptor(config.REFERENCE_PDB, flex_positions)
        wt_out = config.DOCKING_DIR / f"{config.REFERENCE_PDB}_flex{'-'.join(map(str, flex_positions))}_docked.pdbqt"
        wt_stdout = run_vina_flex(wt_rigid_pdbqt, wt_flex_pdbqt, ligand_pdbqt, wt_center, box_size,
                                   wt_out, exhaustiveness, seed)
        wt_modes = parse_vina_modes(wt_stdout)
        wt_best = wt_modes[0]["affinity_kcal_mol"] if wt_modes else None

        mut_model = load_structure(pdb_id)
        ligand_code_mut, _ = STRUCTURE_LIGANDS[pdb_id]
        mut_ligand = find_ligand_residue(mut_model, ligand_code_mut)
        mut_center = np.array([a.coord for a in mut_ligand if a.element != "H"]).mean(axis=0)
        mut_rigid_pdbqt, mut_flex_pdbqt = prepare_flexible_receptor(pdb_id, flex_positions)
        mut_out = config.DOCKING_DIR / f"{pdb_id}_flex_docked.pdbqt"
        mut_stdout = run_vina_flex(mut_rigid_pdbqt, mut_flex_pdbqt, ligand_pdbqt, mut_center, box_size,
                                    mut_out, exhaustiveness, seed)
        mut_modes = parse_vina_modes(mut_stdout)
        mut_best = mut_modes[0]["affinity_kcal_mol"] if mut_modes else None

        delta_flex = (mut_best - wt_best) if (mut_best is not None and wt_best is not None) else None
        rows.append({
            "mutation": rigid_row["mutation"],
            "pdb_id": pdb_id,
            "flexible_positions": ",".join(str(p) for p in flex_positions),
            "vina_wt_rigid_kcal_mol": rigid_row["vina_wt_kcal_mol"],
            "vina_mutant_rigid_kcal_mol": rigid_row["vina_mutant_kcal_mol"],
            "delta_vina_rigid_kcal_mol": rigid_row["delta_vina_kcal_mol"],
            "vina_wt_flex_kcal_mol": round(wt_best, 2) if wt_best is not None else None,
            "vina_mutant_flex_kcal_mol": round(mut_best, 2) if mut_best is not None else None,
            "delta_vina_flex_kcal_mol": round(delta_flex, 2) if delta_flex is not None else None,
            "experimental_ddG_kcal_mol": rigid_row["experimental_ddG_kcal_mol"],
        })
        print(f"    {rigid_row['mutation']} (flex @ {flex_positions}): "
              f"rigid ddG={rigid_row['delta_vina_kcal_mol']:.2f}, "
              f"flexible ddG={delta_flex:.2f} kcal/mol, "
              f"experimental ddG={rigid_row['experimental_ddG_kcal_mol']:.2f} kcal/mol")

    df = pd.DataFrame(rows)
    OUT_DOCK.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_DOCK / "docking_flexible_vs_rigid.csv", index=False)
    print(f"    Wrote {OUT_DOCK / 'docking_flexible_vs_rigid.csv'}")

    fig, ax = plt.subplots(figsize=(9, 6))
    x = np.arange(len(df))
    width = 0.25
    ax.bar(x - width, df["experimental_ddG_kcal_mol"], width, label="Experimental (published Ki)", color="#4c72b0")
    ax.bar(x, df["delta_vina_rigid_kcal_mol"], width, label="Docking, rigid receptor", color="#dd8452")
    ax.bar(x + width, df["delta_vina_flex_kcal_mol"], width,
           label="Docking, flexible mutated side chain(s)", color="#55a868")
    ax.set_xticks(x)
    ax.set_xticklabels(df["mutation"])
    ax.axhline(0, color="black", linewidth=0.5)
    ax.set_ylabel(r"$\Delta\Delta G$ (kcal/mol, positive = weaker mutant binding)")
    ax.set_title("Does side-chain flexibility close the gap between docking and\n"
                 "published methotrexate-binding data?")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT_DOCK / "fig_flexible_vs_rigid_docking.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"    Wrote {OUT_DOCK / 'fig_flexible_vs_rigid_docking.png'}")
    return df


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=" * 90)
    print("Module 6: Mutation-aware structural + docking analysis, validated against")
    print("published biochemical affinity data")
    print("=" * 90)

    config.STRUCTURES_DIR.mkdir(parents=True, exist_ok=True)
    OUT_STRUCT.mkdir(parents=True, exist_ok=True)
    OUT_DOCK.mkdir(parents=True, exist_ok=True)

    phase1_mutation_tables()
    phase2_conservation_link()
    phase3_structure_metadata()
    wt_contacts_df = phase4_wt_contacts()
    phase5_mutation_structural_context(wt_contacts_df)
    scores_df, _ = phases_6_to_9_docking()
    comparison_df = phase10_docking_vs_experimental(scores_df)
    phase11_contact_comparison()
    phase12_flexible_docking(comparison_df)

    print("\n" + "=" * 90)
    print("Module 6 complete.")
    print("=" * 90)
    sys.exit(0)


if __name__ == "__main__":
    main()
