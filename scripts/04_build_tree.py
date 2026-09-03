"""
Module 4: Phylogenetic tree construction.

Builds a tree from the Module 3 alignment (results/alignment/dhfr_aligned.fasta)
using FastTree/VeryFastTree (approximate ML) via subprocess, falling back to
Biopython's own neighbor-joining (Bio.Phylo.TreeConstruction) if neither CLI
tool is installed. Roots the tree on E. coli (the most divergent taxon --
confirmed in Modules 2/3) as outgroup, saves Newick + a rendered PNG, and
checks the resulting topology against textbook species relationships.
"""

import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import config  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from Bio import AlignIO, Phylo  # noqa: E402
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor  # noqa: E402

ALIGNMENT_FASTA = config.RESULTS_ALIGNMENT_DIR / "dhfr_aligned.fasta"
RAW_TREE_OUT = config.TREE_DIR / "dhfr_tree_raw.newick"
FINAL_TREE_OUT = config.RESULTS_TREE_DIR / "dhfr_tree.newick"
TREE_IMAGE_OUT = config.RESULTS_TREE_DIR / "phylo_tree.png"
TOPOLOGY_NOTES_OUT = config.RESULTS_TREE_DIR / "topology_notes.txt"

OUTGROUP_SPECIES = "Escherichia coli"


def find_tree_tool():
    for exe in ("fasttree", "veryfasttree"):
        path = shutil.which(exe)
        if path:
            return exe, path
    return None, None


def run_fasttree(exe, alignment_fasta, output_newick):
    result = subprocess.run(
        [exe, "-out", str(output_newick), "-quiet", str(alignment_fasta)],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(f"{exe} exited {result.returncode}: {result.stderr}")


def run_biopython_nj(alignment_fasta, output_newick):
    alignment = AlignIO.read(alignment_fasta, "fasta")
    calculator = DistanceCalculator("blosum62")
    dm = calculator.get_distance(alignment)
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm)
    Phylo.write(tree, output_newick, "newick")


def build_tree():
    config.TREE_DIR.mkdir(parents=True, exist_ok=True)
    config.RESULTS_TREE_DIR.mkdir(parents=True, exist_ok=True)

    exe, path = find_tree_tool()
    if exe:
        label = "FastTree" if exe == "fasttree" else "VeryFastTree (FastTree-compatible substitute)"
        print(f"Using {label}: {path}")
        run_fasttree(exe, ALIGNMENT_FASTA, RAW_TREE_OUT)
        method = label
    else:
        print("Neither fasttree nor veryfasttree found on PATH -- "
              "falling back to Biopython neighbor-joining (BLOSUM62 distances)...")
        run_biopython_nj(ALIGNMENT_FASTA, RAW_TREE_OUT)
        method = "Biopython neighbor-joining (BLOSUM62 distance)"

    return method


def clean_and_root_tree():
    tree = Phylo.read(RAW_TREE_OUT, "newick")

    # Raw tip names are "Species_Name|Accession" (from the alignment headers);
    # keep the full label in the raw tool output but simplify to species-only
    # names for the cleaned, final Newick/PNG deliverables.
    for tip in tree.get_terminals():
        tip.name = tip.name.split("|")[0]

    outgroup_name = OUTGROUP_SPECIES.replace(" ", "_")
    outgroup_clade = tree.find_any(name=outgroup_name)
    if outgroup_clade is not None:
        tree.root_with_outgroup(outgroup_clade)
    else:
        print(f"WARNING: could not find '{outgroup_name}' to root as outgroup; leaving tree as-is.")

    Phylo.write(tree, FINAL_TREE_OUT, "newick")
    return tree


def render_tree_image(tree):
    fig, ax = plt.subplots(figsize=(9, 4.5))
    display_tree = tree
    for clade in display_tree.find_clades():
        if clade.name:
            clade.name = clade.name.replace("_", " ")
    Phylo.draw(display_tree, axes=ax, do_show=False)
    ax.set_title("DHFR phylogenetic tree (rooted on E. coli outgroup)")
    fig.tight_layout()
    fig.savefig(TREE_IMAGE_OUT, dpi=200)
    plt.close(fig)


def check_topology(tree):
    """Compares the constructed gene tree against textbook organismal
    relationships for this species set and writes a short discussion."""

    def monophyletic(names):
        targets = [tree.find_any(name=n) for n in names]
        if any(t is None for t in targets):
            return None
        return tree.is_monophyletic(targets) is not False

    checks = [
        (
            "Mammals form a clade",
            ["Homo sapiens", "Mus musculus", "Rattus norvegicus"],
            "Human, mouse, and rat are all placental mammals and are expected to "
            "share a more recent common ancestor with each other than with fish, "
            "yeast, or bacteria.",
        ),
        (
            "Rodents are sister taxa",
            ["Mus musculus", "Rattus norvegicus"],
            "Mouse and rat are both murid rodents and are expected to be each "
            "other's closest relative in the tree.",
        ),
        (
            "Vertebrates (mammals + fish) form a clade excluding yeast/bacteria",
            ["Homo sapiens", "Mus musculus", "Rattus norvegicus", "Danio rerio"],
            "Zebrafish is a vertebrate outgroup to mammals, but all four should "
            "still group together to the exclusion of the fungal (yeast) and "
            "bacterial (E. coli) sequences.",
        ),
        (
            "Eukaryotes form a clade excluding the bacterial outgroup",
            ["Homo sapiens", "Mus musculus", "Rattus norvegicus", "Danio rerio", "Saccharomyces cerevisiae"],
            "All eukaryotic sequences (animals + yeast) should exclude E. coli, "
            "which is rooted as the outgroup precisely because bacteria diverged "
            "from eukaryotes earliest of this set.",
        ),
    ]

    lines = [
        "DHFR phylogenetic tree -- topology vs. known species relationships",
        "=" * 70,
        f"Rooted on: {OUTGROUP_SPECIES} (outgroup; most divergent DHFR sequence per "
        "Modules 2-3 BLAST identity / BLOSUM62 similarity)",
        "",
    ]
    for title, names, rationale in checks:
        result = monophyletic(names)
        status = "MATCHES expected topology" if result else (
            "does NOT match (single-gene tree disagrees with species tree here)"
            if result is False else "could not evaluate (a taxon name was not found in the tree)"
        )
        lines.append(f"* {title}: {status}")
        lines.append(f"  Expectation: {rationale}")
        lines.append("")

    lines.append(
        "Note: a tree built from one gene (DHFR) is a *gene tree*, not the "
        "species tree itself -- they usually agree closely for a slowly-"
        "evolving, single-copy, non-hybridizing gene like DHFR, but rate "
        "variation or incomplete lineage sorting can occasionally cause "
        "local disagreements. Any 'does NOT match' result above is worth "
        "discussing rather than treating as an error."
    )

    TOPOLOGY_NOTES_OUT.write_text("\n".join(lines))
    return lines


def main():
    print("=" * 90)
    print("Module 4: Phylogenetic tree construction")
    print("=" * 90)

    if not ALIGNMENT_FASTA.exists():
        print(f"ERROR: {ALIGNMENT_FASTA} not found. Run scripts/03_align_mafft.py first.")
        sys.exit(1)

    method = build_tree()
    print(f"Wrote raw tree: {RAW_TREE_OUT}")

    tree = clean_and_root_tree()
    print(f"Wrote final rooted tree: {FINAL_TREE_OUT}")

    render_tree_image(tree)
    print(f"Wrote tree image: {TREE_IMAGE_OUT}")

    print("\nChecking topology against known species relationships...")
    notes_lines = check_topology(tree)
    print(f"Wrote topology notes: {TOPOLOGY_NOTES_OUT}\n")
    for line in notes_lines:
        print(line)

    print("\n" + "=" * 90)
    print("SUMMARY")
    print("=" * 90)
    print(f"Method: {method}")
    print(f"Taxa: {len(tree.get_terminals())}")
    print(f"Rooted on: {OUTGROUP_SPECIES}")

    sys.exit(0)


if __name__ == "__main__":
    main()
