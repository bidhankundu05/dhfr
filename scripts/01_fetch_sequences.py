"""
Module 1: Sequence retrieval.

Fetches the DHFR protein sequence for human (direct UniProt accession) and
each other species in config.SPECIES_LIST (NCBI Entrez RefSeq search, with a
UniProt REST search fallback). Writes a combined FASTA and a per-species
retrieval log.
"""

import csv
import sys
import textwrap
import time
from io import StringIO
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import config  # noqa: E402

import requests  # noqa: E402
from Bio import Entrez, SeqIO  # noqa: E402

UNIPROT_ACCESSION_URL = "https://rest.uniprot.org/uniprotkb/{accession}.fasta"
UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
REQUEST_TIMEOUT = 30
FASTA_WRAP = 70


def configure_entrez():
    if not config.ENTREZ_EMAIL or config.ENTREZ_EMAIL == "your_email@example.com":
        print(
            "ERROR: config.ENTREZ_EMAIL is not set (still the placeholder). "
            "NCBI's usage policy requires a real email on every Entrez request. "
            "Set ENTREZ_EMAIL in config.py or as an environment variable and re-run."
        )
        sys.exit(1)
    Entrez.email = config.ENTREZ_EMAIL
    if config.ENTREZ_API_KEY:
        Entrez.api_key = config.ENTREZ_API_KEY
    # NCBI rate limits: 3 req/s without an API key, 10 req/s with one.
    return 0.11 if config.ENTREZ_API_KEY else 0.34


def fetch_uniprot_by_accession(accession):
    url = UNIPROT_ACCESSION_URL.format(accession=accession)
    resp = requests.get(url, timeout=REQUEST_TIMEOUT)
    resp.raise_for_status()
    if not resp.text.strip():
        raise ValueError(f"empty response for UniProt accession {accession}")
    record = next(SeqIO.parse(StringIO(resp.text), "fasta"))
    parts = record.id.split("|")
    resolved_accession = parts[1] if len(parts) >= 2 else record.id
    return resolved_accession, str(record.seq)


def gene_symbols_for(species):
    """GENE_SYMBOL first, then any known species-specific aliases."""
    return [config.GENE_SYMBOL] + config.GENE_SYMBOL_ALIASES.get(species, [])


def pick_curated_id(id_list, summaries_by_id):
    """Prefer curated RefSeq (NP_/YP_/AP_) over computationally predicted
    models (XP_/XM_, NCBI's Gnomon pipeline), which can be partial or
    alternate-isoform calls. Falls back to the search engine's own top hit
    if no curated accession is present."""
    for uid in id_list:
        summary = summaries_by_id.get(uid)
        if summary and str(summary.get("AccessionVersion", "")).startswith(("NP_", "YP_", "AP_")):
            return uid
    return id_list[0]


def fetch_entrez_refseq(species, gene_symbol, request_delay):
    term = f"{gene_symbol}[Gene Name] AND {species}[Organism] AND refseq[filter]"
    handle = Entrez.esearch(db="protein", term=term, retmax=20)
    search_record = Entrez.read(handle)
    handle.close()
    time.sleep(request_delay)

    id_list = search_record.get("IdList", [])
    if not id_list:
        return None

    if len(id_list) > 1:
        handle = Entrez.esummary(db="protein", id=",".join(id_list))
        summaries = Entrez.read(handle)
        handle.close()
        time.sleep(request_delay)
        summaries_by_id = {str(s.get("Id")): s for s in summaries}
        best_id = pick_curated_id(id_list, summaries_by_id)
    else:
        best_id = id_list[0]

    handle = Entrez.efetch(db="protein", id=best_id, rettype="fasta", retmode="text")
    fasta_text = handle.read()
    handle.close()
    time.sleep(request_delay)

    if not fasta_text.strip():
        return None

    record = next(SeqIO.parse(StringIO(fasta_text), "fasta"))
    return record.id, str(record.seq)


def fetch_uniprot_search_fallback(species, gene_symbol):
    query = f'gene:{gene_symbol} AND organism_name:"{species}"'
    params = {"query": query, "format": "fasta", "size": 10}
    resp = requests.get(UNIPROT_SEARCH_URL, params=params, timeout=REQUEST_TIMEOUT)
    resp.raise_for_status()
    if not resp.text.strip():
        return None

    records = list(SeqIO.parse(StringIO(resp.text), "fasta"))
    if not records:
        return None

    # Prefer reviewed (sp|, Swiss-Prot) over unreviewed (tr|, TrEMBL) hits,
    # and among reviewed hits skip fragment/partial entries when a
    # non-fragment alternative exists.
    reviewed = [r for r in records if r.id.startswith("sp|")]
    pool = reviewed if reviewed else records
    non_fragment = [r for r in pool if "Fragment" not in r.description]
    chosen = non_fragment[0] if non_fragment else pool[0]

    parts = chosen.id.split("|")
    accession = parts[1] if len(parts) >= 2 else chosen.id
    return accession, str(chosen.seq)


def fetch_species_sequence(species, request_delay):
    """Returns (accession, sequence, source_database) or (None, None, None).

    Tries NCBI Entrez RefSeq across every known gene symbol for this species
    (GENE_SYMBOL, then any GENE_SYMBOL_ALIASES) before ever touching the
    UniProt fallback. This matters: e.g. for E. coli, a UniProt search for
    "DHFR" alone turns up dfrA12, an unrelated plasmid-borne trimethoprim-
    resistance gene, well before the pipeline would try the correct
    chromosomal alias "folA" -- so all Entrez/RefSeq attempts (higher-
    curation, alias-aware) are exhausted first.
    """
    symbols = gene_symbols_for(species)

    for gene_symbol in symbols:
        label = gene_symbol if gene_symbol == config.GENE_SYMBOL else f"{gene_symbol} (alias)"
        try:
            entrez_result = fetch_entrez_refseq(species, gene_symbol, request_delay)
        except Exception as exc:  # noqa: BLE001
            print(f"    NCBI Entrez lookup ({label}) failed: {exc}")
            entrez_result = None

        if entrez_result:
            accession, sequence = entrez_result
            return accession, sequence, "NCBI RefSeq (Entrez)"
        print(f"    No Entrez RefSeq hit for gene symbol '{label}'")

    print("    Entrez exhausted, falling back to UniProt search...")
    for gene_symbol in symbols:
        label = gene_symbol if gene_symbol == config.GENE_SYMBOL else f"{gene_symbol} (alias)"
        try:
            uniprot_result = fetch_uniprot_search_fallback(species, gene_symbol)
        except Exception as exc:  # noqa: BLE001
            print(f"    UniProt fallback search ({label}) failed: {exc}")
            uniprot_result = None

        if uniprot_result:
            accession, sequence = uniprot_result
            return accession, sequence, "UniProt (search fallback)"
        print(f"    No UniProt hit for gene symbol '{label}' either")

    return None, None, None


def write_fasta(records, out_path):
    with open(out_path, "w") as fh:
        for rec in records:
            if rec["sequence"] is None:
                continue
            header_species = rec["species"].replace(" ", "_")
            fh.write(f">{header_species}|{rec['accession']}\n")
            for line in textwrap.wrap(rec["sequence"], FASTA_WRAP):
                fh.write(line + "\n")


def write_log(records, out_path):
    with open(out_path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["species", "accession", "source_database", "sequence_length", "status"])
        for rec in records:
            writer.writerow(
                [
                    rec["species"],
                    rec["accession"] or "",
                    rec["source"] or "",
                    rec["length"] if rec["length"] is not None else "",
                    "OK" if rec["sequence"] is not None else "FAILED",
                ]
            )


def print_summary(records):
    print("\n" + "=" * 78)
    print("SUMMARY")
    print("=" * 78)
    species_width = max([len(rec["species"]) for rec in records] + [len("Species")]) + 2
    source_width = max([len(rec["source"] or "-") for rec in records] + [len("Source")]) + 2
    header = f"{'Species':<{species_width}}{'Accession':<16}{'Source':<{source_width}}{'Length':<8}{'Status'}"
    print(header)
    print("-" * len(header))
    n_ok = 0
    for rec in records:
        status = "OK" if rec["sequence"] is not None else "FAILED"
        n_ok += status == "OK"
        print(
            f"{rec['species']:<{species_width}}"
            f"{(rec['accession'] or '-'):<16}"
            f"{(rec['source'] or '-'):<{source_width}}"
            f"{(str(rec['length']) if rec['length'] is not None else '-'):<8}"
            f"{status}"
        )
    print("-" * len(header))
    print(f"{n_ok}/{len(records)} species retrieved successfully")


def main():
    request_delay = configure_entrez()

    config.SEQUENCES_DIR.mkdir(parents=True, exist_ok=True)
    config.RESULTS_SEQUENCES_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 78)
    print(f"Module 1: Fetching {config.GENE_SYMBOL} sequences for {len(config.SPECIES_LIST)} species")
    print("=" * 78)

    records = []
    for species in config.SPECIES_LIST:
        print(f"\n[{species}]")
        if species == "Homo sapiens":
            try:
                accession, sequence = fetch_uniprot_by_accession(config.HUMAN_UNIPROT)
                source = "UniProt (direct accession)"
                print(f"    Fetched {accession} directly from UniProt ({len(sequence)} aa)")
            except Exception as exc:  # noqa: BLE001
                print(f"    ERROR fetching human sequence from UniProt: {exc}")
                accession, sequence, source = None, None, None
        else:
            accession, sequence, source = fetch_species_sequence(species, request_delay)
            if sequence:
                print(f"    Fetched {accession} via {source} ({len(sequence)} aa)")
            else:
                print("    FAILED to retrieve a sequence for this species")

        records.append(
            {
                "species": species,
                "accession": accession,
                "sequence": sequence,
                "source": source,
                "length": len(sequence) if sequence else None,
            }
        )

    fasta_out = config.SEQUENCES_DIR / "dhfr_multispecies.fasta"
    write_fasta(records, fasta_out)
    print(f"\nWrote combined FASTA: {fasta_out}")

    log_out = config.RESULTS_SEQUENCES_DIR / "retrieval_log.csv"
    write_log(records, log_out)
    print(f"Wrote retrieval log: {log_out}")

    print_summary(records)

    n_failed = sum(1 for rec in records if rec["sequence"] is None)
    sys.exit(1 if n_failed > 0 else 0)


if __name__ == "__main__":
    main()
