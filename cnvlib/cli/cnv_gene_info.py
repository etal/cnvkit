#!/usr/bin/env python
"""Build a CNVkit gene-info table from an Ensembl BioMart export.

The RNA copy-number pipeline (``cnvkit.py import-rna``) needs a per-gene
metadata table. CNVkit bundles one for the human transcriptome (hg38) as
``data/ensembl-gene-info.hg38.tsv``, loaded by ``cnvlib.rna.load_gene_info``.
This script normalizes a BioMart TSV export for *any* reference genome into the
same canonical format, so the RNA pipeline can be run on hg19, mouse
(mm10/mm39), or other assemblies.

Export the following attributes from Ensembl BioMart (Gene attributes; column
order does not matter -- columns are matched by name):

    Gene stable ID
    Gene % GC content
    Chromosome/scaffold name
    Gene start (bp)
    Gene end (bp)
    Gene name
    NCBI gene ID
    Transcript length (including UTRs and CDS)
    Transcript support level (TSL)

Of these, ``Gene name``, ``NCBI gene ID``, and ``Transcript support level`` are
optional and filled with defaults when absent; the rest are required.

Columns are matched case-insensitively against known BioMart display names and
short aliases. Use ``--rename "Source header=canonical_key"`` to map any column
the script does not recognize (keys: gene_id, gc, chromosome, start, end, gene,
entrez_id, tx_length, tx_support).
"""

import argparse
import logging
import os
import sys
from contextlib import nullcontext

import pandas as pd

from skgenome.chromsort import sorter_chrom

# Canonical output columns, in order. The first element is the short key used
# internally and by ``cnvlib.rna.load_gene_info``; the second is the BioMart
# display name written to the header line of the output file.
CANONICAL_COLUMNS = [
    ("gene_id", "Gene stable ID"),
    ("gc", "Gene % GC content"),
    ("chromosome", "Chromosome/scaffold name"),
    ("start", "Gene start (bp)"),
    ("end", "Gene end (bp)"),
    ("gene", "Gene name"),
    ("entrez_id", "NCBI gene ID"),
    ("tx_length", "Transcript length (including UTRs and CDS)"),
    ("tx_support", "Transcript support level (TSL)"),
]
CANONICAL_KEYS = [key for key, _display in CANONICAL_COLUMNS]
# Columns CNVkit cannot function without (coordinates, GC for bias correction,
# transcript length for depth normalization). The rest get sensible defaults.
REQUIRED_KEYS = {"gene_id", "gc", "chromosome", "start", "end", "tx_length"}
DEFAULTS = {"gene": "", "entrez_id": "", "tx_support": "tslNA"}

# Map a normalized (lowercased, whitespace-collapsed) input header to its
# canonical key. Covers BioMart display-name variants across Ensembl releases
# plus common short aliases.
COLUMN_ALIASES = {
    # gene_id
    "gene stable id": "gene_id",
    "gene stable id version": "gene_id",
    "ensembl gene id": "gene_id",
    "gene id": "gene_id",
    "gene_id": "gene_id",
    "geneid": "gene_id",
    # gc
    "gene % gc content": "gc",
    "% gc content": "gc",
    "gc content": "gc",
    "percentage gc content": "gc",
    "gc": "gc",
    # chromosome
    "chromosome/scaffold name": "chromosome",
    "chromosome name": "chromosome",
    "chromosome": "chromosome",
    "chrom": "chromosome",
    "chr": "chromosome",
    "seqname": "chromosome",
    # start
    "gene start (bp)": "start",
    "gene start": "start",
    "start (bp)": "start",
    "start": "start",
    # end
    "gene end (bp)": "end",
    "gene end": "end",
    "end (bp)": "end",
    "end": "end",
    # gene name
    "gene name": "gene",
    "external gene name": "gene",
    "hgnc symbol": "gene",
    "gene_name": "gene",
    "symbol": "gene",
    "name": "gene",
    # entrez_id
    "ncbi gene id": "entrez_id",
    "ncbi gene (formerly entrezgene) id": "entrez_id",
    "entrezgene id": "entrez_id",
    "entrez gene id": "entrez_id",
    "entrezgene": "entrez_id",
    "entrez id": "entrez_id",
    "entrez_id": "entrez_id",
    "ncbi": "entrez_id",
    # tx_length
    "transcript length (including utrs and cds)": "tx_length",
    "transcript length": "tx_length",
    "tx length": "tx_length",
    "tx_length": "tx_length",
    "txlen": "tx_length",
    # tx_support
    "transcript support level (tsl)": "tx_support",
    "transcript support level": "tx_support",
    "tx support": "tx_support",
    "tx_support": "tx_support",
    "tsl": "tx_support",
}


def normalize_header(name):
    """Lowercase and collapse internal whitespace for alias matching."""
    return " ".join(str(name).strip().lower().split())


def resolve_columns(headers, renames):
    """Match input headers to canonical keys; raise if a required key is unmapped.

    ``renames`` maps an exact input header to a canonical key (from --rename) and
    takes precedence over the built-in alias table. The first input column
    matching a key wins.
    """
    found = {}
    for hdr in headers:
        key = renames.get(hdr) or COLUMN_ALIASES.get(normalize_header(hdr))
        if key and key not in found:
            found[key] = hdr
    missing = REQUIRED_KEYS - set(found)
    if missing:
        raise SystemExit(
            "Could not find required column(s) in the input: "
            + ", ".join(sorted(missing))
            + ".\nInput headers were: "
            + ", ".join(map(str, headers))
            + '.\nUse --rename "<input header>=<key>" to map them '
            + "(keys: "
            + ", ".join(CANONICAL_KEYS)
            + ")."
        )
    return found


def build_gene_info(in_fname, renames=None, chromosomes=None):
    """Normalize a BioMart export into a canonical CNVkit gene-info DataFrame.

    Columns are selected by name, reordered, lightly validated, and sorted by
    genomic position. Returns a DataFrame whose columns are the canonical short
    keys in canonical order.
    """
    # ``comment="#"`` skips a leading provenance comment line, so a CNVkit-format
    # gene-info file can be re-normalized as well as a raw BioMart export.
    raw = pd.read_csv(in_fname, sep="\t", dtype=str, na_filter=False, comment="#")
    found = resolve_columns(list(raw.columns), renames or {})
    n_input = len(raw)

    out = pd.DataFrame()
    for key in CANONICAL_KEYS:
        if key in found:
            out[key] = raw[found[key]].astype(str).str.strip()
        else:
            out[key] = DEFAULTS[key]
            logging.info(
                "Column %r not found; filling with default %r", key, DEFAULTS[key]
            )

    # gene_id is the table's key; rows without one are unusable.
    out = out[out["gene_id"] != ""]

    # Required numeric columns must parse, or the RNA pipeline crashes on load
    # (gc -> float/100) or on depth normalization (divide by tx_length).
    gc = pd.to_numeric(out["gc"], errors="coerce")
    start = pd.to_numeric(out["start"], errors="coerce")
    end = pd.to_numeric(out["end"], errors="coerce")
    txlen = pd.to_numeric(out["tx_length"], errors="coerce")
    valid = gc.notna() & start.notna() & end.notna() & txlen.notna()
    if (~valid).any():
        logging.warning(
            "Dropping %d genes with non-numeric gc/start/end/tx_length",
            int((~valid).sum()),
        )
        out, gc, start, end, txlen = (s[valid] for s in (out, gc, start, end, txlen))

    # Transcript length must be positive; align_gene_info_to_samples filters
    # these too, but a clean asset should not ship them.
    positive = txlen > 0
    if (~positive).any():
        logging.warning(
            "Dropping %d genes with transcript length <= 0", int((~positive).sum())
        )
        out, gc, start, end = (s[positive] for s in (out, gc, start, end))

    # Sanity warnings (non-fatal): catch a fraction-vs-percentage mix-up and
    # reversed coordinates without silently mangling the data.
    if ((gc < 0) | (gc > 100)).any():
        logging.warning(
            "%d GC values are outside [0, 100]; the GC column should be a "
            "percentage, not a fraction",
            int(((gc < 0) | (gc > 100)).sum()),
        )
    if (start > end).any():
        logging.warning("%d genes have start > end", int((start > end).sum()))
    n_dupes = int(out["gene_id"].duplicated().sum())
    if n_dupes:
        logging.warning(
            "%d duplicate gene IDs present (cnvlib.rna deduplicates these on load)",
            n_dupes,
        )

    if chromosomes:
        keep = {c.strip() for c in chromosomes}
        before = len(out)
        out = out[out["chromosome"].isin(keep)]
        logging.info(
            "Kept %d/%d genes on chromosomes: %s",
            len(out),
            before,
            ", ".join(sorted(keep)),
        )

    # Sort by genomic position, using CNVkit's canonical chromosome ordering
    # (1..22 < X < Y < MT < scaffolds), for human-readable, stable output.
    sort_idx = (
        out.assign(
            _chrkey=out["chromosome"].map(sorter_chrom),
            _start=pd.to_numeric(out["start"], errors="coerce"),
        )
        .sort_values(["_chrkey", "_start"])
        .index
    )
    out = out.loc[sort_idx]

    logging.info("Built gene info for %d genes (from %d input rows)", len(out), n_input)
    return out


def write_gene_info(out_df, out_fname, genome=None, source=None):
    """Write the gene-info table with a provenance comment + BioMart header.

    The leading comment line is required by ``cnvlib.rna.load_gene_info``, which
    reads with ``header=1`` (i.e. it skips the first line before the column
    header). Emitting it here keeps every gene; a file whose first line is the
    column header instead would lose its first data row on load.
    """
    n = len(out_df)
    comment = f"# CNVkit gene info | genome={genome or 'unknown'}"
    if source:
        comment += f" | source={os.path.basename(source)}"
    comment += f" | {n} genes"

    display = out_df.copy()
    display.columns = [disp for _key, disp in CANONICAL_COLUMNS]

    to_stdout = out_fname in (None, "-")
    ctx = nullcontext(sys.stdout) if to_stdout else open(out_fname, "w")
    with ctx as handle:
        handle.write(comment + "\n")
        display.to_csv(handle, sep="\t", index=False)
    if not to_stdout:
        logging.info("Wrote %s (%d genes)", out_fname, n)


def parse_renames(items):
    """Parse repeated ``--rename HEADER=KEY`` arguments into a dict."""
    renames = {}
    valid_keys = set(CANONICAL_KEYS)
    for item in items:
        if "=" not in item:
            raise SystemExit(f"--rename expects HEADER=KEY, got: {item!r}")
        hdr, key = item.split("=", 1)
        key = key.strip()
        if key not in valid_keys:
            raise SystemExit(
                f"--rename target {key!r} is not a canonical key; choose from: "
                + ", ".join(CANONICAL_KEYS)
            )
        renames[hdr.strip()] = key
    return renames


def argument_parsing():
    AP = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    AP.add_argument(
        "biomart_tsv",
        metavar="BIOMART_TSV",
        help="Ensembl BioMart TSV export (tab-delimited, with a header row).",
    )
    AP.add_argument(
        "-o",
        "--output",
        metavar="FILE",
        help="Output gene-info table (default: stdout).",
    )
    AP.add_argument(
        "-g",
        "--genome",
        help="Reference genome label recorded in the output header comment "
        "(e.g. hg19, mm10, mm39).",
    )
    AP.add_argument(
        "--rename",
        metavar="HEADER=KEY",
        action="append",
        default=[],
        help="Map an unrecognized input column header to a canonical key "
        "(" + ", ".join(CANONICAL_KEYS) + "). Repeatable.",
    )
    AP.add_argument(
        "--chromosomes",
        help="Comma-separated chromosome/scaffold names to keep "
        "(default: keep all). Example: 1,2,3,...,22,X,Y,MT",
    )
    return AP.parse_args()


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    args = argument_parsing()
    renames = parse_renames(args.rename)
    chromosomes = args.chromosomes.split(",") if args.chromosomes else None
    out_df = build_gene_info(args.biomart_tsv, renames=renames, chromosomes=chromosomes)
    write_gene_info(out_df, args.output, genome=args.genome, source=args.biomart_tsv)


if __name__ == "__main__":
    main()
