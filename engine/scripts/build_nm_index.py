#!/usr/bin/env python3
"""Build a species' NM/NR accession typeahead index (`data/index/<species>.json`).

Every curated `NM_` (mRNA) and `NR_` (non-coding RNA) transcript accession of the species
with its gene symbol — the two classes the engine analyzes (app/ncbi.py REFSEQ_PREFIXES).
The result powers the `/suggest` endpoint (see app/ncbi.py `suggest_accessions`).

Two sources:

  --from-ncbi            NCBI Datasets v2, the species' whole product report as a
                         3-column TSV (gene id, symbol, transcript accession), 1000 genes a
                         page. No local data needed; a few minutes per species.
  --cache-dir DIR        a directory of cached `product_report` JSON payloads (e.g. the
                         sibling `ncbi exon/cache`, human only).

Usage:
  python scripts/build_nm_index.py --species mouse --from-ncbi
  python scripts/build_nm_index.py --species human --cache-dir "/path/to/ncbi exon/cache"
  python scripts/build_nm_index.py --all --from-ncbi          # every species but human

The file is two parallel arrays, sorted by accession for bisect:
  {"species": "mouse", "accessions": ["NM_000…", …], "genes": ["Gapdh", …]}
Re-run to refresh when RefSeq updates.
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from app import species as species_mod  # noqa: E402

DATASETS = "https://api.ncbi.nlm.nih.gov/datasets/v2"
PREFIXES = ("NM_", "NR_")
# Everything that can carry a curated transcript. Biological regions and tRNAs cannot, and
# are most of what a taxon-wide listing would otherwise page through (40k regions in mouse).
GENE_TYPES = ("PROTEIN_CODING", "ncRNA", "snoRNA", "snRNA", "scRNA", "rRNA", "miscRNA",
              "PSEUDO", "OTHER")


def from_cache(cache_dir: str, slug: str) -> dict[str, str]:
    index: dict[str, str] = {}
    for f in glob.glob(os.path.join(cache_dir, f"*__{slug}.json")):
        try:
            d = json.load(open(f))
        except Exception:
            continue
        reports = d.get("reports") or []
        if not reports:
            continue
        product = reports[0].get("product") or {}
        symbol = product.get("symbol", "")
        for t in product.get("transcripts") or []:
            acc = (t.get("accession_version") or "").strip()
            if acc.startswith(PREFIXES):
                index[acc] = symbol
    return index


def _get(url: str, api_key: str) -> tuple[str, str]:
    """One paced GET -> (body, next page token). Retries 429/5xx with backoff."""
    headers = {"Accept": "text/tab-separated-values"}
    if api_key:
        headers["api-key"] = api_key
    for attempt in range(6):
        time.sleep(0.12 if api_key else 0.4)
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=headers), timeout=120) as r:
                return r.read().decode(), r.headers.get("x-ncbi-next-page-token", "") or ""
        except urllib.error.HTTPError as e:
            if e.code != 429 and e.code < 500:
                raise
        except (urllib.error.URLError, TimeoutError):
            pass
        time.sleep(2 ** attempt)
    raise RuntimeError(f"NCBI kept failing: {url}")


def from_ncbi(tax_id: str, api_key: str = "") -> dict[str, str]:
    index: dict[str, str] = {}
    params = [("page_size", "1000"), ("table_fields", "gene-id"), ("table_fields", "symbol"),
              ("table_fields", "transcript-accession")] + [("types", t) for t in GENE_TYPES]
    token, page = "", 0
    while True:
        q = params + ([("page_token", token)] if token else [])
        body, token = _get(f"{DATASETS}/gene/taxon/{tax_id}/product_report?"
                           + urllib.parse.urlencode(q), api_key)
        page += 1
        rows = body.splitlines()[1:]                       # header first
        for row in rows:
            cols = row.split("\t")
            if len(cols) >= 3 and cols[2].startswith(PREFIXES):
                index[cols[2].strip()] = cols[1].strip()
        print(f"  page {page}: {len(rows)} rows, {len(index)} accessions so far", flush=True)
        if not token or not rows:
            return index


def write(slug: str, index: dict[str, str], out_dir: str) -> None:
    accessions = sorted(index, key=str.upper)
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, f"{slug}.json")
    json.dump({"species": slug, "accessions": accessions,
               "genes": [index[a] for a in accessions]},
              open(path, "w"), separators=(",", ":"))
    nr = sum(a.startswith("NR_") for a in accessions)
    print(f"{slug}: NM {len(accessions) - nr} | NR {nr} | genes {len(set(index.values()))} "
          f"| wrote {path}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", help="species slug (see app/species.py)")
    ap.add_argument("--all", action="store_true", help="every species except human")
    ap.add_argument("--from-ncbi", action="store_true", help="page the species from NCBI Datasets")
    ap.add_argument("--cache-dir", help="dir of product_report *__<species>.json files")
    ap.add_argument("--out-dir", default=os.path.join(os.path.dirname(__file__), "..", "data", "index"))
    args = ap.parse_args()
    if not (args.all or args.species) or not (args.from_ncbi or args.cache_dir):
        ap.error("give --species or --all, and --from-ncbi or --cache-dir")

    todo = ([s for s in species_mod.SPECIES if s is not species_mod.HUMAN] if args.all
            else [species_mod.get(args.species)])
    for sp in todo:
        print(f"{sp.common} ({sp.scientific}, taxid {sp.tax_id})", flush=True)
        index = (from_ncbi(sp.tax_id, os.environ.get("NCBI_API_KEY", "")) if args.from_ncbi
                 else from_cache(args.cache_dir, sp.slug))
        if not index:
            sys.exit(f"{sp.slug}: nothing found — refusing to write an empty index")
        write(sp.slug, index, args.out_dir)


if __name__ == "__main__":
    main()
