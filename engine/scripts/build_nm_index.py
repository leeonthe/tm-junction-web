#!/usr/bin/env python3
"""Build the human NM accession typeahead index (`data/nm_index.json`).

Reads NCBI Datasets v2 `product_report` JSON payloads from a cache directory and
extracts every curated `NM_` transcript accession with its gene symbol. The result
powers the `/suggest` endpoint (see app/ncbi.py `suggest_accessions`).

Usage:
  python scripts/build_nm_index.py --cache-dir /path/to/product_report/cache
      [--out data/nm_index.json]

Any directory of `*.json` product-report payloads works (e.g. the sibling
`ncbi exon/cache`). Re-run to refresh when RefSeq updates.
"""
from __future__ import annotations

import argparse
import glob
import json
import os


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache-dir", required=True, help="dir of product_report *.json files")
    ap.add_argument("--out", default=os.path.join(os.path.dirname(__file__), "..", "data", "nm_index.json"))
    args = ap.parse_args()

    index: dict[str, str] = {}   # accession -> gene symbol
    genes = 0
    for f in glob.glob(os.path.join(args.cache_dir, "*.json")):
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
            if acc.startswith("NM_"):
                index[acc] = symbol
        genes += 1

    rows = sorted(({"accession": a, "gene": g} for a, g in index.items()),
                  key=lambda r: r["accession"])
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    json.dump(rows, open(args.out, "w"))
    print(f"genes read: {genes} | NM accessions: {len(rows)} | wrote {args.out}")


if __name__ == "__main__":
    main()
