"""NCBI data access — cache-first, live fallback.

Reads seeded fixtures under `data/cache/` first (so the engine + tests run fully
offline for GAPDH/MYC); on a miss, calls NCBI Datasets v2 / E-utilities live and
persists the result in the same cache shape.

Never logs the API key. See vault `02 Data/NCBI Data Source and API.md`.
"""

from __future__ import annotations

import json
import os
import re
from pathlib import Path
from urllib.parse import urlencode

CACHE_DIR = Path(os.environ.get("TMJ_CACHE_DIR", Path(__file__).resolve().parent.parent / "data" / "cache"))
GRCH38_PREFIX = "NC_0000"
API_KEY = os.environ.get("NCBI_API_KEY", "")
DATASETS = "https://api.ncbi.nlm.nih.gov/datasets/v2"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

Interval = tuple[int, int]


class NotFound(Exception):
    pass


def _read_json(p: Path):
    return json.loads(p.read_text()) if p.exists() else None


def _write_json(p: Path, obj) -> None:
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(obj))


def _base_accession(acc: str) -> str:
    return acc.split(".")[0].strip().upper()


# ---------------------------------------------------------------- live helpers

def _http_get_json(url: str, headers: dict | None = None) -> dict:
    import httpx
    h = {"Accept": "application/json"}
    if API_KEY:
        h["api-key"] = API_KEY
    if headers:
        h.update(headers)
    with httpx.Client(timeout=30) as c:
        r = c.get(url, headers=h)
        r.raise_for_status()
        return r.json()


def _http_get_text(url: str) -> str:
    import httpx
    with httpx.Client(timeout=30) as c:
        r = c.get(url)
        r.raise_for_status()
        return r.text


# ---------------------------------------------------------------- public API

def resolve_accession(accession: str) -> tuple[str, str]:
    """accession -> (gene_symbol, gene_id). Cache-first."""
    base = _base_accession(accession)
    cached = _read_json(CACHE_DIR / "resolve" / f"{base}.json")
    if cached:
        return cached["symbol"], cached.get("gene_id", "")

    data = _http_get_json(f"{DATASETS}/gene/accession/{base}")
    reports = data.get("reports") or data.get("gene", {}).get("reports") or []
    if not reports:
        raise NotFound(f"No gene found for accession {accession}")
    gene = reports[0].get("gene") or reports[0].get("product") or reports[0]
    symbol = str(gene.get("symbol", ""))
    gene_id = str(gene.get("gene_id", ""))
    if not symbol:
        raise NotFound(f"Could not resolve symbol for {accession}")
    _write_json(CACHE_DIR / "resolve" / f"{base}.json",
                {"symbol": symbol, "gene_id": gene_id, "resolved_accession": accession})
    return symbol, gene_id


def get_product_report(symbol: str) -> dict:
    """product_report payload for a gene symbol (human). Cache-first."""
    cached = _read_json(CACHE_DIR / "product_report" / f"{symbol}__human.json")
    if cached:
        return cached
    data = _http_get_json(f"{DATASETS}/gene/symbol/{symbol}/taxon/human/product_report")
    _write_json(CACHE_DIR / "product_report" / f"{symbol}__human.json", data)
    return data


def get_sequence(accession: str) -> str:
    """mRNA nucleotide sequence (5'->3', uppercase) for an accession. Cache-first."""
    cached = _read_json(CACHE_DIR / "sequence" / f"{accession}.json")
    if cached:
        return cached["seq"].upper()
    fasta = _http_get_text(
        f"{EUTILS}/efetch.fcgi?db=nuccore&id={accession}&rettype=fasta&retmode=text"
        + (f"&api_key={API_KEY}" if API_KEY else "")
    )
    lines = fasta.splitlines()
    seq = "".join(l.strip() for l in lines if not l.startswith(">")).upper()
    # The defline names the isoform ("... (GAPDH), transcript variant 1, mRNA"). Kept because
    # it is already in hand: it is the fallback when the gene's product report omits `name`,
    # and re-fetching a whole record later just to read its title would be absurd.
    title = next((l[1:].strip() for l in lines if l.startswith(">")), "")
    _write_json(CACHE_DIR / "sequence" / f"{accession}.json",
                {"accession": accession, "seq": seq, "title": title})
    return seq


_VARIANT_IN_TITLE = re.compile(r"\btranscript variant\s+([^,;]+)", re.IGNORECASE)


def variant_from_title(accession: str) -> str | None:
    """NCBI's variant designation read off the sequence record's defline, or None.

    Second source for the same fact the product report's `name` carries, for the case where
    that field is absent. Cache-only: a record whose title was never stored is not worth a
    network round trip to label a row, and the title arrives free the next time the sequence
    itself is fetched.
    """
    cached = _read_json(CACHE_DIR / "sequence" / f"{accession}.json") or {}
    m = _VARIANT_IN_TITLE.search(cached.get("title") or "")
    return f"transcript variant {m.group(1).strip()}" if m else None


# ---------------------------------------------------------------- parsing

def grch38_exons(transcript: dict) -> list[Interval] | None:
    """GRCh38 (begin,end) exon intervals, 1-based inclusive, or None if unavailable."""
    for loc in transcript.get("genomic_locations") or []:
        if not str(loc.get("genomic_accession_version", "")).startswith(GRCH38_PREFIX):
            continue
        out: list[Interval] = []
        for e in loc.get("exons") or []:
            try:
                b, t = int(e["begin"]), int(e["end"])
            except (KeyError, TypeError, ValueError):
                continue
            out.append((min(b, t), max(b, t)))
        return out or None
    return None


def grch38_strand(transcript: dict) -> str:
    """Genomic orientation of the transcript on GRCh38: "+", "-", or "" if unstated.

    Taken from NCBI's explicit `orientation` rather than inferred from exon order, so a
    single-exon transcript — where there is no exon order to read a direction from — still
    reports its strand.
    """
    for loc in transcript.get("genomic_locations") or []:
        if not str(loc.get("genomic_accession_version", "")).startswith(GRCH38_PREFIX):
            continue
        o = str((loc.get("genomic_range") or {}).get("orientation", "")).lower()
        if not o:
            ex = loc.get("exons") or []
            o = str(ex[0].get("orientation", "")).lower() if ex else ""
        return {"plus": "+", "minus": "-"}.get(o, "")
    return ""


def _chromosome_from_acc(gav: str) -> str:
    """GRCh38 genomic accession -> chromosome label. NC_000012 -> '12', 23 -> 'X', 24 -> 'Y'."""
    m = re.match(r"NC_0*(\d+)", gav or "")
    if not m:
        return ""
    n = int(m.group(1))
    if 1 <= n <= 22:
        return str(n)
    return {23: "X", 24: "Y", 12920: "MT"}.get(n, "")


def _grch38_chromosome(transcript: dict) -> str:
    for loc in transcript.get("genomic_locations") or []:
        gav = str(loc.get("genomic_accession_version", ""))
        if gav.startswith(GRCH38_PREFIX):
            return _chromosome_from_acc(gav)
    return ""


def nm_transcripts(product_report: dict) -> tuple[str, str, str, str, str, list[dict]]:
    """Return (gene_id, symbol, description, chromosome, strand, [ {accession, variant,
    is_mane, exons, strand} ... ]) for NM only. Gene strand = the first transcript's — every NM of a
    gene is transcribed from the same strand."""
    reports = product_report.get("reports") or []
    if not reports:
        raise NotFound("Empty product report")
    product = reports[0].get("product") or {}
    gene_id = str(product.get("gene_id", ""))
    symbol = str(product.get("symbol", ""))
    description = str(product.get("description", ""))
    chromosome = ""
    strand = ""
    out: list[dict] = []
    for t in product.get("transcripts") or []:
        acc = (t.get("accession_version") or "").strip()
        if not acc.startswith("NM_"):
            continue
        exons = grch38_exons(t)
        if not exons:
            continue
        if not chromosome:
            chromosome = _grch38_chromosome(t)
        tx_strand = grch38_strand(t)
        if not strand:
            strand = tx_strand
        out.append({
            "accession": acc,
            # NCBI's own designation for the isoform, e.g. "transcript variant 5". Absent on
            # a gene with a single NM — there is no variant to number — which is what the
            # clients render as "mono-isoform".
            "variant": (t.get("name") or "").strip() or None,
            "is_mane": t.get("select_category") == "MANE_SELECT",
            "exons": exons,
            "strand": tx_strand,
            "cds": _cds_range(t),   # (begin, end) in 1-based transcript coords, or None
        })
    return gene_id, symbol, description, chromosome, strand, out


def _cds_range(transcript: dict) -> tuple[int, int] | None:
    """CDS start/end in 1-based transcript (mRNA) coordinates, or None."""
    cds = transcript.get("cds") or {}
    rng = cds.get("range") or []
    if not rng:
        return None
    try:
        b, e = int(rng[0]["begin"]), int(rng[0]["end"])
    except (KeyError, TypeError, ValueError):
        return None
    return (min(b, e), max(b, e))


ACCESSION_RE = re.compile(r"^NM_\d+(\.\d+)?$", re.IGNORECASE)


def is_valid_nm(accession: str) -> bool:
    return bool(ACCESSION_RE.match(accession.strip()))


# ---------------------------------------------------------------- typeahead

_SYMBOL_BEFORE_COMMA = re.compile(r"\(([^)]+)\)\s*,")
_SYMBOL_ANY = re.compile(r"\(([A-Za-z0-9_./-]{2,})\)")


def _symbol_from_title(title: str) -> str:
    """Pull the gene symbol out of an nuccore title, e.g.
    'Homo sapiens ... (GAPDH), transcript variant 3, mRNA' -> 'GAPDH'."""
    m = _SYMBOL_BEFORE_COMMA.search(title or "")
    if m:
        return m.group(1)
    m2 = _SYMBOL_ANY.search(title or "")
    return m2.group(1) if m2 else ""


def _eutils_json(endpoint: str, params: dict) -> dict:
    if API_KEY:
        params = {**params, "api_key": API_KEY}
    return _http_get_json(f"{EUTILS}/{endpoint}?{urlencode(params)}")


INDEX_PATH = Path(os.environ.get("TMJ_NM_INDEX",
                                 Path(__file__).resolve().parent.parent / "data" / "nm_index.json"))
_NM_INDEX: list[dict] | None = None
_NM_KEYS: list[str] | None = None


def _load_index() -> tuple[list[dict], list[str]]:
    """Lazy-load the human NM accession index (accession -> gene). Sorted for bisect."""
    global _NM_INDEX, _NM_KEYS
    if _NM_INDEX is None:
        rows = _read_json(INDEX_PATH) or []
        rows.sort(key=lambda r: r["accession"].upper())
        _NM_INDEX = rows
        _NM_KEYS = [r["accession"].upper() for r in rows]
    return _NM_INDEX, _NM_KEYS or []


def suggest_accessions(q: str, limit: int = 8) -> list[dict]:
    """Typeahead: human NM mRNAs whose accession starts with `q`.

    Primary source is a local index of ~70k real human NM accessions (built from NCBI
    data). Reliable, instant, comprehensive. Falls back to a live NCBI lookup only if the
    index has no match (e.g. an accession newer than the index). Returns [{accession, gene}].
    See vault `02 Data/NCBI Data Source and API.md`.
    """
    base = q.strip().upper().split(".")[0]
    if not base.startswith("NM_") or len(base) < 5:
        return []

    from bisect import bisect_left
    rows, keys = _load_index()
    out: list[dict] = []
    i = bisect_left(keys, base)
    while i < len(keys) and keys[i].startswith(base):
        out.append(rows[i])
        if len(out) >= limit:
            break
        i += 1
    if out:
        return out
    return _suggest_live(base, limit)


def suggest_genes(q: str, limit: int = 8) -> list[dict]:
    """Typeahead: human protein-coding genes whose official symbol starts with `q`.

    Live NCBI E-utilities (db=gene), per-prefix cached. The protein-coding filter drops
    pseudogenes / antisense RNAs that would have no NM transcript to analyze. Best-effort —
    any error (incl. E-utilities rate limiting without an API key) yields []. Non-empty
    results only are cached, so a transient failure doesn't poison the prefix. Returns
    [{symbol, description}]."""
    base = q.strip().upper()
    if len(base) < 2:
        return []
    cached = _read_json(CACHE_DIR / "suggest_gene" / f"{base}.json")
    if cached is not None:
        return cached[:limit]
    try:
        term = (f"{base}*[Preferred Symbol] AND Homo sapiens[Organism] "
                'AND alive[prop] AND "genetype protein coding"[Properties]')
        es = _eutils_json("esearch.fcgi",
                          {"db": "gene", "term": term, "retmax": limit,
                           "retmode": "json", "sort": "relevance"})
        ids = ((es.get("esearchresult") or {}).get("idlist")) or []
        out: list[dict] = []
        if ids:
            summ = _eutils_json("esummary.fcgi",
                                {"db": "gene", "id": ",".join(ids), "retmode": "json"})
            res = summ.get("result") or {}
            for uid in res.get("uids", []):
                r = res.get(uid) or {}
                sym = r.get("name") or ""
                if sym.upper().startswith(base):   # drop Entrez wildcard noise
                    out.append({"symbol": sym, "description": r.get("description", "")})
        if out:
            _write_json(CACHE_DIR / "suggest_gene" / f"{base}.json", out[:limit])
        return out[:limit]
    except Exception:
        return []


def _suggest_live(base: str, limit: int) -> list[dict]:
    """Fallback: NCBI E-utilities lookup (per-prefix cached). Post-filtered to true prefix
    matches (Entrez's accession wildcard is noisy). Best-effort — errors yield []."""
    cached = _read_json(CACHE_DIR / "suggest" / f"{base}.json")
    if cached is not None:
        return cached[:limit]
    try:
        term = f"{base}*[ACCN] AND biomol_mrna[PROP] AND srcdb_refseq[PROP] AND txid9606[ORGN]"
        es = _eutils_json("esearch.fcgi",
                          {"db": "nuccore", "term": term, "retmax": 20, "retmode": "json"})
        ids = ((es.get("esearchresult") or {}).get("idlist")) or []
        out: list[dict] = []
        if ids:
            summ = _eutils_json("esummary.fcgi",
                                {"db": "nuccore", "id": ",".join(ids), "retmode": "json"})
            res = summ.get("result") or {}
            for uid in res.get("uids", []):
                r = res.get(uid) or {}
                acc = r.get("accessionversion") or r.get("caption") or ""
                if acc.upper().startswith(base):   # drop Entrez wildcard noise
                    out.append({"accession": acc, "gene": _symbol_from_title(r.get("title", ""))})
        _write_json(CACHE_DIR / "suggest" / f"{base}.json", out[:limit])
        return out[:limit]
    except Exception:
        return []
