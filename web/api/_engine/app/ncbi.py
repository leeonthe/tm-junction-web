"""NCBI data access — cache-first, live fallback.

Reads seeded fixtures under `data/cache/` first (so the engine + tests run fully
offline for the fixture genes); on a miss, calls NCBI Datasets v2 / E-utilities live and
persists the result in the same cache shape.

Species: every query is scoped by the taxonomy id of a `species.Species`, human unless the
caller says otherwise — see species.py for why ids and not names.

Never logs the API key. See vault `02 Data/NCBI Data Source and API.md`.
"""

from __future__ import annotations

import json
import os
import re
from dataclasses import dataclass, field
from pathlib import Path
from urllib.parse import quote, urlencode

from . import species as species_mod
from .species import HUMAN, Species

CACHE_DIR = Path(os.environ.get("TMJ_CACHE_DIR", Path(__file__).resolve().parent.parent / "data" / "cache"))
API_KEY = os.environ.get("NCBI_API_KEY", "")
DATASETS = "https://api.ncbi.nlm.nih.gov/datasets/v2"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

Interval = tuple[int, int]


class NotFound(Exception):
    pass


def _read_json(p: Path):
    return json.loads(p.read_text()) if p.exists() else None


def _write_json(p: Path, obj) -> None:
    # Best-effort: on a read-only filesystem (a serverless bundle) the cache is a seed to
    # read from, and a failed write must cost a re-fetch later, not the request now.
    try:
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(json.dumps(obj))
    except OSError:
        pass


def _base_accession(acc: str) -> str:
    return acc.split(".")[0].strip().upper()


# ---------------------------------------------------------------- live helpers

class RateLimited(Exception):
    """NCBI said 429 and kept saying it through every retry.

    Kept distinct from an ordinary HTTP failure because the answer to the caller is
    different: nothing is wrong with the request, the shared egress IP is simply over
    NCBI's per-IP budget this second — retrying shortly will work. main.py maps this to
    a 503 with that message instead of a 500 with a stack trace, which is what a single
    unretried 429 used to become.
    """


# NCBI allows ~3 requests/second per IP without an API key (10/s with one). An analyze of
# an uncached gene fires one fetch per isoform back to back, and on serverless hosting the
# egress IP is shared with strangers — so the engine paces itself rather than betting the
# whole request on never being the one over the line.
_MIN_INTERVAL = 0.11 if API_KEY else 0.35
_last_call = 0.0


def _paced_get(url: str, headers: dict | None = None):
    """One NCBI GET: paced under the per-IP budget, retried with backoff on 429/5xx."""
    import time

    import httpx

    global _last_call
    delays = (0.5, 1.0, 2.0)                    # after attempts 1..3; 4 attempts total
    for attempt in range(len(delays) + 1):
        wait = _MIN_INTERVAL - (time.monotonic() - _last_call)
        if wait > 0:
            time.sleep(wait)
        _last_call = time.monotonic()
        with httpx.Client(timeout=30) as c:
            r = c.get(url, headers=headers)
        if r.status_code == 429 or r.status_code >= 500:
            if attempt < len(delays):
                # Honor Retry-After when NCBI names a delay; else back off blind.
                try:
                    retry_after = float(r.headers.get("retry-after", ""))
                except ValueError:
                    retry_after = 0.0
                time.sleep(max(delays[attempt], min(retry_after, 10.0)))
                continue
            if r.status_code == 429:
                raise RateLimited(url)
        r.raise_for_status()
        return r
    raise RateLimited(url)                      # unreachable; keeps type-checkers honest


def _http_get_json(url: str, headers: dict | None = None) -> dict:
    h = {"Accept": "application/json"}
    if API_KEY:
        h["api-key"] = API_KEY
    if headers:
        h.update(headers)
    return _paced_get(url, h).json()


def _http_get_text(url: str) -> str:
    return _paced_get(url).text


# ---------------------------------------------------------------- public API

def resolve_accession(accession: str, refresh: bool = False) -> tuple[str, str, str]:
    """accession -> (gene_symbol, gene_id, tax_id). Cache-first.

    The accession alone names the species — a RefSeq transcript belongs to exactly one
    organism — so an accession search needs no species from the user. tax_id is "" only for
    an entry cached before species support; those are taken as human, and `refresh`
    re-resolves one that turns out not to be (see analyze._gene_of).
    """
    base = _base_accession(accession)
    cached = None if refresh else _read_json(CACHE_DIR / "resolve" / f"{base}.json")
    if cached:
        return cached["symbol"], cached.get("gene_id", ""), str(cached.get("tax_id", ""))

    data = _http_get_json(f"{DATASETS}/gene/accession/{base}")
    reports = data.get("reports") or data.get("gene", {}).get("reports") or []
    if reports:
        gene = reports[0].get("gene") or reports[0].get("product") or reports[0]
    else:
        # Datasets indexes accessions through the genome annotation, so a record curated
        # since the last annotation run (NR_201105.1, HTRA1-AS1) resolves to nothing there.
        # The record itself still names its gene.
        rec = record_structures([accession.strip().upper()]).get(accession.strip().upper()) or {}
        if not rec.get("gene"):
            raise NotFound(f"No gene found for accession {accession}")
        gene = {"symbol": rec["gene"], "gene_id": rec.get("gene_id", ""),
                "tax_id": rec.get("tax_id", "")}
    symbol = str(gene.get("symbol", ""))
    gene_id = str(gene.get("gene_id", ""))
    tax_id = str(gene.get("tax_id", ""))
    if not symbol:
        raise NotFound(f"Could not resolve symbol for {accession}")
    _write_json(CACHE_DIR / "resolve" / f"{base}.json",
                {"symbol": symbol, "gene_id": gene_id, "tax_id": tax_id,
                 "resolved_accession": accession})
    return symbol, gene_id, tax_id


def _file_safe(name: str) -> str:
    """A gene symbol as a file name. Human symbols pass through unchanged; fly's carry
    punctuation — l(2)gl, mt:CoI — and a symbol is user input, so it never reaches a path
    or a URL unescaped."""
    return quote(name, safe="")


def _product_report_path(symbol: str, sp: Species) -> Path:
    return CACHE_DIR / "product_report" / f"{_file_safe(symbol)}__{sp.slug}.json"


def get_product_report(symbol: str, sp: Species = HUMAN) -> dict:
    """product_report payload for a gene symbol in a species. Cache-first.

    NCBI matches symbols case-insensitively and answers with the official spelling, so the
    report is cached under THAT, and looked up under the capitalizations a person would
    type: a mouse "gapdh" finds the cached "Gapdh" instead of costing a request.
    """
    sym = symbol.strip()
    for spelling in dict.fromkeys([sym, sym.upper(), sym.capitalize(), sym.lower()]):
        cached = _read_json(_product_report_path(spelling, sp))
        if cached and cached.get("reports"):
            return cached
    data = _http_get_json(
        f"{DATASETS}/gene/symbol/{_file_safe(sym)}/taxon/{sp.tax_id}/product_report")
    reports = data.get("reports") or []
    if reports:
        # More than one gene can answer to a symbol; the exact spelling wins, then the
        # case-insensitive one. Everything downstream reads reports[0].
        def rank(r: dict) -> int:
            got = str((r.get("product") or {}).get("symbol", ""))
            return 0 if got == sym else 1 if got.lower() == sym.lower() else 2
        reports.sort(key=rank)
        official = str((reports[0].get("product") or {}).get("symbol", "")) or sym
        _write_json(_product_report_path(official, sp), data)
    return data


def get_gene_report(gene_id: str) -> dict:
    """The gene-level Datasets record — chromosome and annotated assemblies. Cache-first,
    and best-effort: it only ever supplies labels (see _labels_from_gene_report)."""
    if not gene_id:
        return {}
    cached = _read_json(CACHE_DIR / "gene_report" / f"{gene_id}.json")
    if cached:
        return cached
    try:
        data = _http_get_json(f"{DATASETS}/gene/id/{gene_id}")
    except RateLimited:
        raise
    except Exception:
        return {}
    gene = ((data.get("reports") or [{}])[0].get("gene")) or {}
    if gene:
        _write_json(CACHE_DIR / "gene_report" / f"{gene_id}.json", gene)
    return gene


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


def get_sequences(accessions: list[str]) -> dict[str, str]:
    """mRNA sequences for many accessions — cache-first, ONE efetch for all the misses.

    The per-accession loop this replaces was the engine's biggest NCBI spender: an
    uncached 25-isoform gene (NRXN1) made 25 serial efetch calls, and under NCBI's
    per-IP request budget — shared, on serverless hosting, with strangers — that is how
    a plain gene search turns into 429s. efetch takes a comma-joined id list, so all
    the misses travel in a single request; the multi-FASTA comes back split by header,
    each record keyed by the accession its defline starts with. Titles are kept, as in
    get_sequence, because they are the fallback source for variant designations.
    """
    out: dict[str, str] = {}
    missing: list[str] = []
    for a in accessions:
        cached = _read_json(CACHE_DIR / "sequence" / f"{a}.json")
        if cached:
            out[a] = cached["seq"].upper()
        else:
            missing.append(a)
    for i in range(0, len(missing), 100):        # eutils is comfortable at 100 ids/call
        chunk = missing[i:i + 100]
        fasta = _http_get_text(
            f"{EUTILS}/efetch.fcgi?db=nuccore&id={','.join(chunk)}&rettype=fasta&retmode=text"
            + (f"&api_key={API_KEY}" if API_KEY else ""))
        for rec in fasta.split("\n>"):
            rec = rec.lstrip(">")
            if not rec.strip():
                continue
            header, _, body = rec.partition("\n")
            acc = header.split()[0].strip()
            seq = "".join(body.split()).upper()
            if not seq:
                continue
            out[acc] = seq
            _write_json(CACHE_DIR / "sequence" / f"{acc}.json",
                        {"accession": acc, "seq": seq, "title": header.strip()})
    # Anything efetch did not return (retired id, odd defline) still gets its own call —
    # correctness beats saving one request when the batch and the answer disagree.
    for a in accessions:
        if a not in out:
            out[a] = get_sequence(a)
    return out


_VARIANT_IN_TITLE = re.compile(r"\btranscript variant\s+([^,;]+)", re.IGNORECASE)


def variant_in(title: str | None) -> str | None:
    """"...(GAPDH), transcript variant 1, mRNA" -> "transcript variant 1".

    FlyBase-derived titles can put the gene symbol AFTER the designation —
    "…dehydrogenase 1, transcript variant A (Gapdh1), mRNA" — so a trailing parenthetical
    is the symbol, not part of the variant's name."""
    m = _VARIANT_IN_TITLE.search(title or "")
    if not m:
        return None
    name = re.sub(r"\s*\([^)]*\)\s*$", "", m.group(1)).strip()
    return f"transcript variant {name}" if name else None


def record_titles(accessions: list[str]) -> dict[str, str]:
    """RefSeq record titles by accession — the record you get by searching the ID itself.

    The gene's product report names each isoform ("transcript variant 1"), but that is one
    surface of NCBI and it can be silent. The record's own title carries the same
    designation and is the surface a user checks by hand, so it is the second source when
    the first says nothing. Cache-first, and batched: one esummary call covers every
    accession still missing, rather than one call per row.
    """
    out: dict[str, str] = {}
    missing: list[str] = []
    for a in accessions:
        cached = _read_json(CACHE_DIR / "title" / f"{a}.json")
        if cached and cached.get("title"):
            out[a] = cached["title"]
            continue
        # Free if the sequence was fetched by this build: its defline is the same title.
        seq_cached = _read_json(CACHE_DIR / "sequence" / f"{a}.json") or {}
        if seq_cached.get("title"):
            out[a] = seq_cached["title"]
            _write_json(CACHE_DIR / "title" / f"{a}.json", {"accession": a, "title": out[a]})
            continue
        missing.append(a)
    if not missing:
        return out
    try:
        raw = _http_get_text(
            f"{EUTILS}/esummary.fcgi?db=nuccore&id={','.join(missing)}&retmode=json"
            + (f"&api_key={API_KEY}" if API_KEY else ""))
        result = (json.loads(raw) or {}).get("result") or {}
    except Exception:
        return out                      # a label is not worth failing an analysis over
    for uid in result.get("uids") or []:
        rec = result.get(uid) or {}
        acc = rec.get("accessionversion") or ""
        title = rec.get("title") or ""
        if acc and title:
            out[acc] = title
            _write_json(CACHE_DIR / "title" / f"{acc}.json", {"accession": acc, "title": title})
    return out


def fill_variants(transcripts: list[dict]) -> None:
    """Give every transcript NCBI's variant designation, in place.

    Only for a gene with isoforms to tell apart: NCBI numbers no variant when there is a
    single transcript, and that absence is the mono-isoform case rather than a gap to fill.
    """
    if len(transcripts) < 2:
        return
    missing = [t["accession"] for t in transcripts if not t.get("variant")]
    if not missing:
        return
    titles = record_titles(missing)
    for t in transcripts:
        if not t.get("variant"):
            t["variant"] = variant_in(titles.get(t["accession"]))


# ---------------------------------------------------------------- unplaced records

_GB_EXON = re.compile(r"^ {5}exon\s+<?(\d+)\.\.>?(\d+)\s*$", re.MULTILINE)
_GB_VERSION = re.compile(r"^VERSION\s+(\S+)", re.MULTILINE)
_GB_ACCESSION = re.compile(r"^ACCESSION\s+(.+(?:\n {12}.+)*)", re.MULTILINE)
_GB_GENE = re.compile(r'^ {21}/gene="([^"]+)"', re.MULTILINE)
_GB_GENE_ID = re.compile(r'/db_xref="GeneID:(\d+)"')
_GB_TAXON = re.compile(r'/db_xref="taxon:(\d+)"')


def _parse_genbank_structure(record: str) -> dict | None:
    """One GenBank flat-file record -> {accession, gene, gene_id, tax_id, exon_lengths,
    replaces}.

    `replaces` is the record's secondary accessions — the records it superseded, which is
    how a curated NR names the XR model it was promoted from."""
    v = _GB_VERSION.search(record)
    if not v:
        return None
    acc = v.group(1)
    names = (_GB_ACCESSION.search(record).group(1).split()
             if _GB_ACCESSION.search(record) else [])
    gene, gene_id = _GB_GENE.search(record), _GB_GENE_ID.search(record)
    taxon = _GB_TAXON.search(record)
    return {
        "accession": acc,
        "gene": gene.group(1) if gene else "",
        "gene_id": gene_id.group(1) if gene_id else "",
        "tax_id": taxon.group(1) if taxon else "",
        "exon_lengths": [int(e) - int(b) + 1 for b, e in _GB_EXON.findall(record)],
        "replaces": [n for n in names if n != _base_accession(acc)],
    }


def record_structures(accessions: list[str]) -> dict[str, dict]:
    """Gene, exon lengths and superseded accessions from each record's own GenBank entry.

    Keyed by the accession as ASKED (an unversioned id answers with its current version).
    Cache-first, one efetch for the misses. Best-effort by design: this only ever feeds
    the fallbacks for records NCBI's annotation has not caught up with, and a transcript
    they cannot help is left out exactly as it was before they existed — never a failed
    request.
    """
    out: dict[str, dict] = {}
    missing: list[str] = []
    for a in accessions:
        cached = _read_json(CACHE_DIR / "record" / f"{a}.json")
        if cached:
            out[a] = cached
        else:
            missing.append(a)
    for i in range(0, len(missing), 50):
        chunk = missing[i:i + 50]
        try:
            text = _http_get_text(
                f"{EUTILS}/efetch.fcgi?db=nuccore&id={','.join(chunk)}&rettype=gb&retmode=text"
                + (f"&api_key={API_KEY}" if API_KEY else ""))
        except Exception:
            continue
        asked = {_base_accession(a): a for a in chunk}
        for rec in text.split("\n//"):
            info = _parse_genbank_structure(rec)
            a = asked.get(_base_accession(info["accession"])) if info else None
            if a:
                out[a] = info
                _write_json(CACHE_DIR / "record" / f"{a}.json", info)
    return out


def _place_unplaced_nr(transcripts: list[dict]) -> dict[str, str]:
    """Give reference exons to curated NR transcripts the annotation has not placed yet.

    A freshly curated NR record is listed in the gene's product report with no genomic
    placement until NCBI's next annotation run — HTRA1-AS1's only curated transcript,
    NR_201105.1, is one: the report places the XR models and leaves the NR bare. But such
    a record names the model it was promoted from (its secondary accession, XR_946382),
    and that model IS placed. The placement is borrowed only when it is demonstrably the
    same structure: the superseded model is in this gene's report, placed on the reference, and
    its exons match the record's own exon features one for one, length for length.
    Anything less and the transcript stays out, as before.

    Mutates the matching transcripts' `genomic_locations` in place and returns
    {NR accession: model accession whose placement it took}.
    """
    bare = [t for t in transcripts
            if is_noncoding(t.get("accession_version") or "") and not reference_exons(t)]
    if not bare:
        return {}
    placed = {_base_accession(t.get("accession_version") or ""): t
              for t in transcripts if reference_exons(t)}
    info = record_structures([t["accession_version"].strip() for t in bare])
    took: dict[str, str] = {}
    for t in bare:
        rec = info.get(t["accession_version"].strip())
        if not rec or not rec["exon_lengths"]:
            continue
        for old in rec["replaces"]:
            model = placed.get(old)
            if model and [e - b + 1 for b, e in reference_exons(model)] == rec["exon_lengths"]:
                t["genomic_locations"] = model["genomic_locations"]
                took[t["accession_version"].strip()] = model["accession_version"].strip()
                break
    return took


# ---------------------------------------------------------------- parsing

# "Chromosome 12 Reference GRCh38.p14 Primary Assembly" -> ("12", "GRCh38.p14")
_REFERENCE_NAME = re.compile(r"^Chromosome (\S+) Reference (\S+)")


def reference_location(transcript: dict) -> dict | None:
    """The transcript's placement on its species' REFERENCE assembly, chromosome-level.

    A product report places a transcript on everything NCBI annotates: the reference
    assembly, its patches and alternate loci, and whole alternate assemblies (human
    T2T-CHM13, rat SHRSP_T2T, zebrafish GRCz12tu). One coordinate system has to be chosen,
    and it is read from NCBI's own naming rather than from a table of accessions, which
    would go stale the day a species is re-annotated on a new assembly (rat and zebrafish
    both were, recently):

      * a placement named "Chromosome N Reference <assembly> …" on an NC_ record — the
        chromosome itself. Patches and alternate loci say "Reference" too, but sit on
        NT_/NW_ scaffolds; alternate assemblies say "Alternate".
      * or, where the assembly names nothing (fly, yeast: one assembly, no alternates), the
        placement there is.

    On human this selects exactly what the old "NC_0000" prefix did — checked over every
    curated transcript of ~20k genes (tests/test_species.py pins it on the fixtures).
    """
    locs = [loc for loc in transcript.get("genomic_locations") or []
            if loc.get("genomic_accession_version") and loc.get("exons")]
    for loc in locs:
        if (_REFERENCE_NAME.match(str(loc.get("sequence_name") or ""))
                and str(loc["genomic_accession_version"]).startswith("NC_")):
            return loc
    if locs and not any(loc.get("sequence_name") for loc in locs):
        return locs[0]
    return None


def reference_exons(transcript: dict) -> list[Interval] | None:
    """Reference-assembly (begin,end) exon intervals, 1-based inclusive, or None."""
    loc = reference_location(transcript)
    if not loc:
        return None
    out: list[Interval] = []
    for e in loc.get("exons") or []:
        try:
            b, t = int(e["begin"]), int(e["end"])
        except (KeyError, TypeError, ValueError):
            continue
        out.append((min(b, t), max(b, t)))
    return out or None


def reference_strand(transcript: dict) -> str:
    """Genomic orientation of the transcript on the reference: "+", "-", or "" if unstated.

    Taken from NCBI's explicit `orientation` rather than inferred from exon order, so a
    single-exon transcript — where there is no exon order to read a direction from — still
    reports its strand.
    """
    loc = reference_location(transcript)
    if not loc:
        return ""
    o = str((loc.get("genomic_range") or {}).get("orientation", "")).lower()
    if not o:
        ex = loc.get("exons") or []
        o = str(ex[0].get("orientation", "")).lower() if ex else ""
    return {"plus": "+", "minus": "-"}.get(o, "")


def _assembly_label(name: str) -> str:
    """"GRCh38.p14" -> "GRCh38": a patch release moves no primary-assembly coordinate."""
    return re.sub(r"\.p\d+$", "", name or "")


def _labels_from_gene_report(gene_id: str, genomic_accession: str) -> tuple[str, str]:
    """(chromosome, assembly) for a placement the product report left unnamed — fly and
    yeast. The gene-level record names both: chromosome "2R" on "Release 6 plus ISO1 MT"."""
    gene = get_gene_report(gene_id)
    chromosome = str((gene.get("chromosomes") or [""])[0])
    assembly = ""
    for ann in gene.get("annotations") or []:
        for loc in ann.get("genomic_locations") or []:
            if loc.get("genomic_accession_version") == genomic_accession:
                assembly = str(ann.get("assembly_name") or "")
                chromosome = str(loc.get("sequence_name") or chromosome)
    return chromosome, _assembly_label(assembly)


def _reference_labels(transcript: dict, gene_id: str) -> tuple[str, str]:
    """(chromosome, assembly) of the transcript's reference placement."""
    loc = reference_location(transcript)
    if not loc:
        return "", ""
    m = _REFERENCE_NAME.match(str(loc.get("sequence_name") or ""))
    if m:
        return m.group(1), _assembly_label(m.group(2))
    return _labels_from_gene_report(gene_id, str(loc.get("genomic_accession_version") or ""))


# The curated RefSeq RNA classes the tool analyzes: NM_ (mRNA) and NR_ (non-coding RNA).
# Both are real molecules in the cDNA a primer pair meets — GAPDH's NR_152150 is amplified
# by any pair that fits it, whether or not it encodes anything — so a gene's NR transcripts
# are siblings of its NM ones, and a gene with only NR transcripts (HTRA1-AS1) is still a
# gene to design for. Model transcripts (XM_/XR_) stay out: predictions, not curated records.
REFSEQ_PREFIXES = ("NM_", "NR_")


def is_noncoding(accession: str) -> bool:
    """True for an NR_ (non-coding RNA) accession."""
    return accession.strip().upper().startswith("NR_")


@dataclass
class GeneTranscripts:
    """A gene's curated transcripts as the engine analyzes them."""
    gene_id: str
    symbol: str
    description: str
    species: Species
    chromosome: str = ""
    # Every transcript of a gene is transcribed from the same strand; this is the first's.
    strand: str = ""
    assembly: str = ""               # the reference assembly the exons are on, e.g. "GRCm39"
    # [{accession, variant, is_mane, exons, strand, cds, placed_via}, ...]
    transcripts: list[dict] = field(default_factory=list)
    # Curated accessions NCBI lists for the gene WITHOUT a reference placement. Named so an
    # empty result can say why: rat Gapdh's only record, NM_017008.5, replaced .4 two weeks
    # after the annotation run that would have placed it, and is listed bare.
    unplaced: list[str] = field(default_factory=list)


def refseq_transcripts(product_report: dict, sp: Species = HUMAN) -> GeneTranscripts:
    """The gene's curated transcripts, NM and NR alike, on the species' reference assembly."""
    reports = product_report.get("reports") or []
    if not reports:
        raise NotFound("Empty product report")
    product = reports[0].get("product") or {}
    gene = GeneTranscripts(
        gene_id=str(product.get("gene_id", "")),
        symbol=str(product.get("symbol", "")),
        description=str(product.get("description", "")),
        # The report says whose gene it is; the caller's species is only what was asked for.
        species=species_mod.by_tax_id(product.get("tax_id")) or sp,
    )
    placed_via = _place_unplaced_nr(product.get("transcripts") or [])
    for t in product.get("transcripts") or []:
        acc = (t.get("accession_version") or "").strip()
        if not acc.startswith(REFSEQ_PREFIXES):
            continue
        exons = reference_exons(t)
        if not exons:
            gene.unplaced.append(acc)
            continue
        if not gene.transcripts:
            gene.chromosome, gene.assembly = _reference_labels(t, gene.gene_id)
        tx_strand = reference_strand(t)
        if not gene.strand:
            gene.strand = tx_strand
        gene.transcripts.append({
            "accession": acc,
            # NCBI's own designation for the isoform, e.g. "transcript variant 5". Absent on
            # a gene with a single transcript — there is no variant to number — which is
            # what the clients render as "mono-isoform".
            "variant": (t.get("name") or "").strip() or None,
            # MANE Select is a human designation; no other species has one.
            "is_mane": t.get("select_category") == "MANE_SELECT",
            "exons": exons,
            "strand": tx_strand,
            "cds": _cds_range(t),   # (begin, end) in 1-based transcript coords, or None
            # The model whose placement this record took — see _place_unplaced_nr.
            "placed_via": placed_via.get(acc),
        })
    return gene


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


ACCESSION_RE = re.compile(r"^N[MR]_\d+(\.\d+)?$", re.IGNORECASE)


def is_valid_refseq(accession: str) -> bool:
    """A curated RefSeq RNA accession: NM_ (mRNA) or NR_ (non-coding RNA)."""
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


# One accession index per species: data/index/<slug>.json, as two parallel arrays —
# {"accessions": [sorted, upper-case order], "genes": [...]} — because a quarter of a million
# {"accession":…, "gene":…} objects cost several times the bytes and the memory.
INDEX_DIR = Path(os.environ.get("TMJ_INDEX_DIR",
                                Path(__file__).resolve().parent.parent / "data" / "index"))
_INDEX: dict[str, tuple[list[str], list[str]]] = {}


def _load_index(sp: Species) -> tuple[list[str], list[str]]:
    """Lazy-load one species' NM/NR accession index: (accessions sorted for bisect, genes)."""
    if sp.slug not in _INDEX:
        data = _read_json(INDEX_DIR / f"{sp.slug}.json") or {}
        _INDEX[sp.slug] = (data.get("accessions") or [], data.get("genes") or [])
    return _INDEX[sp.slug]


def suggest_accessions(q: str, limit: int = 8) -> list[dict]:
    """Typeahead: NM/NR transcripts, of any supported species, whose accession starts with `q`.

    An accession belongs to exactly one organism, so the box that takes one needs no species
    selector: every species' index is searched and each row says whose it is. Returns
    [{accession, gene, species}].

    Primary source is the local indexes of real accessions (built from NCBI data).
    Reliable, instant, comprehensive. Falls back to a live NCBI lookup only if no index has a
    match (e.g. an accession newer than the indexes).
    See vault `02 Data/NCBI Data Source and API.md`.

    The bare prefix is enough to open the list: typing "NR_" is how a user asks whether NR
    accessions are searchable at all, and an empty dropdown answers "no". The index serves
    that for free; only the live fallback waits for 5 characters, since a 3-character
    Entrez wildcard matches everything and means nothing.
    """
    base = q.strip().upper().split(".")[0]
    if not base.startswith(REFSEQ_PREFIXES):
        return []

    from bisect import bisect_left
    out: list[dict] = []
    for sp in species_mod.SPECIES:
        accessions, genes = _load_index(sp)
        i = bisect_left(accessions, base)
        stop = min(len(accessions), i + limit)       # no species can contribute more
        while i < stop and accessions[i].startswith(base):
            out.append({"accession": accessions[i], "gene": genes[i], "species": sp.slug})
            i += 1
    if out or len(base) < 5:
        out.sort(key=lambda r: r["accession"])
        return out[:limit]
    return _suggest_live(base, limit)


_GENE_TYPES = " OR ".join(
    f'"genetype {g}"[Properties]'
    for g in ("protein coding", "ncrna", "snorna", "snrna", "scrna", "rrna"))
# Not `suggest_gene`: the prefixes cached there were answered under the protein-coding-only
# filter, and would keep hiding every non-coding gene they predate. One folder per species.
_GENE_SUGGEST_DIR = CACHE_DIR / "suggest_gene_rna"


def suggest_genes(q: str, sp: Species = HUMAN, limit: int = 8) -> list[dict]:
    """Typeahead: a species' RNA-producing genes whose official symbol starts with `q`.

    Live NCBI E-utilities (db=gene), per-prefix cached. Protein-coding genes and the
    non-coding RNA gene types (lncRNA/antisense/miRNA as ncRNA, plus sno/sn/sc/rRNA) — the
    genes that carry NM or NR transcripts to analyze. Pseudogenes stay out of the
    typeahead: most have no transcript at all, and the few transcribed ones still resolve
    when their symbol is searched in full. Scoped by taxonomy id without subtree expansion,
    so it offers exactly the genes a lookup under that id will find. Best-effort — any
    error (incl. E-utilities rate limiting without an API key) yields []. Non-empty results
    only are cached, so a transient failure doesn't poison the prefix. Symbols come back in
    NCBI's official capitalization. Returns [{symbol, description}]."""
    base = q.strip().upper()
    if len(base) < 2:
        return []
    cache_path = _GENE_SUGGEST_DIR / sp.slug / f"{_file_safe(base)}.json"
    cached = _read_json(cache_path)
    if cached is not None:
        return cached[:limit]
    try:
        term = (f"{base}*[Preferred Symbol] AND txid{sp.tax_id}[Organism:noexp] "
                f"AND alive[prop] AND ({_GENE_TYPES})")
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
            _write_json(cache_path, out[:limit])
        return out[:limit]
    except Exception:
        return []


def _suggest_live(base: str, limit: int) -> list[dict]:
    """Fallback: NCBI E-utilities lookup (per-prefix cached), across the supported species.
    Post-filtered to true prefix matches (Entrez's accession wildcard is noisy).
    Best-effort — errors yield []."""
    # `suggest_refseq`, not `suggest`: what was cached there was asked of human alone.
    cache_path = CACHE_DIR / "suggest_refseq" / f"{base}.json"
    cached = _read_json(cache_path)
    if cached is not None:
        return cached[:limit]
    try:
        # No biomol filter: the NM_/NR_ prefix already names the molecule class, and the
        # post-filter below keeps only true prefix matches.
        taxa = " OR ".join(f"txid{sp.tax_id}[ORGN]" for sp in species_mod.SPECIES)
        term = f"{base}*[ACCN] AND srcdb_refseq[PROP] AND ({taxa})"
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
                sp = species_mod.by_tax_id(r.get("taxid"))
                if acc.upper().startswith(base) and sp:   # drop Entrez wildcard noise
                    out.append({"accession": acc, "gene": _symbol_from_title(r.get("title", "")),
                                "species": sp.slug})
        _write_json(cache_path, out[:limit])
        return out[:limit]
    except Exception:
        return []
