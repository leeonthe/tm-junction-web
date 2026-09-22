"""FastAPI app — the TmJunction engine HTTP surface."""

from __future__ import annotations

import json
import os

from fastapi import FastAPI, Query
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse, StreamingResponse
from pydantic import BaseModel

from . import ncbi, primers
from . import species as species_mod
from .analyze import AnalysisError, analyze, analyze_events, lookup_gene
from .models import FEATURES, AnalyzeResponse, GeneLookupResponse

app = FastAPI(title="TmJunction Engine", version="0.1.0")

# CORS — allow the web frontend origin(s). Configure via TMJ_CORS_ORIGINS:
#   "*"                     -> allow any origin (fine here: no credentials/cookies)
#   "https://app.example"   -> comma-separated explicit list (production)
# Local dev always accepts http(s)://localhost|127.0.0.1 on any port via a regex, so
# the localhost-vs-127.0.0.1 and Vite-port mismatches don't trip the preflight.
_env = os.environ.get("TMJ_CORS_ORIGINS", "").strip()
_cors: dict = {"allow_origin_regex": r"https?://(localhost|127\.0\.0\.1)(:\d+)?"}
if _env == "*":
    _cors = {"allow_origins": ["*"]}
elif _env:
    _cors["allow_origins"] = [o.strip() for o in _env.split(",") if o.strip()]

app.add_middleware(
    CORSMiddleware,
    allow_methods=["*"],
    allow_headers=["*"],
    **_cors,
)


class AnalyzeRequest(BaseModel):
    accession: str
    k: int = 20
    # Transcripts to leave out of the comparison (see analyze._apply_exclusion).
    exclude: list[str] = []


def _exclude_list(exclude: str | None) -> list[str]:
    """`exclude=NM_1,NM_2` on a GET — comma-separated, blanks dropped."""
    return [a for a in (exclude or "").split(",") if a.strip()]


@app.get("/health")
def health() -> dict:
    """Liveness plus what this build implements — see models.FEATURES.

    A deploy that never happened looks exactly like a healthy engine until you ask it for
    something it does not have, so /health answers both questions at once.
    """
    return {"status": "ok", "features": FEATURES}


_QC_MAX_OLIGOS = 50
_QC_MAX_LEN = 60


@app.get("/qc")
def qc(seq: list[str] = Query(default=[])) -> dict:
    """Structure Tm for oligos the browser designed itself.

    The junction designer's second-primer options are computed client-side, and the browser
    can judge Tm, GC, the 3' clamp and homopolymer runs on its own — but hairpin and
    self-dimer stability are primer3's, and live here. This returns just those two numbers
    per oligo, plus the thresholds the engine's own QC gate (primers._evaluate) applies, so
    the page can print the same QC-passed / QC-relaxed verdict the whole-transcript card
    prints, judged by the same rule. Stateless and cheap; unknown characters give 0.0.
    """
    out = []
    for s in seq[:_QC_MAX_OLIGOS]:
        s = (s or "").strip().upper()[:_QC_MAX_LEN]
        if not s or any(c not in "ACGT" for c in s):
            out.append({"seq": s, "hairpin_tm": 0.0, "homodimer_tm": 0.0})
            continue
        out.append({"seq": s, "hairpin_tm": primers.hairpin_tm(s),
                    "homodimer_tm": primers.homodimer_tm(s)})
    return {
        "results": out,
        "thresholds": {
            "tm_min": primers.TM_MIN, "tm_max": primers.TM_MAX,
            "gc_min": primers.GC_MIN, "gc_max": primers.GC_MAX,
            "struct_tm_max": primers.STRUCT_TM_MAX, "poly_max": 5,
        },
    }


_SEQ_MAX_ACCESSIONS = 500


@app.get("/sequences")
def sequences(acc: list[str] = Query(default=[])):
    """Transcript sequences for a set of NM/NR accessions — the whole-transcript designer's siblings.

    The analysis response carries ONE sequence, the target's, because the browser only
    designs against that one. The whole-transcript designer is the exception: it places a
    pair on the target and then has to show which of the gene's other transcripts that
    pair amplifies, at what size — a claim the engine makes from the sequences (see
    panvariant._coverage), and the browser must make the same way or not at all. Sending
    every sibling's mRNA with every analysis would cost a megabyte on a gene like BRCA1
    for a tab most searches never open, so they are fetched here, on demand.

    Cache-first through ncbi.get_sequences, so after an analysis this is a disk read; a
    cold serverless instance pays one efetch for the lot. Accessions that are not NM/NR
    RefSeq ids are a 400, an id NCBI does not know is a 404 — the browser sent
    accessions the engine itself reported, so either is a bug, not a user error.
    """
    accs: list[str] = []
    for a in acc[:_SEQ_MAX_ACCESSIONS]:
        a = (a or "").strip().upper()
        if not a:
            continue
        if not ncbi.is_valid_refseq(a):
            return JSONResponse(status_code=400, content={
                "error": "NOT_NM", "message": f"{a} is not an NM or NR RefSeq accession."})
        if a not in accs:
            accs.append(a)
    if not accs:
        return {"sequences": {}}
    try:
        seqs = ncbi.get_sequences(accs)
    except ncbi.NotFound as e:
        return JSONResponse(status_code=404, content={"error": "NOT_FOUND", "message": str(e)})
    except ncbi.RateLimited:
        return JSONResponse(status_code=503, content=_RATE_LIMIT,
                            headers={"Retry-After": "5"})
    return {"sequences": {a: seqs[a] for a in accs if seqs.get(a)}}


@app.get("/species")
def species_list() -> dict:
    """The species the engine analyzes, in display order — slug, names, taxonomy id."""
    return {"species": [
        {"slug": s.slug, "tax_id": s.tax_id, "scientific": s.scientific,
         "common": s.common, "example": s.example} for s in species_mod.SPECIES]}


def _bad_species(e: Exception) -> JSONResponse:
    return JSONResponse(status_code=400, content={"error": "BAD_SPECIES", "message": str(e)})


@app.get("/suggest")
def suggest(q: str, limit: int = 8) -> dict:
    """Typeahead for NM/NR accessions matching prefix `q`, across every supported species
    (an accession names its own organism; each row says which)."""
    return {"suggestions": ncbi.suggest_accessions(q, limit)}


@app.get("/suggest_genes")
def suggest_genes(q: str, limit: int = 8, species: str | None = None):
    """Live NCBI typeahead for a species' gene symbols (protein-coding + ncRNA) matching
    prefix `q`. Human when `species` is omitted."""
    try:
        sp = species_mod.get(species)
    except species_mod.UnknownSpecies as e:
        return _bad_species(e)
    return {"suggestions": ncbi.suggest_genes(q, sp, limit), "species": sp.slug}


# A 429 from NCBI that survived the client's own retries: the request was fine, the
# shared egress IP is over NCBI's per-IP budget this second. 503 + the honest message,
# because the 500-with-stack-trace this used to become reads as "the tool is broken".
_RATE_LIMIT = {"error": "NCBI_RATE_LIMIT", "message": "NCBI is rate-limiting requests from this server right now — nothing is wrong with your query. Try again in a few seconds."}


@app.get("/gene/{symbol:path}", response_model=GeneLookupResponse)
def gene_lookup(symbol: str, species: str | None = None):
    """Gene name (+ `species`, human when omitted) -> its NM/NR transcripts + exon
    alignment (no classification). A reference step so users can pick which variant to
    analyze. `:path` because a symbol is NCBI's to spell, and fly's carry punctuation."""
    try:
        return lookup_gene(symbol, species)
    except AnalysisError as e:
        status = 404 if e.code in ("NOT_FOUND", "NOT_PLACED") else 400
        return JSONResponse(status_code=status, content={"error": e.code, "message": e.message})
    except ncbi.RateLimited:
        return JSONResponse(status_code=503, content=_RATE_LIMIT,
                            headers={"Retry-After": "5"})


@app.post("/analyze", response_model=AnalyzeResponse)
def analyze_endpoint(req: AnalyzeRequest):
    try:
        return analyze(req.accession, k=req.k, exclude=req.exclude)
    except AnalysisError as e:
        return JSONResponse(status_code=400, content={"error": e.code, "message": e.message})
    except ncbi.RateLimited:
        return JSONResponse(status_code=503, content=_RATE_LIMIT,
                            headers={"Retry-After": "5"})


@app.get("/analyze/stream")
def analyze_stream(accession: str, k: int = 20, exclude: str | None = None):
    """Server-Sent Events: real progress (`progress`) then the full result (`result`).
    Analysis errors are sent as a `failed` event (named so it doesn't clash with the
    browser EventSource's built-in `error`)."""
    def gen():
        try:
            for ev in analyze_events(accession, k, _exclude_list(exclude)):
                if ev.get("type") == "result":
                    yield f"event: result\ndata: {json.dumps(ev['result'].model_dump())}\n\n"
                else:
                    yield ("event: progress\n"
                           f"data: {json.dumps({'pct': ev['pct'], 'detail': ev['detail']})}\n\n")
        except AnalysisError as e:
            yield f"event: failed\ndata: {json.dumps({'error': e.code, 'message': e.message})}\n\n"
        except ncbi.RateLimited:
            yield f"event: failed\ndata: {json.dumps(_RATE_LIMIT)}\n\n"
        except Exception as e:  # pragma: no cover - upstream/network failure
            yield f"event: failed\ndata: {json.dumps({'error': 'UPSTREAM_ERROR', 'message': str(e)})}\n\n"

    return StreamingResponse(gen(), media_type="text/event-stream",
                             headers={"Cache-Control": "no-cache", "X-Accel-Buffering": "no"})


@app.get("/analyze/{accession}", response_model=AnalyzeResponse)
def analyze_get(accession: str, k: int = 20, exclude: str | None = None):
    try:
        return analyze(accession, k=k, exclude=_exclude_list(exclude))
    except AnalysisError as e:
        return JSONResponse(status_code=400, content={"error": e.code, "message": e.message})
    except ncbi.RateLimited:
        return JSONResponse(status_code=503, content=_RATE_LIMIT,
                            headers={"Retry-After": "5"})
