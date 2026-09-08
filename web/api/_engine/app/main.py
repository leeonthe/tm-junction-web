"""FastAPI app — the TmJunction engine HTTP surface."""

from __future__ import annotations

import json
import os

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse, StreamingResponse
from pydantic import BaseModel

from . import ncbi
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


@app.get("/health")
def health() -> dict:
    """Liveness plus what this build implements — see models.FEATURES.

    A deploy that never happened looks exactly like a healthy engine until you ask it for
    something it does not have, so /health answers both questions at once.
    """
    return {"status": "ok", "features": FEATURES}


@app.get("/suggest")
def suggest(q: str, limit: int = 8) -> dict:
    """Live NCBI typeahead for NM accessions matching prefix `q`."""
    return {"suggestions": ncbi.suggest_accessions(q, limit)}


@app.get("/suggest_genes")
def suggest_genes(q: str, limit: int = 8) -> dict:
    """Live NCBI typeahead for human protein-coding gene symbols matching prefix `q`."""
    return {"suggestions": ncbi.suggest_genes(q, limit)}


# A 429 from NCBI that survived the client's own retries: the request was fine, the
# shared egress IP is over NCBI's per-IP budget this second. 503 + the honest message,
# because the 500-with-stack-trace this used to become reads as "the tool is broken".
_RATE_LIMIT = {"error": "NCBI_RATE_LIMIT", "message": "NCBI is rate-limiting requests from this server right now — nothing is wrong with your query. Try again in a few seconds."}


@app.get("/gene/{symbol}", response_model=GeneLookupResponse)
def gene_lookup(symbol: str):
    """Human gene name -> its NM transcripts + exon alignment (no classification).
    A reference step so users can pick which variant to analyze."""
    try:
        return lookup_gene(symbol)
    except AnalysisError as e:
        status = 404 if e.code == "NOT_FOUND" else 400
        return JSONResponse(status_code=status, content={"error": e.code, "message": e.message})
    except ncbi.RateLimited:
        return JSONResponse(status_code=503, content=_RATE_LIMIT,
                            headers={"Retry-After": "5"})


@app.post("/analyze", response_model=AnalyzeResponse)
def analyze_endpoint(req: AnalyzeRequest):
    try:
        return analyze(req.accession, k=req.k)
    except AnalysisError as e:
        return JSONResponse(status_code=400, content={"error": e.code, "message": e.message})
    except ncbi.RateLimited:
        return JSONResponse(status_code=503, content=_RATE_LIMIT,
                            headers={"Retry-After": "5"})


@app.get("/analyze/stream")
def analyze_stream(accession: str, k: int = 20):
    """Server-Sent Events: real progress (`progress`) then the full result (`result`).
    Analysis errors are sent as a `failed` event (named so it doesn't clash with the
    browser EventSource's built-in `error`)."""
    def gen():
        try:
            for ev in analyze_events(accession, k):
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
def analyze_get(accession: str, k: int = 20):
    try:
        return analyze(accession, k=k)
    except AnalysisError as e:
        return JSONResponse(status_code=400, content={"error": e.code, "message": e.message})
    except ncbi.RateLimited:
        return JSONResponse(status_code=503, content=_RATE_LIMIT,
                            headers={"Retry-After": "5"})
