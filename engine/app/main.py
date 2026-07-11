"""FastAPI app — the TmJunction engine HTTP surface."""

from __future__ import annotations

import os

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse
from pydantic import BaseModel

from . import ncbi
from .analyze import AnalysisError, analyze
from .models import AnalyzeResponse

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
    return {"status": "ok"}


@app.get("/suggest")
def suggest(q: str, limit: int = 8) -> dict:
    """Live NCBI typeahead for NM accessions matching prefix `q`."""
    return {"suggestions": ncbi.suggest_accessions(q, limit)}


@app.post("/analyze", response_model=AnalyzeResponse)
def analyze_endpoint(req: AnalyzeRequest):
    try:
        return analyze(req.accession, k=req.k)
    except AnalysisError as e:
        return JSONResponse(status_code=400, content={"error": e.code, "message": e.message})


@app.get("/analyze/{accession}", response_model=AnalyzeResponse)
def analyze_get(accession: str, k: int = 20):
    try:
        return analyze(accession, k=k)
    except AnalysisError as e:
        return JSONResponse(status_code=400, content={"error": e.code, "message": e.message})
