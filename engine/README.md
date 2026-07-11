# TmJunction Engine (Python / FastAPI)

The scientific core: resolve a RefSeq **NM** accession → gene → sibling NM isoforms,
run the **sequence (ATGC) amplifiability** test, and design **Tm-guided** primers
(conventional or exon–exon junction). Deployable independently of the web frontend.

## Layout
```
app/
  overlap.py   coordinate non-unique rule + junction finder (structural / graph layer)
  amplify.py   ⭐ sequence sliding-window uniqueness → 3-tier verdict (authoritative)
               (also exposes unique-window positions so primer specificity is provable)
  primers.py   ⭐ primer3 QC hybrid: Tm + hairpin/dimer gates, multi-length sweep,
               selection ranked by QC + ΔTm specificity (junction primers ΔTm-weighted)
  ncbi.py      cache-first NCBI Datasets v2 / E-utilities client
  analyze.py   orchestration → AnalyzeResponse
  models.py    pydantic response schemas
  main.py      FastAPI app  (POST /analyze, GET /analyze/{acc}, GET /health)
data/cache/    seeded GAPDH + MYC fixtures → runs fully offline for the demo genes
tests/         regression tests locked to the validated GAPDH/MYC matrix
```

## Run locally
```bash
python3 -m venv .venv && . .venv/bin/activate
pip install -r requirements.txt
uvicorn app.main:app --reload --port 8000
# → http://localhost:8000/health
#   curl -X POST localhost:8000/analyze -H 'content-type: application/json' -d '{"accession":"NM_001256799.3"}'
```

## Test (offline, from fixtures)
```bash
.venv/bin/python -m pytest -q      # 8 tests, ~0.4s
```

## Config (env)
| var | meaning |
|---|---|
| `NCBI_API_KEY` | optional; raises NCBI rate limit 5→10 rps (never logged) |
| `TMJ_CORS_ORIGINS` | comma-separated allowed web origins (default `http://localhost:5173`) |
| `TMJ_CACHE_DIR` | override cache dir (default `./data/cache`) |

For GAPDH/MYC everything is cached, so no key/network needed. Any other gene triggers a
live NCBI fetch (product report + mRNA sequences) that is then cached.

## Deploy
Container host (Render / Railway / Fly.io) or any host that runs `uvicorn`. Set
`TMJ_CORS_ORIGINS` to your deployed web URL. Persist `data/cache/` (a volume) so live
lookups are cached across restarts.

## Contract
`POST /analyze {accession, k?}` → `AnalyzeResponse` (see `app/models.py`). Errors:
`{error, message}` with `NOT_NM` / `NOT_FOUND`. Mirrors the web `src/lib/types.ts`.
