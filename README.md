# TmJunction

Web tool for the paper *"Tm-guided exon–exon junction RT-PCR enables specific detection
of RNA variants lacking easily distinguishable exonic regions."* Paste a RefSeq **NM**
transcript accession → learn whether it can be amplified distinctly and get **Tm-guided,
isoform-specific RT-PCR primers** (conventional or exon–exon junction).

## Two independently deployable directories

| dir | stack | role |
|---|---|---|
| [`engine/`](engine/) | Python · FastAPI | the science: sequence amplifiability + Tm-guided primer design + NCBI data |
| [`web/`](web/) | React · Vite · TypeScript | the UI: search → verdict → primers → exon-track graph |

The web app calls the engine over HTTP (`POST /analyze`). Deploy the engine to a
container host and the web app to any static host; wire them with `VITE_API_URL`
(web) and `TMJ_CORS_ORIGINS` (engine).

## Quick start (local)
```bash
# 1. engine
cd engine && python3 -m venv .venv && . .venv/bin/activate
pip install -r requirements.txt && uvicorn app.main:app --reload --port 8000

# 2. web (new terminal)
cd web && npm install && cp .env.example .env && npm run dev
# open http://localhost:5173  — the GAPDH demo loads automatically
```
GAPDH and MYC run fully offline from seeded fixtures; other genes fetch live from NCBI.

## How it works (the science)
The verdict is **sequence-based**: slide a primer-length window along the target mRNA;
a window absent from every sibling isoform is a target-specific primer site — either
**exon-internal** (conventional primer) or **junction-spanning** (EEJ primer). Three tiers:
- 🔵 **Conventional** — has a unique region → ordinary primer.
- 🔴 **Needs EEJ** — no unique region, but a unique junction → junction-spanning primer.
- ⚫ **Hard case** — neither → needs a junction-combination strategy.

Full rationale, algorithms, and validation live in the Obsidian vault
[`tm-junction-web-obsidian/`](tm-junction-web-obsidian/) (start at `Home.md`).

## Design
`design-mocks/result-screen.html` is the standalone design reference. Tokens from the
project's Figma system (Inter, `#237AF2`, bordered cards) — see the vault's UI-UX plan.

## Status
Engine: feature-complete for v1 (8 passing regression tests locked to the validated
GAPDH/MYC matrix). Web: builds clean, renders the full result flow. Next: `primer3`
QC pass, genome-wide specificity check, analytics DB, deploy.
