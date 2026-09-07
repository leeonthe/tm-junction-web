# Exon Junction Primer

Primer design tool for transcript-specific PCR/qPCR.

Most human genes produce several transcript variants, and for many variants there is no
stretch of exonic sequence that belongs to that variant alone — an ordinary primer pair
will amplify the siblings too. This tool tells you, for any RefSeq NM transcript, whether
it can be amplified distinctly from every other NM isoform of its gene, and designs the
primers: a conventional pair when a unique exonic region exists, or an exon–exon junction
(EEJ) primer that spans the one splice junction no sibling carries.

Companion tool to the paper *"Tm-guided exon–exon junction RT-PCR enables specific
detection of RNA variants lacking easily distinguishable exonic regions."*

Live at **https://tm-junction-web.vercel.app**

## What it does

Search by gene symbol, RefSeq accession, or paste your own junction sequence. Every NM
isoform of the gene is fetched from NCBI (sequence + exon structure, GRCh38) and each
transcript is classified by a sequence rule: slide a 20-nt window along the target mRNA —
a window found in no sibling is a target-specific primer site.

- **EEJ-independent** — some unique window sits inside a single exon; an ordinary primer
  pair works.
- **EEJ-dependent** — no exon-internal window is unique, but one exon–exon junction is;
  the primer must straddle that splice.
- **Infeasible** — neither exists; only a junction-combination strategy can isolate it.

For EEJ-dependent transcripts the junction designer applies the paper's Tm rule: the whole
primer must melt inside your chosen range while each arm alone stays ≥ 15 °C below it, so
the oligo primes only across that exact junction — not on pre-mRNA or on isoforms carrying
just one side. Tm is nearest-neighbor (SantaLucia 1998) with salt/Mg²⁺/dNTP corrections you
can edit; arms under 10 nt fall back to the Wallace rule. You can drag-select any window on
the junction and watch it re-evaluate live.

Some details that matter in practice:

- RefSeq sometimes issues several accessions for one molecule (TP53 has 25 NM accessions
  for 13 distinct mRNAs). Accessions with identical sequence *and* identical exon structure
  are treated as one transcript — otherwise every one would be compared against its own
  twin and nothing would look designable. The folded accessions are listed under the row
  that represents them.
- Every amplicon must span at least two exons, so a product can't be confused with one
  amplified off contaminating genomic DNA.
- Besides isoform-specific design, an **all-variant pair** is computed per gene: one primer
  pair that amplifies as many variants as possible at a single product size, for total-
  expression assays.
- NCBI's variant designation ("transcript variant 5") is shown under each accession.

The in-app **Method** page documents the classification rule, the Tm formula, and the
constants.

## Repository layout

    engine/   Python + FastAPI. NCBI client (cache-first), amplifiability classification,
              primer design (primer3 for QC thermodynamics), pan-variant search.
    web/      React + Vite + TypeScript. Search, verdict tables, exon-track graphs,
              the interactive junction designer.
    web/api/  The engine again, as a Vercel Python function — a committed copy of
              engine/ (see web/scripts/vendor-engine.mjs) so one push deploys both.

The web app talks to whichever engine answers `/health` with the capabilities it needs:
the URL in `VITE_API_URL` when that engine is current, else its own same-origin `/api`.

## Running locally

    # engine
    cd engine
    python3 -m venv .venv && . .venv/bin/activate
    pip install -r requirements.txt
    uvicorn app.main:app --reload --port 8000

    # web, in another terminal
    cd web
    npm install && cp .env.example .env
    npm run dev            # http://localhost:5173

Commonly used genes (GAPDH, TP53, ACTB, EGFR, …) are seeded in `engine/data/cache` and
work offline; anything else fetches live from NCBI E-utilities/Datasets. Set
`NCBI_API_KEY` to raise the NCBI rate limit.

Tests: `pytest` in `engine/` (regression suite pinned to validated gene matrices),
`npx vitest run` in `web/`. The web suite includes a guard that fails when `web/api/_engine`
drifts from `engine/` — run `node scripts/vendor-engine.mjs` after engine changes.

## Engine API

    GET  /health                    liveness + implemented capabilities
    GET  /gene/{symbol}             NM transcripts + exon alignment (pick-a-variant step)
    POST /analyze {"accession"}     full classification + primer design for one transcript
    GET  /analyze/stream?accession= same, as SSE with progress events

## License

MIT — see [LICENSE](LICENSE). Depends on
[primer3-py](https://github.com/libnano/primer3-py) (GPL-2.0) as a pip package for primer
QC thermodynamics. Transcript data comes from NCBI RefSeq.
