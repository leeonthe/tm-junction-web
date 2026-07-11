# TmJunction Web (React + Vite + TypeScript)

The frontend: a single-page tool that calls the Python engine and renders the verdict,
Tm-guided primers, an exon-track graph, and the per-isoform table — in the Figma design
system (Inter, `#237AF2` pill buttons, bordered cards). Deployable independently as a
static site.

## Layout
```
src/
  App.tsx                 state + result composition (tabs)
  lib/  api.ts            fetch → engine /analyze
        types.ts          mirrors engine AnalyzeResponse
        tier.ts           tier → color / label / class
  components/
        Nav.tsx  Hero.tsx        nav + search
        VerdictBanner.tsx        3-tier verdict hero
        PrimerCard.tsx           designed primers (seq, Tm, ΔTm, copy)
        ExonTrackGraph.tsx       data-driven SVG, colored by tier, ▾ junction carets
        VerdictTable.tsx         per-isoform table
  styles/app.css          design tokens (light + dark) + component styles
```

## Run locally
```bash
npm install
cp .env.example .env        # VITE_API_URL=http://localhost:8000
npm run dev                 # → http://localhost:5173  (start the engine first)
```

## Build / check
```bash
npm run build       # tsc -b && vite build  → dist/
npm run typecheck
```

## Config
`VITE_API_URL` — URL of the Python engine. Local default `http://localhost:8000`.

## Deploy
Static host (Vercel / Netlify / Cloudflare Pages / S3). `npm run build` → serve `dist/`.
Set `VITE_API_URL` to your deployed engine URL at build time, and add that web origin to
the engine's `TMJ_CORS_ORIGINS`.

## Notes
- Fonts: install **Inter** via the host (or `@fontsource/inter`) for pixel-exact type;
  falls back to system-ui otherwise.
- Theme: light (faithful to Figma) + dark; toggle in the nav (◐).
