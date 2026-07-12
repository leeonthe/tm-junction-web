# Deploying TmJunction

Two services: the Python **engine** on **Render**, the React **web** app on **Vercel**.
Both auto-deploy on push to `main`.

## 1. Engine → Render
1. **render.com → New + → Blueprint** → connect this repo (`leeonthe/tm-junction-web`).
   Render reads [`render.yaml`](render.yaml) and creates the `tmjunction-engine` web service
   (root dir `engine/`, `pip install -r requirements.txt`,
   `uvicorn app.main:app --host 0.0.0.0 --port $PORT`, health check `/health`).
   *Manual alt:* New + → Web Service → root dir `engine`, same build/start commands.
2. In the service → **Environment**, set **`NCBI_API_KEY`** to your NCBI key (optional but
   faster). `TMJ_CORS_ORIGINS` defaults to `*`; set it to your Vercel URL to restrict.
3. Deploy, then note the URL, e.g. `https://tmjunction-engine.onrender.com`. Verify `/health`.
   > ⚠️ The **free** plan sleeps after ~15 min idle → the first request is slow (cold start).
   > Upgrade to a paid instance for always-on.

## 2. Web → Vercel
1. **vercel.com → Add New → Project** → import this repo.
2. **Root Directory = `web`** (Vercel auto-detects Vite; build `npm run build`, output `dist`).
3. Add env var **`VITE_API_URL`** = your Render engine URL from step 1.3 (Production + Preview).
4. Deploy. [`web/vercel.json`](web/vercel.json) handles SPA routing (so `/test` etc. work).

## 3. Wire them together
- Web `VITE_API_URL` → the Render engine URL.
- Engine `TMJ_CORS_ORIGINS` → your Vercel domain (or leave `*`).
- **Redeploy the web app after setting `VITE_API_URL`** — Vite bakes env vars at build time.

## Notes
- `engine/data/nm_index.json` + the GAPDH/MYC cache fixtures are committed, so typeahead and
  the demo genes work immediately with no NCBI calls.
- `primer3-py` installs from a prebuilt wheel — no compiler or Docker needed.
- **Persistent cache (optional):** attach a Render Disk mounted at `engine/data/cache` (paid)
  to keep live-fetched genes cached across restarts; otherwise they re-fetch on a cache miss.
