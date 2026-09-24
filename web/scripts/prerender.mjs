// Write crawlable HTML for "/", "/guide" and "/method" into dist/ — run after `vite build`.
//
// Vite's build emits one index.html whose <body> is an empty root. This renders each page
// with React on the server (src/prerender.tsx, built as an SSR bundle) and writes it into
// that shell with the page's own <title>, description, canonical URL, Open Graph tags and
// the site's JSON-LD — so a search engine indexes the text of the page instead of nothing.
// The app still mounts on load and takes over the root.

import { execSync } from "node:child_process";
import { mkdirSync, readFileSync, rmSync, writeFileSync } from "node:fs";
import { dirname, join } from "node:path";
import { fileURLToPath, pathToFileURL } from "node:url";

const web = join(dirname(fileURLToPath(import.meta.url)), "..");
const dist = join(web, "dist");
const ssrDir = join(web, "dist-ssr");

execSync("npx vite build --ssr src/prerender.tsx --outDir dist-ssr --logLevel warn", { cwd: web, stdio: "inherit" });
const mod = await import(pathToFileURL(join(ssrDir, "prerender.js")).href);
const { render, PAGES, structuredData } = mod;
const SITE_URL = "https://tm-junction-web.vercel.app";

const shell = readFileSync(join(dist, "index.html"), "utf8");
const esc = (s) => s.replace(/&/g, "&amp;").replace(/"/g, "&quot;").replace(/</g, "&lt;");

for (const page of PAGES) {
  const body = render(page.path);
  const url = SITE_URL + (page.path === "/" ? "/" : page.path);
  const head = [
    `<title>${esc(page.title)}</title>`,
    `<meta name="description" content="${esc(page.description)}" />`,
    `<link rel="canonical" href="${url}" />`,
    `<meta property="og:type" content="website" />`,
    `<meta property="og:site_name" content="Exon Junction Primer" />`,
    `<meta property="og:title" content="${esc(page.title)}" />`,
    `<meta property="og:description" content="${esc(page.description)}" />`,
    `<meta property="og:url" content="${url}" />`,
    `<meta name="twitter:card" content="summary" />`,
    `<script type="application/ld+json">${structuredData()}</script>`,
  ].join("\n    ");
  let html = shell.replace(/<title>[^<]*<\/title>/, head);
  html = html.replace('<div id="root"></div>', `<div id="root">${body}</div>`);
  const file = page.path === "/" ? "index.html" : page.path.slice(1) + ".html";
  writeFileSync(join(dist, file), html);
  console.log(`prerendered ${page.path} -> dist/${file} (${(body.length / 1024).toFixed(0)} kB of content)`);
}
rmSync(ssrDir, { recursive: true, force: true });
mkdirSync(dist, { recursive: true });
