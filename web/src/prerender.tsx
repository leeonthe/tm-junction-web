// Server-side render of the three pages, for search engines.
//
// The app is a single page drawn by JavaScript, so a crawler fetching "/guide" got
// `<div id="root"></div>` — not one of the words the page is about. scripts/prerender.mjs
// runs this after `vite build` and writes real HTML for "/", "/guide" and "/method" into
// dist/, with each page's title, description and canonical URL in its <head>. The browser
// still loads the app, which takes over the root on mount; the static HTML is what is
// there before it does, and what a crawler indexes.
//
// The landing page is rendered here as static copy rather than through <App>: App reads
// the address bar and probes the engine on mount, and what a crawler needs from the front
// page is a plain statement of what the tool does, with the way in.

import { renderToString } from "react-dom/server";
import Guide from "./components/Guide";
import Method from "./components/Method";
import { PAGES, PAPER, SITE_NAME } from "./seo";

const noop = () => {};

function Landing() {
  return (
    <>
      <header className="hero">
        <div className="wrap hero-in">
          <span className="eyebrow"><b />Tm-guided exon–exon junction RT-PCR</span>
          <h1>{SITE_NAME}</h1>
          <p className="lede">
            Primer design tool for transcript-specific PCR/qPCR.<br />
            Search a gene or enter a RefSeq accession<br />
            to find the unique exon or junction and design the primers.
          </p>
          <nav className="examples" aria-label="Pages">
            <a className="chip-ex" href="/?g=GAPDH">Try GAPDH</a>
            <a className="chip-ex" href="/guide">How to use</a>
            <a className="chip-ex" href="/method">Method</a>
          </nav>
        </div>
      </header>
      <main className="wrap">
        <section className="card mth-card">
          <h2>Exon junction primer design for isoform-specific RT-PCR</h2>
          <p>
            Most genes produce several transcript variants, and for many variants there is no
            stretch of exonic sequence that belongs to that variant alone — an ordinary primer
            pair amplifies the siblings too. Exon Junction Primer tells you, for any curated
            RefSeq transcript (NM mRNA or NR non-coding RNA), whether it can be amplified
            distinctly from every other isoform of its gene, and designs the primers: a
            conventional pair where a unique exonic region exists, or a <b>Tm-guided exon–exon
            junction (EEJ) primer</b> that spans the one splice junction no sibling carries.
          </p>
          <ul>
            <li><b>EEJ-independent</b> — a unique window sits inside one exon; a conventional pair works.</li>
            <li><b>EEJ-dependent</b> — no exon is unique but one exon–exon junction is; the primer must
              straddle that splice, and the Tm junction rule sets how much of it lies on each side.</li>
            <li><b>Infeasible</b> — neither exists; only a junction-combination strategy can isolate it.</li>
          </ul>
          <p>
            It also designs <b>one primer pair for the whole gene</b> — every variant at one product
            size — for total-expression assays, and lets you exclude transcripts from the comparison.
            Six species are covered from NCBI RefSeq: human (GRCh38), mouse (GRCm39), rat (GRCr8),
            <i> Drosophila melanogaster</i>, <i>Saccharomyces cerevisiae</i> and zebrafish. Every
            verdict is sequence-based: a primer site counts as specific only if it is absent from
            every other isoform's mRNA.
          </p>
          <p>
            Companion tool to the paper <i>“{PAPER.title}.”</i> Read <a href="/guide">how to use it</a>{" "}
            or the <a href="/method">method</a> — every formula and constant behind the numbers.
          </p>
        </section>
      </main>
    </>
  );
}

/** The static HTML for one page's #root, or "" for a path this does not pre-render. */
export function render(path: string): string {
  switch (path) {
    case "/": return renderToString(<Landing />);
    case "/guide": return renderToString(<Guide onBack={noop} onMethod={noop} backLabel="Back to search" />);
    case "/method": return renderToString(<Method onBack={noop} backLabel="Back to search" />);
    default: return "";
  }
}

export { PAGES };
export { structuredData } from "./seo";
