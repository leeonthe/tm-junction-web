import { useMemo, useState } from "react";
import { ArrowRight } from "./icons";
import { JunctionWorkbench, useJunctionSettings, type JunctionGeom } from "./JunctionWorkbench";
import { BracketGlyph, TargetGlyph, buildCells } from "./GraphCard";

/**
 * The drag demo's junction is real: GAPDH NM_001289745.3, exon 1 | exon 2 — the junction
 * the designer itself recommends for that transcript — 48 nt each side, taken from the
 * seeded RefSeq sequence. A made-up sequence would behave differently from anything the
 * user later meets; this one behaves exactly like the page they are being taught.
 */
const DEMO_ARM5 = "AAGGCGGCAGGGGCGGGCGCAGGCCGGATGTGTTCGCGCCGCTGCGGG";
const DEMO_ARM3 = "CCGAGCCACATCGCTCAGACACCATGGGGAAGGTGAAGGTCGGAGTCA";

/** The junction designer, live, exactly as it runs after a real analysis. */
function DragDemo() {
  const s = useJunctionSettings();
  const geom = useMemo((): JunctionGeom => {
    const seq = DEMO_ARM5 + DEMO_ARM3;
    const jx = DEMO_ARM5.length;
    return {
      seq, jx, leftBound: 0, rightBound: seq.length, winStart: 0, winEnd: seq.length,
      classOf: (i) => (i < jx ? "ex-a" : "ex-b"),
      leftLabel: "5′ arm", rightLabel: "3′ arm",
    };
  }, []);
  return <JunctionWorkbench geom={geom} s={s} reseedKey="guide-demo" />;
}

/** The graph key's colour picker, live: pick a colour, the mini track follows. */
const DEMO_TOKENS = [
  { token: "--eej-single", label: "Single EEJ primer", glyph: <BracketGlyph /> },
  { token: "--eej-combo", label: "Double EEJ primer pair", glyph: <BracketGlyph combo /> },
  { token: "--amp-pair", label: "Primer target site", glyph: <TargetGlyph /> },
];

function ColorDemo() {
  const [overrides, setOverrides] = useState<Record<string, string>>({});
  const [open, setOpen] = useState<string | null>(null);
  const [defaults, setDefaults] = useState<Record<string, string>>({});

  const openFor = (token: string) => {
    setDefaults(() => {
      const cs = getComputedStyle(document.documentElement);
      const out: Record<string, string> = {};
      for (const { token: k } of DEMO_TOKENS) out[k] = cs.getPropertyValue(k).trim().toUpperCase();
      return out;
    });
    setOpen((o) => (o === token ? null : token));
  };
  const effective = { ...defaults, ...Object.fromEntries(
    Object.entries(overrides).map(([k, v]) => [k, v.toUpperCase()])) };

  return (
    <div className="guide-colordemo" style={overrides as React.CSSProperties}>
      {/* A miniature of the exon graph, drawn with the same tokens the picker changes —
          the point of the demo is watching the mark follow the choice. */}
      <svg width="100%" height="56" viewBox="0 0 560 56" aria-label="Colour demo track">
        <line x1={10} y1={34} x2={550} y2={34} stroke="var(--border-2)" strokeWidth={1.5} />
        {[[10, 60], [110, 40], [190, 46], [300, 40], [380, 56], [480, 70]].map(([x, w], i) => (
          <rect key={i} x={x} y={26} width={w} height={16} rx={3}
            fill={i === 2 ? "var(--amp-pair)" : "var(--eej)"} />
        ))}
        <path d="M72 22 L72 12 L108 12 L108 22" fill="none" stroke="var(--eej-single)"
          strokeWidth={2.5} strokeLinecap="round" strokeLinejoin="round" />
        <path d="M342 22 L342 12 L378 12 L378 22" fill="none" stroke="var(--eej-combo)"
          strokeWidth={2.5} strokeLinecap="round" strokeLinejoin="round" />
      </svg>
      <p className="g-note" style={{ margin: "6px 0 0" }}>
        {DEMO_TOKENS.map(({ token, label, glyph }, i) => (
          <span key={token}>
            {i > 0 && <span className="g-sep">·</span>}
            <span className="g-key">
              <button type="button" className="g-key-btn" aria-expanded={open === token}
                aria-label={`Change ${label} color`} onClick={() => openFor(token)}>
                {glyph} {label}
              </button>
              {open === token && (
                <div className="sw-pop up" role="dialog" aria-label={`${label} color`}>
                  <div className="sw-grid">
                    {buildCells(defaults, effective).map((c) => (
                      <button type="button" key={c} aria-label={c}
                        className={`sw-cell${effective[token] === c ? " on" : ""}`}
                        style={{ background: c }}
                        onClick={() => { setOverrides((o) => ({ ...o, [token]: c })); setOpen(null); }} />
                    ))}
                  </div>
                </div>
              )}
            </span>
          </span>
        ))}
      </p>
    </div>
  );
}

/**
 * The how-to page. The Method page answers "why are these numbers right"; this one answers
 * "what do I click". Same layout shell (.mth) so the two read as siblings, but everything
 * here names things as they appear on screen, in the order a first analysis meets them.
 */
export default function Guide({ onBack, onMethod, backLabel }: {
  onBack: () => void; onMethod: () => void; backLabel: string;
}) {
  return (
    <div className="mth">
      <button type="button" className="mth-back" onClick={onBack}>
        <span className="rev"><ArrowRight /></span>{backLabel}
      </button>

      <header className="mth-head">
        <p className="card-label" style={{ margin: 0 }}>How to use</p>
        <h1>From a transcript to its primers</h1>
        <p className="mth-lede">
          The tool answers one question: can this transcript be amplified without also
          amplifying its sibling isoforms — and if so, with which primers. A full run takes
          under a minute.
        </p>
      </header>

      <section className="card mth-card">
        <p className="card-label">1 · Enter your target</p>
        <p>
          Search by <b>gene symbol</b> (e.g. GAPDH) to see every NM transcript of the gene
          drawn against the genome, then click the variant you want. If you already know the
          transcript, switch to <b>NCBI ID</b> and paste the RefSeq accession
          (e.g. NM_002046.7) — it goes straight to analysis. <b>Custom sequence</b> skips
          NCBI entirely: paste the two sides of a junction and design against your own
          sequence.
        </p>
        <p className="mth-note">
          Only curated NM transcripts are used, on GRCh38. If a gene lists fewer transcripts
          than you expect, look under the accessions: RefSeq IDs that share an identical
          sequence and exon structure are one transcript here, and the duplicates are named
          under the row that represents them.
        </p>
      </section>

      <section className="card mth-card">
        <p className="card-label">2 · Read the verdict</p>
        <p>Every isoform of the gene is compared against the others and lands in one of three tiers:</p>
        <ul className="mth-list">
          <li><b className="t-conv-ink">EEJ-independent</b> — part of an exon is unique to
            this transcript. An ordinary primer pair placed there is specific.</li>
          <li><b className="t-eej-ink">EEJ-dependent</b> — no exonic stretch is unique, but
            one exon–exon junction is. A primer spanning that junction is specific.</li>
          <li><b className="t-hard-ink">Infeasible</b> — no single unique region or junction
            exists. One primer pair cannot isolate it; it needs two junction primers.</li>
        </ul>
        <p className="mth-note">
          The yellow marks on the exon graph are primer target sites — the regions that
          actually distinguish your transcript. Hover any exon or bracket for its details.
        </p>
        <p>
          Every colour on the graph can be changed — click a swatch in the legend, or an
          entry in the key under the graph. Try it here; the marks follow your choice:
        </p>
        <ColorDemo />
      </section>

      <section className="card mth-card">
        <p className="card-label">3 · Take the primers</p>
        <p>
          For an <b>EEJ-independent</b> target, the <i>Primer pair options</i> card lists
          ready Tm-matched pairs. Adjust the amplicon window or Tm range if your assay needs
          it, pick a pair, and use <i>Copy pair</i>. <i>View on cDNA</i> shows exactly where
          each primer binds.
        </p>
        <p>
          For an <b>EEJ-dependent</b> target, the junction designer opens on a suggested
          primer (highlighted on the sequence). Drag across the sequence to move or resize
          it — every Tm below recalculates as you go. A primer is valid when its
          whole-primer Tm sits in your range while each arm alone stays at least
          15&nbsp;°C below it; that gap is what stops the primer firing on isoforms that
          carry only one side of the junction. The <i>Second primer</i> panel then finds a
          Tm-matched partner for a normal amplicon, and <i>Copy pair</i> gives both oligos.
        </p>
        <p>
          This one is live — GAPDH's exon&nbsp;1–exon&nbsp;2 junction, the real designer.
          Drag across the letters and watch the verdict and the three Tm cards react;
          drag a window that sits all on one side and see <em>why</em> it fails:
        </p>
        <DragDemo />
        <p className="mth-note">
          Tm values follow the reaction conditions in the settings panel (Na⁺/K⁺, Mg²⁺,
          dNTPs, primer concentration). Set them to match your master mix before trusting
          the numbers.
        </p>
      </section>

      <section className="card mth-card">
        <p className="card-label">4 · Gene-wide tools</p>
        <p>
          The <i>Gene classification</i> tab shows the whole gene: every isoform's exon
          structure, its tier, and its per-isoform verdict — click any row to design for
          that transcript instead. The <i>All-variant pair</i> card is the opposite job:
          one primer pair that amplifies as many of the gene's variants as possible at a
          single product size, for measuring total expression rather than one isoform.
        </p>
      </section>

      <section className="card mth-card">
        <p className="card-label">5 · Before you order</p>
        <p>
          Every amplicon here spans at least two exons, so a product off contaminating
          genomic DNA either fails or runs visibly longer on a gel. Specificity is
          established within the gene's NM isoform set — for genome-wide uniqueness, run
          the pair through Primer-BLAST as usual. The formulas and constants behind every
          number are on the{" "}
          <button type="button" className="linkish" onClick={onMethod}>Method</button> page.
        </p>
      </section>

      <button type="button" className="mth-back bottom" onClick={onBack}>
        <span className="rev"><ArrowRight /></span>{backLabel}
      </button>
    </div>
  );
}
