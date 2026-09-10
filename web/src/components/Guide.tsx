import { useMemo, useState } from "react";
import { ArrowRight } from "./icons";
import { JunctionWorkbench, useJunctionSettings, type JunctionGeom } from "./JunctionWorkbench";
import { BracketGlyph, TargetGlyph, buildCells } from "./GraphCard";

/**
 * The drag demo's junction is real: GAPDH NM_001357943.2, exon 3 | exon 4 — the junction
 * the designer itself recommends for that transcript — 48 nt each side, taken from the
 * seeded RefSeq sequence. A made-up sequence would behave differently from anything the
 * user later meets; this one behaves exactly like the page they are being taught. This
 * junction replaced the earlier exon 1 | 2 one, whose GC cliff auto-picked a primer with
 * uncomfortably short arms: here the pick is a textbook 22-mer with 12 / 10 nt arms.
 */
const DEMO_ARM5 = "GTGGATATTGTTGCCATCAATGACCCCTTCATTGACCTCAACTACATG";
const DEMO_ARM3 = "GCTGAGAACGGGAAGCTTGTCATCAATGGAAATCCCATCACCATCTTC";

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

/**
 * One primer-design strategy, drawn the way the exon graph draws it: tier-coloured exons,
 * yellow primer target sites, a red bracket for a single EEJ, magenta for a double. The
 * five figures below are the full vocabulary of the graph's annotations — a reader who can
 * parse these five can parse any transcript the tool shows them.
 */
function StrategyFig({ yellows = [], brackets = [], tier = "--eej" }: {
  yellows?: number[];
  brackets?: { from: number; to: number; combo?: boolean }[];
  tier?: string;
}) {
  const X = [12, 92, 172, 252, 332, 412];      // 6 exons, fixed layout
  const W = [56, 36, 42, 36, 42, 66];
  const mid = (i: number) => X[i] + W[i] / 2;
  return (
    <svg className="strat-fig" viewBox="0 0 490 46" aria-hidden="true">
      <line x1={8} y1={30} x2={484} y2={30} stroke="var(--border-2)" strokeWidth={1.5} />
      {X.map((x, i) => (
        <rect key={i} x={x} y={22} width={W[i]} height={16} rx={3}
          fill={yellows.includes(i) ? "var(--amp-pair)" : `var(${tier})`} />
      ))}
      {brackets.map((b, k) => (
        <path key={k}
          d={`M${mid(b.from)} 18 L${mid(b.from)} 9 L${mid(b.to)} 9 L${mid(b.to)} 18`}
          fill="none" stroke={b.combo ? "var(--eej-combo)" : "var(--eej-single)"}
          strokeWidth={2.5} strokeLinecap="round" strokeLinejoin="round" />
      ))}
    </svg>
  );
}

/** The five annotations, in the order of how much the design asks of the transcript. */
const STRATEGIES: { name: string; fig: React.ReactNode; what: string }[] = [
  {
    name: "Single primer target site",
    fig: <StrategyFig tier="--conv" yellows={[2]} />,
    what: "One exon region (yellow) is unique to this transcript. Either primer of an "
      + "ordinary pair covers it; the other sits wherever the amplicon needs.",
  },
  {
    name: "Double primer target site",
    fig: <StrategyFig tier="--conv" yellows={[1, 4]} />,
    what: "No single exon is unique, but this combination of two exons is — no sibling "
      + "carries both. Forward goes in one, reverse in the other, and only this "
      + "transcript can form the product.",
  },
  {
    name: "Single EEJ",
    fig: <StrategyFig brackets={[{ from: 2, to: 3 }]} />,
    what: "One exon–exon junction (red bracket) is unique. A primer spanning it — the "
      + "junction designer's job — fires only on this transcript; its partner is an "
      + "ordinary primer nearby.",
  },
  {
    name: "Single EEJ + primer target site",
    fig: <StrategyFig brackets={[{ from: 1, to: 2 }]} yellows={[4]} />,
    what: "Neither the junction nor the exon region is unique alone, but no sibling has "
      + "both. The EEJ primer spans the bracket and its partner must sit in the yellow "
      + "region — the pair is the specificity.",
  },
  {
    name: "Double EEJ",
    fig: <StrategyFig brackets={[{ from: 0, to: 1, combo: true }, { from: 3, to: 4, combo: true }]} />,
    what: "Only a combination of two junctions (magenta brackets) isolates this "
      + "transcript. Both primers must be junction-spanning — one across each.",
  },
];

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
        <p className="card-label">3 · The five design strategies</p>
        <p>
          The verdict decides <em>which kind</em> of design isolates your transcript, and the
          exon graph annotates it. These five marks are the whole vocabulary:
        </p>
        <div className="strat-list">
          {STRATEGIES.map((s) => (
            <div className="strat" key={s.name}>
              {s.fig}
              <p><b>{s.name}</b> — {s.what}</p>
            </div>
          ))}
        </div>
        <p className="mth-note">
          Yellow always means "a primer must cover this"; a bracket always means "a primer
          must span this junction". Red marks a single required junction, magenta a
          two-junction combination.
        </p>
      </section>

      <section className="card mth-card">
        <p className="card-label">4 · Take the primers</p>
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
          This one is live — GAPDH's exon&nbsp;3–exon&nbsp;4 junction (NM_001357943.2),
          the real designer.
          Drag across the letters and watch the verdict and the three Tm cards react;
          drag a window that sits all on one side and see <em>why</em> it fails:
        </p>
        <DragDemo />
        <p className="mth-note">
          Both the <b>Primer Tm</b> range and the <b>Reaction conditions</b> (Na⁺/K⁺, Mg²⁺,
          dNTPs, primer concentration) are yours to change, in every designer's settings
          panel — every Tm on the page recalculates from what you set. Match them to your
          own assay and master mix before trusting the numbers; Reset returns the defaults.
        </p>
      </section>

      <section className="card mth-card">
        <p className="card-label">5 · Gene-wide tools</p>
        <p>
          The <i>Gene classification</i> tab shows the whole gene: every isoform's exon
          structure, its tier, and its per-isoform verdict — click any row to design for
          that transcript instead. The <i>Whole transcript amplification</i> tab is the opposite job:
          one primer pair that amplifies as many of the gene's variants as possible at a
          single product size, for measuring total expression rather than one isoform.
        </p>
      </section>

      <section className="card mth-card">
        <p className="card-label">6 · Before you order</p>
        <p>
          Every amplicon here spans at least two exons, so a product off contaminating
          genomic DNA either fails or runs visibly longer on a gel. Specificity is
          established within the gene's NM isoform set — for genome-wide uniqueness, run
          the pair through NCBI Primer-BLAST or UCSC's BLAT (Genome Browser → Tools →
          Blat) as usual. The formulas and constants behind every
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
