import { useEffect, useRef, useState } from "react";
import type { GeneLookupResponse } from "../lib/types";

/**
 * Gene-name result: a reference alignment of every NM transcript for the gene, so the
 * user can choose which variant to analyze. Deliberately tier-less — no Conventional /
 * Needs-EEJ / Hard-case coloring here; that comes only after a variant is analyzed.
 */
export default function GeneTranscriptPicker({
  data, onSelect, busy,
}: {
  data: GeneLookupResponse;
  onSelect: (accession: string) => void;
  busy?: boolean;
}) {
  const { gene, transcripts } = data;
  const [hover, setHover] = useState<number | null>(null);

  const scrollRef = useRef<HTMLDivElement>(null);
  const MIN_W = 720;
  const [W, setW] = useState(1000);
  const [scrollable, setScrollable] = useState(false);
  useEffect(() => {
    const el = scrollRef.current;
    if (!el) return;
    const check = () => {
      setW(Math.max(MIN_W, Math.round(el.clientWidth)));
      setScrollable(el.clientWidth < MIN_W);
    };
    check();
    const ro = new ResizeObserver(check);
    ro.observe(el);
    return () => ro.disconnect();
  }, []);

  const padL = 196, padR = 150, rowH = 46, top = 16, exH = 16;

  let gmin = Infinity, gmax = -Infinity;
  for (const t of transcripts) for (const e of t.exons) {
    if (e.begin < gmin) gmin = e.begin;
    if (e.end > gmax) gmax = e.end;
  }
  const span = gmax - gmin || 1;
  const x = (p: number) => padL + ((p - gmin) / span) * (W - padL - padR);
  const axisY = top + transcripts.length * rowH + 8;
  const height = axisY + 26;
  const ticks = [gmin, (gmin + gmax) / 2, gmax];

  return (
    <section className="card">
      <div className="card-head">
        <div>
          <h3 className="card-title">
            {gene.symbol} <span className="gp-sub">· {transcripts.length} NM transcripts</span>
          </h3>
          <p className="sub">{gene.description} · GRCh38 chr{gene.chromosome || "?"} · Gene {gene.gene_id}</p>
        </div>
        <span className="gp-hint">Reference only — pick a variant to analyze</span>
      </div>

      <div ref={scrollRef} className="scroll-x" style={{ cursor: scrollable ? "grab" : "default" }}>
        <svg width={W} height={height} role="img"
          aria-label={`Exon alignment of ${transcripts.length} ${gene.symbol} transcripts`}>
          {transcripts.map((t, i) => {
            const cy = top + i * rowH + rowH / 2;
            const first = t.exons[0], last = t.exons[t.exons.length - 1];
            const on = hover === i;
            return (
              <g key={t.accession} className="gp-row" role="button" tabIndex={0}
                aria-label={`Analyze ${t.accession}`}
                onMouseEnter={() => setHover(i)} onMouseLeave={() => setHover((h) => (h === i ? null : h))}
                onClick={() => !busy && onSelect(t.accession)}
                onKeyDown={(e) => { if ((e.key === "Enter" || e.key === " ") && !busy) { e.preventDefault(); onSelect(t.accession); } }}
                style={{ cursor: busy ? "wait" : "pointer" }}>
                {/* row hit-area + hover highlight */}
                <rect x={8} y={top + i * rowH + 2} width={W - 16} height={rowH - 4} rx={10}
                  fill={on ? "var(--brand-tint)" : "transparent"}
                  stroke={on ? "color-mix(in srgb,var(--brand) 32%,transparent)" : "transparent"} />
                <text x={20} y={cy + 4} fontFamily="var(--mono)" fontSize={12.5}
                  fontWeight={on ? 700 : 500} fill={on ? "var(--brand-ink)" : "var(--text)"}>
                  {t.accession}
                </text>
                {t.is_mane && (
                  <>
                    <rect x={20 + t.accession.length * 7.1 + 6} y={cy - 9} width={38} height={15} rx={4} fill="var(--brand-tint)" />
                    <text x={20 + t.accession.length * 7.1 + 9} y={cy + 2} fontSize={9.5} fontWeight={700} fill="var(--brand-ink)">MANE</text>
                  </>
                )}
                {/* intron connector + neutral exon blocks */}
                <line x1={x(first.begin)} y1={cy} x2={x(last.end)} y2={cy} stroke="var(--border-2)" strokeWidth={1.5} />
                {t.exons.map((e, k) => (
                  <rect key={k} x={x(e.begin)} y={cy - exH / 2}
                    width={Math.max(2, x(e.end) - x(e.begin))} height={exH} rx={2.5}
                    fill={on ? "var(--pick-exon-on)" : "var(--pick-exon)"} />
                ))}
                {/* right meta (idle) swaps to the analyze affordance on hover — same anchor, no overlap */}
                {on
                  ? <text x={W - 16} y={cy + 4} fontSize={11.5} fontWeight={700} textAnchor="end"
                      fill="var(--brand-ink)">Analyze →</text>
                  : <text x={W - 16} y={cy + 4} fontSize={11.5} fontFamily="var(--mono)"
                      fill="var(--faint)" textAnchor="end">
                      {t.length.toLocaleString()} nt · {t.exon_count} ex
                    </text>}
              </g>
            );
          })}
          <line x1={padL} y1={axisY} x2={x(gmax)} y2={axisY} stroke="var(--border)" strokeWidth={1} />
          {ticks.map((p, i) => (
            <g key={i}>
              <line x1={x(p)} y1={axisY} x2={x(p)} y2={axisY + 4} stroke="var(--faint)" />
              <text x={x(p)} y={axisY + 17} fontSize={10.5} fontFamily="var(--mono)"
                textAnchor="middle" fill="var(--faint)">{(p / 1e6).toFixed(3)} Mb</text>
            </g>
          ))}
        </svg>
      </div>
      <p className="g-note">Exons drawn to GRCh38 scale and aligned across isoforms. Click a transcript to design its primers.</p>
    </section>
  );
}
