import { useEffect, useMemo, useRef, useState } from "react";
import type { Exon, TranscriptVerdict } from "../lib/types";
import { ARM_GAP, autoPick, evalWindow, type WindowEval } from "../lib/tm";
import { Copy } from "./icons";

/**
 * Interactive Tm-guided EEJ primer designer.
 *
 * The user sets a whole-primer Tm range; the tool renders the cDNA around the
 * discriminating exon–exon junction, shades the "warm zone" (where a junction-
 * spanning window melts in range), and drops a green auto-picked primer. The
 * user can drag-select any window across the junction to see it validated live:
 * whole Tm in [min,max] AND each arm ≤ max − 15 °C → green, else red with the
 * exact numbers so they can see where they are. All Tm is computed client-side
 * (lib/tm.ts) — the same formula the engine uses.
 */
export default function JunctionDesigner({ mrna, verdict }: {
  mrna: string;
  verdict: TranscriptVerdict;
}) {
  const exons = verdict.exons;
  const junction = verdict.recommended_junction;

  const [tmMin, setTmMin] = useState(60);
  const [tmMax, setTmMax] = useState(65);
  // Raw input text, so a mid-edit value (empty, single digit) isn't clamped on every
  // keystroke — we only reconcile against the other bound on blur/Enter. The numeric
  // tmMin/tmMax update live while what's typed is already a valid, in-range number.
  const [minStr, setMinStr] = useState("60");
  const [maxStr, setMaxStr] = useState("65");
  const [sel, setSel] = useState<{ s: number; e: number } | null>(null);
  const dragging = useRef(false);
  const anchor = useRef(0);

  const geom = useMemo(() => {
    if (!junction) return null;
    const donor = exons.find((e) => e.order === junction.donor_order);
    const acceptor = exons.find((e) => e.order === junction.acceptor_order);
    if (!donor || !acceptor) return null;
    return {
      donor, acceptor,
      jx: donor.tx_end,               // 0-based index of the first acceptor base
      leftBound: donor.tx_begin - 1,  // 0-based start of the donor exon (arm floor)
      rightBound: acceptor.tx_end,    // 0-based exclusive end of the acceptor exon
    };
  }, [junction, exons]);

  const pick = useMemo(() => {
    if (!geom) return null;
    return autoPick(mrna, geom.jx, geom.leftBound, geom.rightBound, tmMin, tmMax);
  }, [geom, mrna, tmMin, tmMax]);

  // 0-based mRNA index → exon index, for text coloring
  const exonAt = useMemo(() => {
    const arr = new Int16Array(mrna.length).fill(-1);
    exons.forEach((e: Exon, i) => { for (let p = e.tx_begin - 1; p < e.tx_end; p++) arr[p] = i; });
    return arr;
  }, [exons, mrna.length]);

  // reset the selection to the auto-pick whenever the range / junction changes
  useEffect(() => {
    if (pick?.best) setSel({ s: pick.best.s, e: pick.best.e });
    else setSel(null);
  }, [pick]);

  useEffect(() => {
    const up = () => { dragging.current = false; };
    window.addEventListener("mouseup", up);
    return () => window.removeEventListener("mouseup", up);
  }, []);

  if (!geom || !junction) {
    return (
      <section className="card">
        <p className="card-label" style={{ marginBottom: 12 }}>Tm-guided junction designer</p>
        <p className="sub">
          This designer applies to junction-spanning (EEJ) primers. The current target has
          no single discriminating exon–exon junction, so there is nothing to tune here.
        </p>
      </section>
    );
  }

  const { jx, donor, acceptor } = geom;

  const ev: WindowEval | null = sel ? evalWindow(mrna, jx, sel.s, sel.e, tmMin, tmMax) : null;

  // render window centered on the junction, comfortably covering the arms
  const flank = 46;
  const winStart = Math.max(0, jx - flank);
  const winEnd = Math.min(mrna.length, jx + flank);

  function onDown(i: number) {
    dragging.current = true;
    anchor.current = i;
    setSel({ s: i, e: i + 1 });
  }
  function onEnter(i: number) {
    if (!dragging.current) return;
    const a = anchor.current;
    setSel({ s: Math.min(a, i), e: Math.max(a, i) + 1 });
  }

  // Tm range inputs: type freely; live-update only when the value is already valid,
  // and reconcile against the other bound + hard limits [30,95] on blur/Enter.
  function editMin(raw: string) {
    setMinStr(raw);
    const v = parseInt(raw, 10);
    if (!Number.isNaN(v) && v >= 30 && v < tmMax) setTmMin(v);
  }
  function commitMin() {
    const v = parseInt(minStr, 10);
    const c = Number.isNaN(v) ? tmMin : Math.max(30, Math.min(v, tmMax - 1));
    setTmMin(c);
    setMinStr(String(c));
  }
  function editMax(raw: string) {
    setMaxStr(raw);
    const v = parseInt(raw, 10);
    if (!Number.isNaN(v) && v > tmMin && v <= 95) setTmMax(v);
  }
  function commitMax() {
    const v = parseInt(maxStr, 10);
    const c = Number.isNaN(v) ? tmMax : Math.min(95, Math.max(v, tmMin + 1));
    setTmMax(c);
    setMaxStr(String(c));
  }

  const bases = [];
  for (let i = winStart; i < winEnd; i++) {
    const ex = exonAt[i];
    const inSel = sel != null && i >= sel.s && i < sel.e;
    const side = inSel ? (i < jx ? " sel-l" : " sel-r") : "";
    const warm = pick?.warm.has(i) && !inSel ? " warm" : "";
    const jxm = i === jx ? " jx" : "";
    bases.push(
      <span
        key={i}
        className={`b${side}${warm}${jxm}`}
        style={{ color: `var(--exon-${(ex < 0 ? 0 : ex) % 5})` }}
        onMouseDown={(e) => { e.preventDefault(); onDown(i); }}
        onMouseEnter={() => onEnter(i)}
      >{mrna[i]}</span>
    );
  }

  const valid = ev?.valid ?? false;

  return (
    <section className="card elevated jd">
      <div className="card-head">
        <div style={{ display: "flex", alignItems: "center", gap: 12 }}>
          <p className="card-label" style={{ margin: 0 }}>Tm-guided junction designer</p>
          <span className="jd-badge">exon {donor.order}–{acceptor.order} junction</span>
        </div>
        <div className="jd-range">
          <label>Whole-primer Tm
            <input type="number" value={minStr} min={30} max={tmMax - 1} inputMode="numeric"
              onChange={(e) => editMin(e.target.value)} onBlur={commitMin}
              onKeyDown={(e) => { if (e.key === "Enter") e.currentTarget.blur(); }} />
            <span className="dash">–</span>
            <input type="number" value={maxStr} min={tmMin + 1} max={95} inputMode="numeric"
              onChange={(e) => editMax(e.target.value)} onBlur={commitMax}
              onKeyDown={(e) => { if (e.key === "Enter") e.currentTarget.blur(); }} />
            <span className="unit">°C</span>
          </label>
          <span className="jd-cap">each arm ≤ <b>whole-primer Tm − {ARM_GAP} °C</b></span>
        </div>
      </div>

      <p className="sub" style={{ marginBottom: 14 }}>
        Drag across the junction to select a primer. The 5′ arm sits on exon {donor.order},
        the 3′ arm on exon {acceptor.order}. Green = the whole primer melts in range and
        neither arm alone is stable enough to prime — so it fires only on this exact splice.
      </p>

      <div className={`jd-seq mono ${valid ? "valid" : "invalid"}`} onMouseLeave={() => { dragging.current = false; }}>
        {winStart > 0 && <span className="cdna-ellipsis">…{winStart} nt </span>}
        {bases}
        {winEnd < mrna.length && <span className="cdna-ellipsis"> {mrna.length - winEnd} nt…</span>}
      </div>

      <div className="jd-legend">
        <span className="lg"><span className="sw warm" />whole-primer Tm in range (warm zone)</span>
        <span className="lg"><span className="sw sel-ok" />selected — valid</span>
        <span className="lg"><span className="sw sel-bad" />selected — invalid</span>
        <span className="lg"><span className="jx-mark" />exon–exon junction</span>
        {pick && (
          <button className="btn btn-ghost jd-auto"
            onClick={() => pick.best && setSel({ s: pick.best.s, e: pick.best.e })}
            disabled={!pick.best}>Reset to auto-pick</button>
        )}
      </div>

      {ev && <Readout ev={ev} tmMin={tmMin} tmMax={tmMax} />}
    </section>
  );
}

function Readout({ ev, tmMin, tmMax }: {
  ev: WindowEval; tmMin: number; tmMax: number;
}) {
  function copy() { navigator.clipboard?.writeText(ev.whole.seq); }
  // cap tracks the selection's actual whole-primer Tm (ev.armCap = whole Tm − ARM_GAP)
  const cap = ev.armCap.toFixed(1);
  return (
    <div className="jd-readout">
      <div className={`jd-verdict ${ev.valid ? "ok" : "bad"}`}>
        {ev.valid
          ? <span><b>✓ Valid EEJ primer</b> — the whole primer melts in range, but each arm
              alone stays below the cap, so it primes only across this exact junction.</span>
          : <div><b>✗ Not valid yet</b><ul>{ev.reasons.map((r) => <li key={r}>{r}</li>)}</ul></div>}
      </div>

      {ev.valid && ev.notes.length > 0 && (
        <div className="jd-notes">
          {ev.notes.map((n) => <div key={n} className="jd-note-line">⚠ {n}</div>)}
        </div>
      )}

      <div className="jd-metrics">
        <Metric label="Whole primer" seq={ev.whole.seq} tmv={ev.whole.tm}
          note={`target ${tmMin}–${tmMax} °C · GC ${ev.whole.gc.toFixed(0)}% · ${ev.whole.len} nt`}
          pass={ev.whole.pass} />
        <Metric label="5′ arm (donor)" seq={ev.left.seq} tmv={ev.left.tm}
          note={`cap ≤ ${cap} °C · ${ev.left.len} nt`} pass={ev.left.pass} />
        <Metric label="3′ arm (acceptor)" seq={ev.right.seq} tmv={ev.right.tm}
          note={`cap ≤ ${cap} °C · ${ev.right.len} nt`} pass={ev.right.pass} />
      </div>

      <div className="jd-foot">
        <div className="jd-primer mono">
          5′-{ev.left.seq}<span className="jd-split" />{ev.right.seq}-3′
        </div>
        <button className="btn btn-ghost" onClick={copy}><Copy /> Copy primer</button>
      </div>
    </div>
  );
}

function Metric({ label, seq, tmv, note, pass }: {
  label: string; seq: string; tmv: number; note: string; pass: boolean;
}) {
  return (
    <div className={`jd-metric ${pass ? "pass" : "fail"}`}>
      <div className="jd-metric-head">
        <span className="k">{label}</span>
        <span className={`tag ${pass ? "ok" : "bad"}`}>{pass ? "pass" : "fail"}</span>
      </div>
      <div className="jd-tm mono">{seq ? `${tmv.toFixed(1)} °C` : "—"}</div>
      <div className="jd-note">{note}</div>
    </div>
  );
}
