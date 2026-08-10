import {
  type KeyboardEvent, type ReactNode, type RefObject, useEffect, useMemo, useRef, useState,
} from "react";
import type { Exon, TranscriptVerdict } from "../lib/types";
import {
  ARM_GAP, DEFAULT_CONDITIONS, DNTP_MAX, DNTP_MIN, MG_MAX, MG_MIN,
  PRIMER_MAX, PRIMER_MIN, SALT_MAX, SALT_MIN, WALLACE_MAX,
  armTmText, autoPick, evalWindow, isDefaultConditions,
  type TmConditions, type WindowEval,
} from "../lib/tm";
import { numStr } from "../lib/format";
import { Copy } from "./icons";

type CondKey = keyof TmConditions;

const COND_BOUNDS: Record<CondKey, [number, number]> = {
  saltMM: [SALT_MIN, SALT_MAX],
  primerUM: [PRIMER_MIN, PRIMER_MAX],
  mgMM: [MG_MIN, MG_MAX],
  dntpMM: [DNTP_MIN, DNTP_MAX],
};

/** The four editable buffer terms, in the order they appear under the designer. */
const COND_FIELDS: { k: CondKey; label: ReactNode; unit: string; step: number }[] = [
  { k: "saltMM", label: <>Monovalent salt <i>[Na⁺]+[K⁺]</i></>, unit: "mM", step: 5 },
  { k: "mgMM", label: <>Magnesium <i>[Mg²⁺]</i></>, unit: "mM", step: 0.5 },
  { k: "dntpMM", label: <>dNTPs <i>total</i></>, unit: "mM", step: 0.1 },
  { k: "primerUM", label: <>Primer conc. <i>C<sub>T</sub></i></>, unit: "µM", step: 0.05 },
];

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
export default function JunctionDesigner({ mrna, verdict, onMethod }: {
  mrna: string;
  verdict: TranscriptVerdict;
  onMethod?: () => void;
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
  // Reaction conditions feeding the Tm formula. The defaults are the project's
  // calibration point, so the designer opens with exactly the numbers it always had.
  const [cond, setCond] = useState<TmConditions>(DEFAULT_CONDITIONS);
  const [condStr, setCondStr] = useState<Record<CondKey, string>>(() => ({
    saltMM: String(DEFAULT_CONDITIONS.saltMM),
    primerUM: String(DEFAULT_CONDITIONS.primerUM),
    mgMM: String(DEFAULT_CONDITIONS.mgMM),
    dntpMM: String(DEFAULT_CONDITIONS.dntpMM),
  }));
  const [sel, setSel] = useState<{ s: number; e: number } | null>(null);
  const [copied, setCopied] = useState(false);
  const dragging = useRef(false);
  const anchor = useRef(0);
  const seqRef = useRef<HTMLDivElement>(null);
  const saltRef = useRef<HTMLInputElement>(null);

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
    return autoPick(mrna, geom.jx, geom.leftBound, geom.rightBound, tmMin, tmMax, cond);
  }, [geom, mrna, tmMin, tmMax, cond]);

  // 0-based mRNA index → exon index, for text coloring
  const exonAt = useMemo(() => {
    const arr = new Int16Array(mrna.length).fill(-1);
    exons.forEach((e: Exon, i) => { for (let p = e.tx_begin - 1; p < e.tx_end; p++) arr[p] = i; });
    return arr;
  }, [exons, mrna.length]);

  // Re-seed the selection from the auto-pick when the target or the Tm window changes.
  // Deliberately NOT on a conditions change: `pick` is already recomputed by then, but
  // keeping the user's window lets them watch the same primer move as they retune salt
  // or primer concentration.
  useEffect(() => {
    if (pick?.best) setSel({ s: pick.best.s, e: pick.best.e });
    else setSel(null);
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [geom, mrna, tmMin, tmMax]);

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
  // The two exons flanking THIS junction get dedicated colors (exon A = 5′/donor = magenta,
  // exon B = 3′/acceptor = dark green); any other exon that peeks into the window is muted.
  const donorIdx = exons.findIndex((e) => e.order === donor.order);
  const acceptorIdx = exons.findIndex((e) => e.order === acceptor.order);

  const ev: WindowEval | null = sel ? evalWindow(mrna, jx, sel.s, sel.e, tmMin, tmMax, cond) : null;

  // render window centered on the junction, comfortably covering the arms
  const flank = 46;
  const winStart = Math.max(0, jx - flank);
  const winEnd = Math.min(mrna.length, jx + flank);

  function onDown(i: number) {
    dragging.current = true;
    anchor.current = i;
    setSel({ s: i, e: i + 1 });
    seqRef.current?.focus();   // preventDefault on the span blocks auto-focus, so focus here
  }
  function onEnter(i: number) {
    if (!dragging.current) return;
    const a = anchor.current;
    setSel({ s: Math.min(a, i), e: Math.max(a, i) + 1 });
  }
  // Copy the selected primer (5′→3′) on ⌘C / Ctrl+C while the sequence box is focused —
  // the box has user-select:none (so dragging designs a primer), so native copy won't work.
  function copySelected() {
    if (!ev) return;
    navigator.clipboard?.writeText(ev.whole.seq);
    setCopied(true);
    window.setTimeout(() => setCopied(false), 1400);
  }
  function onSeqKeyDown(e: KeyboardEvent) {
    if ((e.metaKey || e.ctrlKey) && (e.key === "c" || e.key === "C")) {
      e.preventDefault();
      copySelected();
    }
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

  // Reaction-condition inputs: same type-freely / clamp-on-commit contract as the Tm
  // range. All four accept decimals (0.25 µM primer, 1.5 mM Mg²⁺, 62.5 mM salt are
  // all real setups).
  function editCond(k: CondKey, raw: string) {
    setCondStr((s) => ({ ...s, [k]: raw }));
    const v = parseFloat(raw);
    const [lo, hi] = COND_BOUNDS[k];
    if (Number.isFinite(v) && v >= lo && v <= hi) setCond((c) => ({ ...c, [k]: v }));
  }
  function commitCond(k: CondKey) {
    const v = parseFloat(condStr[k]);
    const [lo, hi] = COND_BOUNDS[k];
    const c = Number.isFinite(v) ? Math.min(hi, Math.max(lo, v)) : cond[k];
    setCond((p) => ({ ...p, [k]: c }));
    setCondStr((s) => ({ ...s, [k]: numStr(c) }));
  }
  function resetConditions() {
    setCond(DEFAULT_CONDITIONS);
    setCondStr({
      saltMM: String(DEFAULT_CONDITIONS.saltMM),
      primerUM: String(DEFAULT_CONDITIONS.primerUM),
      mgMM: String(DEFAULT_CONDITIONS.mgMM),
      dntpMM: String(DEFAULT_CONDITIONS.dntpMM),
    });
  }
  function focusConditions() {
    saltRef.current?.scrollIntoView({ block: "center", behavior: "smooth" });
    saltRef.current?.focus();
    saltRef.current?.select();
  }

  const bases = [];
  for (let i = winStart; i < winEnd; i++) {
    const ex = exonAt[i];
    const exc = ex === donorIdx ? " ex-a" : ex === acceptorIdx ? " ex-b" : " ex-o";
    const inSel = sel != null && i >= sel.s && i < sel.e;
    const side = inSel ? (i < jx ? " sel-l" : " sel-r") : "";
    const warm = pick?.warm.has(i) && !inSel ? " warm" : "";
    const jxm = i === jx ? " jx" : "";
    bases.push(
      <span
        key={i}
        className={`b${exc}${side}${warm}${jxm}`}
        onMouseDown={(e) => { e.preventDefault(); onDown(i); }}
        onMouseEnter={() => onEnter(i)}
      >{mrna[i]}</span>
    );
  }

  const valid = ev?.valid ?? false;
  // The pristine auto-pick is the "suggested" primer (yellow); once the user drags to a
  // different window it becomes their own selection — blue if valid, orange if not.
  const isSuggested = !!(pick?.best && sel && sel.s === pick.best.s && sel.e === pick.best.e);
  const selState = isSuggested ? "suggested" : valid ? "valid" : "invalid";

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
          <span className="jd-cap">each arm Tm ≤ <b>whole-primer Tm − {ARM_GAP} °C</b></span>
          <button type="button" onClick={focusConditions}
            className={`jd-cond-chip ${isDefaultConditions(cond) ? "" : "mod"}`}
            title="Edit the buffer and primer concentration used by the Tm formula">
            <b>{numStr(cond.saltMM)}</b> mM salt · <b>{numStr(cond.mgMM)}</b> mM Mg²⁺ ·{" "}
            <b>{numStr(cond.primerUM)}</b> µM primer
            {!isDefaultConditions(cond) && <span className="dotmark" aria-hidden="true" />}
          </button>
        </div>
      </div>

      <p className="sub" style={{ marginBottom: 14 }}>
        Drag across the junction to select a primer. The 5′ arm sits on <b className="jd-exa-t">exon {donor.order}</b>,
        the 3′ arm on <b className="jd-exb-t">exon {acceptor.order}</b>. A valid primer means the whole primer melts
        in range while neither arm alone is stable enough to prime — so it fires only on this exact splice.
      </p>

      <div ref={seqRef} tabIndex={0} role="textbox" aria-label="Drag to select a primer; press Cmd or Ctrl + C to copy"
        className={`jd-seq mono ${selState}`}
        onKeyDown={onSeqKeyDown} onMouseLeave={() => { dragging.current = false; }}>
        {winStart > 0 && <span className="cdna-ellipsis">…{winStart} nt </span>}
        {bases}
        {winEnd < mrna.length && <span className="cdna-ellipsis"> {mrna.length - winEnd} nt…</span>}
      </div>

      <div className="jd-copyhint">
        {copied
          ? <span className="ok">✓ Copied primer to clipboard</span>
          : <>Drag to select · press <kbd>⌘C</kbd> / <kbd>Ctrl+C</kbd> to copy the selected primer</>}
      </div>

      <div className="jd-legend">
        <span className="lg"><span className="sw warm" />selectable region (each arm ≤ cap)</span>
        <span className="lg"><span className="sw sel-sug" />suggested primer</span>
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

      <Conditions
        cond={cond} saltRef={saltRef} condStr={condStr}
        editCond={editCond} commitCond={commitCond}
        reset={resetConditions} onMethod={onMethod}
      />
    </section>
  );
}

/**
 * The reaction conditions that feed the Tm formula. Editing any of them re-runs the whole
 * designer — warm zone, auto-pick, arm caps, every metric above. The formula itself is
 * documented on the Method page; this card only exposes the knobs.
 */
function Conditions({ cond, saltRef, condStr, editCond, commitCond, reset, onMethod }: {
  cond: TmConditions;
  saltRef: RefObject<HTMLInputElement>;
  condStr: Record<CondKey, string>;
  editCond: (k: CondKey, v: string) => void;
  commitCond: (k: CondKey) => void;
  reset: () => void;
  onMethod?: () => void;
}) {
  const isDefault = isDefaultConditions(cond);
  const enterBlur = (e: KeyboardEvent<HTMLInputElement>) => {
    if (e.key === "Enter") e.currentTarget.blur();
  };

  return (
    <div className="jd-model">
      <div className="jd-model-head">
        <div>
          <p className="card-label" style={{ margin: 0 }}>Tm formula · reaction conditions</p>
          <p className="jd-model-sub">
            Set your buffer and primer concentration — the whole-primer Tm is recomputed from
            them. Mg²⁺ and the monovalent cations compete for the DNA backbone, so both matter,
            and dNTPs chelate Mg²⁺ so only the surplus counts. Arms under {WALLACE_MAX} nt use
            the Wallace rule, which has no salt or concentration term.
            {onMethod && <> The formula is documented in{" "}
              <button type="button" className="linkish" onClick={onMethod}>Method</button>.</>}
          </p>
        </div>
        <div className="jd-cond">
          {COND_FIELDS.map(({ k, label, unit, step }) => (
            <label key={k}>
              <span className="ck">{label}</span>
              <span className="cin">
                <input ref={k === "saltMM" ? saltRef : undefined}
                  type="number" value={condStr[k]} inputMode="decimal"
                  min={COND_BOUNDS[k][0]} max={COND_BOUNDS[k][1]} step={step}
                  onChange={(e) => editCond(k, e.target.value)}
                  onBlur={() => commitCond(k)} onKeyDown={enterBlur} />
                <span className="unit">{unit}</span>
              </span>
            </label>
          ))}
          <button className="btn btn-ghost jd-reset" onClick={reset} disabled={isDefault}>
            Reset{isDefault ? "" : " to default"}
          </button>
        </div>
      </div>
    </div>
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
        <Metric label="Whole primer" seq={ev.whole.seq} tmText={`${ev.whole.tm.toFixed(1)} °C`}
          note={`target ${tmMin}–${tmMax} °C · GC ${ev.whole.gc.toFixed(0)}% · ${ev.whole.len} nt`}
          pass={ev.whole.pass} />
        <Metric label="5′ arm (donor)" seq={ev.left.seq} tmText={armTmText(ev.left.tm)}
          note={`cap ≤ ${cap} °C · ${ev.left.len} nt`} pass={ev.left.pass} />
        <Metric label="3′ arm (acceptor)" seq={ev.right.seq} tmText={armTmText(ev.right.tm)}
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

function Metric({ label, seq, tmText, note, pass }: {
  label: string; seq: string; tmText: string; note: string; pass: boolean;
}) {
  return (
    <div className={`jd-metric ${pass ? "pass" : "fail"}`}>
      <div className="jd-metric-head">
        <span className="k">{label}</span>
        <span className={`tag ${pass ? "ok" : "bad"}`}>{pass ? "pass" : "fail"}</span>
      </div>
      <div className="jd-tm mono">{seq ? tmText : "—"}</div>
      <div className="jd-note">{note}</div>
    </div>
  );
}
