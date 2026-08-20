import {
  type KeyboardEvent, type ReactNode, useEffect, useMemo, useRef, useState,
} from "react";
import {
  ARM_GAP, DEFAULT_CONDITIONS, DNTP_MAX, DNTP_MIN, MG_MAX, MG_MIN,
  PRIMER_MAX, PRIMER_MIN, SALT_MAX, SALT_MIN, WALLACE_MAX,
  armTmText, autoPick, evalWindow, isDefaultConditions,
  type TmConditions, type WindowEval,
} from "../lib/tm";
import { revComp } from "../lib/partner";
import { numStr } from "../lib/format";
import { Copy } from "./icons";

/**
 * The parts of the Tm-guided designer that do not care WHERE the junction came from.
 *
 * Two callers share this: JunctionDesigner (junction taken from a RefSeq transcript's
 * exon structure) and CustomJunction (junction taken from two arms the user pasted).
 * Keeping the drag-select strip, the verdict and the metric cards here means both
 * modes run the identical evalWindow/autoPick path — "same calculation" is guaranteed
 * by construction rather than by two implementations that can drift apart.
 */

export type CondKey = keyof TmConditions;

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

const TM_FLOOR = 30, TM_CEIL = 95;

/** Everything the workbench needs to know about the junction it is rendering. */
export interface JunctionGeom {
  /** The sequence to design against (a transcript mRNA, or arm5 + arm3 concatenated). */
  seq: string;
  /** 0-based index of the first base on the 3′ side — the cut sits just before it. */
  jx: number;
  /** 0-based floor/ceiling for arm extension (the flanking exons, or the arms' own ends). */
  leftBound: number;
  rightBound: number;
  /** Slice of `seq` to render (transcripts show a window around the junction). */
  winStart: number;
  winEnd: number;
  /** Per-base text colour class: the 5′ side, the 3′ side, or neither. */
  classOf: (i: number) => "ex-a" | "ex-b" | "ex-o";
  /** Column headings for the two arm metric cards, e.g. "5′ arm (donor)". */
  leftLabel: string;
  rightLabel: string;
}

/**
 * Tm range + reaction conditions, with the type-freely / clamp-on-commit contract both
 * modes use: the numeric value updates live while what is typed is already valid, and
 * is only reconciled against the other bound and the hard limits on blur/Enter.
 */
export function useJunctionSettings() {
  const [tmMin, setTmMin] = useState(60);
  const [tmMax, setTmMax] = useState(65);
  const [minStr, setMinStr] = useState("60");
  const [maxStr, setMaxStr] = useState("65");
  const [cond, setCond] = useState<TmConditions>(DEFAULT_CONDITIONS);
  const [condStr, setCondStr] = useState<Record<CondKey, string>>(() => ({
    saltMM: String(DEFAULT_CONDITIONS.saltMM),
    primerUM: String(DEFAULT_CONDITIONS.primerUM),
    mgMM: String(DEFAULT_CONDITIONS.mgMM),
    dntpMM: String(DEFAULT_CONDITIONS.dntpMM),
  }));
  const saltRef = useRef<HTMLInputElement>(null);

  function editMin(raw: string) {
    setMinStr(raw);
    const v = parseInt(raw, 10);
    if (!Number.isNaN(v) && v >= TM_FLOOR && v < tmMax) setTmMin(v);
  }
  function commitMin() {
    const v = parseInt(minStr, 10);
    const c = Number.isNaN(v) ? tmMin : Math.max(TM_FLOOR, Math.min(v, tmMax - 1));
    setTmMin(c); setMinStr(String(c));
  }
  function editMax(raw: string) {
    setMaxStr(raw);
    const v = parseInt(raw, 10);
    if (!Number.isNaN(v) && v > tmMin && v <= TM_CEIL) setTmMax(v);
  }
  function commitMax() {
    const v = parseInt(maxStr, 10);
    const c = Number.isNaN(v) ? tmMax : Math.min(TM_CEIL, Math.max(v, tmMin + 1));
    setTmMax(c); setMaxStr(String(c));
  }
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

  return {
    tmMin, tmMax, minStr, maxStr, editMin, commitMin, editMax, commitMax,
    cond, condStr, editCond, commitCond, resetConditions, saltRef, focusConditions,
  };
}
export type JunctionSettings = ReturnType<typeof useJunctionSettings>;

const enterBlur = (e: KeyboardEvent<HTMLInputElement>) => {
  if (e.key === "Enter") e.currentTarget.blur();
};

/** The whole-primer Tm range, the arm-cap reminder, and the buffer summary chip. */
/**
 * `showArmCap` off for a CONVENTIONAL pair: the arm rule is what makes a junction primer
 * junction-specific, and a primer that spans no junction has no arms. Showing it there
 * would state a constraint the design does not have.
 */
export function TmRangeControls({ s, showChip = true, showArmCap = true }: {
  s: JunctionSettings; showChip?: boolean; showArmCap?: boolean;
}) {
  return (
    <div className="jd-range">
      <label>{showArmCap ? "Whole-primer Tm" : "Primer Tm"}
        <input type="number" value={s.minStr} min={TM_FLOOR} max={s.tmMax - 1} inputMode="numeric"
          onChange={(e) => s.editMin(e.target.value)} onBlur={s.commitMin} onKeyDown={enterBlur} />
        <span className="dash">–</span>
        <input type="number" value={s.maxStr} min={s.tmMin + 1} max={TM_CEIL} inputMode="numeric"
          onChange={(e) => s.editMax(e.target.value)} onBlur={s.commitMax} onKeyDown={enterBlur} />
        <span className="unit">°C</span>
      </label>
      {showArmCap && (
        <span className="jd-cap">each arm Tm ≤ <b>whole-primer Tm − {ARM_GAP} °C</b></span>
      )}
      {showChip && (
        <button type="button" onClick={s.focusConditions}
          className={`jd-cond-chip ${isDefaultConditions(s.cond) ? "" : "mod"}`}
          title="Edit the buffer and primer concentration used by the Tm formula">
          <b>{numStr(s.cond.saltMM)}</b> mM salt · <b>{numStr(s.cond.mgMM)}</b> mM Mg²⁺ ·{" "}
          <b>{numStr(s.cond.primerUM)}</b> µM primer
          {!isDefaultConditions(s.cond) && <span className="dotmark" aria-hidden="true" />}
        </button>
      )}
    </div>
  );
}

/**
 * The reaction conditions that feed the Tm formula. Editing any of them re-runs the whole
 * designer — warm zone, auto-pick, arm caps, every metric above. The formula itself is
 * documented on the Method page; this card only exposes the knobs.
 */
export function Conditions({ s, onMethod }: { s: JunctionSettings; onMethod?: () => void }) {
  const isDefault = isDefaultConditions(s.cond);
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
                <input ref={k === "saltMM" ? s.saltRef : undefined}
                  type="number" value={s.condStr[k]} inputMode="decimal"
                  min={COND_BOUNDS[k][0]} max={COND_BOUNDS[k][1]} step={step}
                  onChange={(e) => s.editCond(k, e.target.value)}
                  onBlur={() => s.commitCond(k)} onKeyDown={enterBlur} />
                <span className="unit">{unit}</span>
              </span>
            </label>
          ))}
          <button className="btn btn-ghost jd-reset" onClick={s.resetConditions} disabled={isDefault}>
            Reset{isDefault ? "" : " to default"}
          </button>
        </div>
      </div>
    </div>
  );
}

/**
 * The drag-to-select sequence strip plus the live verdict for whatever is selected.
 * `reseedKey` re-seeds the selection from the auto-pick when it changes — callers pass
 * whatever identifies "a different junction" (the transcript, or the typed arms).
 * Deliberately NOT keyed on the reaction conditions: keeping the user's window lets them
 * watch the same primer move as they retune salt or primer concentration.
 */
export function JunctionWorkbench({ geom, s, reseedKey, intro, legendExtra, onEval, role = null }: {
  geom: JunctionGeom;
  s: JunctionSettings;
  reseedKey: string;
  intro?: ReactNode;
  /** Extra control(s) rendered at the end of the legend row (e.g. a full-view toggle). */
  legendExtra?: ReactNode;
  /** Live report of the current selection's evaluation — lets a parent design a partner. */
  onEval?: (ev: WindowEval | null) => void;
  /**
   * Which primer of a pair this junction's oligo is, when that is already decided (a
   * two-junction combo fixes both). Set it and the ordered sequence is labelled FORWARD or
   * REVERSE, so nobody has to infer a direction from an unmarked string of bases.
   *
   * "reverse" additionally flips the ordered oligo to the reverse complement of the selected
   * sense window. Only that sequence flips — the strip, the arms and every Tm stay
   * sense-oriented, and the numbers are unaffected because a duplex melts at the same
   * temperature read from either strand.
   *
   * null when the direction is not fixed here (a single-junction EEJ takes its direction
   * from whichever side its partner primer lands on, which the partner panel states).
   */
  role?: "forward" | "reverse" | null;
}) {
  const reverse = role === "reverse";
  const { seq, jx, leftBound, rightBound, winStart, winEnd } = geom;
  const [sel, setSel] = useState<{ s: number; e: number } | null>(null);
  const [copied, setCopied] = useState(false);
  const dragging = useRef(false);
  const anchor = useRef(0);
  const seqRef = useRef<HTMLDivElement>(null);

  const pick = useMemo(
    () => autoPick(seq, jx, leftBound, rightBound, s.tmMin, s.tmMax, s.cond),
    [seq, jx, leftBound, rightBound, s.tmMin, s.tmMax, s.cond],
  );

  useEffect(() => {
    if (pick?.best) setSel({ s: pick.best.s, e: pick.best.e });
    else setSel(null);
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [reseedKey, s.tmMin, s.tmMax]);

  useEffect(() => {
    const up = () => { dragging.current = false; };
    window.addEventListener("mouseup", up);
    return () => window.removeEventListener("mouseup", up);
  }, []);

  const ev: WindowEval | null = sel
    ? evalWindow(seq, jx, sel.s, sel.e, s.tmMin, s.tmMax, s.cond)
    : null;

  // Report the live evaluation upward whenever the selection or the Tm inputs move.
  useEffect(() => {
    onEval?.(ev);
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [sel?.s, sel?.e, s.tmMin, s.tmMax, s.cond, seq, jx]);

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
    navigator.clipboard?.writeText(orderedOligo(ev, reverse));
    setCopied(true);
    window.setTimeout(() => setCopied(false), 1400);
  }
  function onSeqKeyDown(e: KeyboardEvent) {
    if ((e.metaKey || e.ctrlKey) && (e.key === "c" || e.key === "C")) {
      e.preventDefault();
      copySelected();
    }
  }

  const bases = [];
  for (let i = winStart; i < winEnd; i++) {
    const inSel = sel != null && i >= sel.s && i < sel.e;
    const side = inSel ? (i < jx ? " sel-l" : " sel-r") : "";
    const warm = pick?.warm.has(i) && !inSel ? " warm" : "";
    bases.push(
      <span
        key={i}
        className={`b ${geom.classOf(i)}${side}${warm}${i === jx ? " jx" : ""}`}
        onMouseDown={(e) => { e.preventDefault(); onDown(i); }}
        onMouseEnter={() => onEnter(i)}
      >{seq[i]}</span>
    );
  }

  // The pristine auto-pick is the "suggested" primer (yellow); once the user drags to a
  // different window it becomes their own selection — blue if valid, orange if not.
  const isSuggested = !!(pick?.best && sel && sel.s === pick.best.s && sel.e === pick.best.e);
  const selState = isSuggested ? "suggested" : (ev?.valid ?? false) ? "valid" : "invalid";

  return (
    <>
      {intro}
      <div ref={seqRef} tabIndex={0} role="textbox"
        aria-label="Drag to select a primer; press Cmd or Ctrl + C to copy"
        className={`jd-seq mono ${selState}`}
        onKeyDown={onSeqKeyDown} onMouseLeave={() => { dragging.current = false; }}>
        {/* Strand orientation: the strip always reads sense 5′ (left) → 3′ (right). */}
        <div className="jd-ends" aria-hidden="true"><span>5′</span><span>3′</span></div>
        {winStart > 0 && <span className="cdna-ellipsis">…{winStart} nt </span>}
        {bases}
        {winEnd < seq.length && <span className="cdna-ellipsis"> {seq.length - winEnd} nt…</span>}
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
        <span className="lg"><span className="jx-mark" />junction</span>
        <button className="btn btn-ghost jd-auto"
          onClick={() => pick.best && setSel({ s: pick.best.s, e: pick.best.e })}
          disabled={!pick.best}>Reset to auto-pick</button>
        {legendExtra}
      </div>

      {ev && <Readout ev={ev} tmMin={s.tmMin} tmMax={s.tmMax}
        leftLabel={geom.leftLabel} rightLabel={geom.rightLabel} role={role} />}
    </>
  );
}

/**
 * The oligo as ordered, 5′→3′. For a reverse primer that is the reverse complement of the
 * selected sense window — the sequence highlighted on the strip is the template it binds,
 * not the thing you buy.
 */
export function orderedOligo(ev: WindowEval, reverse: boolean): string {
  return reverse ? revComp(ev.whole.seq) : ev.whole.seq;
}

function Readout({ ev, tmMin, tmMax, leftLabel, rightLabel, role }: {
  ev: WindowEval; tmMin: number; tmMax: number; leftLabel: string; rightLabel: string;
  role: "forward" | "reverse" | null;
}) {
  const reverse = role === "reverse";
  function copy() { navigator.clipboard?.writeText(orderedOligo(ev, reverse)); }
  // Reverse-complementing swaps which arm leads: revComp(left+right) = revComp(right)+revComp(left),
  // so the acceptor arm becomes the oligo's 5′ end and the junction mark moves with it.
  const oligo5 = reverse ? revComp(ev.right.seq) : ev.left.seq;
  const oligo3 = reverse ? revComp(ev.left.seq) : ev.right.seq;
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
        <Metric label={leftLabel} seq={ev.left.seq} tmText={armTmText(ev.left.tm)}
          note={`cap ≤ ${cap} °C · ${ev.left.len} nt`} pass={ev.left.pass} />
        <Metric label={rightLabel} seq={ev.right.seq} tmText={armTmText(ev.right.tm)}
          note={`cap ≤ ${cap} °C · ${ev.right.len} nt`} pass={ev.right.pass} />
      </div>

      <div className="jd-foot">
        <div className="jd-primer mono">
          {role && <span className={`jd-role-tag ${role}`}>{role}</span>}
          5′-{oligo5}<span className="jd-split" />{oligo3}-3′
        </div>
        <button className="btn btn-ghost" onClick={copy}><Copy /> Copy primer</button>
      </div>
      {reverse && (
        <p className="sub jd-rev-note">
          This primer runs <b>reverse</b>, so the oligo to order is the <b>reverse complement</b>
          {" "}of the window highlighted above — that window is the template it binds, not the
          sequence you buy. The Tm figures are unchanged: a duplex melts at the same
          temperature read from either strand.
        </p>
      )}
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
