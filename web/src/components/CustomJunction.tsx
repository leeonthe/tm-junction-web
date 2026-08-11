import { useMemo } from "react";
import { MIN_ARM, gcPercent } from "../lib/tm";
import {
  Conditions, JunctionWorkbench, TmRangeControls, useJunctionSettings,
  type JunctionGeom,
} from "./JunctionWorkbench";

/** Shortest arm we accept as input. Below this there is not enough sequence to design against. */
export const MIN_INPUT_ARM = 10;
/** What we nudge toward — long enough that the arm cap, not the sequence, is the constraint. */
const IDEAL_ARM = 20;

/** GAPDH exon 1–2, the same junction the transcript designer opens on. */
export const ARM_EXAMPLE = { five: "TTTTGCGTCGCCAG", three: "CCGAGCCACATCGCTCAGAC" };

/** The two arms as the user typed them, before cleaning. */
export interface Arms { five: string; three: string }
export const EMPTY_ARMS: Arms = { five: "", three: "" };

/**
 * Normalise pasted sequence: drop anything that is not a letter (spaces, line breaks,
 * FASTA position numbers), upper-case it, and read RNA as DNA so a pasted mRNA works.
 */
function clean(raw: string): string {
  return raw.toUpperCase().replace(/[^A-Z]/g, "").replace(/U/g, "T");
}

interface ArmState {
  seq: string;
  /** Letters that survived cleaning but are not bases — N, R, Y and friends. */
  bad: string[];
  tooShort: boolean;
  ok: boolean;
}

export function readArm(raw: string): ArmState {
  const seq = clean(raw);
  const bad = [...new Set(seq.replace(/[ACGT]/g, ""))];
  const tooShort = seq.length > 0 && seq.length < MIN_INPUT_ARM;
  return { seq, bad, tooShort, ok: seq.length >= MIN_INPUT_ARM && bad.length === 0 };
}

/**
 * The paste-two-arms input, shown in the hero in place of the search bar when the user
 * picks the "Custom sequence" mode. Validation lives here, next to each field; the design
 * itself renders below in CustomJunctionResult.
 */
export function CustomArmInputs({ arms, onChange, loading }: {
  arms: Arms;
  onChange: (a: Arms) => void;
  loading?: boolean;
}) {
  const five = readArm(arms.five);
  const three = readArm(arms.three);
  return (
    <div className="cj-entry">
      <div className="cj-inputs">
        <ArmField label="5′ arm" hint="upstream / donor side" arm={five} value={arms.five}
          onChange={(v) => onChange({ ...arms, five: v })} disabled={loading} />
        <span className="cj-join" aria-hidden="true">
          <span className="cj-join-mark" />
          <span className="cj-join-cap">junction</span>
        </span>
        <ArmField label="3′ arm" hint="downstream / acceptor side" arm={three} value={arms.three}
          onChange={(v) => onChange({ ...arms, three: v })} disabled={loading} />
      </div>
      <div className="cj-actions">
        <button type="button" className="btn btn-ghost" onClick={() => onChange(ARM_EXAMPLE)}>
          Load example
        </button>
        <button type="button" className="btn btn-ghost" onClick={() => onChange(EMPTY_ARMS)}
          disabled={!arms.five && !arms.three}>Clear</button>
        <span className="cj-note">
          ~{IDEAL_ARM} nt per arm works best; {MIN_INPUT_ARM} nt is the minimum. A primer still
          needs ≥ {MIN_ARM} nt on each side of the cut to count as junction-spanning.
        </span>
      </div>
    </div>
  );
}

/**
 * Tm-guided EEJ designer for a junction the user typed in, rather than one derived from a
 * RefSeq transcript. The two arms are concatenated into a synthetic "mRNA" whose junction
 * sits exactly where they meet, then handed to the same JunctionWorkbench the transcript
 * designer uses — so the Tm rule, the warm zone and the arm caps are identical, and no
 * second implementation exists to drift.
 */
export default function CustomJunctionResult({ arms, onMethod }: {
  arms: Arms;
  onMethod?: () => void;
}) {
  const s = useJunctionSettings();
  const five = readArm(arms.five);
  const three = readArm(arms.three);
  const ready = five.ok && three.ok;

  const geom = useMemo((): JunctionGeom | null => {
    if (!ready) return null;
    const seq = five.seq + three.seq;
    const jx = five.seq.length;
    return {
      seq, jx,
      leftBound: 0,                 // the pasted arms ARE the bounds — nothing beyond them
      rightBound: seq.length,
      winStart: 0,
      winEnd: seq.length,
      classOf: (i) => (i < jx ? "ex-a" : "ex-b"),
      leftLabel: "5′ arm",
      rightLabel: "3′ arm",
    };
  }, [ready, five.seq, three.seq]);

  if (!geom) {
    const started = five.seq || three.seq;
    return (
      <section className="card cj-idle">
        <p className="card-label" style={{ marginBottom: 10 }}>Tm-guided junction designer</p>
        <p className="sub">
          {started
            ? "Both arms need to be valid before the designer can run — see the notes under each box."
            : `Paste a 5′ and a 3′ arm above (${MIN_INPUT_ARM} nt or more each). The whole primer must melt in your range while each arm alone stays at least 15 °C below it, so it primes only across that exact junction.`}
        </p>
      </section>
    );
  }

  return (
    <section className="card elevated jd">
      <div className="card-head">
        <div style={{ display: "flex", alignItems: "center", gap: 12 }}>
          <p className="card-label" style={{ margin: 0 }}>Tm-guided junction designer</p>
          <span className="jd-badge neutral">
            {five.seq.length} + {three.seq.length} nt · GC {gcPercent(geom.seq).toFixed(0)}%
          </span>
        </div>
        <TmRangeControls s={s} />
      </div>

      <JunctionWorkbench
        geom={geom} s={s} reseedKey={geom.seq}
        intro={
          <p className="sub jd-intro">
            Drag across the junction to select a primer. The 5′ arm is{" "}
            <b className="jd-exa-t">magenta</b>, the 3′ arm <b className="jd-exb-t">green</b>.
            Shaded bases are where an arm still sits under its cap — a selection there can be
            valid. This is the same Tm rule the transcript designer applies.
          </p>
        }
      />

      <Conditions s={s} onMethod={onMethod} />
    </section>
  );
}

function ArmField({ label, hint, arm, value, onChange, disabled }: {
  label: string; hint: string; arm: ArmState;
  value: string; onChange: (v: string) => void; disabled?: boolean;
}) {
  const state = arm.bad.length || arm.tooShort ? "bad" : arm.ok ? "ok" : "";
  return (
    <label className={`cj-arm ${state}`}>
      <span className="cj-arm-head">
        <b>{label}</b><span className="cj-arm-hint">{hint}</span>
        {arm.seq.length > 0 && (
          <span className="cj-arm-len mono">
            {arm.seq.length} nt · GC {gcPercent(arm.seq).toFixed(0)}%
          </span>
        )}
      </span>
      <textarea className="mono" value={value} rows={2} spellCheck={false} disabled={disabled}
        placeholder={label === "5′ arm" ? "e.g. TTTTGCGTCGCCAG" : "e.g. CCGAGCCACATCGCTCAGAC"}
        onChange={(e) => onChange(e.target.value)} />
      {arm.bad.length > 0
        ? <span className="cj-arm-msg bad">Not a DNA base: {arm.bad.join(", ")}</span>
        : arm.tooShort
          ? <span className="cj-arm-msg bad">{arm.seq.length} nt — needs at least {MIN_INPUT_ARM}.</span>
          : <span className="cj-arm-msg">Spaces, line breaks and numbers are ignored; U reads as T.</span>}
    </label>
  );
}
