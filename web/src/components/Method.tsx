import type { ReactNode } from "react";
import {
  ARM_GAP, DEFAULT_CONDITIONS, DNTP_MAX, DNTP_MIN, LEN_MAX, LEN_MIN, MG_DNTP_KA, MG_MAX,
  MG_MIN, MIN_ARM, PRIMER_MAX, PRIMER_MIN, RATIO_MIXED_MAX, RATIO_MONO_MAX, R_GAS,
  SALT_MAX, SALT_MIN, WALLACE_MAX, ZERO_C, primerMolar, saltMolar, tmParts,
  wallaceTm,
} from "../lib/tm";
import { fixed, molarStr, numStr, rounded, signed } from "../lib/format";
import { ArrowRight } from "./icons";

/** The oligo the new formula was validated against — see the worked example below. */
const EXAMPLE = "AACTACATGGCTGAGAAC";
/** An arm-length oligo for the Wallace worked example — also recomputed live. */
const ARM_EXAMPLE = "CTCGCGA";

/**
 * Method page: how the verdict is reached, what the Tm-guided EEJ rule is, and the exact
 * melting-temperature formula behind every number in the designer. Reachable from the nav.
 * Everything numeric here is computed by lib/tm.ts at the default conditions, so the page
 * cannot drift from the engine of the designer.
 */
export default function Method({ onBack, backLabel }: { onBack: () => void; backLabel: string }) {
  const ex = tmParts(EXAMPLE, DEFAULT_CONDITIONS)!;
  const armGC = [...ARM_EXAMPLE].filter((c) => c === "G" || c === "C").length;
  const armAT = ARM_EXAMPLE.length - armGC;

  return (
    <div className="mth">
      <button type="button" className="mth-back" onClick={onBack}>
        <span className="rev"><ArrowRight /></span>{backLabel}
      </button>

      <header className="mth-head">
        <p className="card-label" style={{ margin: 0 }}>Method</p>
        <h1>How TmJunction decides, designs, and computes Tm</h1>
        <p className="mth-lede">
          Every verdict on this site is sequence-based: the target mRNA is compared against the
          other RefSeq <b>NM</b> isoforms of the same gene, and a primer site counts as specific
          only if it is absent from all of them. Nothing here is a homology guess — it is exact
          substring evidence over the isoform set.
        </p>
      </header>

      <section className="card mth-card">
        <p className="card-label">1 · The verdict</p>
        <p>
          A primer-length window is slid along the target mRNA. A window that appears in no
          sibling isoform is a target-specific primer site, which lands the transcript in one of
          three tiers:
        </p>
        <ul className="mth-list">
          <li>
            <b className="t-conv-ink">Conventional</b> — some window sits inside a single exon.
            An ordinary primer already discriminates the target, so no junction trick is needed.
          </li>
          <li>
            <b className="t-eej-ink">Needs EEJ</b> — no exon-internal window is unique, but one
            exon–exon junction is. The primer must straddle that splice; designing it is what the
            junction designer does.
          </li>
          <li>
            <b className="t-hard-ink">Hard case</b> — neither exists. No single primer pair can
            separate this isoform, and a junction-combination strategy is required.
          </li>
        </ul>
        <p className="mth-note">
          Specificity is established within the gene's NM isoform set, not genome-wide.
        </p>
      </section>

      <section className="card mth-card">
        <p className="card-label">2 · The Tm-guided junction rule</p>
        <p>
          A junction-spanning primer has two arms: the 5′ arm on the donor exon, the 3′ arm on
          the acceptor exon. The unspliced pre-mRNA and any isoform that skips the junction carry
          only <em>one</em> of those arms. So specificity is a thermodynamic question, not a
          sequence one — the primer must be stable across the junction and unstable on either arm
          alone:
        </p>
        <div className="mth-rule">
          <div className="mth-rule-row">
            <span className="k">whole primer</span>
            <span className="v mono">Tm within the range you set (default 60–65 °C)</span>
          </div>
          <div className="mth-rule-row">
            <span className="k">each arm</span>
            <span className="v mono">Tm ≤ whole-primer Tm − {ARM_GAP} °C</span>
          </div>
          <div className="mth-rule-row">
            <span className="k">geometry</span>
            <span className="v mono">
              ≥ {MIN_ARM} nt per arm · {LEN_MIN}–{LEN_MAX} nt whole primer
            </span>
          </div>
        </div>
        <p>
          The {ARM_GAP} °C gap is fixed, matching the worked example in the paper, and is measured
          from the primer's <em>own</em> whole-primer Tm rather than from the top of your range.
          The bound is upper-only: a colder arm is more specific, so it is never rejected for
          being too cold. Both arms below the cap means neither half can prime on its own, so the
          primer fires only where the two exons are actually joined.
        </p>
      </section>

      <section className="card mth-card">
        <p className="card-label">3 · The melting-temperature formula</p>
        <p>
          Tm is computed in two steps. Nearest-neighbour thermodynamics — SantaLucia (1998)
          unified stacking parameters — gives the melting temperature the duplex would have in
          1 M Na⁺. A salt term then moves it to your actual buffer. Because the salt effect is
          linear in <em>1/Tm</em> rather than in Tm, the correction is applied to the reciprocal:
        </p>

        <div className="mth-eq-wrap">
          <div className="mth-eq mono">
            <span>T<sub>m</sub><sup>1M</sup></span><span className="op">=</span>
            <span className="frac">
              <span className="num">ΔH<sub>total</sub></span>
              <span className="den">ΔS<sub>total</sub> + R · ln(C<sub>T</sub>)</span>
            </span>
          </div>
        </div>
        <div className="mth-eq-wrap">
          <div className="mth-eq mono">
            <span className="frac">
              <span className="num">1</span><span className="den">T<sub>m</sub></span>
            </span>
            <span className="op">=</span>
            <span className="frac">
              <span className="num">1</span><span className="den">T<sub>m</sub><sup>1M</sup></span>
            </span>
            {/* operator + operand kept in one span so a narrow screen never wraps them apart */}
            <span><span className="op">+ </span>Δ<sub>salt</sub></span>
            <span><span className="op">, </span>T<sub>m</sub>(°C) = T<sub>m</sub> − {ZERO_C}</span>
          </div>
        </div>

        <div className="mth-terms">
          <Term k={<>ΔH<sub>total</sub></>} unit="kcal/mol (cal/mol in the equation)"
            note="Sum of the adjacent-doublet enthalpies plus the initiation term, which depends on whether each end is a G·C or an A·T pair." />
          <Term k={<>ΔS<sub>total</sub></>} unit="cal/(mol·K)"
            note="Sum of the adjacent-doublet entropies plus the same terminal-dependent initiation." />
          <Term k="R" unit={`${R_GAS} cal/(K·mol)`}
            note="Ideal gas constant." />
          <Term k={<>C<sub>T</sub></>} unit="mol/L"
            note={<>Primer concentration — no /4 term, because in PCR the primer is in vast
              excess over the template it anneals to. You enter it in µM; 1 µM = 10⁻⁶ M.</>} />
          <Term k={<>Δ<sub>salt</sub></>} unit="K⁻¹"
            note={<>The Owczarzy (2008) salt term, set by monovalent salt, free Mg²⁺, GC
              fraction and length — see below.</>} />
          <Term k={<>{ZERO_C}</>} unit="K"
            note="Converts the Kelvin result to °C." />
        </div>

        <h3 className="mth-h3">The salt term: monovalent and divalent cations compete</h3>
        <p>
          K⁺/Na⁺ and Mg²⁺ both stabilise the duplex by screening the phosphate backbone, and they
          compete for it — so neither can be folded into the other. Which model applies is decided
          by their ratio:
        </p>
        <div className="mth-eq-wrap">
          <div className="mth-eq mono">
            <span>R</span><span className="op">=</span>
            <span className="frac">
              <span className="num">√[Mg²⁺]<sub>free</sub></span>
              <span className="den">[Mon⁺]</span>
            </span>
          </div>
        </div>
        <div className="mth-rule">
          <div className="mth-rule-row">
            <span className="k mono">R &lt; {RATIO_MONO_MAX}</span>
            <span className="v">Monovalent dominates — the Owczarzy (2004) [Na⁺]-only equation.</span>
          </div>
          <div className="mth-rule-row">
            <span className="k mono">{RATIO_MONO_MAX} ≤ R &lt; {RATIO_MIXED_MAX}</span>
            <span className="v">
              They compete — the full seven-coefficient equation, with three of the coefficients
              themselves refitted as functions of ln[Mon⁺]. <b>An ordinary PCR or qPCR buffer sits
              here</b> ({DEFAULT_CONDITIONS.saltMM} mM K⁺ with {DEFAULT_CONDITIONS.mgMM} mM Mg²⁺
              gives R ≈ {fixed(ex.salt.ratio, 2)}).
            </span>
          </div>
          <div className="mth-rule-row">
            <span className="k mono">R ≥ {RATIO_MIXED_MAX}</span>
            <span className="v">Mg²⁺ dominates — the same equation with its base coefficients.</span>
          </div>
        </div>
        <p>
          In every branch Δ<sub>salt</sub> depends on the GC fraction and the oligo length as well
          as the ion concentrations, which is why two primers in the same tube get different
          corrections. dNTPs are handled first: they chelate Mg²⁺ with an association constant of{" "}
          <span className="mono">{MG_DNTP_KA.toExponential(0)} M⁻¹</span>, so the model uses the
          <b> free</b> Mg²⁺ left over, not what you pipetted. At the defaults that is{" "}
          <b className="mono">{fixed(ex.salt.mgFree * 1e3, 2)} mM</b> free of{" "}
          {numStr(DEFAULT_CONDITIONS.mgMM)} mM total.
        </p>

        <h3 className="mth-h3">Reaction conditions</h3>
        <p>
          All four inputs are editable in the junction designer: <b>monovalent salt</b> [Na⁺]+[K⁺]
          in mM ({SALT_MIN}–{SALT_MAX}), <b>magnesium</b> [Mg²⁺] in mM ({MG_MIN}–{MG_MAX}),{" "}
          <b>dNTPs</b> in mM ({DNTP_MIN}–{DNTP_MAX}), and <b>primer concentration</b> C
          <sub>T</sub> in µM ({PRIMER_MIN}–{PRIMER_MAX}). The designer opens on a standard
          qPCR/RT-PCR buffer — <b className="mono">{numStr(DEFAULT_CONDITIONS.saltMM)} mM</b> salt
          ({molarStr(saltMolar(DEFAULT_CONDITIONS.saltMM))}),{" "}
          <b className="mono">{numStr(DEFAULT_CONDITIONS.mgMM)} mM</b> Mg²⁺,{" "}
          <b className="mono">{numStr(DEFAULT_CONDITIONS.dntpMM)} mM</b> dNTPs and{" "}
          <b className="mono">{numStr(DEFAULT_CONDITIONS.primerUM)} µM</b> primer (
          {molarStr(primerMolar(DEFAULT_CONDITIONS.primerUM))}) — so the numbers you see are the
          ones your reaction will actually run at.
        </p>

        <h3 className="mth-h3">Worked example</h3>
        <p>
          At the default conditions, <span className="mono">{EXAMPLE}</span> evaluates as:
        </p>
        <div className="mth-eq-wrap">
          <div className="mth-eq mono">
            <span>T<sub>m</sub><sup>1M</sup></span><span className="op">=</span>
            <span className="frac">
              <span className="num">{rounded(ex.dh * 1000)}</span>
              <span className="den">{fixed(ex.ds, 1)} + ({fixed(ex.dsTerm, 2)})</span>
            </span>
            <span><span className="op">= </span>{fixed(ex.tm1M + ZERO_C, 1)} K</span>
            <span><span className="op">= </span>{fixed(ex.tm1M, 1)} °C</span>
          </div>
        </div>
        <div className="mth-eq-wrap">
          <div className="mth-eq mono">
            <span>T<sub>m</sub></span><span className="op">=</span>
            <span>(1 / {fixed(ex.tm1M + ZERO_C, 1)} {signed(ex.salt.delta * 1e5)}·10⁻⁵)⁻¹</span>
            <span><span className="op">− </span>{ZERO_C}</span>
            <span><span className="op">= </span>
              <span className="res">{fixed(ex.tm, 1)} °C</span>
            </span>
          </div>
        </div>
        <p className="mth-note">
          ΔH<sub>total</sub> = {fixed(ex.dh, 1)} kcal/mol and ΔS<sub>total</sub> ={" "}
          {fixed(ex.ds, 1)} cal/(mol·K) for this {EXAMPLE.length}-mer; {fixed(ex.dsTerm, 2)} is the
          concentration term at {numStr(DEFAULT_CONDITIONS.primerUM)} µM. The salt term is
          evaluated in the <b>{ex.salt.regime}</b> branch and pulls Tm by{" "}
          {signed(ex.saltShift, 1)} °C. IDT's OligoAnalyzer reports 56 °C for this oligo in qPCR
          mode and 49 °C in monovalent-only mode; this engine gives{" "}
          <b className="mono">{fixed(ex.tm, 1)}</b> and{" "}
          <b className="mono">{fixed(tmParts(EXAMPLE, { ...DEFAULT_CONDITIONS, mgMM: 0, dntpMM: 0 })!.tm, 1)}</b> °C.
        </p>

        <h3 className="mth-h3">Short arms: the Wallace rule</h3>
        <p>
          Nearest-neighbour thermodynamics is not valid on very short oligos: the initiation and
          concentration terms stop being small next to the stacking sum, so a {MIN_ARM}–8 nt arm
          comes out implausibly cold and can even go negative. Arms are short by design here — the
          arm cap is what makes them short — so any arm below{" "}
          <b className="mono">{WALLACE_MAX} nt</b> is scored with the Wallace rule instead:
        </p>
        <div className="mth-eq-wrap">
          <div className="mth-eq mono">
            <span>T<sub>m</sub></span><span className="op">=</span>
            <span>2 °C · (A + T)</span>
            <span><span className="op">+ </span>4 °C · (G + C)</span>
          </div>
        </div>
        <p className="mth-note">
          A + T and G + C are base counts. For the arm{" "}
          <span className="mono">{ARM_EXAMPLE}</span>: 2 · {armAT} + 4 · {armGC} ={" "}
          <b className="mono">{wallaceTm(ARM_EXAMPLE)} °C</b>. The rule has no salt or
          concentration term, so the reaction-condition fields do not move arm Tm below{" "}
          {WALLACE_MAX} nt. Arms of {WALLACE_MAX} nt or longer stay on the nearest-neighbour
          formula above, and the <b>whole primer always does</b>, whatever its length. An arm
          value that still comes out below zero is shown as <span className="mono">NA</span>{" "}
          rather than as a number.
        </p>

        <h3 className="mth-h3">Scope of this number</h3>
        <p className="mth-note">
          Nearest-neighbour with the Owczarzy salt correction is what the junction designer uses
          for every whole-primer number, with the Wallace rule covering short arms, so its figures
          are internally consistent and are what your Tm range should be read against. It is not
          the same estimator as the primer-QC figures on the conventional primer cards, which come
          from primer3 when it is available. The model is fitted to fully complementary B-form DNA
          duplexes; it does not account for secondary structure, mismatches, or modified bases.
        </p>
        <p className="mth-note">
          References: SantaLucia (1998) <i>PNAS</i> 95:1460 for the nearest-neighbour parameters;
          Owczarzy et al. (2004) <i>Biochemistry</i> 43:3537 for the monovalent correction; and
          Owczarzy et al. (2008) <i>Biochemistry</i> 47:5336 for the magnesium and mixed-salt
          model. This is the same combination IDT's OligoAnalyzer uses.
        </p>
      </section>

      <section className="card mth-card">
        <p className="card-label">4 · Data</p>
        <p className="mth-note">
          Transcript sequences, exon coordinates, and the isoform set come from NCBI RefSeq via
          the Datasets v2 API, for the assembly shown with each result.
        </p>
      </section>

      <button type="button" className="mth-back bottom" onClick={onBack}>
        <span className="rev"><ArrowRight /></span>{backLabel}
      </button>
    </div>
  );
}

function Term({ k, unit, note }: { k: ReactNode; unit: string; note: ReactNode }) {
  return (
    <div className="mth-term">
      <div className="tk mono">{k}</div>
      <div className="tu">{unit}</div>
      <div className="tn">{note}</div>
    </div>
  );
}
