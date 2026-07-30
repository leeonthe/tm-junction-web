"""Tm-guided primer design — primer3 QC hybrid.

Design flow:
  1. amplify.py finds WHERE specificity lives (unique exon-internal region / unique
     junction) and the exact positions of every target-unique k-mer window.
  2. Here we sweep candidate oligos (lengths 18–25) over that region/junction, keep only
     those that fully CONTAIN a unique window — which guarantees the oligo is absent from
     every sibling (if the oligo appeared in a sibling, so would its window). No per-
     candidate sequence search needed; specificity is structural.
  3. primer3 provides the QC: thermodynamic Tm (SantaLucia), hairpin, self-dimer and
     cross-dimer ΔG/Tm. We gate candidates on Tm/GC/3'-clamp/structure and pick the best.
  4. ΔTm specificity margin = Tm(full match) − Tm(best off-target contiguous match).

primer3 handles primer *quality*; it does not know about isoform specificity — that stays
our layer. See vault `01 Science/Primer Design Strategy.md`.
"""

from __future__ import annotations

import math
from bisect import bisect_left
from dataclasses import dataclass, field

from .amplify import AmplifyResult, Interval, cumulative_exon_ends

try:
    import primer3  # type: ignore
    HAS_PRIMER3 = True
except Exception:  # pragma: no cover - primer3 is a declared dependency
    HAS_PRIMER3 = False

# ---- design constants (see Open Questions Q14) ----
LEN_MIN, LEN_MAX, LEN_OPT = 18, 25, 20
TM_MIN, TM_MAX, TM_OPT = 57.0, 63.0, 60.0
TM_MIN_JUNCTION = 52.0                 # junction-spanning primers run cooler; allow lower
GC_MIN, GC_MAX = 40.0, 60.0
STRUCT_TM_MAX = 45.0                   # reject hairpin / dimer melting above this
AMPLICON_MIN, AMPLICON_MAX, AMPLICON_OPT = 70, 220, 140
MIN_JX_ARM = 5                         # min bases each side of a junction
DELTA_TM_SAFE = 10.0
SHORTLIST = 60                         # structure-check only the best N by prelim score

_COMP = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}

# nearest-neighbor fallback (only if primer3 is unavailable) -------------------
_NN_H = {"AA": -7.9, "TT": -7.9, "AT": -7.2, "TA": -7.2, "CA": -8.5, "TG": -8.5, "GT": -8.4,
         "AC": -8.4, "CT": -7.8, "AG": -7.8, "GA": -8.2, "TC": -8.2, "CG": -10.6, "GC": -9.8,
         "GG": -8.0, "CC": -8.0}
_NN_S = {"AA": -22.2, "TT": -22.2, "AT": -20.4, "TA": -21.3, "CA": -22.7, "TG": -22.7,
         "GT": -22.4, "AC": -22.4, "CT": -21.0, "AG": -21.0, "GA": -22.2, "TC": -22.2,
         "CG": -27.2, "GC": -24.4, "GG": -19.9, "CC": -19.9}


def revcomp(s: str) -> str:
    return "".join(_COMP.get(c, "N") for c in reversed(s))


def gc_percent(s: str) -> float:
    return 100.0 * sum(c in "GC" for c in s) / len(s) if s else 0.0


def _nn_tm(p: str) -> float:
    if len(p) < 2 or any(c not in "ACGT" for c in p):
        return 0.0
    dh, ds = 0.2, -5.7
    for i in range(len(p) - 1):
        dh += _NN_H[p[i:i + 2]]
        ds += _NN_S[p[i:i + 2]]
    return (dh * 1000) / (ds + 1.987 * math.log(0.5e-6 / 4)) - 273.15 + 16.6 * math.log10(0.05)


def tm(p: str) -> float:
    if not p or any(c not in "ACGT" for c in p):
        return 0.0
    return round(primer3.calc_tm(p), 1) if HAS_PRIMER3 else round(_nn_tm(p), 1)


def hairpin_tm(p: str) -> float:
    return round(primer3.calc_hairpin(p).tm, 1) if HAS_PRIMER3 else 0.0


def homodimer_tm(p: str) -> float:
    return round(primer3.calc_homodimer(p).tm, 1) if HAS_PRIMER3 else 0.0


def heterodimer_tm(a: str, b: str) -> float:
    return round(primer3.calc_heterodimer(a, b).tm, 1) if HAS_PRIMER3 else 0.0


def _homopolymer(p: str) -> bool:
    return any(b * 5 in p for b in "ACGT")


def _lcs_substring(a: str, b: str) -> str:
    """Longest common contiguous substring (the strongest spurious binding stretch)."""
    if not a or not b:
        return ""
    prev = [0] * (len(b) + 1)
    best_len, best_end = 0, 0
    for i in range(1, len(a) + 1):
        cur = [0] * (len(b) + 1)
        ai = a[i - 1]
        for j in range(1, len(b) + 1):
            if ai == b[j - 1]:
                cur[j] = prev[j - 1] + 1
                if cur[j] > best_len:
                    best_len, best_end = cur[j], i
        prev = cur
    return a[best_end - best_len:best_end]


def best_offtarget_tm(primer: str, sibling_seqs: list[str]) -> float:
    """Tm of the longest stretch of the primer also present in any sibling mRNA."""
    best = 0.0
    for s in sibling_seqs:
        sub = _lcs_substring(primer, s.upper())
        if len(sub) >= 8:
            best = max(best, tm(sub))
    return round(best, 1)


# ---- candidate evaluation ---------------------------------------------------

@dataclass
class _Eval:
    tm: float
    gc: float
    hairpin: float
    homodimer: float
    ok: bool
    quality: float


def _evaluate(seq: str, tm_min: float) -> _Eval:
    t, g = tm(seq), gc_percent(seq)
    clamp = seq[-1] in "GC"
    hp, hd = hairpin_tm(seq), homodimer_tm(seq)
    ok = (tm_min <= t <= TM_MAX and GC_MIN <= g <= GC_MAX and clamp
          and hp < STRUCT_TM_MAX and hd < STRUCT_TM_MAX and not _homopolymer(seq))
    quality = (-abs(t - TM_OPT) - max(0.0, g - 55) * 0.3 - max(0.0, 45 - g) * 0.3
               - max(0.0, hp - STRUCT_TM_MAX) * 0.6 - max(0.0, hd - STRUCT_TM_MAX) * 0.6
               + (1.0 if clamp else -2.0))
    return _Eval(t, g, hp, hd, ok, quality)


def _prelim(seq: str) -> float:
    return -abs(tm(seq) - TM_OPT) - abs(gc_percent(seq) - 50) * 0.2


def _specific(sorted_starts: list[int], k: int, s: int, length: int) -> bool:
    """True iff oligo [s, s+length) fully contains at least one unique window [p, p+k)."""
    hi = s + length - k
    if hi < s:
        return False
    idx = bisect_left(sorted_starts, s)
    return idx < len(sorted_starts) and sorted_starts[idx] <= hi


@dataclass
class Primer:
    seq: str
    tm: float
    gc: float
    length: int
    kind: str                 # "conventional" | "junction_spanning"
    role: str                 # "forward" | "reverse"
    anchor: str
    hairpin_tm: float = 0.0
    homodimer_tm: float = 0.0
    qc_pass: bool = True
    tx_start: int = 0         # 0-based start of the primer's binding site in the mRNA


@dataclass
class PrimerDesign:
    tier: str
    mechanism: str
    forward: Primer | None = None
    reverse: Primer | None = None
    amplicon_len: int | None = None
    delta_tm: float | None = None
    confidence: str = "n/a"
    excluded_siblings: list[str] = field(default_factory=list)
    flags: list[str] = field(default_factory=list)
    tm_method: str = "primer3" if HAS_PRIMER3 else "nearest-neighbor"
    pair_dimer_tm: float | None = None


def _best_candidate(cands: list[tuple[int, str]], tm_min: float, sibs: list[str],
                    delta_weight: float) -> tuple[int, str, _Eval, float] | None:
    """cands = [(start, seq)]. Shortlist by prelim QC, then rank by QC + specificity.

    `delta_weight` scales how much the ΔTm specificity margin drives the pick — small for
    conventional primers (already highly specific), large for junction primers where the
    discriminating split is the whole point. Returns (start, seq, eval, delta_tm).
    """
    if not cands:
        return None
    shortlist = sorted(cands, key=lambda c: _prelim(c[1]), reverse=True)[:SHORTLIST]
    scored = []
    for s, w in shortlist:
        ev = _evaluate(w, tm_min)
        d = round(tm(w) - best_offtarget_tm(w, sibs), 1)
        score = ev.quality + delta_weight * min(d, 20.0)
        scored.append((ev.ok, score, s, w, ev, d))
    scored.sort(key=lambda x: (x[0], x[1]), reverse=True)
    ok, score, s, w, ev, d = scored[0]
    return s, w, ev, d


def _specific_forward(seq: str, starts: list[int], k: int) -> list[tuple[int, str]]:
    """Oligos (start, seq) that contain a unique window — one per (start,length)."""
    if not starts:
        return []
    ss = sorted(starts)
    lo = max(0, ss[0] + k - LEN_MAX)
    hi = min(len(seq) - LEN_MIN, ss[-1])
    out: list[tuple[int, str]] = []
    for s in range(lo, hi + 1):
        for L in range(LEN_MIN, LEN_MAX + 1):
            if s + L > len(seq):
                break
            if _specific(ss, k, s, L):
                out.append((s, seq[s:s + L]))
    return out


def _partner(seq: str, lo: int, hi: int) -> list[tuple[int, str]]:
    """Non-specific partner-primer windows (start, seq) in [lo, hi]."""
    out: list[tuple[int, str]] = []
    hi = min(hi, len(seq) - LEN_MIN)
    for r in range(max(0, lo), hi + 1):
        for L in range(LEN_MIN, LEN_MAX + 1):
            if r + L > len(seq):
                break
            out.append((r, seq[r:r + L]))
    return out


def design(target_exons: list[Interval], target_seq: str,
           sibling_seqs: dict[str, str], amp: AmplifyResult) -> PrimerDesign:
    seq = target_seq.upper()
    sibs = [s.upper() for s in sibling_seqs.values()]
    excluded = list(sibling_seqs.keys())
    cum = cumulative_exon_ends(target_exons)
    k = amp.k

    if amp.tier == "NO_SINGLE_UNIQUE_JUNCTION":
        return PrimerDesign(
            tier=amp.tier,
            mechanism="No single unique region or junction — needs a junction-combination "
                      "(dual-junction) strategy.",
            flags=["NO_SINGLE_UNIQUE_JUNCTION"],
        )

    # Conventional by structure (rule 7c: not a trimmed copy of any sibling) but with no
    # single unique exon window — a single primer can't be placed. Its unique *exon
    # combination* is amplifiable by a conventional primer PAIR spanning it; designing that
    # pair is a follow-up. Report honestly rather than crash.
    if amp.tier == "CONVENTIONAL" and not amp.unique_regions:
        return PrimerDesign(
            tier=amp.tier,
            mechanism="Structurally unique (not a trimmed copy of any sibling), but no single "
                      "unique exon window — amplifiable by a conventional primer pair spanning "
                      "its unique exon combination (pair design is a follow-up).",
            flags=["NO_UNIQUE_WINDOW"],
        )

    # NEEDS_EEJ via a combination (junction+exon or two junctions) — no single unique junction,
    # so the single-EEJ designer below doesn't apply. Report the combination honestly.
    if amp.tier == "NEEDS_EEJ" and not amp.unique_junctions:
        if amp.combo_je:
            d, a, e = amp.combo_je
            mech = (f"Amplifiable by a combination: an EEJ primer across exon {d}–exon {a} plus a "
                    f"conventional primer in exon {e} — no other isoform has both (combo design is a follow-up).")
        elif amp.combo_jj:
            (d1, a1), (d2, a2) = amp.combo_jj
            mech = (f"Amplifiable by two EEJ primers — across exon {d1}–exon {a1} and exon {d2}–exon {a2} — "
                    f"which together isolate this isoform (combo design is a follow-up).")
        else:
            mech = "Amplifiable only by a junction combination (combo design is a follow-up)."
        return PrimerDesign(tier=amp.tier, mechanism=mech, flags=["COMBO_EEJ"])

    flags: list[str] = []

    if amp.tier == "CONVENTIONAL":
        region = max(amp.unique_regions, key=lambda r: r.window_count)
        starts = amp.internal_starts.get(region.exon_order, [])
        spec_cands = _specific_forward(seq, starts, k)
        chosen = _best_candidate(spec_cands, TM_MIN, sibs, delta_weight=0.3)
        if chosen is None:
            return _fallback(amp, seq, region, starts, k, sibs, excluded)
        s_start, s_seq, s_eval, dtm = chosen
        specific_is_forward = region.side != "reverse"
        specific = _mk_primer(s_seq, s_eval, "conventional",
                              "forward" if specific_is_forward else "reverse",
                              f"exon {region.exon_order} (unique)",
                              as_reverse=not specific_is_forward, pos=s_start)
        mech = f"Conventional primer in exon {region.exon_order} ({region.window_count} unique sites)."
    else:
        j = amp.unique_junctions[0]
        jstarts = amp.junction_starts.get((j.donor_order, j.acceptor_order), [])
        spec_cands = _specific_forward(seq, jstarts, k)  # junction windows cross the boundary
        # junction primers: specificity margin dominates the pick
        chosen = _best_candidate(spec_cands, TM_MIN_JUNCTION, sibs, delta_weight=2.0)
        if chosen is None:
            boundary = cum[j.donor_order - 1]
            s_start = max(0, boundary - LEN_OPT // 2)
            s_seq = seq[s_start:s_start + LEN_OPT]
            s_eval = _evaluate(s_seq, TM_MIN_JUNCTION)
            dtm = round(tm(s_seq) - best_offtarget_tm(s_seq, sibs), 1)
            flags.append("LOW_QC")
        else:
            s_start, s_seq, s_eval, dtm = chosen
        specific_is_forward = True
        specific = _mk_primer(s_seq, s_eval, "junction_spanning", "forward",
                              f"exon {j.donor_order}–exon {j.acceptor_order} junction", pos=s_start)
        mech = f"Exon–exon junction primer across exon {j.donor_order}–exon {j.acceptor_order}."

    if not s_eval.ok:
        flags.append("LOW_QC")

    # partner primer on the opposite side, forming a sensible amplicon
    partner = _choose_partner(seq, s_start, len(s_seq), specific_is_forward)
    fwd, rev = (specific, partner and partner[0]) if specific_is_forward else (partner and partner[0], specific)

    amplicon = None
    pair_dimer = None
    if fwd and rev and partner:
        p_start = partner[1]
        f_start = s_start if specific_is_forward else p_start
        r_start = p_start if specific_is_forward else s_start
        amplicon = (r_start + rev.length) - f_start
        pair_dimer = heterodimer_tm(fwd.seq, rev.seq)
        if pair_dimer >= STRUCT_TM_MAX:
            flags.append("PAIR_DIMER")

    if dtm < DELTA_TM_SAFE:
        flags.append("LOW_DELTA_TM")

    return PrimerDesign(
        tier=amp.tier, mechanism=mech, forward=fwd, reverse=rev,
        amplicon_len=amplicon, delta_tm=dtm,
        confidence="high" if (dtm >= DELTA_TM_SAFE and s_eval.ok) else "low",
        excluded_siblings=excluded, flags=sorted(set(flags)), pair_dimer_tm=pair_dimer,
    )


def _choose_partner(seq: str, spec_start: int, spec_len: int,
                    spec_is_forward: bool) -> tuple[Primer, int] | None:
    """Pick the partner primer (downstream if specific is forward, else upstream).
    Returns (Primer, start_position) or None. Start is the oligo's 5'-most template index."""
    if spec_is_forward:
        lo = max(spec_start + spec_len, spec_start + AMPLICON_MIN - LEN_MAX)
        cands = _partner(seq, lo, spec_start + AMPLICON_MAX)
        best, best_key = None, None
        for r, win in cands:
            amp_len = (r + len(win)) - spec_start
            if not (AMPLICON_MIN <= amp_len <= AMPLICON_MAX):
                continue
            ev = _evaluate(revcomp(win), TM_MIN)
            key = (ev.ok, ev.quality - abs(amp_len - AMPLICON_OPT) * 0.05)
            if best_key is None or key > best_key:
                best_key, best = key, (r, revcomp(win), ev)
        if not best:
            return None
        r, rev_seq, ev = best
        return _mk_primer(rev_seq, ev, "conventional", "reverse", "shared exon (downstream)", pos=r), r
    else:
        cands = _partner(seq, max(0, spec_start - AMPLICON_MAX), max(0, spec_start - AMPLICON_MIN))
        best, best_key = None, None
        for f, win in cands:
            amp_len = (spec_start + spec_len) - f
            if not (AMPLICON_MIN <= amp_len <= AMPLICON_MAX):
                continue
            ev = _evaluate(win, TM_MIN)
            key = (ev.ok, ev.quality - abs(amp_len - AMPLICON_OPT) * 0.05)
            if best_key is None or key > best_key:
                best_key, best = key, (f, win, ev)
        if not best:
            return None
        f, win, ev = best
        return _mk_primer(win, ev, "conventional", "forward", "shared exon (upstream)", pos=f), f


def _mk_primer(win_seq: str, ev: _Eval, kind: str, role: str, anchor: str,
               as_reverse: bool = False, pos: int = 0) -> Primer:
    oligo = revcomp(win_seq) if as_reverse else win_seq
    # Tm/structure computed on the actual oligo sequence. `pos` = 0-based start of the
    # binding site (the sense window) in the mRNA — same for forward and reverse.
    t = tm(oligo)
    return Primer(seq=oligo, tm=t, gc=round(gc_percent(oligo), 0), length=len(oligo),
                  kind=kind, role=role, anchor=anchor,
                  hairpin_tm=hairpin_tm(oligo), homodimer_tm=homodimer_tm(oligo),
                  qc_pass=ev.ok, tx_start=pos)


def _fallback(amp, seq, region, starts, k, sibs, excluded) -> PrimerDesign:
    """No QC-passing conventional candidate — return best-effort with LOW_QC."""
    ss = sorted(starts)
    s_start = ss[0] if ss else region.tx_start
    s_seq = seq[s_start:s_start + LEN_OPT]
    ev = _evaluate(s_seq, TM_MIN)
    fwd = _mk_primer(s_seq, ev, "conventional", "forward", f"exon {region.exon_order} (unique)", pos=s_start)
    partner = _choose_partner(seq, s_start, len(s_seq), True)
    dtm = round(tm(s_seq) - best_offtarget_tm(s_seq, sibs), 1)
    return PrimerDesign(
        tier=amp.tier,
        mechanism=f"Conventional primer in exon {region.exon_order} (QC-relaxed).",
        forward=fwd, reverse=partner[0] if partner else None,
        amplicon_len=((partner[1] + partner[0].length) - s_start) if partner else None,
        delta_tm=dtm, confidence="low", excluded_siblings=excluded, flags=["LOW_QC"],
    )
