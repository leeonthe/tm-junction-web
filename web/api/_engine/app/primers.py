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
    # single unique exon window — no single primer is specific. The specificity is the
    # exon COMBINATION: no sibling carries both exons of amp.exon_pair, so a primer pair
    # spanning them makes an amplicon only this transcript can form. Design that pair.
    if amp.tier == "CONVENTIONAL" and not amp.unique_regions:
        if amp.exon_pair:
            return _design_exon_pair(amp, seq, cum, excluded)
        return PrimerDesign(
            tier=amp.tier,
            mechanism="Structurally unique (not a trimmed copy of any sibling), but no single "
                      "unique exon window and no isolating 2-exon pair.",
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
        # Only candidates a junction-spanning partner can pair with; see _can_span.
        spannable = [c for c in spec_cands if _can_span(cum, c[0], len(c[1]))]
        chosen = _best_candidate(spannable or spec_cands, TM_MIN, sibs, delta_weight=0.3)
        if chosen is None:
            return _fallback(amp, seq, region, starts, k, sibs, excluded, cum)
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
    partner = _choose_partner(seq, s_start, len(s_seq), specific_is_forward, cum)
    if partner is None:
        # The preferred side cannot reach a junction — a specific primer in the terminal
        # exon has nothing downstream of it. Try the other orientation before giving up:
        # the specific oligo is the same sequence either way, only its role changes.
        flipped = _choose_partner(seq, s_start, len(s_seq), not specific_is_forward, cum)
        if flipped is not None:
            specific_is_forward = not specific_is_forward
            specific = _mk_primer(s_seq, s_eval, specific.kind,
                                  "forward" if specific_is_forward else "reverse",
                                  specific.anchor, as_reverse=not specific_is_forward,
                                  pos=s_start)
            partner = flipped
    fwd, rev = (specific, partner and partner[0]) if specific_is_forward else (partner and partner[0], specific)

    # Terminal check of the intron-spanning rule. _choose_partner already refuses a partner
    # that keeps the product inside one exon, so this fires only if some future path reaches
    # here another way — the rule is a correctness requirement (a single-exon product cannot
    # be told from one amplified off contaminating gDNA), so it is enforced where the pair is
    # emitted rather than trusted to every branch that builds one.
    if fwd and rev and partner and not _spans_junction(
            cum, min(s_start, partner[1]),
            max(s_start + len(s_seq), partner[1] + partner[0].length)):
        # Keep the specific oligo — it is still the right primer — and drop the partner that
        # would have made an unusable product, rather than emit the pair.
        partner = None
        fwd, rev = (specific, None) if specific_is_forward else (None, specific)
        flags.append("NO_SPANNING_PAIR")

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


def _design_exon_pair(amp: AmplifyResult, seq: str, cum: list[int],
                      excluded: list[str]) -> PrimerDesign:
    """7c-Blue design: forward in exon FA, reverse in exon FB (amp.exon_pair).

    The PAIR is the specificity — no single sibling carries both exons, so the amplicon
    forms only on the target. The individual oligos need not be unique, which is why the
    unique-window sweep above does not apply. Instead primer3's own pair designer (the
    engine behind Primer3Plus) picks a Tm-matched, structure-checked pair, confined to
    the two exons via SEQUENCE_PRIMER_PAIR_OK_REGION_LIST; a plain QC sweep is the
    fallback if primer3 is unavailable or finds nothing under its constraints.
    """
    fa, fb = amp.exon_pair                                 # 1-based, fa < fb
    a = ((cum[fa - 2] if fa > 1 else 0), cum[fa - 1] - 1)  # 0-based inclusive tx spans
    b = ((cum[fb - 2] if fb > 1 else 0), cum[fb - 1] - 1)
    mech = (f"Conventional primer pair spanning the unique exon combination — forward in "
            f"exon {fa}, reverse in exon {fb}. No other isoform carries both exons, so "
            f"the amplicon forms only on this transcript.")
    pick = _p3_pair(seq, a, b) if HAS_PRIMER3 else None
    if pick is None:
        pick = _sweep_pair(seq, a, b)
    if pick is None:                                       # exons too short for any oligo
        return PrimerDesign(tier=amp.tier, mechanism=mech,
                            excluded_siblings=excluded, flags=["NO_PAIR"])
    (f_start, f_seq, f_ev), (r_start, r_seq, r_ev) = pick
    # Two different exons put a junction between them by construction, but the pair is
    # emitted here, so the rule is checked here too rather than assumed.
    if not _spans_junction(cum, f_start, r_start + len(r_seq)):
        return PrimerDesign(tier=amp.tier, mechanism=mech,
                            excluded_siblings=excluded, flags=["NO_SPANNING_PAIR"])
    fwd = _mk_primer(f_seq, f_ev, "conventional", "forward", f"exon {fa} (pair)", pos=f_start)
    rev = _mk_primer(r_seq, r_ev, "conventional", "reverse", f"exon {fb} (pair)", pos=r_start)
    pair_dimer = heterodimer_tm(fwd.seq, rev.seq)
    flags = []
    if not (f_ev.ok and r_ev.ok):
        flags.append("LOW_QC")
    if pair_dimer >= STRUCT_TM_MAX:
        flags.append("PAIR_DIMER")
    return PrimerDesign(
        tier=amp.tier, mechanism=mech, forward=fwd, reverse=rev,
        amplicon_len=(r_start + len(r_seq)) - f_start,
        delta_tm=None,   # specificity is combinatorial here, not a per-oligo Tm margin
        confidence="high" if (f_ev.ok and r_ev.ok) else "low",
        excluded_siblings=excluded, flags=sorted(set(flags)), pair_dimer_tm=pair_dimer,
    )


def _p3_ranked(seq: str, a: tuple[int, int], b: tuple[int, int]):
    """primer3 pair design confined to sense-strand spans a (left) and b (right).

    Returns a ranked list of ((f_start, f_oligo, f_eval), (r_start, r_oligo, r_eval)) with
    starts = 0-based binding-site starts on the mRNA, oligos 5'->3' as ordered. Pairs that
    pass OUR QC gate come first (each class in primer3's own penalty order); [] on failure.
    """
    (a_lo, a_hi), (b_lo, b_hi) = a, b
    lo = max(AMPLICON_MIN, b_lo - a_hi + 2 * LEN_MIN - 1)  # shortest reachable product
    hi = b_hi - a_lo + 1
    if lo > hi:
        return []
    try:
        res = primer3.bindings.design_primers(
            {
                "SEQUENCE_ID": "target",
                "SEQUENCE_TEMPLATE": seq,
                "SEQUENCE_PRIMER_PAIR_OK_REGION_LIST":
                    [[a_lo, a_hi - a_lo + 1, b_lo, b_hi - b_lo + 1]],
            },
            {
                "PRIMER_TASK": "generic",
                "PRIMER_PICK_LEFT_PRIMER": 1, "PRIMER_PICK_RIGHT_PRIMER": 1,
                "PRIMER_PICK_INTERNAL_OLIGO": 0,
                "PRIMER_NUM_RETURN": 10,
                "PRIMER_MIN_SIZE": LEN_MIN, "PRIMER_OPT_SIZE": LEN_OPT,
                "PRIMER_MAX_SIZE": LEN_MAX,
                "PRIMER_MIN_TM": TM_MIN, "PRIMER_OPT_TM": TM_OPT, "PRIMER_MAX_TM": TM_MAX,
                "PRIMER_MIN_GC": GC_MIN, "PRIMER_MAX_GC": GC_MAX,
                "PRIMER_GC_CLAMP": 1,            # match _evaluate's 3'-clamp gate
                "PRIMER_MAX_POLY_X": 4,          # match _homopolymer (5+ run rejected)
                "PRIMER_PRODUCT_SIZE_RANGE": [[lo, hi]],
                # nudge toward compact products; quality still dominates the penalty
                "PRIMER_PRODUCT_OPT_SIZE": min(max(AMPLICON_OPT, lo), hi),
                "PRIMER_PAIR_WT_PRODUCT_SIZE_GT": 0.05,
                "PRIMER_PAIR_WT_PRODUCT_SIZE_LT": 0.05,
            })
    except Exception:
        return []
    # QC-passing pairs first (primer3's own penalty order within each class), then the
    # best-effort ones — a caller wanting one pair takes [0]; the whole-transcript
    # designer takes several, so one placement can feed a real options list.
    passing, fallback, seen = [], [], set()
    for i in range(res.get("PRIMER_PAIR_NUM_RETURNED", 0)):
        f_seq = res[f"PRIMER_LEFT_{i}_SEQUENCE"].upper()
        r_seq = res[f"PRIMER_RIGHT_{i}_SEQUENCE"].upper()
        if (f_seq, r_seq) in seen:
            continue
        seen.add((f_seq, r_seq))
        f_start = res[f"PRIMER_LEFT_{i}"][0]
        r_pos, r_len = res[f"PRIMER_RIGHT_{i}"]   # r_pos = 3'-most template index
        f_ev, r_ev = _evaluate(f_seq, TM_MIN), _evaluate(r_seq, TM_MIN)
        cand = ((f_start, f_seq, f_ev), (r_pos - r_len + 1, r_seq, r_ev))
        (passing if f_ev.ok and r_ev.ok else fallback).append(cand)
    return passing + fallback


def _region_best(seq: str, lo: int, hi: int, as_reverse: bool):
    """Best QC oligo whose binding site fits inside [lo, hi] (0-based inclusive)."""
    best, best_key = None, None
    for s in range(max(0, lo), hi - LEN_MIN + 2):
        for L in range(LEN_MIN, LEN_MAX + 1):
            if s + L - 1 > hi or s + L > len(seq):
                break
            oligo = revcomp(seq[s:s + L]) if as_reverse else seq[s:s + L]
            ev = _evaluate(oligo, TM_MIN)
            key = (ev.ok, ev.quality)
            if best_key is None or key > best_key:
                best_key, best = key, (s, oligo, ev)
    return best


def _p3_pair(seq: str, a: tuple[int, int], b: tuple[int, int]):
    """The single best pair from _p3_ranked, or None — the one-pair callers' entry."""
    ranked = _p3_ranked(seq, a, b)
    return ranked[0] if ranked else None


def _sweep_pair(seq: str, a: tuple[int, int], b: tuple[int, int]):
    """Fallback pair pick: best-QC oligo per exon, no pair-level optimization."""
    f = _region_best(seq, a[0], a[1], as_reverse=False)
    r = _region_best(seq, b[0], b[1], as_reverse=True)
    return (f, r) if f and r else None


def _exon_of(cum: list[int], pos: int) -> int:
    """0-based mRNA position -> index of the exon holding it (cum = exclusive exon ends)."""
    for i, c in enumerate(cum):
        if pos < c:
            return i
    return len(cum) - 1


def _spans_junction(cum: list[int], start: int, end: int) -> bool:
    """Does the amplicon [start, end) cross at least one exon-exon junction?

    A product contained in ONE exon is indistinguishable from one amplified off
    contaminating genomic DNA -- the same primer sites sit uninterrupted in the genome.
    Spanning a junction makes gDNA either fail or give a visibly longer band, so this is a
    correctness requirement of the design rather than a preference. Unknown exon structure
    passes: better to return the pair than to reject everything on a claim we cannot check.
    """
    if not cum:
        return True
    return _exon_of(cum, start) != _exon_of(cum, end - 1)


def _can_span(cum: list[int], s: int, length: int) -> bool:
    """Could SOME partner put a junction inside this primer's amplicon?

    Pure geometry, checked BEFORE the specific primer is chosen. Without it the search picks
    the best-QC specific primer first and only then discovers no partner can reach a
    junction -- which is how ACTB, whose unique region is a 744-nt terminal exon, ended up
    with a forward primer and no reverse at all. A downstream partner needs a boundary
    within AMPLICON_MAX of the primer's start; an upstream one, within AMPLICON_MAX of its
    end. Either orientation counts here; _choose_partner decides which is actually used.
    """
    if not cum:
        return True
    e = s + length
    return any(s < b < s + AMPLICON_MAX or e - AMPLICON_MAX < b < e
               for b in cum[:-1])          # cum[-1] is the transcript end, not a junction


def _choose_partner(seq: str, spec_start: int, spec_len: int,
                    spec_is_forward: bool, cum: list[int] | None = None) -> tuple[Primer, int] | None:
    """Pick the partner primer (downstream if specific is forward, else upstream).
    Returns (Primer, start_position) or None. Start is the oligo's 5'-most template index.

    The pair must span a junction (see _spans_junction), so a partner that would keep the
    whole product inside one exon is rejected however good its Tm."""
    cum = cum or []
    if spec_is_forward:
        lo = max(spec_start + spec_len, spec_start + AMPLICON_MIN - LEN_MAX)
        cands = _partner(seq, lo, spec_start + AMPLICON_MAX)
        best, best_key = None, None
        for r, win in cands:
            amp_len = (r + len(win)) - spec_start
            if not (AMPLICON_MIN <= amp_len <= AMPLICON_MAX):
                continue
            if not _spans_junction(cum, spec_start, r + len(win)):
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
            if not _spans_junction(cum, f, spec_start + spec_len):
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


def _fallback(amp, seq, region, starts, k, sibs, excluded, cum) -> PrimerDesign:
    """No QC-passing conventional candidate — return best-effort with LOW_QC."""
    ss = sorted(starts)
    s_start = ss[0] if ss else region.tx_start
    s_seq = seq[s_start:s_start + LEN_OPT]
    ev = _evaluate(s_seq, TM_MIN)
    fwd = _mk_primer(s_seq, ev, "conventional", "forward", f"exon {region.exon_order} (unique)", pos=s_start)
    partner = _choose_partner(seq, s_start, len(s_seq), True, cum)
    dtm = round(tm(s_seq) - best_offtarget_tm(s_seq, sibs), 1)
    return PrimerDesign(
        tier=amp.tier,
        mechanism=f"Conventional primer in exon {region.exon_order} (QC-relaxed).",
        forward=fwd, reverse=partner[0] if partner else None,
        amplicon_len=((partner[1] + partner[0].length) - s_start) if partner else None,
        delta_tm=dtm, confidence="low", excluded_siblings=excluded, flags=["LOW_QC"],
    )
