"""Sequence (ATGC) layer — the AUTHORITATIVE amplifiability test.

Slides a primer-length window along the target mRNA; a window absent from every
sibling mRNA is a target-specific primer site. Classifies each unique window as
exon-internal (conventional primer) or junction-spanning (EEJ primer), and assigns
the three-tier verdict.

Validated on GAPDH + MYC — see vault `01 Science/Sequence-Based Amplifiability (ATGC).md`.
"""

from __future__ import annotations

from dataclasses import dataclass, field

Interval = tuple[int, int]
DEFAULT_K = 20


@dataclass
class UniqueRegion:
    exon_order: int          # 1-based exon number (transcript 5'->3')
    tx_start: int            # 0-based inclusive start in the mRNA
    tx_end: int              # 0-based inclusive end
    window_count: int        # number of unique k-mers landing in this exon
    side: str                # "forward" | "reverse" | "either"
    # The DISCRIMINATING sequence itself: the largest genomic (begin,end) sub-span of this
    # exon that no sibling carries. This is the actual isoform difference, and it is what
    # the UI reports and highlights.
    #
    # It is much smaller than [tx_start, tx_end], which is the k-mer PLACEMENT ENVELOPE —
    # a window is target-specific as soon as it OVERLAPS the difference, so the envelope
    # runs up to (k-1) nt wider on each side. Reporting the envelope as "the unique
    # sequence" overstated a 5-nt alternative donor as a 23-nt unique stretch (APEX1
    # NM_001641.4 exon 1: envelope mRNA 147-169, real difference mRNA 165-169).
    # None only when sibling exon structures were not supplied.
    uniq_span: tuple[int, int] | None = None


@dataclass
class UniqueJunction:
    donor_order: int         # 1-based exon number
    acceptor_order: int
    window_count: int


@dataclass
class AmplifyResult:
    tier: str                # "CONVENTIONAL" | "NEEDS_EEJ" | "NO_SINGLE_UNIQUE_JUNCTION"
    amplifiable: bool
    needs_eej: bool
    unique_regions: list[UniqueRegion] = field(default_factory=list)
    unique_junctions: list[UniqueJunction] = field(default_factory=list)
    len_ok: bool = True      # sum(exon lengths) == len(mRNA)
    k: int = DEFAULT_K
    # positions (0-based mRNA start) of each target-unique k-mer window, so the primer
    # designer can guarantee specificity: any oligo that fully contains one of these
    # windows is itself absent from every sibling. See primers.py.
    internal_starts: dict[int, list[int]] = field(default_factory=dict)          # exon_order -> [starts]
    junction_starts: dict[tuple[int, int], list[int]] = field(default_factory=dict)  # (donor,acceptor) -> [starts]
    # for a 7c-Blue transcript (unique exon combination, no single unique window): the
    # exon pair (1-based forward, reverse) to target with a conventional primer pair.
    exon_pair: tuple[int, int] | None = None
    # combination rescues for would-be hard cases (no unique region/junction, no exon pair):
    # combo_je = (donor, acceptor, exon) — an EEJ across that junction + a conventional primer
    #   in that exon isolate the transcript (no sibling has both). NEEDS_EEJ.
    combo_je: tuple[int, int, int] | None = None
    # combo_jj = ((d1,a1),(d2,a2)) — two EEJs that together isolate it (no sibling has both).
    #   NEEDS_EEJ; the two junctions are the discriminating pair (shown magenta).
    combo_jj: tuple[tuple[int, int], tuple[int, int]] | None = None
    # For a junction+exon combo (combo_je): the GENOMIC (begin,end) sub-span of the combo exon
    # that actually distinguishes the transcript — the exon minus the parts the junction-holding
    # siblings still carry (the conventional primer must sit here, so only this part is the yellow
    # target site; the overlapped rest stays the tier color). None if the whole exon distinguishes.
    combo_exon_region: tuple[int, int] | None = None


def cumulative_exon_ends(exons: list[Interval]) -> list[int]:
    """Transcript-coordinate END (0-based exclusive) of each exon = running sum of lengths."""
    ends: list[int] = []
    s = 0
    for b, e in exons:
        s += (e - b + 1)
        ends.append(s)
    return ends


def _exon_at(cum: list[int], pos: int) -> int:
    """0-based transcript position -> 1-based exon number."""
    for i, c in enumerate(cum):
        if pos < c:
            return i + 1
    return len(cum)


def _exon_seqs(exons: list[Interval], seq: str) -> list[tuple[int, int, str]]:
    """(genomic_begin, genomic_end, exon_sequence) per exon, in transcript 5'->3' order."""
    cum = cumulative_exon_ends(exons)
    out: list[tuple[int, int, str]] = []
    for i, (b, e) in enumerate(exons):
        tx0 = cum[i - 1] if i else 0
        out.append((b, e, seq[tx0:cum[i]].upper()))
    return out


def _genomic_overlap(a: Interval, b: Interval) -> bool:
    lo1, hi1 = (a[0], a[1]) if a[0] <= a[1] else (a[1], a[0])
    lo2, hi2 = (b[0], b[1]) if b[0] <= b[1] else (b[1], b[0])
    return max(lo1, lo2) <= min(hi1, hi2)


def is_exon_subset(target_exons: list[Interval], target_seq: str,
                   sib_exons: list[Interval], sib_seq: str) -> bool:
    """True iff the target is a 'trimmed version' of the sibling — EVERY exon of the
    target is included in the sibling. An exon E1 is 'included' iff some sibling exon E2
    overlaps E1 at the same genomic locus AND E1's sequence is fully contained in E2's
    (per the agreed definition). If so, no single-exon (conventional) primer can tell the
    target from that sibling — only a junction can."""
    s_ex = _exon_seqs(sib_exons, sib_seq)
    for (tb, te, tseq) in _exon_seqs(target_exons, target_seq):
        if not any(_genomic_overlap((tb, te), (sb, se)) and tseq in sseq
                   for (sb, se, sseq) in s_ex):
            return False
    return True


def exon_subset_of_any(target_exons: list[Interval], target_seq: str,
                       siblings: list[tuple[list[Interval], str]]) -> bool:
    """True iff the target is an exon-subset of at least one sibling."""
    return any(is_exon_subset(target_exons, target_seq, se, ss) for (se, ss) in siblings)


def discriminating_exon_pair(
    target_exons: list[Interval], target_seq: str,
    siblings: list[tuple[list[Interval], str]],
) -> tuple[int, int] | None:
    """For a 7c transcript (unique exon combination but no single unique window), find the
    exon PAIR to put a conventional primer pair in: 1-based (forward_exon, reverse_exon),
    forward upstream of reverse, such that NO single sibling contains BOTH exons — so no
    sibling can produce that amplicon and the pair is transcript-specific. A conventional
    primer pair IS two primer sites (two exons); there is no 3+-exon single pair. Returns
    None if no 2-exon pair isolates the transcript — then it is not amplifiable by a single
    conventional primer pair. Prefers the most robust pair (fewest siblings sharing either
    exon), then the smallest exon span."""
    t_ex = _exon_seqs(target_exons, target_seq)
    sib_ex = [_exon_seqs(se, ss) for (se, ss) in siblings]

    def holders(exon: tuple[int, int, str]) -> set[int]:
        tb, te, tseq = exon
        return {si for si, s_ex in enumerate(sib_ex)
                if any(_genomic_overlap((tb, te), (sb, se)) and tseq in sseq
                       for (sb, se, sseq) in s_ex)}

    hold = [holders(e) for e in t_ex]
    best: tuple[tuple[int, int], int, int] | None = None
    n = len(t_ex)
    for i in range(n):
        for j in range(i + 1, n):
            if hold[i] & hold[j]:            # some sibling has both -> not specific
                continue
            score = (len(hold[i]) + len(hold[j]), j - i)   # robust first, then compact
            if best is None or score < best[0]:
                best = (score, i + 1, j + 1)
    return (best[1], best[2]) if best else None


def _exon_holders(target_exons: list[Interval], target_seq: str,
                  siblings: list[tuple[list[Interval], str]]) -> list[set[int]]:
    """holders[exon_index] = sibling indices whose sequence contains that exon (genomic
    overlap + sequence containment)."""
    t_ex = _exon_seqs(target_exons, target_seq)
    sib_ex = [_exon_seqs(se, ss) for (se, ss) in siblings]
    out = []
    for (tb, te, tseq) in t_ex:
        out.append({si for si, s_ex in enumerate(sib_ex)
                    if any(_genomic_overlap((tb, te), (sb, se)) and tseq in sseq
                           for (sb, se, sseq) in s_ex)})
    return out


def _junction_holders(target_exons: list[Interval], target_seq: str,
                      sibling_seqs: list[str], k: int) -> dict[tuple[int, int], set[int]]:
    """holders[(donor,acceptor)] = sibling indices whose mRNA contains the junction-spanning
    k-mer (centered on the boundary) — i.e. siblings that share that exact splice."""
    cum = cumulative_exon_ends(target_exons)
    others = [s.upper() for s in sibling_seqs]
    out: dict[tuple[int, int], set[int]] = {}
    for i in range(len(target_exons) - 1):
        b = cum[i]                                   # 0-based first base of exon i+2
        w = target_seq[b - k // 2:b - k // 2 + k]
        if len(w) < k:
            continue
        out[(i + 1, i + 2)] = {si for si, o in enumerate(others) if w in o}
    return out


def uncovered_spans(
    target_exons: list[Interval], exon_order: int, other_exons: list[list[Interval]],
) -> list[tuple[int, int]]:
    """Genomic sub-spans of the target's exon that NONE of `other_exons` carries.

    Isoforms of one gene are aligned to one genome, so "the sibling has this sequence" and
    "the sibling has an exon at this locus" are the same statement — which makes subtracting
    the siblings' exon coverage an exact way to find the sequence that distinguishes this
    exon, with no k-mer window length in the answer. Returns the gaps 5'->3' in GENOMIC
    order, or [(exon)] when nothing overlaps at all.
    """
    b, e = target_exons[exon_order - 1]
    covered = []
    for sib_exons in other_exons:
        for (sb, se) in sib_exons:
            lo, hi = max(sb, b), min(se, e)
            if lo <= hi:
                covered.append((lo, hi))
    if not covered:
        return [(b, e)]
    covered.sort()
    merged = [list(covered[0])]
    for lo, hi in covered[1:]:
        if lo <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], hi)
        else:
            merged.append([lo, hi])
    gaps, cur = [], b
    for lo, hi in merged:
        if lo > cur:
            gaps.append((cur, lo - 1))
        cur = max(cur, hi + 1)
    if cur <= e:
        gaps.append((cur, e))
    return gaps


def _largest_uncovered(
    target_exons: list[Interval], exon_order: int, other_exons: list[list[Interval]],
) -> tuple[int, int] | None:
    """The widest uncovered sub-span, or None when the exon is fully covered. A primer only
    has to OVERLAP unique sequence, so any gap works; the widest is the one worth naming."""
    gaps = uncovered_spans(target_exons, exon_order, other_exons)
    return max(gaps, key=lambda g: g[1] - g[0]) if gaps else None


def _combo_exon_region(
    target_exons: list[Interval], exon_order: int, hold_exons: list[list[Interval]],
) -> tuple[int, int] | None:
    """Distinguishing genomic sub-span of a combo exon: the exon MINUS the parts covered by the
    junction-holding siblings' exons at that locus (the siblings the conventional primer must
    still exclude — the others are already excluded by the junction)."""
    return _largest_uncovered(target_exons, exon_order, hold_exons)


def discriminating_junction_exon(
    target_exons: list[Interval], target_seq: str,
    siblings: list[tuple[list[Interval], str]], k: int = DEFAULT_K,
) -> tuple[int, int, int] | None:
    """Rescue a would-be hard case with a junction + exon combination: an EEJ across junction
    (donor,acceptor) plus a conventional primer in exon E, where NO single sibling has BOTH
    the junction and the exon — so the amplicon forms only in this transcript. Returns
    (donor, acceptor, exon_order) or None. Prefers the most discriminating, compact combo."""
    seq = target_seq.upper()
    e_hold = _exon_holders(target_exons, seq, siblings)
    j_hold = _junction_holders(target_exons, seq, [s for (_, s) in siblings], k)
    best: tuple[tuple[int, int], int, int, int] | None = None
    for (d, a), hj in sorted(j_hold.items(), key=lambda kv: len(kv[1])):
        for ei in range(len(e_hold)):
            eo = ei + 1
            if eo in (d, a) or (hj & e_hold[ei]):
                continue
            score = (len(hj) + len(e_hold[ei]), abs(eo - a))
            if best is None or score < best[0]:
                best = (score, d, a, eo)
    return (best[1], best[2], best[3]) if best else None


def discriminating_two_junctions(
    target_exons: list[Interval], target_seq: str,
    siblings: list[tuple[list[Interval], str]], k: int = DEFAULT_K,
) -> tuple[tuple[int, int], tuple[int, int]] | None:
    """Rescue a would-be hard case with two EEJs: junctions J1 (upstream) and J2 (downstream)
    where NO single sibling has BOTH — so an amplicon spanning J1..J2 forms only in this
    transcript. Returns ((d1,a1),(d2,a2)) or None. Most discriminating, compact first."""
    j_hold = _junction_holders(target_exons, target_seq.upper(),
                               [s for (_, s) in siblings], k)
    items = sorted(j_hold.items())
    best: tuple[tuple[int, int], tuple[int, int], tuple[int, int]] | None = None
    for x in range(len(items)):
        for y in range(x + 1, len(items)):
            (j1, h1), (j2, h2) = items[x], items[y]
            if h1 & h2:
                continue
            score = (len(h1) + len(h2), j2[0] - j1[1])
            if best is None or score < best[0]:
                best = (score, j1, j2)
    return (best[1], best[2]) if best else None


def analyze_amplifiability(
    target_exons: list[Interval],
    target_seq: str,
    sibling_seqs: list[str],
    k: int = DEFAULT_K,
    sibling_exons: list[list[Interval]] | None = None,
) -> AmplifyResult:
    """Core sequence test. `target_seq` and `sibling_seqs` are mRNA sequences (5'->3').

    A window is unique iff it is a substring of no sibling. Tiering (with `sibling_exons`
    supplied, paired index-wise with `sibling_seqs`):
      - **CONVENTIONAL** if it has a unique exon-internal window (original rule),
      - else **NEEDS_EEJ** if it has a junction-spanning unique window (a real unique
        junction always keeps EEJ — a junction primer is the right tool),
      - else **CONVENTIONAL** if the target is NOT an exon-subset of any sibling (rule 7c)
        **and** a discriminating 2-exon pair exists (no sibling carries both), so a
        conventional primer pair isolates it despite no single unique window,
      - else **NO_SINGLE_UNIQUE_JUNCTION** — including a not-a-subset transcript with no
        isolating 2-exon pair: it is not amplifiable by a single conventional primer pair
        (a 3+-exon combination is not one PCR), so it is a hard case, not Blue.
    Without `sibling_exons` (legacy) it falls back to the window-only tiering.
    """
    seq = target_seq.upper()
    cum = cumulative_exon_ends(target_exons)
    len_ok = (cum[-1] == len(seq)) if cum else False

    others = [s.upper() for s in sibling_seqs]
    internal: dict[int, list[int]] = {}   # exon# -> list of window start positions
    junction: dict[tuple[int, int], list[int]] = {}  # (donor,acceptor) -> window starts

    n = len(seq)
    for i in range(n - k + 1):
        w = seq[i:i + k]
        if all(w not in o for o in others):
            e1 = _exon_at(cum, i)
            e2 = _exon_at(cum, i + k - 1)
            if e1 == e2:
                internal.setdefault(e1, []).append(i)
            else:
                junction.setdefault((e1, e2), []).append(i)

    unique_regions: list[UniqueRegion] = []
    for exon_order, starts in sorted(internal.items()):
        first = min(starts)
        last = max(starts) + k - 1
        frac = (first / n) if n else 0.0
        side = "forward" if frac < 0.40 else ("reverse" if frac > 0.60 else "either")
        # The sequence that actually distinguishes this exon, independent of k. An exon can
        # only carry unique windows if some of it is uncovered, so this is normally present;
        # it is None only without sibling exon structures to subtract.
        uniq = (_largest_uncovered(target_exons, exon_order, sibling_exons)
                if sibling_exons is not None else None)
        unique_regions.append(
            UniqueRegion(exon_order=exon_order, tx_start=first, tx_end=last,
                         window_count=len(starts), side=side, uniq_span=uniq)
        )

    unique_junctions = [
        UniqueJunction(donor_order=d, acceptor_order=a, window_count=len(starts))
        for (d, a), starts in sorted(junction.items())
    ]

    not_subset = False
    if sibling_exons is not None:
        siblings = list(zip(sibling_exons, sibling_seqs))
        not_subset = not exon_subset_of_any(target_exons, seq, siblings)

    # Order matters: a genuine unique junction stays EEJ (a junction primer is the right
    # tool). The structural rule 7c only rescues what would otherwise be a hard case — a
    # transcript with no unique window at all that still isn't a trimmed copy of any
    # sibling, so a conventional primer *pair* spanning its unique exon combination works.
    exon_pair: tuple[int, int] | None = None
    combo_je: tuple[int, int, int] | None = None
    combo_jj: tuple[tuple[int, int], tuple[int, int]] | None = None
    have_sibs = sibling_exons is not None
    if unique_regions:
        tier = "CONVENTIONAL"
    elif unique_junctions:
        tier = "NEEDS_EEJ"
    elif not_subset and (exon_pair := discriminating_exon_pair(target_exons, seq, siblings)):
        # 7c: a 2-exon conventional pair isolates it (no sibling carries both). CONVENTIONAL.
        tier = "CONVENTIONAL"
    elif have_sibs and (combo_je := discriminating_junction_exon(target_exons, seq, siblings, k)):
        # rescue: an EEJ + a conventional exon primer isolate it. NEEDS_EEJ.
        tier = "NEEDS_EEJ"
    elif have_sibs and (combo_jj := discriminating_two_junctions(target_exons, seq, siblings, k)):
        # rescue: two EEJs together isolate it. NEEDS_EEJ (the two junctions shown magenta).
        tier = "NEEDS_EEJ"
    else:
        tier = "NO_SINGLE_UNIQUE_JUNCTION"

    # For a junction+exon combo, pin down which part of the combo exon actually distinguishes
    # the transcript (the part the junction-holding siblings lack) — only that is the target site.
    combo_exon_region: tuple[int, int] | None = None
    if combo_je:
        d, a, eo = combo_je
        jh = _junction_holders(target_exons, seq, [s for (_, s) in siblings], k).get((d, a), set())
        combo_exon_region = _combo_exon_region(target_exons, eo, [siblings[si][0] for si in jh])

    return AmplifyResult(
        tier=tier,
        amplifiable=tier in ("CONVENTIONAL", "NEEDS_EEJ"),
        needs_eej=(tier == "NEEDS_EEJ"),
        unique_regions=unique_regions,
        unique_junctions=unique_junctions,
        len_ok=len_ok,
        k=k,
        internal_starts=internal,
        junction_starts=junction,
        exon_pair=exon_pair,
        combo_je=combo_je,
        combo_jj=combo_jj,
        combo_exon_region=combo_exon_region,
    )
