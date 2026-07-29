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
        unique_regions.append(
            UniqueRegion(exon_order=exon_order, tx_start=first, tx_end=last,
                         window_count=len(starts), side=side)
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
    if unique_regions:
        tier = "CONVENTIONAL"
    elif unique_junctions:
        tier = "NEEDS_EEJ"
    elif not_subset:
        # 7c: it has a unique exon combination, but it is only CONVENTIONAL if a *2-exon*
        # pair actually isolates it (no sibling carries both). A conventional primer pair is
        # two primer sites — a 3+-exon combination is not a single PCR, so if no 2-exon pair
        # exists the transcript is NOT conventionally/singly amplifiable → hard case, NOT Blue.
        exon_pair = discriminating_exon_pair(target_exons, seq, siblings)
        tier = "CONVENTIONAL" if exon_pair else "NO_SINGLE_UNIQUE_JUNCTION"
    else:
        tier = "NO_SINGLE_UNIQUE_JUNCTION"

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
    )
