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


def analyze_amplifiability(
    target_exons: list[Interval],
    target_seq: str,
    sibling_seqs: list[str],
    k: int = DEFAULT_K,
) -> AmplifyResult:
    """Core sequence test. `target_seq` and `sibling_seqs` are mRNA sequences (5'->3').

    A window is unique iff it is a substring of no sibling. Exon-internal unique
    windows -> CONVENTIONAL; else junction-spanning unique windows -> NEEDS_EEJ;
    else NO_SINGLE_UNIQUE_JUNCTION.
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

    if unique_regions:
        tier = "CONVENTIONAL"
    elif unique_junctions:
        tier = "NEEDS_EEJ"
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
    )
