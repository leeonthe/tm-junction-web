"""Coordinate (structural) layer — the exon-overlap structural-containment rule and
the coordinate junction finder.

Ported from the sibling `ncbi exon` project (`ncbi_exon/parser.py`). This layer is
the FAST STRUCTURAL PROXY used for the exon-track graph and as a cross-check. The
authoritative amplifiability answer comes from the sequence layer (`amplify.py`).

See vault: `01 Science/Structurally Contained Transcript Definition.md`,
`01 Science/Amplifiability and Primer Logic.md`.
"""

from __future__ import annotations

Interval = tuple[int, int]


def intervals_overlap(a: Interval, b: Interval) -> bool:
    """Definition #1 — 1-based inclusive interval overlap: share >= 1 base."""
    return a[0] <= b[1] and b[0] <= a[1]


def is_structurally_contained_in(t_exons: list[Interval], other_exons: list[Interval]) -> bool:
    """t is structurally contained in `other` iff EVERY exon of t overlaps >= 1 exon of other
    (the term was "non-unique" until 2026-10-23)."""
    if not t_exons or not other_exons:
        return False
    return all(any(intervals_overlap(ea, eb) for eb in other_exons) for ea in t_exons)


def containing_partners(
    target_exons: list[Interval],
    siblings: dict[str, list[Interval]],
) -> list[str]:
    """Return the sibling accessions that structurally contain `target` (pairwise, one-vs-one).

    `siblings` excludes the target. Never takes the union of all siblings.
    """
    return [
        acc
        for acc, ex in siblings.items()
        if is_structurally_contained_in(target_exons, ex)
    ]


def coord_unique_junctions(
    target_exons: list[Interval],
    siblings: dict[str, list[Interval]],
) -> list[tuple[int, int]]:
    """Coordinate junction finder: return (donor_exon#, acceptor_exon#) 1-based pairs
    whose (donor.end, acceptor.begin) coordinate pair is reproduced by no sibling.

    Cross-checks the sequence junction finder — validated to agree on GAPDH/MYC.
    """
    def junctions(ex: list[Interval]) -> set[tuple[int, int]]:
        return {(ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1)}

    sibling_junctions: set[tuple[int, int]] = set()
    for ex in siblings.values():
        sibling_junctions |= junctions(ex)

    out: list[tuple[int, int]] = []
    for i in range(len(target_exons) - 1):
        pair = (target_exons[i][1], target_exons[i + 1][0])
        if pair not in sibling_junctions:
            out.append((i + 1, i + 2))
    return out
