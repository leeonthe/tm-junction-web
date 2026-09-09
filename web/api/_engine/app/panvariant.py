"""One primer pair for the whole gene — the opposite question to the rest of the tool.

Everywhere else the job is to amplify ONE isoform and not its siblings. Here it is to
amplify them ALL: a total-expression assay measures a gene, not a variant, so the pair must
land where every transcript agrees and must give one band, not a ladder.

Three requirements, and the middle one is what makes this more than "primers in a shared
exon":

  * Common region — both primer sites must exist in every variant covered.
  * ONE product size — the amplicon must be the same length in every variant covered.
    Two transcripts can share both primer sites and still differ between them (a cassette
    exon spliced in by one and out by another), which gives two bands from one pair and
    makes the assay unquantifiable.
  * At least two exons — the product must cross a junction, or it cannot be told from one
    amplified off contaminating genomic DNA (the same rule as every other design here).

Structure proposes, sequence decides. Candidates are placed inside a RUN of consecutive
exons — identical exons in the same order with nothing spliced in between means the
distance between any two points in that run is the same in every variant carrying it, which
is where a same-size product is most likely to be found. But the run is a search heuristic,
not the verdict: a pair is credited only with the variants it is SHOWN to amplify, by
locating both sites in each transcript and measuring the product. That is what lets it
cover variants whose exon boundaries differ outside the primer sites — MYC's two
transcripts share no identical run at all, yet one pair measures both at 138 bp.

When no run reaches every variant, the best run by COVERAGE wins — the most transcripts one
pair can measure — and the variants left out are named rather than quietly dropped.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .amplify import Interval, cumulative_exon_ends
from .primers import (
    AMPLICON_MAX, AMPLICON_MIN, AMPLICON_OPT, LEN_MIN, Primer, _mk_primer, _p3_ranked,
    _sweep_pair,
    heterodimer_tm, revcomp, STRUCT_TM_MAX,
)

# How many (run, junction) placements to hand to primer3 before settling for the best found.
# Runs are tried best-coverage first, so this bounds work without bounding coverage.
MAX_PLACEMENTS = 40


@dataclass
class PanVariant:
    """A pair that amplifies `covered` with one product length, or None if none exists."""
    forward: Primer
    reverse: Primer
    amplicon_len: int                       # identical in every covered variant
    covered: list[str] = field(default_factory=list)
    uncovered: list[str] = field(default_factory=list)
    reference: str = ""                     # whose exon numbering `exons` refers to
    exons: list[int] = field(default_factory=list)   # 1-based, the two the primers sit in
    flags: list[str] = field(default_factory=list)


def _consecutive_index(hay: list[Interval], needle: tuple[Interval, ...]) -> int | None:
    """Where `needle` sits in `hay` as a CONSECUTIVE run, or None.

    Consecutive matters more than present: a variant carrying the same exons with another
    one spliced in between makes a longer product from the same primer sites, which is the
    second band this design exists to avoid.
    """
    k = len(needle)
    for s in range(len(hay) - k + 1):
        if tuple(hay[s:s + k]) == needle:
            return s
    return None


def _amplifies(seq: str, fwd: str, rev: str) -> tuple[int, int] | None:
    """Where this pair amplifies `seq`, as (product_start, product_end), or None.

    Both sites must occur exactly ONCE: a primer with a second binding site in the same
    transcript makes the product ambiguous, which is the same problem as a second band.
    """
    f = seq.find(fwd)
    if f < 0 or seq.find(fwd, f + 1) >= 0:
        return None
    site = revcomp(rev)
    r = seq.find(site)
    if r < 0 or seq.find(site, r + 1) >= 0:
        return None
    if r + len(site) <= f + len(fwd):        # the reverse primer must sit downstream
        return None
    return f, r + len(site)


def _coverage(pair: tuple[str, str], seqs: dict[str, str],
              accs: list[str]) -> tuple[list[str], list[str], int | None]:
    """Which variants this pair amplifies AT ONE SIZE — the sequence-level verdict.

    A variant that amplifies at a different length is NOT covered: it would show up as a
    second band. The size that carries the most variants wins, so one outlier cannot veto
    a pair that measures everything else cleanly.
    """
    fwd, rev = pair
    hits: dict[str, int] = {}
    for a in accs:
        got = _amplifies(seqs[a], fwd, rev)
        if got:
            hits[a] = got[1] - got[0]
    if not hits:
        return [], list(accs), None
    sizes: dict[int, list[str]] = {}
    for a, n in hits.items():
        sizes.setdefault(n, []).append(a)
    size = max(sizes, key=lambda n: (len(sizes[n]), -n))
    covered = [a for a in accs if a in sizes[size]]
    return covered, [a for a in accs if a not in sizes[size]], size


# How many options the whole-transcript designer offers, and how many distinct pairs a
# single junction placement may contribute — several placements is the interesting axis
# (different junctions, different products), so one junction may not fill the list alone.
MAX_OPTIONS = 5
PER_PLACEMENT = 3


def design(transcripts: dict[str, dict], seqs: dict[str, str],
           order: list[str] | None = None) -> PanVariant | None:
    """The single best pair — see design_options; kept for callers wanting one answer."""
    got = design_options(transcripts, seqs, order)
    return got[0] if got else None


def design_options(transcripts: dict[str, dict], seqs: dict[str, str],
                   order: list[str] | None = None) -> list[PanVariant]:
    """Ranked pair options for the gene: most variants covered, one product size, two exons.

    A list, not a winner, for the same reason the EEJ-independent designer offers one: the
    best pair by this ranking is not always the best pair for someone's assay (a probe to
    fit, a size to match an old gel, a primer already in the freezer). Best first — most
    variants covered, then clean over flagged, then the tidier product.

    `transcripts` maps accession -> {"exons": [...]} in transcript order; `order` fixes the
    reporting order (defaults to the mapping's).
    """
    accs = list(order or transcripts)
    if not accs:
        return []
    exons_by_acc = {a: [tuple(e) for e in transcripts[a]["exons"]] for a in accs}

    # The search space is which JUNCTION to straddle, not which region to sit in: the pair
    # must cross one, and only the AMPLICON_MAX bases either side of it can hold a primer
    # that makes a legal product. Each junction is identified by its two flanking exons, so
    # the same junction proposed by several variants is one placement, not several.
    placements: dict[tuple, tuple[list[str], str, int]] = {}
    for ref in accs:
        ex = exons_by_acc[ref]
        for jx in range(len(ex) - 1):
            key = (ex[jx], ex[jx + 1])
            if key in placements:
                continue
            # How many variants splice these two exons straight together. Where they do, the
            # distance across the junction is identical, so a product spanning it is one
            # size in all of them — the reason to try the best-carried junctions first.
            carriers = [a for a in accs
                        if _consecutive_index(exons_by_acc[a], key) is not None]
            placements[key] = (carriers, ref, jx)
    if not placements:
        return []

    ranked = sorted(placements.items(),
                    key=lambda kv: (-len(kv[1][0]), kv[1][1], kv[1][2]))

    options: list[PanVariant] = []
    tried = 0
    for _key, (carriers, ref, jx) in ranked:
        # No early exit on carrier count: a junction carried by one variant can still yield a
        # pair that amplifies several, when what differs between them lies outside the primer
        # sites. Only the verified coverage below is allowed to rank a candidate.
        if tried >= MAX_PLACEMENTS:
            break
        seq = seqs[ref]
        cum = cumulative_exon_ends(exons_by_acc[ref])
        boundary = cum[jx]                                  # 0-based first base after it
        ex_lo = cum[jx - 1] if jx else 0                    # start of the donor exon
        ex_hi = cum[jx + 1] - 1                             # last base of the acceptor exon
        lo = max(ex_lo, boundary - AMPLICON_MAX)
        hi = min(ex_hi, boundary + AMPLICON_MAX - 1)
        if boundary - lo < LEN_MIN or hi - boundary + 1 < LEN_MIN:
            continue
        tried += 1
        # primer3 gets the window, not the transcript: the sites are confined to it either
        # way, and handing over a 5 kb template to place a 140 bp product inside 400 of them
        # costs seconds per call. Coordinates come back window-relative.
        sub = seq[lo:hi + 1]
        picks = _p3_ranked(sub, (0, boundary - 1 - lo), (boundary - lo, hi - lo))[:PER_PLACEMENT]
        if not picks:
            one = _sweep_pair(sub, (0, boundary - 1 - lo), (boundary - lo, hi - lo))
            picks = [one] if one else []
        for pick in picks:
            (f_rel, f_seq, f_ev), (r_rel, r_seq, r_ev) = pick
            f_start, r_start = f_rel + lo, r_rel + lo
            covered, uncovered, size = _coverage((f_seq, r_seq), seqs, accs)
            if size is None or not covered or not (AMPLICON_MIN <= size <= AMPLICON_MAX):
                continue
            flags: list[str] = []
            if not (f_ev.ok and r_ev.ok):
                flags.append("LOW_QC")
            if heterodimer_tm(f_seq, r_seq) >= STRUCT_TM_MAX:
                flags.append("PAIR_DIMER")
            options.append(PanVariant(
                forward=_mk_primer(f_seq, f_ev, "conventional", "forward",
                                   f"exon {_exon_order(cum, f_start)} (all-variant)", pos=f_start),
                reverse=_mk_primer(r_seq, r_ev, "conventional", "reverse",
                                   f"exon {_exon_order(cum, r_start)} (all-variant)", pos=r_start),
                amplicon_len=size, covered=covered, uncovered=uncovered, reference=ref,
                exons=[_exon_order(cum, f_start), _exon_order(cum, r_start)],
                flags=sorted(set(flags)),
            ))
        full_clean = [o for o in options
                      if len(o.covered) == len(accs) and not o.flags]
        if len(full_clean) >= MAX_OPTIONS:
            break                                           # enough of the best class
    options.sort(key=lambda o: (-len(o.covered), bool(o.flags),
                                abs(o.amplicon_len - AMPLICON_OPT)))
    # Dedup by oligo pair — the same primers can fall out of two overlapping placements.
    seen: set[tuple[str, str]] = set()
    out: list[PanVariant] = []
    for o in options:
        key = (o.forward.seq, o.reverse.seq)
        if key in seen:
            continue
        seen.add(key)
        out.append(o)
        if len(out) >= MAX_OPTIONS:
            break
    return out


def _exon_order(cum: list[int], pos: int) -> int:
    for i, c in enumerate(cum):
        if pos < c:
            return i + 1
    return len(cum)



