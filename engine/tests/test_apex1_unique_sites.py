"""APEX1 — what "unique sequence" means for a near-identical exon 1.

NM_001641.4 and NM_001244249.2 differ by only 5 nt in exon 1, yet the tool used to
report a 23-nt unique span (mRNA 147-169) and "4 target-specific 20-nt primer
sites". Both numbers were true, but neither answered what a reader was asking:

  * Both exon 1s start at chr14:20,455,226. NM_001641.4 uses a 5-nt-later donor, so
    its exon 1 is 169 nt vs the sibling's 164 nt -- genomic 20,455,390-20,455,394,
    mRNA 165-169, is the sequence no sibling carries. THAT is the unique sequence,
    and it is what the API reports and the graph highlights.
  * mRNA 147-169 is the PLACEMENT ENVELOPE: a k-mer is target-specific as soon as it
    OVERLAPS the difference, so the positions a 20-nt primer may start at span up to
    (k-1) nt wider. Reporting it as "unique sequence" overstated a 5-nt difference
    more than fourfold. The envelope is still computed -- primers.py places primers
    from it -- but it stays internal (AmplifyResult), off the API.
  * The window COUNT (4) is a number of candidate primer positions, not a length,
    and is no longer offered as a measure of uniqueness anywhere in the UI.

One deliberate subtlety, pinned below so it is not "corrected" later: as raw mRNA
strings the two transcripts first differ at position 166, not 165, because the
sibling's base 165 (the first base of ITS exon 2) happens to be the same 'G'. That
is a coincidence of transcript concatenation. At the level the graph draws -- which
part of THIS EXON no sibling has -- the answer is the full 5 nt, 165-169.
"""

import random

from app.amplify import DEFAULT_K, analyze_amplifiability, uncovered_spans
from app.analyze import analyze

TARGET = "NM_001641.4"
NEAR_SIBLING = "NM_001244249.2"


def _by_acc(resp):
    return {t.accession: t for t in resp.transcripts}


def _exon1_region(resp):
    regions = [u for u in resp.target_verdict.unique_regions if u.exon_order == 1]
    assert len(regions) == 1
    return regions[0]


def test_exon1_geometry_is_a_five_nt_alternative_donor():
    """Same start coordinate, 5-nt-longer exon 1 on the target."""
    r = analyze(TARGET)
    v = _by_acc(r)
    t1 = v[TARGET].exons[0]
    s1 = v[NEAR_SIBLING].exons[0]
    assert (t1.tx_begin, t1.tx_end) == (1, 169)
    assert (s1.tx_begin, s1.tx_end) == (1, 164)
    assert t1.begin == s1.begin == 20455226      # shared 5' end
    assert t1.end - s1.end == 5                  # alternative donor, 5 nt later


def test_reported_unique_sequence_is_the_five_nt_difference():
    """The ticket's requirement: report the difference, 165-169 (5 nt)."""
    u = _exon1_region(analyze(TARGET))
    assert (u.tx_begin, u.tx_end) == (165, 169)
    assert u.uniq_len == 5
    assert (u.begin, u.end) == (20455390, 20455394)      # genomic, for the highlight
    # Not the old envelope, and not the placement count dressed up as a length.
    assert (u.tx_begin, u.tx_end) != (147, 169)
    assert u.uniq_len != u.window_count


def test_no_sibling_carries_any_of_the_reported_region():
    """Ground truth from the exon structures: the reported span is uncovered by ALL
    siblings, and the base just before it IS covered — so the boundary is tight."""
    r = analyze(TARGET)
    u = _exon1_region(r)
    sib_exons = [[(e.begin, e.end) for e in v.exons]
                 for v in r.transcripts if v.accession != TARGET]
    covered = {p for ex in sib_exons for (b, e) in ex for p in range(b, e + 1)}
    assert not (set(range(u.begin, u.end + 1)) & covered)   # nothing shared
    assert u.begin - 1 in covered                           # one base earlier is shared


def test_placement_envelope_still_exists_but_stays_internal():
    """primers.py needs the envelope to place oligos; it just is not the API's answer."""
    r = analyze(TARGET)
    amp = _amp_for_target(r)
    internal = [u for u in amp.unique_regions if u.exon_order == 1][0]
    assert (internal.tx_start + 1, internal.tx_end + 1) == (147, 169)   # the envelope
    assert internal.window_count == 4
    assert _by_acc(r)[TARGET].exons[0].unique_sites == 4
    # ...and it is genuinely wider than the difference it advertises.
    envelope = internal.tx_end - internal.tx_start + 1
    assert envelope == 23 > _exon1_region(r).uniq_len == 5


def test_the_first_differing_base_is_166_and_that_does_not_shrink_the_region():
    """The coincidence documented above: as strings they first differ at 166, because
    the sibling's base 165 is the same 'G'. The exon-level answer is still 5 nt."""
    r = analyze(TARGET)
    seq = r.target_mrna.upper()
    sib = analyze(NEAR_SIBLING).target_mrna.upper()
    first_diff = next(i for i in range(min(len(seq), len(sib))) if seq[i] != sib[i])
    assert first_diff == 165                      # 0-based -> mRNA 166
    assert seq[164] == sib[164] == "G"            # why it is 166 and not 165
    # The sibling's exon 1 nevertheless stops at 164, so its base 165 belongs to a
    # DIFFERENT exon -- the target's 165 is not shared exonic sequence.
    assert _by_acc(r)[NEAR_SIBLING].exons[0].tx_end == 164
    assert _exon1_region(r).tx_begin == 165


def test_every_reported_window_is_genuinely_absent_from_all_siblings():
    """The envelope's own ground truth, recomputed from raw sequence."""
    r = analyze(TARGET)
    seq = r.target_mrna.upper()
    sibs = [analyze(a).target_mrna.upper() for a in _by_acc(r) if a != TARGET]
    k, exon1_end = DEFAULT_K, 169

    unique_in_exon1 = [
        i for i in range(len(seq) - k + 1)
        if i + k <= exon1_end and all(seq[i:i + k] not in o for o in sibs)
    ]
    assert unique_in_exon1 == [146, 147, 148, 149]          # 0-based starts
    internal = [u for u in _amp_for_target(r).unique_regions if u.exon_order == 1][0]
    assert internal.tx_start == min(unique_in_exon1)
    assert internal.tx_end == max(unique_in_exon1) + k - 1
    # Every one of them overlaps the reported 5-nt difference — which is exactly why
    # a primer covering the difference is what specificity requires.
    u = _exon1_region(r)
    for s in unique_in_exon1:
        assert s + 1 <= u.tx_end and s + k >= u.tx_begin

    # The window just outside the envelope really is shared -- the boundary is tight.
    assert any(seq[145:145 + k] in o for o in sibs)


def test_near_sibling_exon1_has_no_unique_site():
    """Its exon 1 is a strict prefix of the target's, so nothing there is specific."""
    r = analyze(NEAR_SIBLING)
    assert _by_acc(r)[NEAR_SIBLING].exons[0].unique_sites == 0
    assert not [u for u in r.target_verdict.unique_regions if u.exon_order == 1]


def test_uncovered_spans_subtracts_every_sibling_not_just_the_nearest():
    """A closer sibling must be able to shrink the answer; the widest coverage wins."""
    target = [(100, 200)]
    # One sibling reaching 150, another reaching 180 -> only 181-200 is left.
    assert uncovered_spans(target, 1, [[(100, 150)], [(100, 180)]]) == [(181, 200)]
    # A sibling strictly inside leaves a gap on both sides.
    assert uncovered_spans(target, 1, [[(120, 180)]]) == [(100, 119), (181, 200)]
    # No overlapping sibling exon at all -> the whole exon is unique.
    assert uncovered_spans(target, 1, [[(900, 950)]]) == [(100, 200)]


def test_a_one_nt_difference_yields_a_2k_minus_1_envelope():
    """Synthetic control isolating the envelope rule from any APEX1 specifics.

    A 1-nt difference with room on BOTH sides is overlapped by k windows, covering
    2k-1 nt. APEX1's envelope is the narrower handle + k - 1 only because its
    difference is flush with the exon 3' end, so windows running past it are
    junction-spanning and leave exon 1's tally.
    """
    rng = random.Random(7)                          # non-repeating: no k-mer collisions
    head = "".join(rng.choice("ACGT") for _ in range(60))
    tail = "".join(rng.choice("ACGT") for _ in range(60))
    target, sibling = head + "A" + tail, head + "C" + tail
    kmers = [target[i:i + DEFAULT_K] for i in range(len(target) - DEFAULT_K + 1)]
    assert len(set(kmers)) == len(kmers), "fixture must not repeat a k-mer"

    exons = [(1, len(target))]                      # single exon: no junction windows
    u = analyze_amplifiability(exons, target, [sibling], k=DEFAULT_K).unique_regions[0]
    assert u.window_count == DEFAULT_K              # k windows overlap a 1-nt handle
    assert u.tx_end - u.tx_start + 1 == 2 * DEFAULT_K - 1 == 39
    # Without sibling exon structures there is nothing to subtract, so no region is
    # claimed rather than a wrong one being invented.
    assert u.uniq_span is None


def _amp_for_target(resp):
    """Re-run the sequence layer for the analyzed target (the API hides it)."""
    from app import ncbi
    seqs = {v.accession: ncbi.get_sequence(v.accession) for v in resp.transcripts}
    exons = {v.accession: [(e.begin, e.end) for e in v.exons] for v in resp.transcripts}
    acc = resp.target_accession
    sibs = [a for a in seqs if a != acc]
    return analyze_amplifiability(
        exons[acc], seqs[acc], [seqs[o] for o in sibs], k=DEFAULT_K,
        sibling_exons=[exons[o] for o in sibs])
