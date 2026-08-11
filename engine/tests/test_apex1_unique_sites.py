"""APEX1 — the unique-site span reported for a near-identical exon 1.

Filed as a suspected recognition error: NM_001641.4 and NM_001244249.2 differ by
only 5 nt in exon 1, yet the tool reports a 23-nt unique span (mRNA 147-169).

It is not an error, and these tests pin why, so the numbers cannot drift or be
"simplified" back:

  * Both exon 1s start at chr14:20,455,226. NM_001641.4 uses a 5-nt-later donor,
    so its exon 1 is 169 nt vs the sibling's 164 nt.
  * The transcripts nevertheless first differ at mRNA 166, not 165 -- the target's
    first extra base is 'G' and the sibling's base 165 (the first base of ITS
    exon 2) is also 'G'. The real discriminating handle is 4 nt, not 5.
  * A k-mer is target-specific as soon as it OVERLAPS that handle, so the reported
    span is the PLACEMENT ENVELOPE for a k-nt primer -- up to (k-1) nt wider than
    the difference. 4 windows fit inside exon 1; they cover mRNA 147-169.

The span is therefore NOT a run of unique bases, which is what the original
tooltip wording implied and what prompted the ticket.
"""

import random

from app.amplify import DEFAULT_K, analyze_amplifiability
from app.analyze import analyze

TARGET = "NM_001641.4"
NEAR_SIBLING = "NM_001244249.2"


def _by_acc(resp):
    return {t.accession: t for t in resp.transcripts}


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


def test_reported_unique_span_is_147_169_with_four_windows():
    """The exact numbers from the ticket, and the count that explains them."""
    r = analyze(TARGET)
    regions = [u for u in r.target_verdict.unique_regions if u.exon_order == 1]
    assert len(regions) == 1
    u = regions[0]
    assert (u.tx_begin, u.tx_end) == (147, 169)
    # 23 nt of span, but only 4 placeable primer sites -- the span is an envelope.
    assert u.tx_end - u.tx_begin + 1 == 23
    assert _by_acc(r)[TARGET].exons[0].unique_sites == 4


def test_span_is_an_envelope_not_a_run_of_unique_bases():
    """The invariant behind the ticket: envelope width = handle + (k-1), not handle."""
    r = analyze(TARGET)
    u = [x for x in r.target_verdict.unique_regions if x.exon_order == 1][0]
    seq = r.target_mrna.upper()
    sib = analyze(NEAR_SIBLING).target_mrna.upper()

    first_diff = next(i for i in range(min(len(seq), len(sib))) if seq[i] != sib[i])
    assert first_diff == 165                      # 0-based -> mRNA 166
    assert seq[164] == sib[164] == "G"            # why it is 166 and not 165

    handle = 169 - (first_diff + 1) + 1           # mRNA 166..169 = 4 nt
    assert handle == 4
    # The whole point: the span is far wider than the sequence difference itself.
    # The handle is flush with the exon 3' end (it ends at 169), so windows running
    # past it become junction-spanning and drop out -> handle + k - 1, not 2k - 1.
    span = u.tx_end - u.tx_begin + 1
    assert span == handle + DEFAULT_K - 1 == 23


def test_every_reported_window_is_genuinely_absent_from_all_siblings():
    """Ground truth, recomputed from raw sequence rather than trusting the engine."""
    r = analyze(TARGET)
    seq = r.target_mrna.upper()
    sibs = [analyze(a).target_mrna.upper() for a in _by_acc(r) if a != TARGET]
    k, exon1_end = DEFAULT_K, 169

    unique_in_exon1 = [
        i for i in range(len(seq) - k + 1)
        if i + k <= exon1_end and all(seq[i:i + k] not in o for o in sibs)
    ]
    assert unique_in_exon1 == [146, 147, 148, 149]          # 0-based starts
    assert [i + 1 for i in unique_in_exon1] == [147, 148, 149, 150]
    u = [x for x in r.target_verdict.unique_regions if x.exon_order == 1][0]
    assert u.tx_begin == min(unique_in_exon1) + 1
    assert u.tx_end == max(unique_in_exon1) + k

    # The windows just outside the span really are shared -- the boundary is tight.
    assert any(seq[145:145 + k] in o for o in sibs)


def test_near_sibling_exon1_has_no_unique_site():
    """Its exon 1 is a strict prefix of the target's, so nothing there is specific."""
    r = analyze(NEAR_SIBLING)
    assert _by_acc(r)[NEAR_SIBLING].exons[0].unique_sites == 0
    assert not [u for u in r.target_verdict.unique_regions if u.exon_order == 1]


def test_a_one_nt_difference_yields_a_2k_minus_1_envelope():
    """Synthetic control isolating the rule from any APEX1 specifics.

    A 1-nt difference with room on BOTH sides is overlapped by k windows, covering
    2k-1 nt. APEX1 reports the narrower handle + k - 1 only because its difference
    is flush with the exon 3' end, so windows running past it are junction-spanning
    and leave exon 1's tally.
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
