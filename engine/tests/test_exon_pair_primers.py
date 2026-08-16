"""7c-Blue (unique exon combination) — the primer PAIR is now designed, not deferred.

A transcript can have no unique window at all and still be CONVENTIONAL: rule 7c
rescues it when some 2-exon combination exists that no single sibling carries. The
specificity then lives in the PAIR — forward in one exon, reverse in the other —
because the amplicon needs both primer sites on one molecule and only the target
has them. primers.design used to return "pair design is a follow-up" with no
primers for this tier (ticket 12); these tests pin the design it now produces,
via primer3's own pair engine (the one behind Primer3Plus) confined to the two
exons with SEQUENCE_PRIMER_PAIR_OK_REGION_LIST.

Scenario: target = exons A+B+C. Sibling S1 = A+B, sibling S2 = B+C (same genomic
loci, identical sequences). Every internal window and both splice junctions of
the target exist in a sibling, so no single unique window — but no sibling has
BOTH A and C, so the pair (exon 1, exon 3) isolates the target.
"""

import random

from app.amplify import analyze_amplifiability
from app import primers

_R = random.Random(7)
_A = "".join(_R.choice("ACGT") for _ in range(150))
_B = "".join(_R.choice("ACGT") for _ in range(120))
_C = "".join(_R.choice("ACGT") for _ in range(180))

TARGET_EXONS = [(1000, 1149), (2000, 2119), (3000, 3179)]
TARGET_SEQ = _A + _B + _C
S1_EXONS = [(1000, 1149), (2000, 2119)]          # A + B
S1_SEQ = _A + _B
S2_EXONS = [(2000, 2119), (3000, 3179)]          # B + C
S2_SEQ = _B + _C

A_SPAN = (0, 149)                                 # 0-based inclusive tx spans
C_SPAN = (270, 449)


def _amp():
    return analyze_amplifiability(
        TARGET_EXONS, TARGET_SEQ, [S1_SEQ, S2_SEQ], k=20,
        sibling_exons=[S1_EXONS, S2_EXONS])


def _design():
    return primers.design(TARGET_EXONS, TARGET_SEQ,
                          {"S1": S1_SEQ, "S2": S2_SEQ}, _amp())


def test_scenario_is_7c_blue_with_the_expected_exon_pair():
    amp = _amp()
    assert amp.tier == "CONVENTIONAL"
    assert not amp.unique_regions and not amp.unique_junctions
    assert amp.exon_pair == (1, 3)


def test_pair_is_designed_not_deferred():
    d = _design()
    assert d.forward is not None and d.reverse is not None
    assert "NO_UNIQUE_WINDOW" not in d.flags
    assert d.forward.role == "forward" and d.reverse.role == "reverse"
    assert d.forward.kind == d.reverse.kind == "conventional"
    assert "exon 1" in d.forward.anchor and "exon 3" in d.reverse.anchor
    assert "pair" in d.mechanism


def test_primers_sit_inside_their_exons_and_amplicon_spans_them():
    d = _design()
    f, r = d.forward, d.reverse
    assert A_SPAN[0] <= f.tx_start and f.tx_start + f.length - 1 <= A_SPAN[1]
    assert C_SPAN[0] <= r.tx_start and r.tx_start + r.length - 1 <= C_SPAN[1]
    assert d.amplicon_len == (r.tx_start + r.length) - f.tx_start
    # the ordered oligos really bind the target: forward is a sense window, the
    # reverse oligo is the reverse complement of one
    assert TARGET_SEQ[f.tx_start:f.tx_start + f.length] == f.seq
    assert TARGET_SEQ[r.tx_start:r.tx_start + r.length] == primers.revcomp(r.seq)


def test_specificity_is_the_combination_no_sibling_has_both_sites():
    d = _design()
    f_site = d.forward.seq
    r_site = primers.revcomp(d.reverse.seq)
    for sib in (S1_SEQ, S2_SEQ):
        assert not (f_site in sib and r_site in sib)
    # and combinatorial specificity is not reported as a per-oligo Tm margin
    assert d.delta_tm is None
    assert "LOW_DELTA_TM" not in d.flags
