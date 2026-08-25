"""Every designed amplicon must cross at least one exon-exon junction.

A product contained inside ONE exon is indistinguishable from one amplified off
contaminating genomic DNA: the same primer sites sit uninterrupted in the genome, so a
gDNA template gives the same band. Spanning a junction makes gDNA either fail to amplify
or give a visibly longer product, which is why intron-spanning design is standard practice
for RT-PCR -- and why it is a correctness requirement here, not a preference.

Ticket 22: the conventional unique-region path could place both primers in one exon,
because the partner search only ever checked amplicon LENGTH.
"""

from app.amplify import cumulative_exon_ends
from app.analyze import analyze
from app.primers import _exon_of, _spans_junction

CUM = [150, 330, 520, 700, 950, 1200]


def test_spans_junction_reads_the_exon_boundaries():
    assert _spans_junction(CUM, 10, 140) is False      # wholly inside exon 1
    assert _spans_junction(CUM, 160, 300) is False     # wholly inside exon 2
    assert _spans_junction(CUM, 140, 160) is True      # straddles 1|2
    assert _spans_junction(CUM, 10, 900) is True       # several exons


def test_a_product_ending_on_a_boundary_is_still_one_exon():
    """[0,150) is exon 1 entire; the junction is crossed only once base 150 is included."""
    assert _spans_junction(CUM, 0, 150) is False
    assert _spans_junction(CUM, 0, 151) is True
    assert _exon_of(CUM, 149) == 0 and _exon_of(CUM, 150) == 1


def test_unknown_exon_structure_passes_rather_than_rejecting_everything():
    """Nothing to check against -> offer the pair; silently returning no design would be
    a worse failure than an unverified one."""
    assert _spans_junction([], 10, 40) is True


def _check(accession: str):
    r = analyze(accession)
    d = r.primer_design
    if not (d.forward and d.reverse):
        return None
    exons = [(e.begin, e.end) for e in r.target_verdict.exons]
    cum = cumulative_exon_ends(exons)
    start = d.forward.tx_start
    end = d.reverse.tx_start + d.reverse.length
    return _spans_junction(cum, start, end), start, end, cum


def test_real_designs_span_a_junction():
    """The genes the suite already exercises, checked against their own exon structures."""
    for acc in ("NM_001641.4", "NM_001256799.3", "NM_002046.7", "NM_001101.5"):
        got = _check(acc)
        if got is None:
            continue                       # no pair designed for this tier — nothing to check
        spans, start, end, cum = got
        assert spans, f"{acc}: amplicon {start}-{end} sits inside one exon (cum={cum})"


def test_actb_the_mono_isoform_case_still_gets_a_pair():
    """ACTB has ONE NM transcript and a 1.8 kb final exon — exactly where an intra-exon
    product was easiest to produce. It must still yield a design, and that design must span."""
    got = _check("NM_001101.5")
    assert got is not None, "ACTB lost its design entirely"
    assert got[0]
