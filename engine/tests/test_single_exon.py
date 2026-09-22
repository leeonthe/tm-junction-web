"""Single-exon transcripts get primer pairs.

Every product here used to have to cross an exon-exon junction, because a product inside
one exon cannot be told from one amplified off contaminating genomic DNA. That rule
presumes a junction exists. An intronless transcript has none — the rule in yeast, and not
rare elsewhere (human JUN, the histones; 101 cached human genes mix single- and multi-exon
transcripts) — so it was left with a forward primer, no reverse, and no assay at all.

Now: a single-exon transcript gets a pair inside its one exon, still specific where it has
siblings, and flagged SAME_EXON so the page can state the cost (DNase + a no-RT control).
And a gene with a single-exon transcript gets whole-transcript pairs in the stretch its
transcripts share, since no junction-crossing product can include a transcript with no
junction. Transcripts WITH junctions are held to the rule exactly as before.
"""

import pytest

from app import ncbi, panvariant
from app.analyze import analyze
from app.panvariant import _amplifies, _shared_stretches
from app.primers import (
    AMPLICON_MAX, AMPLICON_MIN, MIN_JX_ARM, _can_span, _product_ok, _spans_junction, _usable_exons,
    revcomp,
)


def _inside_one_exon(r):
    """Both primer sites, and everything between them, lie in the target's single exon."""
    d, v = r.primer_design, r.target_verdict
    span = _amplifies(r.target_mrna, d.forward.seq, d.reverse.seq)
    assert span is not None
    assert len(v.exons) == 1 and 0 <= span[0] < span[1] <= v.exons[0].tx_end
    return span


# ---------------------------------------------------------------- the rule

def test_the_junction_rule_binds_only_where_there_is_a_junction():
    one_exon, three_exons = [999], [150, 450, 1000]
    assert _can_span(one_exon, 100, 20) is True             # nothing to reach
    assert _product_ok(three_exons, 10, 140) is False       # inside exon 1: refused
    assert _product_ok(three_exons, 100, 240) is True       # crosses 1|2
    # Same-exon is something the CALLER establishes, never a default.
    assert _product_ok(one_exon, 100, 240) is False
    assert _product_ok(one_exon, 100, 240, same_exon_ok=True) is True
    assert (_usable_exons([999]), _usable_exons([10, 1128]), _usable_exons(three_exons)) == (1, 1, 3)


# ---------------------------------------------------------------- transcript-specific

@pytest.mark.parametrize("accession, symbol", [
    ("NM_001181321.3", "TDH3"),      # baker's yeast — intronless, like most of its genes
    ("NM_002228.4", "JUN"),          # human — the same case
    ("NM_019845.3", "RPRM"),         # human, minus strand — the transcript the ticket names
])
def test_a_sole_single_exon_transcript_gets_a_pair(accession, symbol):
    r = analyze(accession)
    d = r.primer_design
    assert r.gene.symbol == symbol and len(r.target_verdict.exons) == 1
    assert d.forward is not None and d.reverse is not None       # it used to stop at the forward
    assert AMPLICON_MIN <= d.amplicon_len <= AMPLICON_MAX
    start, end = _inside_one_exon(r)
    assert end - start == d.amplicon_len
    # The cost is stated, in the flag the page reads and in the sentence it prints.
    assert "SAME_EXON" in d.flags and "NO_SPANNING_PAIR" not in d.flags
    assert "no-RT control" in d.mechanism


def test_a_single_exon_transcript_with_a_sibling_is_still_specific():
    """Fly Gapdh1: two transcripts, one exon each. The pair for variant A must not amplify
    variant B — same-exon is a concession on gDNA, never on specificity."""
    r = analyze("NM_080369.3")
    d = r.primer_design
    assert [len(t.exons) for t in r.transcripts] == [1, 1]
    assert d.forward is not None and d.reverse is not None and "SAME_EXON" in d.flags
    _inside_one_exon(r)
    sibling = ncbi.get_sequence("NM_001038847.2")
    assert _amplifies(sibling, d.forward.seq, d.reverse.seq) is None
    assert d.forward.seq not in sibling or revcomp(d.reverse.seq) not in sibling


def test_transcripts_with_junctions_are_held_to_the_rule_as_before():
    r = analyze("NM_001256799.3")                                # GAPDH variant 2, 8 exons
    d = r.primer_design
    assert "SAME_EXON" not in d.flags
    cum = [e.tx_end for e in r.target_verdict.exons]
    start, end = _amplifies(r.target_mrna, d.forward.seq, d.reverse.seq)
    assert _product_ok(cum, start, end) and len(cum) > 1
    assert all("SAME_EXON" not in o.flags for o in r.pan_variant_options)


# ---------------------------------------------------------------- whole-transcript

def test_a_sole_single_exon_transcript_gets_a_whole_transcript_pair():
    r = analyze("NM_001181321.3")
    assert r.pan_variant is not None and r.pan_variant_options
    for o in r.pan_variant_options:
        assert o.covered == ["NM_001181321.3"] and not o.uncovered
        assert o.exons == [1, 1] and "SAME_EXON" in o.flags
        assert AMPLICON_MIN <= o.amplicon_len <= AMPLICON_MAX


def test_single_exon_siblings_are_measured_together_in_what_they_share():
    """Fly Gapdh1's two single-exon variants differ at their ends; one pair in the stretch
    between covers both at one size."""
    r = analyze("NM_080369.3")
    best = r.pan_variant_options[0]
    assert sorted(best.covered) == ["NM_001038847.2", "NM_080369.3"] and not best.uncovered
    assert "SAME_EXON" in best.flags
    sizes = {a: _amplifies(ncbi.get_sequence(a), best.forward.seq, best.reverse.seq) for a in best.covered}
    assert {e - s for s, e in sizes.values()} == {best.amplicon_len}


def test_a_single_exon_sibling_is_covered_where_no_junction_crossing_pair_can_reach_it():
    """GHSR: NM_198407.2 has two exons; NM_004122.2 (GHSR1b) is ONE, reading on into the
    intron. No product across the 1|2 junction can include it — so the pair that measures
    the gene sits in the exon-1 stretch both carry, and outranks the junction-crossing pair,
    which reaches only one of the two."""
    r = analyze("NM_198407.2")
    assert sorted(len(t.exons) for t in r.transcripts) == [1, 2]
    best = r.pan_variant_options[0]
    assert sorted(best.covered) == ["NM_004122.2", "NM_198407.2"] and "SAME_EXON" in best.flags
    crossing = [o for o in r.pan_variant_options if "SAME_EXON" not in o.flags]
    assert crossing and all(o.uncovered == ["NM_004122.2"] for o in crossing)
    assert r.pan_variant_options.index(crossing[0]) > 0          # coverage ranks first
    # Verified from sequence, one size in both.
    for a in best.covered:
        s, e = _amplifies(ncbi.get_sequence(a), best.forward.seq, best.reverse.seq)
        assert e - s == best.amplicon_len
    # And the two-exon transcript's OWN assay still crosses its junction.
    assert "SAME_EXON" not in r.primer_design.flags


# ---------------------------------------------------------------- the shared stretch

def _stretch_of(ref, exons, strands):
    """The proposal made from `ref` — every single-exon transcript makes its own."""
    return next(r for r in _shared_stretches(exons, list(exons), strands) if r[0] == ref)


def test_the_shared_stretch_is_what_every_carrier_covers():
    exons = {"ONE": [(1000, 1999)], "TWO": [(900, 1499), (3000, 3400)], "FAR": [(9000, 9500)]}
    _ref, lo, hi, carriers = _stretch_of("ONE", exons, {"ONE": "+"})
    assert carriers == ["ONE", "TWO"]                            # FAR shares nothing
    assert (lo, hi) == (0, 500)                                  # genomic 1000-1499, on ONE's mRNA


def test_the_shared_stretch_is_read_from_the_5_prime_end_on_the_minus_strand():
    exons = {"ONE": [(1000, 1999)], "TWO": [(900, 1499)]}
    _ref, lo, hi, _c = _stretch_of("ONE", exons, {"ONE": "-"})
    assert (lo, hi) == (500, 1000)                               # the mRNA starts at 1999


def test_a_sliver_of_overlap_does_not_shrink_the_stretch_to_nothing():
    """A transcript sharing less than a product's worth is left out, not allowed to narrow
    the stretch until no pair fits for anyone."""
    exons = {"ONE": [(1000, 1999)], "BIG": [(1000, 1800)], "TINY": [(1990, 2500)]}
    _ref, lo, hi, carriers = _stretch_of("ONE", exons, {"ONE": "+"})
    assert carriers == ["ONE", "BIG"] and (lo, hi) == (0, 801)
    # No single-exon transcript, no same-exon proposal at all.
    assert _shared_stretches({"A": [(1, 500), (900, 1400)]}, ["A"], {}) == []


def test_same_exon_pairs_are_proposed_only_for_genes_with_a_single_exon_transcript():
    ts = ncbi.refseq_transcripts(ncbi.get_product_report("GAPDH")).transcripts
    tmap = {t["accession"]: t for t in ts}
    seqs = ncbi.get_sequences(list(tmap))
    assert all("SAME_EXON" not in o.flags for o in panvariant.design_options(tmap, seqs))


# ---------------------------------------------------------------- strand, with no exon order

def test_a_single_exon_transcript_takes_its_strand_from_ncbi_not_from_exon_order():
    """Strand was inferred from exon order, which one exon does not have, and defaulted to
    plus — mirroring every position inside the exon on a minus-strand gene. Fly Gapdh1's
    variant A is unique over mRNA 1-150 (its 5' end, the exon's HIGH coordinates) and was
    reported as unique over 1331-1480. The designs were right; the sentence was not."""
    from app.analyze import _is_plus, _region_tx, _uniq_genomic
    one = [(1000, 1999)]
    assert _is_plus(one, "+") and _is_plus(one, "") and not _is_plus(one, "-")
    assert _is_plus([(5000, 5100), (1000, 1100)], "+") is False     # two exons: order decides
    assert _region_tx(one, 1, 1850, 1999, "-") == (1, 150)
    assert _region_tx(one, 1, 1850, 1999, "+") == (851, 1000)
    assert _uniq_genomic(one, 0, 149, 1, "-") == (1850, 1999)

    r = analyze("NM_080369.3")
    assert r.gene.strand == "-" and len(r.target_verdict.exons) == 1
    u = r.target_verdict.unique_regions[0]
    assert (u.tx_begin, u.tx_end) == (1, 150)
    assert u.window_starts[0][0] == 1                 # where the specific primers really start
    exon = r.target_verdict.exons[0]
    assert (u.begin, u.end) == (exon.end - 149, exon.end)


def test_the_capability_is_advertised():
    from app.models import FEATURES
    assert "single_exon_pairs" in FEATURES


# ---------------------------------------------------------------- an exon too short for a primer

def test_a_too_short_exon_is_a_partial_binding_site_so_the_product_still_crosses():
    """Yeast ACT1: a 10-nt first exon, then 1118. No primer fits in 10 nt — but one can START
    there and run on into exon 2, and then the product crosses the junction and genomic DNA
    cannot give it. That beats a same-exon pair, so it is what ACT1 gets, unflagged."""
    r = analyze("NM_001179927.1")
    d, cum = r.primer_design, [e.tx_end for e in r.target_verdict.exons]
    assert cum == [10, 1128]
    assert d.forward is not None and d.reverse is not None       # it used to stop at the forward
    start, end = _amplifies(r.target_mrna, d.forward.seq, d.reverse.seq)
    assert _spans_junction(cum, start, end) and "SAME_EXON" not in d.flags
    assert AMPLICON_MIN <= end - start <= AMPLICON_MAX
    # The forward primer has a real foothold on each side of the junction.
    assert cum[0] - start >= MIN_JX_ARM and start + d.forward.length - cum[0] >= MIN_JX_ARM
    # The whole-transcript pairs do the same.
    best = r.pan_variant_options[0]
    assert best.exons == [1, 2] and "SAME_EXON" not in best.flags and best.covered == [r.target_accession]
    assert best.forward.tx_start <= cum[0] - MIN_JX_ARM


def test_a_junction_in_reach_means_a_legal_product_is_in_reach():
    """What kept ACT1 without a pair: a specific primer 45 nt in can only make a 66 bp
    product with a partner upstream of it, under the 70 bp floor — yet the pre-check called
    the junction reachable, so that primer was chosen and no partner existed."""
    cum = [10, 1128]
    assert _can_span(cum, 45, 21) is False            # 66 bp at most: not a product
    assert _can_span(cum, 60, 20) is True             # 80 bp with a partner at base 0
    assert _can_span(cum, 400, 20) is False           # 220 bp cannot stretch back to base 5
    assert _can_span([150, 450, 1000], 100, 20) is True           # ordinary exons, as before


def test_both_primers_go_in_the_long_exon_only_when_nothing_crosses():
    """The fallback, built so that no junction-crossing pair can exist: a 4-nt exon cannot
    give a primer its MIN_JX_ARM foothold. Then — and only then — the pair sits in the long
    exon and says so."""
    from app.amplify import analyze_amplifiability
    from app import primers
    from tests.test_pan_variant import _bases
    seq = _bases(11, 600)
    exons = [(1000, 1003), (2000, 2595)]                           # 4 nt, then 596
    amp = analyze_amplifiability(exons, seq, [], sibling_exons=[])
    d = primers.design(exons, seq, {}, amp)
    assert d.forward is not None and d.reverse is not None
    start, end = _amplifies(seq, d.forward.seq, d.reverse.seq)
    assert not _spans_junction([4, 600], start, end)
    assert "SAME_EXON" in d.flags and "too short to hold a primer" in d.mechanism
