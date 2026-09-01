"""One pair for the whole gene — ticket 25.

The opposite of every other design here: it must amplify as MANY of a gene's transcripts as
possible rather than exactly one. Three requirements, and the second is the one that makes
this hard:

  * both primer sites present in every transcript claimed,
  * ONE product size across all of them — a second size is a second band, and a band per
    isoform makes a total-expression assay unquantifiable,
  * at least two exons, so the product cannot be confused with one off genomic DNA.

Coverage is credited from SEQUENCE — both sites located in each transcript and the product
measured — not from the exon bookkeeping that proposes where to look.
"""

from app import ncbi, panvariant
from app.analyze import analyze
from app.panvariant import _amplifies, _consecutive_index, _coverage
from app.primers import revcomp


def _gene(symbol):
    _, _, _, _, _, ts = ncbi.nm_transcripts(ncbi.get_product_report(symbol))
    seqs = {t["accession"]: ncbi.get_sequence(t["accession"]) for t in ts}
    return {t["accession"]: t for t in ts}, seqs


def test_amplifies_locates_the_product_or_refuses():
    seq = "AAAA" + "ACGTACGTACGTACGT" + "TTTT" + "GGGGCCCCGGGGCCCC" + "AAAA"
    fwd = "ACGTACGTACGTACGT"
    rev = revcomp("GGGGCCCCGGGGCCCC")
    assert _amplifies(seq, fwd, rev) == (4, 40)
    # The reverse primer must sit downstream of the forward, not before it.
    assert _amplifies(seq, revcomp("GGGGCCCCGGGGCCCC"), revcomp(fwd)) is None
    # A site absent from the transcript is not amplified.
    assert _amplifies(seq, "CCCCCCCCCCCCCCCC", rev) is None


def test_a_second_binding_site_disqualifies_the_transcript():
    """Two sites make the product ambiguous — the same problem as a second band."""
    repeated = "ACGTACGTACGTACGT"
    seq = "AA" + repeated + "TT" + repeated + "GG" + "GGGGCCCCGGGGCCCC"
    assert _amplifies(seq, repeated, revcomp("GGGGCCCCGGGGCCCC")) is None


def test_coverage_counts_only_transcripts_that_amplify_at_ONE_size():
    """The size carrying the most transcripts wins; the odd length out is not covered."""
    fwd, rev_site = "ACGTACGTACGTACGT", "GGGGCCCCGGGGCCCC"
    rev = revcomp(rev_site)
    short = fwd + "T" * 40 + rev_site           # 72 nt product
    long_ = fwd + "T" * 90 + rev_site           # 122 nt product
    seqs = {"A": "CC" + short, "B": "GG" + short, "C": short, "D": "AA" + long_}
    covered, uncovered, size = _coverage((fwd, rev), seqs, ["A", "B", "C", "D"])
    assert covered == ["A", "B", "C"]
    assert uncovered == ["D"]
    assert size == 72


def test_a_run_must_be_consecutive_not_merely_present():
    """An exon spliced in between makes a longer product from the same sites."""
    hay = [(1, 10), (20, 30), (40, 50), (60, 70)]
    assert _consecutive_index(hay, ((20, 30), (40, 50))) == 1
    assert _consecutive_index(hay, ((20, 30), (60, 70))) is None      # exon 3 sits between
    assert _consecutive_index(hay, ((99, 100),)) is None


def test_gapdh_gets_one_pair_for_every_variant():
    tmap, seqs = _gene("GAPDH")
    d = panvariant.design(tmap, seqs)
    assert d is not None
    assert len(d.covered) == len(tmap) == 5
    assert not d.uncovered
    # Every claim re-checked from the sequences, independently of how it was designed.
    sizes = {a: _amplifies(seqs[a], d.forward.seq, d.reverse.seq) for a in tmap}
    assert all(sizes.values())
    assert {e - b for b, e in sizes.values()} == {d.amplicon_len}


def test_tp53_all_25_accessions_from_one_pair():
    """The gene the fold came from: 25 accessions, 13 transcripts, one pair for all."""
    tmap, seqs = _gene("TP53")
    d = panvariant.design(tmap, seqs)
    assert len(d.covered) == 25
    assert not d.uncovered


def test_the_product_always_crosses_a_junction():
    """The two-exon rule, on the gene-wide pair as much as the isoform-specific one."""
    from app.amplify import cumulative_exon_ends
    from app.primers import _spans_junction
    for symbol in ("GAPDH", "TP53", "MYC"):
        tmap, seqs = _gene(symbol)
        d = panvariant.design(tmap, seqs)
        cum = cumulative_exon_ends([tuple(e) for e in tmap[d.reference]["exons"]])
        span = _amplifies(seqs[d.reference], d.forward.seq, d.reverse.seq)
        assert span is not None
        assert _spans_junction(cum, span[0], span[1]), symbol
        assert d.exons[0] != d.exons[1], symbol


def test_coverage_is_sequence_evidence_not_exon_bookkeeping():
    """MYC's two transcripts share no identical exon run, yet one pair measures both.

    Their difference lies outside the primer sites, so exon bookkeeping alone would credit
    the pair with one transcript and hand the user a worse assay than exists.
    """
    tmap, seqs = _gene("MYC")
    exons = {a: [tuple(e) for e in t["exons"]] for a, t in tmap.items()}
    (a1, a2) = list(exons)
    shared_runs = [k for k in (tuple(exons[a1][i:j + 1])
                               for i in range(len(exons[a1]))
                               for j in range(i + 1, len(exons[a1])))
                   if _consecutive_index(exons[a2], k) is not None]
    assert not shared_runs                       # no identical run between the two
    d = panvariant.design(tmap, seqs)
    assert len(d.covered) == 2                   # ...and yet


def test_the_response_carries_it():
    r = analyze("NM_002046.7")
    p = r.pan_variant
    assert p is not None and p.forward and p.reverse
    assert len(p.covered) == len(r.transcripts)
    assert p.amplicon_len > 0 and p.exons[0] != p.exons[1]


# ---------------------------------------------------------------- constructed genes
#
# Every gene in the cache admits a pair covering all of its variants, so the two cases the
# ticket cares about most — a cassette exon that changes the SIZE, and falling back to the
# largest coverable set — are built here rather than hoped for.

def _bases(seed: int, n: int) -> str:
    """Deterministic, GC-balanced filler primer3 can actually design in."""
    out, x = [], seed
    for _ in range(n):
        x = (x * 1103515245 + 12345) & 0x7FFFFFFF
        out.append("ACGT"[(x >> 16) % 4])
    return "".join(out)


EX1, EX2, EX3 = _bases(1, 300), _bases(2, 300), _bases(3, 300)
CASSETTE = _bases(4, 120)


def _tx(exons):
    return {"exons": [list(e) for e in exons], "is_mane": False}


def test_a_cassette_exon_between_the_primers_is_not_covered():
    """Same two primer sites, longer product — a second band, so it does NOT count.

    This is the requirement that makes the feature more than "primers in a shared exon".
    The third transcript carries BOTH sites; a coverage test that only asked "are the sites
    present?" would claim it and hand the user an assay that cannot be quantified. Two
    exons only, so the sole junction available is the one the cassette sits in — the case
    cannot be dodged by placing the pair elsewhere.
    """
    a, c = (1000, 1299), (1800, 2099)
    cassette = (1500, 1619)
    tmap = {
        "NM_100001.1": _tx([a, c]),
        "NM_100002.1": _tx([a, c]),
        "NM_100003.1": _tx([a, cassette, c]),         # +120 nt between the primer sites
    }
    seqs = {
        "NM_100001.1": EX1 + EX3,
        "NM_100002.1": EX1 + EX3,
        "NM_100003.1": EX1 + CASSETTE + EX3,
    }
    d = panvariant.design(tmap, seqs)
    assert d is not None
    assert d.covered == ["NM_100001.1", "NM_100002.1"]
    assert d.uncovered == ["NM_100003.1"]
    # Excluded on SIZE, not on absence: both sites are there, 120 nt further apart.
    third = _amplifies(seqs["NM_100003.1"], d.forward.seq, d.reverse.seq)
    assert third is not None
    assert third[1] - third[0] == d.amplicon_len + 120


def test_falls_back_to_the_most_common_region():
    """No region reaches all three, so the pair covering two is the right answer."""
    a, b = (1000, 1299), (1400, 1699)
    other = (5000, 5299), (5400, 5699)
    tmap = {
        "NM_200001.1": _tx([a, b]),
        "NM_200002.1": _tx([a, b]),
        "NM_200003.1": _tx(list(other)),             # shares nothing with the other two
    }
    seqs = {
        "NM_200001.1": EX1 + EX2,
        "NM_200002.1": EX1 + EX2,
        "NM_200003.1": _bases(9, 300) + _bases(10, 300),
    }
    d = panvariant.design(tmap, seqs)
    assert d is not None
    assert d.covered == ["NM_200001.1", "NM_200002.1"]
    assert d.uncovered == ["NM_200003.1"]
    assert d.exons[0] != d.exons[1]                  # still two exons


def test_no_pair_at_all_is_reported_as_none():
    """A single-exon transcript has no junction to cross, so there is nothing to offer."""
    tmap = {"NM_300001.1": _tx([(1000, 1299)])}
    assert panvariant.design(tmap, {"NM_300001.1": EX1}) is None
