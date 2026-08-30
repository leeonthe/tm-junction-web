"""Several accessions, one molecule — TP53's duplicated NM records.

RefSeq mints more than one NM accession for the same mRNA: TP53 has 25 NM accessions
covering 13 distinct sequences, e.g. NM_001126115.2 and NM_001276697.3, both annotated
"transcript variant 5" and byte-identical. Treated as separate isoforms, each transcript
is compared against its own twin, no primer can tell the two apart, and the tool calls a
designable transcript EEJ-infeasible — 24 of 25 TP53 isoforms used to land there.

Ticket 27: fold accessions by SEQUENCE (and by exon structure where the gene picker has
no sequences), and name the folded accessions on the surviving verdict.
"""

from app.analyze import _pick_representative, _same_sequence_groups, _same_structure_groups


def acc(name, mane=False, exons=((1, 100), (200, 300))):
    return {"accession": name, "is_mane": mane, "exons": [list(e) for e in exons]}


def accs(*ts):
    return {t["accession"]: t for t in ts}


def test_identical_sequences_collapse_to_one_transcript():
    a = accs(acc("NM_000001.1"), acc("NM_000002.1"), acc("NM_000003.1"))
    seqs = {"NM_000001.1": "ACGT", "NM_000002.1": "ACGT", "NM_000003.1": "TTTT"}
    reps, same = _same_sequence_groups(a, seqs)
    assert reps == ["NM_000001.1", "NM_000003.1"]
    assert same["NM_000001.1"] == ["NM_000002.1"]
    assert same["NM_000003.1"] == []


def test_the_analyzed_accession_represents_its_own_group():
    """The user asked about their accession, not about a synonym of it."""
    a = accs(acc("NM_000001.1", mane=True), acc("NM_000002.1"))
    seqs = {"NM_000001.1": "ACGT", "NM_000002.1": "ACGT"}
    reps, same = _same_sequence_groups(a, seqs, target_acc="NM_000002.1")
    assert reps == ["NM_000002.1"]
    assert same["NM_000002.1"] == ["NM_000001.1"]


def test_mane_speaks_for_the_group_when_the_target_is_elsewhere():
    a = accs(acc("NM_000009.1"), acc("NM_000002.1", mane=True))
    seqs = {"NM_000009.1": "ACGT", "NM_000002.1": "ACGT"}
    reps, _ = _same_sequence_groups(a, seqs, target_acc="NM_999999.9")
    assert reps == ["NM_000002.1"]
    # With no MANE either, the choice still has to be stable between runs.
    assert _pick_representative(["NM_000009.1", "NM_000002.1"], a, None) == "NM_000002.1"


def test_sequences_decide_it_accessions_do_not():
    """Same variant NAME, different sequence: still two transcripts."""
    a = accs(acc("NM_000001.1"), acc("NM_000002.1"))
    seqs = {"NM_000001.1": "ACGT", "NM_000002.1": "ACGTA"}
    reps, same = _same_sequence_groups(a, seqs)
    assert reps == ["NM_000001.1", "NM_000002.1"]
    assert same == {"NM_000001.1": [], "NM_000002.1": []}


def test_the_gene_picker_folds_on_exon_structure():
    """/gene fetches no sequences; identical exon structures splice to identical mRNA."""
    same_struct = ((1, 100), (200, 300))
    other = ((1, 100), (250, 300))
    reps, same = _same_structure_groups([
        acc("NM_000001.1", exons=same_struct),
        acc("NM_000002.1", exons=same_struct),
        acc("NM_000003.1", exons=other),
    ])
    assert [t["accession"] for t in reps] == ["NM_000001.1", "NM_000003.1"]
    assert same["NM_000001.1"] == ["NM_000002.1"]


def test_tp53_variant_5_is_one_transcript_with_a_design():
    """The reported case, end to end: NM_001126115.2 == NM_001276697.3."""
    from app.analyze import analyze
    r = analyze("NM_001126115.2")
    assert r.target_verdict.same_sequence_accessions == ["NM_001276697.3"]
    # 25 accessions, 13 sequences — and the table lists each sequence once.
    assert r.summary.nm_count == 13
    assert r.summary.merged_accession_count == 12
    assert len(r.transcripts) == 13
    assert len({v.accession for v in r.transcripts}) == 13
    # Its twin is no longer a sibling to be told apart, so the transcript is designable.
    assert r.target_verdict.tier == "NEEDS_EEJ"
    # Every folded accession is named somewhere, so nothing silently disappears.
    listed = {v.accession for v in r.transcripts}
    listed |= {a for v in r.transcripts for a in v.same_sequence_accessions}
    assert len(listed) == 25


def test_every_transcript_carries_ncbis_variant_designation():
    """Ticket 27.1: the accession alone does not say WHICH isoform it is.

    NCBI numbers a variant only when the gene has more than one to tell apart, so an
    absent name is the mono-isoform case rather than missing data — the clients render
    that as "mono-isoform".
    """
    from app.analyze import analyze, lookup_gene
    r = analyze("NM_001126115.2")
    by_acc = {v.accession: v for v in r.transcripts}
    # The reported pair: one row, one variant, both accessions named on it.
    assert by_acc["NM_001126115.2"].variant == "transcript variant 5"
    assert by_acc["NM_001126115.2"].same_sequence_accessions == ["NM_001276697.3"]
    assert by_acc["NM_000546.6"].variant == "transcript variant 1"
    assert all(v.variant for v in r.transcripts)      # TP53 numbers every isoform
    # The gene picker carries it too, so a variant is named before anything is analyzed.
    assert all(t.variant for t in lookup_gene("TP53").transcripts)


def test_a_sole_isoform_has_no_variant_number():
    """ACTB has one NM, so NCBI names no variant — there is nothing to number it against."""
    from app.analyze import analyze
    r = analyze("NM_001101.5")
    assert len(r.transcripts) == 1
    assert r.transcripts[0].variant is None
