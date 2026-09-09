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


def test_three_or_more_accessions_leave_exactly_one_row():
    """The ticket's "if more than two" case, which no cached gene actually exhibits.

    Keep one, remove ALL the others, and name every one of them under the kept transcript —
    so the row still accounts for each accession the gene has.
    """
    same = ((1, 100), (200, 300))
    a = accs(acc("NM_000004.1", exons=same), acc("NM_000002.1", exons=same, mane=True),
             acc("NM_000003.1", exons=same), acc("NM_000009.1", exons=((1, 100), (250, 350))))
    seqs = {"NM_000004.1": "ACGT", "NM_000002.1": "ACGT", "NM_000003.1": "ACGT",
            "NM_000009.1": "TTTT"}
    reps, folded = _same_sequence_groups(a, seqs)
    assert reps == ["NM_000002.1", "NM_000009.1"]                 # MANE speaks for the trio
    assert folded["NM_000002.1"] == ["NM_000003.1", "NM_000004.1"]
    # Nothing lost: every accession is either a row or named under one.
    assert set(reps) | {x for v in folded.values() for x in v} == set(a)


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


def test_identical_sequence_is_not_enough_without_the_structure():
    """Both criteria, because both are what the claim rests on.

    Same molecule means the same mRNA AND the same exons at the same coordinates. A
    coincidence of sequence at a different structure is not the same transcript, and
    folding it away would hide a real isoform behind another's accession.
    """
    a = accs(acc("NM_000001.1", exons=((1, 100), (200, 300))),
             acc("NM_000002.1", exons=((1, 100), (250, 350))))
    seqs = {"NM_000001.1": "ACGT", "NM_000002.1": "ACGT"}
    reps, same = _same_sequence_groups(a, seqs)
    assert reps == ["NM_000001.1", "NM_000002.1"]
    assert same == {"NM_000001.1": [], "NM_000002.1": []}


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


def test_the_variant_is_read_off_the_record_when_the_gene_report_omits_it():
    """Ticket 27.1: "can be found in NCBI database when Refseq ID is searched".

    The gene's product report names each isoform, but that is one surface of NCBI and it
    can be silent. The record you get by searching the accession carries the same
    designation in its title, so it is the second source — and the picker needs it as much
    as the verdict table, since that is where the numbers were missing.
    """
    import copy
    from app import ncbi
    from app.analyze import analyze, lookup_gene

    real = ncbi.get_product_report

    def without_names(symbol):
        d = copy.deepcopy(real(symbol))
        for r in d.get("reports") or []:
            for tr in (r.get("product") or {}).get("transcripts") or []:
                tr.pop("name", None)
        return d

    ncbi.get_product_report = without_names
    try:
        picker = {t_.accession: t_.variant for t_ in lookup_gene("GAPDH").transcripts}
        table = {v.accession: v.variant for v in analyze("NM_002046.7").transcripts}
    finally:
        ncbi.get_product_report = real

    # The ticket's own examples, recovered with the report saying nothing.
    for got in (picker, table):
        assert got["NM_002046.7"] == "transcript variant 1"
        assert got["NM_001289745.3"] == "transcript variant 3"
        assert all(got.values())


def test_a_title_names_its_variant_or_it_does_not():
    from app.ncbi import variant_in
    assert variant_in("Homo sapiens ... (GAPDH), transcript variant 1, mRNA") == "transcript variant 1"
    assert variant_in("Homo sapiens ... (TP53), transcript variant X2, mRNA") == "transcript variant X2"
    assert variant_in("Homo sapiens beta-actin (ACTB), mRNA") is None
    assert variant_in(None) is None


def test_a_sole_transcript_is_left_alone():
    """No variant to number, so nothing is fetched and nothing is invented."""
    from app import ncbi
    one = [{"accession": "NM_001101.5"}]
    ncbi.fill_variants(one)
    assert one[0].get("variant") is None


def test_the_fold_is_name_blind_and_keeps_both_names():
    """Ticket 32: one molecule, two accessions, two DIFFERENT variant names.

    The grouping keys on sequence + exon structure and never reads the name, so a
    beta2/beta3-style pair folds regardless — and the folded accession's own designation
    rides along, aligned, so the row can show both. (The ticket's live example,
    NM_001330091.2/NM_001330092.2, is NOT such a pair on current NCBI records: their
    exon 6 acceptors differ by 9 nt, so they are distinguishable and correctly kept
    apart. The rule still has to hold for pairs that ARE identical.)
    """
    same = ((1, 100), (200, 300))
    a = accs(acc("NM_000001.1", exons=same), acc("NM_000002.1", exons=same))
    a["NM_000001.1"]["variant"] = "transcript variant beta2"
    a["NM_000002.1"]["variant"] = "transcript variant beta3"
    seqs = {"NM_000001.1": "ACGT", "NM_000002.1": "ACGT"}
    reps, folded = _same_sequence_groups(a, seqs)
    assert reps == ["NM_000001.1"]
    assert folded["NM_000001.1"] == ["NM_000002.1"]
    # What analyze passes through to the response, aligned with the accession list.
    assert [a[x]["variant"] for x in folded["NM_000001.1"]] == ["transcript variant beta3"]


def test_the_ticketed_nrxn1_pair_is_genuinely_different_and_stays_apart():
    """Guards the premise check: beta2/beta3 differ at exon 6's acceptor (9 nt), so the
    fold must NOT combine them — each has junction k-mers the other lacks."""
    from app import ncbi
    _, _, _, _, _, ts = ncbi.nm_transcripts(ncbi.get_product_report("NRXN1"))
    a = next(t_ for t_ in ts if t_["accession"] == "NM_001330091.2")
    b = next(t_ for t_ in ts if t_["accession"] == "NM_001330092.2")
    assert a["exons"] != b["exons"]
    assert ncbi.get_sequence(a["accession"]) != ncbi.get_sequence(b["accession"])
    reps, _f = _same_sequence_groups(
        {t_["accession"]: t_ for t_ in ts},
        {t_["accession"]: ncbi.get_sequence(t_["accession"]) for t_ in ts})
    assert "NM_001330091.2" in reps and "NM_001330092.2" in reps
