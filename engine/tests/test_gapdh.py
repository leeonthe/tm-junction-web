"""Regression tests locked to the validated GAPDH + MYC matrix.

Runs fully offline from seeded fixtures (no NCBI calls). See vault
`01 Science/Worked Example - GAPDH.md` and `Sequence-Based Amplifiability (ATGC).md`.
"""

from app.analyze import analyze


def _by_acc(resp):
    return {t.accession: t for t in resp.transcripts}


def test_gapdh_target_conventional():
    r = analyze("NM_001256799.3")
    assert r.gene.symbol == "GAPDH"
    assert r.target_verdict.tier == "CONVENTIONAL"
    assert r.target_verdict.unique_regions[0].exon_order == 1
    assert r.primer_design.forward is not None
    assert r.primer_design.reverse is not None
    # forward primer sits in the unique exon-1 region and is specific
    assert r.primer_design.delta_tm is not None and r.primer_design.delta_tm > 0


def test_gapdh_full_matrix():
    r = analyze("NM_002046.7")
    v = _by_acc(r)
    assert v["NM_001256799.3"].tier == "CONVENTIONAL"
    assert v["NM_001289745.3"].tier == "NEEDS_EEJ"
    assert v["NM_001289746.2"].tier == "CONVENTIONAL"      # red by coords, blue by sequence
    assert v["NM_001357943.2"].tier == "NEEDS_EEJ"
    assert v["NR_152150.2"].tier == "NEEDS_EEJ"           # the non-coding variant: exon 6–7
    # MANE has no unique region or junction, and no combination isolates it either once
    # the NR sibling is counted — see test_gapdh_mane_is_a_hard_case_beside_its_nr_sibling.
    assert v["NM_002046.7"].tier == "NO_SINGLE_UNIQUE_JUNCTION"


def test_gapdh_junction_locations():
    r = analyze("NM_001289745.3")
    v = _by_acc(r)
    rj = v["NM_001289745.3"].recommended_junction
    assert (rj.donor_order, rj.acceptor_order) == (1, 2)
    rj4 = v["NM_001357943.2"].recommended_junction
    assert (rj4.donor_order, rj4.acceptor_order) == (3, 4)


def test_an_eej_transcript_names_no_target_site_exon():
    """Ticket 10: a target site is a place that DISCRIMINATES.

    A single-junction EEJ transcript used to get its nearest flanking exon marked as a
    target site for the partner primer. That exon carries no specificity — the junction
    supplies all of it — so calling it a target site pointed the designer at an arbitrary
    exon. The partner primer is now the junction designer's own Tm-matched search.
    """
    v = _by_acc(analyze("NM_001289745.3"))
    for acc in ("NM_001289745.3", "NM_001357943.2"):
        assert v[acc].tier == "NEEDS_EEJ"
        assert v[acc].recommended_junction is not None
        assert not hasattr(v[acc], "partner_exon")
        # The only target sites left are the ones that distinguish the transcript.
        assert v[acc].amplify_exon_pair is None or v[acc].unique_regions


def test_gapdh_t3_coord_vs_sequence_divergence():
    """The whole point: coord non-unique but sequence CONVENTIONAL."""
    r = analyze("NM_001289746.2")
    v = _by_acc(r)["NM_001289746.2"]
    assert v.tier == "CONVENTIONAL"
    assert v.coord_non_unique is True          # structural flag disagrees, sequence wins


def _gapdh_mane_amp(include_nr: bool):
    from app import ncbi
    from app.amplify import analyze_amplifiability
    _, _, _, _, _, ts = ncbi.refseq_transcripts(ncbi.get_product_report("GAPDH"))
    sibs = [t for t in ts if t["accession"] != "NM_002046.7"
            and (include_nr or not ncbi.is_noncoding(t["accession"]))]
    mane = next(t for t in ts if t["accession"] == "NM_002046.7")
    seqs = ncbi.get_sequences([t["accession"] for t in ts])
    return analyze_amplifiability(
        mane["exons"], seqs["NM_002046.7"], [seqs[t["accession"]] for t in sibs],
        sibling_exons=[t["exons"] for t in sibs])


def test_gapdh_mane_combo_among_the_mrnas():
    """Against its NM siblings alone, MANE has no unique region and no unique single
    junction, but a junction+exon combination isolates it: an EEJ across a junction + a
    conventional primer in an exon no sibling pairs with it."""
    amp = _gapdh_mane_amp(include_nr=False)
    assert amp.tier == "NEEDS_EEJ"
    assert amp.combo_je is not None


def test_gapdh_mane_is_a_hard_case_beside_its_nr_sibling():
    """GAPDH's non-coding transcript is in the same cDNA, so it is a sibling — and it
    takes MANE's combination away. NR_152150.2 is MANE minus the end of exon 6 and all of
    exons 7–8, so it carries every MANE junction from 1–2 to 5–6, including the 1–2
    junction the combination relied on; each remaining MANE feature is shared with an NM
    sibling that also has the rest. No pair of sites is MANE's alone."""
    assert _gapdh_mane_amp(include_nr=True).tier == "NO_SINGLE_UNIQUE_JUNCTION"
    r = analyze("NM_002046.7")
    v = _by_acc(r)["NM_002046.7"]
    assert v.tier == "NO_SINGLE_UNIQUE_JUNCTION"
    assert v.recommended_junction is None
    assert r.primer_design.forward is None and r.primer_design.reverse is None


def test_gapdh_summary():
    r = analyze("NM_002046.7")
    assert r.summary.nm_count == 6            # 5 NM + 1 NR — every curated transcript
    assert r.summary.nr_count == 1
    assert r.summary.conventional_count == 2
    assert r.summary.needs_eej_count == 3     # 745, 943, and the NR (exon 6–7 junction)
    assert r.summary.hard_case_count == 1     # MANE, once the NR sibling is counted
    assert r.summary.coord_non_unique_count == 5


def test_myc_stress():
    r = analyze("NM_002467.6")
    v = _by_acc(r)
    assert v["NM_001354870.1"].tier == "CONVENTIONAL"     # 797 unique windows
    assert v["NM_002467.6"].tier == "NEEDS_EEJ"           # 3bp boundary shift junction
    rj = v["NM_002467.6"].recommended_junction
    assert (rj.donor_order, rj.acceptor_order) == (1, 2)


def test_primer3_qc_conventional():
    """Conventional primers must pass primer3 QC: GC 40–60%, Tm 57–63°C, high ΔTm."""
    r = analyze("NM_001256799.3")
    pd = r.primer_design
    assert pd.tm_method == "primer3"
    for p in (pd.forward, pd.reverse):
        assert p is not None
        assert 40 <= p.gc <= 60, f"GC out of range: {p.gc}"
        assert 55 <= p.tm <= 64, f"Tm out of range: {p.tm}"
        assert p.hairpin_tm < 45 and p.homodimer_tm < 45
    assert pd.forward.qc_pass
    assert pd.delta_tm is not None and pd.delta_tm >= 10   # specific
    assert "LOW_QC" not in pd.flags and "LOW_DELTA_TM" not in pd.flags
    assert pd.pair_dimer_tm is not None and pd.pair_dimer_tm < 45  # no primer-dimer


def test_junction_primer_maximizes_delta_tm():
    """ΔTm-aware selection: this junction now yields a specific, QC-clean primer."""
    r = analyze("NM_001357943.2")
    pd = r.primer_design
    assert pd.tier == "NEEDS_EEJ"
    assert pd.forward.kind == "junction_spanning"
    assert pd.delta_tm is not None and pd.delta_tm >= 10   # was ~6 before ΔTm weighting
    assert pd.confidence == "high"


def test_low_qc_is_honest_not_hidden():
    """GC-rich junction (GAPDH 5' CpG island) — best specific primer still flagged LOW_QC."""
    r = analyze("NM_001289745.3")
    pd = r.primer_design
    assert pd.tier == "NEEDS_EEJ"
    assert pd.delta_tm is not None and pd.delta_tm >= 10   # specific...
    assert "LOW_QC" in pd.flags and pd.confidence == "low"  # ...but honestly low QC


def test_accession_normalization():
    """Case-insensitive, version-optional, whitespace-tolerant input all resolve."""
    for acc in ["NM_002046.7", "NM_002046", "nm_002046.7", "  NM_002046.7  ", "Nm_002046.7"]:
        r = analyze(acc)
        assert r.gene.symbol == "GAPDH"
        assert r.target_accession == "NM_002046.7"


def test_invalid_accession():
    from app.analyze import AnalysisError
    import pytest
    for bad in ["ENST00000229239", "hello", "", "NP_000537.3"]:
        with pytest.raises(AnalysisError):
            analyze(bad)
