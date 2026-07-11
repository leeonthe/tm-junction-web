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
    assert v["NM_002046.7"].tier == "NO_SINGLE_UNIQUE_JUNCTION"


def test_gapdh_junction_locations():
    r = analyze("NM_001289745.3")
    v = _by_acc(r)
    rj = v["NM_001289745.3"].recommended_junction
    assert (rj.donor_order, rj.acceptor_order) == (1, 2)
    rj4 = v["NM_001357943.2"].recommended_junction
    assert (rj4.donor_order, rj4.acceptor_order) == (3, 4)


def test_gapdh_t3_coord_vs_sequence_divergence():
    """The whole point: coord non-unique but sequence CONVENTIONAL."""
    r = analyze("NM_001289746.2")
    v = _by_acc(r)["NM_001289746.2"]
    assert v.tier == "CONVENTIONAL"
    assert v.coord_non_unique is True          # structural flag disagrees, sequence wins


def test_gapdh_mane_hard_case():
    r = analyze("NM_002046.7")
    v = _by_acc(r)["NM_002046.7"]
    assert v.tier == "NO_SINGLE_UNIQUE_JUNCTION"
    assert r.primer_design.tier == "NO_SINGLE_UNIQUE_JUNCTION"
    assert "NO_SINGLE_UNIQUE_JUNCTION" in r.primer_design.flags


def test_gapdh_summary():
    r = analyze("NM_002046.7")
    assert r.summary.nm_count == 5
    assert r.summary.conventional_count == 2
    assert r.summary.needs_eej_count == 2
    assert r.summary.hard_case_count == 1
    assert r.summary.coord_non_unique_count == 4


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
