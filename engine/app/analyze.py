"""Orchestration: accession -> gene -> siblings -> per-transcript verdict + primers.

Ties together ncbi (data), overlap (structural), amplify (sequence, authoritative),
and primers (Tm design) into the AnalyzeResponse contract.
"""

from __future__ import annotations

from . import ncbi, overlap, primers
from .amplify import AmplifyResult, analyze_amplifiability, cumulative_exon_ends
from .models import (
    AnalyzeResponse, Exon, GeneInfo, GeneSummary, JunctionOut, PrimerDesignOut,
    PrimerOut, TranscriptVerdict, UniqueRegionOut,
)


class AnalysisError(Exception):
    def __init__(self, code: str, message: str):
        self.code = code
        self.message = message
        super().__init__(message)


def _cds_status(tx_begin: int, tx_end: int, cds: tuple[int, int] | None) -> str:
    if not cds:
        return "noncoding"
    cb, ce = cds
    if tx_end < cb:
        return "5utr"
    if tx_begin > ce:
        return "3utr"
    return "cds"


def _exons_out(exons: list[tuple[int, int]], seq: str, cds: tuple[int, int] | None,
               amp: AmplifyResult) -> list[Exon]:
    cum = cumulative_exon_ends(exons)   # 0-based exclusive running mRNA ends
    out: list[Exon] = []
    for i, (b, e) in enumerate(exons):
        tx_begin = (cum[i - 1] if i else 0) + 1     # 1-based mRNA position
        tx_end = cum[i]
        sub = seq[tx_begin - 1:tx_end]
        gc = round(100.0 * sum(c in "GC" for c in sub) / len(sub), 0) if sub else 0.0
        out.append(Exon(
            order=i + 1, begin=b, end=e, length=e - b + 1,
            tx_begin=tx_begin, tx_end=tx_end,
            cds=_cds_status(tx_begin, tx_end, cds), gc=gc,
            unique_sites=len(amp.internal_starts.get(i + 1, [])),
        ))
    return out


def _junction_label(donor: int, acceptor: int) -> str:
    return f"exon {donor}–exon {acceptor}"


def _amp_to_verdict(acc: str, is_mane: bool, exons: list[tuple[int, int]],
                    seq: str, cds: tuple[int, int] | None,
                    amp: AmplifyResult, coord_non_unique: bool) -> TranscriptVerdict:
    junctions = [
        JunctionOut(donor_order=j.donor_order, acceptor_order=j.acceptor_order,
                    label=_junction_label(j.donor_order, j.acceptor_order))
        for j in amp.unique_junctions
    ]
    recommended = junctions[0] if (amp.needs_eej and junctions) else None
    return TranscriptVerdict(
        accession=acc,
        is_mane=is_mane,
        tier=amp.tier,
        amplifiable=amp.amplifiable,
        needs_eej=amp.needs_eej,
        unique_regions=[UniqueRegionOut(exon_order=r.exon_order, window_count=r.window_count,
                                        side=r.side) for r in amp.unique_regions],
        unique_junctions=junctions,
        recommended_junction=recommended,
        coord_non_unique=coord_non_unique,
        exons=_exons_out(exons, seq, cds, amp),
        amplify_exon_pair=list(amp.exon_pair) if amp.exon_pair else None,
    )


def _primer_out(p) -> PrimerOut | None:
    if p is None:
        return None
    return PrimerOut(seq=p.seq, tm=p.tm, gc=p.gc, length=p.length, kind=p.kind,
                     role=p.role, anchor=p.anchor, hairpin_tm=p.hairpin_tm,
                     homodimer_tm=p.homodimer_tm, qc_pass=p.qc_pass, tx_start=p.tx_start)


def analyze(accession: str, k: int = 20) -> AnalyzeResponse:
    """Run the full analysis and return the response (non-streaming wrapper)."""
    result: AnalyzeResponse | None = None
    for ev in analyze_events(accession, k):
        if ev.get("type") == "result":
            result = ev["result"]
    assert result is not None
    return result


def analyze_events(accession: str, k: int = 20):
    """Generator yielding real progress events, then the result.

    Events: {"type":"progress","pct":int,"detail":str} … {"type":"result","result":AnalyzeResponse}.
    The sequence fetches (one per NM isoform) are the bulk of the wall-clock on a cold cache,
    so progress is driven mostly by those. Raises AnalysisError for bad input / not-found.
    """
    accession = accession.strip().upper()   # RefSeq accessions are uppercase
    if not ncbi.is_valid_nm(accession):
        raise AnalysisError("NOT_NM", "Enter a curated NM RefSeq accession, e.g. NM_002046.7.")

    yield {"type": "progress", "pct": 3, "detail": "Resolving gene…"}
    try:
        symbol, _gid = ncbi.resolve_accession(accession)
        report = ncbi.get_product_report(symbol)
        gene_id, symbol, description, chromosome, transcripts = ncbi.nm_transcripts(report)
    except ncbi.NotFound as e:
        raise AnalysisError("NOT_FOUND", str(e))

    accs = {t["accession"]: t for t in transcripts}
    target_acc = accession
    if target_acc not in accs:
        base = accession.split(".")[0]
        match = next((a for a in accs if a.split(".")[0] == base), None)
        if match:
            target_acc = match
    if target_acc not in accs:
        raise AnalysisError("NOT_FOUND", f"{accession} is not an NM transcript of {symbol}.")

    n = len(accs)
    yield {"type": "progress", "pct": 10, "detail": f"{symbol}: {n} NM isoform{'s' if n != 1 else ''}"}

    # fetch sequences (cache-first) — the main cost; emit progress per transcript
    seqs: dict[str, str] = {}
    for i, a in enumerate(accs):
        seqs[a] = ncbi.get_sequence(a)
        yield {"type": "progress", "pct": 10 + round(68 * (i + 1) / n),
               "detail": f"Fetching mRNA sequences {i + 1}/{n}"}

    exons_by_acc = {a: t["exons"] for a, t in accs.items()}

    yield {"type": "progress", "pct": 82, "detail": "Classifying isoforms…"}
    verdicts: list[TranscriptVerdict] = []
    for a, t in accs.items():
        sib_accs = [o for o in accs if o != a]
        amp = analyze_amplifiability(
            t["exons"], seqs[a], [seqs[o] for o in sib_accs], k=k,
            sibling_exons=[exons_by_acc[o] for o in sib_accs])
        sib_exons = {o: exons_by_acc[o] for o in sib_accs}
        coord_nu = bool(overlap.non_unique_partners(t["exons"], sib_exons))
        verdicts.append(_amp_to_verdict(a, t["is_mane"], t["exons"], seqs[a],
                                        t.get("cds"), amp, coord_nu))

    yield {"type": "progress", "pct": 94, "detail": "Designing Tm-guided primers…"}
    tgt = accs[target_acc]
    tgt_sib_accs = [o for o in accs if o != target_acc]
    tgt_sibs = {o: seqs[o] for o in tgt_sib_accs}
    tgt_amp = analyze_amplifiability(
        tgt["exons"], seqs[target_acc], [seqs[o] for o in tgt_sib_accs], k=k,
        sibling_exons=[exons_by_acc[o] for o in tgt_sib_accs])
    design = primers.design(tgt["exons"], seqs[target_acc], tgt_sibs, tgt_amp)
    tgt_verdict = next(v for v in verdicts if v.accession == target_acc)

    summary = GeneSummary(
        nm_count=len(accs),
        conventional_count=sum(v.tier == "CONVENTIONAL" for v in verdicts),
        needs_eej_count=sum(v.tier == "NEEDS_EEJ" for v in verdicts),
        hard_case_count=sum(v.tier == "NO_SINGLE_UNIQUE_JUNCTION" for v in verdicts),
        coord_non_unique_count=sum(v.coord_non_unique for v in verdicts),
    )
    verdicts.sort(key=lambda v: (v.accession != target_acc, not v.is_mane, v.accession))

    response = AnalyzeResponse(
        target_accession=target_acc,
        gene=GeneInfo(gene_id=gene_id, symbol=symbol, description=description, chromosome=chromosome),
        target_verdict=tgt_verdict,
        target_mrna=seqs[target_acc],
        primer_design=PrimerDesignOut(
            tier=design.tier, mechanism=design.mechanism,
            forward=_primer_out(design.forward), reverse=_primer_out(design.reverse),
            amplicon_len=design.amplicon_len, delta_tm=design.delta_tm,
            confidence=design.confidence, excluded_siblings=design.excluded_siblings,
            flags=design.flags, tm_method=design.tm_method, pair_dimer_tm=design.pair_dimer_tm,
        ),
        transcripts=verdicts,
        summary=summary,
        meta={"assembly": "GRCh38", "k": k},
    )
    yield {"type": "progress", "pct": 100, "detail": "Done"}
    yield {"type": "result", "result": response}
