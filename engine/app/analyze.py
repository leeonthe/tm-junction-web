"""Orchestration: accession -> gene -> siblings -> per-transcript verdict + primers.

Ties together ncbi (data), overlap (structural), amplify (sequence, authoritative),
and primers (Tm design) into the AnalyzeResponse contract.
"""

from __future__ import annotations

from . import ncbi, overlap, primers
from .amplify import AmplifyResult, analyze_amplifiability, cumulative_exon_ends
from .models import (
    AnalyzeResponse, Exon, GeneExonOut, GeneInfo, GeneLookupResponse, GeneSummary,
    GeneTranscriptOut, JunctionOut, PrimerDesignOut, PrimerOut, TranscriptVerdict,
    UniqueRegionOut,
)


class AnalysisError(Exception):
    def __init__(self, code: str, message: str):
        self.code = code
        self.message = message
        super().__init__(message)


def lookup_gene(symbol: str) -> GeneLookupResponse:
    """Human gene symbol -> its NM transcripts and exon alignment (no classification).

    A lightweight reference step: no mRNA sequences are fetched and no amplifiability is
    computed — just the structural picture so the user can choose a variant to analyze.
    """
    name = (symbol or "").strip()
    if not name:
        raise AnalysisError("BAD_REQUEST", "Enter a gene name, e.g. GAPDH.")
    try:
        report = ncbi.get_product_report(name)
        gene_id, sym, description, chromosome, strand, transcripts = ncbi.nm_transcripts(report)
    except ncbi.NotFound:
        raise AnalysisError("NOT_FOUND", f"No human gene found for “{name}”.")
    if not transcripts:
        raise AnalysisError("NOT_FOUND", f"“{name}” has no NM (mRNA) RefSeq transcripts on GRCh38.")

    out: list[GeneTranscriptOut] = []
    for t in transcripts:
        exons = t["exons"]                      # genomic (begin,end), ascending
        cds = t.get("cds")
        out.append(GeneTranscriptOut(
            accession=t["accession"],
            is_mane=t["is_mane"],
            exon_count=len(exons),
            length=sum(e - b + 1 for b, e in exons),
            cds_begin=cds[0] if cds else None,
            cds_end=cds[1] if cds else None,
            exons=[GeneExonOut(order=i + 1, begin=b, end=e) for i, (b, e) in enumerate(exons)],
        ))
    # MANE first, then longest transcript to shortest — a stable, useful reading order.
    out.sort(key=lambda x: (not x.is_mane, -x.length))
    gene = GeneInfo(gene_id=gene_id, symbol=sym, description=description,
                    chromosome=chromosome, strand=strand)
    return GeneLookupResponse(gene=gene, transcripts=out)


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


def _runs(positions: list[int]) -> list[list[int]]:
    """Sorted 0-based positions -> 1-based inclusive [lo, hi] runs of consecutive values."""
    out: list[list[int]] = []
    for p in sorted(positions):
        if out and p == out[-1][1]:      # out[-1][1] holds the 1-based end == next 0-based
            out[-1][1] = p + 1
        else:
            out.append([p + 1, p + 1])
    return out


def _uniq_genomic(exons: list[tuple[int, int]], tx_start: int, tx_end: int,
                  exon_order: int) -> tuple[int, int]:
    """Map a unique window span [tx_start, tx_end] (0-based mRNA) to its genomic (begin,end)
    within its exon, honoring strand (inferred from transcript vs genomic exon order)."""
    cum = cumulative_exon_ends(exons)                 # 0-based exclusive mRNA ends
    ei = exon_order - 1
    gb, ge = exons[ei]                                # genomic begin <= end
    exon_tx0 = cum[ei - 1] if ei else 0               # 0-based mRNA start of the exon
    off_s, off_e = tx_start - exon_tx0, tx_end - exon_tx0
    plus = len(exons) < 2 or exons[0][0] <= exons[-1][0]
    b, e = (gb + off_s, gb + off_e) if plus else (ge - off_e, ge - off_s)
    return (min(b, e), max(b, e))


def _region_tx(exons: list[tuple[int, int]], exon_order: int,
               gb: int, ge: int) -> tuple[int, int]:
    """Inverse of _uniq_genomic for a sub-span: map a genomic (gb,ge) inside exon `exon_order`
    to 1-based mRNA (tx_begin, tx_end), honoring strand."""
    cum = cumulative_exon_ends(exons)
    ei = exon_order - 1
    b, e = exons[ei]
    exon_tx1 = (cum[ei - 1] if ei else 0) + 1          # 1-based mRNA start of the exon
    plus = len(exons) < 2 or exons[0][0] <= exons[-1][0]
    if plus:
        return (exon_tx1 + (gb - b), exon_tx1 + (ge - b))
    return (exon_tx1 + (e - ge), exon_tx1 + (e - gb))


def _partner_exon(exons: list[tuple[int, int]], donor: int, acceptor: int) -> int | None:
    """Nearest exon to the exon-exon junction (donor|acceptor), for the partner primer.

    Candidates are the exon before the donor and the exon after the acceptor. The junction
    sits between the donor and acceptor exons; the traversal distance to the upstream candidate
    is the donor exon's length, and to the downstream candidate the acceptor exon's length —
    so pick the flanking exon that is shorter. Ties go upstream. Returns a 1-based exon order,
    or None when the junction is terminal (no exon on the chosen side).
    """
    n = len(exons)
    up = donor - 1        # exon order before the donor
    down = acceptor + 1   # exon order after the acceptor
    up_ok, down_ok = up >= 1, down <= n
    if not up_ok and not down_ok:
        return None
    if up_ok and down_ok:
        len_donor = exons[donor - 1][1] - exons[donor - 1][0]
        len_acceptor = exons[acceptor - 1][1] - exons[acceptor - 1][0]
        return up if len_donor <= len_acceptor else down
    return up if up_ok else down


def _amp_to_verdict(acc: str, is_mane: bool, exons: list[tuple[int, int]],
                    seq: str, cds: tuple[int, int] | None,
                    amp: AmplifyResult, coord_non_unique: bool) -> TranscriptVerdict:
    junctions = [
        JunctionOut(donor_order=j.donor_order, acceptor_order=j.acceptor_order,
                    label=_junction_label(j.donor_order, j.acceptor_order))
        for j in amp.unique_junctions
    ]
    # Recommended EEJ + orange exons depend on how the transcript is amplifiable:
    #   - single unique junction: recommended = that junction (red caret).
    #   - junction+exon combo:    recommended = the junction (red caret), orange = the exon.
    #   - two-junction combo:     no single recommended; the two junctions are magenta.
    #   - 7c conventional pair:   orange = the exon pair, no junction.
    recommended = junctions[0] if (amp.needs_eej and junctions) else None
    exon_pair_out = list(amp.exon_pair) if amp.exon_pair else None
    combo_junctions = None
    if amp.combo_je:
        d, a, e = amp.combo_je
        recommended = JunctionOut(donor_order=d, acceptor_order=a, label=_junction_label(d, a))
        exon_pair_out = [e]
    elif amp.combo_jj:
        (d1, a1), (d2, a2) = amp.combo_jj
        combo_junctions = [[d1, a1], [d2, a2]]
    # Original single-unique-junction EEJ: mark the nearest flanking exon as the target site for
    # the conventional partner primer. Of exon (donor-1) and exon (acceptor+1), pick whichever is
    # closer to the junction — i.e. the shorter of the two flanking exons the junction sits between.
    partner_exon = None
    if recommended is not None and exon_pair_out is None and combo_junctions is None:
        partner_exon = _partner_exon(exons, recommended.donor_order, recommended.acceptor_order)
    uniq_out = []
    for r in amp.unique_regions:
        # Report the DISCRIMINATING sequence (the exon minus what siblings carry), not the
        # k-mer placement envelope — see UniqueRegionOut. Fall back to the envelope only
        # when no sibling exon structures were available to subtract.
        if r.uniq_span:
            gb, ge = r.uniq_span
            tb, te = _region_tx(exons, r.exon_order, gb, ge)
        else:
            gb, ge = _uniq_genomic(exons, r.tx_start, r.tx_end, r.exon_order)
            tb, te = r.tx_start + 1, r.tx_end + 1
        uniq_out.append(UniqueRegionOut(
            exon_order=r.exon_order, window_count=r.window_count, side=r.side,
            begin=gb, end=ge, tx_begin=tb, tx_end=te, uniq_len=te - tb + 1,
            window_starts=_runs(amp.internal_starts.get(r.exon_order, []))))
    # Junction+exon combo: only the distinguishing part of the combo exon is the target site.
    # Emit it as a sub-span (genomic + mRNA) so the graph paints just that slice yellow and
    # leaves the overlapped remainder the tier color (red). exon_pair_out keeps the exon for the
    # verdict text; the sub-span makes the fill precise.
    if amp.combo_je and amp.combo_exon_region:
        _, _, ce = amp.combo_je
        gb, ge = amp.combo_exon_region
        tb, te = _region_tx(exons, ce, gb, ge)
        uniq_out.append(UniqueRegionOut(exon_order=ce, window_count=0, side="either",
                                        begin=gb, end=ge, tx_begin=tb, tx_end=te,
                                        uniq_len=te - tb + 1))
    return TranscriptVerdict(
        accession=acc,
        is_mane=is_mane,
        tier=amp.tier,
        amplifiable=amp.amplifiable,
        needs_eej=amp.needs_eej,
        unique_regions=uniq_out,
        unique_junctions=junctions,
        recommended_junction=recommended,
        coord_non_unique=coord_non_unique,
        exons=_exons_out(exons, seq, cds, amp),
        amplify_exon_pair=exon_pair_out,
        combo_junctions=combo_junctions,
        partner_exon=partner_exon,
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
        gene_id, symbol, description, chromosome, strand, transcripts = ncbi.nm_transcripts(report)
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
        gene=GeneInfo(gene_id=gene_id, symbol=symbol, description=description,
                      chromosome=chromosome, strand=strand),
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
