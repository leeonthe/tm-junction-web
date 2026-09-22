"""Orchestration: accession -> gene -> siblings -> per-transcript verdict + primers.

Ties together ncbi (data), overlap (structural), amplify (sequence, authoritative),
and primers (Tm design) into the AnalyzeResponse contract.
"""

from __future__ import annotations

from . import ncbi, overlap, panvariant, primers
from . import species as species_mod
from .amplify import AmplifyResult, analyze_amplifiability, cumulative_exon_ends
from .models import (
    FEATURES, AnalyzeResponse, Exon, GeneExonOut, GeneInfo, GeneLookupResponse, GeneSummary,
    GeneTranscriptOut, JunctionOut, PanVariantOut, PrimerDesignOut, PrimerOut,
    TranscriptVerdict, UniqueRegionOut,
)


class AnalysisError(Exception):
    def __init__(self, code: str, message: str):
        self.code = code
        self.message = message
        super().__init__(message)


def _gene_info(g: ncbi.GeneTranscripts) -> GeneInfo:
    sp = g.species
    return GeneInfo(gene_id=g.gene_id, symbol=g.symbol, description=g.description,
                    assembly=g.assembly, chromosome=g.chromosome, strand=g.strand,
                    species=sp.slug, organism=sp.scientific, common_name=sp.common,
                    tax_id=sp.tax_id)


def _no_transcripts(g: ncbi.GeneTranscripts) -> AnalysisError:
    """Why a gene NCBI knows has nothing to analyze — the two reasons are different facts.

    Either NCBI has curated no transcript for it (only XM/XR models, or none), or it has
    and has not placed them on the reference assembly yet: a record revised since the last
    annotation run is listed bare until the next one. Rat Gapdh is the second kind, and
    "no transcripts" would be false of it — there is one, NM_017008.5, without coordinates.
    """
    who = f"{g.species.common} {g.symbol}"
    if g.unplaced:
        n = len(g.unplaced)
        listed = ", ".join(g.unplaced[:4]) + (f" and {n - 4} more" if n > 4 else "")
        return AnalysisError(
            "NOT_PLACED",
            f"NCBI lists {n} curated transcript{'s' if n != 1 else ''} for {who} ({listed}) but "
            f"has not placed {'them' if n != 1 else 'it'} on the {g.species.scientific} reference "
            "assembly — usually a record revised since NCBI's last annotation run. Without "
            "genomic coordinates there is no exon structure to design against; NCBI adds "
            "them at its next annotation release for this species.")
    return AnalysisError(
        "NOT_FOUND",
        f"{who} has no curated RefSeq transcripts (NM mRNA or NR non-coding RNA) — NCBI "
        "has only predicted models (XM/XR) for it, or none.")


def lookup_gene(symbol: str, species: str | None = None) -> GeneLookupResponse:
    """Gene symbol (+ species, human by default) -> its NM/NR transcripts and exon
    alignment (no classification).

    A lightweight reference step: no mRNA sequences are fetched and no amplifiability is
    computed — just the structural picture so the user can choose a variant to analyze.
    """
    name = (symbol or "").strip()
    if not name:
        raise AnalysisError("BAD_REQUEST", "Enter a gene name, e.g. GAPDH.")
    try:
        sp = species_mod.get(species)
    except species_mod.UnknownSpecies as e:
        raise AnalysisError("BAD_SPECIES", str(e))
    try:
        gene = ncbi.refseq_transcripts(ncbi.get_product_report(name, sp), sp)
    except ncbi.NotFound:
        raise AnalysisError("NOT_FOUND", f"No {sp.common.lower()} ({sp.scientific}) gene found "
                                         f"for “{name}”.")
    transcripts = gene.transcripts
    if not transcripts:
        raise _no_transcripts(gene)

    # Every isoform gets its designation, from the product report or from the record's own
    # title — the picker names variants as much as the verdict table does.
    ncbi.fill_variants(transcripts)
    # One entry per distinct molecule — several accessions for one exon structure are one
    # transcript, and listing them separately invites picking a "different" variant that is
    # the same sequence. See _same_structure_groups.
    reps, same_struct = _same_structure_groups(transcripts)
    by_acc = {t["accession"]: t for t in transcripts}
    out: list[GeneTranscriptOut] = []
    for t in reps:
        exons = t["exons"]                      # genomic (begin,end), ascending
        cds = t.get("cds")
        folded = same_struct.get(t["accession"], [])
        out.append(GeneTranscriptOut(
            accession=t["accession"],
            variant=t.get("variant"),
            same_structure_accessions=folded,
            same_structure_variants=[by_acc[n].get("variant") for n in folded],
            placed_via=t.get("placed_via"),
            is_mane=t["is_mane"],
            exon_count=len(exons),
            length=sum(e - b + 1 for b, e in exons),
            cds_begin=cds[0] if cds else None,
            cds_end=cds[1] if cds else None,
            exons=[GeneExonOut(order=i + 1, begin=b, end=e) for i, (b, e) in enumerate(exons)],
        ))
    # MANE first, then longest transcript to shortest — a stable, useful reading order.
    out.sort(key=lambda x: (not x.is_mane, -x.length))
    return GeneLookupResponse(gene=_gene_info(gene), transcripts=out)


def _pick_representative(members: list[str], accs: dict, target_acc: str | None) -> str:
    """Which accession speaks for a set of accessions that are the same molecule.

    The one the user asked about, if it is in the set — the analysis is about their
    accession, not a synonym of it. Otherwise MANE, then the lowest accession, so the
    choice is stable between runs.
    """
    if target_acc in members:
        return target_acc
    mane = sorted(m for m in members if accs[m]["is_mane"])
    return mane[0] if mane else sorted(members)[0]


def _same_sequence_groups(accs: dict, seqs: dict[str, str],
                          target_acc: str | None = None) -> tuple[list[str], dict[str, list[str]]]:
    """Fold accessions that are the same molecule into one transcript each.

    Same molecule means: byte-identical mRNA AND the identical exon structure at the
    identical GRCh38 coordinates. Both are checked because both are what the claim rests
    on — the accession number is not evidence of anything. (Across every gene cached here
    the two criteria never disagree; requiring both keeps a coincidence of sequence at a
    different locus from being folded away silently.)

    RefSeq issues several accessions for one molecule — TP53 alone has 25 NM accessions
    for 13 distinct sequences, e.g. NM_001126115.2 and NM_001276697.3, both "transcript
    variant 5". Treating those as separate isoforms is not a cosmetic duplication: each is
    then compared against an identical sibling, no primer can tell them apart, and the tool
    calls a perfectly designable transcript EEJ-infeasible.

    Returns (representatives in the input order, {representative: other accessions}).
    """
    by_seq: dict[tuple, list[str]] = {}
    for a in accs:
        key = (seqs[a], tuple(map(tuple, accs[a]["exons"])))
        by_seq.setdefault(key, []).append(a)
    order = {a: i for i, a in enumerate(accs)}
    reps: list[str] = []
    same: dict[str, list[str]] = {}
    for members in by_seq.values():
        rep = _pick_representative(members, accs, target_acc)
        reps.append(rep)
        same[rep] = sorted(m for m in members if m != rep)
    reps.sort(key=lambda a: order[a])
    return reps, same


def _same_structure_groups(transcripts: list[dict]) -> tuple[list[dict], dict[str, list[str]]]:
    """The same fold for the gene picker, which has exon structures but no sequences.

    Identical exon structures on one genome splice to the identical mRNA, so this reaches
    the same grouping as _same_sequence_groups without paying for a sequence fetch per
    transcript. Kept separate because the evidence differs: a RefSeq sequence may carry an
    annotated difference from the genome, which only the sequences themselves would show.
    """
    by_struct: dict[tuple, list[dict]] = {}
    for t in transcripts:
        by_struct.setdefault(tuple(map(tuple, t["exons"])), []).append(t)
    accs = {t["accession"]: t for t in transcripts}
    reps: list[dict] = []
    same: dict[str, list[str]] = {}
    for members in by_struct.values():
        names = [m["accession"] for m in members]
        rep = _pick_representative(names, accs, None)
        reps.append(accs[rep])
        same[rep] = sorted(n for n in names if n != rep)
    return reps, same


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


def _is_plus(exons: list[tuple[int, int]], strand: str = "") -> bool:
    """Is the transcript read left-to-right on the genome?

    Exon order says so for any transcript with two exons. A single-exon transcript HAS no
    exon order, and assuming "plus" for it mirrored every position inside the exon on a
    minus-strand gene: fly Gapdh1's variant A is unique over mRNA 1–150 and was reported as
    unique over 1331–1480, the other end. So one exon takes NCBI's stated strand.
    """
    if len(exons) >= 2:
        return exons[0][0] <= exons[-1][0]
    return strand != "-"


def _uniq_genomic(exons: list[tuple[int, int]], tx_start: int, tx_end: int,
                  exon_order: int, strand: str = "") -> tuple[int, int]:
    """Map a unique window span [tx_start, tx_end] (0-based mRNA) to its genomic (begin,end)
    within its exon, honoring strand (inferred from transcript vs genomic exon order)."""
    cum = cumulative_exon_ends(exons)                 # 0-based exclusive mRNA ends
    ei = exon_order - 1
    gb, ge = exons[ei]                                # genomic begin <= end
    exon_tx0 = cum[ei - 1] if ei else 0               # 0-based mRNA start of the exon
    off_s, off_e = tx_start - exon_tx0, tx_end - exon_tx0
    plus = _is_plus(exons, strand)
    b, e = (gb + off_s, gb + off_e) if plus else (ge - off_e, ge - off_s)
    return (min(b, e), max(b, e))


def _region_tx(exons: list[tuple[int, int]], exon_order: int,
               gb: int, ge: int, strand: str = "") -> tuple[int, int]:
    """Inverse of _uniq_genomic for a sub-span: map a genomic (gb,ge) inside exon `exon_order`
    to 1-based mRNA (tx_begin, tx_end), honoring strand."""
    cum = cumulative_exon_ends(exons)
    ei = exon_order - 1
    b, e = exons[ei]
    exon_tx1 = (cum[ei - 1] if ei else 0) + 1          # 1-based mRNA start of the exon
    plus = _is_plus(exons, strand)
    if plus:
        return (exon_tx1 + (gb - b), exon_tx1 + (ge - b))
    return (exon_tx1 + (e - ge), exon_tx1 + (e - gb))


def _amp_to_verdict(acc: str, is_mane: bool, exons: list[tuple[int, int]],
                    seq: str, cds: tuple[int, int] | None,
                    amp: AmplifyResult, coord_non_unique: bool,
                    same_sequence: list[str] | None = None,
                    variant: str | None = None,
                    same_variants: list[str | None] | None = None,
                    placed_via: str | None = None, strand: str = "") -> TranscriptVerdict:
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
    uniq_out = []
    for r in amp.unique_regions:
        # Report the DISCRIMINATING sequence (the exon minus what siblings carry), not the
        # k-mer placement envelope — see UniqueRegionOut. Fall back to the envelope only
        # when no sibling exon structures were available to subtract.
        if r.uniq_span:
            gb, ge = r.uniq_span
            tb, te = _region_tx(exons, r.exon_order, gb, ge, strand)
        else:
            gb, ge = _uniq_genomic(exons, r.tx_start, r.tx_end, r.exon_order, strand)
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
        tb, te = _region_tx(exons, ce, gb, ge, strand)
        uniq_out.append(UniqueRegionOut(exon_order=ce, window_count=0, side="either",
                                        begin=gb, end=ge, tx_begin=tb, tx_end=te,
                                        uniq_len=te - tb + 1))
    return TranscriptVerdict(
        accession=acc,
        variant=variant,
        same_sequence_accessions=list(same_sequence or []),
        same_sequence_variants=list(same_variants or []),
        placed_via=placed_via,
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
    )


def _pan_out(pan) -> PanVariantOut:
    return PanVariantOut(
        forward=_primer_out(pan.forward), reverse=_primer_out(pan.reverse),
        amplicon_len=pan.amplicon_len, covered=pan.covered, uncovered=pan.uncovered,
        reference=pan.reference, exons=pan.exons, flags=pan.flags,
    )


def _primer_out(p) -> PrimerOut | None:
    if p is None:
        return None
    return PrimerOut(seq=p.seq, tm=p.tm, gc=p.gc, length=p.length, kind=p.kind,
                     role=p.role, anchor=p.anchor, hairpin_tm=p.hairpin_tm,
                     homodimer_tm=p.homodimer_tm, qc_pass=p.qc_pass, tx_start=p.tx_start)


def _match_version(accession: str, accs) -> str | None:
    """The gene's accession an input names: exactly, or by its base when the version given
    is absent or not the current one."""
    if accession in accs:
        return accession
    base = accession.split(".")[0]
    return next((a for a in accs if a.split(".")[0] == base), None)


def _gene_of(accession: str) -> tuple[ncbi.GeneTranscripts, str]:
    """An accession -> (its gene's curated transcripts, the accession as the gene lists it).

    The species is never asked for: the accession resolves to its own gene, and the gene's
    record says whose it is. A species the engine does not cover is an error that says so.

    An entry resolved before species support carries no tax id and is read as human. If
    the gene that yields does not contain the accession, the entry may be a non-human
    record from back when that could not be told apart — so it is resolved afresh, once.
    """
    symbol, _gid, tax_id = ncbi.resolve_accession(accession)
    for attempt in (0, 1):
        sp = species_mod.by_tax_id(tax_id)
        if sp is None:
            raise AnalysisError(
                "UNSUPPORTED_SPECIES",
                f"{accession} belongs to a species this tool does not cover (NCBI taxonomy "
                f"id {tax_id}). Supported: {species_mod.supported_names()}.")
        gene = ncbi.refseq_transcripts(ncbi.get_product_report(symbol, sp), sp)
        target = _match_version(accession, [t["accession"] for t in gene.transcripts])
        if target or tax_id or attempt:
            return gene, target or accession
        try:
            symbol, _gid, tax_id = ncbi.resolve_accession(accession, refresh=True)
        except ncbi.RateLimited:
            raise
        except Exception:
            return gene, accession          # offline or unknown: the first answer stands
    raise AssertionError("unreachable")


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
    The sequence fetches (one per isoform) are the bulk of the wall-clock on a cold cache,
    so progress is driven mostly by those. Raises AnalysisError for bad input / not-found.
    """
    accession = accession.strip().upper()   # RefSeq accessions are uppercase
    if not ncbi.is_valid_refseq(accession):
        # The code keeps its NOT_NM name — clients match on it — though NR is now valid too.
        raise AnalysisError("NOT_NM", "Enter a curated RefSeq transcript accession — NM (mRNA) "
                                      "or NR (non-coding RNA), e.g. NM_002046.7 or NR_152150.2.")

    yield {"type": "progress", "pct": 3, "detail": "Resolving gene…"}
    try:
        gene, target_acc = _gene_of(accession)
    except ncbi.NotFound as e:
        raise AnalysisError("NOT_FOUND", str(e))
    symbol = gene.symbol

    # NM and NR together: a gene's non-coding transcripts are in the same cDNA as its
    # mRNAs, so they are siblings to be told apart from (or co-amplified with) like any
    # other — GAPDH's NR_152150 can take a primer pair as readily as NM_002046 can.
    accs = {t["accession"]: t for t in gene.transcripts}
    if target_acc not in accs:
        if accession.split(".")[0] in {a.split(".")[0] for a in gene.unplaced}:
            raise _no_transcripts(gene)
        raise AnalysisError("NOT_FOUND", f"{accession} is not a curated RefSeq transcript of "
                                         f"{gene.species.common.lower()} {symbol} on "
                                         f"{gene.assembly or 'the reference assembly'}.")

    n = len(accs)
    n_nr = sum(ncbi.is_noncoding(a) for a in accs)
    kinds = f" ({n - n_nr} NM + {n_nr} NR)" if n_nr and n_nr < n else " (NR)" if n_nr else ""
    who = symbol if gene.species is species_mod.HUMAN else f"{gene.species.common} {symbol}"
    yield {"type": "progress", "pct": 10,
           "detail": f"{who}: {n} isoform{'s' if n != 1 else ''}{kinds}"}

    # Fetch sequences — the main cost. Cache-first, and every miss travels in ONE efetch
    # rather than a call per isoform: the per-transcript loop this replaces is what turned
    # a big uncached gene into a burst of NCBI requests and, on a shared egress IP, 429s.
    yield {"type": "progress", "pct": 30, "detail": f"Fetching transcript sequences ({n})…"}
    seqs = ncbi.get_sequences(list(accs))
    yield {"type": "progress", "pct": 78, "detail": "Sequences fetched"}

    exons_by_acc = {a: t["exons"] for a, t in accs.items()}

    # A gene with more than one transcript has isoforms to tell apart, so each one has a variant
    # designation; where the product report omitted it, take it from the record's own title.
    ncbi.fill_variants(list(accs.values()))

    # One transcript per SEQUENCE, not per accession: a byte-identical twin is the same
    # molecule, and comparing a transcript against it would find nothing that tells the two
    # apart — see _same_sequence_groups. The target always represents its own group.
    reps, same_seq = _same_sequence_groups(accs, seqs, target_acc)
    if len(reps) < len(accs):
        merged = len(accs) - len(reps)
        yield {"type": "progress", "pct": 80,
               "detail": f"{merged} accession{'s' if merged != 1 else ''} share a sequence "
                         f"with another — {len(reps)} distinct isoforms"}

    yield {"type": "progress", "pct": 82, "detail": "Classifying isoforms…"}
    verdicts: list[TranscriptVerdict] = []
    for a in reps:
        t = accs[a]
        sib_accs = [o for o in reps if o != a]
        amp = analyze_amplifiability(
            t["exons"], seqs[a], [seqs[o] for o in sib_accs], k=k,
            sibling_exons=[exons_by_acc[o] for o in sib_accs])
        sib_exons = {o: exons_by_acc[o] for o in sib_accs}
        coord_nu = bool(overlap.non_unique_partners(t["exons"], sib_exons))
        verdicts.append(_amp_to_verdict(a, t["is_mane"], t["exons"], seqs[a],
                                        t.get("cds"), amp, coord_nu,
                                        same_seq.get(a, []), t.get("variant"),
                                        [accs[x].get("variant") for x in same_seq.get(a, [])],
                                        t.get("placed_via"), t.get("strand") or gene.strand))

    yield {"type": "progress", "pct": 94, "detail": "Designing Tm-guided primers…"}
    tgt = accs[target_acc]
    tgt_sib_accs = [o for o in reps if o != target_acc]
    tgt_sibs = {o: seqs[o] for o in tgt_sib_accs}
    tgt_amp = analyze_amplifiability(
        tgt["exons"], seqs[target_acc], [seqs[o] for o in tgt_sib_accs], k=k,
        sibling_exons=[exons_by_acc[o] for o in tgt_sib_accs])
    design = primers.design(tgt["exons"], seqs[target_acc], tgt_sibs, tgt_amp)
    tgt_verdict = next(v for v in verdicts if v.accession == target_acc)

    # The other question a user can ask of a gene: one pair for ALL of its variants. Design
    # it over the folded set, so coverage counts distinct transcripts rather than accessions.
    yield {"type": "progress", "pct": 97, "detail": "Designing whole-transcript pairs…"}
    pan_opts = panvariant.design_options({a: accs[a] for a in reps}, seqs, order=reps)

    summary = GeneSummary(
        nm_count=len(reps),
        nr_count=sum(ncbi.is_noncoding(a) for a in reps),
        merged_accession_count=len(accs) - len(reps),
        conventional_count=sum(v.tier == "CONVENTIONAL" for v in verdicts),
        needs_eej_count=sum(v.tier == "NEEDS_EEJ" for v in verdicts),
        hard_case_count=sum(v.tier == "NO_SINGLE_UNIQUE_JUNCTION" for v in verdicts),
        coord_non_unique_count=sum(v.coord_non_unique for v in verdicts),
    )
    verdicts.sort(key=lambda v: (v.accession != target_acc, not v.is_mane, v.accession))

    response = AnalyzeResponse(
        target_accession=target_acc,
        gene=_gene_info(gene),
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
        pan_variant=_pan_out(pan_opts[0]) if pan_opts else None,
        pan_variant_options=[_pan_out(o) for o in pan_opts],
        meta={"assembly": gene.assembly, "species": gene.species.slug, "k": k,
              "features": FEATURES},
    )
    yield {"type": "progress", "pct": 100, "detail": "Done"}
    yield {"type": "result", "result": response}
