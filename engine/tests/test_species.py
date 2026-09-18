"""Other species: mouse, rat, fruit fly, baker's yeast, zebrafish — as well as human.

Human was pinned in four places: the Datasets URL (`/taxon/human/`), the "NC_0000" prefix
that picked the GRCh38 placement, `Homo sapiens[Organism]` / `txid9606` in the two
typeaheads, and a human-only accession index. Each is now a property of a `Species` or is
read from NCBI's own data. What these tests pin is the evidence gathered before writing a
line of it — NCBI does NOT behave the same for every species:

  * Saccharomyces cerevisiae's genes are filed under the S288C strain (559292); a lookup
    under the species name or id (4932) finds no gene at all.
  * Symbol lookup is case-insensitive and answers with the official spelling, which is
    all-capitals in human and yeast, first-letter in mouse and rat, all lower case in
    zebrafish and mixed in fly. So capitalization is NCBI's answer, never the engine's rule.
  * Rat, zebrafish and human are annotated on a second, "Alternate" assembly as well.
    Fly and yeast placements carry no name to tell a reference by.
  * A record NCBI has revised since its last annotation run is listed with no placement
    (rat Gapdh, mouse Trp53) — rare (0-3% of genes sampled per species), and not "no
    transcripts".

All offline, from fixtures in data/cache.
"""

import glob
import json
import re

import pytest
from fastapi.testclient import TestClient

from app import ncbi
from app import species as sm
from app.analyze import AnalysisError, analyze, lookup_gene
from app.main import app
from app.models import FEATURES

client = TestClient(app)


# ---------------------------------------------------------------- the registry

def test_every_query_is_scoped_by_taxonomy_id():
    assert [s.slug for s in sm.SPECIES] == ["human", "mouse", "rat", "fly", "yeast", "zebrafish"]
    assert {s.slug: s.tax_id for s in sm.SPECIES} == {
        "human": "9606", "mouse": "10090", "rat": "10116", "fly": "7227",
        "yeast": "559292",                     # S288C — NOT 4932, where NCBI has no genes
        "zebrafish": "7955"}


def test_a_species_can_be_named_any_reasonable_way():
    mouse = sm.get("mouse")
    for name in ("Mouse", "MOUSE", "10090", "Mus musculus", "mus musculus", " mouse "):
        assert sm.get(name) is mouse
    assert sm.get(None) is sm.HUMAN and sm.get("") is sm.HUMAN      # every pre-species caller
    yeast = sm.get("yeast")
    for name in ("4932", "559292", "Saccharomyces cerevisiae", "Baker's yeast"):
        assert sm.get(name) is yeast
    with pytest.raises(sm.UnknownSpecies) as e:
        sm.get("unicorn")
    assert "Zebrafish" in str(e.value)                              # says what IS supported


def test_a_record_names_its_species_by_tax_id():
    assert sm.by_tax_id("10090").slug == "mouse"
    assert sm.by_tax_id(7955).slug == "zebrafish"
    assert sm.by_tax_id("") is sm.HUMAN          # cached before species support: human
    assert sm.by_tax_id("9598") is None          # chimpanzee — a real record, not covered


# ---------------------------------------------------------------- which placement

def _loc(acc, name=None, exons=((100, 200),), orientation="plus"):
    loc = {"genomic_accession_version": acc,
           "genomic_range": {"orientation": orientation},
           "exons": [{"begin": str(b), "end": str(e)} for b, e in exons]}
    if name:
        loc["sequence_name"] = name
    return loc


def test_the_reference_chromosome_wins_over_every_other_placement():
    ref = _loc("NC_086030.1", "Chromosome 12 Reference GRCr8")
    # Listed LAST, behind an alternate assembly, a patch and an alternate locus.
    t = {"genomic_locations": [
        _loc("NC_141545.1", "Chromosome 12 Alternate SHRSP_T2T Primary Assembly"),
        _loc("NW_009646202.1", "Chromosome 12 Reference GRCr8 PATCHES"),
        _loc("NT_167244.2", "Chromosome 12 Reference GRCr8 ALT_REF_LOCI_1"),
        ref]}
    assert ncbi.reference_location(t) is ref
    assert ncbi._reference_labels(t, "0") == ("12", "GRCr8")


def test_a_transcript_placed_only_off_the_reference_is_not_placed():
    for only in (_loc("NC_060936.1", "Chromosome 12 Alternate T2T-CHM13v2.0"),
                 _loc("NT_167244.2", "Chromosome 6 Reference GRCh38.p14 ALT_REF_LOCI_1"),
                 _loc("NW_009646202.1", "Chromosome 1 Reference GRCh38.p14 PATCHES"),
                 _loc("NT_187383.1", "NT_187383 Chromosome 16 Unlocalized Scaffold Reference "
                                     "GRCh38.p14 Primary Assembly"),
                 _loc("NG_000017.2", "genomic region")):
        assert ncbi.reference_location({"genomic_locations": [only]}) is None


def test_an_assembly_that_names_nothing_has_one_placement():
    """Fly and yeast: a single assembly, and the product report names no sequence."""
    arm = _loc("NT_033778.4", orientation="minus")        # 2R — an NT_ record, and the chromosome
    t = {"genomic_locations": [arm]}
    assert ncbi.reference_location(t) is arm
    assert ncbi.reference_strand(t) == "-"
    # An empty placement is no placement (18 human records carry `[{}]`).
    assert ncbi.reference_location({"genomic_locations": [{}]}) is None
    assert ncbi.reference_exons({"genomic_locations": [{}]}) is None


def test_a_patch_release_is_not_a_different_assembly():
    assert ncbi._assembly_label("GRCh38.p14") == "GRCh38"
    assert ncbi._assembly_label("GRCm39") == "GRCm39"
    assert ncbi._assembly_label("Release 6 plus ISO1 MT") == "Release 6 plus ISO1 MT"


def test_on_human_the_rule_selects_exactly_what_the_NC_0000_prefix_did():
    """The regression that licenses the change: for every curated transcript of every human
    fixture, the species-neutral rule picks the placement the hardcoded prefix picked, and
    reads the same chromosome off its name that the accession number used to encode."""
    def old_chromosome(acc):
        n = int(re.match(r"NC_0*(\d+)", acc).group(1))
        return str(n) if 1 <= n <= 22 else {23: "X", 24: "Y"}[n]

    checked = 0
    for path in glob.glob(str(ncbi.CACHE_DIR / "product_report" / "*__human.json")):
        reports = json.load(open(path)).get("reports") or []
        product = reports[0]["product"] if reports else {}
        for t in product.get("transcripts") or []:
            if not (t.get("accession_version") or "").startswith(ncbi.REFSEQ_PREFIXES):
                continue
            old = next((loc for loc in t.get("genomic_locations") or []
                        if str(loc.get("genomic_accession_version", "")).startswith("NC_0000")), None)
            new = ncbi.reference_location(t)
            if not (old and old.get("exons")):
                assert new is None, t["accession_version"]
                continue
            assert new is old, t["accession_version"]
            acc = old["genomic_accession_version"]
            assert ncbi._reference_labels(t, product["gene_id"]) == (old_chromosome(acc), "GRCh38")
            checked += 1
    assert checked > 500


# ---------------------------------------------------------------- gene lookup

def test_mouse_gapdh_the_ticketed_search():
    g = lookup_gene("Gapdh", "mouse")
    assert (g.gene.symbol, g.gene.gene_id) == ("Gapdh", "14433")
    assert (g.gene.species, g.gene.organism, g.gene.common_name, g.gene.tax_id) == (
        "mouse", "Mus musculus", "Mouse", "10090")
    assert (g.gene.assembly, g.gene.chromosome, g.gene.strand) == ("GRCm39", "6", "-")
    assert len(g.transcripts) == 8
    assert all(t.accession.startswith("NM_") for t in g.transcripts)
    assert not any(t.is_mane for t in g.transcripts)      # MANE is a human designation


def test_capitalization_is_ncbis_answer_not_a_rule():
    """Typed any way, the gene is the same and comes back in its official spelling."""
    for typed in ("Gapdh", "gapdh", "GAPDH"):
        assert lookup_gene(typed, "mouse").gene.symbol == "Gapdh"      # first letter only
    assert lookup_gene("GAPDH", "zebrafish").gene.symbol == "gapdh"    # all lower case
    assert lookup_gene("tdh3", "yeast").gene.symbol == "TDH3"          # all capitals
    assert lookup_gene("DPP", "fly").gene.symbol == "dpp"
    assert lookup_gene("gapdh").gene.symbol == "GAPDH"                 # and human, as ever


def test_the_same_symbol_is_a_different_gene_in_each_species():
    ids = {sp: lookup_gene("gapdh", sp).gene.gene_id for sp in ("human", "mouse", "zebrafish")}
    assert ids == {"human": "2597", "mouse": "14433", "zebrafish": "317743"}


def test_rat_ignores_its_second_assembly():
    g = lookup_gene("Tp53", "rat")
    assert (g.gene.assembly, g.gene.chromosome, g.gene.strand) == ("GRCr8", "10", "+")
    assert len(g.transcripts) == 5
    # Every exon is on GRCr8's chromosome 10 (~54.5 Mb), none on SHRSP_T2T's.
    report = ncbi.get_product_report("Tp53", sm.get("rat"))
    for t in report["reports"][0]["product"]["transcripts"]:
        if t["accession_version"].startswith("NM_"):
            assert len(t["genomic_locations"]) == 2                        # both are offered
            assert ncbi.reference_location(t)["genomic_accession_version"] == "NC_086028.1"


def test_zebrafish_reads_its_reference_not_the_alternate():
    g = lookup_gene("tp53", "zebrafish")
    assert (g.gene.assembly, g.gene.chromosome) == ("GRCz12ab", "5")
    assert [t.variant for t in g.transcripts if t.variant] != []       # named from the records


def test_fly_names_its_chromosome_arm_and_assembly_from_the_gene_record():
    g = lookup_gene("dpp", "fly")
    assert (g.gene.assembly, g.gene.chromosome, g.gene.strand) == ("Release 6 plus ISO1 MT", "2L", "+")
    assert sorted(t.variant for t in g.transcripts) == [
        "transcript variant A", "transcript variant B", "transcript variant C", "transcript variant E"]


def test_fly_titles_can_put_the_symbol_after_the_variant():
    assert ncbi.variant_in("Drosophila melanogaster glyceraldehyde 3 phosphate dehydrogenase 1, "
                           "transcript variant A (Gapdh1), mRNA") == "transcript variant A"
    assert ncbi.variant_in("Drosophila melanogaster decapentaplegic (dpp), transcript variant C, "
                           "mRNA") == "transcript variant C"
    assert {t.variant for t in lookup_gene("Gapdh1", "fly").transcripts} == {
        "transcript variant A", "transcript variant B"}


def test_yeast_is_found_under_its_reference_strain():
    g = lookup_gene("ACT1", "yeast")
    assert (g.gene.symbol, g.gene.tax_id, g.gene.organism) == ("ACT1", "559292", "Saccharomyces cerevisiae")
    assert (g.gene.assembly, g.gene.chromosome, g.gene.strand) == ("R64", "VI", "-")
    assert [(t.exon_count, t.length) for t in g.transcripts] == [(2, 1128)]   # its one intron
    assert lookup_gene("TDH3", "yeast").transcripts[0].exon_count == 1


def test_an_unplaced_record_is_not_a_gene_without_transcripts():
    """Rat Gapdh: NM_017008.5 replaced .4 on Sep 1, 2026, after the annotation run that
    would have placed it, so NCBI lists it bare. The answer names it and says why."""
    with pytest.raises(AnalysisError) as e:
        lookup_gene("Gapdh", "rat")
    assert e.value.code == "NOT_PLACED"
    assert "NM_017008.5" in e.value.message and "Rat Gapdh" in e.value.message
    # Reached by accession, the same fact — not "not a transcript of Gapdh".
    with pytest.raises(AnalysisError) as e:
        analyze("NM_017008.5")
    assert e.value.code == "NOT_PLACED"
    # A gene NCBI has curated nothing for says THAT instead.
    with pytest.raises(AnalysisError) as e:
        lookup_gene("mt:CoI", "fly")
    assert e.value.code == "NOT_FOUND" and "no curated RefSeq transcripts" in e.value.message


def test_errors_name_the_species():
    with pytest.raises(AnalysisError) as e:
        lookup_gene("Gapdh", "unicorn")
    assert e.value.code == "BAD_SPECIES"


# ---------------------------------------------------------------- analysis by accession

def test_an_accession_names_its_own_species():
    """No species is asked for: NM_008084 is a mouse transcript because NCBI says so."""
    r = analyze("nm_008084")                       # any case, no version — as for human
    assert r.target_accession == "NM_008084.4"
    assert (r.gene.symbol, r.gene.species, r.gene.assembly) == ("Gapdh", "mouse", "GRCm39")
    assert r.meta["species"] == "mouse" and r.meta["assembly"] == "GRCm39"
    assert r.summary.nm_count == 8
    assert {t.accession for t in r.transcripts} == {t.accession for t in lookup_gene("Gapdh", "mouse").transcripts}
    # One pair for all eight, verified from their sequences like any human gene's.
    assert len(r.pan_variant_options[0].covered) == 8


def test_primers_designed_for_another_species_are_specific_within_it():
    r = analyze("NM_001411845.1")                  # mouse Gapdh variant 8: unique first exon
    assert r.target_verdict.tier == "CONVENTIONAL"
    fwd, rev = r.primer_design.forward, r.primer_design.reverse
    assert fwd is not None and rev is not None
    sibs = ncbi.get_sequences([t.accession for t in r.transcripts if t.accession != r.target_accession])
    assert fwd.seq in r.target_mrna
    assert all(fwd.seq not in s for s in sibs.values())


@pytest.mark.parametrize("accession, species, symbol, assembly, isoforms", [
    ("NM_001429993.1", "rat", "Tp53", "GRCr8", 5),
    ("NM_057963.5", "fly", "dpp", "Release 6 plus ISO1 MT", 4),
    ("NM_131327.2", "zebrafish", "tp53", "GRCz12ab", 4),
    ("NM_001179927.1", "yeast", "ACT1", "R64", 1),
])
def test_every_species_analyzes_end_to_end(accession, species, symbol, assembly, isoforms):
    r = analyze(accession)
    assert (r.gene.species, r.gene.symbol, r.gene.assembly) == (species, symbol, assembly)
    assert r.summary.nm_count == isoforms == len(r.transcripts)
    # The exon structure fits the record's own sequence, whatever assembly it came from:
    # exon for exon, up to the poly(A) tail a record may carry and no genome contains.
    placed = sum(e.length for e in r.target_verdict.exons)
    assert placed <= len(r.target_mrna) and set(r.target_mrna[placed:]) <= {"A"}
    assert r.target_verdict.amplifiable


def test_a_species_the_engine_does_not_cover_says_so(monkeypatch):
    monkeypatch.setattr(ncbi, "resolve_accession", lambda acc, refresh=False: ("GAPDH", "449484", "9598"))
    with pytest.raises(AnalysisError) as e:
        analyze("NM_001009133.1")
    assert e.value.code == "UNSUPPORTED_SPECIES"
    assert "9598" in e.value.message and "Mouse" in e.value.message


def test_an_entry_resolved_before_species_support_heals_itself(monkeypatch):
    """Cached without a tax id, an accession is read as human. A mouse record resolved back
    then yields human GAPDH (the symbol matches case-insensitively), which does not contain
    it — so it is resolved afresh, once, rather than failing forever on a stale cache."""
    calls = []

    def resolve(acc, refresh=False):
        calls.append(refresh)
        return ("Gapdh", "14433", "10090" if refresh else "")

    monkeypatch.setattr(ncbi, "resolve_accession", resolve)
    r = analyze("NM_008084.4")
    assert r.gene.species == "mouse" and calls == [False, True]


# ---------------------------------------------------------------- files and URLs

def test_a_symbol_never_reaches_a_path_or_a_url_unescaped():
    fly = sm.get("fly")
    assert ncbi._product_report_path("GAPDH", sm.HUMAN).name == "GAPDH__human.json"   # unchanged
    assert ncbi._product_report_path("l(2)gl", fly).name == "l%282%29gl__fly.json"
    assert ncbi._product_report_path("mt:CoI", fly).name == "mt%3ACoI__fly.json"
    hostile = ncbi._product_report_path("../../etc/passwd", fly)
    assert hostile.parent == ncbi.CACHE_DIR / "product_report" and "/" not in hostile.name


# ---------------------------------------------------------------- typeahead

def test_each_species_has_an_accession_index():
    for sp in sm.SPECIES:
        accessions, genes = ncbi._load_index(sp)
        assert len(accessions) == len(genes) > 5000, sp.slug
        assert accessions == sorted(accessions) and len(set(accessions)) == len(accessions)
        assert all(a.startswith(ncbi.REFSEQ_PREFIXES) for a in accessions[:50] + accessions[-50:])


def test_the_accession_box_searches_every_species_at_once():
    """An accession belongs to one organism, so the box needs no species — each row says."""
    for prefix, want in [("NM_008084", ("NM_008084.4", "Gapdh", "mouse")),
                         ("NM_031144", ("NM_031144.3", "Actb", "rat")),
                         ("NM_057963", ("NM_057963.5", "dpp", "fly")),
                         ("NM_001181321", ("NM_001181321.3", "TDH3", "yeast")),
                         ("NM_131327", ("NM_131327.2", "tp53", "zebrafish")),
                         ("NM_002046", ("NM_002046.7", "GAPDH", "human"))]:
        rows = ncbi.suggest_accessions(prefix)
        assert (rows[0]["accession"], rows[0]["gene"], rows[0]["species"]) == want
    merged = ncbi.suggest_accessions("NM_0010")
    assert len(merged) == 8 and merged == sorted(merged, key=lambda r: r["accession"])
    assert len({r["species"] for r in ncbi.suggest_accessions("NM_00100", limit=50)}) > 1


def test_gene_typeahead_is_per_species():
    mouse = ncbi.suggest_genes("Gapd", sm.get("mouse"))
    assert [g["symbol"] for g in mouse][:2] == ["Gapdh", "Gapdhs"]      # official spelling
    assert [g["symbol"] for g in ncbi.suggest_genes("gapd", sm.get("zebrafish"))][:1] == ["gapdh"]
    assert "GAPDH" in [g["symbol"] for g in ncbi.suggest_genes("GAPD")]


# ---------------------------------------------------------------- HTTP

def test_the_feature_is_advertised_and_the_species_listed():
    assert "multi_species" in FEATURES
    assert "multi_species" in client.get("/health").json()["features"]
    listed = client.get("/species").json()["species"]
    assert [s["slug"] for s in listed] == [s.slug for s in sm.SPECIES]
    assert listed[1] == {"slug": "mouse", "tax_id": "10090", "scientific": "Mus musculus",
                         "common": "Mouse", "example": "Gapdh"}


def test_gene_endpoint_takes_a_species_and_defaults_to_human():
    r = client.get("/gene/Gapdh", params={"species": "mouse"})
    assert r.status_code == 200 and r.json()["gene"]["species"] == "mouse"
    assert client.get("/gene/Gapdh", params={"species": "Mus musculus"}).json()["gene"]["gene_id"] == "14433"
    human = client.get("/gene/GAPDH").json()["gene"]
    assert (human["species"], human["assembly"], human["gene_id"]) == ("human", "GRCh38", "2597")
    bad = client.get("/gene/Gapdh", params={"species": "unicorn"})
    assert bad.status_code == 400 and bad.json()["error"] == "BAD_SPECIES"
    gone = client.get("/gene/Gapdh", params={"species": "rat"})
    assert gone.status_code == 404 and gone.json()["error"] == "NOT_PLACED"


def test_a_punctuated_fly_symbol_survives_the_url():
    r = client.get("/gene/Su(var)3-9", params={"species": "fly"})
    assert r.status_code == 200 and r.json()["gene"]["symbol"] == "Su(var)3-9"


def test_typeahead_endpoints():
    r = client.get("/suggest_genes", params={"q": "Gapd", "species": "mouse"}).json()
    assert r["species"] == "mouse" and r["suggestions"][0]["symbol"] == "Gapdh"
    assert client.get("/suggest_genes", params={"q": "Gapd", "species": "unicorn"}).status_code == 400
    rows = client.get("/suggest", params={"q": "NM_008084"}).json()["suggestions"]
    assert rows[0] == {"accession": "NM_008084.4", "gene": "Gapdh", "species": "mouse"}


def test_analysis_over_http_carries_the_species():
    body = client.get("/analyze/NM_001411845.1").json()
    assert (body["gene"]["species"], body["gene"]["organism"]) == ("mouse", "Mus musculus")
    assert body["meta"]["features"].count("multi_species") == 1
