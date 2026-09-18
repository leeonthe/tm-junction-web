"""The species the engine analyzes, and how each is named to NCBI.

Every NCBI query is scoped by TAXONOMY ID, never by a species name. Names look equivalent
and are not: NCBI files Saccharomyces cerevisiae's genes under the reference strain, S288C
(559292), so a Datasets lookup under "Saccharomyces cerevisiae" (4932) finds no gene at
all, while the same lookup under 559292 finds every one. The id is the species as NCBI
indexes it; the names here are for people.

Gene-symbol capitalization is deliberately NOT encoded here. The conventions differ by
community — human and yeast all capitals (GAPDH, TDH3), mouse and rat first letter only
(Gapdh), zebrafish all lower case (gapdh), fly mixed and meaningful (w, N, Gapdh1) — but
NCBI's symbol lookup is case-insensitive and answers with the official spelling, so the
engine asks as typed and shows what NCBI returns rather than guessing a convention.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Species:
    slug: str          # the API's `species` value and the cache-file suffix
    tax_id: str        # NCBI Taxonomy id — what every NCBI query is scoped by
    scientific: str    # "Mus musculus"
    common: str        # "Mouse"
    example: str       # a gene symbol in the species' own convention


HUMAN = Species("human", "9606", "Homo sapiens", "Human", "GAPDH")

# Order is the order the UI lists them in.
SPECIES: tuple[Species, ...] = (
    HUMAN,
    Species("mouse", "10090", "Mus musculus", "Mouse", "Gapdh"),
    Species("rat", "10116", "Rattus norvegicus", "Rat", "Actb"),
    Species("fly", "7227", "Drosophila melanogaster", "Fruit fly", "Gapdh1"),
    Species("yeast", "559292", "Saccharomyces cerevisiae", "Baker's yeast", "ACT1"),
    Species("zebrafish", "7955", "Danio rerio", "Zebrafish", "gapdh"),
)

_BY_SLUG = {s.slug: s for s in SPECIES}
_BY_TAX = {s.tax_id: s for s in SPECIES}
# What a person might reasonably send instead of the slug. Yeast's species-level id and
# its strain-qualified name both mean the one yeast NCBI has genes for.
_ALIASES = {
    **{s.scientific.lower(): s for s in SPECIES},
    **{s.common.lower(): s for s in SPECIES},
    "4932": _BY_SLUG["yeast"],
    "saccharomyces cerevisiae s288c": _BY_SLUG["yeast"],
    "fruit fly": _BY_SLUG["fly"], "drosophila": _BY_SLUG["fly"],
    "baker's yeast": _BY_SLUG["yeast"], "bakers yeast": _BY_SLUG["yeast"],
}


class UnknownSpecies(ValueError):
    pass


def supported_names() -> str:
    """"Human, Mouse, Rat, …" — for messages that have to say what IS supported."""
    return ", ".join(s.common for s in SPECIES)


def get(name: str | None) -> Species:
    """Slug, taxonomy id, scientific or common name -> Species. Empty means human, so every
    caller that predates species support keeps its meaning."""
    key = (name or "").strip()
    if not key:
        return HUMAN
    found = _BY_SLUG.get(key.lower()) or _BY_TAX.get(key) or _ALIASES.get(key.lower())
    if not found:
        raise UnknownSpecies(f"Unsupported species “{key}”. Supported: {supported_names()}.")
    return found


def by_tax_id(tax_id: str | int | None) -> Species | None:
    """The supported species an NCBI record belongs to, or None. A missing id means human:
    every record cached before species support was human by construction."""
    key = str(tax_id or "").strip()
    if not key:
        return HUMAN
    return _BY_TAX.get(key) or _ALIASES.get(key)
