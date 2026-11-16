from __future__ import annotations

import hashlib
import re
from dataclasses import asdict, dataclass, field


def digest(value: str) -> str:
    return hashlib.sha256(value.encode()).hexdigest()


def publication_id(value: str) -> str:
    """Normalize publication IDs, never application numbers or ambiguous bare digits."""
    result = re.sub(r"[\s.,\-/]", "", value).upper()
    if not re.fullmatch(r"[A-Z]{2}\d{4,12}(?:[A-Z]\d?)?", result):
        raise ValueError(f"Expected a publication number such as US7662379B2, got {value!r}")
    return result


def publication_parts(value: str) -> tuple[str, str, str]:
    match = re.fullmatch(r"([A-Z]{2})(\d+)([A-Z]\d?)?", publication_id(value))
    assert match
    return match[1], match[2], match[3] or ""


def same_publication(a: str, b: str) -> bool:
    ac, an, ak = publication_parts(a)
    bc, bn, bk = publication_parts(b)
    return (ac, an) == (bc, bn) and (not ak or not bk or ak == bk)


def protein_sequence(value: str) -> str:
    result = re.sub(r"\s+", "", value).upper()
    if not result or not re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWYBXZJUO]+", result):
        raise ValueError("Invalid amino-acid sequence; residues are never silently repaired")
    return result


@dataclass
class TargetProfile:
    canonical_name: str
    synonyms: list[str]
    confounders: list[str]
    conditional_terms: list[str] = field(default_factory=list)
    target_species: str = "human"
    notes: str = ""
    mocked: bool = True

    def __post_init__(self):
        if self.target_species != "human":
            raise ValueError("This version requires human target binding")
        if not self.canonical_name or not self.synonyms:
            raise ValueError("A profile needs a canonical name and at least one synonym")
        for term in self.synonyms + self.confounders + self.conditional_terms:
            if not isinstance(term, str) or not term.strip() or len(term) > 200:
                raise ValueError("Profile terms must be nonempty strings of at most 200 characters")


def mentions(text: str, terms: list[str]) -> list[tuple[int, int, str]]:
    result = []
    for term in terms:
        tokens = re.findall(r"[a-z]+|\d+", term.casefold())
        if not tokens:
            continue
        pattern = r"(?<![a-z0-9])" + r"[\s\-‐‑‒–—_/]*".join(map(re.escape, tokens)) + r"(?![a-z0-9])"
        result.extend((m.start(), m.end(), term) for m in re.finditer(pattern, text, flags=re.I))
    return result


def target_mentions(text: str, profile: TargetProfile) -> dict:
    negatives = mentions(text, profile.confounders)
    positives = [m for m in mentions(text, profile.synonyms)
                 if not any(m[0] >= n[0] and m[1] <= n[1] for n in negatives)]
    conditional = mentions(text, profile.conditional_terms)
    return {"matches": sorted({m[2] for m in positives}),
            "confounders": sorted({m[2] for m in negatives}),
            "conditional_terms": sorted({m[2] for m in conditional}),
            "possible_match": bool(positives or conditional)}


@dataclass
class Patent:
    publication: str
    title: str = ""
    abstract: str = ""
    cpc: list[str] = field(default_factory=list)
    applicants: list[str] = field(default_factory=list)
    assignees_original: list[str] = field(default_factory=list)
    assignees_current: list[str] = field(default_factory=list)
    publication_date: str = ""
    filing_date: str = ""
    priority_date: str = ""
    family_id: str = ""
    family: list[str] = field(default_factory=list)
    backward_citations: list[str] = field(default_factory=list)
    sections: list[dict] = field(default_factory=list)
    listing_urls: list[str] = field(default_factory=list)
    provenance: list[dict] = field(default_factory=list)
    gaps: list[str] = field(default_factory=list)

    @property
    def antibody_cpc(self) -> bool:
        return any(re.match(r"^C07K16/\d+$", re.sub(r"\s+", "", c).upper()) for c in self.cpc)

    @property
    def url(self) -> str:
        return f"https://patents.google.com/patent/{self.publication}/en"


@dataclass
class SequenceRecord:
    publication: str
    sequence: str
    accession: str = ""
    seq_id: str = ""
    title: str = ""
    source_url: str = ""
    source: str = "ncbi"
    created: str = ""
    updated: str = ""
    reference: str = ""

    @property
    def id(self) -> str:
        return digest("|".join([self.publication, self.accession, self.seq_id, self.sequence, self.source]))


@dataclass
class Domain:
    id: str
    record_id: str
    publication: str
    accession: str
    seq_id: str
    chain: str
    sequence: str
    cdrh3: str
    start: int
    end: int
    score: float
    complete: bool
    numbering: list[list]


def canonical_key(vh: str, vl: str) -> str:
    return "ab_" + digest(vh + "|" + vl)


def serializable(value):
    return asdict(value)
