from __future__ import annotations

import csv
import json
from dataclasses import asdict
from pathlib import Path

from .models import Domain, Patent, SequenceRecord, TargetProfile, canonical_key, same_publication

COLUMNS = ["canonical_key", "vh_seq", "vl_seq", "cdrh3", "target", "sources", "structures", "cluster_id", "pairing_status"]


def antibody_row(vh: Domain | None, vl: Domain | None, target: str, sources: list[str]) -> dict:
    heavy, light = vh.sequence if vh else "", vl.sequence if vl else ""
    return {"canonical_key": canonical_key(heavy, light), "vh_seq": heavy, "vl_seq": light,
        "cdrh3": vh.cdrh3 if vh else "", "target": target, "sources": sorted(set(sources)),
        "structures": [], "cluster_id": "", "pairing_status": "paired" if vh and vl else "unpaired"}


def deduplicate_rows(rows: list[dict]) -> list[dict]:
    result = {}
    for row in rows:
        key = (row["canonical_key"], row["target"])
        if key not in result:
            result[key] = {**row, "sources": list(row["sources"])}
        else:
            result[key]["sources"] = sorted(set(result[key]["sources"] + row["sources"]))
    return sorted(result.values(), key=lambda row: (row["target"], row["canonical_key"]))


def validate_review(review: dict, profile: TargetProfile, patents: list[Patent], domains: list[Domain], *, allow_mock=False):
    """Structural/evidence checks; binding interpretation belongs to the reviewer, not a regex."""
    if review.get("mocked") and not allow_mock:
        return [], ["Mock review cannot confirm live antibodies"]
    if review.get("status") != "confirmed":
        return [], []
    catalogue = {d.id: d for d in domains}
    sections = {(p.publication, s["section"], s["locator"]): s["text"] for p in patents for s in p.sections}
    rows, errors = [], []
    for assertion in review.get("assertions", []):
        try:
            vh = catalogue[assertion["vh_domain_id"]] if assertion.get("vh_domain_id") else None
            vl = catalogue[assertion["vl_domain_id"]] if assertion.get("vl_domain_id") else None
            if not (vh or vl):
                raise ValueError("Assertion has no domains")
            if (vh and vh.chain != "H") or (vl and vl.chain not in ("K", "L")):
                raise ValueError("Wrong heavy/light domain type")
            if any(not d.complete for d in (vh, vl) if d):
                raise ValueError("Incomplete variable domain cannot be exported as confirmed")
            if any(set(d.sequence) - set("ACDEFGHIKLMNPQRSTVWY") for d in (vh, vl) if d):
                raise ValueError("Ambiguous/nonstandard residues require review")
            expected_pairing = "paired" if vh and vl else "unpaired"
            if assertion.get("pairing_status") != expected_pairing:
                raise ValueError("Pairing status contradicts selected domains")
            if assertion.get("human_target_binding") is not True:
                raise ValueError("Human target binding is not supported")
            evidence = assertion.get("evidence", [])
            if not evidence:
                raise ValueError("No source evidence")
            for item in evidence:
                quote = item.get("quote", "")
                body = sections.get((item["publication"], item["section"], item["locator"]), "")
                if len(quote) < 20 or quote not in body:
                    raise ValueError("Evidence quote/anchor not found in supplied source text")
            # Cross-publication sequence mappings must be explicit, not inferred from family membership.
            pubs = {e["publication"] for e in evidence}
            if any(not any(same_publication(d.publication, pub) for pub in pubs) for d in (vh, vl) if d):
                raise ValueError("Evidence does not map SEQ IDs in the domains' source publication")
            rows.append(antibody_row(vh, vl, profile.canonical_name,
                [f"https://patents.google.com/patent/{pub}/en" for pub in pubs]))
        except (KeyError, TypeError, ValueError) as error:
            errors.append(str(error))
    return deduplicate_rows(rows), errors


def write_tsv(path: Path, rows: list[dict], columns: list[str] | None = None):
    columns = columns or COLUMNS
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: json.dumps(v, ensure_ascii=False, separators=(",", ":")) if isinstance(v, (list, dict)) else v
                             for k, v in row.items() if k in columns})


def write_jsonl(path: Path, records: list):
    with path.open("w") as handle:
        for record in records:
            handle.write(json.dumps(record, ensure_ascii=False) + "\n")
