"""IMGT variable-domain identity, evidence-preserving aggregation, and H3 clusters."""
from __future__ import annotations

from collections import defaultdict
import copy
import hashlib
import importlib
from importlib.metadata import PackageNotFoundError, version
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys

from common import atomic_write, stable_json, value

AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")
HMM_ALPHABET = AMINO_ACIDS | set("BZXOU")
REGIONS = {"fr1": (1, 26), "cdr1": (27, 38), "fr2": (39, 55), "cdr2": (56, 65),
           "fr3": (66, 104), "cdr3": (105, 117), "fr4": (118, 128)}
TARGET_EVIDENCE = {"accession", "curated_name", "structural_contact", "exact_sequence"}


def clean_sequence(sequence):
    sequence = value(sequence)
    return re.sub(r"\s+", "", str(sequence)).upper() if sequence is not None else None


def canonical_key(heavy, light):
    if not heavy and not light:
        return None
    # Explicit roles prevent H-only and L-only records from sharing an identity.
    return hashlib.sha256(f"IMGT-v1|H:{heavy or ''}|L:{light or ''}".encode()).hexdigest()


class Numberer:
    def __init__(self, cache_dir):
        package = importlib.import_module("anarci")
        root = Path(package.__file__).parent
        binary = Path(sys.executable).parent / "hmmscan"
        if not binary.exists():
            binary = Path(shutil.which("hmmscan") or "/missing/hmmscan")
        if not binary.exists():
            raise RuntimeError("HMMER hmmscan is missing; use the supplied Conda environment")
        self.hmmerpath = str(binary.parent)
        hmmer = subprocess.run([str(binary), "-h"], capture_output=True, text=True, check=True).stdout
        try:
            package_version = version("anarci")
        except PackageNotFoundError:
            package_version = "unknown"
        fingerprints = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in (
            root / "anarci.py", root / "schemes.py", root / "dat" / "HMMs" / "ALL.hmm")}
        self.engine = {"anarci_version": package_version, "hmmer_version": next(
            (line for line in hmmer.splitlines() if "HMMER" in line), "unknown"),
            "fingerprints": fingerprints, "scheme": "imgt", "bit_score_threshold": 80,
            "allowed_species": None, "normalization_version": 1}
        digest = hashlib.sha256(stable_json(self.engine).encode()).hexdigest()
        self.cache_path = Path(cache_dir) / "numbering" / (digest + ".json")
        self.cache = json.loads(self.cache_path.read_text()) if self.cache_path.exists() else {}
        self.results = {}
        self.anarci = package.anarci

    def number_many(self, sequences):
        sequences = sorted({clean_sequence(s) for s in sequences if clean_sequence(s)})
        pending = []
        for seq in sequences:
            key = hashlib.sha256(seq.encode()).hexdigest()
            if key in self.cache:
                self.results[seq] = self.cache[key]
            elif set(seq) - HMM_ALPHABET:
                self.results[seq] = {"domains": [], "error": "invalid_amino_acid_characters"}
            else:
                pending.append((key, seq))
        if pending:
            numbered, details, _ = self.anarci(pending, scheme="imgt", allow={"H", "K", "L"},
                hmmerpath=self.hmmerpath, ncpu=1, allowed_species=None, assign_germline=False)
            if len(numbered) != len(pending) or len(details) != len(pending):
                raise RuntimeError("ANARCI returned a different number of results than input sequences")
            for (key, seq), domains, annotations in zip(pending, numbered, details):
                result = {"domains": [], "error": None if domains else "no_antibody_domain_numbered"}
                for domain, annotation in zip(domains or [], annotations or []):
                    positions, start, end = domain
                    numbering = [{"position": int(pos), "insertion": insertion.strip(), "residue": residue}
                                 for (pos, insertion), residue in positions]
                    ungapped = "".join(p["residue"] for p in numbering if p["residue"] != "-")
                    regions = {name: "".join(p["residue"] for p in numbering
                               if left <= p["position"] <= right and p["residue"] != "-")
                               for name, (left, right) in REGIONS.items()}
                    occupied = [p["position"] for p in numbering if p["residue"] != "-"]
                    result["domains"].append({"chain_type": annotation["chain_type"],
                        "sequence": ungapped, "numbering": numbering, "regions": regions,
                        "raw_start": int(start), "raw_end_inclusive": int(end),
                        "complete_variable_domain": bool(occupied and min(occupied) <= 1 and max(occupied) >= 127),
                        "quality_flags": ["ambiguous_variable_residues"] if set(ungapped) - AMINO_ACIDS else []})
                self.results[seq] = result
                self.cache[key] = result
            atomic_write(self.cache_path, stable_json(self.cache).encode())
        return self.results


def merge_records(records):
    """A feedback query can rediscover an observation; it cannot duplicate that assay."""
    merged = {}
    for rec in records:
        key = (rec["source"], rec["source_id"])
        if key not in merged:
            merged[key] = copy.deepcopy(rec)
        else:
            previous = merged[key]
            if any(previous.get(k) != rec.get(k) for k in ("vh_seq", "vl_seq", "raw", "structures")):
                raise ValueError(f"Conflicting source record ID: {key}; refusing to discard a variant")
            previous["query_terms"] = unique_objects(previous.get("query_terms", []) + rec.get("query_terms", []))
    return [merged[k] for k in sorted(merged)]


def unique_objects(objects):
    return [obj for _, obj in sorted({stable_json(o): o for o in objects}.items())]


def assign_domains(rec, numbered):
    result, reasons = {"H": None, "L": None}, []
    raw_h, raw_l = clean_sequence(rec.get("vh_seq")), clean_sequence(rec.get("vl_seq"))
    extras = [clean_sequence(s) for s in rec.get("unassigned_sequences", []) if clean_sequence(s)]
    is_scfv = "scfv" in str(rec.get("antibody_format") or "").lower()
    for role, seq in (("H", raw_h), ("L", raw_l)):
        if not seq:
            continue
        output = numbered[seq]
        if output["error"]:
            reasons.append("numbering_failed:" + output["error"])
            continue
        domains = output["domains"]
        # A documented single Fv can encode H and L on the same polypeptide.
        scfv_pair = is_scfv and role == "H" and not raw_l and len(domains) == 2 and sum(d["chain_type"] == "H" for d in domains) == 1
        if len(domains) != 1 and not scfv_pair:
            reasons.append("multiple_variable_domains_without_arm_mapping")
            continue
        for domain in domains:
            found_role = "H" if domain["chain_type"] == "H" else "L"
            if (found_role != role and not scfv_pair) or result[found_role] is not None:
                reasons.append("chain_role_mismatch")
                continue
            result[found_role] = domain
            reasons += domain["quality_flags"]
    if extras:
        reasons.append("ambiguous_chain_assignment")
    if result["H"] and result["L"] and not rec.get("pairing") and not is_scfv:
        reasons.append("pairing_not_documented")
    return result, sorted(set(reasons))


def observed_target_context(rec, target):
    matched = sorted(set(rec.get("target_accessions", [])) & set(target["accessions"]))
    if len(matched) == 1:
        rec["target_subunit"] = target["components"][matched[0]]["gene"]
        rec["target_entity"] = rec["target_subunit"]
    raw = rec.get("raw", {})
    text = stable_json({"target": rec.get("target"), "antigen": raw.get("curated_source_antigen"),
                        "pdb_name": raw.get("pdb_name")}).lower()
    for label, pattern in (("IL23", r"\bil[- ]?23\b|interleukin[- ]23"),
                           ("IL35", r"\bil[- ]?35\b|interleukin[- ]35"),
                           ("IL12p70", r"\bil[- ]?12\s*p70\b")):
        if re.search(pattern, text):
            rec["complex_state"] = label
            rec["target_entity"] = label
            if target["name"].lower() == "il12p70" and label in ("IL23", "IL35"):
                rec["candidate_reasons"].append("related_complex_without_requested_complex_binding")
            break


def normalize_records(records, target, numberer):
    records = merge_records(records)
    sequences = [s for r in records for s in [r.get("vh_seq"), r.get("vl_seq"), *r.get("unassigned_sequences", [])] if value(s) is not None]
    numbered = numberer.number_many(sequences)
    observations = []
    for rec in records:
        observed_target_context(rec, target)
        domains, reasons = assign_domains(rec, numbered)
        rec["numbering"] = domains
        rec["normalized_vh_seq"] = domains["H"]["sequence"] if domains["H"] else None
        rec["normalized_vl_seq"] = domains["L"]["sequence"] if domains["L"] else None
        rec["observation_id"] = f"{rec['source']}:{rec['source_id']}"
        rec["requested_target"] = target["name"]
        rec["candidate_reasons"] += reasons
        if rec["target_evidence"] not in TARGET_EVIDENCE:
            rec["candidate_reasons"].append("target_binding_not_established")
        rec["candidate_reasons"] = sorted(set(rec["candidate_reasons"]))
        rec["canonical_key"] = canonical_key(
            domains["H"]["sequence"] if domains["H"] else None,
            domains["L"]["sequence"] if domains["L"] else None) if not reasons else None
        rec["accepted_for_target"] = not rec["candidate_reasons"]
        rec["identity_resolution"] = "sequence" if rec["canonical_key"] else None
        observations.append(rec)

    # An entry with no inline sequence may be linked through a PDB containing one
    # unique antibody. Multiple different antibodies/variants never get guessed.
    structural_keys = defaultdict(set)
    for rec in observations:
        if rec["source"] == "SAbDab2" and rec["canonical_key"] and rec["accepted_for_target"]:
            for structure in rec["structures"]:
                structural_keys[structure["pdb_id"]].add(rec["canonical_key"])
    for rec in observations:
        if rec["canonical_key"] or rec.get("vh_seq") or rec.get("vl_seq") or not rec["accepted_for_target"] or not rec.get("name"):
            continue
        keys = {key for s in rec["structures"] for key in structural_keys.get(s["pdb_id"], set())}
        if len(keys) == 1:
            rec["canonical_key"] = next(iter(keys))
            rec["identity_resolution"] = "unambiguous_structural_cross_reference"

    groups = defaultdict(list)
    for rec in observations:
        if rec["accepted_for_target"] and not rec.get("entity_only") and (rec["canonical_key"] or rec.get("name")):
            groups[rec["canonical_key"] or rec["observation_id"]].append(rec)
    by_key = {r["canonical_key"]: r for r in observations if r["canonical_key"] and any(r["numbering"].values()) and r["accepted_for_target"]}
    antibodies = []
    for identifier, evidence in sorted(groups.items()):
        key = evidence[0]["canonical_key"]
        representative = by_key[key] if key else evidence[0]
        h, l = representative["numbering"]["H"], representative["numbering"]["L"]
        formats = sorted({r["antibody_format"] for r in evidence if r["antibody_format"]})
        single_domain = any(re.search(r"\bvhh\b|single.?domain|^sd[-_]?[hl]$", f, re.I) for f in formats)
        complete = bool((h and l and h["complete_variable_domain"] and l["complete_variable_domain"]) or
                        (single_domain and (h or l) and (h or l)["complete_variable_domain"]))
        antibodies.append({"antibody_id": identifier, "canonical_key": key,
            "vh_seq": h["sequence"] if h else None, "vl_seq": l["sequence"] if l else None,
            "imgt_vh": h["numbering"] if h else None, "imgt_vl": l["numbering"] if l else None,
            **{("cdrh" if role == "H" else "cdrl") + str(n): d["regions"]["cdr" + str(n)] if d else None
               for role, d in (("H", h), ("L", l)) for n in (1, 2, 3)},
            "frameworks": {"VH": {k: v for k, v in h["regions"].items() if k.startswith("fr")} if h else None,
                           "VL": {k: v for k, v in l["regions"].items() if k.startswith("fr")} if l else None},
            "names": sorted({r["name"] for r in evidence if r["name"]}),
            "aliases": sorted({a for r in evidence for a in r["aliases"] if value(a)}),
            "evidence_bucket": "paired" if h and l else "heavy_only" if h else "light_only" if l else "named_no_sequence",
            "antibody_formats": formats, "is_complete": complete, "target": target["name"],
            "target_entities": sorted({r["target_entity"] for r in evidence if r["target_entity"]}),
            "target_species": sorted({r["target_species"] for r in evidence if r["target_species"] is not None}),
            "sources": sorted({r["source"] for r in evidence}),
            "observation_ids": sorted(r["observation_id"] for r in evidence),
            "structures": unique_objects(s for r in evidence for s in r["structures"]),
            "cluster_id": None, "representative_member": None, "cluster_size": 0})
    cluster_h3(antibodies)
    candidates = [r for r in observations if not r["accepted_for_target"]]
    return antibodies, observations, candidates, {"unique_input_sequences": len(numbered),
        "numbered_input_sequences": sum(bool(r["domains"]) for r in numbered.values()),
        "failed_input_sequences": sum(bool(r["error"]) for r in numbered.values()),
        "canonical_antibodies": sum(a["canonical_key"] is not None for a in antibodies),
        "named_no_sequence": sum(a["canonical_key"] is None for a in antibodies),
        "candidate_observations": len(candidates), "observations": len(observations),
        "numbering_engine": numberer.engine}


def global_identity(a, b):
    if a == b:
        return 1.0
    if not a or not b:
        return 0.0
    from Bio.Align import PairwiseAligner

    aligner = PairwiseAligner(mode="global", match_score=1, mismatch_score=-1,
                              open_gap_score=-1, extend_gap_score=-1)
    alignment = aligner.align(a, b)[0]
    left, right = str(alignment[0]), str(alignment[1])
    return sum(x == y and x != "-" for x, y in zip(left, right)) / len(left)


def cluster_h3(antibodies):
    representatives, members = [], defaultdict(list)
    for antibody in sorted(antibodies, key=lambda a: a["canonical_key"] or a["antibody_id"]):
        h3 = antibody.get("cdrh3")
        if not h3:
            continue
        chosen = None
        for rep in representatives:
            other = rep["cdrh3"]
            if min(len(h3), len(other)) / max(len(h3), len(other)) >= 0.95 and global_identity(h3, other) >= 0.95:
                chosen = rep
                break
        if chosen is None:
            representatives.append(antibody)
            chosen = antibody
        members[chosen["canonical_key"]].append(antibody)
    for representative, group in members.items():
        for antibody in group:
            antibody.update(cluster_id="cdrh3:" + representative,
                             representative_member=representative, cluster_size=len(group))
