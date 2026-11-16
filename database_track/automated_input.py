#!/usr/bin/env python3
"""Aggregate target-specific antibody evidence. See README.md for formats and caveats."""
from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
from datetime import datetime, timezone
import hashlib
import io
import json
import logging
from pathlib import Path
import re
import time

from common import HttpCache, atomic_write, resolve_target, stable_json
from normalize import Numberer, merge_records, normalize_records
import sources

ROOT = Path(__file__).resolve().parent
LOG = logging.getLogger("zymo")
TSV_COLUMNS = ["antibody_id", "canonical_key", "vh_seq", "vl_seq", "cdrh1", "cdrh2", "cdrh3", "cdrl1", "cdrl2", "cdrl3",
               "names", "aliases", "evidence_bucket", "antibody_formats", "is_complete", "target", "target_entities",
               "target_species", "sources", "structures", "observation_ids", "cluster_id", "representative_member", "cluster_size"]


def retrieve(name, function, *args):
    started = time.monotonic()
    LOG.info("Retrieving %s", name)
    try:
        records, metadata = function(*args)
        status = {"status": "partial" if metadata.get("errors") else "success", "records": len(records),
                  "seconds": round(time.monotonic() - started, 3), **metadata}
        LOG.info("%s: %d records (%s)", name, len(records), status["status"])
        return records, status
    except Exception as exc:
        LOG.exception("%s failed", name)
        return [], {"status": "failed", "records": 0, "seconds": round(time.monotonic() - started, 3),
                    "errors": [{"type": type(exc).__name__, "error": str(exc)}]}


def feedback_names(records, target):
    terms = set()
    supported_keys = {r["canonical_key"] for r in records if r.get("accepted_for_target") and r.get("canonical_key")}
    existing = {a.casefold() for a in target["aliases"] + target["near_positive_terms"]}
    generic = {"antibody", "fab", "heavy chain", "light chain", "human", "mouse", "nanobody", "igg", "immunoglobulin"}
    for rec in records:
        if rec["source"] not in ("PLAbDab", "Thera-SAbDab", "SAbDab2"):
            continue
        if not rec.get("accepted_for_target") and rec.get("canonical_key") not in supported_keys:
            continue
        if rec["query_terms"] and all(m["query_class"] == "NEAR_POSITIVE" for m in rec["query_terms"]):
            continue
        # PLAbDab aliases include accession IDs; these are not clone names.
        names = [rec.get("name")]
        if rec["source"] == "Thera-SAbDab":
            names += rec.get("aliases", [])
        for name in names:
            if name and 3 <= len(name) <= 80 and re.search(r"[a-z]", name, re.I) and name.casefold() not in existing | generic:
                terms.add(name)
    return sorted(terms)


def write_outputs(directory, antibodies, observations, candidates, manifest):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    artifacts = {"antibodies.jsonl": antibodies, "observations.jsonl": observations, "candidates.jsonl": candidates}
    checksums = {}
    for filename, rows in artifacts.items():
        # Full IMGT maps belong in the antibody table. Repeating them on every
        # assay and keyword hit needlessly multiplies the exported data volume.
        export = rows if filename == "antibodies.jsonl" else (
            {k: v for k, v in r.items() if k != "numbering"} for r in rows)
        data = "".join(stable_json(r) + "\n" for r in export).encode()
        atomic_write(directory / filename, data)
        checksums[filename] = {"rows": len(rows), "sha256": hashlib.sha256(data).hexdigest()}
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=TSV_COLUMNS, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for antibody in antibodies:
        writer.writerow({k: stable_json(v) if isinstance(v, (dict, list)) else v for k, v in antibody.items() if k in TSV_COLUMNS})
    data = buffer.getvalue().encode()
    atomic_write(directory / "antibodies.tsv", data)
    checksums["antibodies.tsv"] = {"rows": len(antibodies), "sha256": hashlib.sha256(data).hexdigest()}
    manifest["artifacts"] = checksums
    atomic_write(directory / "manifest.json", (json.dumps(manifest, indent=2, ensure_ascii=False, allow_nan=False) + "\n").encode())


def pipeline_build(target, output_dir, *, asd_path=None, aacdb_enabled=False, offline=False,
                   cache_dir=None, refresh=False):
    """Return a manifest; incomplete enabled sources set manifest['status'] to 'partial'.

    Configuration/target-resolution errors raise rather than manufacturing an empty dataset.
    ASD is read only from a local release the caller is permitted to use.
    """
    started = time.monotonic()
    client = HttpCache(cache_dir or ROOT / ".cache", offline=offline, refresh=refresh)
    target = resolve_target(target, client)
    statuses, records = {}, []
    tasks = {"ABCD": (sources.fetch_abcd, client, target), "IEDB": (sources.fetch_iedb, client, target),
             "PLAbDab": (sources.fetch_plabdab, client, target), "Thera-SAbDab": (sources.fetch_thera, client, target),
             "RCSB": (sources.fetch_rcsb, client, target)}
    if asd_path:
        tasks["ASD"] = (sources.fetch_asd, asd_path, target)
    else:
        statuses["ASD"] = {"status": "disabled", "reason": "No permitted local ASD release supplied"}
    with ThreadPoolExecutor(max_workers=3) as pool:
        jobs = {pool.submit(retrieve, name, *args): name for name, args in tasks.items()}
        for future in as_completed(jobs):
            name = jobs[future]
            found, statuses[name] = future.result()
            records.extend(found)
    structural, statuses["SAbDab2"] = retrieve("SAbDab2", sources.fetch_sabdab, client, target, records)
    if statuses["RCSB"]["status"] != "success":
        statuses["SAbDab2"]["status"] = "partial"
        statuses["SAbDab2"].setdefault("errors", []).append({"error": "RCSB retrieval was incomplete; structural expansion coverage is incomplete"})
    records.extend(structural)
    if aacdb_enabled:
        annotations, statuses["AACDB"] = retrieve("AACDB", sources.fetch_aacdb, client, records)
        if isinstance(annotations, dict):
            for rec in records:
                if rec["source_id"] in annotations and rec["source"] == "SAbDab2":
                    rec["annotations"] = {"AACDB": annotations[rec["source_id"]]}
    else:
        statuses["AACDB"] = {"status": "disabled", "reason": "Optional annotation enrichment disabled"}

    normalization_error, terms = None, []
    try:
        numberer = Numberer(client.directory)
        LOG.info("Numbering initial antibody sequences")
        antibodies, observations, candidates, normalization = normalize_records(records, target, numberer)
        terms = feedback_names(observations, target)
        if terms:
            LOG.info("One feedback pass: %d names/aliases, Channel B only", len(terms))
            feedback_tasks = {"PLAbDab": sources.fetch_plabdab, "Thera-SAbDab": sources.fetch_thera}
            with ThreadPoolExecutor(max_workers=2) as pool:
                jobs = {pool.submit(retrieve, name + " feedback", fn, client, target, terms): name
                        for name, fn in feedback_tasks.items()}
                for future in as_completed(jobs):
                    name = jobs[future]
                    found, status = future.result()
                    records.extend(found)
                    statuses[name]["feedback"] = status
                    if status["status"] != "success":
                        statuses[name]["status"] = "partial"
            LOG.info("Normalizing final records")
            antibodies, observations, candidates, normalization = normalize_records(records, target, numberer)
    except Exception as exc:
        LOG.exception("Normalization failed; retaining retrieved evidence")
        normalization_error = {"type": type(exc).__name__, "error": str(exc)}
        antibodies, observations = [], merge_records(records)
        for rec in observations:
            rec.update(canonical_key=None, accepted_for_target=False)
            rec["candidate_reasons"] = sorted(set(rec["candidate_reasons"] + ["normalization_unavailable"]))
        candidates, normalization = observations, {"error": normalization_error}

    manifest = {"schema_version": 1, "status": "failed" if normalization_error else "partial" if any(
        s["status"] not in ("success", "disabled") for s in statuses.values()) else "complete",
        "created_at": datetime.now(timezone.utc).isoformat(), "seconds": round(time.monotonic() - started, 3),
        "options": {"offline": offline, "refresh": refresh, "asd_enabled": bool(asd_path), "aacdb_enabled": aacdb_enabled},
        "target": target, "sources": dict(sorted(statuses.items())), "normalization": normalization,
        "feedback": {"passes": int(bool(terms)), "channels": ["PLAbDab", "Thera-SAbDab"], "terms": terms},
        "evidence_buckets": dict(Counter(a["evidence_bucket"] for a in antibodies)),
        "cache": [client.used[k] for k in sorted(client.used)],
        "coverage_notes": ["Records exposed by these databases are not exhaustive literature coverage.",
            "Title-only, inferred-target, incompatible-species and unresolved-arm hits remain candidates.",
            "Unknown source species remain unknown; an HMM species prediction is not species evidence.",
            "Variable-domain identity does not imply identical Fc, payload, therapeutic product or binding measurements.",
            "H3 clustering is an annotation, not proof of clonal lineage or an independent train/test split.",
            "Negative assay outcomes retain their epitope, method and experimental context."]}
    write_outputs(output_dir, antibodies, observations, candidates, manifest)
    LOG.info("%s: %d antibody rows, %d observations, %d candidates (%s)", target["name"],
             len(antibodies), len(observations), len(candidates), manifest["status"])
    return manifest


def smoke_check(manifest, directory):
    checks = []
    antibodies = [json.loads(line) for line in (Path(directory) / "antibodies.jsonl").read_text().splitlines()]
    observations = [json.loads(line) for line in (Path(directory) / "observations.jsonl").read_text().splitlines()]
    keys = [a["canonical_key"] for a in antibodies if a["canonical_key"]]
    checks.append(("enabled_sources_complete", manifest["status"] == "complete"))
    checks.append(("canonical_antibodies_found", bool(keys)))
    checks.append(("canonical_keys_unique", len(keys) == len(set(keys))))
    checks.append(("at_most_one_feedback_pass", manifest["feedback"]["passes"] <= 1))
    checks.append(("uncertain_observations_excluded", all(not r["candidate_reasons"] for r in observations if r.get("accepted_for_target"))))
    name = manifest["target"]["name"]
    if name == "IL10":
        found = [a for a in antibodies if any("9d7" in n.lower() for n in a["names"])]
        checks.append(("9D7_cross_source_reconciliation", any({"ABCD", "IEDB", "SAbDab2"} <= set(a["sources"]) for a in found)))
    if name == "IL12p70":
        components = manifest["target"]["components"]
        checks.append(("correct_IL12_accessions", components["P29459"]["gene"] == "IL12A" and components["P29460"]["gene"] == "IL12B"))
        checks.append(("related_complexes_not_promoted", all(not r["accepted_for_target"] for r in observations if r.get("complex_state") in ("IL23", "IL35"))))
    if name == "CSF2":
        names = {n.lower() for a in antibodies for n in a["names"]}
        checks.append(("GMCSF_therapeutic_found", bool(names & {"namilumab", "otilimab", "gimsilumab"})))
        checks.append(("GMCSF_receptor_therapeutics_excluded", not names & {"mavrilimumab", "birinkibart", "trabikibart"}))
    return {"target": name, "passed": all(ok for _, ok in checks), "checks": dict(checks),
            "manifest": str(Path(directory) / "manifest.json"), "antibodies": len(antibodies),
            "observations": len(observations), "candidates": manifest["artifacts"]["candidates.jsonl"]["rows"]}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    for command in ("run", "smoke"):
        p = subparsers.add_parser(command)
        p.add_argument("--out", type=Path, required=True)
        p.add_argument("--cache", type=Path, default=ROOT / ".cache")
        mode = p.add_mutually_exclusive_group()
        mode.add_argument("--offline", action="store_true")
        mode.add_argument("--refresh", action="store_true")
        p.add_argument("--asd", type=Path, help="Permitted local ASD Delta table root")
        p.add_argument("--aacdb", action="store_true", help="Enable annotation enrichment only")
        if command == "run":
            p.add_argument("--target", type=Path, required=True)
        else:
            p.add_argument("--live", action="store_true", help="Explicitly allow network access for the smoke suite")
    args = parser.parse_args(argv)
    if args.command == "smoke" and not (args.live or args.offline):
        parser.error("smoke requires --live or --offline")
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s", datefmt="%H:%M:%S")
    try:
        options = dict(asd_path=args.asd, aacdb_enabled=args.aacdb, offline=args.offline, cache_dir=args.cache, refresh=args.refresh)
        if args.command == "run":
            manifest = pipeline_build(args.target, args.out, **options)
            return 0 if manifest["status"] == "complete" else 1
        reports = []
        for preset in ("il10", "il12p70", "csf2"):
            directory = args.out / preset
            manifest = pipeline_build(ROOT / "targets" / (preset + ".json"), directory, **options)
            reports.append(smoke_check(manifest, directory))
        report = {"passed": all(r["passed"] for r in reports), "targets": reports}
        atomic_write(args.out / "smoke_report.json", (json.dumps(report, indent=2) + "\n").encode())
        print(json.dumps(report, indent=2))
        return 0 if report["passed"] else 1
    except Exception as exc:
        LOG.exception("Pipeline could not start: %s", exc)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
