from __future__ import annotations

import importlib.metadata
import json
import platform
from dataclasses import asdict
from pathlib import Path

from .domains import number_domains, preflight
from .export import COLUMNS, antibody_row, deduplicate_rows, validate_review, write_jsonl, write_tsv
from .http import RetrievalError
from .llm import profile_prompt, review_patent, review_request
from .models import Patent, SequenceRecord, TargetProfile, target_mentions
from .sources import FIXTURES, LiveSources, merge_patents
from .telemetry import Telemetry


def family_order(publications: list[str]) -> list[str]:
    return sorted(publications, key=lambda p: ({"EP": 0, "WO": 1, "US": 2}.get(p[:2], 3), p))


class Pipeline:
    def __init__(self, source, telemetry: Telemetry, cache: Path,
                 max_patents: int | None = None, max_citations: int | None = None):
        self.source, self.t, self.cache = source, telemetry, cache
        self.max_patents, self.max_citations = max_patents, max_citations
        self.patents: dict[str, Patent] = {}
        self.records: dict[str, SequenceRecord] = {}
        self.domains = []
        self.decisions = []
        self.reviews, self.requests = [], []
        self.confirmed, self.mock_rows = [], []
        self.survivors: list[Patent] = []
        self.status = "running"

    def process(self, ids: list[str], profile: TargetProfile, prefix: str, limit: int | None):
        selected = []
        inspected = 0
        with self.t.stage(prefix + ": metadata and text filter") as step:
            step["input_count"] = len(ids)
            for pub in dict.fromkeys(ids):
                if pub in self.patents:
                    continue
                if limit is not None and inspected >= limit:
                    self.t.issue("limit", f"{prefix} processing capped at {limit} candidate publications (at most {limit} families)")
                    break
                inspected += 1
                try:
                    patent = self.source.patent(pub)
                    self.patents[pub] = patent
                    if not patent.cpc and isinstance(self.source, LiveSources):
                        try:
                            merge_patents(patent, self.source.google(pub))
                        except (RetrievalError, ValueError) as error:
                            self.t.issue("metadata_gap", str(error), pub)
                    decision = {"publication": pub, "route": prefix, "cpc": patent.cpc}
                    if not patent.antibody_cpc:
                        decision["status"] = "missing_cpc" if not patent.cpc else "non_antibody_cpc"
                        self.decisions.append(decision)
                        continue
                    match = target_mentions(patent.title + "\n" + patent.abstract, profile)
                    if not match["possible_match"]:
                        self.source.fulltext(patent)
                        match = target_mentions("\n".join(s["text"] for s in patent.sections if s["section"] == "description"), profile)
                        decision["match_field"] = "description"
                    else:
                        decision["match_field"] = "title_or_abstract"
                    decision.update(match)
                    if not match["possible_match"]:
                        decision["status"] = "no_target_match" if patent.sections else "target_unresolved_missing_text"
                        self.decisions.append(decision)
                        continue
                    decision["status"] = "sequence_retrieval"
                    self.decisions.append(decision)
                    selected.append(patent)
                except (RetrievalError, ValueError) as error:
                    self.t.issue("metadata_error", str(error), pub)
            step["output_count"] = len(selected)
        group_records = []
        owners: dict[str, set[str]] = {}
        with self.t.stage(prefix + ": sequence recovery") as step:
            step["input_count"] = len(selected)
            for patent in selected:
                records = []
                try:
                    records = self.source.sequences(patent.publication)
                except (RetrievalError, ValueError) as error:
                    self.t.issue("sequence_gap", str(error), patent.publication)
                records.extend(self.source.listings(patent))
                if not records:
                    # The original might be pre-grant or lack indexed proteins. Recover from the first useful family member.
                    for family_pub in family_order(self.source.family(patent)):
                        try:
                            records = self.source.sequences(family_pub)
                            family_patent = self.source.patent(family_pub)
                            self.patents[family_pub] = family_patent
                            if not records:
                                self.source.google_text(family_patent)
                            records.extend(self.source.listings(family_patent))
                            if records:
                                self.t.issue("family_recovery", f"Recovered {len(records)} records from {family_pub}; SEQ IDs retain that publication", patent.publication)
                                break
                        except (RetrievalError, ValueError) as error:
                            self.t.issue("sequence_gap", str(error), family_pub)
                if not records:
                    self.source.google_text(patent)
                    records.extend(self.source.listings(patent))
                if not records:
                    self.t.issue("sequence_gap", "No retrievable protein sequences; not evidence of absence", patent.publication)
                for record in records:
                    self.records[record.id] = record
                    group_records.append(record)
                    owners.setdefault(record.id, set()).add(patent.publication)
            group_records = list({r.id: r for r in group_records}.values())
            step["output_count"] = len(group_records)
        with self.t.stage(prefix + ": ANARCI") as step:
            step["input_count"] = len(group_records)
            domains = number_domains(group_records, self.cache, self.t)
            self.domains.extend(domains)
            step["output_count"] = len(domains)
        passing = {pub for domain in domains for pub in owners.get(domain.record_id, set())}
        missing = [p for p in selected if p.publication not in passing]
        # Indexed proteins can be antigens/CDRs only. Look for accessible listings before rejecting the patent.
        if missing and isinstance(self.source, LiveSources):
            rescued = []
            with self.t.stage(prefix + ": missing-domain listing rescue") as step:
                step["input_count"] = len(missing)
                for patent in missing:
                    self.source.google_text(patent)
                    for record in self.source.listings(patent):
                        if record.id not in self.records:
                            self.records[record.id] = record
                            rescued.append(record)
                        owners.setdefault(record.id, set()).add(patent.publication)
                step["output_count"] = len(rescued)
            if rescued:
                with self.t.stage(prefix + ": ANARCI rescued listings") as step:
                    step["input_count"] = len(rescued)
                    extra_domains = number_domains(rescued, self.cache, self.t)
                    domains.extend(extra_domains)
                    self.domains.extend(extra_domains)
                    step["output_count"] = len(extra_domains)
            passing = {pub for domain in domains for pub in owners.get(domain.record_id, set())}
        survivors = [p for p in selected if p.publication in passing]
        for decision in self.decisions:
            if decision["publication"] in {p.publication for p in selected}:
                decision["status"] = "has_antibody_domains" if decision["publication"] in passing else "no_antibody_domains"
        return survivors

    def run(self, profile: TargetProfile, seeds: list[str]):
        out = self.t.out
        (out / "target_profile.json").write_text(json.dumps(asdict(profile), indent=2) + "\n")
        (out / "profile_prompt.txt").write_text(profile_prompt(profile.canonical_name))
        try:
            with self.t.stage("preflight") as step:
                preflight()
                try:
                    self.source.preflight()
                except RetrievalError:
                    self.status = "blocked"
                    raise
                step["output_count"] = 1
            with self.t.stage("discovery") as step:
                ids, discovered_records = self.source.discover(profile, self.max_patents)
                # Discovery sequences locate publications; every accepted publication is independently retrieved in full.
                step["input_count"] = 1
                step["output_count"] = len(ids)
                step["protein_discovery_records"] = len(discovered_records)
                ids = list(dict.fromkeys(seeds + ids))
            roots = self.process(ids, profile, "seeds", self.max_patents)
            with self.t.stage("one-hop citation expansion") as step:
                step["input_count"] = len(roots)
                citations = []
                for patent in roots:
                    citations.extend(self.source.citations(patent, self.max_citations))
                citations = [p for p in dict.fromkeys(citations) if p not in self.patents]
                step["output_count"] = len(citations)
            cited = self.process(citations, profile, "citations", self.max_citations) if citations else []
            self.survivors = roots + cited
            with self.t.stage("full patent text and evidence review") as step:
                step["input_count"] = len(self.survivors)
                for patent in self.survivors:
                    self.source.fulltext(patent)
                    review_patents = [patent]
                    if not {"claims", "description"} <= {s["section"] for s in patent.sections}:
                        for pub in family_order(self.source.family(patent)):
                            if not pub.startswith(("EP", "WO")):
                                continue
                            try:
                                sibling = self.source.fulltext(self.source.patent(pub))
                                self.patents[pub] = sibling
                                if sibling.sections:
                                    review_patents.append(sibling)
                                if {"claims", "description"} <= {s["section"] for s in sibling.sections}:
                                    break
                            except (RetrievalError, ValueError) as error:
                                self.t.issue("text_gap", str(error), pub)
                        self.source.google_text(patent)
                    # The full page can expose listing URLs absent from the bibliographic response.
                    for document in review_patents:
                        if document.listing_urls:
                            new_records = [r for r in self.source.listings(document) if r.id not in self.records]
                            if new_records:
                                self.records.update((r.id, r) for r in new_records)
                                self.domains.extend(number_domains(new_records, self.cache, self.t))
                    for section in ("claims", "description"):
                        if not any(s["section"] == section for p in review_patents for s in p.sections):
                            self.t.issue("text_gap", f"Missing {section} across retrieved source and family text", patent.publication)
                    allowed_pubs = {p.publication for p in review_patents} | set(patent.family)
                    relevant_domains = [d for d in self.domains if d.publication in allowed_pubs]
                    request = review_request(profile, review_patents, relevant_domains)
                    self.requests.append(request)
                    review = review_patent(profile, review_patents, relevant_domains)
                    self.reviews.append({"publication": patent.publication, **review})
                    rows, errors = validate_review(review, profile, review_patents, relevant_domains)
                    self.confirmed.extend(rows)
                    for error in errors:
                        if not review.get("mocked"):
                            self.t.issue("review_error", error, patent.publication)
                    if self.source.fixture:
                        self.mock_demo(profile, review_patents, relevant_domains)
                step["output_count"] = len(self.confirmed)
                step["pending_patents"] = sum(r["status"] == "pending_review" for r in self.reviews)
            self.status = "fixture_complete" if self.source.fixture else "pending_review"
        except Exception as error:
            if self.status != "blocked":
                self.status = "failed"
            self.t.issue("run_error", str(error))
            raise
        finally:
            self.export(profile)
        return self.status

    def mock_demo(self, profile: TargetProfile, patents: list[Patent], domains):
        path = FIXTURES / "mock_reviews.json"
        if not path.exists():
            return
        mock = json.loads(path.read_text()).get(profile.canonical_name)
        if not mock:
            return
        def find(seq_id, chain):
            return next((d for d in domains if d.seq_id == str(seq_id) and d.chain in chain and d.complete), None)
        vh, vl = find(mock["vh_seq_id"], "H"), find(mock["vl_seq_id"], "KL")
        if not vh or not vl:
            raise ValueError("Fixture review references missing/incomplete variable domains")
        assertion = {"vh_domain_id": vh.id, "vl_domain_id": vl.id, "pairing_status": "paired", "human_target_binding": True,
                     "evidence": mock["evidence"]}
        review = {"status": "confirmed", "mocked": True, "assertions": [assertion]}
        rows, errors = validate_review(review, profile, patents, domains, allow_mock=True)
        if errors:
            raise ValueError("Invalid mock review fixture: " + "; ".join(errors))
        self.mock_rows.extend(rows)

    def export(self, profile: TargetProfile):
        out = self.t.out
        with self.t.stage("export") as step:
            candidates = []
            for domain in self.domains:
                record = self.records[domain.record_id]
                row = antibody_row(domain if domain.chain == "H" else None, domain if domain.chain != "H" else None,
                                   profile.canonical_name, [record.source_url, f"https://patents.google.com/patent/{record.publication}/en"])
                row.update(review_status="pending_review", complete_domain=domain.complete)
                candidates.append(row)
            candidates = deduplicate_rows(candidates)
            write_tsv(out / "confirmed.tsv", deduplicate_rows(self.confirmed))
            write_tsv(out / "candidates.tsv", candidates, COLUMNS + ["review_status", "complete_domain"])
            if self.source.fixture:
                write_tsv(out / "mock_confirmed.tsv", deduplicate_rows(self.mock_rows))
            else:
                (out / "mock_confirmed.tsv").unlink(missing_ok=True)
            write_jsonl(out / "evidence.jsonl", [{"domain": asdict(d), "sequence_record": asdict(self.records[d.record_id])} for d in self.domains])
            write_jsonl(out / "patents.jsonl", [asdict(p) for p in self.patents.values()])
            write_jsonl(out / "filter_decisions.jsonl", self.decisions)
            write_jsonl(out / "review_requests.jsonl", self.requests)
            write_jsonl(out / "reviews.jsonl", self.reviews)
            environment = {}
            for name in ["anarci", "biopython", "requests", "beautifulsoup4"]:
                try:
                    environment[name] = importlib.metadata.version(name)
                except importlib.metadata.PackageNotFoundError:
                    environment[name] = "missing"
            summary = {"status": self.status, "fixture": self.source.fixture, "llm_mocked": True,
                "target": profile.canonical_name, "patents_inspected": len(self.patents), "surviving_patents": len(self.survivors),
                "profile_notes": profile.notes,
                "sequences": len(self.records), "unique_sequences": len({r.sequence for r in self.records.values()}),
                "antibody_domains": len(self.domains), "complete_domains": sum(d.complete for d in self.domains),
                "unique_antibody_domains": len({(d.chain, d.sequence) for d in self.domains}),
                "unique_complete_domains": len({(d.chain, d.sequence) for d in self.domains if d.complete}),
                "confirmed_antibodies": len(deduplicate_rows(self.confirmed)), "candidate_rows": len(candidates),
                "mock_demonstration_rows": len(deduplicate_rows(self.mock_rows)),
                "review_payload_bytes": sum(len(json.dumps(request, ensure_ascii=False).encode()) for request in self.requests),
                "limits": {"seed_publications": self.max_patents, "citation_publications": self.max_citations},
                "python": platform.python_version(), "dependencies": environment}
            (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
            step["output_count"] = len(candidates)
        self.t.write()
        self.report(summary)

    def report(self, summary: dict):
        stats = self.t.snapshot()
        lines = [f"# {summary['target']} patent retrieval", "", f"Status: **{summary['status']}**. ZymoLLM calls are mocked.", "",
            "Live confirmed.tsv stays empty until evidence review is implemented. Candidate domains are not confirmed target binders, and an unpaired row does not establish a naturally single-domain antibody.", "",
            "## Counts", "", "| Measure | Count |", "|---|---:|"]
        if summary["profile_notes"]:
            lines[4:4] = ["Target profile: " + summary["profile_notes"], ""]
        for key in ("patents_inspected", "surviving_patents", "sequences", "unique_sequences", "antibody_domains", "unique_antibody_domains", "unique_complete_domains", "candidate_rows", "confirmed_antibodies", "mock_demonstration_rows"):
            lines.append(f"| {key.replace('_', ' ')} | {summary[key]} |")
        lines += ["", "Antibody domains counts source observations; unique antibody domains deduplicates identical chain sequences. Complete-domain counts use IMGT framework/CDR3 coverage, not target-binding evidence."]
        lines += ["", "## Retrieved patent content", "", "| Publication | Abstract | CPC codes | Assignee fields | Claims blocks | Description blocks | Listing links |", "|---|---:|---:|---:|---:|---:|---:|"]
        for p in self.patents.values():
            lines.append(f"| {p.publication} | {int(bool(p.abstract))} | {len(p.cpc)} | {int(bool(p.assignees_original or p.assignees_current))} | {sum(s['section']=='claims' for s in p.sections)} | {sum(s['section']=='description' for s in p.sections)} | {len(p.listing_urls)} |")
        lines += ["", "## Timings", "", f"Measured total: {stats['total_seconds']:.3f} seconds. Offline fixture measurements exclude live service latency.", "",
            "| Stage | Inputs | Outputs | Seconds |", "|---|---:|---:|---:|"]
        for stage in stats["stages"]:
            lines.append(f"| {stage['stage']} | {stage['input_count']} | {stage['output_count']} | {stage['seconds']:.4f} |")
        longest = max(stats["stages"], key=lambda s: s["seconds"], default=None)
        if longest:
            lines += ["", f"Longest measured stage: **{longest['stage']}**, {longest['seconds']:.3f}s."]
        req = stats["requests"]
        lines += ["", f"HTTP requests/attempts: {sum(r['attempts'] for r in req)}; cache hits: {sum(r['cache_hit'] for r in req)}; response bytes (including cached reads): {sum(r['bytes'] for r in req)}.",
            f"Network: {sum(r['network_seconds'] for r in req):.3f}s; rate-limit waiting: {sum(r['wait_seconds'] for r in req):.3f}s; retry waiting: {sum(r['retry_seconds'] for r in req):.3f}s.", "",
            "## Gaps and inefficiencies", "",
            "- Title/abstract discovery is primary. Description-only patents not reached through NCBI, seeds or citations may be missed.",
            "- CPC C07K16 descendants are mandatory. Missing CPC is unresolved, not proof that a patent lacks antibodies.",
            "- Family text and sequence IDs retain their own publications. Family membership never establishes pairing or identical claims.",
            "- Protein recovery stops at the first useful family member when the original has none. Other family versions may contain additional sequences.",
            "- NCBI title/reference matches can be antigen proteins, receptors, short CDRs, or germline controls. ANARCI identifies domains, not target specificity.",
            "- Accessible ST.25/ST.26 protein listings are parsed; inaccessible attachments, unannotated DNA and OCR-only sequences remain gaps.",
            "- If indexed proteins have no antibody domains, the known patent page is checked for listing links before rejection. Full-page retrieval may also expose additional listings after initial filtering.",
            "- Citation expansion is one hop. Both directions can add unrelated or already-seen families; search rate limits can dominate runtime.",
            "- Raw HTTP and ANARCI results are cached. --refresh re-fetches HTTP; historical raw cache is retained for provenance.",
            f"- Active candidate-publication limits (a conservative bound on family counts): {summary['limits']}. Inspect limit/incomplete events below; bounded smoke runs are not exhaustive searches.",
            "- structures and cluster_id are reserved; no structural lookup or similarity clustering is performed.", ""]
        lines += [f"Prepared review payloads: {summary['review_payload_bytes']:,} UTF-8 bytes total. They include complete retrieved patent sections; numbering arrays are excluded from model payloads. API latency, token use, and any need for chunking remain unmeasured while ZymoLLM is mocked.", ""]
        if self.source.fixture:
            lines += ["Offline data are bounded real patent/NCBI responses. Re-encoded listing fixtures and mock reviews are labelled demonstrations. Live discovery, OPS authentication and citation retrieval are not validated by this run.", ""]
        for issue in stats["issues"]:
            lines.append(f"- **{issue['kind']}** {issue['publication']}: {issue['message']}")
        for patent in self.patents.values():
            lines.extend(f"- **patent gap** {patent.publication}: {gap}" for gap in sorted(set(patent.gaps)))
        lines += ["", "## Estimates", "",
            "Planning ranges: offline smoke 10–60 seconds for all targets; bounded live smoke 3–15 minutes per target, 10–45 minutes total. Uncapped runs can take tens of minutes to several hours. These are estimates, not benchmarks; counts, throttling and retries determine actual runtime.", "",
            "Detailed per-request timings and source retrieval timestamps are in metrics.json and the HTTP cache metadata. Publication, filing and priority dates are separate from NCBI creation/update dates and cache retrieval dates.", "",
            "Source coverage and release/update information: see README.md. NCBI disclaimer: https://www.ncbi.nlm.nih.gov/About/disclaimer.html", ""]
        (self.t.out / "report.md").write_text("\n".join(lines))
