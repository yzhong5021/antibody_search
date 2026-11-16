"""Small live public-source probes, independent of OPS credentials. Not a full live run."""
from __future__ import annotations

import argparse
import json
import os
import re
from pathlib import Path

from bs4 import BeautifulSoup

from patent_pipeline.http import HTTP, RetrievalError
from patent_pipeline.parsers import parse_google
from patent_pipeline.sources import NCBI, OPS
from patent_pipeline.telemetry import Telemetry


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("runs/source_probe"))
    parser.add_argument("--cache", type=Path, default=Path(".cache"))
    args = parser.parse_args()
    t = Telemetry(args.out)
    client, rows = HTTP(args.cache, t), []
    ncbi = NCBI(client)
    for target, term, pub in [("IL-10", "interleukin-10", "US7662379B2"),
                              ("GM-CSF", "granulocyte macrophage colony stimulating factor", "US12291567B2"),
                              ("IL-12p70", "interleukin-12", "US7700739B2")]:
        row = {"target": target, "publication": pub}
        with t.stage(target + " NCBI target probe") as stage:
            begin = len(t.requests)
            records = ncbi.query(term + " AND gbdiv_pat[PROP]", 2)
            row["broad_protein_count"] = t.requests[begin].get("result_count")
            row["sample"] = [{"publication": r.publication, "accession": r.accession, "seq_id": r.seq_id,
                              "length": len(r.sequence), "title": r.title} for r in records]
            stage["output_count"] = len(records)
        with t.stage(target + " NCBI publication probe") as stage:
            records = ncbi.for_publication(pub)
            row["publication_protein_count"] = len(records)
            row["source_publications"] = sorted(set(r.publication for r in records))
            stage["output_count"] = len(records)
        with t.stage(target + " Google publication probe") as stage:
            url = f"https://patents.google.com/patent/{pub}/en"
            try:
                p = parse_google(client.get(url, group="google"), pub, url)
                row["google"] = {"title": bool(p.title), "abstract": bool(p.abstract), "cpc": len(p.cpc),
                    "claims": sum(s["section"] == "claims" for s in p.sections),
                    "description_blocks": sum(s["section"] == "description" for s in p.sections),
                    "original_assignees": p.assignees_original, "current_assignees": p.assignees_current,
                    "publication_date": p.publication_date, "family_publications": len(p.family),
                    "backward_citations": len(p.backward_citations), "linked_sequence_listings": len(p.listing_urls)}
                stage["output_count"] = 1
            except (RetrievalError, ValueError) as error:
                row["google"] = {"error": str(error)}
        rows.append(row)
    with t.stage("EP patent protein reference check") as stage:
        ep_records = ncbi.for_publication("EP2192132A2", 2)
        non_us = {"publication": "EP2192132A2", "records": [{"publication": r.publication,
            "accession": r.accession, "seq_id": r.seq_id, "journal": r.reference} for r in ep_records],
            "notice": "Verified patent references, not a country keyword match. Empty SEQ IDs are retained as unknown; the publication date is never used as a SEQ ID."}
        stage["output_count"] = len(ep_records)
    directories = []
    for name, url, pattern in [("GenBank PAT nucleotide", "https://ftp.ncbi.nlm.nih.gov/genbank/", r"gbpat\d+\.seq\.gz"),
        ("GenBank daily nucleotide updates", "https://ftp.ncbi.nlm.nih.gov/genbank/daily-nc/", r".*\.gz"),
        ("BLAST patent protein/nucleotide", "https://ftp.ncbi.nlm.nih.gov/blast/db/", r"(?:pataa|patnt)[^/]*")]:
        with t.stage(name) as stage:
            try:
                raw = client.get(url, group="ncbi", max_bytes=3_000_000).decode()
                links = re.findall(r'href="([^"]+)"', raw)
                matched = [link for link in links if re.fullmatch(pattern, link)]
                directories.append({"source": name, "url": url, "files": len(matched), "examples": matched[:6]})
                stage["output_count"] = len(matched)
            except (RetrievalError, ValueError) as error:
                directories.append({"source": name, "url": url, "error": str(error)})
    ep_url = "https://data.epo.org/publication-server/rest/v1.2/patents/EP3827016NWA1"
    with t.stage("EPO publication server") as stage:
        try:
            soup = BeautifulSoup(client.get(ep_url, group="epo_publication"), "html.parser")
            eps = {"url": ep_url, "formats": [{"label": a.get_text(" ", strip=True), "href": a["href"]} for a in soup.select("a[href]")]}
            stage["output_count"] = len(eps["formats"])
        except (RetrievalError, ValueError) as error:
            eps = {"url": ep_url, "error": str(error)}
    if os.getenv("EPO_OPS_KEY") and os.getenv("EPO_OPS_SECRET"):
        with t.stage("OPS authenticated metadata probe") as stage:
            try:
                p = OPS(client).biblio("EP1907421A2")
                ops = {"status": "ok", "publication": p.publication, "title": p.title, "cpc": p.cpc}
                stage["output_count"] = 1
            except RetrievalError as error:
                ops = {"status": "failed", "reason": str(error)}
    else:
        ops = {"status": "blocked", "reason": "EPO_OPS_KEY and EPO_OPS_SECRET are not configured; no authenticated OPS requests made"}
    result = {"kind": "bounded live source probes, not a complete pipeline smoke test", "targets": rows,
              "non_us_reference_control": non_us, "ftp_directories": directories, "epo_publication_server": eps, "epo_ops": ops}
    (args.out / "probe.json").write_text(json.dumps(result, indent=2) + "\n")
    t.write()
    lines = ["# Live public-source probes", "", result["kind"], "", "| Target | Broad protein hits | Proteins in control patent | Google text/metadata |",
        "|---|---:|---:|---|"]
    for row in rows:
        g = row["google"]
        lines.append(f"| {row['target']} | {row['broad_protein_count']} | {row['publication_protein_count']} | " +
            (f"{g['claims']} claim blocks; {g['description_blocks']} description blocks; {g['cpc']} CPC codes" if "claims" in g else "retrieval failed") + " |")
    lines += ["", "Broad protein counts include unrelated sequences and receptor/subunit matches; these are not antibody counts.",
              "", f"Verified non-US control: {len(ep_records)} proteins with EP2192132A2 references. Their accession metadata are preserved even when the source lacks a SEQ ID. A country keyword query alone does not establish patent authority.",
              "", f"OPS: **{ops['status']}**. {ops.get('reason', '')}", "",
              "See probe.json for metadata/sequence availability and FTP directory observations; metrics.json contains measured timings and request provenance. These probes do not validate authenticated OPS discovery, full-text retrieval, or citation expansion.", ""]
    (args.out / "report.md").write_text("\n".join(lines))
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
