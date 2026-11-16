"""Explicit fixture refresh; never imported or run by tests. Downloads small public samples."""
from __future__ import annotations

import gzip
import json
import time
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parents[1] / "patent_pipeline" / "fixtures"
TARGETS = {"IL-10": "US7662379B2", "GM-CSF": "US12291567B2", "IL-12p70": "US7700739B2"}


def main():
    ROOT.mkdir(parents=True, exist_ok=True)
    manifest = {"notice": "Real public responses, bounded for offline regression tests; not an exhaustive search.", "targets": {}}
    for target, pub in TARGETS.items():
        base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
        # Kind codes may contain digits; parse rather than treating all digits as the number.
        from patent_pipeline.models import publication_parts
        number = publication_parts(pub)[1]
        def get(url, params=None):
            with urlopen(url + ("?" + urlencode(params) if params else ""), timeout=30) as response:
                return response.read()
        search = json.loads(get(base + "esearch.fcgi", {"db": "protein", "term": number + " AND gbdiv_pat[PROP]",
            "retmode": "json", "retmax": 250, "tool": "ZymoPatentPipeline"}))["esearchresult"]
        time.sleep(.4)
        proteins = get(base + "efetch.fcgi", {"db": "protein", "id": ",".join(search["idlist"]),
            "rettype": "gp", "retmode": "xml", "tool": "ZymoPatentPipeline"})
        (ROOT / f"{pub}.xml.gz").write_bytes(gzip.compress(proteins, mtime=0))
        time.sleep(.4)
        url = f"https://patents.google.com/patent/{pub}/en"
        patent = get(url)
        (ROOT / f"{pub}.html.gz").write_bytes(gzip.compress(patent, mtime=0))
        manifest["targets"][target] = {"publication": pub, "patent_url": url,
            "protein_query": search.get("querytranslation"), "protein_count": int(search["count"]),
            "fixture_record_count": len(search["idlist"])}
        print(target, pub, search["count"], len(patent), flush=True)
        time.sleep(1.1)
    (ROOT / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
