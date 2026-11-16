# Patent antibody retrieval

A small, sequential Python workflow:

`target profile → title/abstract discovery → antibody CPC + target mentions → protein/listing retrieval → ANARCI → one-hop citations → patent text → review → TSV`

**ZymoLLM is currently mocked.** Live `confirmed.tsv` stays empty. Recovered variable
domains are candidates, not confirmed binders. Offline smoke tests additionally
write `mock_confirmed.tsv` with three labelled, source-grounded examples.

## Run

From this directory:

```bash
conda env create --prefix ./.conda --file environment.yml
.conda/bin/python -m pip install -e '.[test]'
.conda/bin/python -m patent_pipeline smoke --out runs/smoke
.conda/bin/python -m pytest -q
```

The environment uses Python
3.11, ANARCI's packaged models, and HMMER 3.4. HMMER 3.3.2 was unavailable for
this ARM Mac. Running via `.conda/bin/python` works without activating Conda.

Live discovery requires free [EPO OPS registration](https://www.epo.org/en/searching-for-patents/data/web-services/ops):

```bash
export EPO_OPS_KEY='your-consumer-key'
export EPO_OPS_SECRET='your-consumer-secret'
export NCBI_EMAIL='your-contact-address'
# Optional: export NCBI_API_KEY='your-key'
.conda/bin/python -m patent_pipeline smoke --live --out runs/live_smoke
.conda/bin/python -m patent_pipeline run --target 'IL-10' --out runs/il10
.conda/bin/python -m patent_pipeline run --target 'GM-CSF' --publication US12291567B2 --out runs/gmcsf
```

Missing OPS credentials cause exit code **2** and a blocked report. Failed runs
return **1**. Retrieval completion with pending LLM review returns **0**; that is
not confirmation of antibody binding. An empty successful search is different
from a failed or blocked request. `.env.example` documents variables; `.env` is
not loaded automatically and credentials are never written to the cache.

Smoke uses at most ten seed candidate publications and ten additional citation
candidate publications per target, conservatively bounding family counts. It also
caps each NCBI search at 500 protein records. Limits and provider truncation are
reported. Normal `run` is uncapped unless `--max-patents`, `--max-citations`, or
`--max-proteins` is supplied. Each patent's protein search is paginated separately
from discovery, so a target search's first 500 hits do not masquerade as all
sequences from that patent.

```bash
# Optional bounded public-source probes; does not require OPS credentials.
.conda/bin/python scripts/probe_sources.py --out runs/source_probe
# A single offline target:
.conda/bin/python -m patent_pipeline run --target 'IL-12p70' --offline --out runs/il12_offline
```

## Target and evidence rules

`propose_target_profile` and `review_patent` in
[llm.py](patent_pipeline/llm.py) are the only two LLM integration points. They have
targeted prompts and predictable data contracts. Replace their bodies with your
ZymoLLM calls; no provider SDK is coupled to the rest of the pipeline.

The profile call returns `TargetProfile` (canonical name, synonyms, confounders,
conditional terms, human target species). Known profiles cover IL-10, GM-CSF,
and IL-12p70. Unknown names use a literal-only stub and explicitly report that
synonym expansion is absent. Pass `--profile profile.json` to supply a profile:

```json
{
  "canonical_name": "IL-10",
  "synonyms": ["IL-10", "IL10", "interleukin-10"],
  "confounders": ["IL-10R", "IL10RA", "IL10RB", "interleukin-10 receptor"],
  "conditional_terms": [],
  "target_species": "human"
}
```

Profiles become deterministic EPO CQL and NCBI Entrez queries, rather than
LLM-generated SQL. CPC `C07K16/00` and its descendants are mandatory. Missing CPC
remains unresolved. IPC and antibody-related therapeutic classes are not silently
substituted for CPC. Confounders suppress overlapping mentions, not whole patents:
an IL-10 patent can legitimately discuss the IL-10 receptor.

Human **target** binding is required; antibody origin need not be human. ANARCI
species predictions are not used to establish target species. IL-12p35/p40 terms
are conditional, and a binder needs evidence of assembled human IL-12p70 binding;
binding IL-23 as well is permitted. A cytokine fused to an unrelated antibody is
not that antibody's target.

The review function receives the target, source-anchored claims/descriptions, and
numbered domains with accession/publication/SEQ ID mappings. Its response contains
`status`, `mocked`, and `assertions`. Each assertion references existing
`vh_domain_id` and/or `vl_domain_id`, `pairing_status`, `human_target_binding`, and
`evidence` entries (`publication`, `section`, `locator`, `quote`). Evidence anchors
and exact quotes are validated before export. Incomplete/ambiguous domains cannot
be confirmed. Family membership alone never maps sequence identifiers. A patent's
list of heavy and light chains never creates a Cartesian product of pairs.

The stub always returns `pending_review`. Saved `review_requests.jsonl` contains
the complete future API payload; `reviews.jsonl` preserves review results.

## Outputs and provenance

Each output directory contains:

| File | Meaning |
|---|---|
| `confirmed.tsv` | Supported antibodies after a real reviewer; currently header only |
| `candidates.tsv` | Deduplicated, unpaired domain candidates, completeness and pending status |
| `mock_confirmed.tsv` | Offline-only demonstration pairs; never live confirmations |
| `evidence.jsonl` | Every domain observation, original protein, accession, SEQ ID, source and coordinates |
| `patents.jsonl` | Metadata and section-anchored full text actually retrieved |
| `filter_decisions.jsonl` | CPC and target-match decisions, including unresolved cases |
| `target_profile.json`, `profile_prompt.txt` | Exact profile and synonym prompt |
| `review_requests.jsonl`, `reviews.jsonl` | Review input and response, with explicit mock status |
| `summary.json`, `metrics.json`, `report.md` | Counts, timings, source gaps and run status |

Final TSV columns are:

```text
canonical_key  vh_seq  vl_seq  cdrh3  target  sources  structures  cluster_id  pairing_status
```

The actual delimiter is a tab. Arrays are compact JSON in cells. Keys hash the
uppercase VH/VL sequence pair, independently of target or source. Missing chains
are empty cells. CDRH3 uses IMGT positions 105–117, including insertions, with
alignment gaps removed. Domain coordinates in `evidence.jsonl` are zero-based,
half-open coordinates into the original accession/listing sequence.

Candidate rows contain one chain each. `unpaired` means pairing is unsupported;
it does not establish a naturally single-domain antibody. `structures=[]` and
`cluster_id=""` are reserved for later work. Exact sequence deduplication retains
every source observation in `evidence.jsonl`. Reported domain observations can
exceed distinct candidate domains because full chains, variable fragments, and
listings can disclose the same domain.

Raw HTTP bodies and retrieval metadata are retained in `.cache/http`. ANARCI
results are cached by sequence, version and numbering options in `.cache/domains`.
`--refresh` re-fetches HTTP responses. Interrupted or failed requests are not
cached as successful empty results. If a cached NCBI history handle expires before
its protein batches are fetched, use `--refresh` to obtain a new handle.

## Source assessment

| Source | Coverage and useful content | Filters, freshness and limitations |
|---|---|---|
| [EPO OPS](https://www.epo.org/en/searching-for-patents/data/web-services/ops) | Worldwide bibliographic metadata, titles/abstracts, CPC, applicants, families and citations; claims/descriptions for supported authorities | Required primary discovery source. English title/abstract CQL + CPC. Free registration, OAuth and a 4 GB/week free threshold. Search pages max 100; max 2,000 retrieved results per query, so oversized searches are date-partitioned. |
| [OPS full-text coverage](https://www.epo.org/en/service-support/faq/searching-patents/open-patent-services/coverage/what-coverage-full-text) | EP/WO and several national collections | US full text is absent from this coverage list. Try English EP/WO family text, retaining that publication's claims and SEQ IDs; then the original Google page. Missing or non-English text remains a gap. |
| [NCBI Protein / GenBank PAT](https://www.ncbi.nlm.nih.gov/education/patent_and_ip_faqs/) | Protein sequences, accession versions, patent title/reference, usually SEQ ID, record creation/update dates; includes patent data beyond the US through partner sources | Entrez text and `gbdiv_pat[PROP]`, plus verified publication-number retrieval. No patent CPC, abstract, complete claims, ownership history or reliable pairing. NCBI documents exclusion of US pre-grant sequences and some algorithm-generated/invalid sequences. Publication-to-index lag varies. |
| [GenBank flat files and daily updates](https://www.ncbi.nlm.nih.gov/genbank/genbank/) | `gbpat*.seq.gz` contains **nucleotide** PAT records. `/ncbi-asn1/` is an alternative serialization, not extra patent text. | Full releases every two months; INSDC exchange and incremental updates daily. No CPC filter. Bulk download/indexing is unnecessary for these bounded targets; CDS annotations do not guarantee translations for every patent sequence. |
| [NCBI BLAST databases](https://ftp.ncbi.nlm.nih.gov/blast/db/) | `pataa` protein and `patnt` nucleotide patent sequence databases; useful for local similarity searches at larger scale | Directory/manifest access can establish available versions, not complete legal/patent coverage. Antigen similarity identifies an antigen-containing patent, not an antigen-binding antibody. Local BLAST still has output/configuration limits and indexing costs. Deferred. |
| [European Publication Server](https://www.epo.org/en/searching-for-patents/data/web-services/publication-server) | EP A/B publications as XML/HTML/PDF/ZIP, and accessible listing links when exposed | Publication-number retrieval; updated weekly. Probe and follow explicitly labelled machine-readable listing links. Main patent text is not assumed to contain a complete sequence listing. |
| [Google Patents](https://support.google.com/faqs/answer/7049585) | Broad international metadata and translated full text, original/current assignees, family links, citations, and some listing links | Known-publication HTML fallback only. Coverage is not guaranteed; automated access and page structure can fail. Primary assignees and CPC are parsed separately from cited patents and IPC. Titles/abstracts remain the discovery focus. |
| [ENA / DDBJ](https://www.ddbj.nig.ac.jp/ddbj/patent-data-e.html) | Partner patent sequence sources; useful for data not present in NCBI, with authority- and sequence-type-specific coverage | Additional adapters could improve sequence recovery but would add another reconciliation path. Assessed here; not implemented in this first version. |
| BigQuery / Lens PatSeq | Potential broader structured discovery/sequence data | Excluded under the stated access constraint. |

The small live probes write measured source availability to
`runs/source_probe/report.md` and `probe.json`; they are separate from the blocked
authenticated OPS smoke test. Actual availability is recorded per publication.
No global source is treated as complete. Original/current assignee fields retain
provider assertions and are not legal ownership determinations.

### Concrete findings

- The supplied `us|7662379` example expands to `us OR 7662379` in Entrez and yielded
  more than 42 million hits in the exploratory probe. A patent-number query yielded
  32 proteins. Country/number/kind verification prevents number collisions.
- A broad IL-10 query returned an IL-10 receptor-2 patent and an eight-residue
  peptide. Target text alone does not identify a full antibody variable domain.
- The IL-10 control contains full chains, variable domains and incomplete germline
  controls. ANARCI filtering and completeness flags are both useful.
- Google HTML and NCBI protein XML successfully supplied the three regression
  controls. One HTTP client received Google 503 responses while the standard
  Python HTTP client served the same public pages; fallback access remains a
  measured source dependency.
- Additional EP2192132A2 records demonstrate non-US patent references. These EP
  references lack SEQ IDs: a publication day such as `02-JUN-2010` must not become
  SEQ ID `02`. Accessions and sequences are retained, with missing identifiers
  reported. A `JP` keyword probe returned US references mentioning Japan, so that
  probe does not validate JP-specific protein coverage.

## Timings, estimates and remaining gaps

Every stage records wall-clock duration and input/output counts. Every HTTP event
records source, retrieval timestamp, bytes, cache use, attempts, network time,
rate-limit wait and retry wait. The report separates cached reads from network
traffic. Source publication/filing/priority dates are distinct from NCBI record
dates and retrieval dates.

Planning estimates were **16–24 engineering hours** for implementation;
**10–60 seconds** for offline smoke; **3–15 minutes per target** for bounded live
smoke, or **10–45 minutes** for all three. Uncapped searches may take tens of
minutes to several hours. These are planning ranges. Actual measurements appear
in each run's report, including faster cached/offline runs.

Retrieval is deliberately sequential: NCBI stays below three requests/second;
OPS searches use a 6.1-second interval, other OPS/family calls at least 1.05
seconds; Google/EPS requests at least 1.05 seconds. Retries are bounded at three
attempts and honor `Retry-After`; deferrals above 60 seconds are reported for a
later rerun instead of retrying early. OPS token expiry is refreshed
once after an authorization failure. These defaults follow the
[OPS restrictions](https://www.epo.org/en/service-support/faq/searching-patents/open-patent-services/general-information/are-there-any)
and [NCBI usage guidance](https://eutilities.github.io/site/API_Key/usageandkey/).

Remaining gaps are explicit: description-only discovery, missing CPC, inaccessible
listings, proteins missing from NCBI, unannotated DNA, scanned sequences, and
unreviewed target/chain mappings. Original-patent recovery tries family proteins
only when no proteins/listings were found and stops at the first useful family
member. Indexed proteins with no antibody domains trigger a further check of the
known patent page for accessible listing links before rejection. Additional
versions can contain additional sequences. Citation expansion
is one hop, not an exhaustive patent graph. There is no OCR, DNA translation,
bulk indexing, structure lookup or similarity clustering in this version.

NCBI data use: [NCBI disclaimer and copyright notice](https://www.ncbi.nlm.nih.gov/About/disclaimer.html).
