from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .http import HTTP, RetrievalError
from .llm import propose_target_profile
from .models import TargetProfile, publication_id
from .pipeline import Pipeline
from .sources import FIXTURES, FixtureSources, LiveSources
from .telemetry import Telemetry


def positive_int(value):
    number = int(value)
    if number <= 0:
        raise argparse.ArgumentTypeError("must be positive")
    return number


def run(args, target: str, out: Path, seed: list[str], fixture: bool) -> dict:
    telemetry = Telemetry(out)
    with telemetry.stage("target profile") as step:
        step["input_count"] = 1
        profile = TargetProfile(**json.loads(Path(args.profile).read_text())) if args.profile else propose_target_profile(target)
        step["output_count"] = 1
    cache = Path(args.cache)
    source = FixtureSources(telemetry) if fixture else LiveSources(HTTP(cache, telemetry, args.refresh), args.max_proteins)
    pipeline = Pipeline(source, telemetry, cache, args.max_patents, args.max_citations)
    try:
        pipeline.run(profile, seed)
    except (RetrievalError, RuntimeError, ValueError, OSError) as error:
        print(f"{target}: {error}", file=sys.stderr)
    summary_path = out / "summary.json"
    if not summary_path.exists():
        telemetry.write()
        return {"status": "failed", "target": target, "error": "Run failed before export; inspect metrics.json"}
    return json.loads(summary_path.read_text())


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Patent antibody retrieval; ZymoLLM calls are explicit stubs.")
    sub = parser.add_subparsers(dest="command", required=True)
    for name in ("run", "smoke"):
        p = sub.add_parser(name)
        p.add_argument("--out", type=Path, default=Path("runs") / name)
        p.add_argument("--cache", default=".cache")
        p.add_argument("--refresh", action="store_true", help="Re-fetch HTTP responses; keep deterministic domain cache")
        p.add_argument("--profile", help="JSON target profile (run only)")
        p.add_argument("--max-patents", type=positive_int, default=10 if name == "smoke" else None)
        p.add_argument("--max-citations", type=positive_int, default=10 if name == "smoke" else None)
        p.add_argument("--max-proteins", type=positive_int, default=500 if name == "smoke" else None,
                       help="Optional per-search protein cap; reported as incomplete coverage")
        if name == "run":
            p.add_argument("--target", required=True)
            p.add_argument("--publication", action="append", default=[], type=publication_id)
            p.add_argument("--offline", action="store_true", help="Use bounded offline fixtures")
        else:
            p.add_argument("--live", action="store_true", help="Requires EPO OPS credentials and real ANARCI")
    args = parser.parse_args(argv)
    if args.command == "run":
        summary = run(args, args.target, args.out, args.publication, args.offline)
        print(json.dumps(summary, indent=2))
        return 2 if summary["status"] == "blocked" else int(summary["status"] == "failed")
    if args.profile:
        parser.error("--profile is supported by run; smoke always uses the three defined targets")
    manifest = json.loads((FIXTURES / "manifest.json").read_text())
    summaries = []
    for target, item in manifest["targets"].items():
        out = args.out / target.lower().replace("-", "").replace(" ", "")
        summary = run(args, target, out, [item["publication"]], not args.live)
        if summary["status"] not in ("blocked", "failed") and not summary["antibody_domains"]:
            summary["status"] = "failed"
            summary["smoke_error"] = "No antibody domains reached real ANARCI numbering"
            (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
            with (out / "report.md").open("a") as handle:
                handle.write("\nSmoke validation failed: no antibody domains reached real ANARCI numbering.\n")
        summaries.append(summary)
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "smoke_summary.json").write_text(json.dumps(summaries, indent=2) + "\n")
    lines = ["# Three-target smoke results", "", "Mode: " + ("live retrieval" if args.live else "bounded offline fixtures with real ANARCI") + ". ZymoLLM is mocked.", "",
             "| Target | Status | Unique proteins | Distinct candidate domains | Confirmed | Mock examples | Seconds |", "|---|---|---:|---:|---:|---:|---:|"]
    for summary in summaries:
        folder = summary["target"].lower().replace("-", "").replace(" ", "")
        metrics = json.loads((args.out / folder / "metrics.json").read_text())
        lines.append(f"| [{summary['target']}]({folder}/report.md) | {summary['status']} | {summary['unique_sequences']} | {summary['candidate_rows']} | {summary['confirmed_antibodies']} | {summary['mock_demonstration_rows']} | {metrics['total_seconds']:.3f} |")
    lines += ["", "Offline results exercise cached public patent/NCBI examples and format fixtures; they do not validate live OPS discovery. Mock examples are confined to mock_confirmed.tsv. Candidate domains are not confirmed binders. Live runs require EPO_OPS_KEY and EPO_OPS_SECRET.", "",
              "Planning estimates after setup: offline smoke 10–60 seconds; bounded live smoke 3–15 minutes per target, 10–45 minutes total. The table reports actual per-target measurements, including faster cached/offline processing.", ""]
    (args.out / "report.md").write_text("\n".join(lines))
    print(json.dumps(summaries, indent=2))
    return 2 if any(s["status"] == "blocked" for s in summaries) else int(any(s["status"] == "failed" for s in summaries))


if __name__ == "__main__":
    raise SystemExit(main())
