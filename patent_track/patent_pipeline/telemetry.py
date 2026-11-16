from __future__ import annotations

import json
import time
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path


def utcnow() -> str:
    return datetime.now(timezone.utc).isoformat()


class Telemetry:
    def __init__(self, out: Path):
        self.out = out
        out.mkdir(parents=True, exist_ok=True)
        self.started = time.perf_counter()
        self.started_at = utcnow()
        self.stages: list[dict] = []
        self.requests: list[dict] = []
        self.issues: list[dict] = []
        self.current_stage = "setup"

    @contextmanager
    def stage(self, name: str):
        previous = self.current_stage
        self.current_stage = name
        record = {"stage": name, "input_count": 0, "output_count": 0}
        start = time.perf_counter()
        print(f"[{name}] starting", flush=True)
        try:
            yield record
        except Exception:
            record["failed"] = True
            raise
        finally:
            record["seconds"] = round(time.perf_counter() - start, 4)
            self.stages.append(record)
            self.current_stage = previous
            print(f"[{name}] {record['output_count']} outputs in {record['seconds']:.2f}s", flush=True)

    def issue(self, kind: str, message: str, publication: str = ""):
        value = {"kind": kind, "message": message, "publication": publication,
                 "stage": self.current_stage}
        self.issues.append(value)
        print(f"[{kind}] {publication} {message}", flush=True)

    def snapshot(self) -> dict:
        return {"started_at": self.started_at, "finished_at": utcnow(),
                "total_seconds": round(time.perf_counter() - self.started, 4),
                "stages": self.stages, "requests": self.requests, "issues": self.issues}

    def write(self):
        (self.out / "metrics.json").write_text(json.dumps(self.snapshot(), indent=2) + "\n")
