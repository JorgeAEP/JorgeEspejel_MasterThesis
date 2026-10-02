#!/usr/bin/env python3
"""Run the RAVEN Bipolaris phosphate-gradient MICOM checks.

This focused panel fills in intermediate phosphate bounds for the Results
phosphate-limitation curve. Bounds 1000 and 10 already exist in the core panel,
so the default run only adds bounds 100 and 30.
"""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from datetime import datetime
from pathlib import Path


THIS_SCRIPT = Path(__file__).resolve()
PACKAGE = THIS_SCRIPT.parents[1]
ROOT = PACKAGE.parent

RUNNER = ROOT / "Bipolaris__sorokiniana/02_py_scripts/306_run_phaseA_kegg_medium_emma_style.py"
PYTHON = ROOT / "Bipolaris__sorokiniana/.venv-report/bin/python"

REPORT_ROOT = PACKAGE / "reports/116_raven_phosphate_gradient_micom"

MANIFESTS = [
    "fungus_alone__raven_independent",
    "hvsc1_plus_bipolaris__raven_independent",
    "hvsc2_plus_bipolaris__raven_independent",
]

SCENARIOS = {
    "phosphate_bound_100": PACKAGE / "media/sensitivity_panel/phosphate_bound_100.csv",
    "phosphate_bound_30": PACKAGE / "media/sensitivity_panel/phosphate_bound_30.csv",
    "phosphate_bound_3": PACKAGE / "media/sensitivity_panel/phosphate_bound_3.csv",
}


def output_summary_path(manifest: str, scenario: str) -> Path:
    return (
        REPORT_ROOT
        / manifest
        / scenario
        / f"{manifest}_{scenario}_emma_style_community_growth_summary.csv"
    )


def run_job(
    manifest: str,
    scenario: str,
    medium_csv: Path,
    solver: str,
    fraction: float,
    debug: bool,
    force: bool,
) -> dict[str, object]:
    summary_path = output_summary_path(manifest, scenario)
    if summary_path.exists() and not force:
        return {
            "manifest": manifest,
            "scenario": scenario,
            "medium_csv": str(medium_csv),
            "returncode": 0,
            "status": "skipped_existing",
            "summary_path": str(summary_path),
            "stdout_log": "",
            "stderr_log": "",
        }

    cmd = [
        str(PYTHON),
        str(RUNNER),
        "--package-root",
        str(PACKAGE),
        "--report-root",
        str(REPORT_ROOT),
        "--manifest-name",
        manifest,
        "--medium-name",
        "full",
        "--medium-csv",
        str(medium_csv),
        "--medium-label",
        scenario,
        "--solver",
        solver,
        "--fraction",
        str(fraction),
    ]
    if debug:
        cmd.append("--debug")

    log_dir = REPORT_ROOT / "_logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    stdout_log = log_dir / f"{manifest}__{scenario}.stdout.log"
    stderr_log = log_dir / f"{manifest}__{scenario}.stderr.log"

    with stdout_log.open("w", encoding="utf-8") as out, stderr_log.open("w", encoding="utf-8") as err:
        completed = subprocess.run(cmd, stdout=out, stderr=err, text=True)

    return {
        "manifest": manifest,
        "scenario": scenario,
        "medium_csv": str(medium_csv),
        "returncode": completed.returncode,
        "status": "ok" if completed.returncode == 0 else "failed",
        "summary_path": str(summary_path),
        "stdout_log": str(stdout_log),
        "stderr_log": str(stderr_log),
        "cmd": " ".join(cmd),
    }


def read_growth_summary(row: dict[str, object]) -> dict[str, object]:
    path = Path(str(row["summary_path"]))
    out = {
        "manifest": row["manifest"],
        "scenario": row["scenario"],
        "medium_csv": row["medium_csv"],
        "run_status": row["status"],
        "returncode": row["returncode"],
        "summary_path": str(path),
    }
    if not path.exists():
        out.update({"status": "missing", "community_growth": None, "objective_value": None})
        return out
    with path.open(newline="", encoding="utf-8") as handle:
        summary = next(csv.DictReader(handle), {})
    out.update(summary)
    return out


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fieldnames = sorted({key for row in rows for key in row})
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--solver", default="cplex")
    parser.add_argument("--fraction", type=float, default=1.0)
    parser.add_argument("--debug", action="store_true")
    parser.add_argument("--force", action="store_true", help="Rerun jobs even if summaries already exist.")
    parser.add_argument("--manifest", action="append", choices=MANIFESTS)
    parser.add_argument("--scenario", action="append", choices=sorted(SCENARIOS))
    args = parser.parse_args()

    REPORT_ROOT.mkdir(parents=True, exist_ok=True)
    if not PYTHON.exists():
        raise FileNotFoundError(PYTHON)
    if not RUNNER.exists():
        raise FileNotFoundError(RUNNER)

    manifests = args.manifest or MANIFESTS
    scenarios = {name: SCENARIOS[name] for name in (args.scenario or SCENARIOS)}
    for manifest in manifests:
        manifest_path = PACKAGE / "manifests" / f"{manifest}.csv"
        if not manifest_path.exists():
            raise FileNotFoundError(manifest_path)
    for medium_csv in scenarios.values():
        if not medium_csv.exists():
            raise FileNotFoundError(medium_csv)

    print("Script 116 RAVEN phosphate-gradient MICOM")
    print(f"started: {datetime.now().isoformat(timespec='seconds')}")
    print(f"package_root: {PACKAGE}")
    print(f"report_root: {REPORT_ROOT}")
    print(f"manifests: {';'.join(manifests)}")
    print(f"scenarios: {';'.join(scenarios)}")

    job_rows = []
    total = len(manifests) * len(scenarios)
    counter = 0
    for scenario, medium_csv in scenarios.items():
        for manifest in manifests:
            counter += 1
            print(f"[{counter}/{total}] {manifest} :: {scenario}", flush=True)
            row = run_job(
                manifest=manifest,
                scenario=scenario,
                medium_csv=medium_csv,
                solver=args.solver,
                fraction=args.fraction,
                debug=args.debug,
                force=args.force,
            )
            job_rows.append(row)
            print(f"  {row['status']} ({row['returncode']})", flush=True)
            write_csv(REPORT_ROOT / "job_status_live.csv", job_rows)

    write_csv(REPORT_ROOT / "job_status.csv", job_rows)
    write_csv(
        REPORT_ROOT / "phosphate_gradient_growth_summary.csv",
        [read_growth_summary(row) for row in job_rows],
    )

    failed_jobs = sum(1 for row in job_rows if row["returncode"] != 0)

    summary = {
        "package_root": str(PACKAGE),
        "report_root": str(REPORT_ROOT),
        "solver": args.solver,
        "fraction": args.fraction,
        "manifests": manifests,
        "scenarios": list(scenarios),
        "job_count": int(len(job_rows)),
        "failed_jobs": failed_jobs,
        "all_jobs_ok": failed_jobs == 0,
        "finished": datetime.now().isoformat(timespec="seconds"),
    }
    (REPORT_ROOT / "run_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print("Script 116 compact summary")
    print(f"job_count: {summary['job_count']}")
    print(f"failed_jobs: {summary['failed_jobs']}")
    print(f"all_jobs_ok: {summary['all_jobs_ok']}")
    print(f"growth_summary: {REPORT_ROOT / 'phosphate_gradient_growth_summary.csv'}")
    print(f"outputs: {REPORT_ROOT}")

    if not summary["all_jobs_ok"]:
        sys.exit(1)


if __name__ == "__main__":
    main()
