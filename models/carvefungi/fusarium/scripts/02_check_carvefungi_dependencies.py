#!/usr/bin/env python3
"""Check the validated shared CarveFungi toolchain without changing it."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import os
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PRECEDENT = ROOT.parent / "Bipolaris_sorokiniana_carveme"
SHARED_ENV = PRECEDENT / "02_environment"
INPUT_FASTA = ROOT / "01_Input/Fusarium_graminearum_proteins.faa"
REPORT_DIR = ROOT / "06_reports/02_carvefungi_dependency_check"
REPORT_PATH = REPORT_DIR / "dependency_check.json"

VENV = SHARED_ENV / ".venv-carveme"
PYTHON = VENV / "bin/python"
CARVE = VENV / "bin/carve"
DIAMOND = SHARED_ENV / "tools/diamond_bin/usr/bin/diamond"
EGGNOG = SHARED_ENV / "eggnog-mapper/emapper.py"
EGGNOG_DATA = SHARED_ENV / "eggnog-mapper/data"
PSIPRED_RUNSEC = SHARED_ENV / "psipred/runSec.sh"
CARVEFUNGI_SOURCE = SHARED_ENV / "CarveFungi"
CREATE_MODEL = CARVEFUNGI_SOURCE / "bin/CreateModelEggNogPool.py"
COMPARTMENT_PREDICT = CARVEFUNGI_SOURCE / "bin/compartmentPrediction/predict.py"
DEEP_MODELS = PRECEDENT / "01_Input/Models_trained/models-trained"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def module_check() -> dict[str, bool]:
    modules = ["Bio", "cobra", "carveme", "keras", "numpy", "pandas", "tensorflow"]
    if not PYTHON.exists():
        return {module: False for module in modules}
    code = (
        "import importlib.util,json; "
        f"m={modules!r}; print(json.dumps({{x: bool(importlib.util.find_spec(x)) for x in m}}))"
    )
    completed = subprocess.run([str(PYTHON), "-c", code], capture_output=True, text=True)
    if completed.returncode:
        return {module: False for module in modules}
    return json.loads(completed.stdout)


def collect_status() -> dict[str, object]:
    eggnog_databases = [
        EGGNOG_DATA / "eggnog.db",
        EGGNOG_DATA / "eggnog.taxa.db",
        EGGNOG_DATA / "eggnog_proteins.dmnd",
    ]
    required = {
        "input_fasta": INPUT_FASTA.is_file(),
        "shared_venv": VENV.is_dir(),
        "shared_python": PYTHON.is_file(),
        "carve_command": CARVE.is_file(),
        "diamond": DIAMOND.is_file(),
        "eggnog_mapper": EGGNOG.is_file(),
        "eggnog_databases": all(path.is_file() and path.stat().st_size > 0 for path in eggnog_databases),
        "psipred_runSec": PSIPRED_RUNSEC.is_file(),
        "carvefungi_source": CARVEFUNGI_SOURCE.is_dir(),
        "CreateModelEggNogPool": CREATE_MODEL.is_file(),
        "compartment_prediction": COMPARTMENT_PREDICT.is_file(),
        "compartment_models": DEEP_MODELS.is_dir(),
    }
    modules = module_check()
    disk = shutil.disk_usage(ROOT)
    return {
        "checked_at_utc": datetime.now(timezone.utc).isoformat(),
        "project_root": str(ROOT),
        "precedent_project": str(PRECEDENT),
        "shared_environment_read_only_intent": True,
        "required": required,
        "python_modules": modules,
        "input_fasta_sha256": sha256(INPUT_FASTA) if INPUT_FASTA.is_file() else None,
        "tool_provenance_sha256": {
            "emapper.py": sha256(EGGNOG) if EGGNOG.is_file() else None,
            "CreateModelEggNogPool.py": sha256(CREATE_MODEL) if CREATE_MODEL.is_file() else None,
            "compartment_predict.py": sha256(COMPARTMENT_PREDICT)
            if COMPARTMENT_PREDICT.is_file() else None,
        },
        "free_disk_bytes": disk.free,
        "ready": all(required.values()) and all(modules.values()),
        "environment": {
            "EGGNOG_DATA_DIR": str(EGGNOG_DATA),
            "EGGNOG_CPUS_default": os.environ.get("EGGNOG_CPUS", "4"),
            "COMPARTMENT_BATCH_SIZE_default": os.environ.get("COMPARTMENT_BATCH_SIZE", "100"),
        },
    }


def main() -> int:
    status = collect_status()
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_PATH.write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
    print(f"Input FASTA SHA256: {status['input_fasta_sha256']}")
    print(f"Required paths:     {sum(status['required'].values())}/{len(status['required'])}")
    print(f"Python modules:     {sum(status['python_modules'].values())}/{len(status['python_modules'])}")
    print(f"Reconstruction ready: {'YES' if status['ready'] else 'NO'}")
    print(f"Report: {REPORT_PATH}")
    return 0 if status["ready"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
