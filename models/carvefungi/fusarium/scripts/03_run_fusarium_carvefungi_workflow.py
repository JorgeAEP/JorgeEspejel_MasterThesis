#!/usr/bin/env python3
"""Run the Fusarium CarveFungi workflow; default behavior is check-only."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import importlib.util


ROOT = Path(__file__).resolve().parents[1]
CHECKER_PATH = Path(__file__).with_name("02_check_carvefungi_dependencies.py")
spec = importlib.util.spec_from_file_location("dependency_check", CHECKER_PATH)
if spec is None or spec.loader is None:
    raise RuntimeError(f"Cannot load dependency checker: {CHECKER_PATH}")
dependency_check = importlib.util.module_from_spec(spec)
spec.loader.exec_module(dependency_check)

INPUT_FASTA = ROOT / "01_Input/Fusarium_graminearum_proteins.faa"
FUNCTIONAL = ROOT / "01_Input/Fusarium_graminearum_proteins.faa.emapper.annotations"
FUNCTIONAL_CF = ROOT / "01_Input/Fusarium_graminearum_proteins.carvefungi_compat.annotations"
LOCALIZATION = Path(str(FUNCTIONAL_CF) + ".loc_pred")
PSIPRED_INPUTS = ROOT / "01_Input/psipred_inputs"
RUNTIME = ROOT / "02_environment/CarveFungi_runtime"
NOSPACE_BIN = Path("/tmp/fusarium_carvefungi_bin")
MODEL_OUT = ROOT / "04_models/Fusarium_graminearum_carvefungi.xml"
REPORT_DIR = ROOT / "06_reports/03_carvefungi_reconstruction"
CHECK_REPORT = REPORT_DIR / "workflow_readiness.json"
EGGNOG_CPUS = os.environ.get("EGGNOG_CPUS", "4")
BATCH_SIZE = int(os.environ.get("COMPARTMENT_BATCH_SIZE", "100"))
PSIPRED_MAX_LENGTH = 10_000
CPLEX_LIB = Path(os.environ.get(
    "CPLEX_STUDIO_LIB",
    "/home/tlu/IBM/ILOG/CPLEX_Studio222/cplex/bin/x86-64_linux/libcplex2220.so",
))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def ensure_nospace_tools() -> None:
    """Expose tools through paths safe for eggNOG's shell command builder."""
    NOSPACE_BIN.mkdir(parents=True, exist_ok=True)
    diamond_link = NOSPACE_BIN / "diamond"
    expected_target = dependency_check.DIAMOND.resolve()
    if diamond_link.is_symlink() and diamond_link.resolve() != expected_target:
        diamond_link.unlink()
    if not diamond_link.exists():
        diamond_link.symlink_to(expected_target)
    if diamond_link.resolve() != expected_target:
        raise RuntimeError(f"Unexpected DIAMOND shim target: {diamond_link}")


def run_command(command: list[str], log, cwd: Path = ROOT) -> None:
    print("COMMAND: " + " ".join(command), file=log, flush=True)
    environment = os.environ.copy()
    environment["EGGNOG_DATA_DIR"] = str(dependency_check.EGGNOG_DATA)
    environment.setdefault("TF_USE_LEGACY_KERAS", "1")
    ensure_nospace_tools()
    environment["PATH"] = str(NOSPACE_BIN) + os.pathsep + environment.get("PATH", "")
    if CPLEX_LIB.is_file():
        existing = environment.get("LD_PRELOAD", "")
        environment["LD_PRELOAD"] = str(CPLEX_LIB) + (os.pathsep + existing if existing else "")
    completed = subprocess.run(
        command, cwd=cwd, env=environment, stdout=log, stderr=subprocess.STDOUT, text=True
    )
    if completed.returncode:
        raise RuntimeError(f"Command failed ({completed.returncode}): {' '.join(command)}")


def read_single_fasta(path: Path) -> tuple[str, str]:
    header: str | None = None
    sequence_parts: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    raise RuntimeError(f"Expected one FASTA record in {path}")
                header = line[1:]
            else:
                if header is None:
                    raise RuntimeError(f"Sequence found before FASTA header in {path}")
                sequence_parts.append(line)
    if header is None or not sequence_parts:
        raise RuntimeError(f"Missing FASTA record or sequence in {path}")
    return header, "".join(sequence_parts)


def write_long_sequence_fallback_report(records: list[dict[str, object]]) -> None:
    report_path = REPORT_DIR / "psipred_long_sequence_fallbacks.json"
    report_path.write_text(
        json.dumps(
            {
                "method": (
                    "For proteins above PSIPRED's 10,000-residue input limit, run PSIPRED on "
                    "the same 500 edge residues retained by CarveFungi localization "
                    "(original residues 2-251 plus the final 250), then expand the secondary-"
                    "structure prediction to original length before CarveFungi applies its slice."
                ),
                "carvefungi_long_protein_slice": "sequence[1:251] + sequence[-250:]",
                "records": records,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )


def run_psipred(fasta: Path, log, fallback_records: list[dict[str, object]]) -> None:
    """Run PSIPRED, preserving CarveFungi's edge representation for huge proteins."""
    horizontal = fasta.with_suffix(".horiz")
    if horizontal.is_file() and horizontal.stat().st_size:
        return
    header, sequence = read_single_fasta(fasta)
    if len(sequence) <= PSIPRED_MAX_LENGTH:
        run_command([str(dependency_check.PSIPRED_RUNSEC), str(fasta.with_suffix(""))], log)
        return

    # predict.py later reduces every long protein to these exact positions.
    edge_sequence = sequence[1:251] + sequence[-250:]
    if len(edge_sequence) != 500:
        raise RuntimeError(f"Unexpected edge representation length for {fasta}: {len(edge_sequence)}")
    with tempfile.TemporaryDirectory(prefix="fusarium_psipred_edges_") as temporary_directory:
        edge_fasta = Path(temporary_directory) / f"{fasta.stem}.fasta"
        edge_fasta.write_text(f">{header}\n{edge_sequence}\n", encoding="utf-8")
        run_command(
            [str(dependency_check.PSIPRED_RUNSEC), str(edge_fasta.with_suffix(""))],
            log,
        )
        edge_horizontal = edge_fasta.with_suffix(".horiz")
        predicted = "".join(
            line.split("Pred:", 1)[1].replace(" ", "").strip()
            for line in edge_horizontal.read_text(encoding="utf-8").splitlines()
            if "Pred:" in line
        )
    if len(predicted) != 500:
        raise RuntimeError(
            f"Expected 500 PSIPRED edge predictions for {fasta.stem}, found {len(predicted)}"
        )

    # CarveFungi allocates the original length, then retains [1:251] and [-250:].
    expanded = ["C"] * len(sequence)
    expanded[1:251] = predicted[:250]
    expanded[-250:] = predicted[250:]
    horizontal.write_text("Pred: " + "".join(expanded) + "\n", encoding="utf-8")
    record = {
        "protein_id": fasta.stem,
        "original_length": len(sequence),
        "psipred_edge_input_length": len(edge_sequence),
        "expanded_prediction_length": len(expanded),
        "original_fasta_sha256": sha256(fasta),
        "horizontal_prediction": str(horizontal),
    }
    fallback_records.append(record)
    write_long_sequence_fallback_report(fallback_records)
    print(
        f"PSIPRED long-sequence fallback: {fasta.stem} ({len(sequence)} aa)",
        file=log,
        flush=True,
    )


def prepare_isolated_runtime() -> None:
    if not RUNTIME.exists():
        shutil.copytree(
            dependency_check.CARVEFUNGI_SOURCE,
            RUNTIME,
            ignore=shutil.ignore_patterns(".git", "results", "annotations", "__pycache__"),
        )
    for directory in (
        RUNTIME / "results",
        RUNTIME / "data/annotations/fasta",
        RUNTIME / "data/annotations/functional",
        RUNTIME / "data/annotations/localization",
    ):
        directory.mkdir(parents=True, exist_ok=True)
    selected_models = RUNTIME / "bin/compartmentPrediction/selected_models"
    if not selected_models.exists():
        selected_models.symlink_to(dependency_check.DEEP_MODELS, target_is_directory=True)


def write_compatibility_table(log) -> None:
    code = f"""
from pathlib import Path
import pandas as pd
src=Path({str(FUNCTIONAL)!r}); dst=Path({str(FUNCTIONAL_CF)!r})
with src.open() as fh:
    header_line=next(i for i,line in enumerate(fh) if line.startswith('#query'))
df=pd.read_csv(src, skiprows=header_line, sep='\\t')
query_col='#query_name' if '#query_name' in df.columns else '#query'
df=df[df['EC'].notna() & df['EC'].astype(str).ne('-')].copy()
out=pd.DataFrame({{
0:df[query_col].astype(str), 1:df.get('seed_ortholog','-'), 2:df.get('evalue',''),
3:df.get('score',0), 4:df.get('eggNOG_OGs',''), 5:df.get('tax_scope',''),
6:df.get('COG_category',''),
7:df['EC'].astype(str).str.replace('ec:','',regex=False).str.replace(';',',',regex=False)}})
with dst.open('w') as fh:
    fh.write('## CarveFungi-compatible eggNOG table generated from eggNOG-mapper\\n')
    fh.write('## Columns: query seed evalue score OGs tax COG EC\\n##\\n##\\n')
out.to_csv(dst, sep='\\t', header=False, index=False, mode='a')
print('EC-annotated proteins:',len(out))
"""
    run_command([str(dependency_check.PYTHON), "-c", code], log)


def split_fastas(log) -> None:
    code = f"""
from Bio import SeqIO
import pandas as pd
from pathlib import Path
src=Path({str(FUNCTIONAL_CF)!r}); out=Path({str(PSIPRED_INPUTS)!r}); out.mkdir(parents=True,exist_ok=True)
df=pd.read_csv(src,skiprows=4,sep='\\t',header=None); ids=set(df[0].astype(str))
n=0
for record in SeqIO.parse({str(INPUT_FASTA)!r},'fasta'):
    if record.id in ids:
        SeqIO.write(record,out/(record.id+'.fasta'),'fasta'); n+=1
print('PSIPRED FASTAs:',n)
"""
    run_command([str(dependency_check.PYTHON), "-c", code], log)


def predict_localization(log) -> None:
    functional_dir = RUNTIME / "data/annotations/functional"
    predictor_dir = RUNTIME / "bin/compartmentPrediction"
    code = f"""
from pathlib import Path
import pandas as pd
src=Path({str(FUNCTIONAL_CF)!r}); out=Path({str(functional_dir)!r}); size={BATCH_SIZE}
df=pd.read_csv(src,skiprows=4,sep='\\t',header=None)
for i in range(0,len(df),size):
    p=out/f'fusarium.batch{{i//size:04d}}.annotations'
    with p.open('w') as fh: fh.write('## batch\\n##\\n##\\n##\\n')
    df.iloc[i:i+size].to_csv(p,sep='\\t',header=False,index=False,mode='a')
    print(p.name)
"""
    completed = subprocess.run(
        [str(dependency_check.PYTHON), "-c", code], capture_output=True, text=True
    )
    if completed.returncode:
        raise RuntimeError(completed.stderr)
    batch_names = [name for name in completed.stdout.splitlines() if name.strip()]
    predictions: list[Path] = []
    for name in batch_names:
        prediction = FUNCTIONAL_CF.parent / f"{name}.loc_pred"
        predictions.append(prediction)
        if not prediction.is_file() or not prediction.stat().st_size:
            predict_code = (
                "import sys; sys.path.insert(0,{!r}); import predict; "
                "predict.predict({!r},{!r},{!r})"
            ).format(str(predictor_dir), str(PSIPRED_INPUTS) + "/", name, str(FUNCTIONAL_CF.parent) + "/")
            run_command([str(dependency_check.PYTHON), "-c", predict_code], log, predictor_dir)
    merge_code = f"""
from pathlib import Path
import pandas as pd
parts=[Path(x) for x in {repr([str(path) for path in predictions])}]
pd.concat([pd.read_csv(p) for p in parts],ignore_index=True).to_csv(Path({str(LOCALIZATION)!r}),index=False)
for p in parts: p.unlink(missing_ok=True)
"""
    run_command([str(dependency_check.PYTHON), "-c", merge_code], log)


def workflow_readiness() -> dict[str, object]:
    status = dependency_check.collect_status()
    return {
        **status,
        "mode": "check-only",
        "isolated_runtime": str(RUNTIME),
        "no_space_tool_directory": str(NOSPACE_BIN),
        "planned_model_output": str(MODEL_OUT),
        "biomass_provenance_warning": (
            "CarveFungi CreateModelEggNogPool injects its built-in artificial fungal biomass; "
            "this is reconstruction provenance, not a validated Fusarium biomass assertion."
        ),
        "existing_model_will_not_be_overwritten": not MODEL_OUT.exists(),
    }


def execute() -> None:
    readiness = workflow_readiness()
    if not readiness["ready"]:
        raise SystemExit("Dependencies are not ready; inspect workflow_readiness.json")
    if MODEL_OUT.exists():
        raise SystemExit(f"Refusing to overwrite existing model: {MODEL_OUT}")
    input_before = sha256(INPUT_FASTA)
    prepare_isolated_runtime()
    fallback_records: list[dict[str, object]] = []
    log_path = ROOT / "05_logs" / f"fusarium_carvefungi_{datetime.now():%Y-%m-%d_%H%M%S}.log"
    with log_path.open("w", encoding="utf-8") as log:
        if not FUNCTIONAL.exists():
            run_command([
                str(dependency_check.PYTHON), str(dependency_check.EGGNOG),
                "-i", str(INPUT_FASTA), "--itype", "proteins", "-m", "diamond",
                "--tax_scope", "Fungi", "--output", "Fusarium_graminearum_proteins.faa",
                "--output_dir", str(INPUT_FASTA.parent), "--override", "--cpu", EGGNOG_CPUS,
                "--annot_batch_size", "500",
            ], log)
        write_compatibility_table(log)
        shutil.copy2(FUNCTIONAL_CF, RUNTIME / "data/annotations/functional" / FUNCTIONAL_CF.name)
        shutil.copy2(INPUT_FASTA, RUNTIME / "data/annotations/fasta" / INPUT_FASTA.name)
        split_fastas(log)
        for fasta in sorted(PSIPRED_INPUTS.glob("*.fasta")):
            run_psipred(fasta, log, fallback_records)
        if not LOCALIZATION.exists():
            predict_localization(log)
        shutil.copy2(LOCALIZATION, RUNTIME / "data/annotations/localization" / LOCALIZATION.name)
        create_model = RUNTIME / "bin/CreateModelEggNogPool.py"
        code = (
            "import sys; sys.path.insert(0,{!r}); import CreateModelEggNogPool; "
            "CreateModelEggNogPool.carveFungi_pipeline('Fusarium_graminearum',{!r},{!r})"
        ).format(str(create_model.parent), str(FUNCTIONAL_CF), str(LOCALIZATION))
        run_command([str(dependency_check.PYTHON), "-c", code], log, create_model.parent)
        candidates = sorted((RUNTIME / "results").glob("Fusarium_graminearum*.sbml"))
        if not candidates:
            raise RuntimeError("CarveFungi completed without a Fusarium SBML candidate")
        shutil.copy2(candidates[-1], MODEL_OUT)
    input_after = sha256(INPUT_FASTA)
    if input_before != input_after:
        raise RuntimeError("Input FASTA checksum changed during reconstruction")
    result = {
        "completed_at_utc": datetime.now(timezone.utc).isoformat(),
        "input_sha256_before": input_before,
        "input_sha256_after": input_after,
        "input_unchanged": input_before == input_after,
        "mother_model": str(MODEL_OUT),
        "mother_model_sha256": sha256(MODEL_OUT),
        "log": str(log_path),
    }
    (REPORT_DIR / "reconstruction_result.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", action="store_true", help="execute the expensive reconstruction")
    args = parser.parse_args()
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    readiness = workflow_readiness()
    CHECK_REPORT.write_text(json.dumps(readiness, indent=2) + "\n", encoding="utf-8")
    print(f"Dependencies ready: {'YES' if readiness['ready'] else 'NO'}")
    print(f"Check report: {CHECK_REPORT}")
    if args.run:
        execute()
    else:
        print("Check-only mode. Use --run to start reconstruction.")
    return 0 if readiness["ready"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
