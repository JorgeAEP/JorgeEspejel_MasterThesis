#!/usr/bin/env python3
"""Reproducible Bipolaris sorokiniana CarveFungi/CarveMe workflow.

Default mode is a dependency/input check only. Use --run to execute the
annotation, localization and reconstruction steps after dependencies are ready.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path


ROOT = Path("/home/tlu/Desktop/Jorge/VS CODE/Bipolaris_sorokiniana_carveme")
CURATED = Path("/home/tlu/Desktop/Jorge/VS CODE/Bipolaris__sorokiniana")
PROJECT_PY = Path("/home/tlu/Desktop/Jorge/VS CODE/.venv/bin/python")
CARVEME_VENV = ROOT / "02_environment/.venv-carveme"
INPUT_FASTA = ROOT / "01_Input/protein.faa"
PREP_FASTA = ROOT / "01_Input/Bipolaris_sorokiniana_proteins.faa"
FUNCTIONAL = ROOT / "01_Input/Bipolaris_sorokiniana_proteins.faa.emapper.annotations"
FUNCTIONAL_CF = ROOT / "01_Input/Bipolaris_sorokiniana_proteins.carvefungi_compat.annotations"
LOCALIZATION = FUNCTIONAL_CF.parent / (FUNCTIONAL_CF.name + ".loc_pred")
MODEL_OUT = ROOT / "04_models/Bipolaris_sorokiniana_carvefungi.xml"
REPORT = ROOT / "06_reports/Bipolaris_carvefungi_workflow_report.json"
BIOMASS_CSV = ROOT / "06_reports/curated_Bipolaris_biomass_stoichiometry.csv"
LOG = ROOT / "05_logs" / f"bipolaris_carvefungi_{datetime.now():%Y-%m-%d_%H%M%S}.log"
NOSPACE_BIN = Path("/tmp/bipolaris_carvefungi_bin")
PSIPRED_INPUTS = ROOT / "01_Input/psipred_inputs"
EGGNOG_CPUS = os.environ.get("EGGNOG_CPUS", "4")
COMPARTMENT_BATCH_SIZE = int(os.environ.get("COMPARTMENT_BATCH_SIZE", "100"))
CPLEX_STUDIO_LIB = Path(os.environ.get(
    "CPLEX_STUDIO_LIB",
    "/home/tlu/IBM/ILOG/CPLEX_Studio222/cplex/bin/x86-64_linux/libcplex2220.so",
))

# CarveFungi tutorial-specific files. Set these with command-line options or env vars.
CARVEFUNGI_BIN = Path(os.environ.get("CARVEFUNGI_BIN", ROOT / "02_environment/CarveFungi/bin"))
EGGNOG = Path(os.environ.get("EMAPPER", ROOT / "02_environment/eggnog-mapper/emapper.py"))
EGGNOG_DATA = Path(os.environ.get("EGGNOG_DATA_DIR", ROOT / "02_environment/eggnog-mapper/data"))
PSIPRED_RUNSEC = Path(os.environ.get("PSIPRED_RUNSEC", ROOT / "02_environment/psipred/runSec.sh"))
COMPARTMENT_PRED = Path(os.environ.get("COMPARTMENT_PRED", CARVEFUNGI_BIN / "compartmentPrediction/predict.py"))
DEEP_MODELS = Path(os.environ.get("DEEP_MODELS", ROOT / "01_Input/Models_trained/models-trained"))
CREATE_MODEL = Path(os.environ.get("CREATE_MODEL", CARVEFUNGI_BIN / "CreateModelEggNogPool.py"))


def exists_path(path: Path) -> bool:
    return str(path) not in {"", "."} and path.exists()


def run(cmd: list[str], log, cwd: Path = ROOT) -> None:
    print(" ".join(cmd), file=log, flush=True)
    env = os.environ.copy()
    env.setdefault("TF_USE_LEGACY_KERAS", "1")
    env.setdefault("EGGNOG_DATA_DIR", str(EGGNOG_DATA))
    if CPLEX_STUDIO_LIB.exists():
        prior = env.get("LD_PRELOAD", "")
        env["LD_PRELOAD"] = str(CPLEX_STUDIO_LIB) + (os.pathsep + prior if prior else "")
    env["PATH"] = str(NOSPACE_BIN) + os.pathsep + env.get("PATH", "")
    p = subprocess.run(cmd, cwd=cwd, env=env, text=True, stdout=log, stderr=subprocess.STDOUT)
    if p.returncode:
        raise SystemExit(f"Command failed ({p.returncode}): {' '.join(cmd)}")


def latest_curated_model() -> Path | None:
    models = sorted((CURATED / "03_intermediate_models").glob("*/*_post*.xml"))
    if not models:
        models = sorted((CURATED / "03_intermediate_models").glob("*/*.xml"))
    return models[-1] if models else None


def extract_biomass_reference() -> dict:
    """Use the real Bipolaris biomass from the curated model, not the tutorial yeast biomass."""
    model_path = latest_curated_model()
    if not model_path or not PROJECT_PY.exists():
        return {"ok": False, "reason": "missing curated model or project Python", "model": str(model_path)}
    code = f"""
import cobra, csv
model = cobra.io.read_sbml_model({str(model_path)!r})
rid = 'BIOMASS_BIPOLARIS'
if not model.reactions.has_id(rid):
    raise SystemExit('BIOMASS_BIPOLARIS not found')
r = model.reactions.get_by_id(rid)
with open({str(BIOMASS_CSV)!r}, 'w', newline='') as fh:
    w = csv.writer(fh)
    w.writerow(['reaction_id','reaction_name','metabolite_id','metabolite_name','coefficient','formula','charge','compartment'])
    for m,c in sorted(r.metabolites.items(), key=lambda x: x[0].id):
        w.writerow([r.id, r.name, m.id, m.name, c, m.formula, m.charge, m.compartment])
print(len(r.metabolites))
"""
    p = subprocess.run([str(PROJECT_PY), "-c", code], text=True, capture_output=True)
    return {
        "ok": p.returncode == 0,
        "model": str(model_path),
        "biomass_csv": str(BIOMASS_CSV),
        "biomass_metabolites": p.stdout.strip(),
        "error": p.stderr.strip(),
    }


def check() -> dict:
    ensure_nospace_tools()
    diamond = ROOT / "02_environment/tools/diamond_bin/usr/bin/diamond"
    diamond_nospace = NOSPACE_BIN / "diamond"
    carve = CARVEME_VENV / "bin/carve"
    py = CARVEME_VENV / "bin/python"
    eggnog_db_files = list(EGGNOG_DATA.glob("*.dmnd")) + list(EGGNOG_DATA.glob("*.db")) + list(EGGNOG_DATA.glob("*.sqlite"))
    required = {
        "input_fasta": INPUT_FASTA.exists(),
        "carveme_venv": CARVEME_VENV.exists(),
        "carve_command": carve.exists(),
        "diamond": diamond_nospace.exists() or diamond.exists() or bool(shutil.which("diamond")),
        "eggnog_mapper": exists_path(EGGNOG),
        "eggnog_database_files": EGGNOG_DATA.exists() and bool(eggnog_db_files),
        "psipred_runSec": exists_path(PSIPRED_RUNSEC),
        "compartment_predict_py": exists_path(COMPARTMENT_PRED),
        "compartment_deep_models": exists_path(DEEP_MODELS),
        "CreateModelEggNogPool_py": exists_path(CREATE_MODEL),
    }
    py_mods = {}
    if py.exists():
        code = "import importlib.util; mods=['pandas','numpy','Bio','tensorflow','keras','carveme']; print({m: bool(importlib.util.find_spec(m)) for m in mods})"
        p = subprocess.run([str(py), "-c", code], text=True, capture_output=True)
        py_mods["raw"] = p.stdout.strip()
    biomass = extract_biomass_reference()
    return {"required": required, "carveme_python_modules": py_mods, "biomass_reference": biomass}


def ensure_nospace_tools() -> None:
    NOSPACE_BIN.mkdir(parents=True, exist_ok=True)
    diamond_src = ROOT / "02_environment/tools/diamond_bin/usr/bin/diamond"
    diamond_dst = NOSPACE_BIN / "diamond"
    if diamond_src.exists() and not diamond_dst.exists():
        diamond_dst.symlink_to(diamond_src)


def prepare_fasta() -> None:
    if not INPUT_FASTA.exists():
        raise SystemExit(f"Missing input FASTA: {INPUT_FASTA}")
    if INPUT_FASTA.resolve() != PREP_FASTA.resolve():
        shutil.copyfile(INPUT_FASTA, PREP_FASTA)


def split_annotated_fastas(log) -> None:
    code = f"""
from Bio import SeqIO
import pandas as pd
from pathlib import Path
ann = Path({str(FUNCTIONAL)!r})
with ann.open() as fh:
    header_line = next(i for i, line in enumerate(fh) if line.startswith('#query'))
functional = pd.read_csv(ann, skiprows=header_line, sep='\\t')
query_col = '#query_name' if '#query_name' in functional.columns else '#query'
functional = functional[functional['EC'].notna() & (functional['EC'].astype(str) != '-')]
ids = set(functional[query_col].astype(str))
out = Path({str(PSIPRED_INPUTS)!r})
out.mkdir(parents=True, exist_ok=True)
count = 0
for record in SeqIO.parse({str(PREP_FASTA)!r}, 'fasta'):
    if record.id in ids:
        SeqIO.write(record, out / (record.id + '.fasta'), 'fasta')
        count += 1
print(count)
"""
    run([str(CARVEME_VENV / "bin/python"), "-c", code], log)


def write_carvefungi_eggnog_compat(log) -> None:
    """CarveFungi expects older eggNOG columns with EC numbers at column index 7."""
    code = f"""
from pathlib import Path
import pandas as pd
src = Path({str(FUNCTIONAL)!r})
dst = Path({str(FUNCTIONAL_CF)!r})
with src.open() as fh:
    header_line = next(i for i, line in enumerate(fh) if line.startswith('#query'))
df = pd.read_csv(src, skiprows=header_line, sep='\\t')
query_col = '#query_name' if '#query_name' in df.columns else '#query'
ec = df['EC'].fillna('-').astype(str)
df = df[ec.ne('-')]
clean_ec = df['EC'].astype(str).str.replace('ec:', '', regex=False).str.replace(';', ',', regex=False)
out = pd.DataFrame({{
    0: df[query_col].astype(str),
    1: df.get('seed_ortholog', '-'),
    2: df.get('evalue', ''),
    3: df.get('score', 0),
    4: df.get('eggNOG_OGs', ''),
    5: df.get('tax_ceiling', ''),
    6: df.get('COG_category', ''),
    7: clean_ec,
}})
with dst.open('w') as fh:
    fh.write('## CarveFungi-compatible eggNOG table generated from eggNOG v3\\n')
    fh.write('## Columns: query seed evalue score OGs tax COG EC\\n')
    fh.write('## Source: ' + str(src) + '\\n')
    fh.write('##\\n')
out.to_csv(dst, sep='\\t', header=False, index=False, mode='a')
print(len(out))
"""
    run([str(CARVEME_VENV / "bin/python"), "-c", code], log)


def ensure_compartment_models_link() -> None:
    selected = COMPARTMENT_PRED.parent / "selected_models"
    if not selected.exists():
        selected.symlink_to(DEEP_MODELS)


def run_compartment_prediction_batched(log, cf_func: Path) -> None:
    """Avoid OOM by running CarveFungi localization prediction in small chunks."""
    code = f"""
from pathlib import Path
import pandas as pd
src = Path({str(FUNCTIONAL_CF)!r})
cf_func = Path({str(cf_func)!r})
out_file = Path({str(LOCALIZATION)!r})
batch_size = {COMPARTMENT_BATCH_SIZE}
df = pd.read_csv(src, skiprows=4, sep='\\t', header=None)
batch_names = []
for i in range(0, len(df), batch_size):
    batch = df.iloc[i:i+batch_size]
    name = f"{{src.stem}}.batch{{i//batch_size:04d}}.annotations"
    path = cf_func / name
    with path.open('w') as fh:
        fh.write('## CarveFungi-compatible eggNOG batch\\n##\\n##\\n##\\n')
    batch.to_csv(path, sep='\\t', header=False, index=False, mode='a')
    batch_names.append(name)
print('\\n'.join(batch_names))
"""
    p = subprocess.run([str(CARVEME_VENV / "bin/python"), "-c", code], text=True, capture_output=True)
    if p.returncode:
        print(p.stdout, file=log)
        print(p.stderr, file=log)
        raise SystemExit("Could not prepare localization batches")
    batch_names = [x for x in p.stdout.splitlines() if x.strip()]
    pred_files = []
    for name in batch_names:
        pred = LOCALIZATION.parent / (name + ".loc_pred")
        pred_files.append(pred)
        if not pred.exists() or not pred.stat().st_size:
            code = (
                "import sys; sys.path.insert(0, {!r}); import predict; "
                "predict.predict({!r}, {!r}, {!r})"
            ).format(
                str(COMPARTMENT_PRED.parent),
                str(PSIPRED_INPUTS) + "/",
                name,
                str(LOCALIZATION.parent) + "/",
            )
            run([str(CARVEME_VENV / "bin/python"), "-c", code], log, cwd=COMPARTMENT_PRED.parent)
    code = f"""
from pathlib import Path
import pandas as pd
parts = [Path(x) for x in {repr([str(p) for p in pred_files])}]
out = Path({str(LOCALIZATION)!r})
df = pd.concat([pd.read_csv(p) for p in parts], ignore_index=True)
df.to_csv(out, index=False)
print(len(df))
for p in parts:
    p.unlink(missing_ok=True)
for p in Path({str(cf_func)!r}).glob({(FUNCTIONAL_CF.stem + '.batch*.annotations')!r}):
    p.unlink(missing_ok=True)
"""
    run([str(CARVEME_VENV / "bin/python"), "-c", code], log)


def run_workflow() -> None:
    ensure_nospace_tools()
    status = check()
    missing = [k for k, ok in status["required"].items() if not ok]
    if missing:
        raise SystemExit("Missing required dependencies/files: " + ", ".join(missing))
    prepare_fasta()
    ensure_compartment_models_link()
    cf_data = CREATE_MODEL.parent.parent / "data/annotations"
    cf_func = cf_data / "functional"
    cf_fasta = cf_data / "fasta"
    cf_loc = cf_data / "localization"
    for d in [cf_func, cf_fasta, cf_loc, CREATE_MODEL.parent.parent / "results"]:
        d.mkdir(parents=True, exist_ok=True)
    with LOG.open("w") as log:
        if FUNCTIONAL.exists() and FUNCTIONAL.stat().st_size:
            print(f"Reusing existing eggNOG annotations: {FUNCTIONAL}", file=log, flush=True)
        else:
            run([
                str(CARVEME_VENV / "bin/python"), str(EGGNOG),
                "-i", str(PREP_FASTA), "--itype", "proteins",
                "-m", "diamond", "--tax_scope", "Fungi",
                "--output", "Bipolaris_sorokiniana_proteins.faa",
                "--output_dir", str(FUNCTIONAL.parent),
                "--override", "--cpu", EGGNOG_CPUS, "--annot_batch_size", "500",
            ], log)
        write_carvefungi_eggnog_compat(log)
        shutil.copyfile(FUNCTIONAL_CF, cf_func / FUNCTIONAL_CF.name)
        shutil.copyfile(PREP_FASTA, cf_fasta / PREP_FASTA.name)
        split_annotated_fastas(log)
        for fasta in sorted(PSIPRED_INPUTS.glob("*.fasta")):
            if not fasta.with_suffix(".horiz").exists():
                run([str(PSIPRED_RUNSEC), str(fasta.with_suffix(""))], log)
        if LOCALIZATION.exists() and LOCALIZATION.stat().st_size:
            print(f"Reusing existing localization predictions: {LOCALIZATION}", file=log, flush=True)
        else:
            run_compartment_prediction_batched(log, cf_func)
        if LOCALIZATION.exists():
            shutil.copyfile(LOCALIZATION, cf_loc / LOCALIZATION.name)
        code = (
            "import sys; sys.path.insert(0, {!r}); import CreateModelEggNogPool; "
            "CreateModelEggNogPool.carveFungi_pipeline('Bipolaris_sorokiniana', {!r}, {!r})"
        ).format(str(CREATE_MODEL.parent), str(FUNCTIONAL_CF), str(LOCALIZATION))
        run([str(CARVEME_VENV / "bin/python"), "-c", code], log, cwd=CREATE_MODEL.parent)
        merged = sorted((CREATE_MODEL.parent.parent / "results").glob("Bipolaris_sorokiniana--*.sbml"))
        candidates = merged or sorted((CREATE_MODEL.parent.parent / "results").glob("Bipolaris_sorokiniana*.sbml"))
        if candidates:
            shutil.copyfile(candidates[-1], MODEL_OUT)
        if MODEL_OUT.exists():
            stats = "import cobra, json; m=cobra.io.read_sbml_model({!r}); s=m.optimize(); print(json.dumps({{'reactions':len(m.reactions),'metabolites':len(m.metabolites),'genes':len(m.genes),'status':s.status,'growth':float(s.objective_value) if s.objective_value is not None else None}}))".format(str(MODEL_OUT))
            run([str(PROJECT_PY), "-c", stats], log)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", action="store_true", help="execute workflow; default only checks readiness")
    args = parser.parse_args()
    ROOT.joinpath("04_models").mkdir(exist_ok=True)
    ROOT.joinpath("05_logs").mkdir(exist_ok=True)
    ROOT.joinpath("06_reports").mkdir(exist_ok=True)
    status = check()
    REPORT.write_text(json.dumps(status, indent=2))
    print(json.dumps(status, indent=2))
    if args.run:
        run_workflow()
    else:
        print(f"Check-only report written: {REPORT}")
        print("Use --run only after missing dependencies/files are resolved.")


if __name__ == "__main__":
    main()
