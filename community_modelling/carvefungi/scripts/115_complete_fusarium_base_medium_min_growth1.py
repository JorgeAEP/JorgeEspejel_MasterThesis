#!/usr/bin/env python3
"""Complete the base apoplastic medium for Fusarium growth >= 1.

This is a focused copy-only diagnostic. It starts from the pre-Fusarium
shared/apoplastic medium translated to Fusarium CarveFungi model exchange IDs,
then asks MICOM/COBRA's ``complete_medium`` to add the smallest number of
model exchange components required to reach a Fusarium growth rate of 1.

Two Fusarium biomass formulations are tested:

1. shared 52-component biomass
2. Fusarium fujikuroi-derived species-specific biomass

No source model, manifest, or medium file is overwritten.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
from datetime import datetime
from pathlib import Path

import cobra
from cobra.exceptions import OptimizationError
import pandas as pd
from micom.media import complete_medium


THIS_SCRIPT = Path(__file__).resolve()
PACKAGE = THIS_SCRIPT.parents[1]
WORKSPACE = PACKAGE.parent
OUT = PACKAGE / "reports/115_fusarium_base_medium_completion_min_growth1"

SCRIPT_06 = (
    WORKSPACE
    / "Fusarium_graminearum_carveme/03_scripts/"
    "06_test_bipolaris_52_component_biomass_and_shared_medium_copy_only.py"
)

MODELS = {
    "fusarium_52_biomass": PACKAGE
    / "fungal_models/biomass_variant_panel/Fusarium_graminearum_carvefungi_52biomass.xml",
    "fusarium_specific_biomass": PACKAGE
    / "fungal_models/biomass_variant_panel/Fusarium_graminearum_carvefungi_fujikuroi_biomass.xml",
}

ZERO = 1e-9


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


S06 = load_module(SCRIPT_06, "fusarium_script06")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    fields = sorted({key for row in rows for key in row}) if rows else ["empty"]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def optimize_with_medium(model: cobra.Model, medium: dict[str, float]) -> dict[str, object]:
    test = model.copy()
    test.medium = {
        reaction_id: float(flux)
        for reaction_id, flux in medium.items()
        if test.reactions.has_id(reaction_id) and abs(float(flux)) > ZERO
    }
    solution = test.optimize()
    value = float(solution.objective_value) if solution.objective_value is not None else 0.0
    return {
        "status": solution.status,
        "growth": value,
        "applied_medium_rows": len(test.medium),
    }


def exchange_metadata(model: cobra.Model, reaction_id: str) -> dict[str, object]:
    reaction = model.reactions.get_by_id(reaction_id)
    mets = list(reaction.metabolites)
    met = mets[0] if mets else None
    return {
        "reaction": reaction_id,
        "reaction_name": reaction.name,
        "metabolite_id": met.id if met is not None else "",
        "metabolite_name": met.name if met is not None else "",
        "formula": met.formula if met is not None else "",
        "charge": met.charge if met is not None else "",
        "has_carbon": bool(S06.exchange_has_carbon(reaction)),
    }


def build_base_medium_on_model(model: cobra.Model) -> tuple[dict[str, float], list[dict[str, object]]]:
    """Map the shared base medium to Fusarium model-level exchange IDs.

    Uses the same mapping table and temporary shared-interface additions that
    were used in the earlier Fusarium shared-medium diagnostics.
    """

    equivalents = S06.read_shared_medium_equivalents()
    shared_existing, mapping_rows = S06.build_shared_medium(model, equivalents, allow_missing=False)

    all_shared_medium = {
        str(row["carvefungi_exchange"]): float(row["full_flux"])
        for row in equivalents
    }
    missing_interfaces = sorted(set(all_shared_medium).difference(shared_existing))
    interface_rows = S06.add_all_missing_temporary_interfaces(model, missing_interfaces)

    base_medium = {
        reaction_id: flux
        for reaction_id, flux in all_shared_medium.items()
        if model.reactions.has_id(reaction_id) and abs(float(flux)) > ZERO
    }

    rows = []
    for row in mapping_rows:
        rid = str(row["carvefungi_exchange"])
        rows.append(
            {
                **row,
                "present_after_temporary_interfaces": model.reactions.has_id(rid),
                "included_in_base_medium_for_completion": rid in base_medium,
            }
        )
    return base_medium, rows + interface_rows


def completion_for_model(
    model_key: str,
    model_path: Path,
    growth: float,
    max_import: float,
    solver: str,
    minimize_components: bool,
) -> dict[str, object]:
    model = cobra.io.read_sbml_model(str(model_path))
    model.id = f"{model.id}_completion_copy"
    try:
        model.solver = solver
    except Exception:
        pass

    objective_ids = [r.id for r in model.reactions if abs(r.objective_coefficient) > ZERO]
    base_medium, mapping_rows = build_base_medium_on_model(model)
    base_result = optimize_with_medium(model, base_medium)

    medium_series = pd.Series(base_medium, dtype=float)
    try:
        completed = complete_medium(
            model,
            medium_series,
            growth=growth,
            max_import=max_import,
            minimize_components=minimize_components,
            strict=list(medium_series.index),
        )
    except OptimizationError as err:
        out_dir = OUT / model_key
        out_dir.mkdir(parents=True, exist_ok=True)
        write_csv(out_dir / "base_medium_mapping_audit.csv", mapping_rows)
        summary = {
            "model_key": model_key,
            "model_path": str(model_path),
            "model_sha256": sha256(model_path),
            "objective_reactions": objective_ids,
            "growth_target": growth,
            "max_import_for_added_components": max_import,
            "minimize_components": minimize_components,
            "base_medium_rows": len(base_medium),
            "base_growth_status": base_result["status"],
            "base_growth": base_result["growth"],
            "completed_medium_rows": None,
            "added_component_count": None,
            "completed_growth_status": "medium_completion_failed",
            "completed_growth": None,
            "added_component_ids": [],
            "error": str(err),
            "out_dir": str(out_dir),
        }
        (out_dir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
        return summary
    completed_medium = {str(k): float(v) for k, v in completed.items() if float(v) > ZERO}
    completed_result = optimize_with_medium(model, completed_medium)

    added_ids = sorted(set(completed_medium).difference(base_medium))
    added_rows = []
    for reaction_id in added_ids:
        added_rows.append(
            {
                "model_key": model_key,
                **exchange_metadata(model, reaction_id),
                "completed_flux": completed_medium[reaction_id],
                "in_native_model_medium": reaction_id in model.medium,
                "native_model_medium_flux": float(model.medium.get(reaction_id, 0.0)),
            }
        )

    full_rows = []
    for reaction_id, flux in sorted(completed_medium.items()):
        row = exchange_metadata(model, reaction_id)
        row.update(
            {
                "model_key": model_key,
                "completed_flux": flux,
                "source": "base_medium" if reaction_id in base_medium else "completion_added",
            }
        )
        full_rows.append(row)

    out_dir = OUT / model_key
    out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(out_dir / "base_medium_mapping_audit.csv", mapping_rows)
    write_csv(out_dir / "completion_added_components.csv", added_rows)
    write_csv(out_dir / "completed_medium_model_exchange_rows.csv", full_rows)
    pd.Series(completed_medium, name="flux").to_csv(out_dir / "completed_medium_series.csv", header=True)

    summary = {
        "model_key": model_key,
        "model_path": str(model_path),
        "model_sha256": sha256(model_path),
        "objective_reactions": objective_ids,
        "growth_target": growth,
        "max_import_for_added_components": max_import,
        "minimize_components": minimize_components,
        "base_medium_rows": len(base_medium),
        "base_growth_status": base_result["status"],
        "base_growth": base_result["growth"],
        "completed_medium_rows": len(completed_medium),
        "added_component_count": len(added_ids),
        "completed_growth_status": completed_result["status"],
        "completed_growth": completed_result["growth"],
        "added_component_ids": added_ids,
        "out_dir": str(out_dir),
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--growth", type=float, default=1.0)
    parser.add_argument("--max-import", type=float, default=100.0)
    parser.add_argument("--solver", default="cplex")
    parser.add_argument(
        "--minimize-components",
        action="store_true",
        default=True,
        help="Minimize number of added medium components. Enabled by default.",
    )
    parser.add_argument(
        "--minimize-total-import",
        action="store_true",
        help="Minimize total import flux instead of number of components.",
    )
    args = parser.parse_args()

    minimize_components = not args.minimize_total_import

    OUT.mkdir(parents=True, exist_ok=True)
    print("Script 115 Fusarium base-medium completion")
    print(f"started: {datetime.now().isoformat(timespec='seconds')}")
    print(f"package: {PACKAGE}")
    print(f"report_root: {OUT}")
    print(f"growth_target: {args.growth}")
    print(f"max_import_for_added_components: {args.max_import}")
    print(f"minimize_components: {minimize_components}")

    summaries = []
    total = len(MODELS)
    for idx, (model_key, model_path) in enumerate(MODELS.items(), start=1):
        print(f"[{idx}/{total}] {model_key}", flush=True)
        summary = completion_for_model(
            model_key=model_key,
            model_path=model_path,
            growth=args.growth,
            max_import=args.max_import,
            solver=args.solver,
            minimize_components=minimize_components,
        )
        summaries.append(summary)
        print(
            "  "
            f"base_growth={summary['base_growth']} | "
            f"completed_growth={summary['completed_growth']} | "
            f"added={summary['added_component_count']}",
            flush=True,
        )

    pd.DataFrame(summaries).to_csv(OUT / "completion_summary.csv", index=False)
    (OUT / "run_summary.json").write_text(
        json.dumps(
            {
                "created": datetime.now().isoformat(timespec="seconds"),
                "growth_target": args.growth,
                "max_import": args.max_import,
                "minimize_components": minimize_components,
                "reports": summaries,
            },
            indent=2,
        ),
        encoding="utf-8",
    )

    print("Script 115 compact summary")
    for summary in summaries:
        print(
            f"{summary['model_key']}: "
            f"base_growth={summary['base_growth']} "
            f"completed_growth={summary['completed_growth']} "
            f"added_components={summary['added_component_count']}"
        )
    print(f"outputs: {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
