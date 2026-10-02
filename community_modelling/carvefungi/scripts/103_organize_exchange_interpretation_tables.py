#!/usr/bin/env python3
"""Organize MICOM exchange fluxes into interpretation-ready tables.

This script does not run MICOM. It reads existing exchange flux outputs from:

- Script 96 fungal-species × SynCom factorial panel
- Script 100 Fusarium-active panel

and writes smaller tables for:

- medium uptake/secretion;
- fungal exchange fluxes;
- bacterial exchange fluxes;
- candidate bacteria-to-fungus cross-feeding;
- candidate fungus-to-bacteria cross-feeding;
- resource competition;
- HvSC1-vs-HvSC2 exchange differences;
- top exchanges by run.
"""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path


PACKAGE_ROOT = Path(__file__).resolve().parents[1]
REPORTS = PACKAGE_ROOT / "reports"

PANEL96_SUMMARIES = REPORTS / "96_fungal_species_syncom_factorial_panel" / "summaries"
PANEL100_RAW = REPORTS / "100_fusarium_active_growth_micom_panel" / "raw_micom_runs"

OUTDIR = REPORTS / "103_exchange_interpretation_tables"

EPS = 1e-6
TOP_N_PER_RUN = 50


LONG_COLUMNS = [
    "source_panel",
    "run_id",
    "manifest",
    "scenario",
    "biomass_mode",
    "syncom",
    "fungus_context",
    "community_type",
    "organism",
    "role",
    "fungus",
    "exchange",
    "metabolite_key",
    "flux",
    "flux_direction",
]


SUMMARY_COLUMNS = [
    "run_id",
    "manifest",
    "scenario",
    "biomass_mode",
    "syncom",
    "fungus_context",
    "community_type",
    "exchange",
    "metabolite_key",
    "medium_flux",
    "medium_uptake",
    "medium_secretion",
    "fungal_net_flux",
    "fungal_uptake_total",
    "fungal_secretion_total",
    "bacterial_net_flux",
    "bacterial_uptake_total",
    "bacterial_secretion_total",
    "fungal_producer_count",
    "fungal_consumer_count",
    "bacterial_producer_count",
    "bacterial_consumer_count",
    "fungal_producers",
    "fungal_consumers",
    "bacterial_producers",
    "bacterial_consumers",
    "candidate_bacteria_to_fungus",
    "candidate_fungus_to_bacteria",
    "candidate_resource_competition",
]


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({col: row.get(col, "") for col in columns})


def as_float(value: object) -> float | None:
    try:
        if value is None or value == "":
            return None
        return float(value)
    except (TypeError, ValueError):
        return None


def fmt(value: float | None) -> str:
    if value is None:
        return ""
    return f"{value:.12g}"


def metabolite_key(exchange: str) -> str:
    """Collapse common exchange namespace suffixes for comparison."""
    if exchange.endswith("_e") or exchange.endswith("_m"):
        return exchange[:-2]
    return exchange


def flux_direction(flux: float) -> str:
    if flux > EPS:
        return "secretion_positive"
    if flux < -EPS:
        return "uptake_negative"
    return "near_zero"


def infer_role(organism: str) -> tuple[str, str]:
    lowered = organism.lower()
    if organism == "medium":
        return "medium", "none"
    if "bipolaris" in lowered:
        return "fungus", "Bs"
    if "fusarium" in lowered:
        return "fungus", "Fg"
    return "bacteria", "none"


def infer_metadata_from_manifest(manifest: str, scenario: str) -> dict[str, str]:
    if "bipolaris" in manifest and "fusarium" in manifest:
        fungus_context = "Bs+Fg"
    elif "fungi_pair" in manifest:
        fungus_context = "Bs+Fg"
    elif "bipolaris" in manifest:
        fungus_context = "Bs"
    elif "fusarium" in manifest:
        fungus_context = "Fg"
    else:
        fungus_context = "none"

    if "hvsc1" in manifest:
        syncom = "HvSC1"
    elif "hvsc2" in manifest:
        syncom = "HvSC2"
    else:
        syncom = "none"

    if "fungi_pair" in manifest:
        community_type = "dual_fungus_syncom" if syncom != "none" else "fungi_pair"
    elif "fungus_alone" in manifest:
        community_type = "fungus_alone"
    elif "bacteria_only" in manifest:
        community_type = "bacteria_only"
    elif fungus_context in {"Bs", "Fg"} and syncom != "none":
        community_type = "single_fungus_syncom"
    else:
        community_type = ""

    biomass_mode = "species_specific" if "specific" in manifest else "shared_52"

    return {
        "manifest": manifest,
        "scenario": scenario,
        "biomass_mode": biomass_mode,
        "syncom": syncom,
        "fungus_context": fungus_context,
        "community_type": community_type,
    }


def load_panel96_long() -> list[dict[str, object]]:
    path = PANEL96_SUMMARIES / "exchange_flux_long_nonzero.csv"
    rows = []
    for row in read_csv(path):
        flux = as_float(row.get("flux"))
        if flux is None or abs(flux) <= EPS:
            continue
        exchange = row.get("exchange", "")
        rows.append(
            {
                "source_panel": "96_fungal_species_syncom_factorial_panel",
                "run_id": row.get("run_id", ""),
                "manifest": row.get("manifest", ""),
                "scenario": row.get("scenario", ""),
                "biomass_mode": row.get("biomass_mode", ""),
                "syncom": row.get("syncom", ""),
                "fungus_context": row.get("fungus_context", ""),
                "community_type": row.get("community_type", ""),
                "organism": row.get("organism", ""),
                "role": row.get("role", ""),
                "fungus": row.get("fungus", ""),
                "exchange": exchange,
                "metabolite_key": metabolite_key(exchange),
                "flux": flux,
                "flux_direction": flux_direction(flux),
            }
        )
    return rows


def load_panel100_long() -> list[dict[str, object]]:
    rows = []
    for path in sorted(PANEL100_RAW.glob("*/*/*_emma_style_fluxes.csv")):
        scenario = path.parent.name
        manifest = path.parent.parent.name
        meta = infer_metadata_from_manifest(manifest, scenario)
        run_id = f"{manifest}__{scenario}"
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle)
            for raw in reader:
                organism = raw.get("compartment", "")
                role, fungus = infer_role(organism)
                for exchange, value in raw.items():
                    if not exchange.startswith("EX_"):
                        continue
                    flux = as_float(value)
                    if flux is None or abs(flux) <= EPS:
                        continue
                    rows.append(
                        {
                            "source_panel": "100_fusarium_active_growth_micom_panel",
                            "run_id": run_id,
                            **meta,
                            "organism": organism,
                            "role": role,
                            "fungus": fungus,
                            "exchange": exchange,
                            "metabolite_key": metabolite_key(exchange),
                            "flux": flux,
                            "flux_direction": flux_direction(flux),
                        }
                    )
    return rows


def rows_for_csv(rows: list[dict[str, object]], columns: list[str]) -> list[dict[str, object]]:
    out = []
    for row in rows:
        new = dict(row)
        for key, value in list(new.items()):
            if isinstance(value, float):
                new[key] = fmt(value)
        out.append(new)
    return out


def summarize_by_run_exchange(long_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in long_rows:
        grouped[(str(row["run_id"]), str(row["exchange"]))].append(row)

    summaries = []
    for (_run_id, _exchange), group in grouped.items():
        first = group[0]

        medium_flux = 0.0
        fungal_net = 0.0
        bacterial_net = 0.0
        fungal_uptake = 0.0
        fungal_secretion = 0.0
        bacterial_uptake = 0.0
        bacterial_secretion = 0.0
        fungal_producers: set[str] = set()
        fungal_consumers: set[str] = set()
        bacterial_producers: set[str] = set()
        bacterial_consumers: set[str] = set()

        for row in group:
            flux = float(row["flux"])
            role = row.get("role", "")
            organism = str(row.get("organism", ""))
            if role == "medium":
                medium_flux += flux
            elif role == "fungus":
                fungal_net += flux
                if flux > EPS:
                    fungal_secretion += flux
                    fungal_producers.add(organism)
                elif flux < -EPS:
                    fungal_uptake += -flux
                    fungal_consumers.add(organism)
            elif role == "bacteria":
                bacterial_net += flux
                if flux > EPS:
                    bacterial_secretion += flux
                    bacterial_producers.add(organism)
                elif flux < -EPS:
                    bacterial_uptake += -flux
                    bacterial_consumers.add(organism)

        summary = {
            "run_id": first.get("run_id", ""),
            "manifest": first.get("manifest", ""),
            "scenario": first.get("scenario", ""),
            "biomass_mode": first.get("biomass_mode", ""),
            "syncom": first.get("syncom", ""),
            "fungus_context": first.get("fungus_context", ""),
            "community_type": first.get("community_type", ""),
            "exchange": first.get("exchange", ""),
            "metabolite_key": first.get("metabolite_key", ""),
            "medium_flux": medium_flux,
            "medium_uptake": -medium_flux if medium_flux < -EPS else 0.0,
            "medium_secretion": medium_flux if medium_flux > EPS else 0.0,
            "fungal_net_flux": fungal_net,
            "fungal_uptake_total": fungal_uptake,
            "fungal_secretion_total": fungal_secretion,
            "bacterial_net_flux": bacterial_net,
            "bacterial_uptake_total": bacterial_uptake,
            "bacterial_secretion_total": bacterial_secretion,
            "fungal_producer_count": len(fungal_producers),
            "fungal_consumer_count": len(fungal_consumers),
            "bacterial_producer_count": len(bacterial_producers),
            "bacterial_consumer_count": len(bacterial_consumers),
            "fungal_producers": ";".join(sorted(fungal_producers)),
            "fungal_consumers": ";".join(sorted(fungal_consumers)),
            "bacterial_producers": ";".join(sorted(bacterial_producers)),
            "bacterial_consumers": ";".join(sorted(bacterial_consumers)),
            "candidate_bacteria_to_fungus": bool(bacterial_producers and fungal_consumers),
            "candidate_fungus_to_bacteria": bool(fungal_producers and bacterial_consumers),
            "candidate_resource_competition": bool(fungal_consumers and bacterial_consumers),
        }
        summaries.append(summary)

    return sorted(
        summaries,
        key=lambda row: (
            str(row["run_id"]),
            -max(
                abs(float(row["medium_flux"])),
                abs(float(row["fungal_net_flux"])),
                abs(float(row["bacterial_net_flux"])),
            ),
        ),
    )


def filter_rows(rows: list[dict[str, object]], predicate) -> list[dict[str, object]]:
    return [row for row in rows if predicate(row)]


def top_exchange_summary(summary_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    grouped: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in summary_rows:
        grouped[str(row["run_id"])].append(row)
    out = []
    for run_id, rows in grouped.items():
        ranked = sorted(
            rows,
            key=lambda row: max(
                abs(float(row["medium_flux"])),
                abs(float(row["fungal_net_flux"])),
                abs(float(row["bacterial_net_flux"])),
                float(row["fungal_uptake_total"]),
                float(row["fungal_secretion_total"]),
                float(row["bacterial_uptake_total"]),
                float(row["bacterial_secretion_total"]),
            ),
            reverse=True,
        )
        for rank, row in enumerate(ranked[:TOP_N_PER_RUN], start=1):
            new = dict(row)
            new["rank_within_run"] = rank
            out.append(new)
    return out


def hvsc1_vs_hvsc2(summary_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    keyed: dict[tuple[str, str, str, str, str, str], dict[str, dict[str, object]]] = defaultdict(dict)
    for row in summary_rows:
        syncom = str(row.get("syncom", ""))
        if syncom not in {"HvSC1", "HvSC2"}:
            continue
        key = (
            str(row.get("scenario", "")),
            str(row.get("biomass_mode", "")),
            str(row.get("community_type", "")),
            str(row.get("fungus_context", "")),
            str(row.get("exchange", "")),
            str(row.get("metabolite_key", "")),
        )
        keyed[key][syncom] = row

    out = []
    for key, pair in keyed.items():
        if "HvSC1" not in pair or "HvSC2" not in pair:
            continue
        h1 = pair["HvSC1"]
        h2 = pair["HvSC2"]
        row = {
            "scenario": key[0],
            "biomass_mode": key[1],
            "community_type": key[2],
            "fungus_context": key[3],
            "exchange": key[4],
            "metabolite_key": key[5],
            "hvsc1_medium_flux": h1["medium_flux"],
            "hvsc2_medium_flux": h2["medium_flux"],
            "delta_medium_flux_hvsc1_minus_hvsc2": float(h1["medium_flux"])
            - float(h2["medium_flux"]),
            "hvsc1_fungal_net_flux": h1["fungal_net_flux"],
            "hvsc2_fungal_net_flux": h2["fungal_net_flux"],
            "delta_fungal_net_flux_hvsc1_minus_hvsc2": float(h1["fungal_net_flux"])
            - float(h2["fungal_net_flux"]),
            "hvsc1_bacterial_net_flux": h1["bacterial_net_flux"],
            "hvsc2_bacterial_net_flux": h2["bacterial_net_flux"],
            "delta_bacterial_net_flux_hvsc1_minus_hvsc2": float(h1["bacterial_net_flux"])
            - float(h2["bacterial_net_flux"]),
            "hvsc1_candidate_bacteria_to_fungus": h1["candidate_bacteria_to_fungus"],
            "hvsc2_candidate_bacteria_to_fungus": h2["candidate_bacteria_to_fungus"],
            "hvsc1_candidate_fungus_to_bacteria": h1["candidate_fungus_to_bacteria"],
            "hvsc2_candidate_fungus_to_bacteria": h2["candidate_fungus_to_bacteria"],
            "hvsc1_resource_competition": h1["candidate_resource_competition"],
            "hvsc2_resource_competition": h2["candidate_resource_competition"],
        }
        out.append(row)

    return sorted(
        out,
        key=lambda row: max(
            abs(float(row["delta_medium_flux_hvsc1_minus_hvsc2"])),
            abs(float(row["delta_fungal_net_flux_hvsc1_minus_hvsc2"])),
            abs(float(row["delta_bacterial_net_flux_hvsc1_minus_hvsc2"])),
        ),
        reverse=True,
    )


def write_readme(stats: dict[str, int]) -> None:
    text = f"""# Exchange interpretation tables

This folder reorganizes the large MICOM exchange-flux output into smaller,
question-oriented tables. No MICOM simulations were run by this step.

## Generated tables

- `01_exchange_flux_long_nonzero_combined.csv`
  - Combined non-zero exchange fluxes from Script 96 and Script 100.
  - One row = one organism × exchange × run.

- `02_exchange_summary_by_run_exchange.csv`
  - One row = one run × exchange.
  - Aggregates medium flux, fungal uptake/secretion, bacterial uptake/secretion,
    producer/consumer counts, and candidate interaction flags.

- `03_medium_uptake_by_run.csv`
  - Medium components consumed from the community medium.

- `04_fungal_exchange_fluxes_by_run.csv`
  - Fungal uptake and secretion flux rows.

- `05_bacterial_exchange_fluxes_by_run.csv`
  - Bacterial uptake and secretion flux rows.

- `06_candidate_bacteria_to_fungus_crossfeeding.csv`
  - Exchanges where at least one bacterium secretes and at least one fungus
    consumes the same exchange.

- `07_candidate_fungus_to_bacteria_crossfeeding.csv`
  - Exchanges where at least one fungus secretes and at least one bacterium
    consumes the same exchange.

- `08_resource_competition_candidates.csv`
  - Exchanges where fungi and bacteria both consume the same exchange.

- `09_hvsc1_vs_hvsc2_exchange_differences.csv`
  - Exchange-level HvSC1 minus HvSC2 differences for matched contexts.

- `10_top_exchange_summary_by_run.csv`
  - Top {TOP_N_PER_RUN} largest exchange signals per run.

## Row counts

```text
combined_long_rows: {stats['combined_long_rows']}
summary_run_exchange_rows: {stats['summary_run_exchange_rows']}
medium_uptake_rows: {stats['medium_uptake_rows']}
fungal_exchange_rows: {stats['fungal_exchange_rows']}
bacterial_exchange_rows: {stats['bacterial_exchange_rows']}
bacteria_to_fungus_candidate_rows: {stats['bacteria_to_fungus_candidate_rows']}
fungus_to_bacteria_candidate_rows: {stats['fungus_to_bacteria_candidate_rows']}
resource_competition_rows: {stats['resource_competition_rows']}
hvsc1_vs_hvsc2_rows: {stats['hvsc1_vs_hvsc2_rows']}
top_exchange_rows: {stats['top_exchange_rows']}
```

## How to interpret flux signs

For organism rows:

- positive flux = secretion/export;
- negative flux = uptake/import.

For medium rows:

- negative medium flux indicates medium uptake by the community;
- positive medium flux indicates net secretion to the medium.

Candidate cross-feeding and competition rows are hypotheses from flux patterns,
not proof of direct transfer. They should be used to prioritize follow-up
inspection or targeted blocking/supplementation tests.
"""
    (OUTDIR / "README.md").write_text(text, encoding="utf-8")


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)

    long_rows = load_panel96_long() + load_panel100_long()
    summary_rows = summarize_by_run_exchange(long_rows)

    medium_uptake = filter_rows(
        summary_rows, lambda row: float(row["medium_uptake"]) > EPS
    )
    fungal_fluxes = filter_rows(long_rows, lambda row: row.get("role") == "fungus")
    bacterial_fluxes = filter_rows(long_rows, lambda row: row.get("role") == "bacteria")
    b2f = filter_rows(summary_rows, lambda row: bool(row["candidate_bacteria_to_fungus"]))
    f2b = filter_rows(summary_rows, lambda row: bool(row["candidate_fungus_to_bacteria"]))
    competition = filter_rows(summary_rows, lambda row: bool(row["candidate_resource_competition"]))
    hvsc_diff = hvsc1_vs_hvsc2(summary_rows)
    top_rows = top_exchange_summary(summary_rows)

    write_csv(
        OUTDIR / "01_exchange_flux_long_nonzero_combined.csv",
        rows_for_csv(long_rows, LONG_COLUMNS),
        LONG_COLUMNS,
    )
    write_csv(
        OUTDIR / "02_exchange_summary_by_run_exchange.csv",
        rows_for_csv(summary_rows, SUMMARY_COLUMNS),
        SUMMARY_COLUMNS,
    )
    write_csv(
        OUTDIR / "03_medium_uptake_by_run.csv",
        rows_for_csv(medium_uptake, SUMMARY_COLUMNS),
        SUMMARY_COLUMNS,
    )
    write_csv(
        OUTDIR / "04_fungal_exchange_fluxes_by_run.csv",
        rows_for_csv(fungal_fluxes, LONG_COLUMNS),
        LONG_COLUMNS,
    )
    write_csv(
        OUTDIR / "05_bacterial_exchange_fluxes_by_run.csv",
        rows_for_csv(bacterial_fluxes, LONG_COLUMNS),
        LONG_COLUMNS,
    )
    write_csv(
        OUTDIR / "06_candidate_bacteria_to_fungus_crossfeeding.csv",
        rows_for_csv(b2f, SUMMARY_COLUMNS),
        SUMMARY_COLUMNS,
    )
    write_csv(
        OUTDIR / "07_candidate_fungus_to_bacteria_crossfeeding.csv",
        rows_for_csv(f2b, SUMMARY_COLUMNS),
        SUMMARY_COLUMNS,
    )
    write_csv(
        OUTDIR / "08_resource_competition_candidates.csv",
        rows_for_csv(competition, SUMMARY_COLUMNS),
        SUMMARY_COLUMNS,
    )

    hvsc_columns = [
        "scenario",
        "biomass_mode",
        "community_type",
        "fungus_context",
        "exchange",
        "metabolite_key",
        "hvsc1_medium_flux",
        "hvsc2_medium_flux",
        "delta_medium_flux_hvsc1_minus_hvsc2",
        "hvsc1_fungal_net_flux",
        "hvsc2_fungal_net_flux",
        "delta_fungal_net_flux_hvsc1_minus_hvsc2",
        "hvsc1_bacterial_net_flux",
        "hvsc2_bacterial_net_flux",
        "delta_bacterial_net_flux_hvsc1_minus_hvsc2",
        "hvsc1_candidate_bacteria_to_fungus",
        "hvsc2_candidate_bacteria_to_fungus",
        "hvsc1_candidate_fungus_to_bacteria",
        "hvsc2_candidate_fungus_to_bacteria",
        "hvsc1_resource_competition",
        "hvsc2_resource_competition",
    ]
    write_csv(
        OUTDIR / "09_hvsc1_vs_hvsc2_exchange_differences.csv",
        rows_for_csv(hvsc_diff, hvsc_columns),
        hvsc_columns,
    )
    write_csv(
        OUTDIR / "10_top_exchange_summary_by_run.csv",
        rows_for_csv(top_rows, ["rank_within_run"] + SUMMARY_COLUMNS),
        ["rank_within_run"] + SUMMARY_COLUMNS,
    )

    stats = {
        "combined_long_rows": len(long_rows),
        "summary_run_exchange_rows": len(summary_rows),
        "medium_uptake_rows": len(medium_uptake),
        "fungal_exchange_rows": len(fungal_fluxes),
        "bacterial_exchange_rows": len(bacterial_fluxes),
        "bacteria_to_fungus_candidate_rows": len(b2f),
        "fungus_to_bacteria_candidate_rows": len(f2b),
        "resource_competition_rows": len(competition),
        "hvsc1_vs_hvsc2_rows": len(hvsc_diff),
        "top_exchange_rows": len(top_rows),
    }
    write_readme(stats)

    print("Script 103 compact summary")
    for key, value in stats.items():
        print(f"{key}: {value}")
    print(f"outputs: {OUTDIR}")


if __name__ == "__main__":
    main()
