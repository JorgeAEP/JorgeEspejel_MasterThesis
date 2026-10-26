#!/usr/bin/env python3
"""Build the consolidated identifier-harmonisation table for the archive."""

from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
HARMONISATION = ROOT / "community_modelling" / "carvefungi" / "harmonisation"
OUTPUT = HARMONISATION / "consolidated_identifier_harmonisation.csv"
LATEX_OUTPUT = ROOT.parent / "Writing_Master_Thesis" / "appendix_identifier_harmonisation.tex"


FIELDS = [
    "mapping_scope",
    "source_identifier",
    "kegg_id",
    "chemical_name",
    "target_metabolite_id",
    "target_exchange_id",
    "mapping_status",
    "equivalence_basis",
]


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def build_rows() -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []

    for record in read_rows(HARMONISATION / "fungal_curated76_interface_changes.csv"):
        created = any(
            record.get(field) == "True"
            for field in ("created_external_metabolite", "created_exchange", "created_transport")
        )
        rows.append(
            {
                "mapping_scope": "CarveFungi fungal interface",
                "source_identifier": record["kegg_id"],
                "kegg_id": record["kegg_id"],
                "chemical_name": record["metabolite_name"],
                "target_metabolite_id": record["external_id"],
                "target_exchange_id": record["exchange_id"],
                "mapping_status": "created" if created else "retained/reopened",
                "equivalence_basis": "curated KEGG-to-fungal interface mapping",
            }
        )

    for record in read_rows(ROOT / "community_modelling" / "raven" / "harmonisation" / "curated76_bridge_table.csv"):
        rows.append(
            {
                "mapping_scope": "RAVEN bridge",
                "source_identifier": record["kegg_id"],
                "kegg_id": record["kegg_id"],
                "chemical_name": record["metabolite_name"],
                "target_metabolite_id": record["external_id"],
                "target_exchange_id": record["exchange_id"],
                "mapping_status": record["fungal_presence_status"],
                "equivalence_basis": "curated KEGG-to-RAVEN bridge mapping",
            }
        )

    medium_records = read_rows(HARMONISATION / "step88_to_carvefungi_exchange_mapping_by_scenario.csv")
    grouped: dict[tuple[str, ...], dict[str, object]] = {}
    for record in medium_records:
        key = (
            record["step88_community_exchange"],
            record["kegg"],
            record["carvefungi_metabolite"],
            record["carvefungi_exchange"],
        )
        entry = grouped.setdefault(
            key,
            {
                "record": record,
                "scenarios": set(),
                "statuses": set(),
            },
        )
        entry["scenarios"].add(record["scenario"])
        entry["statuses"].add(record["mapping_status"])

    for entry in grouped.values():
        record = entry["record"]
        statuses = sorted(entry["statuses"])
        rows.append(
            {
                "mapping_scope": "CarveFungi medium mapping",
                "source_identifier": record["step88_community_exchange"],
                "kegg_id": record["kegg"],
                "chemical_name": record["carvefungi_metabolite_name"],
                "target_metabolite_id": record["carvefungi_metabolite"],
                "target_exchange_id": record["carvefungi_exchange"],
                "mapping_status": ";".join(statuses),
                "equivalence_basis": "medium exchange mapped through KEGG compound identity",
            }
        )

    return sorted(
        rows,
        key=lambda row: (
            row["mapping_scope"],
            row["kegg_id"],
            row["source_identifier"],
            row["target_exchange_id"],
        ),
    )


def latex_escape(value: str) -> str:
    replacements = {
        "\\": r"\textbackslash{}",
        "&": r"\&",
        "%": r"\%",
        "_": r"\_\allowbreak{}",
        "#": r"\#",
    }
    return "".join(replacements.get(char, char) for char in value)


def write_latex(rows: list[dict[str, str]]) -> None:
    lines = [
        r"% Generated from consolidated_identifier_harmonisation.csv.",
        r"\tiny",
        r"\setlength{\tabcolsep}{2pt}",
        r"\setlength{\emergencystretch}{2em}",
        r"\begin{longtable}{@{}p{1.8cm}p{1.5cm}p{0.9cm}p{2.0cm}p{1.5cm}p{1.7cm}p{1.8cm}p{2.0cm}@{}}",
        r"\caption{Consolidated identifier harmonisation table showing the source identifiers, KEGG chemical equivalents and target model interfaces used for community modelling.}\label{tab:identifier-harmonisation}\\",
        r"\toprule",
        r"Mapping scope & Source identifier & KEGG ID & Chemical equivalent & Target metabolite ID & Target exchange ID & Status & Equivalence basis \\",
        r"\midrule",
        r"\endfirsthead",
        r"\toprule",
        r"Mapping scope & Source identifier & KEGG ID & Chemical equivalent & Target metabolite ID & Target exchange ID & Status & Equivalence basis \\",
        r"\midrule",
        r"\endhead",
        r"\bottomrule",
        r"\endfoot",
    ]
    for row in rows:
        values = [latex_escape(row[field]) for field in FIELDS]
        lines.append(" & ".join(values) + r" \\")
    lines.extend([r"\end{longtable}", r"\normalsize", ""])
    LATEX_OUTPUT.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    rows = build_rows()
    with OUTPUT.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    write_latex(rows)
    print(f"wrote {len(rows)} rows to {OUTPUT}")
    print(f"wrote LaTeX table to {LATEX_OUTPUT}")


if __name__ == "__main__":
    main()
