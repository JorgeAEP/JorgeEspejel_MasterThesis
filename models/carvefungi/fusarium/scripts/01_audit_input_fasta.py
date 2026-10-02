#!/usr/bin/env python3
"""Audit the Fusarium protein FASTA without modifying either input copy."""

from __future__ import annotations

import hashlib
import json
import re
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE_FASTA = Path(
    "/home/tlu/Desktop/Jorge/MATLABWORK/fusarium_graminearum/"
    "Fusarium_graminearum_proteins_clean.faa"
)
INPUT_FASTA = ROOT / "01_Input/Fusarium_graminearum_proteins.faa"
REPORT_DIR = ROOT / "06_reports/01_input_fasta_audit"
REPORT_PATH = REPORT_DIR / "input_fasta_audit.json"

STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")
PERMITTED_EXTRA = set("BXZJUO*")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def audit_fasta(path: Path) -> dict[str, object]:
    sequence_count = 0
    empty_sequences: list[str] = []
    duplicate_ids: list[str] = []
    invalid_characters: Counter[str] = Counter()
    extra_characters: Counter[str] = Counter()
    lengths: list[int] = []
    seen_ids: set[str] = set()
    current_id: str | None = None
    current_length = 0

    def finish_record() -> None:
        nonlocal current_id, current_length
        if current_id is None:
            return
        lengths.append(current_length)
        if current_length == 0:
            empty_sequences.append(current_id)

    with path.open(encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                finish_record()
                sequence_count += 1
                header = line[1:].strip()
                if not header:
                    raise ValueError(f"Empty FASTA header at line {line_number}")
                current_id = header.split()[0]
                if current_id in seen_ids:
                    duplicate_ids.append(current_id)
                seen_ids.add(current_id)
                current_length = 0
                continue
            if current_id is None:
                raise ValueError(f"Sequence data before first header at line {line_number}")
            sequence = re.sub(r"\s+", "", line).upper()
            current_length += len(sequence)
            for character in sequence:
                if character in STANDARD_AA:
                    continue
                if character in PERMITTED_EXTRA:
                    extra_characters[character] += 1
                else:
                    invalid_characters[character] += 1

    finish_record()
    total_residues = sum(lengths)
    return {
        "sequence_count": sequence_count,
        "unique_identifier_count": len(seen_ids),
        "duplicate_identifiers": sorted(set(duplicate_ids)),
        "empty_sequences": empty_sequences,
        "total_residues": total_residues,
        "minimum_length": min(lengths) if lengths else None,
        "maximum_length": max(lengths) if lengths else None,
        "mean_length": total_residues / len(lengths) if lengths else None,
        "permitted_nonstandard_character_counts": dict(sorted(extra_characters.items())),
        "invalid_character_counts": dict(sorted(invalid_characters.items())),
    }


def main() -> int:
    missing = [str(path) for path in (SOURCE_FASTA, INPUT_FASTA) if not path.is_file()]
    if missing:
        print("Missing required FASTA file(s):", file=sys.stderr)
        for path in missing:
            print(f"  {path}", file=sys.stderr)
        return 2

    source_before = sha256(SOURCE_FASTA)
    input_before = sha256(INPUT_FASTA)
    fasta_audit = audit_fasta(INPUT_FASTA)
    source_after = sha256(SOURCE_FASTA)
    input_after = sha256(INPUT_FASTA)

    report = {
        "audit_timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "methodological_precedent": str(ROOT.parent / "Bipolaris_sorokiniana_carveme"),
        "source_fasta": str(SOURCE_FASTA),
        "workspace_fasta": str(INPUT_FASTA),
        "source_sha256_before": source_before,
        "source_sha256_after": source_after,
        "workspace_sha256_before": input_before,
        "workspace_sha256_after": input_after,
        "source_unchanged_during_audit": source_before == source_after,
        "workspace_copy_unchanged_during_audit": input_before == input_after,
        "source_and_workspace_copy_identical": source_after == input_after,
        "file_size_bytes": INPUT_FASTA.stat().st_size,
        "fasta": fasta_audit,
    }
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_PATH.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    passed = (
        report["source_unchanged_during_audit"]
        and report["workspace_copy_unchanged_during_audit"]
        and report["source_and_workspace_copy_identical"]
        and fasta_audit["sequence_count"] > 0
        and not fasta_audit["duplicate_identifiers"]
        and not fasta_audit["empty_sequences"]
        and not fasta_audit["invalid_character_counts"]
    )
    print(f"Source SHA256:    {source_after}")
    print(f"Workspace SHA256: {input_after}")
    print(f"Sequences:        {fasta_audit['sequence_count']}")
    print(f"Residues:         {fasta_audit['total_residues']}")
    print(f"Length min/mean/max: {fasta_audit['minimum_length']}/"
          f"{fasta_audit['mean_length']:.2f}/{fasta_audit['maximum_length']}")
    print(f"Audit status:     {'PASS' if passed else 'FAIL'}")
    print(f"Report:           {REPORT_PATH}")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
