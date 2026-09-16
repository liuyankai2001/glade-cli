"""Controlled real GenBank/manifest fixtures; no network or fake runtime data."""

import hashlib
import json
from pathlib import Path

from Bio import SeqIO

from src.config.run_config import RunConfig
from src.expression_box.expression_burden import (
    calculate_expression_burden,
    expression_burden_summary,
)
from test_plasmid_sequence import expression_record


def make_project(root: Path, *, count=1, enzymes=None):
    config = RunConfig(target_name="C00031")
    config.project_output_path = root / "outputs" / "C00031"
    config.manifest_output_path = config.project_output_path / "design_manifest.json"
    config.project_output_path.mkdir(parents=True)
    genes = [
        {
            "accession": f"protein_{letter}",
            "cds_length_nt": 9,
            "rbs": {"activity_percentile": 50},
            "ostir": {"translation_initiation_rate": 100, "unintended_start_count": 0},
        }
        for letter in "abc"
    ]
    summary = expression_burden_summary(
        calculate_expression_burden(
            [
                {
                    "cassette_index": 1,
                    "promoter": {"activity_percentile": 50},
                    "genes": genes,
                }
            ],
            {},
            fallback_promoter_percentile=50,
        )
    )
    constructs = []
    for design_id in range(1, count + 1):
        path = config.project_output_path / f"expression_{design_id}.gb"
        record = expression_record()
        SeqIO.write(record, path, "genbank")
        constructs.append(
            {
                "parts_design_id": design_id,
                "rank": design_id,
                "path": path.name,
                "length_bp": len(record),
                "sequence_sha256": hashlib.sha256(str(record.seq).encode()).hexdigest(),
                "file_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                "cassette_count": 1,
                "component_count": 5,
                "cassette_ranges": [
                    {"protein_accessions": [g["accession"] for g in genes]}
                ],
            }
        )
    manifest = {
        "schema_version": "design_manifest.v1",
        "target_compound_id": "C00031",
        "revision": 1,
        "cds_selection": {
            "host": {"chassis_key": "ecoli_mg1655", "name": "Escherichia coli MG1655"},
            "restriction_enzymes": ["eCoRI", "hINDiii"] if enzymes is None else enzymes,
        },
        "parts_selection": {
            "schema_version": "parts_selection.v2",
            "status": "selected",
            "selection_fingerprint": "a" * 64,
            "selected_design_ids": list(range(1, count + 1)),
            "design_references": [
                {"design_id": i, "expression_burden": summary}
                for i in range(1, count + 1)
            ],
        },
        "assembled_expression_constructs": {
            "schema_version": "assembled_expression_constructs.v1",
            "status": "assembled",
            "source_parts_selection_fingerprint": "a" * 64,
            "design_count": count,
            "constructs": constructs,
        },
    }
    config.manifest_output_path.write_text(json.dumps(manifest), encoding="utf-8")
    return config
