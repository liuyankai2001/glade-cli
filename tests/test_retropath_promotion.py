from __future__ import annotations

import csv
import hashlib
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from src.pathway_analyze.retropath_analyze import (
    CANDIDATE_ROUTE_COLUMNS,
    CANDIDATE_ROUTES_FILE_NAME,
    CANDIDATE_STEP_COLUMNS,
    CANDIDATE_STEPS_FILE_NAME,
)
from src.pathway_analyze.retropath_gem_validation import (
    HYPOTHESIS_COLUMNS,
    STOICHIOMETRY_HYPOTHESES_FILE_NAME,
    STOICHIOMETRY_TERMS_FILE_NAME,
    SUMMARY_COLUMNS,
    TERM_COLUMNS,
    VALIDATION_MANIFEST_FILE_NAME,
    VALIDATION_SUMMARY_FILE_NAME,
)
from src.pathway_analyze.retropath_promotion import (
    PROMOTION_MANIFEST_FILE_NAME,
    materialize_retropath_solutions,
    verify_retropath_solution_promotion,
)
from src.write_manifest.solution import write_solution


def _write_csv(path: Path, columns: tuple[str, ...], rows: list[dict]) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return hashlib.sha256(path.read_bytes()).hexdigest()


class _FakeKegg:
    def get_compound_name(self, compound_id: str) -> str:
        return {"C12345": "target", "C00010": "substrate"}.get(
            compound_id, compound_id
        )

    def try_get_reaction(self, reaction_id: str):
        return None


class _FakeMnxref:
    def __init__(self) -> None:
        self.template = SimpleNamespace(
            mnxr_id="MNXR1",
            terms=(
                SimpleNamespace(
                    side="left", coefficient=1.0, mnxm_id="MNXM10"
                ),
                SimpleNamespace(
                    side="right", coefficient=1.0, mnxm_id="MNXM11"
                ),
            ),
            reaction_xrefs=("kegg:R12345", "rhea:12345"),
        )

    def __enter__(self):
        return self

    def __exit__(self, *_args):
        return None

    def templates_for_rules(self, _rule_ids):
        return (self.template,)

    def chemicals(self, ids):
        values = {
            "MNXM10": SimpleNamespace(
                charge=0, inchikey="AAAAAAAAAAAAAA-BBBBBBBBBB-C"
            ),
            "MNXM11": SimpleNamespace(
                charge=0, inchikey="CCCCCCCCCCCCCC-DDDDDDDDDD-E"
            ),
        }
        return {value: values[value] for value in ids}


class RetroPathPromotionTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.gap_dir = self.root / "gap" / "depth0"
        self.retropath_dir = self.gap_dir / "retropath"
        self.p8_dir = self.retropath_dir / "gem_validation"
        self.candidate_id = "RP2ROUTE:" + "1" * 64
        self.step_id = "RP2STEP:" + "2" * 64
        self.hypothesis_id = "RP2STOICH:" + "3" * 64
        self.combination_id = "RP2GEM:" + "4" * 64
        self.rules = self.root / "rules.csv"
        self.rules.write_text(
            '"Rule ID","Rule","EC number"\n'
            '"RR-TEST","[C:1]>>[C:1]","1.1.1.1"\n',
            encoding="utf-8",
        )
        self.config = SimpleNamespace(
            target_name="C12345",
            depth=0,
            gap_output_path=self.root / "gap",
            retropath_rules_path=self.rules,
            cache_dir=self.root / "cache",
            data_dir=self.root / "data",
        )
        self._write_inputs()

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def _write_inputs(self) -> None:
        self.gap_dir.mkdir(parents=True)
        _write_csv(
            self.gap_dir / "solutions.csv",
            ("solution_id", "target_compound_id", "eligible_for_recommendation"),
            [{
                "solution_id": 1,
                "target_compound_id": "C12345",
                "eligible_for_recommendation": "true",
            }],
        )
        _write_csv(
            self.gap_dir / "all_solution_steps.csv",
            ("solution_id", "step_index", "reaction_id"),
            [{"solution_id": 1, "step_index": 1, "reaction_id": "R00001"}],
        )
        _write_csv(
            self.gap_dir / "solution_electron_summary.csv",
            ("solution_id", "electron_system_status"),
            [{"solution_id": 1, "electron_system_status": "not_required"}],
        )
        _write_csv(
            self.gap_dir / "route_electron_requirements.csv",
            ("solution_id", "step_index", "reaction_id"),
            [],
        )
        _write_csv(
            self.retropath_dir / CANDIDATE_ROUTES_FILE_NAME,
            CANDIDATE_ROUTE_COLUMNS,
            [{
                "candidate_rank": 1,
                "candidate_id": self.candidate_id,
                "target_compound_id": "C12345",
                "sink_kegg_ids": "C00010",
                "kegg_prefix_steps": 0,
                "retropath_steps": 1,
                "total_steps": 1,
                "maximum_sink_depth": 0,
                "upstream_enumeration_truncated": "false",
                "candidate_top_k_truncated": "false",
            }],
        )
        _write_csv(
            self.retropath_dir / CANDIDATE_STEPS_FILE_NAME,
            CANDIDATE_STEP_COLUMNS,
            [{
                "candidate_id": self.candidate_id,
                "step_index": 1,
                "step_id": self.step_id,
                "step_source": "retropath",
                "status": "predicted",
                "orientation": "biosynthetic",
                "direction": "left_to_right",
                "reaction_option_ids": "RP2:test",
                "reaction_smiles": "[C]>>[C]",
                "substrate_compound_ids": "C00010",
                "product_compound_ids": "C12345",
                "substrate_stoichiometry_json": '[["C00010",1.0]]',
                "product_stoichiometry_json": '[["C12345",1.0]]',
                "rule_ids": "RR-TEST",
                "source_reaction_ids": "MNXR1",
                "source_ec_numbers": "1.1.1.1",
                "source_uniprot_ids": "P12345",
                "expansion_depth": 0,
            }],
        )
        summary_hash = _write_csv(
            self.p8_dir / VALIDATION_SUMMARY_FILE_NAME,
            SUMMARY_COLUMNS,
            [{
                "candidate_rank": 1,
                "candidate_id": self.candidate_id,
                "combination_id": self.combination_id,
                "candidate_status": "PASS_STRICT_HYPOTHESIS_EXISTS",
                "validation_status": "PASS_STRICT_ROUTE_FLUX",
                "stoichiometry_hypothesis_ids": self.hypothesis_id,
                "fba_status": "optimal",
                "target_flux": 1.0,
                "growth_flux": 0.5,
                "pfba_status": "optimal",
                "pfba_target_flux": 1.0,
                "fva_status": "not_run",
                "route_step_count": 1,
                "active_route_step_count": 1,
                "combination_truncated": "false",
            }],
        )
        hypothesis_hash = _write_csv(
            self.p8_dir / STOICHIOMETRY_HYPOTHESES_FILE_NAME,
            HYPOTHESIS_COLUMNS,
            [{
                "candidate_rank": 1,
                "candidate_id": self.candidate_id,
                "step_id": self.step_id,
                "hypothesis_id": self.hypothesis_id,
                "rule_id": "RR-TEST",
                "source_mnxr_id": "MNXR1",
                "source_orientation": "left_to_right",
                "evidence_grade": "template_exact",
                "balance_status": "balanced",
                "cofactor_reconstruction_status": "complete",
            }],
        )
        term_hash = _write_csv(
            self.p8_dir / STOICHIOMETRY_TERMS_FILE_NAME,
            TERM_COLUMNS,
            [
                {
                    "candidate_rank": 1,
                    "candidate_id": self.candidate_id,
                    "step_id": self.step_id,
                    "hypothesis_id": self.hypothesis_id,
                    "side": "left",
                    "coefficient": 1,
                    "role": "core",
                    "compound_id": "C00010",
                    "source_mnxm_id": "MNXM10",
                    "name": "substrate",
                    "formula": "C",
                    "charge": 0,
                    "inchi": "InChI=1S/C",
                    "inchikey": "AAAAAAAAAAAAAA-BBBBBBBBBB-C",
                    "smiles": "[C]",
                },
                {
                    "candidate_rank": 1,
                    "candidate_id": self.candidate_id,
                    "step_id": self.step_id,
                    "hypothesis_id": self.hypothesis_id,
                    "side": "right",
                    "coefficient": 1,
                    "role": "core",
                    "compound_id": "C12345",
                    "source_mnxm_id": "MNXM11",
                    "name": "target",
                    "formula": "C",
                    "charge": 0,
                    "inchi": "InChI=1S/C",
                    "inchikey": "CCCCCCCCCCCCCC-DDDDDDDDDD-E",
                    "smiles": "[C]",
                },
            ],
        )
        manifest = {
            "schema_version": "retropath_gem_validation.v1",
            "target_compound": "C12345",
            "expansion_depth": 0,
            "inputs": {"fixture": True},
            "artifacts": {
                VALIDATION_SUMMARY_FILE_NAME: {
                    "path": str((self.p8_dir / VALIDATION_SUMMARY_FILE_NAME).resolve()),
                    "sha256": summary_hash,
                },
                STOICHIOMETRY_HYPOTHESES_FILE_NAME: {
                    "path": str((self.p8_dir / STOICHIOMETRY_HYPOTHESES_FILE_NAME).resolve()),
                    "sha256": hypothesis_hash,
                },
                STOICHIOMETRY_TERMS_FILE_NAME: {
                    "path": str((self.p8_dir / STOICHIOMETRY_TERMS_FILE_NAME).resolve()),
                    "sha256": term_hash,
                },
            },
        }
        self.validation_manifest = self.p8_dir / VALIDATION_MANIFEST_FILE_NAME
        self.validation_manifest.write_text(
            json.dumps(manifest), encoding="utf-8"
        )

    def test_pass_is_appended_after_kegg_and_committed(self) -> None:
        with (
            patch(
                "src.pathway_analyze.retropath_promotion.MnxrefIndex",
                return_value=_FakeMnxref(),
            ),
            patch(
                "src.pathway_analyze.retropath_promotion.KeggRestClient",
                return_value=_FakeKegg(),
            ),
        ):
            result = materialize_retropath_solutions(
                self.config,
                validation_manifest_path=self.validation_manifest,
            )
        self.assertEqual([2], result["formal_solution_ids"])
        with (self.gap_dir / "solutions.csv").open(
            "r", encoding="utf-8", newline=""
        ) as handle:
            summaries = list(csv.DictReader(handle))
        self.assertEqual(["1", "2"], [row["solution_id"] for row in summaries])
        self.assertEqual("kegg", summaries[0]["solution_source"])
        self.assertEqual("retropath", summaries[1]["solution_source"])
        with (self.gap_dir / "all_solution_steps.csv").open(
            "r", encoding="utf-8", newline=""
        ) as handle:
            steps = list(csv.DictReader(handle))
        predicted = next(row for row in steps if row["solution_id"] == "2")
        self.assertEqual(self.hypothesis_id, predicted["reaction_id"])
        self.assertEqual("true", predicted["formal_mapping_exact"])
        promotion = verify_retropath_solution_promotion(
            gap_dir=self.gap_dir,
            target_compound="C12345",
            expansion_depth=0,
            solution_id=2,
        )
        self.assertEqual(result["promotion_id"], promotion["promotion_id"])
        self.assertTrue(
            (self.retropath_dir / PROMOTION_MANIFEST_FILE_NAME).is_file()
        )

    def test_tampered_formal_csv_fails_closed(self) -> None:
        with (
            patch(
                "src.pathway_analyze.retropath_promotion.MnxrefIndex",
                return_value=_FakeMnxref(),
            ),
            patch(
                "src.pathway_analyze.retropath_promotion.KeggRestClient",
                return_value=_FakeKegg(),
            ),
        ):
            materialize_retropath_solutions(
                self.config,
                validation_manifest_path=self.validation_manifest,
            )
        with (self.gap_dir / "solutions.csv").open("a", encoding="utf-8") as handle:
            handle.write("tampered\n")
        with self.assertRaisesRegex(ValueError, "promotion artifact mismatch"):
            verify_retropath_solution_promotion(
                gap_dir=self.gap_dir,
                target_compound="C12345",
                expansion_depth=0,
                solution_id=2,
            )

    def test_write_solution_uses_common_manifest_section(self) -> None:
        with (
            patch(
                "src.pathway_analyze.retropath_promotion.MnxrefIndex",
                return_value=_FakeMnxref(),
            ),
            patch(
                "src.pathway_analyze.retropath_promotion.KeggRestClient",
                return_value=_FakeKegg(),
            ),
        ):
            materialize_retropath_solutions(
                self.config,
                validation_manifest_path=self.validation_manifest,
            )
        self.config.solution = 2
        self.config.manifest_output_path = self.root / "design_manifest.json"
        result = write_solution(self.config)
        self.assertTrue(result["运行成功"])
        manifest = json.loads(
            self.config.manifest_output_path.read_text(encoding="utf-8")
        )
        solution = manifest["solution"]
        self.assertEqual("retropath", solution["prediction"]["source"])
        self.assertEqual("pending", solution["prediction"]["review_status"])
        self.assertEqual(self.hypothesis_id, solution["steps"][0]["reaction_id"])
        self.assertEqual(
            self.hypothesis_id,
            solution["steps"][0]["retropath_hypothesis_id"],
        )
        self.assertTrue(solution["steps"][0]["formal_mapping_exact"])


if __name__ == "__main__":
    unittest.main()
