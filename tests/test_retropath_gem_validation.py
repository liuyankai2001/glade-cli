from __future__ import annotations

import csv
import hashlib
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import cobra

from src.cli.app import build_parser
from src.pathway_analyze.gem_validation import KeggRestClient, run_validation
from src.pathway_analyze.retropath_analyze import (
    CANDIDATE_ROUTE_COLUMNS,
    CANDIDATE_STEP_COLUMNS,
    REJECTED_ROUTE_COLUMNS,
)
from src.pathway_analyze.retropath_gem_validation import (
    STOICHIOMETRY_TERMS_FILE_NAME,
    _validate_combination,
    validate_retropath_candidates,
)
from src.pathway_analyze.retropath_pipeline import RETROPATH_PIPELINE_SCHEMA
from src.pathway_analyze.retropath_stoichiometry import (
    CompletedReactionHypothesis,
    CompletedTerm,
    CompoundProperty,
)


def _property(compound_id: str, formula: str) -> CompoundProperty:
    return CompoundProperty(
        compound_id=compound_id,
        name=compound_id,
        formula=formula,
        charge=0,
        source="test",
    )


def _hypothesis(
    candidate_id: str,
    step_id: str,
    left: list[CompoundProperty],
    right: list[CompoundProperty],
) -> CompletedReactionHypothesis:
    terms = tuple(
        [CompletedTerm("left", 1.0, item, "predicted_core") for item in left]
        + [CompletedTerm("right", 1.0, item, "predicted_core") for item in right]
    )
    return CompletedReactionHypothesis(
        hypothesis_id="RP2STOICH:test",
        candidate_id=candidate_id,
        step_id=step_id,
        rule_id="RR-TEST",
        source_mnxr_id="MNXR1",
        source_reference="fixture",
        source_equation="fixture",
        source_orientation="left_to_right",
        evidence_grade="rr02_mnxref_v3_template_balanced",
        terms=terms,
        balance_status="balanced",
        cofactor_reconstruction_status="not_applicable",
    )


def _metabolite(
    identifier: str,
    formula: str,
    *,
    kegg_id: str | None = None,
) -> cobra.Metabolite:
    metabolite = cobra.Metabolite(
        identifier,
        formula=formula,
        charge=0,
        compartment="c",
    )
    if kegg_id:
        metabolite.annotation["kegg.compound"] = kegg_id
    return metabolite


def _source(model: cobra.Model, metabolite: cobra.Metabolite, bound: float = 10) -> None:
    reaction = cobra.Reaction(f"SRC_{metabolite.id}", lower_bound=0, upper_bound=bound)
    reaction.add_metabolites({metabolite: 1})
    model.add_reactions([reaction])


class RetroPathStrictGemTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.kegg_client = KeggRestClient(self.root / "kegg", request_sleep_seconds=0)

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def test_balanced_forced_candidate_route_passes_with_growth(self) -> None:
        model = cobra.Model("pass")
        substrate = _metabolite("a_c", "C", kegg_id="C00010")
        target = _metabolite("t_c", "C", kegg_id="C12345")
        model.add_metabolites([substrate, target])
        _source(model, substrate)
        biomass = cobra.Reaction("BIOMASS", lower_bound=0, upper_bound=1000)
        biomass.add_metabolites({substrate: -1})
        model.add_reactions([biomass])
        model.objective = biomass
        candidate_id = "RP2ROUTE:test"
        step_id = "RP2STEP:test"
        substrate_property = _property("C00010", "C")
        target_property = _property("C12345", "C")
        hypothesis = _hypothesis(
            candidate_id,
            step_id,
            [substrate_property],
            [target_property],
        )

        result = _validate_combination(
            base_model=model,
            target="C12345",
            candidate_rank=1,
            candidate_id=candidate_id,
            steps=(
                {
                    "step_index": "1",
                    "step_id": step_id,
                    "step_source": "retropath",
                },
            ),
            hypotheses=(hypothesis,),
            known_properties={"C00010": substrate_property, "C12345": target_property},
            kegg_client=self.kegg_client,
            baseline_growth=10.0,
            required_growth=1.0,
            biomass_reaction_id="BIOMASS",
            minimum_flux=1e-4,
            run_fva=False,
            combination_truncated=False,
        )

        self.assertEqual("PASS_STRICT_ROUTE_FLUX", result.row["validation_status"])
        self.assertEqual(1, result.row["active_route_step_count"])
        self.assertGreaterEqual(float(result.row["target_flux"]), 1e-4)

    def test_target_bypass_cannot_hide_blocked_rp2_cofactor(self) -> None:
        model = cobra.Model("blocked")
        growth_source = _metabolite("g_c", "C")
        substrate = _metabolite("a_c", "C", kegg_id="C00010")
        target = _metabolite("t_c", "CH2", kegg_id="C12345")
        model.add_metabolites([growth_source, substrate, target])
        _source(model, growth_source)
        _source(model, substrate)
        biomass = cobra.Reaction("BIOMASS", lower_bound=0, upper_bound=1000)
        biomass.add_metabolites({growth_source: -1})
        bypass = cobra.Reaction("TARGET_BYPASS", lower_bound=0, upper_bound=1000)
        bypass.add_metabolites({substrate: -1, target: 1})
        model.add_reactions([biomass, bypass])
        model.objective = biomass
        candidate_id = "RP2ROUTE:test"
        step_id = "RP2STEP:test"
        substrate_property = _property("C00010", "C")
        cofactor_property = CompoundProperty(
            compound_id="MNXM999",
            name="missing cofactor",
            formula="H2",
            charge=0,
            source="mnxref:3.0",
            source_mnxm_id="MNXM999",
        )
        target_property = _property("C12345", "CH2")
        hypothesis = _hypothesis(
            candidate_id,
            step_id,
            [substrate_property, cofactor_property],
            [target_property],
        )

        result = _validate_combination(
            base_model=model,
            target="C12345",
            candidate_rank=1,
            candidate_id=candidate_id,
            steps=(
                {
                    "step_index": "1",
                    "step_id": step_id,
                    "step_source": "retropath",
                },
            ),
            hypotheses=(hypothesis,),
            known_properties={"C00010": substrate_property, "C12345": target_property},
            kegg_client=self.kegg_client,
            baseline_growth=10.0,
            required_growth=1.0,
            biomass_reaction_id="BIOMASS",
            minimum_flux=1e-4,
            run_fva=False,
            combination_truncated=False,
        )

        self.assertEqual("FAIL_GEM_INFEASIBLE", result.row["validation_status"])
        self.assertEqual(0, result.row["active_route_step_count"])

    def test_cli_has_independent_candidate_selector(self) -> None:
        parser = build_parser()
        all_candidates = parser.parse_args(
            ["validate", "--input", "example.json", "--retropath-candidates"]
        )
        self.assertEqual([], all_candidates.retropath_candidates)
        selected = parser.parse_args(
            [
                "validate",
                "--input",
                "example.json",
                "--retropath-candidates",
                "1",
                "2",
                "--depth",
                "3",
            ]
        )
        self.assertEqual([1, 2], selected.retropath_candidates)
        self.assertEqual(3, selected.depth)

    def test_run_validation_dispatches_without_touching_kegg_default(self) -> None:
        config = SimpleNamespace(retropath_candidates=[])
        expected = {"ok": True, "search_engine": "retropath"}
        with (
            patch(
                "src.pathway_analyze.retropath_gem_validation.validate_retropath_candidates",
                return_value=expected,
            ) as validate_retropath,
            patch("builtins.print"),
        ):
            self.assertIs(run_validation(config), expected)
        validate_retropath.assert_called_once_with(config)

        config = SimpleNamespace(retropath_candidates=None)
        kegg_expected = {"ok": True, "search_engine": "kegg"}
        with (
            patch(
                "src.pathway_analyze.gem_validation.gem_validate",
                return_value=kegg_expected,
            ) as validate_kegg,
            patch("builtins.print"),
        ):
            self.assertIs(run_validation(config), kegg_expected)
        validate_kegg.assert_called_once_with(config)

    @staticmethod
    def _write_csv(path: Path, columns: tuple[str, ...], rows: list[dict]) -> str:
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns, lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        return hashlib.sha256(path.read_bytes()).hexdigest()

    def test_full_orchestrator_writes_hashed_isolated_outputs(self) -> None:
        target_id = "C12345"
        substrate_id = "C00010"
        candidate_id = "RP2ROUTE:" + "1" * 64
        step_id = "RP2STEP:" + "2" * 64
        retropath_dir = self.root / "gap" / "depth0" / "retropath"
        input_dir = retropath_dir / "input"
        input_dir.mkdir(parents=True)
        route_path = retropath_dir / "candidate_routes.csv"
        step_path = retropath_dir / "candidate_steps.csv"
        rejected_path = retropath_dir / "rejected_routes.csv"
        route_hash = self._write_csv(
            route_path,
            CANDIDATE_ROUTE_COLUMNS,
            [
                {
                    "candidate_rank": 1,
                    "candidate_id": candidate_id,
                    "source_retrosynthetic_path_id": "RP2PATH:test",
                    "target_compound_id": target_id,
                    "sink_kegg_ids": substrate_id,
                    "sink_depths": f"{substrate_id}:0",
                    "sink_inchikeys": "TEST",
                    "kegg_prefix_reaction_ids": "",
                    "retropath_step_ids": step_id,
                    "retropath_reaction_option_ids": "RP2:" + "3" * 64,
                    "kegg_prefix_steps": 0,
                    "retropath_steps": 1,
                    "total_steps": 1,
                    "maximum_sink_depth": 0,
                    "minimum_rule_specificity": 8,
                    "worst_rule_score": 0.1,
                    "score_semantics": "lower_is_better",
                    "contains_auxiliary_fragments": "false",
                    "route_source": "kegg_retropath",
                    "contains_predicted_steps": "true",
                    "validation_status": "raw",
                    "review_required": "true",
                    "upstream_enumeration_truncated": "false",
                    "candidate_top_k_truncated": "false",
                }
            ],
        )
        step_hash = self._write_csv(
            step_path,
            CANDIDATE_STEP_COLUMNS,
            [
                {
                    "candidate_id": candidate_id,
                    "step_index": 1,
                    "step_id": step_id,
                    "step_source": "retropath",
                    "status": "predicted",
                    "orientation": "biosynthetic",
                    "direction": "biosynthetic",
                    "reaction_option_ids": "RP2:" + "3" * 64,
                    "reaction_smiles": "[C]>>[C]",
                    "substrate_compound_ids": substrate_id,
                    "product_compound_ids": target_id,
                    "substrate_stoichiometry_json": f'[["{substrate_id}",1.0]]',
                    "product_stoichiometry_json": f'[["{target_id}",1.0]]',
                    "depends_on_step_ids": "",
                    "source_transformation_ids": "TRS-1",
                    "sink_anchor_kegg_ids": substrate_id,
                    "expansion_depth": 0,
                    "is_endogenous": "",
                    "rule_ids": "RR-TEST",
                    "source_reaction_ids": "MNXR1",
                    "source_ec_numbers": "1.1.1.1",
                    "minimum_rule_specificity": 8,
                    "worst_rule_score": 0.1,
                    "score_semantics": "lower_is_better",
                    "balance_status": "not_checked",
                    "cofactor_reconstruction_status": "not_checked",
                }
            ],
        )
        rejected_hash = self._write_csv(
            rejected_path,
            REJECTED_ROUTE_COLUMNS,
            [],
        )
        mapping_path = input_dir / "compound_mapping.csv"
        mapping_columns = (
            "role",
            "kegg_id",
            "representative_kegg_id",
            "is_representative",
            "minimum_depth",
            "inchi",
            "inchikey",
            "isomeric_smiles",
            "formula",
            "charge",
            "structure_provenance",
        )
        mapping_hash = self._write_csv(
            mapping_path,
            mapping_columns,
            [
                {
                    "role": "target",
                    "kegg_id": target_id,
                    "representative_kegg_id": target_id,
                    "is_representative": "true",
                    "minimum_depth": "",
                    "inchi": "InChI=1S/C",
                    "inchikey": "VNWKTOKETHGBQD-UHFFFAOYSA-N",
                    "isomeric_smiles": "[C]",
                    "formula": "C",
                    "charge": 0,
                    "structure_provenance": "[]",
                },
                {
                    "role": "sink",
                    "kegg_id": substrate_id,
                    "representative_kegg_id": substrate_id,
                    "is_representative": "true",
                    "minimum_depth": 0,
                    "inchi": "InChI=1S/C",
                    "inchikey": "VNWKTOKETHGBQD-UHFFFAOYSA-N",
                    "isomeric_smiles": "[C]",
                    "formula": "C",
                    "charge": 0,
                    "structure_provenance": "[]",
                },
            ],
        )
        pipeline_path = retropath_dir / "pipeline_result.json"
        pipeline = {
            "schema_version": RETROPATH_PIPELINE_SCHEMA,
            "ok": True,
            "status": "retropath_candidates_found",
            "target_compound": target_id,
            "expansion_depth": 0,
            "candidate_count": 1,
            "artifacts": {
                "candidate_routes": {"path": str(route_path), "sha256": route_hash},
                "candidate_steps": {"path": str(step_path), "sha256": step_hash},
                "rejected_routes": {
                    "path": str(rejected_path),
                    "sha256": rejected_hash,
                },
                "compound_mapping": {
                    "path": str(mapping_path),
                    "sha256": mapping_hash,
                },
            },
        }
        pipeline_path.write_text(json.dumps(pipeline), encoding="utf-8")

        model = cobra.Model("integration")
        substrate = _metabolite("a_c", "C", kegg_id=substrate_id)
        target = _metabolite("t_c", "C", kegg_id=target_id)
        model.add_metabolites([substrate, target])
        _source(model, substrate)
        biomass = cobra.Reaction("BIOMASS", lower_bound=0, upper_bound=1000)
        biomass.add_metabolites({substrate: -1})
        model.add_reactions([biomass])
        model.objective = biomass
        model_path = self.root / "model.json"
        cobra.io.save_json_model(model, model_path)
        medium_path = self.root / "medium.json"
        medium_path.write_text(json.dumps({"SRC_a_c": 10}), encoding="utf-8")
        rules_path = self.root / "rules.csv"
        rules_path.write_text("Rule ID\nRR-TEST\n", encoding="utf-8")
        hypothesis = _hypothesis(
            candidate_id,
            step_id,
            [_property(substrate_id, "C")],
            [_property(target_id, "C")],
        )
        reconstruction = SimpleNamespace(
            candidate_id=candidate_id,
            step_id=step_id,
            status="complete",
            hypotheses=(hypothesis,),
            rejections=tuple(),
            truncated=False,
        )
        fake_index = MagicMock()
        fake_index.__enter__.return_value = fake_index
        fake_index.__exit__.return_value = None
        fake_index.manifest = {
            "index_path": str(self.root / "index.sqlite3"),
            "index_sha256": "f" * 64,
        }
        config = SimpleNamespace(
            target_name=target_id,
            depth=0,
            gap_output_path=self.root / "gap",
            model_path=model_path,
            medium_path=medium_path,
            cache_dir=self.root / "cache",
            data_dir=self.root / "data",
            retropath_rules_path=rules_path,
            retropath_candidates=[],
            validation_mode="per",
            validation_cofactor_mode="strict",
            validation_skip_fva=True,
        )

        with (
            patch(
                "src.pathway_analyze.retropath_gem_validation.MnxrefIndex",
                return_value=fake_index,
            ),
            patch(
                "src.pathway_analyze.retropath_gem_validation.reconstruct_retropath_step",
                return_value=reconstruction,
            ),
            patch(
                "src.pathway_analyze.retropath_promotion.materialize_retropath_solutions",
                return_value={
                    "formal_solution_ids": [1],
                    "solution_mappings": [{"solution_id": 1}],
                    "promotion_manifest": str(
                        retropath_dir / "formal_solution_promotion.json"
                    ),
                    "promotion_manifest_sha256": "a" * 64,
                },
            ),
        ):
            result = validate_retropath_candidates(config)

        self.assertEqual(
            "PASS_STRICT_HYPOTHESIS_EXISTS",
            result["candidate_statuses"]["1"],
        )
        validation_dir = Path(result["validation_dir"])
        self.assertTrue((validation_dir / "validation_manifest.json").is_file())
        manifest = json.loads(
            (validation_dir / "validation_manifest.json").read_text(encoding="utf-8")
        )
        self.assertTrue(manifest["formal_promotion_allowed"])
        self.assertEqual([1], result["formal_solution_ids"])
        self.assertEqual(1, result["summary_row_count"])
        with (validation_dir / STOICHIOMETRY_TERMS_FILE_NAME).open(
            "r",
            encoding="utf-8",
            newline="",
        ) as handle:
            terms = list(csv.DictReader(handle))
        self.assertTrue(terms)
        self.assertTrue(all("smiles" in row for row in terms))


if __name__ == "__main__":
    unittest.main()
