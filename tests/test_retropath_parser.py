from __future__ import annotations

import csv
import hashlib
import json
import tempfile
import unittest
from pathlib import Path
from typing import Iterable, Mapping

import rdkit
from rdkit import Chem
from rdkit.Chem import inchi as rd_inchi
from rdkit.Chem import rdMolDescriptors

from src.pathway_analyze.retropath_client import RetroPathClientRun
from src.pathway_analyze.retropath_input import RetroPathInputBundle
from src.pathway_analyze.retropath_models import (
    PredictedCompound,
    RetroPathRunResult,
    RetroPathRuntimeProvenance,
)
from src.pathway_analyze.retropath_parser import (
    RETROPATH_RESULTS_COLUMNS,
    RetroPathParseError,
    parse_retropath_network,
)
from src.pathway_analyze.retropath_routes import parse_and_enumerate_retropath

RULE_COLUMNS = (
    "Rule ID",
    "Rule",
    "EC number",
    "Reaction order",
    "Diameter",
    "Score",
    "Legacy ID",
    "Reaction direction",
    "Rule relative direction",
    "Rule usage",
    "Score normalized",
)


def make_compound(
    kegg_id: str,
    smiles: str,
    *,
    minimum_depth: int | None = None,
) -> PredictedCompound:
    molecule = Chem.MolFromSmiles(smiles)
    assert molecule is not None
    inchi = rd_inchi.MolToInchi(molecule)
    return PredictedCompound.create(
        compound_id=kegg_id,
        inchi=inchi,
        inchikey=rd_inchi.InchiToInchiKey(inchi),
        isomeric_smiles=Chem.MolToSmiles(molecule, isomericSmiles=True),
        formula=rdMolDescriptors.CalcMolFormula(molecule),
        charge=sum(atom.GetFormalCharge() for atom in molecule.GetAtoms()),
        kegg_ids=(kegg_id,),
        minimum_depth=minimum_depth,
        structure_provenance=("test", f"rdkit:{rdkit.__version__}"),
    )


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_csv(
    path: Path,
    fieldnames: Iterable[str],
    rows: Iterable[Mapping[str, object]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=tuple(fieldnames), lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


class RetroPathParserTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.raw_dir = self.root / "run" / "raw"
        self.raw_dir.mkdir(parents=True)
        self.target = make_compound("C90000", "CCO")
        self.sink_methanol = make_compound("C00002", "CO", minimum_depth=0)
        self.sink_methane = make_compound("C00003", "C", minimum_depth=2)
        self.intermediate = make_compound("C90001", "CC=O")
        self.bundle = self.make_bundle(
            self.target,
            (self.sink_methanol, self.sink_methane),
        )
        self.rules_path = self.root / "rules.csv"
        write_csv(
            self.rules_path,
            RULE_COLUMNS,
            (
                self.rule_row(
                    "RULE-A",
                    diameter=8,
                    score=0.4,
                    ec="1.1.1.1",
                    legacy="MNXR100_MNXM1",
                ),
                self.rule_row(
                    "RULE-B", diameter=8, score=0.6, ec="NOEC", legacy="MNXR101_MNXM2"
                ),
                self.rule_row(
                    "RULE-C",
                    diameter=6,
                    score=0.8,
                    ec="2.2.2.-",
                    legacy="MNXR102_MNXM3",
                ),
                self.rule_row("UNUSED", diameter=2, score=2.0, ec="NOEC", legacy=""),
            ),
        )

    def tearDown(self) -> None:
        self.temporary.cleanup()

    @staticmethod
    def rule_row(
        rule_id: str,
        *,
        diameter: int,
        score: float,
        ec: str,
        legacy: str,
    ) -> dict[str, object]:
        return {
            "Rule ID": rule_id,
            "Rule": "[C:1]>>[C:1]",
            "EC number": ec,
            "Reaction order": 1,
            "Diameter": diameter,
            "Score": score,
            "Legacy ID": legacy,
            "Reaction direction": 0,
            "Rule relative direction": -1,
            "Rule usage": "retro",
            "Score normalized": 1.0 / (1.0 + score),
        }

    def make_bundle(
        self,
        target: PredictedCompound,
        sinks: tuple[PredictedCompound, ...],
    ) -> RetroPathInputBundle:
        input_dir = self.root / "input"
        input_dir.mkdir(exist_ok=True)
        source_path = input_dir / "target_source.csv"
        sink_path = input_dir / "chassis_sink.csv"
        source_path.write_text(
            f'Name,InChI\n{target.compound_id},"{target.inchi}"\n',
            encoding="utf-8",
        )
        sink_path.write_text(
            "Name,InChI\n"
            + "".join(f'{item.compound_id},"{item.inchi}"\n' for item in sinks),
            encoding="utf-8",
        )
        return RetroPathInputBundle(
            expansion_depth=2,
            reachable_compound_count=len(sinks),
            target_compound=target,
            sink_compounds=sinks,
            mappings=tuple(),
            rejected_compounds=tuple(),
            target_source_path=source_path,
            chassis_sink_path=sink_path,
            compound_mapping_path=input_dir / "compound_mapping.csv",
            rejected_compounds_path=input_dir / "rejected_compounds.csv",
            target_source_sha256=hashlib.sha256(source_path.read_bytes()).hexdigest(),
            chassis_sink_sha256=hashlib.sha256(sink_path.read_bytes()).hexdigest(),
        )

    def make_run(
        self,
        status: str = "succeeded",
        *,
        rules_sha256: str | None = None,
    ) -> RetroPathClientRun:
        return_code = {
            "succeeded": 0,
            "source_in_sink": 10,
            "no_solution": 11,
            "failed": 1,
            "timed_out": -15,
        }[status]
        result = RetroPathRunResult(
            job_id="rp2-parser-test",
            status=status,
            return_code=return_code,
            provenance=RetroPathRuntimeProvenance(
                wrapper_version="3.9.1",
                wrapper_reported_version="3.9.1",
                workflow_version="r20260212",
                knime_version="4.7.0",
                rdkit_plugin_version="4.9.1",
                rules_version="rr02-rp2-hs",
                rules_sha256=rules_sha256 or sha256_file(self.rules_path),
            ),
            parameters=(("max_steps", 3),),
            artifacts=("raw/target_scope.csv",),
        )
        run_dir = self.raw_dir.parent
        return RetroPathClientRun(
            result=result,
            request_fingerprint="f" * 64,
            output_dir=run_dir,
            raw_dir=self.raw_dir,
            run_manifest_path=run_dir / "run_manifest.json",
            client_state_path=run_dir / "client_state.json",
            cache_hit=False,
        )

    @staticmethod
    def result_row(
        transformation_id: str,
        reaction_smiles: str,
        substrate: PredictedCompound,
        product: PredictedCompound,
        *,
        iteration: int,
        diameter: int,
        rule_ids: str,
        ec_numbers: str,
        score: float,
        in_sink: bool,
        sink_name: str = "None",
    ) -> dict[str, object]:
        return {
            "Initial source": "[target]",
            "Transformation ID": transformation_id,
            "Reaction SMILES": reaction_smiles,
            "Substrate SMILES": substrate.isomeric_smiles,
            "Substrate InChI": substrate.inchi,
            "Product SMILES": product.isomeric_smiles,
            "Product InChI": product.inchi,
            "In Sink": int(in_sink),
            "Sink name": f"[{sink_name}]",
            "Diameter": diameter,
            "Rule ID": f"[{rule_ids}]",
            "EC number": f"[{ec_numbers}]",
            "Score": score,
            "Iteration": iteration,
        }

    def write_complete_scope(self) -> None:
        rows = (
            self.result_row(
                "TRS_0_0_1",
                "CCO>>CC=O.CO",
                self.target,
                self.intermediate,
                iteration=0,
                diameter=8,
                rule_ids="RULE-A, RULE-B",
                ec_numbers="1.1.1.1",
                score=0.4,
                in_sink=False,
            ),
            self.result_row(
                "TRS_0_0_1",
                "CCO>>CC=O.CO",
                self.target,
                self.sink_methanol,
                iteration=0,
                diameter=8,
                rule_ids="RULE-A, RULE-B",
                ec_numbers="1.1.1.1",
                score=0.4,
                in_sink=True,
                sink_name="C00002",
            ),
            self.result_row(
                "TRS_0_1_1",
                "CC=O>>C",
                self.intermediate,
                self.sink_methane,
                iteration=1,
                diameter=6,
                rule_ids="RULE-C",
                ec_numbers="2.2.2.-",
                score=0.8,
                in_sink=True,
                sink_name="C00003",
            ),
        )
        write_csv(self.raw_dir / "target_scope.csv", RETROPATH_RESULTS_COLUMNS, rows)

    def test_parses_transformations_rule_variants_and_exact_sink_matches(self) -> None:
        self.write_complete_scope()

        network = parse_retropath_network(
            self.make_run(),
            self.bundle,
            self.rules_path,
        )

        self.assertEqual(network.status, "succeeded")
        self.assertEqual(len(network.transformations), 2)
        self.assertEqual(len(network.transformations[0].reaction_variants), 2)
        self.assertEqual(network.transformations[0].score_semantics, "lower_is_better")
        self.assertEqual(
            network.transformations[0].reaction_variants[0].rule_specificity,
            8,
        )
        self.assertIn(
            "MNXR100",
            network.transformations[0].reaction_variants[0].source_reaction_ids,
        )
        self.assertEqual(
            {item.representative_kegg_id for item in network.sink_matches},
            {"C00002", "C00003"},
        )
        self.assertTrue(
            all(
                item.compound.compound_id.startswith("RP2CPD:")
                for item in network.compounds
                if not item.compound.kegg_ids
            )
        )
        self.assertFalse(network.rejections)

    def test_json_topology_is_cross_checked_against_csv(self) -> None:
        self.write_complete_scope()
        node_ids = {
            "target": self.target.inchikey,
            "intermediate": self.intermediate.inchikey,
            "methanol": self.sink_methanol.inchikey,
            "methane": self.sink_methane.inchikey,
        }
        payload = {
            "elements": {
                "nodes": [
                    {
                        "data": {
                            "id": node_ids[name],
                            "type": "compound",
                            "InChI": compound.inchi,
                            "SMILES": compound.isomeric_smiles,
                        }
                    }
                    for name, compound in (
                        ("target", self.target),
                        ("intermediate", self.intermediate),
                        ("methanol", self.sink_methanol),
                        ("methane", self.sink_methane),
                    )
                ]
                + [
                    {"data": {"id": "TRS_0_0_1", "type": "reaction"}},
                    {"data": {"id": "TRS_0_1_1", "type": "reaction"}},
                ],
                "edges": [
                    {"data": {"source": node_ids["target"], "target": "TRS_0_0_1"}},
                    {
                        "data": {
                            "source": "TRS_0_0_1",
                            "target": node_ids["intermediate"],
                        }
                    },
                    {"data": {"source": "TRS_0_0_1", "target": node_ids["methanol"]}},
                    {
                        "data": {
                            "source": node_ids["intermediate"],
                            "target": "TRS_0_1_1",
                        }
                    },
                    {"data": {"source": "TRS_0_1_1", "target": node_ids["target"]}},
                ],
            }
        }
        (self.raw_dir / "target_scope.json").write_text(
            json.dumps(payload),
            encoding="utf-8",
        )

        network = parse_retropath_network(
            self.make_run(),
            self.bundle,
            self.rules_path,
        )

        rejected = {
            item.transformation_id: item.reason_code for item in network.rejections
        }
        self.assertEqual(rejected["TRS_0_1_1"], "artifact_inconsistent")
        self.assertEqual(len(network.transformations), 1)

    def test_wrapper_sink_flag_cannot_replace_exact_p2_identity(self) -> None:
        write_csv(
            self.raw_dir / "target_scope.csv",
            RETROPATH_RESULTS_COLUMNS,
            (
                self.result_row(
                    "TRS_0_0_1",
                    "CCO>>CC=O",
                    self.target,
                    self.intermediate,
                    iteration=0,
                    diameter=8,
                    rule_ids="RULE-A",
                    ec_numbers="1.1.1.1",
                    score=0.4,
                    in_sink=True,
                    sink_name="C00002",
                ),
            ),
        )

        network = parse_retropath_network(
            self.make_run(),
            self.bundle,
            self.rules_path,
        )

        self.assertFalse(network.transformations)
        self.assertEqual(network.rejections[0].reason_code, "sink_identity_mismatch")

    def test_rules_checksum_mismatch_is_a_run_level_error(self) -> None:
        self.write_complete_scope()

        with self.assertRaisesRegex(RetroPathParseError, "rules_checksum_mismatch"):
            parse_retropath_network(
                self.make_run(rules_sha256="a" * 64),
                self.bundle,
                self.rules_path,
            )

    def test_parse_and_enumerate_returns_only_the_complete_branched_route(self) -> None:
        self.write_complete_scope()

        result = parse_and_enumerate_retropath(
            self.make_run(),
            self.bundle,
            self.rules_path,
        )

        self.assertEqual(result.complete_path_count, 1)
        self.assertEqual(
            {item.representative_kegg_id for item in result.paths[0].sink_matches},
            {"C00002", "C00003"},
        )

    def test_product_row_order_does_not_change_reaction_or_path_identity(self) -> None:
        self.write_complete_scope()
        first = parse_and_enumerate_retropath(
            self.make_run(),
            self.bundle,
            self.rules_path,
        )
        scope_path = self.raw_dir / "target_scope.csv"
        with scope_path.open("r", encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle))
        write_csv(scope_path, RETROPATH_RESULTS_COLUMNS, reversed(rows))

        second = parse_and_enumerate_retropath(
            self.make_run(),
            self.bundle,
            self.rules_path,
        )

        self.assertEqual(
            first.network.predicted_reactions,
            second.network.predicted_reactions,
        )
        self.assertEqual(first.paths[0].path_id, second.paths[0].path_id)

    def test_non_structural_auxiliary_fragment_is_audited_not_silently_lost(
        self,
    ) -> None:
        sink_row = self.result_row(
            "TRS_0_0_1",
            "CCO>>CO.[H+]",
            self.target,
            self.sink_methanol,
            iteration=0,
            diameter=8,
            rule_ids="RULE-A",
            ec_numbers="1.1.1.1",
            score=0.4,
            in_sink=False,
        )
        proton_row = dict(sink_row)
        proton_row.update(
            {
                "Product SMILES": "[H+]",
                "Product InChI": "InChI=1S/p+1",
            }
        )
        write_csv(
            self.raw_dir / "results.csv",
            RETROPATH_RESULTS_COLUMNS,
            (sink_row, proton_row),
        )

        network = parse_retropath_network(
            self.make_run(),
            self.bundle,
            self.rules_path,
        )

        transformation = network.transformations[0]
        self.assertEqual(len(transformation.auxiliary_fragments), 1)
        self.assertEqual(
            transformation.reaction_variants[0].cofactor_reconstruction_status,
            "incomplete",
        )
        self.assertIn(
            f"sink_flag_disagreement:{self.sink_methanol.compound_id}",
            network.warnings,
        )

    def test_terminal_statuses_do_not_create_fake_predictions(self) -> None:
        no_solution = parse_retropath_network(
            self.make_run("no_solution"),
            self.bundle,
            self.rules_path,
        )
        self.assertFalse(no_solution.transformations)

        target_as_sink = make_compound("C90000", "CCO", minimum_depth=0)
        source_bundle = self.make_bundle(self.target, (target_as_sink,))
        source_in_sink = parse_retropath_network(
            self.make_run("source_in_sink"),
            source_bundle,
            self.rules_path,
        )
        self.assertIsNotNone(source_in_sink.source_in_sink)
        self.assertFalse(source_in_sink.predicted_reactions)

        with self.assertRaisesRegex(RetroPathParseError, "execution_failed"):
            parse_retropath_network(
                self.make_run("failed"),
                self.bundle,
                self.rules_path,
            )


if __name__ == "__main__":
    unittest.main()
