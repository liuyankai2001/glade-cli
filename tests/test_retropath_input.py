from __future__ import annotations

import csv
import hashlib
import tempfile
import unittest
from pathlib import Path
from typing import Optional

from services.retropath.app.validation import validate_compound_csv
from src.pathway_analyze.expand_chassis_metabolites import ExpansionBundle
from src.pathway_analyze.retropath_input import (
    RetroPathInputBuildError,
    build_retropath_inputs,
)
from src.pathway_analyze.retropath_models import PredictedCompound
from src.pathway_analyze.retropath_structure import StructureResolutionError


STRUCTURES = {
    "C00001": {
        "inchi": "InChI=1S/H2O/h1H2",
        "inchikey": "XLYOFNOQVPJJNP-UHFFFAOYSA-N",
        "smiles": "O",
        "formula": "H2O",
        "charge": 0,
    },
    # Deliberate KEGG alias for the same test structure.
    "C00002": {
        "inchi": "InChI=1S/H2O/h1H2",
        "inchikey": "XLYOFNOQVPJJNP-UHFFFAOYSA-N",
        "smiles": "O",
        "formula": "H2O",
        "charge": 0,
    },
    "C00003": {
        "inchi": "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3",
        "inchikey": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
        "smiles": "CCO",
        "formula": "C2H6O",
        "charge": 0,
    },
    "C00004": {
        "inchi": "InChI=1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)",
        "inchikey": "QTBSBXVTEAMEQO-UHFFFAOYSA-N",
        "smiles": "CC(=O)O",
        "formula": "C2H4O2",
        "charge": 0,
    },
    # Deliberately conflicting structure with C00001's InChIKey.
    "C00005": {
        "inchi": "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3",
        "inchikey": "XLYOFNOQVPJJNP-UHFFFAOYSA-N",
        "smiles": "CCO",
        "formula": "C2H6O",
        "charge": 0,
    },
    "C00999": {
        "inchi": "InChI=1S/CH4/h1H4",
        "inchikey": "VNWKTOKETHGBQD-UHFFFAOYSA-N",
        "smiles": "C",
        "formula": "CH4",
        "charge": 0,
    },
}


class FakeStructureProvider:
    def __init__(
        self,
        *,
        failures: Optional[set[str]] = None,
        incomplete: Optional[set[str]] = None,
    ) -> None:
        self.failures = failures or set()
        self.incomplete = incomplete or set()
        self.calls: list[tuple[str, Optional[int]]] = []

    def resolve(
        self,
        compound_id: str,
        *,
        minimum_depth: Optional[int] = None,
    ) -> PredictedCompound:
        self.calls.append((compound_id, minimum_depth))
        if compound_id in self.failures or compound_id not in STRUCTURES:
            raise StructureResolutionError(
                "kegg_structure_not_found",
                compound_id,
                "test structure is unavailable",
            )
        structure = STRUCTURES[compound_id]
        complete = compound_id not in self.incomplete
        return PredictedCompound.create(
            compound_id=compound_id,
            inchi=structure["inchi"],
            inchikey=structure["inchikey"],
            isomeric_smiles=structure["smiles"] if complete else None,
            formula=structure["formula"] if complete else None,
            charge=structure["charge"] if complete else None,
            kegg_ids=(compound_id,),
            minimum_depth=minimum_depth,
            structure_provenance=(f"fake:{compound_id}", "rdkit:test"),
        )


def expansion_bundle(depth_by_compound: dict[str, int], depth: int) -> ExpansionBundle:
    return ExpansionBundle(
        depth=depth,
        base_compounds=frozenset(
            compound_id
            for compound_id, minimum_depth in depth_by_compound.items()
            if minimum_depth == 0
        ),
        reachable_compounds=frozenset(depth_by_compound),
        depth_by_compound=dict(depth_by_compound),
        witnesses_by_product={},
        expanded_file=Path("unused.csv"),
        manifest={},
    )


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


class RetroPathInputTests(unittest.TestCase):
    def test_depth_zero_builds_service_compatible_source_and_sink(self) -> None:
        provider = FakeStructureProvider()
        bundle = expansion_bundle({"C00001": 0, "C00003": 0}, depth=0)

        with tempfile.TemporaryDirectory() as temporary_directory:
            result = build_retropath_inputs(
                "C00999",
                bundle,
                provider,
                temporary_directory,
            )
            source_bytes = result.target_source_path.read_bytes()
            sink_bytes = result.chassis_sink_path.read_bytes()
            source_rows = read_rows(result.target_source_path)
            sink_rows = read_rows(result.chassis_sink_path)

        self.assertEqual(result.expansion_depth, 0)
        self.assertEqual(result.reachable_compound_count, 2)
        self.assertEqual(result.sink_structure_count, 2)
        self.assertEqual(source_rows[0]["Name"], "C00999")
        self.assertEqual({row["Name"] for row in sink_rows}, {"C00001", "C00003"})
        self.assertEqual(validate_compound_csv(source_bytes, kind="source").row_count, 1)
        self.assertEqual(validate_compound_csv(sink_bytes, kind="sink").row_count, 2)
        self.assertEqual(
            result.target_source_sha256,
            hashlib.sha256(source_bytes).hexdigest(),
        )
        self.assertEqual(
            result.chassis_sink_sha256,
            hashlib.sha256(sink_bytes).hexdigest(),
        )
        self.assertTrue(source_bytes.startswith(b"Name,InChI\n"))

    def test_cumulative_sink_deduplicates_structure_and_uses_lowest_depth_alias(self) -> None:
        provider = FakeStructureProvider()
        bundle = expansion_bundle(
            {"C00001": 2, "C00002": 0, "C00003": 1},
            depth=2,
        )

        with tempfile.TemporaryDirectory() as temporary_directory:
            result = build_retropath_inputs(
                "C00999",
                bundle,
                provider,
                temporary_directory,
            )
            sink_rows = read_rows(result.chassis_sink_path)
            mapping_rows = read_rows(result.compound_mapping_path)

        self.assertEqual({row["Name"] for row in sink_rows}, {"C00002", "C00003"})
        water = next(item for item in result.sink_compounds if item.compound_id == "C00002")
        self.assertEqual(water.kegg_ids, ("C00001", "C00002"))
        self.assertEqual(water.minimum_depth, 0)
        sink_mapping = {
            row["kegg_id"]: row
            for row in mapping_rows
            if row["role"] == "sink"
        }
        self.assertEqual(sink_mapping["C00001"]["representative_kegg_id"], "C00002")
        self.assertEqual(sink_mapping["C00001"]["is_representative"], "false")
        self.assertEqual(sink_mapping["C00002"]["is_representative"], "true")

    def test_same_depth_alias_uses_lexical_kegg_id(self) -> None:
        bundle = expansion_bundle({"C00002": 1, "C00001": 1}, depth=1)

        with tempfile.TemporaryDirectory() as temporary_directory:
            result = build_retropath_inputs(
                "C00999",
                bundle,
                FakeStructureProvider(),
                temporary_directory,
            )

        self.assertEqual(result.sink_compounds[0].compound_id, "C00001")

    def test_partial_sink_failure_is_audited_without_blocking_valid_sink(self) -> None:
        bundle = expansion_bundle({"C00001": 0, "C00004": 1}, depth=1)

        with tempfile.TemporaryDirectory() as temporary_directory:
            result = build_retropath_inputs(
                "C00999",
                bundle,
                FakeStructureProvider(failures={"C00004"}),
                temporary_directory,
            )
            rejected_rows = read_rows(result.rejected_compounds_path)

        self.assertEqual(result.sink_structure_count, 1)
        self.assertEqual(result.rejected_compound_count, 1)
        self.assertEqual(rejected_rows[0]["kegg_id"], "C00004")
        self.assertEqual(rejected_rows[0]["reason_code"], "kegg_structure_not_found")

    def test_incomplete_structure_is_rejected(self) -> None:
        bundle = expansion_bundle({"C00001": 0, "C00004": 1}, depth=1)

        with tempfile.TemporaryDirectory() as temporary_directory:
            result = build_retropath_inputs(
                "C00999",
                bundle,
                FakeStructureProvider(incomplete={"C00004"}),
                temporary_directory,
            )

        self.assertEqual(result.rejected_compounds[0].reason_code, "structure_fields_incomplete")

    def test_target_failure_writes_audit_artifacts_and_blocks_build(self) -> None:
        bundle = expansion_bundle({"C00001": 0}, depth=0)

        with tempfile.TemporaryDirectory() as temporary_directory:
            output_dir = Path(temporary_directory)
            with self.assertRaises(RetroPathInputBuildError) as caught:
                build_retropath_inputs(
                    "C00999",
                    bundle,
                    FakeStructureProvider(failures={"C00999"}),
                    output_dir,
                )
            source_rows = read_rows(output_dir / "target_source.csv")
            rejected_rows = read_rows(output_dir / "rejected_compounds.csv")

        self.assertEqual(caught.exception.code, "target_structure_invalid")
        self.assertEqual(source_rows, [])
        self.assertEqual(rejected_rows[0]["role"], "target")

    def test_all_sink_failures_write_audit_and_block_build(self) -> None:
        bundle = expansion_bundle({"C00004": 0}, depth=0)

        with tempfile.TemporaryDirectory() as temporary_directory:
            output_dir = Path(temporary_directory)
            with self.assertRaises(RetroPathInputBuildError) as caught:
                build_retropath_inputs(
                    "C00999",
                    bundle,
                    FakeStructureProvider(failures={"C00004"}),
                    output_dir,
                )
            source_rows = read_rows(output_dir / "target_source.csv")
            sink_rows = read_rows(output_dir / "chassis_sink.csv")

        self.assertEqual(caught.exception.code, "sink_structure_empty")
        self.assertEqual(len(source_rows), 1)
        self.assertEqual(sink_rows, [])

    def test_inchikey_to_multiple_inchis_rejects_conflicting_group(self) -> None:
        bundle = expansion_bundle(
            {"C00001": 0, "C00005": 1, "C00003": 0},
            depth=1,
        )

        with tempfile.TemporaryDirectory() as temporary_directory:
            result = build_retropath_inputs(
                "C00999",
                bundle,
                FakeStructureProvider(),
                temporary_directory,
            )

        self.assertEqual(tuple(item.compound_id for item in result.sink_compounds), ("C00003",))
        conflicts = [
            item
            for item in result.rejected_compounds
            if item.reason_code == "structure_identity_conflict"
        ]
        self.assertEqual({item.kegg_id for item in conflicts}, {"C00001", "C00005"})

    def test_missing_minimum_depth_is_audited(self) -> None:
        bundle = ExpansionBundle(
            depth=1,
            base_compounds=frozenset({"C00003"}),
            reachable_compounds=frozenset({"C00003", "C00004"}),
            depth_by_compound={"C00003": 0},
            witnesses_by_product={},
            expanded_file=Path("unused.csv"),
            manifest={},
        )

        with tempfile.TemporaryDirectory() as temporary_directory:
            result = build_retropath_inputs(
                "C00999",
                bundle,
                FakeStructureProvider(),
                temporary_directory,
            )

        self.assertEqual(result.rejected_compounds[0].reason_code, "missing_minimum_depth")

    def test_depth_beyond_requested_cumulative_layer_is_audited(self) -> None:
        bundle = expansion_bundle({"C00003": 0, "C00004": 2}, depth=1)

        with tempfile.TemporaryDirectory() as temporary_directory:
            result = build_retropath_inputs(
                "C00999",
                bundle,
                FakeStructureProvider(),
                temporary_directory,
            )

        self.assertEqual(
            result.rejected_compounds[0].reason_code,
            "minimum_depth_exceeds_requested_depth",
        )

    def test_target_in_sink_is_preserved_for_source_in_sink_status(self) -> None:
        bundle = expansion_bundle({"C00999": 0, "C00001": 0}, depth=0)

        with tempfile.TemporaryDirectory() as temporary_directory:
            result = build_retropath_inputs(
                "C00999",
                bundle,
                FakeStructureProvider(),
                temporary_directory,
            )
            sink_names = {
                row["Name"] for row in read_rows(result.chassis_sink_path)
            }

        self.assertIn("C00999", sink_names)

    def test_output_bytes_and_hashes_are_stable(self) -> None:
        first_bundle = expansion_bundle(
            {"C00003": 1, "C00002": 0, "C00001": 2},
            depth=2,
        )
        second_bundle = expansion_bundle(
            {"C00001": 2, "C00002": 0, "C00003": 1},
            depth=2,
        )

        with tempfile.TemporaryDirectory() as first_directory:
            first = build_retropath_inputs(
                "C00999",
                first_bundle,
                FakeStructureProvider(),
                first_directory,
            )
            first_bytes = {
                path.name: path.read_bytes()
                for path in (
                    first.target_source_path,
                    first.chassis_sink_path,
                    first.compound_mapping_path,
                    first.rejected_compounds_path,
                )
            }
        with tempfile.TemporaryDirectory() as second_directory:
            second = build_retropath_inputs(
                "C00999",
                second_bundle,
                FakeStructureProvider(),
                second_directory,
            )
            second_bytes = {
                path.name: path.read_bytes()
                for path in (
                    second.target_source_path,
                    second.chassis_sink_path,
                    second.compound_mapping_path,
                    second.rejected_compounds_path,
                )
            }

        self.assertEqual(first_bytes, second_bytes)
        self.assertEqual(first.target_source_sha256, second.target_source_sha256)
        self.assertEqual(first.chassis_sink_sha256, second.chassis_sink_sha256)


if __name__ == "__main__":
    unittest.main()
