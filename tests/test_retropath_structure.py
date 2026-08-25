from __future__ import annotations

import hashlib
import json
import tempfile
import unittest
from pathlib import Path
from urllib.error import HTTPError, URLError

import rdkit
from rdkit import Chem

from src.pathway_analyze.retropath_structure import (
    KeggMolStructureProvider,
    StructureResolutionError,
    compound_from_kegg_mol,
)


def mol_block(smiles: str) -> str:
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise AssertionError(f"invalid test SMILES: {smiles}")
    return Chem.MolToMolBlock(molecule)


class RetroPathStructureTests(unittest.TestCase):
    def test_water_mol_generates_complete_standard_structure(self) -> None:
        compound = compound_from_kegg_mol(
            "c00001",
            mol_block("O"),
            minimum_depth=0,
            source_url="https://example.test/C00001/mol",
        )

        self.assertEqual(compound.compound_id, "C00001")
        self.assertEqual(compound.inchi, "InChI=1S/H2O/h1H2")
        self.assertEqual(compound.inchikey, "XLYOFNOQVPJJNP-UHFFFAOYSA-N")
        self.assertEqual(compound.isomeric_smiles, "O")
        self.assertEqual(compound.formula, "H2O")
        self.assertEqual(compound.charge, 0)
        self.assertEqual(compound.minimum_depth, 0)
        self.assertIn(f"rdkit:{rdkit.__version__}", compound.structure_provenance)
        self.assertTrue(
            any(
                value.startswith("kegg_mol_sha256:")
                for value in compound.structure_provenance
            )
        )

    def test_isomeric_smiles_preserves_stereochemistry(self) -> None:
        compound = compound_from_kegg_mol(
            "C00186",
            mol_block("C[C@H](O)C(=O)O"),
        )

        self.assertIn("@", compound.isomeric_smiles or "")
        self.assertEqual(compound.formula, "C3H6O3")

    def test_line_endings_do_not_change_mol_identity(self) -> None:
        unix_mol = mol_block("CCO").replace("\r\n", "\n")
        windows_mol = unix_mol.replace("\n", "\r\n")

        first = compound_from_kegg_mol("C00001", unix_mol)
        second = compound_from_kegg_mol("C00001", windows_mol)

        self.assertEqual(first, second)

    def test_invalid_mol_and_checksum_are_rejected(self) -> None:
        with self.assertRaisesRegex(StructureResolutionError, "mol_parse_failed"):
            compound_from_kegg_mol("C00001", "not a mol block")
        with self.assertRaisesRegex(StructureResolutionError, "mol_checksum_mismatch"):
            compound_from_kegg_mol(
                "C00001",
                mol_block("O"),
                mol_sha256="0" * 64,
            )
        with self.assertRaisesRegex(
            StructureResolutionError,
            "invalid_kegg_compound_id",
        ):
            compound_from_kegg_mol("not-kegg", mol_block("O"))

    def test_provider_uses_checksummed_cache(self) -> None:
        calls: list[str] = []

        def fetch(url: str, timeout: float) -> str:
            calls.append(f"{url}|{timeout}")
            return mol_block("O")

        with tempfile.TemporaryDirectory() as temporary_directory:
            provider = KeggMolStructureProvider(
                temporary_directory,
                retries=1,
                request_sleep_seconds=0,
                fetch_text=fetch,
            )
            first = provider.resolve("C00001", minimum_depth=0)
            second = provider.resolve("C00001", minimum_depth=2)
            cache_dir = Path(temporary_directory) / "mol"
            metadata = json.loads(
                (cache_dir / "C00001.json").read_text(encoding="utf-8")
            )
            cached_mol = (cache_dir / "C00001.mol").read_text(encoding="utf-8")

        self.assertEqual(len(calls), 1)
        self.assertEqual(first.inchi, second.inchi)
        self.assertEqual(first.minimum_depth, 0)
        self.assertEqual(second.minimum_depth, 2)
        self.assertEqual(metadata["kegg_id"], "C00001")
        normalized = cached_mol.replace("\r\n", "\n").replace("\r", "\n")
        self.assertEqual(
            metadata["mol_sha256"],
            hashlib.sha256(normalized.encode("utf-8")).hexdigest(),
        )

    def test_provider_refetches_when_cache_checksum_is_stale(self) -> None:
        calls = 0

        def fetch(url: str, timeout: float) -> str:
            nonlocal calls
            calls += 1
            return mol_block("O")

        with tempfile.TemporaryDirectory() as temporary_directory:
            provider = KeggMolStructureProvider(
                temporary_directory,
                retries=1,
                request_sleep_seconds=0,
                fetch_text=fetch,
            )
            provider.resolve("C00001")
            (Path(temporary_directory) / "mol" / "C00001.mol").write_text(
                mol_block("CCO"),
                encoding="utf-8",
            )
            refreshed = provider.resolve("C00001")

        self.assertEqual(calls, 2)
        self.assertEqual(refreshed.inchi, "InChI=1S/H2O/h1H2")

    def test_provider_does_not_cache_unparseable_response(self) -> None:
        def invalid(url: str, timeout: float) -> str:
            return "not a mol block"

        with tempfile.TemporaryDirectory() as temporary_directory:
            provider = KeggMolStructureProvider(
                temporary_directory,
                retries=1,
                request_sleep_seconds=0,
                fetch_text=invalid,
            )
            with self.assertRaisesRegex(
                StructureResolutionError,
                "mol_parse_failed",
            ):
                provider.resolve("C00001")
            cache_files = list((Path(temporary_directory) / "mol").glob("*"))

        self.assertEqual(cache_files, [])

    def test_provider_reports_network_and_empty_response_failures(self) -> None:
        def unavailable(url: str, timeout: float) -> str:
            raise URLError("offline")

        def empty(url: str, timeout: float) -> str:
            return ""

        def missing(url: str, timeout: float) -> str:
            raise HTTPError(url, 404, "not found", None, None)

        with tempfile.TemporaryDirectory() as first_directory:
            provider = KeggMolStructureProvider(
                first_directory,
                retries=1,
                request_sleep_seconds=0,
                fetch_text=unavailable,
            )
            with self.assertRaisesRegex(
                StructureResolutionError,
                "kegg_fetch_failed",
            ):
                provider.resolve("C00001")
        with tempfile.TemporaryDirectory() as second_directory:
            provider = KeggMolStructureProvider(
                second_directory,
                retries=1,
                request_sleep_seconds=0,
                fetch_text=empty,
            )
            with self.assertRaisesRegex(
                StructureResolutionError,
                "empty_mol_response",
            ):
                provider.resolve("C00001")
        with tempfile.TemporaryDirectory() as third_directory:
            provider = KeggMolStructureProvider(
                third_directory,
                retries=3,
                request_sleep_seconds=0,
                fetch_text=missing,
            )
            with self.assertRaisesRegex(
                StructureResolutionError,
                "kegg_structure_not_found",
            ):
                provider.resolve("C00001")

    def test_provider_configuration_is_validated(self) -> None:
        invalid_factories = (
            lambda: KeggMolStructureProvider("cache", timeout_seconds=0),
            lambda: KeggMolStructureProvider("cache", retries=0),
            lambda: KeggMolStructureProvider("cache", request_sleep_seconds=-1),
        )
        for factory in invalid_factories:
            with self.subTest(factory=factory), self.assertRaises(ValueError):
                factory()


if __name__ == "__main__":
    unittest.main()
