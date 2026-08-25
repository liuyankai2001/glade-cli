from __future__ import annotations

import csv
import hashlib
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from src.cli.app import build_parser
from src.main_protein_selection.retropath_enzyme_selection import (
    ENZYME_CANDIDATE_COLUMNS,
    RETROPATH_ENZYME_SELECTION_SCHEMA,
    STEP_ENZYME_CANDIDATES_FILE_NAME,
    RetropathEnzymeRequirement,
    _build_requirements,
    _load_p9_inputs,
    _search_requirement,
    _SearchOutcome,
    _ValidatedP9Inputs,
    select_retropath_enzymes,
)
from src.main_protein_selection.select_main_enzymes import (
    run_main_protein_selection,
)
from src.main_protein_selection.selenzyme_retrieval import (
    SelenzymeClient,
    SelenzymeSourceUnavailable,
)
from src.main_protein_selection.uniprot_protein_candidates import ProteinCandidate
from src.pathway_analyze.retropath_analyze import (
    CANDIDATE_ROUTE_COLUMNS,
    CANDIDATE_ROUTES_FILE_NAME,
    CANDIDATE_STEP_COLUMNS,
    CANDIDATE_STEPS_FILE_NAME,
    REJECTED_ROUTE_COLUMNS,
    REJECTED_ROUTES_FILE_NAME,
)
from src.pathway_analyze.retropath_gem_validation import (
    HYPOTHESIS_COLUMNS,
    RETROPATH_GEM_VALIDATION_SCHEMA,
    STOICHIOMETRY_HYPOTHESES_FILE_NAME,
    STOICHIOMETRY_TERMS_FILE_NAME,
    SUMMARY_COLUMNS,
    TERM_COLUMNS,
    VALIDATION_MANIFEST_FILE_NAME,
    VALIDATION_SUMMARY_FILE_NAME,
)
from src.pathway_analyze.retropath_mnxref import (
    MnxrefChemical,
    MnxrefReactionTemplate,
    MnxrefReactionTerm,
)
from src.pathway_analyze.retropath_pipeline import (
    PIPELINE_RESULT_FILE_NAME,
    RETROPATH_PIPELINE_SCHEMA,
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_csv(path: Path, columns: tuple[str, ...], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in columns})


class _Response:
    status_code = 200
    text = "fixture-response"

    def raise_for_status(self) -> None:
        return None

    def json(self) -> dict:
        return {
            "app": "Selenzy",
            "version": "fixture",
            "data": json.dumps([
                {
                    "Seq. ID": "P12345",
                    "Score": 90.0,
                    "Rxn. ID": "MNXR1",
                    "EC Number": "1.1.1.1",
                    "Reaction similarity": 0.9,
                    "sim_RF": 0.8,
                }
            ]),
        }


class _Session:
    def __init__(self) -> None:
        self.calls: list[dict] = []

    def request(self, method: str, url: str, **kwargs):
        self.calls.append({"method": method, "url": url, **kwargs})
        return _Response()


class _FakeIndex:
    def __init__(self) -> None:
        self.template = MnxrefReactionTemplate(
            rule_id="RR-TEST",
            mnxr_id="MNXR1",
            main_mnxm_id="MNXM1",
            reaction_direction="0",
            rule_relative_direction="0",
            rule_usage="both",
            equation="1 MNXM1 = 1 MNXM2",
            balanced=True,
            transport=False,
            reference="fixture",
            parse_status="ok",
            terms=(
                MnxrefReactionTerm("left", 1.0, "MNXM1", "MNXD1", 0),
                MnxrefReactionTerm("right", 1.0, "MNXM2", "MNXD1", 0),
            ),
            reaction_xrefs=("kegg:R00001", "rhea:10000"),
        )
        self.values = {
            "MNXM1": MnxrefChemical(
                "MNXM1",
                "methane",
                "CH4",
                0,
                16.0,
                "InChI=1S/CH4/h1H4",
                "C",
                "fixture:methane",
                "VNWKTOKETHGBQD-UHFFFAOYSA-N",
            ),
            "MNXM2": MnxrefChemical(
                "MNXM2",
                "water",
                "H2O",
                0,
                18.0,
                "InChI=1S/H2O/h1H2",
                "O",
                "fixture:water",
                "XLYOFNOQVPJJNP-UHFFFAOYSA-N",
            ),
        }

    def templates_for_rules(self, _rule_ids):
        return (self.template,)

    def chemicals(self, mnxm_ids):
        return {value: self.values[value] for value in mnxm_ids}


def _protein(similarity: float = 1.0) -> ProteinCandidate:
    candidate = ProteinCandidate(
        accession="P12345",
        entry_name="TEST_ECOLI",
        protein_name="fixture enzyme",
        organism_name="Escherichia coli",
        organism_id=562,
        reviewed=True,
        length=300,
        ec_numbers=["1.1.1.1"],
        score=88.0,
        reasons=["fixture"],
        sequence="M" * 300,
        sequence_sha256="a" * 64,
    )
    candidate.retrieval_strategy = "selenzyme_full_reaction_smiles_exact"
    candidate.retrieval_query_id = "selenzyme_fixture"
    candidate.reaction_confidence = "selenzyme_structural_prediction"
    candidate.selenzyme_reaction_similarity = similarity
    candidate.selenzyme_sim_rf = 0.95
    candidate.selenzyme_matched_reaction_id = "MNXR1"
    return candidate


class RetroPathEnzymeSelectionTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def test_selenzyme_client_posts_structural_query_and_reuses_cache(self) -> None:
        session = _Session()
        client = SelenzymeClient(
            session=session,
            cache_root=self.root / "cache",
            rest_url="http://selenzyme.test/REST",
        )
        first = client.query_reaction_smarts(
            "CCO>>CC=O",
            query_kind="full_reaction_smiles",
            host_taxon_id=511145,
            targets=50,
        )
        second = client.query_reaction_smarts(
            "CCO>>CC=O",
            query_kind="full_reaction_smiles",
            host_taxon_id=511145,
            targets=50,
        )
        self.assertEqual(1, len(session.calls))
        payload = session.calls[0]["json"]
        self.assertEqual("CCO>>CC=O", payload["smarts"])
        self.assertEqual("Morgan", payload["fp"])
        self.assertNotIn("db", payload)
        self.assertNotIn("rxnid", payload)
        self.assertEqual("selenzyme_client.v4", first["cache_schema"])
        self.assertTrue(second["cache_hit"])
        self.assertEqual(0.9, first["rows"][0]["reaction_similarity"])

    def test_selenzyme_null_data_is_not_misreported_as_no_hit(self) -> None:
        class NullResponse(_Response):
            def json(self) -> dict:
                return {"app": "Selenzy", "version": "fixture", "data": None}

        class NullSession(_Session):
            def request(self, method: str, url: str, **kwargs):
                self.calls.append({"method": method, "url": url, **kwargs})
                return NullResponse()

        client = SelenzymeClient(
            session=NullSession(),
            cache_root=self.root / "cache",
            rest_url="http://selenzyme.test/REST",
        )
        with patch(
            "src.main_protein_selection.selenzyme_retrieval."
            "SELENZYME_HTTP_CONFIG",
            SimpleNamespace(retries=1, timeout_seconds=1.0, sleep_seconds=0.0),
        ):
            with self.assertRaisesRegex(
                SelenzymeSourceUnavailable,
                "response data is null",
            ):
                client.query_reaction_smarts(
                    "CCO>>CC=O",
                    query_kind="full_reaction_smiles",
                    host_taxon_id=511145,
                    targets=50,
                )

    def test_build_requirements_separates_hypotheses_and_exact_mapping(self) -> None:
        candidate_id = "RP2ROUTE:" + "1" * 64
        step_id = "RP2STEP:" + "2" * 64
        hypotheses = (
            {
                "candidate_rank": "1",
                "candidate_id": candidate_id,
                "step_id": step_id,
                "hypothesis_id": "RP2STOICH:" + "3" * 64,
                "rule_id": "RR-TEST",
                "source_mnxr_id": "MNXR1",
                "source_orientation": "left_to_right",
            },
            {
                "candidate_rank": "1",
                "candidate_id": candidate_id,
                "step_id": step_id,
                "hypothesis_id": "RP2STOICH:" + "4" * 64,
                "rule_id": "RR-TEST",
                "source_mnxr_id": "MNXR1",
                "source_orientation": "left_to_right",
            },
        )
        base_terms = [
            {
                "candidate_rank": "1",
                "candidate_id": candidate_id,
                "step_id": step_id,
                "side": "left",
                "coefficient": "1",
                "role": "predicted_core",
                "compound_id": "RP2CPD:a",
                "source_mnxm_id": "",
                "name": "methane",
                "formula": "CH4",
                "charge": "0",
                "inchi": "InChI=1S/CH4/h1H4",
                "inchikey": "VNWKTOKETHGBQD-UHFFFAOYSA-N",
                "smiles": "C",
                "xrefs": "",
            },
            {
                "candidate_rank": "1",
                "candidate_id": candidate_id,
                "step_id": step_id,
                "side": "right",
                "coefficient": "1",
                "role": "predicted_core",
                "compound_id": "RP2CPD:b",
                "source_mnxm_id": "",
                "name": "water",
                "formula": "H2O",
                "charge": "0",
                "inchi": "InChI=1S/H2O/h1H2",
                "inchikey": "XLYOFNOQVPJJNP-UHFFFAOYSA-N",
                "smiles": "O",
                "xrefs": "",
            },
        ]
        terms = []
        for hypothesis in hypotheses:
            for row in base_terms:
                terms.append({**row, "hypothesis_id": hypothesis["hypothesis_id"]})
        terms[-2] = {
            **terms[-2],
            "formula": "C2H6",
            "inchi": "",
            "inchikey": "",
            "smiles": "",
        }
        inputs = _ValidatedP9Inputs(
            candidate_inputs=SimpleNamespace(),
            candidate_rank=1,
            candidate_id=candidate_id,
            route={},
            steps=({
                "candidate_id": candidate_id,
                "step_index": "1",
                "step_id": step_id,
                "step_source": "retropath",
                "status": "predicted",
                "reaction_smiles": "C>>O",
                "source_reaction_ids": "MNXR1",
                "source_ec_numbers": "1.1.1.1",
                "source_uniprot_ids": "P12345",
            },),
            p8_dir=self.root,
            p8_manifest_path=self.root / "validation_manifest.json",
            p8_manifest={},
            passing_rows=(
                {
                    "combination_id": "RP2COMB:one",
                    "stoichiometry_hypothesis_ids": hypotheses[0]["hypothesis_id"],
                },
                {
                    "combination_id": "RP2COMB:two",
                    "stoichiometry_hypothesis_ids": hypotheses[1]["hypothesis_id"],
                },
            ),
            hypotheses=hypotheses,
            terms=tuple(terms),
        )
        requirements = _build_requirements(
            inputs,
            _FakeIndex(),
            {"RR-TEST": {"Rule": "[C:1]>>[O:1]", "EC number": "1.1.1.1"}},
        )
        self.assertEqual(2, len(requirements))
        self.assertTrue(requirements[0].formal_mapping_exact)
        self.assertFalse(requirements[1].formal_mapping_exact)
        self.assertNotEqual(
            requirements[0].reaction_signature_sha256,
            requirements[1].reaction_signature_sha256,
        )
        self.assertEqual(("10000",), requirements[1].source_rhea_ids)
        self.assertEqual(tuple(), requirements[1].exact_rhea_ids)
        self.assertEqual("", requirements[1].full_reaction_smiles)
        self.assertEqual("C>>O", requirements[1].core_reaction_smiles)

    def test_structural_similarity_one_remains_manual_review(self) -> None:
        requirement = RetropathEnzymeRequirement(
            candidate_rank=1,
            candidate_id="RP2ROUTE:" + "1" * 64,
            combination_id="RP2COMB:test",
            step_index=1,
            step_id="RP2STEP:" + "2" * 64,
            step_source="retropath",
            step_status="predicted",
            hypothesis_id="RP2STOICH:" + "3" * 64,
            reaction_signature_sha256="4" * 64,
            full_reaction_smiles="CCO>>CC=O",
            core_reaction_smiles="CCO>>CC=O",
            rule_id="RR-TEST",
            rule_smarts="[C:1][O:2]>>[C:1]=[O:2]",
            source_mnxr_id="MNXR1",
            source_reaction_ids=("MNXR1",),
            source_ec_numbers=tuple(),
            source_uniprot_ids=tuple(),
            source_rhea_ids=tuple(),
            exact_kegg_reaction_ids=tuple(),
            exact_rhea_ids=tuple(),
            formal_mapping_exact=False,
        )
        query_result = {
            "query_id": "selenzyme_fixture",
            "query_type": "selenzyme_by_full_reaction_smiles",
            "query_kind": "full_reaction_smiles",
            "query_sha256": "5" * 64,
            "status": "ok",
            "rows": [],
        }
        client = SimpleNamespace(
            query_reaction_smarts=lambda *args, **kwargs: query_result
        )
        with patch(
            "src.main_protein_selection.retropath_enzyme_selection."
            "retrieve_selenzyme_candidates",
            return_value=([_protein(1.0)], [], [], {}),
        ):
            outcome = _search_requirement(
                requirement,
                chassis_key="ecoli_mg1655",
                top_n=5,
                max_results=100,
                allow_transmembrane=False,
                session=SimpleNamespace(),
                rhea_client=SimpleNamespace(),
                selenzyme_client=client,
                selenzyme_circuit_error="",
                entry_cache={},
            )
        self.assertEqual(1, len(outcome.candidates))
        candidate = outcome.candidates[0]
        self.assertEqual("full_reaction_similarity", candidate["evidence_tier"])
        self.assertEqual("manual_review", candidate["fit_status"])
        self.assertEqual("true", candidate["manual_review_required"])

    def _write_p5_p8_fixture(self) -> tuple[SimpleNamespace, Path]:
        target = "C00001"
        candidate_id = "RP2ROUTE:" + "1" * 64
        step_id = "RP2STEP:" + "2" * 64
        hypothesis_id = "RP2STOICH:" + "3" * 64
        gap_dir = self.root / "gap"
        retropath_dir = gap_dir / "depth0" / "retropath"
        route_path = retropath_dir / CANDIDATE_ROUTES_FILE_NAME
        step_path = retropath_dir / CANDIDATE_STEPS_FILE_NAME
        rejected_path = retropath_dir / REJECTED_ROUTES_FILE_NAME
        mapping_path = retropath_dir / "input" / "compound_mapping.csv"
        _write_csv(route_path, CANDIDATE_ROUTE_COLUMNS, [{
            "candidate_rank": 1,
            "candidate_id": candidate_id,
            "target_compound_id": target,
            "total_steps": 1,
        }])
        _write_csv(step_path, CANDIDATE_STEP_COLUMNS, [{
            "candidate_id": candidate_id,
            "step_index": 1,
            "step_id": step_id,
            "step_source": "retropath",
            "status": "predicted",
            "reaction_smiles": "C>>O",
            "rule_ids": "RR-TEST",
            "source_reaction_ids": "MNXR1",
            "source_ec_numbers": "1.1.1.1",
        }])
        _write_csv(rejected_path, REJECTED_ROUTE_COLUMNS, [])
        mapping_path.parent.mkdir(parents=True, exist_ok=True)
        mapping_path.write_text("compound_id\n", encoding="utf-8")
        pipeline_path = retropath_dir / PIPELINE_RESULT_FILE_NAME
        pipeline = {
            "schema_version": RETROPATH_PIPELINE_SCHEMA,
            "ok": True,
            "status": "retropath_candidates_found",
            "target_compound": target,
            "expansion_depth": 0,
            "candidate_count": 1,
            "artifacts": {
                "candidate_routes": {
                    "path": str(route_path.resolve()),
                    "sha256": _sha256(route_path),
                },
                "candidate_steps": {
                    "path": str(step_path.resolve()),
                    "sha256": _sha256(step_path),
                },
                "rejected_routes": {
                    "path": str(rejected_path.resolve()),
                    "sha256": _sha256(rejected_path),
                },
                "compound_mapping": {
                    "path": str(mapping_path.resolve()),
                    "sha256": _sha256(mapping_path),
                },
            },
        }
        pipeline_path.write_text(json.dumps(pipeline), encoding="utf-8")

        p8_dir = retropath_dir / "gem_validation"
        hypothesis_path = p8_dir / STOICHIOMETRY_HYPOTHESES_FILE_NAME
        terms_path = p8_dir / STOICHIOMETRY_TERMS_FILE_NAME
        summary_path = p8_dir / VALIDATION_SUMMARY_FILE_NAME
        _write_csv(hypothesis_path, HYPOTHESIS_COLUMNS, [{
            "candidate_rank": 1,
            "candidate_id": candidate_id,
            "step_id": step_id,
            "hypothesis_id": hypothesis_id,
            "rule_id": "RR-TEST",
            "source_mnxr_id": "MNXR1",
            "source_orientation": "left_to_right",
        }])
        _write_csv(terms_path, TERM_COLUMNS, [{
            "candidate_rank": 1,
            "candidate_id": candidate_id,
            "step_id": step_id,
            "hypothesis_id": hypothesis_id,
            "side": "left",
            "coefficient": 1,
            "charge": 0,
            "inchi": "InChI=1S/CH4/h1H4",
            "inchikey": "VNWKTOKETHGBQD-UHFFFAOYSA-N",
            "smiles": "C",
        }, {
            "candidate_rank": 1,
            "candidate_id": candidate_id,
            "step_id": step_id,
            "hypothesis_id": hypothesis_id,
            "side": "right",
            "coefficient": 1,
            "charge": 0,
            "inchi": "InChI=1S/H2O/h1H2",
            "inchikey": "XLYOFNOQVPJJNP-UHFFFAOYSA-N",
            "smiles": "O",
        }])
        _write_csv(summary_path, SUMMARY_COLUMNS, [{
            "candidate_rank": 1,
            "candidate_id": candidate_id,
            "combination_id": "RP2COMB:test",
            "candidate_status": "PASS_STRICT_HYPOTHESIS_EXISTS",
            "validation_status": "PASS_STRICT_ROUTE_FLUX",
            "stoichiometry_hypothesis_ids": hypothesis_id,
        }])
        manifest_path = p8_dir / VALIDATION_MANIFEST_FILE_NAME
        manifest = {
            "schema_version": RETROPATH_GEM_VALIDATION_SCHEMA,
            "target_compound": target,
            "expansion_depth": 0,
            "candidate_statuses": {"1": "PASS_STRICT_HYPOTHESIS_EXISTS"},
            "inputs": {
                "pipeline_result": {
                    "path": str(pipeline_path.resolve()),
                    "sha256": _sha256(pipeline_path),
                },
                "candidate_routes_sha256": _sha256(route_path),
                "candidate_steps_sha256": _sha256(step_path),
            },
            "artifacts": {
                name: {"path": str(path.resolve()), "sha256": _sha256(path)}
                for name, path in {
                    STOICHIOMETRY_HYPOTHESES_FILE_NAME: hypothesis_path,
                    STOICHIOMETRY_TERMS_FILE_NAME: terms_path,
                    VALIDATION_SUMMARY_FILE_NAME: summary_path,
                }.items()
            },
        }
        manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
        config = SimpleNamespace(
            target_name=target,
            gap_output_path=gap_dir,
            depth=0,
            retropath_candidate=1,
        )
        return config, summary_path

    def test_p8_gate_rejects_tampered_artifact(self) -> None:
        config, summary_path = self._write_p5_p8_fixture()
        inputs = _load_p9_inputs(config)
        self.assertEqual(1, inputs.candidate_rank)
        summary_path.write_text("tampered", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "checksum"):
            _load_p9_inputs(config)

    def test_orchestrator_writes_isolated_hashed_review_artifacts(self) -> None:
        rules_path = self.root / "rules.csv"
        rules_path.write_text(
            '"Rule ID","Rule","EC number"\n'
            '"RR-TEST","C>>O","NOEC"\n',
            encoding="utf-8",
        )
        pipeline_path = self.root / "pipeline.json"
        pipeline_path.write_text("{}", encoding="utf-8")
        p8_path = self.root / "validation_manifest.json"
        p8_path.write_text("{}", encoding="utf-8")
        retropath_dir = self.root / "retropath"
        p8_manifest = {
            "inputs": {"rr02_sha256": _sha256(rules_path)},
        }
        candidate_inputs = SimpleNamespace(
            retropath_dir=retropath_dir,
            target_compound="C00001",
            depth=0,
            pipeline_path=pipeline_path,
            pipeline={
                "artifacts": {
                    "candidate_routes": {"sha256": "1" * 64},
                    "candidate_steps": {"sha256": "2" * 64},
                }
            },
        )
        fake_inputs = _ValidatedP9Inputs(
            candidate_inputs=candidate_inputs,
            candidate_rank=1,
            candidate_id="RP2ROUTE:" + "1" * 64,
            route={},
            steps=tuple(),
            p8_dir=self.root,
            p8_manifest_path=p8_path,
            p8_manifest=p8_manifest,
            passing_rows=tuple(),
            hypotheses=tuple(),
            terms=tuple(),
        )
        requirement = RetropathEnzymeRequirement(
            candidate_rank=1,
            candidate_id=fake_inputs.candidate_id,
            combination_id="RP2COMB:test",
            step_index=1,
            step_id="RP2STEP:" + "2" * 64,
            step_source="retropath",
            step_status="predicted",
            hypothesis_id="RP2STOICH:" + "3" * 64,
            reaction_signature_sha256="4" * 64,
            full_reaction_smiles="C>>O",
            core_reaction_smiles="C>>O",
            rule_id="RR-TEST",
            rule_smarts="C>>O",
            source_mnxr_id="MNXR1",
            source_reaction_ids=("MNXR1",),
            source_ec_numbers=tuple(),
            source_uniprot_ids=tuple(),
            source_rhea_ids=tuple(),
            exact_kegg_reaction_ids=tuple(),
            exact_rhea_ids=tuple(),
            formal_mapping_exact=False,
        )
        record = {
            "protein_candidate_rank": 1,
            "accession": "P12345",
            "reviewed": "true",
            "evidence_tier": "full_reaction_similarity",
            "evidence_tiers": "full_reaction_similarity",
            "fit_status": "manual_review",
            "manual_review_required": "true",
            "selection_status": "selected",
            "reaction_similarity": 1.0,
            "protein_score": 88.0,
            "direction_verdict": "unknown",
        }

        class FakeIndex:
            manifest = {
                "index_path": str(self.root / "index.sqlite3"),
                "index_sha256": "5" * 64,
            }

            def __enter__(inner_self):
                return inner_self

            def __exit__(inner_self, *args):
                return None

        config = SimpleNamespace(
            retropath_rules_path=rules_path,
            data_dir=self.root,
            cache_dir=self.root / "cache",
            chassis_key="ecoli_mg1655",
        )
        with (
            patch(
                "src.main_protein_selection.retropath_enzyme_selection._load_p9_inputs",
                return_value=fake_inputs,
            ),
            patch(
                "src.main_protein_selection.retropath_enzyme_selection.MnxrefIndex",
                return_value=FakeIndex(),
            ),
            patch(
                "src.main_protein_selection.retropath_enzyme_selection._build_requirements",
                return_value=[requirement],
            ),
            patch(
                "src.main_protein_selection.retropath_enzyme_selection._search_requirement",
                return_value=_SearchOutcome([record], [record], [], False),
            ),
        ):
            result = select_retropath_enzymes(
                config,
                selenzyme_client=SimpleNamespace(),
            )
        self.assertTrue(result["ok"])
        selection_path = Path(result["selection_manifest"])
        selection = json.loads(selection_path.read_text(encoding="utf-8"))
        self.assertEqual(RETROPATH_ENZYME_SELECTION_SCHEMA, selection["schema_version"])
        self.assertEqual("ready_for_review", selection["status"])
        self.assertTrue(selection["review_required"])
        self.assertFalse(selection["formal_promotion_allowed"])
        candidate_path = selection_path.parent / STEP_ENZYME_CANDIDATES_FILE_NAME
        with candidate_path.open("r", encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(ENZYME_CANDIDATE_COLUMNS, tuple(rows[0]))
        self.assertEqual("manual_review", rows[0]["fit_status"])
        self.assertEqual("1.0", rows[0]["reaction_similarity"])
        self.assertEqual(
            _sha256(candidate_path),
            selection["artifacts"][STEP_ENZYME_CANDIDATES_FILE_NAME]["sha256"],
        )

    def test_cli_keeps_default_and_adds_retropath_candidate_mode(self) -> None:
        parser = build_parser()
        regular = parser.parse_args(["main-enzyme", "-i", "demo.json"])
        self.assertIsNone(regular.retropath_candidate)
        self.assertEqual(0, regular.depth)
        retropath = parser.parse_args([
            "main-enzyme",
            "-i",
            "demo.json",
            "--retropath-candidate",
            "2",
            "--depth",
            "3",
        ])
        self.assertEqual(2, retropath.retropath_candidate)
        self.assertEqual(3, retropath.depth)

    def test_main_enzyme_dispatches_only_explicit_retropath_mode(self) -> None:
        retropath_config = SimpleNamespace(retropath_candidate=1)
        with patch(
            "src.main_protein_selection.retropath_enzyme_selection."
            "run_retropath_enzyme_selection",
            return_value={"status": "ready_for_review"},
        ) as retropath_run:
            result = run_main_protein_selection(retropath_config)
        self.assertEqual("ready_for_review", result["status"])
        retropath_run.assert_called_once_with(retropath_config)

        regular_config = SimpleNamespace(
            retropath_candidate=None,
            manifest_output_path=self.root / "manifest.json",
            project_output_path=self.root,
            cache_dir=self.root / "cache",
        )
        with patch(
            "src.main_protein_selection.select_main_enzymes.select_main_enzymes",
            return_value={"status": "complete"},
        ) as regular_run:
            result = run_main_protein_selection(regular_config)
        self.assertEqual("complete", result["status"])
        regular_run.assert_called_once()


if __name__ == "__main__":
    unittest.main()
