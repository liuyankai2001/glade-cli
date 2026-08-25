from __future__ import annotations

import json
import re
import unittest
from dataclasses import replace

from src.pathway_analyze.retropath_models import (
    CandidateRoute,
    PredictedCompound,
    PredictedReaction,
    RETROPATH_MODEL_SCHEMA_VERSION,
    RetroPathRunResult,
    RetroPathRuntimeProvenance,
    retropath_result_from_json,
    retropath_result_to_json,
)


WATER_INCHI = "InChI=1S/H2O/h1H2"
WATER_INCHIKEY = "XLYOFNOQVPJJNP-UHFFFAOYSA-N"
ETHANOL_INCHI = "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
ETHANOL_INCHIKEY = "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"


class RetroPathModelTests(unittest.TestCase):
    def make_known_compound(self) -> PredictedCompound:
        return PredictedCompound.create(
            compound_id="c00001",
            name="Water",
            inchi=WATER_INCHI,
            inchikey=WATER_INCHIKEY.lower(),
            isomeric_smiles="O",
            formula="H2O",
            charge=0,
            kegg_ids=("C00001",),
            minimum_depth=0,
            structure_provenance=("RDKit 2025.03", "KEGG MOL 2026-08-24"),
        )

    def make_predicted_compound(self) -> PredictedCompound:
        return PredictedCompound.create(
            name="Predicted ethanol",
            inchi=ETHANOL_INCHI,
            inchikey=ETHANOL_INCHIKEY,
            isomeric_smiles="CCO",
        )

    def make_reaction(self, **overrides: object) -> PredictedReaction:
        predicted = self.make_predicted_compound()
        values: dict[str, object] = {
            "rule_id": "RR-02-1234",
            "reaction_smiles": "O.CCO>>CC=O",
            "substrate_compounds": ("C00001", predicted.compound_id, "C00001"),
            "product_compounds": (predicted.compound_id,),
            "orientation": "biosynthetic",
            "source_reaction_ids": ("MNXR2", "R00001"),
            "source_ec_numbers": ("1.1.1.1",),
            "source_uniprot_ids": ("P12345",),
            "rule_specificity": 8,
            "rule_specificity_semantics": "diameter",
            "rule_score_raw": 0.91,
            "score_semantics": "higher_is_better",
        }
        values.update(overrides)
        return PredictedReaction.create(**values)  # type: ignore[arg-type]

    def make_route(self, **overrides: object) -> CandidateRoute:
        reaction = self.make_reaction()
        values: dict[str, object] = {
            "target_compound_id": self.make_predicted_compound().compound_id,
            "matched_sink_kegg_id": "C00001",
            "matched_sink_depth": 2,
            "kegg_prefix_reaction_ids": ("R00001", "R00002"),
            "retropath_reaction_ids": (reaction.reaction_id,),
            "minimum_rule_specificity": 8,
        }
        values.update(overrides)
        return CandidateRoute.create(**values)  # type: ignore[arg-type]

    def make_provenance(self) -> RetroPathRuntimeProvenance:
        return RetroPathRuntimeProvenance(
            wrapper_version="3.9.1",
            wrapper_reported_version="3.9.0",
            workflow_version="r20260212",
            knime_version="4.7.0",
            rdkit_plugin_version="4.9.1",
            rules_version="rr02-rp2-hs",
            rules_sha256="a" * 64,
        )

    def test_predicted_compound_uses_inchikey_namespace(self) -> None:
        compound = self.make_predicted_compound()

        self.assertEqual(compound.compound_id, f"RP2CPD:{ETHANOL_INCHIKEY}")
        self.assertFalse(compound.compound_id.startswith("C"))

    def test_predicted_compound_without_inchikey_uses_stable_full_hash(self) -> None:
        first = PredictedCompound.create(inchi=ETHANOL_INCHI)
        second = PredictedCompound.create(inchi=f"  {ETHANOL_INCHI}  ")

        self.assertEqual(first.compound_id, second.compound_id)
        self.assertRegex(first.compound_id, r"^RP2CPD:[0-9a-f]{64}$")

    def test_known_compound_keeps_kegg_id_and_normalizes_evidence(self) -> None:
        compound = PredictedCompound.create(
            compound_id="c00001",
            inchi=WATER_INCHI,
            kegg_ids=("C00002", "c00001", "C00002"),
            structure_provenance=("z", "a", "z"),
        )

        self.assertEqual(compound.compound_id, "C00001")
        self.assertEqual(compound.kegg_ids, ("C00001", "C00002"))
        self.assertEqual(compound.structure_provenance, ("a", "z"))

    def test_staged_optional_compound_fields_are_allowed(self) -> None:
        compound = PredictedCompound.create(inchi=ETHANOL_INCHI)

        self.assertIsNone(compound.inchikey)
        self.assertIsNone(compound.isomeric_smiles)
        self.assertIsNone(compound.formula)
        self.assertIsNone(compound.charge)
        self.assertIsNone(compound.minimum_depth)

    def test_reaction_id_is_stable_for_set_like_evidence_reordering(self) -> None:
        first = self.make_reaction(
            source_reaction_ids=("MNXR2", "R00001", "MNXR2"),
            source_ec_numbers=("2.2.2.2", "1.1.1.1"),
        )
        second = self.make_reaction(
            source_reaction_ids=("R00001", "MNXR2"),
            source_ec_numbers=("1.1.1.1", "2.2.2.2"),
            substrate_compounds=(
                "C00001",
                "C00001",
                self.make_predicted_compound().compound_id,
            ),
        )

        self.assertEqual(first.reaction_id, second.reaction_id)
        self.assertEqual(first.source_reaction_ids, second.source_reaction_ids)
        self.assertRegex(first.reaction_id, r"^RP2:[0-9a-f]{64}$")
        self.assertFalse(first.reaction_id.startswith("R0"))

    def test_reaction_compound_multisets_are_sorted_but_not_deduplicated(self) -> None:
        reaction = self.make_reaction()

        self.assertEqual(reaction.substrate_compounds.count("C00001"), 2)
        self.assertEqual(
            reaction.substrate_compounds,
            tuple(sorted(reaction.substrate_compounds)),
        )

    def test_reaction_direction_and_side_reversal_change_identity(self) -> None:
        baseline = self.make_reaction()
        reversed_orientation = self.make_reaction(orientation="retrosynthetic")
        reversed_sides = self.make_reaction(
            substrate_compounds=baseline.product_compounds,
            product_compounds=baseline.substrate_compounds,
        )

        self.assertNotEqual(baseline.reaction_id, reversed_orientation.reaction_id)
        self.assertNotEqual(baseline.reaction_id, reversed_sides.reaction_id)

    def test_route_id_preserves_step_order_and_boundary_identity(self) -> None:
        reaction = self.make_reaction()
        reaction_two = self.make_reaction(rule_id="RR-02-9999")
        baseline = self.make_route(
            retropath_reaction_ids=(reaction.reaction_id, reaction_two.reaction_id)
        )
        reordered = self.make_route(
            retropath_reaction_ids=(reaction_two.reaction_id, reaction.reaction_id)
        )
        other_depth = self.make_route(matched_sink_depth=3)
        other_target = self.make_route(target_compound_id="C00002")

        self.assertNotEqual(baseline.candidate_id, reordered.candidate_id)
        self.assertNotEqual(baseline.candidate_id, other_depth.candidate_id)
        self.assertNotEqual(baseline.candidate_id, other_target.candidate_id)
        self.assertRegex(baseline.candidate_id, r"^RP2ROUTE:[0-9a-f]{64}$")
        self.assertEqual(baseline.kegg_prefix_steps, 2)
        self.assertEqual(baseline.retropath_steps, 2)
        self.assertEqual(baseline.total_steps, 4)
        self.assertEqual(baseline.route_source, "kegg_retropath")
        self.assertTrue(baseline.contains_predicted_steps)

    def test_depth_zero_route_allows_empty_kegg_prefix(self) -> None:
        route = self.make_route(
            matched_sink_depth=0,
            kegg_prefix_reaction_ids=tuple(),
        )

        self.assertEqual(route.kegg_prefix_steps, 0)
        self.assertEqual(route.total_steps, route.retropath_steps)

    def test_nested_run_result_round_trips_through_dict_and_json(self) -> None:
        compound = self.make_known_compound()
        predicted = self.make_predicted_compound()
        reaction = self.make_reaction()
        route = self.make_route()
        result = RetroPathRunResult(
            job_id="rp2-test-job",
            status="succeeded",
            return_code=0,
            provenance=self.make_provenance(),
            parameters=(("timeout", 600), ("max_steps", 5), ("enabled", True)),
            artifacts=("raw/scope.csv", "raw/results.csv"),
            compounds=(compound, predicted),
            reactions=(reaction,),
            candidate_routes=(route,),
        )

        dictionary_round_trip = RetroPathRunResult.from_dict(result.to_dict())
        json_text = retropath_result_to_json(result)
        json_round_trip = retropath_result_from_json(json_text)

        self.assertEqual(result, dictionary_round_trip)
        self.assertEqual(result, json_round_trip)
        self.assertTrue(json_round_trip.is_terminal)
        self.assertTrue(json_round_trip.has_candidates)
        self.assertEqual(json.loads(json_text)["schema_version"], 1)

    def test_run_result_normalizes_parameters_and_artifacts(self) -> None:
        result = RetroPathRunResult(
            job_id="job",
            status="queued",
            parameters=(("z", 1), ("a", "x")),
            artifacts=("b.csv", "a.csv", "b.csv"),
        )

        self.assertEqual(result.parameters, (("a", "x"), ("z", 1)))
        self.assertEqual(result.artifacts, ("a.csv", "b.csv"))
        self.assertFalse(result.is_terminal)
        self.assertFalse(result.has_candidates)

    def test_unsupported_schema_version_is_rejected(self) -> None:
        payload = {
            "schema_version": RETROPATH_MODEL_SCHEMA_VERSION + 1,
            "job_id": "job",
            "status": "queued",
        }

        with self.assertRaisesRegex(ValueError, "unsupported.*schema_version"):
            RetroPathRunResult.from_dict(payload)

    def test_tampered_reaction_and_route_ids_are_rejected(self) -> None:
        reaction = self.make_reaction()
        route = self.make_route()

        with self.assertRaisesRegex(ValueError, "reaction_id does not match"):
            replace(reaction, reaction_id=f"RP2:{'0' * 64}")
        with self.assertRaisesRegex(ValueError, "candidate_id does not match"):
            replace(route, candidate_id=f"RP2ROUTE:{'0' * 64}")

    def test_serialized_derived_route_fields_are_checked(self) -> None:
        payload = self.make_route().to_dict()
        payload["total_steps"] = 999

        with self.assertRaisesRegex(ValueError, "total_steps does not match"):
            CandidateRoute.from_dict(payload)

    def test_invalid_values_are_rejected(self) -> None:
        reaction = self.make_reaction()
        invalid_calls = (
            lambda: PredictedCompound.create(inchi="InChI=1/C2H6O"),
            lambda: PredictedCompound.create(
                compound_id="C00001",
                inchi=WATER_INCHI,
                kegg_ids=("not-kegg",),
            ),
            lambda: PredictedReaction.create(
                rule_id="rule",
                reaction_smiles="a>>b",
                substrate_compounds=tuple(),
                product_compounds=("C00001",),
                orientation="biosynthetic",
            ),
            lambda: PredictedReaction.create(
                rule_id="rule",
                reaction_smiles="a>>b",
                substrate_compounds=("C00001",),
                product_compounds=("C00002",),
                orientation="sideways",
            ),
            lambda: CandidateRoute.create(
                target_compound_id="C00002",
                matched_sink_kegg_id="C00001",
                matched_sink_depth=-1,
                retropath_reaction_ids=(reaction.reaction_id,),
            ),
            lambda: CandidateRoute.create(
                target_compound_id="C00002",
                matched_sink_kegg_id="C00001",
                matched_sink_depth=0,
                retropath_reaction_ids=(reaction.reaction_id,),
                review_required=False,
            ),
            lambda: RetroPathRunResult(job_id="job", status="unknown"),
            lambda: RetroPathRunResult(job_id="job", status="succeeded"),
        )

        for invalid_call in invalid_calls:
            with self.subTest(call=invalid_call), self.assertRaises(ValueError):
                invalid_call()

    def test_all_generated_hash_namespaces_use_full_sha256(self) -> None:
        compound = PredictedCompound.create(inchi=ETHANOL_INCHI)
        reaction = self.make_reaction()
        route = self.make_route()

        patterns = (
            (compound.compound_id, r"^RP2CPD:[0-9a-f]{64}$"),
            (reaction.reaction_id, r"^RP2:[0-9a-f]{64}$"),
            (route.candidate_id, r"^RP2ROUTE:[0-9a-f]{64}$"),
        )
        for identifier, pattern in patterns:
            with self.subTest(identifier=identifier):
                self.assertIsNotNone(re.fullmatch(pattern, identifier))

    def test_schema_v1_identity_golden_values_do_not_change_silently(self) -> None:
        compound = PredictedCompound.create(inchi=ETHANOL_INCHI)
        reaction = self.make_reaction()
        route = self.make_route()

        self.assertEqual(
            compound.compound_id,
            "RP2CPD:e17e6cda77d21f5abf6c176f48b8e506"
            "b14cebf3d01d96bbd160dd7c84e74a6f",
        )
        self.assertEqual(
            reaction.reaction_id,
            "RP2:cee5f537caad7eed358e53f735263a6c"
            "2a092d7683ea0c5750ad25fcf977233d",
        )
        self.assertEqual(
            route.candidate_id,
            "RP2ROUTE:12e551ecdd56474f10cd3a14e5573268"
            "8163bb7608b1bc2d505c1ca64eaa7c3f",
        )


if __name__ == "__main__":
    unittest.main()
