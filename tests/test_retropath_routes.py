from __future__ import annotations

import unittest

from rdkit import Chem
from rdkit.Chem import inchi as rd_inchi
from rdkit.Chem import rdMolDescriptors

from src.pathway_analyze.retropath_models import PredictedCompound, PredictedReaction
from src.pathway_analyze.retropath_parser import (
    ParsedCompoundNode,
    ParsedRetroPathNetwork,
    ParsedTransformation,
    SinkMatch,
)
from src.pathway_analyze.retropath_routes import enumerate_sink_routes


def make_compound(
    smiles: str,
    *,
    kegg_id: str | None = None,
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
        kegg_ids=tuple() if kegg_id is None else (kegg_id,),
        minimum_depth=minimum_depth,
        structure_provenance=("test",),
    )


class RetroPathRouteTests(unittest.TestCase):
    def setUp(self) -> None:
        self.target = make_compound("CCO", kegg_id="C90000")
        self.intermediate = make_compound("CC=O")
        self.sink_one = make_compound("CO", kegg_id="C00002", minimum_depth=0)
        self.sink_two = make_compound("C", kegg_id="C00003", minimum_depth=2)

    @staticmethod
    def transformation(
        transformation_id: str,
        substrate: PredictedCompound,
        products: tuple[PredictedCompound, ...],
        *,
        iteration: int,
        rule_id: str,
        diameter: int = 8,
        score: float = 0.5,
    ) -> ParsedTransformation:
        reaction_smiles = f"{substrate.isomeric_smiles}>>" + ".".join(
            item.isomeric_smiles or "" for item in products
        )
        reaction = PredictedReaction.create(
            rule_id=rule_id,
            reaction_smiles=reaction_smiles,
            substrate_compounds=(substrate.compound_id,),
            product_compounds=tuple(item.compound_id for item in products),
            orientation="retrosynthetic",
            rule_specificity=diameter,
            rule_specificity_semantics="diameter",
            rule_score_raw=score,
            score_semantics="lower_is_better",
        )
        return ParsedTransformation(
            transformation_id=transformation_id,
            substrate_compound_id=substrate.compound_id,
            product_compound_ids=tuple(item.compound_id for item in products),
            reaction_smiles=reaction_smiles,
            iteration=iteration,
            diameter=diameter,
            score_raw=score,
            score_semantics="lower_is_better",
            rule_ids=(rule_id,),
            reported_ec_numbers=tuple(),
            reaction_variants=(reaction,),
            auxiliary_fragments=tuple(),
        )

    @staticmethod
    def sink_match(compound: PredictedCompound) -> SinkMatch:
        assert compound.inchikey is not None
        assert compound.minimum_depth is not None
        return SinkMatch(
            compound_id=compound.compound_id,
            inchikey=compound.inchikey,
            representative_kegg_id=compound.compound_id,
            kegg_ids=compound.kegg_ids,
            minimum_depth=compound.minimum_depth,
            wrapper_in_sink=True,
            wrapper_sink_names=(compound.compound_id,),
        )

    def network(
        self,
        transformations: tuple[ParsedTransformation, ...],
        *,
        status: str = "succeeded",
        sinks: tuple[PredictedCompound, ...] | None = None,
        max_steps: int = 3,
    ) -> ParsedRetroPathNetwork:
        sink_compounds = sinks or (self.sink_one, self.sink_two)
        compounds = {
            compound.compound_id: compound
            for compound in (
                self.target,
                self.intermediate,
                self.sink_one,
                self.sink_two,
                *sink_compounds,
            )
        }
        return ParsedRetroPathNetwork(
            status=status,
            target_compound_id=self.target.compound_id,
            max_steps=max_steps,
            compounds=tuple(
                ParsedCompoundNode(
                    compound=item,
                    is_target=item.compound_id == self.target.compound_id,
                    wrapper_in_sink=item in sink_compounds,
                    wrapper_sink_names=(item.compound_id,)
                    if item in sink_compounds
                    else tuple(),
                )
                for item in sorted(
                    compounds.values(), key=lambda value: value.compound_id
                )
            ),
            transformations=transformations,
            sink_matches=tuple(self.sink_match(item) for item in sink_compounds),
            rule_evidence=tuple(),
            rejections=tuple(),
            warnings=tuple(),
        )

    def complete_branch_network(self) -> ParsedRetroPathNetwork:
        return self.network(
            (
                self.transformation(
                    "TRS_ROOT",
                    self.target,
                    (self.intermediate, self.sink_one),
                    iteration=0,
                    rule_id="RULE-ROOT",
                ),
                self.transformation(
                    "TRS_BRANCH",
                    self.intermediate,
                    (self.sink_two,),
                    iteration=1,
                    rule_id="RULE-BRANCH",
                    diameter=6,
                    score=0.8,
                ),
            )
        )

    def test_complete_and_or_branch_requires_every_leaf_to_hit_sink(self) -> None:
        result = enumerate_sink_routes(self.complete_branch_network())

        self.assertEqual(result.complete_path_count, 1)
        path = result.paths[0]
        self.assertEqual(path.transformation_ids, ("TRS_ROOT", "TRS_BRANCH"))
        self.assertEqual(
            {item.representative_kegg_id for item in path.sink_matches},
            {"C00002", "C00003"},
        )
        self.assertEqual(path.reaction_count, 2)
        self.assertEqual(path.maximum_branch_depth, 2)
        self.assertEqual(path.minimum_rule_specificity, 6)
        self.assertRegex(path.path_id, r"^RP2PATH:[0-9a-f]{64}$")
        self.assertFalse(result.truncated)

    def test_unresolved_branch_rejects_the_whole_candidate(self) -> None:
        root_only = self.network(
            (
                self.transformation(
                    "TRS_ROOT",
                    self.target,
                    (self.intermediate, self.sink_one),
                    iteration=0,
                    rule_id="RULE-ROOT",
                ),
            )
        )

        result = enumerate_sink_routes(root_only)

        self.assertFalse(result.paths)
        self.assertIn(
            "unresolved_non_sink_leaf",
            {item.reason_code for item in result.rejections},
        )

    def test_cycle_is_detected_instead_of_being_linearized(self) -> None:
        cyclic = self.network(
            (
                self.transformation(
                    "TRS_ROOT",
                    self.target,
                    (self.intermediate,),
                    iteration=0,
                    rule_id="RULE-ROOT",
                ),
                self.transformation(
                    "TRS_CYCLE",
                    self.intermediate,
                    (self.target,),
                    iteration=1,
                    rule_id="RULE-CYCLE",
                ),
            )
        )

        result = enumerate_sink_routes(cyclic)

        self.assertFalse(result.paths)
        self.assertIn(
            "cycle_detected", {item.reason_code for item in result.rejections}
        )

    def test_iteration_and_depth_are_hard_path_constraints(self) -> None:
        wrong_iteration = self.network(
            (
                self.transformation(
                    "TRS_ROOT",
                    self.target,
                    (self.intermediate,),
                    iteration=1,
                    rule_id="RULE-ROOT",
                ),
            )
        )
        iteration_result = enumerate_sink_routes(wrong_iteration)
        self.assertIn(
            "reaction_direction_invalid",
            {item.reason_code for item in iteration_result.rejections},
        )

        too_deep = self.network(
            (
                self.transformation(
                    "TRS_ROOT",
                    self.target,
                    (self.intermediate,),
                    iteration=0,
                    rule_id="RULE-ROOT",
                ),
                self.transformation(
                    "TRS_BRANCH",
                    self.intermediate,
                    (self.sink_two,),
                    iteration=1,
                    rule_id="RULE-BRANCH",
                ),
            ),
            max_steps=1,
        )
        depth_result = enumerate_sink_routes(too_deep)
        self.assertIn(
            "depth_exceeded", {item.reason_code for item in depth_result.rejections}
        )

    def test_equivalent_network_transformations_are_deduplicated_by_path_identity(
        self,
    ) -> None:
        first = self.transformation(
            "TRS_Z_HIGH_SPECIFICITY",
            self.target,
            (self.sink_one,),
            iteration=0,
            rule_id="RULE-SAME",
        )
        second = self.transformation(
            "TRS_A_LOW_SPECIFICITY",
            self.target,
            (self.sink_one,),
            iteration=0,
            rule_id="RULE-SAME",
        )

        result = enumerate_sink_routes(self.network((first, second)))

        self.assertEqual(result.complete_path_count, 1)

    def test_route_limit_is_explicit_and_deterministic(self) -> None:
        first = self.transformation(
            "TRS_Z_HIGH_SPECIFICITY",
            self.target,
            (self.sink_one,),
            iteration=0,
            rule_id="RULE-A",
            diameter=8,
        )
        second = self.transformation(
            "TRS_A_LOW_SPECIFICITY",
            self.target,
            (self.sink_two,),
            iteration=0,
            rule_id="RULE-B",
            diameter=6,
        )

        result = enumerate_sink_routes(
            self.network((first, second)),
            max_routes=1,
        )

        self.assertEqual(result.complete_path_count, 1)
        self.assertEqual(
            result.paths[0].transformation_ids,
            ("TRS_Z_HIGH_SPECIFICITY",),
        )
        self.assertTrue(result.truncated)
        self.assertIn(
            "enumeration_limit_reached",
            {item.reason_code for item in result.rejections},
        )

    def test_search_state_limit_preserves_routes_found_before_truncation(self) -> None:
        first = self.transformation(
            "TRS_A",
            self.target,
            (self.sink_one,),
            iteration=0,
            rule_id="RULE-A",
        )
        second = self.transformation(
            "TRS_B",
            self.target,
            (self.sink_two,),
            iteration=0,
            rule_id="RULE-B",
        )

        result = enumerate_sink_routes(
            self.network((first, second)),
            max_search_states=3,
        )

        self.assertEqual(result.complete_path_count, 1)
        self.assertEqual(result.explored_states, 3)
        self.assertTrue(result.truncated)
        self.assertIn(
            "max_search_states=3",
            next(
                item.reason_detail
                for item in result.rejections
                if item.reason_code == "enumeration_limit_reached"
            ),
        )

    def test_non_candidate_terminal_statuses_return_empty_enumeration(self) -> None:
        for status in ("no_solution", "source_in_sink"):
            with self.subTest(status=status):
                result = enumerate_sink_routes(self.network(tuple(), status=status))
                self.assertFalse(result.paths)
                self.assertFalse(result.truncated)


if __name__ == "__main__":
    unittest.main()
