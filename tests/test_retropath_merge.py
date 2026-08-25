from __future__ import annotations

import unittest
from dataclasses import replace

from rdkit import Chem
from rdkit.Chem import inchi as rd_inchi
from rdkit.Chem import rdMolDescriptors

from src.pathway_analyze.expand_chassis_metabolites import (
    ExpansionBundle,
    ForwardExpansionWitness,
)
from src.pathway_analyze.kegg_gap_analyze import ReactionRecord
from src.pathway_analyze.retropath_merge import (
    flip_predicted_reaction,
    merge_retropath_candidates,
)
from src.pathway_analyze.retropath_models import PredictedCompound, PredictedReaction
from src.pathway_analyze.retropath_parser import (
    ParsedRetroPathNetwork,
    ParsedTransformation,
    SinkMatch,
)
from src.pathway_analyze.retropath_routes import (
    RetroPathEnumerationResult,
    RetrosyntheticPath,
    RouteRejection,
)


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


def reaction_record(
    reaction_id: str,
    substrates: tuple[str, ...],
    products: tuple[str, ...],
) -> ReactionRecord:
    return ReactionRecord(
        reaction_id=reaction_id,
        name=reaction_id,
        comment="",
        equation=" + ".join(substrates) + " => " + " + ".join(products),
        annotation_text="",
        left_stoichiometry=tuple((item, 1.0) for item in substrates),
        right_stoichiometry=tuple((item, 1.0) for item in products),
        enzyme_ecs=("1.1.1.1",),
        ko_ids=tuple(),
        pathway_ids=tuple(),
        module_ids=tuple(),
    )


class FakeKeggClient:
    def __init__(self, reactions: dict[str, ReactionRecord]) -> None:
        self.reactions = reactions

    def try_get_reaction(self, reaction_id: str) -> ReactionRecord | None:
        return self.reactions.get(reaction_id)


class RetroPathMergeTests(unittest.TestCase):
    def setUp(self) -> None:
        self.base_id = "C10000"
        self.target = make_compound("CCO", kegg_id="C90000")
        self.intermediate = make_compound("CC=O")
        self.sink_one = make_compound("CO", kegg_id="C10001", minimum_depth=1)
        self.sink_two = make_compound("C", kegg_id="C10002", minimum_depth=2)
        self.reactions = {
            "R10001": reaction_record("R10001", (self.base_id,), ("C10001",)),
            "R10002": reaction_record("R10002", ("C10001",), ("C10002",)),
        }
        self.client = FakeKeggClient(self.reactions)

    @staticmethod
    def witness(
        product: str,
        depth: int,
        reaction_id: str,
        substrates: tuple[str, ...],
        *,
        is_endogenous: bool,
    ) -> ForwardExpansionWitness:
        return ForwardExpansionWitness(
            product_compound=product,
            depth=depth,
            reaction_id=reaction_id,
            direction="left_to_right",
            substrate_compounds=substrates,
            product_compounds=(product,),
            is_endogenous=is_endogenous,
            thermo_direction="favored",
            oxygen_required=False,
            nadph_burden=0.0,
            sam_burden=0.0,
            coa_burden=0.0,
            electron_risk_level="none",
            electron_risk_score=0,
            enzyme_ecs=("1.1.1.1",),
            ko_ids=tuple(),
        )

    def expansion_bundle(
        self,
        *,
        depth: int = 2,
        include_second_witness: bool = True,
        alternate_first_witness: bool = False,
    ) -> ExpansionBundle:
        witnesses: dict[str, tuple[ForwardExpansionWitness, ...]] = {
            "C10001": (
                self.witness(
                    "C10001",
                    1,
                    "R10001",
                    (self.base_id,),
                    is_endogenous=True,
                ),
            )
        }
        if alternate_first_witness:
            self.reactions["R10003"] = reaction_record(
                "R10003",
                (self.base_id,),
                ("C10001",),
            )
            witnesses["C10001"] = (
                *witnesses["C10001"],
                self.witness(
                    "C10001",
                    1,
                    "R10003",
                    (self.base_id,),
                    is_endogenous=True,
                ),
            )
        if include_second_witness:
            witnesses["C10002"] = (
                self.witness(
                    "C10002",
                    2,
                    "R10002",
                    ("C10001",),
                    is_endogenous=False,
                ),
            )
        return ExpansionBundle(
            depth=depth,
            base_compounds=frozenset({self.base_id}),
            reachable_compounds=frozenset({self.base_id, "C10001", "C10002"}),
            depth_by_compound={self.base_id: 0, "C10001": 1, "C10002": 2},
            witnesses_by_product=witnesses,
            expanded_file=None,  # type: ignore[arg-type]
            manifest={},
        )

    @staticmethod
    def sink_match(compound: PredictedCompound, depth: int) -> SinkMatch:
        assert compound.inchikey is not None
        return SinkMatch(
            compound_id=compound.compound_id,
            inchikey=compound.inchikey,
            representative_kegg_id=compound.compound_id,
            kegg_ids=compound.kegg_ids,
            minimum_depth=depth,
            wrapper_in_sink=True,
            wrapper_sink_names=(compound.compound_id,),
        )

    @staticmethod
    def transformation(
        transformation_id: str,
        substrate: PredictedCompound,
        products: tuple[PredictedCompound, ...],
        *,
        iteration: int,
        rules: tuple[str, ...],
        diameter: int,
        score: float,
    ) -> ParsedTransformation:
        reaction_smiles = f"{substrate.isomeric_smiles}>>" + ".".join(
            item.isomeric_smiles or "" for item in products
        )
        variants = tuple(
            PredictedReaction.create(
                rule_id=rule_id,
                reaction_smiles=reaction_smiles,
                substrate_compounds=(substrate.compound_id,),
                product_compounds=tuple(item.compound_id for item in products),
                orientation="retrosynthetic",
                source_reaction_ids=(f"MNXR{index + 1}",),
                source_ec_numbers=(f"1.1.1.{index + 1}",),
                rule_specificity=diameter,
                rule_specificity_semantics="diameter",
                rule_score_raw=score + index * 0.1,
                score_semantics="lower_is_better",
            )
            for index, rule_id in enumerate(rules)
        )
        return ParsedTransformation(
            transformation_id=transformation_id,
            substrate_compound_id=substrate.compound_id,
            product_compound_ids=tuple(sorted(item.compound_id for item in products)),
            reaction_smiles=reaction_smiles,
            iteration=iteration,
            diameter=diameter,
            score_raw=score,
            score_semantics="lower_is_better",
            rule_ids=rules,
            reported_ec_numbers=tuple(),
            reaction_variants=variants,
            auxiliary_fragments=tuple(),
        )

    def enumeration(
        self,
        *,
        sink_depths: tuple[int, int] = (1, 2),
        sink_order_reversed: bool = False,
        status: str = "succeeded",
        p4_rejections: tuple[RouteRejection, ...] = tuple(),
    ) -> RetroPathEnumerationResult:
        root = self.transformation(
            "TRS_ROOT",
            self.target,
            (self.intermediate, self.sink_one),
            iteration=0,
            rules=("RULE-A", "RULE-B"),
            diameter=8,
            score=0.4,
        )
        branch = self.transformation(
            "TRS_BRANCH",
            self.intermediate,
            (self.sink_two,),
            iteration=1,
            rules=("RULE-C",),
            diameter=6,
            score=0.8,
        )
        sink_matches = (
            self.sink_match(self.sink_one, sink_depths[0]),
            self.sink_match(self.sink_two, sink_depths[1]),
        )
        if sink_order_reversed:
            sink_matches = tuple(reversed(sink_matches))
        network = ParsedRetroPathNetwork(
            status=status,
            target_compound_id=self.target.compound_id,
            max_steps=3,
            compounds=tuple(),
            transformations=(root, branch) if status == "succeeded" else tuple(),
            sink_matches=sink_matches if status == "succeeded" else tuple(),
            rule_evidence=tuple(),
            rejections=tuple(),
            warnings=tuple(),
        )
        paths = tuple()
        if status == "succeeded":
            paths = (
                RetrosyntheticPath(
                    path_id=f"RP2PATH:{'a' * 64}",
                    target_compound_id=self.target.compound_id,
                    transformation_ids=("TRS_ROOT", "TRS_BRANCH"),
                    reaction_ids_by_transformation=(
                        (
                            "TRS_ROOT",
                            tuple(item.reaction_id for item in root.reaction_variants),
                        ),
                        (
                            "TRS_BRANCH",
                            tuple(
                                item.reaction_id for item in branch.reaction_variants
                            ),
                        ),
                    ),
                    sink_matches=sink_matches,
                    reaction_count=2,
                    maximum_branch_depth=2,
                    minimum_rule_specificity=6,
                    worst_rule_score=0.8,
                    score_semantics="lower_is_better",
                    contains_auxiliary_fragments=False,
                ),
            )
        return RetroPathEnumerationResult(
            network=network,
            paths=paths,
            rejections=p4_rejections,
            explored_states=10,
            max_routes=1000,
            max_search_states=100_000,
            truncated=False,
        )

    def depth_zero_bundle(self) -> ExpansionBundle:
        return ExpansionBundle(
            depth=0,
            base_compounds=frozenset({"C10001", "C10002"}),
            reachable_compounds=frozenset({"C10001", "C10002"}),
            depth_by_compound={"C10001": 0, "C10002": 0},
            witnesses_by_product={},
            expanded_file=None,  # type: ignore[arg-type]
            manifest={},
        )

    def test_flip_creates_a_new_biosynthetic_reaction_identity(self) -> None:
        reverse = self.enumeration().network.transformations[0].reaction_variants[0]

        forward = flip_predicted_reaction(reverse)

        self.assertEqual(forward.orientation, "biosynthetic")
        self.assertEqual(forward.substrate_compounds, reverse.product_compounds)
        self.assertEqual(forward.product_compounds, reverse.substrate_compounds)
        self.assertEqual(
            forward.reaction_smiles,
            ".".join(
                (
                    self.intermediate.isomeric_smiles or "",
                    self.sink_one.isomeric_smiles or "",
                )
            )
            + f">>{self.target.isomeric_smiles}",
        )
        self.assertNotEqual(forward.reaction_id, reverse.reaction_id)

    def test_multi_sink_witnesses_and_rp2_steps_form_one_biosynthetic_dag(self) -> None:
        result = merge_retropath_candidates(
            self.enumeration(),
            self.expansion_bundle(),
            self.client,  # type: ignore[arg-type]
        )

        self.assertEqual(result.candidate_count, 1)
        candidate = result.candidates[0]
        self.assertEqual(
            [item.reaction_option_ids[0] for item in candidate.steps[:2]],
            ["R10001", "R10002"],
        )
        self.assertEqual(
            [item.source_transformation_ids for item in candidate.steps[2:]],
            [("TRS_BRANCH",), ("TRS_ROOT",)],
        )
        self.assertEqual(candidate.kegg_prefix_steps, 2)
        self.assertEqual(candidate.retropath_steps, 2)
        self.assertEqual(candidate.total_steps, 4)
        self.assertEqual(
            candidate.candidate_id,
            "RP2ROUTE:02ed7d4cea75cd6d42df95016436f4db40f798a3df5f1fb4ca0458a219fbe065",
        )
        self.assertEqual(
            {item.representative_kegg_id for item in candidate.sink_matches},
            {"C10001", "C10002"},
        )
        root_step = next(
            item
            for item in candidate.steps
            if item.source_transformation_ids == ("TRS_ROOT",)
        )
        self.assertEqual(len(root_step.reaction_option_ids), 2)
        self.assertEqual(result.candidate_count, 1)
        self.assertTrue(
            all(
                item.orientation == "biosynthetic"
                for item in result.biosynthetic_reactions
            )
        )

    def test_depth_zero_has_no_kegg_prefix(self) -> None:
        result = merge_retropath_candidates(
            self.enumeration(sink_depths=(0, 0)),
            self.depth_zero_bundle(),
            self.client,  # type: ignore[arg-type]
        )

        candidate = result.candidates[0]
        self.assertEqual(candidate.kegg_prefix_steps, 0)
        self.assertEqual(candidate.retropath_steps, 2)
        self.assertTrue(
            all(item.step_source == "retropath" for item in candidate.steps)
        )

    def test_missing_witness_and_depth_mismatch_are_rejected(self) -> None:
        missing = merge_retropath_candidates(
            self.enumeration(),
            self.expansion_bundle(include_second_witness=False),
            self.client,  # type: ignore[arg-type]
        )
        self.assertFalse(missing.candidates)
        self.assertIn(
            "expansion_witness_missing",
            {item.reason_code for item in missing.rejections},
        )

        mismatch = merge_retropath_candidates(
            self.enumeration(sink_depths=(1, 1)),
            self.expansion_bundle(),
            self.client,  # type: ignore[arg-type]
        )
        self.assertFalse(mismatch.candidates)
        self.assertIn(
            "sink_depth_mismatch",
            {item.reason_code for item in mismatch.rejections},
        )

    def test_combined_step_and_enzyme_limits_include_rp2_steps(self) -> None:
        result = merge_retropath_candidates(
            self.enumeration(),
            self.expansion_bundle(),
            self.client,  # type: ignore[arg-type]
            max_new_enzymes=2,
        )

        self.assertFalse(result.candidates)
        self.assertIn(
            "candidate_limit_exceeded",
            {item.reason_code for item in result.rejections},
        )

    def test_top_k_is_deterministic_and_reports_pruning(self) -> None:
        result = merge_retropath_candidates(
            self.enumeration(),
            self.expansion_bundle(alternate_first_witness=True),
            self.client,  # type: ignore[arg-type]
            max_candidates=1,
        )

        self.assertEqual(result.candidate_count, 1)
        self.assertTrue(result.truncated)
        self.assertEqual(
            result.candidates[0].kegg_prefix_reaction_ids,
            ("R10001", "R10002"),
        )
        self.assertIn("top_k_pruned", {item.reason_code for item in result.rejections})

    def test_sink_input_order_does_not_change_candidate_identity(self) -> None:
        forward = merge_retropath_candidates(
            self.enumeration(),
            self.expansion_bundle(),
            self.client,  # type: ignore[arg-type]
        )
        reversed_order = merge_retropath_candidates(
            self.enumeration(sink_order_reversed=True),
            self.expansion_bundle(),
            self.client,  # type: ignore[arg-type]
        )

        self.assertEqual(
            forward.candidates[0].candidate_id,
            reversed_order.candidates[0].candidate_id,
        )

    def test_p4_rejections_and_non_candidate_statuses_are_preserved(self) -> None:
        rejection = RouteRejection(
            "cycle_detected",
            "test cycle",
            compound_id=self.intermediate.compound_id,
        )
        result = merge_retropath_candidates(
            self.enumeration(p4_rejections=(rejection,)),
            self.expansion_bundle(),
            self.client,  # type: ignore[arg-type]
        )
        self.assertIn(
            ("p4", "cycle_detected"),
            {(item.source_stage, item.reason_code) for item in result.rejections},
        )

        no_solution = merge_retropath_candidates(
            self.enumeration(status="no_solution"),
            self.expansion_bundle(),
            self.client,  # type: ignore[arg-type]
        )
        self.assertFalse(no_solution.candidates)
        self.assertFalse(no_solution.biosynthetic_reactions)

    def test_unused_sink_boundary_is_not_silently_added_to_candidate(self) -> None:
        enumeration = self.enumeration()
        extra_sink = make_compound("O", kegg_id="C10003", minimum_depth=0)
        extra_match = self.sink_match(extra_sink, 0)
        path = replace(
            enumeration.paths[0],
            sink_matches=(*enumeration.paths[0].sink_matches, extra_match),
        )
        network = replace(
            enumeration.network,
            sink_matches=(*enumeration.network.sink_matches, extra_match),
        )
        enumeration = replace(enumeration, network=network, paths=(path,))
        bundle = self.expansion_bundle()
        bundle = replace(
            bundle,
            base_compounds=frozenset({*bundle.base_compounds, "C10003"}),
            reachable_compounds=frozenset({*bundle.reachable_compounds, "C10003"}),
            depth_by_compound={**bundle.depth_by_compound, "C10003": 0},
        )

        result = merge_retropath_candidates(
            enumeration,
            bundle,
            self.client,  # type: ignore[arg-type]
        )

        self.assertFalse(result.candidates)
        self.assertIn(
            "candidate_merge_invalid",
            {item.reason_code for item in result.rejections},
        )


if __name__ == "__main__":
    unittest.main()
