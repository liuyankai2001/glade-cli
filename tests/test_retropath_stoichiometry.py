from __future__ import annotations

import csv
import tempfile
import unittest
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import inchi as rd_inchi

from src.pathway_analyze.retropath_mnxref import (
    MNXREF_FILES,
    MnxrefIndex,
    build_mnxref_subset,
)
from src.pathway_analyze.retropath_stoichiometry import (
    enumerate_candidate_hypotheses,
    is_balanced,
    parse_formula,
    reconstruct_retropath_step,
)


def _inchikey(smiles: str) -> str:
    molecule = Chem.MolFromSmiles(smiles)
    return rd_inchi.InchiToInchiKey(rd_inchi.MolToInchi(molecule)).upper()


class RetroPathStoichiometryTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.rules = self.root / "rules.csv"
        with self.rules.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=(
                    "Rule ID",
                    "Legacy ID",
                    "Reaction direction",
                    "Rule relative direction",
                    "Rule usage",
                ),
            )
            writer.writeheader()
            writer.writerows(
                [
                    {
                        "Rule ID": "RR-TEST",
                        "Legacy ID": "MNXR1_MNXM10",
                        "Reaction direction": "0",
                        "Rule relative direction": "1",
                        "Rule usage": "both",
                    },
                    {
                        "Rule ID": "RR-TEST",
                        "Legacy ID": "MNXR2_MNXM10",
                        "Reaction direction": "0",
                        "Rule relative direction": "1",
                        "Rule usage": "both",
                    },
                ]
            )
        sources = {name: self.root / name for name in MNXREF_FILES}
        sources["reac_prop.tsv"].write_text(
            "MNXR1\t1 MNXM2@MNXD1 + 1 MNXM10@MNXD1 = 1 MNXM11@MNXD1\t"
            "water + propanal = propanediol\ttrue\tfalse\tfixture:r1\n"
            "MNXR2\t1 MNXM20@MNXD1 + 1 MNXM10@MNXD1 = 1 MNXM11@MNXD1\t"
            "water-alt + propanal = propanediol\ttrue\tfalse\tfixture:r2\n",
            encoding="utf-8",
        )
        sources["chem_prop.tsv"].write_text(
            "MNXM2\twater\tH2O\t0\t18\tInChI=1S/H2O/h1H2\tO\t"
            "chebi:15377\tXLYOFNOQVPJJNP-UHFFFAOYSA-N\n"
            "MNXM20\twater-alt\tH2O\t0\t18\tInChI=1S/H2O/h1H2\tO\t"
            "fixture:water\tXLYOFNOQVPJJNP-UHFFFAOYSA-N\n"
            "MNXM10\tpropanal\tC3H6O\t0\t58\t"
            "InChI=1S/C3H6O/c1-2-3-4/h3H,2H2,1H3\tCCC=O\t"
            "fixture:propanal\tNBBJYMSMWIIQGU-UHFFFAOYSA-N\n"
            "MNXM11\tpropanediol\tC3H8O2\t0\t76\t"
            "InChI=1S/C3H8O2/c1-3(5)2-4/h3-5H,2H2,1H3\tCC(O)CO\t"
            "fixture:diol\tPQRJMXQQBMCJPP-UHFFFAOYSA-N\n",
            encoding="utf-8",
        )
        sources["chem_xref.tsv"].write_text(
            "chebi:15377\tMNXM2\tidentity\twater\n",
            encoding="utf-8",
        )
        sources["reac_xref.tsv"].write_text(
            "fixture:R1\tMNXR1\nfixture:R2\tMNXR2\n",
            encoding="utf-8",
        )
        self.index_dir = self.root / "index"
        build_mnxref_subset(
            rules_path=self.rules,
            source_paths=sources,
            output_dir=self.index_dir,
        )
        substrate_id = f"RP2CPD:{_inchikey('CCC=O')}"
        product_id = f"RP2CPD:{_inchikey('CC(O)CO')}"
        self.step = {
            "candidate_id": "RP2ROUTE:test",
            "step_id": "RP2STEP:test",
            "step_source": "retropath",
            "reaction_smiles": "CCC=O>>CC(O)CO",
            "substrate_stoichiometry_json": f'[["{substrate_id}",1.0]]',
            "product_stoichiometry_json": f'[["{product_id}",1.0]]',
            "rule_ids": "RR-TEST",
        }

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def test_all_source_backed_cofactor_hypotheses_are_retained(self) -> None:
        with MnxrefIndex(self.index_dir, self.rules) as index:
            result = reconstruct_retropath_step(self.step, index, {})

        self.assertEqual("complete", result.status)
        self.assertEqual(2, len(result.hypotheses))
        self.assertEqual(
            [("MNXM2",), ("MNXM20",)],
            sorted(item.recovered_compound_ids for item in result.hypotheses),
        )
        self.assertTrue(all(is_balanced(item.terms) for item in result.hypotheses))
        self.assertTrue(
            all(
                item.evidence_grade == "rr02_mnxref_v3_template_balanced"
                for item in result.hypotheses
            )
        )

    def test_unknown_rule_is_not_completed_from_ec_or_balance_guessing(self) -> None:
        self.step["rule_ids"] = "RR-NOT-IN-INDEX"
        with MnxrefIndex(self.index_dir, self.rules) as index:
            result = reconstruct_retropath_step(self.step, index, {})

        self.assertEqual("incomplete", result.status)
        self.assertFalse(result.hypotheses)
        self.assertEqual("source_template_missing", result.rejections[0].reason_code)

    def test_hypothesis_and_combination_caps_are_explicit(self) -> None:
        with MnxrefIndex(self.index_dir, self.rules) as index:
            result = reconstruct_retropath_step(
                self.step,
                index,
                {},
                max_hypotheses=1,
            )
        self.assertTrue(result.truncated)
        combinations, truncated = enumerate_candidate_hypotheses(
            [result, result],
            max_combinations=1,
        )
        self.assertEqual(1, len(combinations))
        self.assertFalse(truncated)

        with MnxrefIndex(self.index_dir, self.rules) as index:
            full = reconstruct_retropath_step(self.step, index, {})
        combinations, truncated = enumerate_candidate_hypotheses(
            [full, full],
            max_combinations=2,
        )
        self.assertEqual(2, len(combinations))
        self.assertTrue(truncated)

    def test_formula_parser_rejects_incomplete_or_symbolic_formulas(self) -> None:
        self.assertEqual({"C": 6.0, "H": 12.0, "O": 6.0}, parse_formula("C6H12O6"))
        for formula in ("", "C2H4R", "C2H4*"):
            with self.subTest(formula=formula):
                with self.assertRaises(ValueError):
                    parse_formula(formula)


if __name__ == "__main__":
    unittest.main()
