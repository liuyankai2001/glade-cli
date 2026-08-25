from __future__ import annotations

import csv
import json
import tempfile
import unittest
from pathlib import Path

from src.pathway_analyze.retropath_mnxref import (
    MNXREF_FILES,
    MnxrefIndex,
    MnxrefIndexError,
    build_mnxref_subset,
    validate_mnxref_index,
)


class MnxrefSubsetTests(unittest.TestCase):
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
        self.sources = {name: self.root / name for name in MNXREF_FILES}
        self.sources["reac_prop.tsv"].write_text(
            "# fixture\n"
            "MNXR1\t1 MNXM2@MNXD1 + 1 MNXM10@MNXD1 = 1 MNXM11@MNXD1\t"
            "water + propanal = propanediol\ttrue\tfalse\tfixture:r1\n"
            "MNXR2\t1 MNXM20@MNXD1 + 1 MNXM10@MNXD1 = 1 MNXM11@MNXD1\t"
            "water-alt + propanal = propanediol\ttrue\tfalse\tfixture:r2\n"
            "MNXR999\t1 MNXM2@MNXD1 = 1 MNXM2@MNXD1\tignored\ttrue\tfalse\tfixture\n",
            encoding="utf-8",
        )
        chemicals = (
            (
                "MNXM2", "water", "H2O", "0", "18", "InChI=1S/H2O/h1H2",
                "O", "chebi:15377", "XLYOFNOQVPJJNP-UHFFFAOYSA-N",
            ),
            (
                "MNXM20", "water-alt", "H2O", "0", "18", "InChI=1S/H2O/h1H2",
                "O", "fixture:water", "XLYOFNOQVPJJNP-UHFFFAOYSA-N",
            ),
            (
                "MNXM10", "propanal", "C3H6O", "0", "58",
                "InChI=1S/C3H6O/c1-2-3-4/h3H,2H2,1H3", "CCC=O",
                "fixture:propanal", "NBBJYMSMWIIQGU-UHFFFAOYSA-N",
            ),
            (
                "MNXM11", "propanediol", "C3H8O2", "0", "76",
                "InChI=1S/C3H8O2/c1-3(5)2-4/h3-5H,2H2,1H3", "CC(O)CO",
                "fixture:diol", "PQRJMXQQBMCJPP-UHFFFAOYSA-N",
            ),
        )
        self.sources["chem_prop.tsv"].write_text(
            "# fixture\n" + "".join("\t".join(row) + "\n" for row in chemicals),
            encoding="utf-8",
        )
        self.sources["chem_xref.tsv"].write_text(
            "# fixture\n"
            "kegg:C00001\tMNXM2\tidentity\twater\n"
            "fixture:propanal\tMNXM10\tidentity\tpropanal\n"
            "fixture:diol\tMNXM11\tidentity\tpropanediol\n",
            encoding="utf-8",
        )
        self.sources["reac_xref.tsv"].write_text(
            "fixture:R1\tMNXR1\nfixture:R2\tMNXR2\n",
            encoding="utf-8",
        )
        self.output = self.root / "index"

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def test_build_validate_and_query_subset(self) -> None:
        manifest = build_mnxref_subset(
            rules_path=self.rules,
            source_paths=self.sources,
            output_dir=self.output,
        )

        self.assertEqual(2, manifest["counts"]["indexed_reactions"])
        self.assertEqual(4, manifest["counts"]["indexed_chemicals"])
        validate_mnxref_index(self.output, self.rules)
        with MnxrefIndex(self.output, self.rules) as index:
            templates = index.templates_for_rules(["RR-TEST"])
            self.assertEqual(["MNXR1", "MNXR2"], [item.mnxr_id for item in templates])
            self.assertTrue(all(item.balanced for item in templates))
            chemicals = index.chemicals(["MNXM2", "MNXM10"])
            self.assertEqual("H2O", chemicals["MNXM2"].formula)
            self.assertIn("kegg:C00001", chemicals["MNXM2"].xrefs)

    def test_tampered_index_and_rules_are_rejected(self) -> None:
        build_mnxref_subset(
            rules_path=self.rules,
            source_paths=self.sources,
            output_dir=self.output,
        )
        index_path = Path(
            json.loads(
                (self.output / "mnxref_rr02_subset_manifest.json").read_text(
                    encoding="utf-8"
                )
            )["index_path"]
        )
        with index_path.open("ab") as handle:
            handle.write(b"tampered")
        with self.assertRaisesRegex(MnxrefIndexError, "checksum"):
            validate_mnxref_index(self.output, self.rules)

        build_mnxref_subset(
            rules_path=self.rules,
            source_paths=self.sources,
            output_dir=self.output,
        )
        self.rules.write_text(
            self.rules.read_text(encoding="utf-8") + "\n",
            encoding="utf-8",
        )
        with self.assertRaisesRegex(MnxrefIndexError, "different RR02"):
            validate_mnxref_index(self.output, self.rules)


if __name__ == "__main__":
    unittest.main()
