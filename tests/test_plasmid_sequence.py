import importlib
import importlib.util
import unittest
from copy import deepcopy
from itertools import permutations
from pathlib import Path
from types import SimpleNamespace

from Bio.Restriction import EcoRI, HindIII
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord


def sequence_module(test):
    try:
        spec = importlib.util.find_spec("src.plasmid_design.sequence")
    except ModuleNotFoundError:
        spec = None
    test.assertIsNotNone(spec, "sequence design service must exist")
    return importlib.import_module("src.plasmid_design.sequence")


def expression_record():
    record = SeqRecord(Seq("ATGAAATAACCATGCCCTAAGGATGGGGTAA"), id="expression")
    record.annotations["molecule_type"] = "DNA"
    record.annotations["topology"] = "linear"
    for start, end, gene in [
        (0, 9, "protein_a"),
        (11, 20, "protein_b"),
        (22, 31, "protein_c"),
    ]:
        record.features.append(
            SeqFeature(
                SimpleLocation(start, end, strand=1),
                type="CDS",
                qualifiers={"gene": [gene]},
            )
        )
    return record


class RestrictionTests(unittest.TestCase):
    def test_enzyme_case_is_normalized_and_two_distinct_names_required(self):
        m = sequence_module(self)
        enzymes = m.resolve_enzymes(["hINDiii", " eCoRI "])
        self.assertEqual([e.name for e in enzymes], ["EcoRI", "HindIII"])
        for names in [[], ["EcoRI"], ["ecori", "EcoRI"], ["EcoRI", "HindIII", "BamHI"]]:
            with self.assertRaises(m.DesignError):
                m.resolve_enzymes(names)

    def test_compatible_ends_and_type_iis_are_rejected(self):
        m = sequence_module(self)
        for names in [["BamHI", "BglII"], ["BsaI", "EcoRI"], ["SmaI", "EcoRI"]]:
            with self.assertRaises(m.DesignError):
                m.resolve_enzymes(names)

    def test_circular_closure_site_is_detected_without_duplicate_allowed_sites(self):
        m = sequence_module(self)
        enzymes = m.resolve_enzymes(["EcoRI", "HindIII"])
        issues = m.audit_sequence("TTCCCGAA", enzymes, circular=True)
        self.assertEqual(
            [(i["enzyme"], i["start_bp"], i["end_bp"]) for i in issues],
            [("EcoRI", 6, 11)],
        )
        self.assertEqual(
            m.audit_sequence(
                "GAATTCACGAAGCTT",
                enzymes,
                circular=True,
                allowed={"EcoRI": {1}, "HindIII": {10}},
            ),
            [],
        )


class MolecularDesignTests(unittest.TestCase):
    def setUp(self):
        self.m = sequence_module(self)
        from src.plasmid_design.catalog import ModuleCatalog

        self.catalog = ModuleCatalog(Path(__file__).resolve().parents[1] / "data")

    def test_all_component_orders_preserve_dna_annotations_and_matching_preparations(
        self,
    ):
        insert = expression_record()
        blocks = {
            "resistance": self.catalog.get_resistance("basic_seva_ap").sequence
            + self.catalog.scaffold["resistance_to_replication"],
            "replication": self.catalog.get_replication("basic_seva_p15a").sequence
            + self.catalog.scaffold["replication_to_t1"]
            + self.catalog.scaffold["t1"],
            "expression": "GAATTC" + str(insert.seq) + "AAGCTT",
        }
        for order in permutations(["resistance", "replication", "expression"]):
            with self.subTest(order=order):
                design = self.m.build_design(
                    self.catalog,
                    insert,
                    ["EcoRI", "HindIII"],
                    "basic_seva_ap",
                    "basic_seva_p15a",
                    1,
                    component_order=list(order),
                )
                self.assertTrue(design.preview["valid"], design.preview["issues"])
                self.assertEqual(
                    str(design.final_record.seq), "".join(blocks[k] for k in order)
                )
                self.assertEqual(design.preview["component_order"], list(order))
                self.assertEqual(
                    [
                        s["kind"]
                        for s in design.preview["segments"]
                        if s["kind"] in blocks
                    ],
                    list(order),
                )
                for feature in insert.features:
                    gene = feature.qualifiers["gene"]
                    mapped = next(
                        f
                        for f in design.final_record.features
                        if f.qualifiers.get("gene") == gene
                    )
                    self.assertEqual(
                        mapped.extract(design.final_record.seq),
                        feature.extract(insert.seq),
                    )
                for key in ("backbone", "insert"):
                    prepared = design.preparation_records[key]
                    self.assertEqual(len(EcoRI.search(prepared.seq)), 1)
                    self.assertEqual(len(HindIII.search(prepared.seq)), 1)
                start = design.plan["target"]["replace_start_bp"] - 1
                end = design.plan["target"]["replace_end_bp"]
                outside = str(design.backbone_record.seq[end:]) + str(
                    design.backbone_record.seq[:start]
                )
                self.assertEqual(
                    str(design.preparation_records["backbone"].seq)[6:-6],
                    "AAGCTT" + outside + "GAATTC",
                )
                for feature in design.backbone_record.features:
                    if feature.qualifiers.get("web_kind") not in (
                        ["resistance"],
                        ["replication"],
                        ["terminator"],
                    ):
                        continue
                    mapped = next(
                        f
                        for f in design.preparation_records["backbone"].features
                        if f.qualifiers.get("label") == feature.qualifiers.get("label")
                    )
                    self.assertEqual(
                        mapped.extract(design.preparation_records["backbone"].seq),
                        feature.extract(design.backbone_record.seq),
                    )

    def test_invalid_component_orders_are_rejected(self):
        for order in (
            [],
            ["resistance", "resistance", "expression"],
            ["resistance", "replication", "other"],
            "resistance",
            None,
        ):
            with self.subTest(order=order), self.assertRaises(ValueError):
                self.m.build_design(
                    self.catalog,
                    expression_record(),
                    ["EcoRI", "HindIII"],
                    "basic_seva_ap",
                    "basic_seva_p15a",
                    1,
                    component_order=order,
                )

    def test_reordering_checks_new_module_connection_sites(self):
        from src.plasmid_design.catalog import Module

        resistance = Module("res", "test resistance", "resistance", "TTCCCCC")
        replication = Module("rep", "test replication", "replication", "CCCC")
        catalog = SimpleNamespace(
            get_resistance=lambda _: resistance,
            get_replication=lambda _: replication,
            scaffold={
                "resistance_to_replication": "CCCC",
                "replication_to_t1": "G",
                "t1": "AA",
                "landing_pad_spacer": "CCCC",
            },
        )
        before = self.m.build_design(
            catalog, expression_record(), ["EcoRI", "HindIII"], "res", "rep", 1
        )
        self.assertTrue(before.preview["valid"], before.preview["issues"])
        after = self.m.build_design(
            catalog,
            expression_record(),
            ["EcoRI", "HindIII"],
            "res",
            "rep",
            1,
            component_order=["replication", "resistance", "expression"],
        )
        self.assertFalse(after.preview["valid"])
        self.assertTrue(
            any(
                i.get("enzyme") == "EcoRI"
                and i.get("module") == "模块连接边界"
                and i["start_bp"] == 5
                and i["end_bp"] == 10
                for i in after.preview["issues"]
            )
        )
        self.assertEqual(after.preparation_records, {})

    def test_actual_modules_insert_and_all_three_cds_are_preserved(self):
        insert = expression_record()
        original = deepcopy(insert)
        design = self.m.build_design(
            self.catalog,
            insert,
            ["EcoRI", "HindIII"],
            "basic_seva_ap",
            "basic_seva_p15a",
            1,
        )
        self.assertTrue(design.preview["valid"], design.preview["issues"])
        sequence = str(design.final_record.seq)
        self.assertEqual(sequence.count(str(insert.seq)), 1)
        self.assertTrue(
            sequence.startswith(self.catalog.get_resistance("basic_seva_ap").sequence)
        )
        self.assertIn(
            self.catalog.get_replication("basic_seva_p15a").sequence, sequence
        )
        self.assertEqual(len(EcoRI.search(Seq(sequence), linear=False)), 1)
        self.assertEqual(len(HindIII.search(Seq(sequence), linear=False)), 1)
        inserted = {
            f.qualifiers["gene"][0]: f
            for f in design.final_record.features
            if f.type == "CDS"
            and f.qualifiers.get("gene", [""])[0].startswith("protein_")
        }
        self.assertEqual(set(inserted), {"protein_a", "protein_b", "protein_c"})
        for feature in original.features:
            mapped = inserted[feature.qualifiers["gene"][0]]
            self.assertEqual(
                str(mapped.extract(design.final_record.seq)),
                str(feature.extract(original.seq)),
            )
            self.assertEqual(
                mapped.extract(design.final_record.seq).translate(),
                feature.extract(original.seq).translate(),
            )
        self.assertEqual(str(insert.seq), str(original.seq))
        self.assertEqual(
            [(str(f.location), f.qualifiers) for f in insert.features],
            [(str(f.location), f.qualifiers) for f in original.features],
        )
        self.assertEqual(design.final_record.annotations["topology"], "circular")

    def test_internal_insert_conflict_reports_position_and_keeps_source(self):
        insert = expression_record()
        insert.seq = Seq("GAATTC" + str(insert.seq))
        before = str(insert.seq)
        design = self.m.build_design(
            self.catalog,
            insert,
            ["EcoRI", "HindIII"],
            "basic_seva_ap",
            "basic_seva_p15a",
            1,
        )
        self.assertFalse(design.preview["valid"])
        issues = [i for i in design.preview["issues"] if i.get("enzyme") == "EcoRI"]
        self.assertTrue(any(i["module"] == "完整表达构建" for i in issues))
        self.assertEqual(str(insert.seq), before)
        self.assertEqual(design.preparation_records, {})

    def test_preparation_fragments_have_guard_bases_and_matching_digest_ends(self):
        design = self.m.build_design(
            self.catalog,
            expression_record(),
            ["EcoRI", "HindIII"],
            "basic_seva_ap",
            "basic_seva_p15a",
            1,
        )
        fragments = design.preparation_records
        for key in ["backbone", "insert"]:
            s = fragments[key].seq
            self.assertEqual(len(EcoRI.search(s, linear=True)), 1)
            self.assertEqual(len(HindIII.search(s, linear=True)), 1)
            self.assertGreater(
                min(EcoRI.search(s, linear=True), HindIII.search(s, linear=True))[0], 6
            )
        backbone = str(fragments["backbone"].seq)
        insert = str(fragments["insert"].seq)
        self.assertLess(backbone.index("AAGCTT"), backbone.index("GAATTC"))
        self.assertLess(insert.index("GAATTC"), insert.index("AAGCTT"))
        self.assertIn(str(expression_record().seq), insert)


if __name__ == "__main__":
    unittest.main()
