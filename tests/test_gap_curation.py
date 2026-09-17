import unittest


class ReviewedGapExtractionTests(unittest.TestCase):
    def test_reviewed_flanks_are_removed_and_cds_coordinates_preserved(self):
        from src.plasmid_design.gap_curation import extract_reviewed_gaps

        features = [
            {
                "type": "CDS",
                "parts": [{"start": 2, "end": 11, "strand": 1}],
                "operator": None,
                "qualifiers": {"gene": ["test"]},
            }
        ]
        gaps = [
            {
                "start_1based": 1,
                "end_1based": 2,
                "evidence": "reviewed upstream boundary",
            },
            {
                "start_1based": 12,
                "end_1based": 13,
                "evidence": "reviewed downstream boundary",
            },
        ]
        clean, mapped, retained = extract_reviewed_gaps("AAATGCCCTAAGG", features, gaps)
        self.assertEqual(clean, "ATGCCCTAA")
        self.assertEqual(mapped[0]["parts"], [{"start": 0, "end": 9, "strand": 1}])
        self.assertEqual(
            retained,
            [
                {
                    "original_start_1based": 3,
                    "original_end_1based": 11,
                    "module_start_1based": 1,
                    "module_end_1based": 9,
                }
            ],
        )
        self.assertEqual(features[0]["parts"][0]["start"], 2)

    def test_ambiguous_or_function_overlapping_deletions_are_rejected(self):
        from src.plasmid_design.gap_curation import extract_reviewed_gaps

        features = [
            {
                "type": "rep_origin",
                "parts": [{"start": 2, "end": 6, "strand": 1}],
                "qualifiers": {},
            }
        ]
        for gaps in (
            [{"start_1based": 1, "end_1based": 2}],
            [
                {
                    "start_1based": 3,
                    "end_1based": 4,
                    "evidence": "unreviewed functional region",
                }
            ],
            [{"start_1based": 0, "end_1based": 2, "evidence": "invalid range"}],
            [
                {"start_1based": 1, "end_1based": 2, "evidence": "one"},
                {"start_1based": 2, "end_1based": 2, "evidence": "overlap"},
            ],
        ):
            with self.subTest(gaps=gaps), self.assertRaises(ValueError):
                extract_reviewed_gaps("AAACGTGG", features, gaps)
