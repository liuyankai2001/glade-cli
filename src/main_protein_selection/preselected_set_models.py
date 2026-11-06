"""V8 set validation: fixed reaction groups plus scored automatic assignments."""
from __future__ import annotations

from math import isclose
from typing import Literal

from pydantic import Field, model_validator

from .models import (
    MAIN_ENZYME_SETS_SCHEMA_VERSION, MAIN_ENZYME_SETS_ALGORITHM_VERSION,
    MainEnzymeSet, MainEnzymeSetsResult, PreselectedReaction, _require_sorted_unique_positive,
)

PRESELECTED_SETS_ALGORITHM_VERSION = "evidence_constrained_assignment.v7_preselected_reactions"


def metadata_metrics(proteins):
    return {"protein_count": len(proteins),
            "organism_count": len({p["organism_name"] for p in proteins if p["organism_name"]}),
            "reviewed_fraction": sum(bool(p["reviewed"]) for p in proteins) / len(proteins)}


class PreselectedEnzymeSet(MainEnzymeSet):
    def automatic_projection(self):
        result = self.model_dump(mode="json")
        indexes = {assignment.step_index for assignment in self.step_assignments}
        accessions = {assignment.accession for assignment in self.step_assignments}
        result["covered_step_indexes"] = sorted(indexes)
        result["review_required_step_indexes"] = sorted(set(self.review_required_step_indexes) & indexes)
        proteins = []
        for protein in result["proteins"]:
            if protein["accession"] not in accessions:
                continue
            protein["assigned_step_indexes"] = sorted(set(protein["assigned_step_indexes"]) & indexes)
            protein["capable_step_indexes"] = sorted(set(protein["capable_step_indexes"]) & indexes)
            proteins.append(protein)
        result["proteins"] = proteins
        result["protein_count"] = len(proteins)
        if proteins:
            result["metrics"].update(metadata_metrics(proteins))
        result["status"] = "review_required" if result["review_required_step_indexes"] or result["warnings"] else "complete"
        return result

    @model_validator(mode="after")
    def validate_consistency(self):
        if not self.step_assignments:
            raise ValueError("all-fixed selections must use already_selected, without a generated set")
        for name in ("covered_step_indexes", "uncovered_step_indexes", "review_required_step_indexes"):
            _require_sorted_unique_positive(getattr(self, name), name)
        if not self.coverage_complete or self.uncovered_step_indexes:
            raise ValueError("generated v8 sets must have complete coverage")
        if not set(self.review_required_step_indexes).issubset(self.covered_step_indexes):
            raise ValueError("review indexes must be covered")
        accessions = [protein.accession for protein in self.proteins]
        if accessions != sorted(set(accessions)) or self.protein_count != len(accessions):
            raise ValueError("v8 set proteins must be sorted, unique and counted in full")
        expected = metadata_metrics([protein.model_dump() for protein in self.proteins])
        for name, value in expected.items():
            if not isclose(getattr(self.metrics, name), value):
                raise ValueError(f"{name} must describe the complete selected protein set")
        if self.protein_count != self.metrics.protein_count:
            raise ValueError("protein count and metrics differ")
        if self.metrics.warning_count != len(self.warnings):
            raise ValueError("warning count differs")
        needs_review = bool(self.review_required_step_indexes or self.warnings)
        if (self.status == "review_required") != needs_review:
            raise ValueError("v8 set review status differs from its review indicators")
        MainEnzymeSet.model_validate(self.automatic_projection())
        return self


class PreselectedEnzymeSetsResult(MainEnzymeSetsResult):
    schema_version: Literal["main_enzyme_sets.v8"] = "main_enzyme_sets.v8"
    algorithm_version: Literal["evidence_constrained_assignment.v7_preselected_reactions"] = PRESELECTED_SETS_ALGORITHM_VERSION
    status: Literal["complete", "review_required", "infeasible", "truncated", "stale_input", "source_unavailable", "already_selected"]
    sets: list[PreselectedEnzymeSet] = Field(default_factory=list)
    preselected_reactions: dict[str, PreselectedReaction]
    preselected_review_step_indexes: list[int] = Field(default_factory=list)

    @model_validator(mode="after")
    def validate_consistency(self):
        if not self.preselected_reactions:
            raise ValueError("v8 set results require preselected reactions")
        if any(not reaction.strip() for reaction in self.preselected_reactions):
            raise ValueError("fixed reaction IDs must not be empty")
        _require_sorted_unique_positive(self.required_step_indexes, "required_step_indexes")
        _require_sorted_unique_positive(self.preselected_review_step_indexes, "preselected_review_step_indexes")
        fixed = {}
        step_list = []
        hashes = {}
        for group in self.preselected_reactions.values():
            step_list.extend(group.step_indexes)
            for protein in group.proteins:
                fixed.setdefault(protein.accession, set()).update(group.step_indexes)
                if hashes.setdefault(protein.accession, protein.sequence_sha256) != protein.sequence_sha256:
                    raise ValueError("conflicting fixed sequences")
        if len(step_list) != len(set(step_list)):
            raise ValueError("fixed reaction groups overlap")
        fixed_steps = set(step_list)
        required = set(self.required_step_indexes)
        if not fixed_steps.issubset(required) or not set(self.preselected_review_step_indexes).issubset(fixed_steps):
            raise ValueError("fixed steps or reviews do not match the required route")
        if self.genome_ec_filter is not None:
            if (required != set(self.genome_ec_filter.external_required_step_indexes) | fixed_steps
                    or self.chassis_key != self.genome_ec_filter.chassis_key):
                raise ValueError("genome coverage does not match v8 requirements")
        automatic = required - fixed_steps
        if self.status == "already_selected":
            if (automatic or self.sets or self.uncovered_step_indexes or not self.ok or not self.search_complete
                    or self.search_nodes != 0 or self.minimum_protein_count != len(fixed)):
                raise ValueError("invalid already_selected result")
            return self
        if not automatic:
            raise ValueError("all-fixed successful input must report already_selected")
        for selected in self.sets:
            if set(selected.covered_step_indexes) != required:
                raise ValueError("v8 set must cover fixed and automatic route steps")
            if {item.step_index for item in selected.step_assignments} != automatic:
                raise ValueError("only remaining steps may have automatic assignments")
            assignments = {accession: set(indexes) for accession, indexes in fixed.items()}
            for item in selected.step_assignments:
                assignments.setdefault(item.accession, set()).add(item.step_index)
            actual = {p.accession: set(p.assigned_step_indexes) for p in selected.proteins}
            if actual != assignments:
                raise ValueError("v8 set must retain every fixed protein and every automatic assignment")
            if any(not set(protein.capable_step_indexes).issubset(required) for protein in selected.proteins):
                raise ValueError("protein capabilities exceed route steps")
            if not set(self.preselected_review_step_indexes).issubset(selected.review_required_step_indexes):
                raise ValueError("fixed review steps must be retained")
        projection = self.model_dump(mode="json", exclude={"preselected_reactions", "preselected_review_step_indexes"})
        projection.update(schema_version=MAIN_ENZYME_SETS_SCHEMA_VERSION, algorithm_version=MAIN_ENZYME_SETS_ALGORITHM_VERSION,
                          genome_ec_filter=None, required_step_indexes=sorted(automatic),
                          sets=[selected.automatic_projection() for selected in self.sets])
        projection["minimum_protein_count"] = min((item["protein_count"] for item in projection["sets"]), default=None)
        if self.status in {"complete", "review_required"} and self.sets:
            if self.sets[0].status != self.status:
                raise ValueError("result status must describe the complete top set")
            projection["status"] = projection["sets"][0]["status"]
        MainEnzymeSetsResult.model_validate(projection)
        if self.sets:
            if self.minimum_protein_count is None or not len(fixed) <= self.minimum_protein_count <= min(s.protein_count for s in self.sets):
                raise ValueError("minimum count must use the complete protein union")
        elif self.minimum_protein_count is not None:
            raise ValueError("failed search must not report a minimum count")
        return self
