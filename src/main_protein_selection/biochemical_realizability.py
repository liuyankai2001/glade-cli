from __future__ import annotations

from typing import Any

from src.main_protein_selection.reaction_direction_verifier import (
    DIRECTION_CONTRADICTED,
    DIRECTION_SUPPORTED,
    DIRECTION_UNKNOWN,
)


REACTION_FIT_VERIFIED = "verified"
REACTION_FIT_VERIFIED_WITH_RISK = "verified_with_risk"
REACTION_FIT_MANUAL_REVIEW = "manual_review"
REACTION_FIT_REJECTED = "rejected"


# These rules describe directionally unsupported uses of otherwise correctly
# annotated EC classes.  They are keyed by biochemical function, never by
# benchmark target or literature-gold accession.
REVERSE_DIRECTION_BLOCK_RULES: dict[str, dict[str, str]] = {
    "1.2.1.67": {
        "rule_id": "reverse_vanillin_dehydrogenase_not_supported",
        "evidence": (
            "EC 1.2.1.67 annotations describe aldehyde oxidation to vanillate; "
            "the reverse carboxylate-to-aldehyde realization requires a different enzyme class"
        ),
        "repair_class": "carboxylate_to_aldehyde_reduction",
        "suggested_enzyme_family": "carboxylic acid reductase",
        "required_cofactors": "ATP;NADPH;4'-phosphopantetheine",
    },
    "2.1.1.341": {
        "rule_id": "reverse_vanillate_demethylase_not_supported",
        "evidence": (
            "EC 2.1.1.341 annotations describe tetrahydrofolate-dependent O-demethylation; "
            "the reverse biosynthetic methylation is not supported by that enzyme class"
        ),
        "repair_class": "biosynthetic_o_methylation",
        "suggested_enzyme_family": "SAM-dependent O-methyltransferase",
        "required_cofactors": "S-adenosyl-L-methionine;Mg2+",
    },
}


def _split(value: Any) -> list[str]:
    if value is None:
        return []
    if isinstance(value, (list, tuple, set)):
        return [str(item).strip() for item in value if str(item).strip()]
    return [
        item.strip()
        for item in str(value).replace("|", ";").split(";")
        if item.strip()
    ]


def _normalized_ecs(value: Any) -> list[str]:
    return [item.strip() for item in _split(value) if item.strip()]


def _requirement_ecs(requirement: dict[str, Any]) -> list[str]:
    ecs = _normalized_ecs(requirement.get("ec_numbers"))
    if ecs:
        return ecs
    return _normalized_ecs(
        requirement.get("locked_ec_numbers")
        or requirement.get("locked_enzyme_ecs")
    )


def ec_status(requirement: dict[str, Any]) -> str:
    explicit = str(requirement.get("ec_status") or "").strip().lower()
    if explicit in {"complete", "partial", "missing"}:
        return explicit
    ecs = _requirement_ecs(requirement)
    if not ecs:
        return "missing"
    if any("-" in ec or len(ec.split(".")) != 4 for ec in ecs):
        return "partial"
    return "complete"


def evaluate_candidate_reaction_fit(
    requirement: dict[str, Any],
    candidate: dict[str, Any],
) -> dict[str, Any]:
    status = ec_status(requirement)
    candidate_rhea = set(_normalized_ecs(candidate.get("matched_rhea_ids")))
    requirement_rhea = set(_normalized_ecs(requirement.get("rhea_ids")))
    requirement_rhea.update(_normalized_ecs(requirement.get("rhea_master_ids")))
    requirement_rhea.update(_normalized_ecs(requirement.get("required_rhea_direction_ids")))
    requirement_rhea.update(_normalized_ecs(requirement.get("rhea_bidirectional_ids")))
    exact_rhea_match = bool(candidate_rhea & requirement_rhea)
    candidate_ko = set(_normalized_ecs(candidate.get("matched_ko_ids")))
    requirement_ko = set(_normalized_ecs(requirement.get("ko_ids")))
    retrieval_strategies = {
        value.lower()
        for value in (
            _split(candidate.get("retrieval_strategy"))
            + _split(candidate.get("retrieval_strategies"))
        )
    }
    exact_ko_match = (
        bool(candidate_ko & requirement_ko)
        and "kegg_ko_exact" in retrieval_strategies
    )
    reaction_confidence = {
        value.lower() for value in _split(candidate.get("reaction_confidence"))
    }
    exact_selenzyme_match = (
        "selenzyme_kegg_exact" in retrieval_strategies
        and "selenzyme_exact" in reaction_confidence
    )
    risk_selenzyme_match = (
        "selenzyme_kegg_risk" in retrieval_strategies
        and "selenzyme_risk" in reaction_confidence
    )
    ec_selenzyme_match = (
        "selenzyme_ec_risk" in retrieval_strategies
        and "selenzyme_ec_risk" in reaction_confidence
    )
    selenzyme_ec_relations = set(_split(candidate.get("selenzyme_risk_status")))
    shared_reaction_ec_overlap = (
        "shared_reaction_ec_overlap" in selenzyme_ec_relations
    )
    literature_grade = (
        "a"
        if "literature_grade_a" in reaction_confidence
        else "b"
        if "literature_grade_b" in reaction_confidence
        else ""
    )
    literature_activity_match = (
        "literature_experimental_activity" in retrieval_strategies
        and bool(literature_grade)
    )
    candidate_ec_value = (
        candidate.get("ec_numbers")
        if "ec_numbers" in candidate
        else candidate.get("ec_number")
    )
    candidate_declared_ecs = set(_normalized_ecs(candidate_ec_value))
    requirement_declared_ecs = set(_requirement_ecs(requirement))
    exact_ec_match = (
        status == "complete"
        and (
            "ec_exact" in retrieval_strategies
            or (
                not retrieval_strategies
                and bool(candidate_declared_ecs & requirement_declared_ecs)
            )
        )
    )
    if (
        status != "complete"
        and not exact_rhea_match
        and not exact_ko_match
        and not literature_activity_match
        and not exact_selenzyme_match
        and not ec_selenzyme_match
        and not risk_selenzyme_match
    ):
        return {
            "status": REACTION_FIT_MANUAL_REVIEW,
            "score": 0.0,
            "rule_ids": [f"{status}_ec_requires_fallback"],
            "evidence": [
                f"Reaction has {status} EC annotation and lacks an exact reaction-backed candidate"
            ],
        }

    direction = str(requirement.get("direction") or "").strip().lower()
    candidate_ecs = set(_normalized_ecs(candidate.get("ec_numbers")))
    matched_rules = [
        rule
        for ec, rule in REVERSE_DIRECTION_BLOCK_RULES.items()
        if direction == "right_to_left" and ec in candidate_ecs
    ]
    if matched_rules:
        return {
            "status": REACTION_FIT_REJECTED,
            "score": 0.0,
            "rule_ids": [rule["rule_id"] for rule in matched_rules],
            "evidence": [rule["evidence"] for rule in matched_rules],
        }

    if (
        status == "complete"
        and (exact_selenzyme_match or ec_selenzyme_match or risk_selenzyme_match)
        and candidate_declared_ecs
        and requirement_declared_ecs
        and not (candidate_declared_ecs & requirement_declared_ecs)
        and not (ec_selenzyme_match and shared_reaction_ec_overlap)
    ):
        return {
            "status": REACTION_FIT_REJECTED,
            "score": 0.0,
            "rule_ids": ["selenzyme_candidate_ec_contradicts_requirement"],
            "evidence": [
                "Selenzyme candidate has an official EC annotation incompatible with the locked complete EC",
                f"candidate ECs: {';'.join(sorted(candidate_declared_ecs))}",
                f"required ECs: {';'.join(sorted(requirement_declared_ecs))}",
            ],
        }

    direction_verdict = str(candidate.get("direction_verdict") or "").strip().lower()
    direction_confidence = str(candidate.get("direction_confidence") or "").strip().lower()
    direction_evidence = _split(candidate.get("direction_evidence"))
    if direction_verdict == DIRECTION_CONTRADICTED:
        return {
            "status": REACTION_FIT_REJECTED,
            "score": 0.0,
            "rule_ids": ["protein_direction_unsupported"],
            "evidence": direction_evidence or [
                "Candidate official evidence contradicts the locked reaction direction"
            ],
        }

    def direction_checked(
        result: dict[str, Any],
        *,
        supported_rule: str,
    ) -> dict[str, Any]:
        if direction_verdict == DIRECTION_UNKNOWN:
            return {
                **result,
                "status": REACTION_FIT_VERIFIED_WITH_RISK,
                "score": min(float(result.get("score") or 0), 69.0),
                "rule_ids": [*result.get("rule_ids", []), "direction_unknown_risk"],
                "evidence": [
                    *result.get("evidence", []),
                    *(direction_evidence or [
                        "Official evidence does not establish the locked reaction direction"
                    ]),
                ],
            }
        if direction_verdict == DIRECTION_SUPPORTED:
            return {
                **result,
                "rule_ids": [*result.get("rule_ids", []), supported_rule],
                "evidence": [
                    *result.get("evidence", []),
                    *(direction_evidence or [
                        f"Direction support established ({direction_confidence or 'unspecified'} confidence)"
                    ]),
                ],
            }
        # Test/legacy records that lack the new context retain their historical
        # classification. Production requirements are enriched before this
        # function is called.
        return result

    if exact_rhea_match:
        return direction_checked({
            "status": REACTION_FIT_VERIFIED,
            "score": 100.0,
            "rule_ids": ["rhea_exact_reaction_match"],
            "evidence": [
                "Candidate is annotated to the Rhea reaction mapped from the locked KEGG reaction",
            ],
        }, supported_rule="rhea_direction_supported")

    if exact_ko_match:
        matched_ko_ids = sorted(candidate_ko & requirement_ko)
        return direction_checked({
            "status": REACTION_FIT_VERIFIED,
            "score": 95.0,
            "rule_ids": ["kegg_ko_exact_annotation"],
            "evidence": [
                "Candidate is mapped from an exact KEGG Orthology annotation",
                f"matched KO IDs: {';'.join(matched_ko_ids)}",
            ],
        }, supported_rule="ko_candidate_direction_supported")

    if exact_ec_match:
        catalytic = str(candidate.get("catalytic_activities") or "").strip()
        rhea_ids = _split(candidate.get("rhea_ids"))
        evidence = ["Exact complete EC with no known direction/mechanism contradiction"]
        score = 70.0
        if catalytic:
            evidence.append("UniProt catalytic activity annotation is present")
            score += 15.0
        if rhea_ids:
            evidence.append("Candidate has a Rhea reaction cross-reference")
            score += 15.0
        return direction_checked({
            "status": REACTION_FIT_VERIFIED,
            "score": min(score, 100.0),
            "rule_ids": ["exact_ec_no_direction_contradiction"],
            "evidence": evidence,
        }, supported_rule="ec_candidate_direction_supported")

    if literature_activity_match:
        level = literature_grade.upper()
        return direction_checked({
            "status": REACTION_FIT_VERIFIED_WITH_RISK,
            "score": 90.0 if literature_grade == "a" else 80.0,
            "rule_ids": [
                f"literature_experimental_activity_grade_{literature_grade}"
            ],
            "evidence": [
                f"Literature Grade {level} experimental activity supports the locked reaction",
                "Literature-derived non-standard activity requires human review",
            ],
        }, supported_rule="literature_candidate_direction_supported")

    if exact_selenzyme_match:
        return direction_checked({
            "status": REACTION_FIT_VERIFIED,
            "score": 100.0,
            "rule_ids": ["selenzyme_exact_reaction_no_direction_contradiction"],
            "evidence": [
                "Candidate has a unit-similarity SelenzymeRF match to the locked KEGG reaction",
            ],
        }, supported_rule="selenzyme_candidate_direction_supported")

    if ec_selenzyme_match:
        current_uniprot_ec_confirmed = bool(
            candidate_declared_ecs & requirement_declared_ecs
        )
        if current_uniprot_ec_confirmed:
            ec_relation = "exact_current_uniprot_ec"
        elif shared_reaction_ec_overlap:
            ec_relation = "shared_reaction_ec_overlap"
        elif not candidate_declared_ecs:
            ec_relation = "unannotated_current_uniprot_ec"
        else:
            return {
                "status": REACTION_FIT_REJECTED,
                "score": 0.0,
                "rule_ids": ["selenzyme_candidate_ec_contradicts_requirement"],
                "evidence": [
                    "Candidate EC has no exact or shared-reaction relationship with the locked EC"
                ],
            }
        score = {
            "exact_current_uniprot_ec": 69.0,
            "shared_reaction_ec_overlap": 60.0,
            "unannotated_current_uniprot_ec": 55.0,
        }[ec_relation]
        relation_evidence = {
            "exact_current_uniprot_ec": (
                "Current UniProt annotation contains the locked EC"
            ),
            "shared_reaction_ec_overlap": (
                "Candidate and locked EC numbers are both attached to the same Selenzyme reaction record"
            ),
            "unannotated_current_uniprot_ec": (
                "Current UniProt annotation has no EC number confirming the Selenzyme association"
            ),
        }[ec_relation]
        return direction_checked({
            "status": REACTION_FIT_VERIFIED_WITH_RISK,
            "score": score,
            "rule_ids": [
                "selenzyme_ec_association_requires_target_reaction_review",
                f"selenzyme_ec_relation_{ec_relation}",
            ],
            "evidence": [
                "SelenzymeRF associates the protein with the locked complete EC",
                "The EC query uses a representative EC reaction and does not establish the locked substrate/product specificity",
                relation_evidence,
            ],
        }, supported_rule="selenzyme_ec_candidate_direction_supported")

    if risk_selenzyme_match:
        try:
            similarity = float(candidate.get("selenzyme_reaction_similarity"))
        except (TypeError, ValueError):
            similarity = -1.0
        if not 0.0 <= similarity < 1.0:
            return {
                "status": REACTION_FIT_REJECTED,
                "score": 0.0,
                "rule_ids": ["selenzyme_invalid_combined_reaction_similarity"],
                "evidence": [
                    "SelenzymeRF risk fallback lacks a valid combined reaction-similarity score"
                ],
            }
        return direction_checked({
            "status": REACTION_FIT_VERIFIED_WITH_RISK,
            "score": round(similarity * 100.0, 6),
            "rule_ids": ["selenzyme_combined_similarity_risk_fallback"],
            "evidence": [
                "Candidate is accepted by the SelenzymeRF risk-fallback policy",
                f"combined reaction similarity: {similarity:.10g}",
            ],
        }, supported_rule="selenzyme_candidate_direction_supported")

    return {
        "status": REACTION_FIT_MANUAL_REVIEW,
        "score": 0.0,
        "rule_ids": ["candidate_lacks_supported_reaction_evidence"],
        "evidence": [
            "Candidate lacks exact EC/Rhea/KO, Grade A/B literature, or Selenzyme reaction evidence"
        ],
    }


def candidate_is_reaction_verified(row: dict[str, Any]) -> bool:
    return str(row.get("reaction_fit_status") or "").strip() in {
        REACTION_FIT_VERIFIED,
        REACTION_FIT_VERIFIED_WITH_RISK,
    }


def route_repair_requests_for_requirement(requirement: dict[str, Any]) -> list[dict[str, Any]]:
    direction = str(requirement.get("direction") or "").strip().lower()
    if direction != "right_to_left":
        return []
    requests: list[dict[str, Any]] = []
    for ec in _requirement_ecs(requirement):
        rule = REVERSE_DIRECTION_BLOCK_RULES.get(ec)
        if not rule:
            continue
        requests.append({
            "step_index": int(requirement.get("step_index") or 0),
            "reaction_id": str(requirement.get("reaction_id") or ""),
            "locked_direction": direction,
            "blocked_ec": ec,
            "blocking_rule_id": rule["rule_id"],
            "evidence": rule["evidence"],
            "repair_class": rule.get("repair_class", "replace_biochemical_reaction"),
            "suggested_enzyme_family": rule.get("suggested_enzyme_family", ""),
            "required_cofactors": _split(rule.get("required_cofactors", "")),
            "status": "proposal_only",
            "recommended_action": "replace reaction stoichiometry and rerun GEM/FBA/FVA before enzyme reselection",
            "requires_gem_revalidation": True,
        })
    return requests
