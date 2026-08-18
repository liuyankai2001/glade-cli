"""Deterministic input validation for the Protein Supply workflow."""

import re
from typing import Any, Protocol

from src.protein_selection.research_context import build_research_context
from src.protein_selection.progress import emit_progress
from src.protein_selection.services.errors import DatabaseLookupError, LookupErrorKind
from src.protein_selection.state import (
    ReactionRecord,
    UniProtRecord,
    ValidationField,
    ValidationIssue,
    ValidationState,
)


UNIPROT_ACCESSION_PATTERN = re.compile(
    r"(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})"
)
KEGG_REACTION_PATTERN = re.compile(r"R[0-9]{5}")


class UniProtLookup(Protocol):
    def get_record(self, accession: str) -> UniProtRecord: ...

    def get_full_record(self, accession: str) -> dict[str, Any]: ...


class ReactionLookup(Protocol):
    def get_reaction(self, reaction_id: str) -> ReactionRecord: ...


def normalize_identifier(value: object) -> str:
    """Normalize a user-provided identifier without guessing its meaning."""

    return value.strip().upper() if isinstance(value, str) else ""


def validate_input(
    state: ValidationState,
    *,
    uniprot_client: UniProtLookup,
    kegg_client: ReactionLookup,
) -> ValidationState:
    """Validate input syntax and confirm both records with official databases."""

    uniprot_id = normalize_identifier(state.get("uniprot_id"))
    reaction_id = normalize_identifier(state.get("reaction_id"))
    issues: list[ValidationIssue] = []
    emit_progress("validation", "started", "正在校验输入 ID 格式")

    _validate_required_and_format(
        value=uniprot_id,
        field="uniprot_id",
        pattern=UNIPROT_ACCESSION_PATTERN,
        expected="a 6- or 10-character UniProtKB accession",
        issues=issues,
    )
    _validate_required_and_format(
        value=reaction_id,
        field="reaction_id",
        pattern=KEGG_REACTION_PATTERN,
        expected="a KEGG reaction ID in the form Rxxxxx",
        issues=issues,
    )

    if issues:
        emit_progress(
            "validation",
            "error",
            f"输入格式无效，共发现 {len(issues)} 个问题",
        )
        return {
            "uniprot_id": uniprot_id,
            "reaction_id": reaction_id,
            "validation_status": "invalid",
            "validation_errors": issues,
            "uniprot_record": None,
            "reaction_record": None,
        }
    emit_progress("validation", "completed", "输入 ID 格式正确")

    uniprot_record: UniProtRecord | None = None
    uniprot_annotation: dict[str, Any] | None = None
    reaction_record: ReactionRecord | None = None
    service_error = False

    try:
        emit_progress("uniprot", "started", f"正在查询 {uniprot_id}")
        combined_lookup = getattr(
            uniprot_client,
            "get_record_and_full",
            None,
        )
        if callable(combined_lookup):
            uniprot_record, uniprot_annotation = combined_lookup(uniprot_id)
        else:
            uniprot_record = uniprot_client.get_record(uniprot_id)
            full_lookup = getattr(uniprot_client, "get_full_record", None)
            if callable(full_lookup):
                uniprot_annotation = full_lookup(
                    uniprot_record["primary_accession"]
                )
            else:
                # Lightweight protocol fakes and legacy integrations may only
                # expose the compact lookup. Production UniProtClient always
                # supplies the complete record through get_record_and_full.
                uniprot_annotation = {
                    "primaryAccession": uniprot_record["primary_accession"],
                    "uniProtkbId": uniprot_record["entry_name"],
                    "entryType": uniprot_record["entry_type"],
                    "organism": {
                        "scientificName": uniprot_record["organism_name"],
                        "taxonId": uniprot_record["taxon_id"],
                    },
                }
        emit_progress(
            "uniprot",
            "completed",
            (
                f"查询成功：{uniprot_record['entry_name']}，"
                f"{uniprot_record['organism_name']}"
            ),
        )
    except DatabaseLookupError as exc:
        issues.append(_lookup_issue("uniprot_id", exc))
        service_error = exc.kind in {
            LookupErrorKind.SERVICE_ERROR,
            LookupErrorKind.INVALID_RESPONSE,
        }
        emit_progress("uniprot", "error", str(exc))

    try:
        emit_progress("kegg", "started", f"正在查询 {reaction_id}")
        reaction_record = kegg_client.get_reaction(reaction_id)
        reaction_detail = []
        if reaction_record.get("enzyme_ids"):
            reaction_detail.append(
                "EC " + ", ".join(reaction_record["enzyme_ids"][:3])
            )
        if reaction_record.get("rhea_ids"):
            reaction_detail.append(
                ", ".join(reaction_record["rhea_ids"][:3])
            )
        suffix = "；" + "，".join(reaction_detail) if reaction_detail else ""
        emit_progress(
            "kegg",
            "completed",
            f"查询成功：{reaction_record['reaction_id']}{suffix}",
        )
    except DatabaseLookupError as exc:
        issues.append(_lookup_issue("reaction_id", exc))
        service_error = service_error or exc.kind in {
            LookupErrorKind.SERVICE_ERROR,
            LookupErrorKind.INVALID_RESPONSE,
        }
        emit_progress("kegg", "error", str(exc))

    if issues:
        emit_progress(
            "validation",
            "error" if service_error else "warning",
            "官方数据库校验未通过，工作流将停止",
        )
        return {
            "uniprot_id": uniprot_id,
            "reaction_id": reaction_id,
            "validation_status": "service_error" if service_error else "invalid",
            "validation_errors": issues,
            "uniprot_record": uniprot_record,
            "reaction_record": reaction_record,
            "uniprot_annotation": uniprot_annotation,
            "research_context": None,
        }

    if uniprot_record is None or reaction_record is None:
        raise RuntimeError("successful validation did not produce both records")
    if uniprot_annotation is None:
        raise RuntimeError("successful validation did not produce UniProt annotation")

    research_context = build_research_context(
        uniprot_annotation,
        reaction_record,
    )
    reaction_match_message = (
        "matched：UniProt 与 KEGG 属于同一 Rhea 反应族"
        if research_context.preliminary_reaction_match == "matched"
        else (
            f"{research_context.preliminary_reaction_match}："
            "缺少精确共同 Rhea 映射，需要后续证据核验"
        )
    )
    emit_progress(
        "reaction_match",
        (
            "completed"
            if research_context.preliminary_reaction_match == "matched"
            else "warning"
        ),
        reaction_match_message,
    )
    emit_progress(
        "reaction_match",
        "info",
        research_context.preliminary_reaction_match_reason,
        verbose_only=True,
    )
    emit_progress("validation", "completed", "官方输入记录校验完成")

    return {
        "uniprot_id": uniprot_record["primary_accession"],
        "reaction_id": reaction_record["reaction_id"],
        "validation_status": "valid",
        "validation_errors": [],
        "uniprot_record": uniprot_record,
        "reaction_record": reaction_record,
        "uniprot_annotation": uniprot_annotation,
        "research_context": research_context.model_dump(mode="json"),
        "preliminary_reaction_match": (
            research_context.preliminary_reaction_match
        ),
        "preliminary_reaction_match_reason": (
            research_context.preliminary_reaction_match_reason
        ),
    }


def _validate_required_and_format(
    *,
    value: str,
    field: ValidationField,
    pattern: re.Pattern[str],
    expected: str,
    issues: list[ValidationIssue],
) -> None:
    if not value:
        issues.append(
            {"code": "missing_input", "field": field, "message": "value is required"}
        )
    elif pattern.fullmatch(value) is None:
        issues.append(
            {
                "code": "invalid_format",
                "field": field,
                "message": f"expected {expected}",
            }
        )


def _lookup_issue(field: ValidationField, exc: DatabaseLookupError) -> ValidationIssue:
    return {"code": exc.kind.value, "field": field, "message": str(exc)}
