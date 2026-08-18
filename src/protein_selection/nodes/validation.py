"""Deterministic input validation for the Protein Supply workflow."""

from collections.abc import Mapping
import re
from typing import Any, Protocol

from src.protein_selection.research_context import (
    build_main_enzyme_research_context,
    build_research_context,
)
from src.protein_selection.progress import emit_progress
from src.protein_selection.services.errors import DatabaseLookupError, LookupErrorKind
from src.protein_selection.state import (
    ReactionRecord,
    ProteinSupplyState,
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
    state: ProteinSupplyState,
    *,
    uniprot_client: UniProtLookup,
    kegg_client: ReactionLookup,
) -> ProteinSupplyState:
    """Validate either a legacy pair or one manifest-derived research unit."""

    if state.get("research_unit") is not None:
        return validate_research_unit(
            state,
            uniprot_client=uniprot_client,
            kegg_client=kegg_client,
        )
    return _validate_legacy_input(
        state,
        uniprot_client=uniprot_client,
        kegg_client=kegg_client,
    )


def _validate_legacy_input(
    state: ValidationState,
    *,
    uniprot_client: UniProtLookup,
    kegg_client: ReactionLookup,
) -> ProteinSupplyState:
    """Preserve the original single-UniProt/single-reaction validation path."""

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


def validate_research_unit(
    state: ProteinSupplyState,
    *,
    uniprot_client: UniProtLookup,
    kegg_client: ReactionLookup,
) -> ProteinSupplyState:
    """Validate one manifest-derived main-enzyme research unit."""

    issues: list[ValidationIssue] = []
    emit_progress(
        "validation",
        "started",
        "正在校验主酶研究单元及全部反应 ID",
    )
    raw_unit = state.get("research_unit")
    if not isinstance(raw_unit, Mapping):
        issues.append(
            {
                "code": "missing_input",
                "field": "uniprot_id",
                "message": "research_unit must be an object",
            }
        )
        return _research_unit_failure(
            {},
            "",
            [],
            issues,
            service_error=False,
        )

    research_unit = dict(raw_unit)
    uniprot_id, reaction_ids = _research_unit_identifiers(
        research_unit,
        issues,
    )
    if issues:
        emit_progress(
            "validation",
            "error",
            f"研究单元格式无效，共发现 {len(issues)} 个问题",
        )
        return _research_unit_failure(
            research_unit,
            uniprot_id,
            reaction_ids,
            issues,
            service_error=False,
        )
    emit_progress(
        "validation",
        "completed",
        f"主酶和 {len(reaction_ids)} 个唯一反应 ID 格式正确",
    )

    uniprot_record: UniProtRecord | None = None
    uniprot_annotation: dict[str, Any] | None = None
    reaction_records: dict[str, ReactionRecord] = {}
    service_error = False

    try:
        emit_progress("uniprot", "started", f"正在查询 {uniprot_id}")
        uniprot_record, uniprot_annotation = _lookup_uniprot(
            uniprot_client,
            uniprot_id,
        )
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

    emit_progress(
        "kegg",
        "started",
        f"正在查询 {len(reaction_ids)} 个唯一 KEGG 反应",
    )
    for reaction_id in reaction_ids:
        try:
            record = kegg_client.get_reaction(reaction_id)
            reaction_records[reaction_id] = record
            emit_progress(
                "kegg",
                "info",
                f"查询成功：{reaction_id}",
                verbose_only=True,
            )
        except DatabaseLookupError as exc:
            issue = _lookup_issue("reaction_id", exc)
            issue["message"] = f"{reaction_id}: {issue['message']}"
            issues.append(issue)
            service_error = service_error or exc.kind in {
                LookupErrorKind.SERVICE_ERROR,
                LookupErrorKind.INVALID_RESPONSE,
            }
            emit_progress("kegg", "error", f"{reaction_id}: {exc}")
    if len(reaction_records) == len(reaction_ids):
        emit_progress(
            "kegg",
            "completed",
            f"全部 {len(reaction_ids)} 个 KEGG 反应查询成功",
        )

    if issues:
        emit_progress(
            "validation",
            "error" if service_error else "warning",
            "官方数据库校验未通过，当前研究单元将停止",
        )
        return _research_unit_failure(
            research_unit,
            uniprot_id,
            reaction_ids,
            issues,
            service_error=service_error,
            uniprot_record=uniprot_record,
            uniprot_annotation=uniprot_annotation,
            reaction_records=reaction_records,
        )

    if uniprot_record is None or uniprot_annotation is None:
        raise RuntimeError(
            "successful research-unit validation omitted UniProt records"
        )
    try:
        main_context = build_main_enzyme_research_context(
            uniprot_annotation,
            research_unit,  # type: ignore[arg-type]
            reaction_records,
        )
    except ValueError as exc:
        message = str(exc)
        field: ValidationField = (
            "uniprot_id"
            if "accession" in message.casefold()
            or "uniprot" in message.casefold()
            else "reaction_id"
        )
        issues.append(
            {
                "code": "invalid_response",
                "field": field,
                "message": message,
            }
        )
        emit_progress("validation", "warning", message)
        return _research_unit_failure(
            research_unit,
            uniprot_id,
            reaction_ids,
            issues,
            service_error=False,
            uniprot_record=uniprot_record,
            uniprot_annotation=uniprot_annotation,
            reaction_records=reaction_records,
        )

    is_single_step = main_context.reaction_scope == "single_step"
    legacy_reaction_record = (
        reaction_records[reaction_ids[0]]
        if is_single_step
        else None
    )
    legacy_context = (
        build_research_context(
            uniprot_annotation,
            legacy_reaction_record,
        )
        if legacy_reaction_record is not None
        else None
    )
    match_status = main_context.preliminary_reaction_match
    emit_progress(
        "reaction_match",
        "completed" if match_status == "matched" else "warning",
        (
            f"{match_status}：主酶研究范围包含 "
            f"{len(main_context.reaction_steps)} 个 Step"
        ),
    )
    emit_progress(
        "reaction_match",
        "info",
        main_context.preliminary_reaction_match_reason,
        verbose_only=True,
    )
    emit_progress("validation", "completed", "主酶研究单元校验完成")

    result: ProteinSupplyState = {
        "research_unit": research_unit,  # type: ignore[typeddict-item]
        "uniprot_id": uniprot_record["primary_accession"],
        "validation_status": "valid",
        "validation_errors": [],
        "uniprot_record": uniprot_record,
        "reaction_record": legacy_reaction_record,
        "reaction_records": reaction_records,
        "uniprot_annotation": uniprot_annotation,
        "research_context": (
            legacy_context.model_dump(mode="json")
            if legacy_context is not None
            else None
        ),
        "main_enzyme_research_context": main_context.model_dump(mode="json"),
        "preliminary_reaction_match": match_status,
        "preliminary_reaction_match_reason": (
            main_context.preliminary_reaction_match_reason
        ),
    }
    if is_single_step:
        result["reaction_id"] = reaction_ids[0]
    return result


def _research_unit_identifiers(
    research_unit: Mapping[str, Any],
    issues: list[ValidationIssue],
) -> tuple[str, list[str]]:
    uniprot_id = normalize_identifier(research_unit.get("accession"))
    _validate_required_and_format(
        value=uniprot_id,
        field="uniprot_id",
        pattern=UNIPROT_ACCESSION_PATTERN,
        expected="a 6- or 10-character UniProtKB accession",
        issues=issues,
    )

    scope = research_unit.get("reaction_scope")
    if scope not in {"single_step", "multi_step"}:
        issues.append(
            {
                "code": "invalid_format",
                "field": "reaction_id",
                "message": "reaction_scope must be single_step or multi_step",
            }
        )
    raw_indexes = research_unit.get("assigned_step_indexes")
    indexes: list[int] = []
    if not isinstance(raw_indexes, list):
        issues.append(
            {
                "code": "invalid_format",
                "field": "reaction_id",
                "message": "assigned_step_indexes must be a list",
            }
        )
    else:
        for raw_index in raw_indexes:
            if isinstance(raw_index, bool):
                indexes = []
                break
            try:
                indexes.append(int(raw_index))
            except (TypeError, ValueError):
                indexes = []
                break
        if (
            not indexes
            or any(index < 1 for index in indexes)
            or indexes != sorted(set(indexes))
        ):
            issues.append(
                {
                    "code": "invalid_format",
                    "field": "reaction_id",
                    "message": (
                        "assigned_step_indexes must contain sorted unique "
                        "positive integers"
                    ),
                }
            )

    raw_steps = research_unit.get("reaction_steps")
    step_indexes: list[int] = []
    reaction_ids: list[str] = []
    if not isinstance(raw_steps, list) or any(
        not isinstance(step, Mapping) for step in raw_steps
    ):
        issues.append(
            {
                "code": "invalid_format",
                "field": "reaction_id",
                "message": "reaction_steps must be an object list",
            }
        )
        return uniprot_id, reaction_ids

    for step in raw_steps:
        raw_step_index = step.get("step_index")
        try:
            step_index = int(raw_step_index)
        except (TypeError, ValueError):
            step_index = 0
        if isinstance(raw_step_index, bool) or step_index < 1:
            issues.append(
                {
                    "code": "invalid_format",
                    "field": "reaction_id",
                    "message": "reaction_steps contains an invalid step_index",
                }
            )
        step_indexes.append(step_index)

        reaction_id = normalize_identifier(step.get("reaction_id"))
        before = len(issues)
        _validate_required_and_format(
            value=reaction_id,
            field="reaction_id",
            pattern=KEGG_REACTION_PATTERN,
            expected="a KEGG reaction ID in the form Rxxxxx",
            issues=issues,
        )
        if len(issues) == before and reaction_id not in reaction_ids:
            reaction_ids.append(reaction_id)

    if step_indexes != indexes:
        issues.append(
            {
                "code": "invalid_format",
                "field": "reaction_id",
                "message": (
                    "assigned_step_indexes and reaction_steps must describe "
                    "the same ordered Steps"
                ),
            }
        )
    if scope == "single_step" and len(indexes) != 1:
        issues.append(
            {
                "code": "invalid_format",
                "field": "reaction_id",
                "message": "single_step research units require one Step",
            }
        )
    if scope == "multi_step" and len(indexes) < 2:
        issues.append(
            {
                "code": "invalid_format",
                "field": "reaction_id",
                "message": "multi_step research units require multiple Steps",
            }
        )
    return uniprot_id, reaction_ids


def _lookup_uniprot(
    uniprot_client: UniProtLookup,
    uniprot_id: str,
) -> tuple[UniProtRecord, dict[str, Any]]:
    combined_lookup = getattr(
        uniprot_client,
        "get_record_and_full",
        None,
    )
    if callable(combined_lookup):
        return combined_lookup(uniprot_id)

    record = uniprot_client.get_record(uniprot_id)
    full_lookup = getattr(uniprot_client, "get_full_record", None)
    if callable(full_lookup):
        return record, full_lookup(record["primary_accession"])
    return record, {
        "primaryAccession": record["primary_accession"],
        "uniProtkbId": record["entry_name"],
        "entryType": record["entry_type"],
        "organism": {
            "scientificName": record["organism_name"],
            "taxonId": record["taxon_id"],
        },
    }


def _research_unit_failure(
    research_unit: Mapping[str, Any],
    uniprot_id: str,
    reaction_ids: list[str],
    issues: list[ValidationIssue],
    *,
    service_error: bool,
    uniprot_record: UniProtRecord | None = None,
    uniprot_annotation: dict[str, Any] | None = None,
    reaction_records: dict[str, ReactionRecord] | None = None,
) -> ProteinSupplyState:
    records = reaction_records or {}
    result: ProteinSupplyState = {
        "research_unit": dict(research_unit),  # type: ignore[typeddict-item]
        "uniprot_id": uniprot_id,
        "validation_status": (
            "service_error" if service_error else "invalid"
        ),
        "validation_errors": issues,
        "uniprot_record": uniprot_record,
        "reaction_record": (
            records.get(reaction_ids[0])
            if len(reaction_ids) == 1
            else None
        ),
        "reaction_records": records,
        "uniprot_annotation": uniprot_annotation,
        "research_context": None,
        "main_enzyme_research_context": None,
    }
    if len(reaction_ids) == 1:
        result["reaction_id"] = reaction_ids[0]
    return result


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
