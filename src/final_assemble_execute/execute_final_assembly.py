"""Execute all selected final-assembly plans as an in-silico output bundle."""

from __future__ import annotations

import json
import shutil
import tempfile
import uuid
from collections import Counter
from collections.abc import Callable, Mapping
from dataclasses import dataclass, field
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from src.final_assemble_execute.common import (
    assemble_sequence,
    build_final_record,
    load_execution_context,
    sha256_file,
    sha256_sequence,
    stable_json_hash,
    validate_written_outputs,
    write_fasta,
    write_genbank,
    write_json_atomic,
)
from src.final_assemble_execute.config import (
    FINAL_ASSEMBLY_DIRNAME,
    FINAL_ASSEMBLY_EXECUTION_ALGORITHM_VERSION,
    FINAL_ASSEMBLY_SCHEMA_VERSION,
    FINAL_DESIGN_REPORT_FILENAME,
    FINAL_DESIGN_REPORT_SCHEMA_VERSION,
    RUN_SUMMARY_FILENAME,
    THEORETICAL_ASSEMBLY_WARNING,
)
from src.final_assemble_execute.export_final_design_report import (
    ReportGenerator,
    build_report_payload,
    generate_final_design_report,
)
from src.final_assemble_execute.models import FinalAssemblyExecutionContext
from src.write_manifest.store import update_design_manifest


DesignExecutor = Callable[
    [FinalAssemblyExecutionContext, Mapping[str, Any], Path],
    dict[str, Any],
]


def _utc_now() -> str:
    return datetime.now(UTC).isoformat().replace("+00:00", "Z")


def _final_relative(filename: str) -> str:
    return f"{FINAL_ASSEMBLY_DIRNAME}/{filename}"


def _staged_file_info(
    path: Path,
    filename: str,
    *,
    file_format: str,
    sequence: str | None = None,
) -> dict[str, Any]:
    result: dict[str, Any] = {
        "path": _final_relative(filename),
        "format": file_format,
        "file_sha256": sha256_file(path),
    }
    if sequence is not None:
        result.update(
            {
                "length_bp": len(sequence),
                "sequence_sha256": sha256_sequence(sequence),
            }
        )
    return result


def _execute_one_design(
    context: FinalAssemblyExecutionContext,
    plan: Mapping[str, Any],
    staging_dir: Path,
) -> dict[str, Any]:
    design_id = int(plan["parts_design_id"])
    construct = next(
        item
        for item in context.source_context.constructs
        if item.design_id == design_id
    )
    insert_record = context.insert_records[design_id]
    assembled = assemble_sequence(
        context.source_context.backbone.sequence,
        construct.sequence,
        plan,
    )
    record_id = f"{context.source_context.target_compound_id}_D{design_id:03d}_FIN"
    final_record, feature_warnings = build_final_record(
        record_id=record_id,
        backbone_record=context.backbone_record,
        insert_record=insert_record,
        result=assembled,
        plan=plan,
    )
    prefix = f"design_{design_id:03d}_final"
    genbank_name = f"{prefix}.gb"
    fasta_name = f"{prefix}.fasta"
    assembly_name = f"design_{design_id:03d}_assembly.json"
    genbank_path = staging_dir / genbank_name
    fasta_path = staging_dir / fasta_name
    assembly_path = staging_dir / assembly_name
    write_genbank(genbank_path, final_record)
    write_fasta(fasta_path, record_id, assembled.sequence)
    validation = validate_written_outputs(
        genbank_path=genbank_path,
        fasta_path=fasta_path,
        expected_sequence=assembled.sequence,
        inserted_start_bp=assembled.inserted_start_bp,
        inserted_end_bp=assembled.inserted_end_bp,
        expected_insert=construct.sequence,
    )
    plan_warnings = plan.get("warnings")
    plan_warnings = plan_warnings if isinstance(plan_warnings, list) else []
    warnings_list = [str(item) for item in [*plan_warnings, *feature_warnings]]
    linearization = plan.get("backbone_linearization")
    linearization = linearization if isinstance(linearization, Mapping) else {}
    files = {
        "genbank": _staged_file_info(
            genbank_path,
            genbank_name,
            file_format="genbank",
            sequence=assembled.sequence,
        ),
        "fasta": _staged_file_info(
            fasta_path,
            fasta_name,
            file_format="fasta",
            sequence=assembled.sequence,
        ),
    }
    report_payload = {
        "schema_version": "final_assembly_design_report.v1",
        "status": "assembled_in_silico",
        "algorithm_version": FINAL_ASSEMBLY_EXECUTION_ALGORITHM_VERSION,
        "target_compound_id": context.source_context.target_compound_id,
        "parts_design_id": design_id,
        "source": {
            "final_assembly_plan_selection_fingerprint": (
                context.plan_selection_fingerprint
            ),
            "plan_fingerprint": plan.get("plan_fingerprint"),
            "backbone_file_sha256": context.source_context.backbone.file_sha256,
            "insert_file_sha256": construct.file_sha256,
        },
        "assembly_method": plan.get("assembly_method"),
        "backbone_linearization": dict(linearization),
        "target": assembled.target_audit,
        "final_construct": {
            "record_id": record_id,
            "length_bp": len(assembled.sequence),
            "topology": "circular",
            "sequence_sha256": sha256_sequence(assembled.sequence),
            "feature_count": len(final_record.features),
            "files": files,
        },
        "validation": validation,
        "warnings": [THEORETICAL_ASSEMBLY_WARNING, *warnings_list],
    }
    write_json_atomic(assembly_path, report_payload)
    files["assembly_report"] = _staged_file_info(
        assembly_path,
        assembly_name,
        file_format="json",
    )
    return {
        "parts_design_id": design_id,
        "status": "assembled_in_silico",
        "assembly_method": plan.get("assembly_method"),
        "enzyme_summary": linearization.get("enzyme_summary"),
        "plan_fingerprint": plan.get("plan_fingerprint"),
        "target": assembled.target_audit,
        "length_bp": len(assembled.sequence),
        "topology": "circular",
        "sequence_sha256": sha256_sequence(assembled.sequence),
        "feature_count": len(final_record.features),
        "files": files,
        "validation": validation,
        "warnings": warnings_list,
    }


def _remove_transaction_path(path: Path, project_root: Path, prefix: str) -> None:
    resolved = path.resolve()
    root = project_root.resolve()
    if resolved.parent != root or not resolved.name.startswith(prefix):
        raise ValueError(f"refusing to remove unsafe transaction path: {resolved}")
    if resolved.exists():
        shutil.rmtree(resolved)


@dataclass(slots=True)
class FinalAssemblyOutputTransaction:
    project_root: Path
    target_dir: Path
    staging_dir: Path
    backup_dir: Path | None = field(default=None, init=False)
    installed: bool = field(default=False, init=False)

    def install(self) -> None:
        root = self.project_root.resolve()
        target = self.target_dir.resolve()
        if target.parent != root or target.name != FINAL_ASSEMBLY_DIRNAME:
            raise ValueError(f"unsafe final assembly target: {target}")
        if self.target_dir.exists():
            self.backup_dir = root / (
                f".{FINAL_ASSEMBLY_DIRNAME}.backup-{uuid.uuid4().hex}"
            )
            self.target_dir.rename(self.backup_dir)
        try:
            self.staging_dir.rename(self.target_dir)
        except Exception:
            if self.backup_dir is not None and self.backup_dir.exists():
                self.backup_dir.rename(self.target_dir)
            raise
        self.installed = True

    def rollback(self) -> None:
        if not self.installed:
            self.cleanup_staging()
            return
        failed = self.project_root / (
            f".{FINAL_ASSEMBLY_DIRNAME}.failed-{uuid.uuid4().hex}"
        )
        if self.target_dir.exists():
            self.target_dir.rename(failed)
        if self.backup_dir is not None and self.backup_dir.exists():
            self.backup_dir.rename(self.target_dir)
        _remove_transaction_path(
            failed,
            self.project_root,
            f".{FINAL_ASSEMBLY_DIRNAME}.failed-",
        )
        self.installed = False
        self.backup_dir = None
        self.cleanup_staging()

    def finalize(self) -> str | None:
        warning: str | None = None
        if self.backup_dir is not None and self.backup_dir.exists():
            try:
                _remove_transaction_path(
                    self.backup_dir,
                    self.project_root,
                    f".{FINAL_ASSEMBLY_DIRNAME}.backup-",
                )
            except OSError as exc:
                warning = f"could not remove previous final-assembly backup: {exc}"
        self.backup_dir = None
        self.cleanup_staging()
        return warning

    def cleanup_staging(self) -> None:
        if self.staging_dir.exists():
            _remove_transaction_path(
                self.staging_dir,
                self.project_root,
                f".{FINAL_ASSEMBLY_DIRNAME}.staging-",
            )


def _status(planned: int, succeeded: int) -> str:
    if succeeded == planned:
        return "complete"
    return "partial" if succeeded else "failed"


def _selection_fingerprint(section: Mapping[str, Any]) -> str:
    unsigned = dict(section)
    unsigned.pop("selection_fingerprint", None)
    return stable_json_hash(unsigned)


def execute_final_assembly(
    config: Any,
    *,
    reporter: ReportGenerator | None = None,
    design_executor: DesignExecutor | None = None,
) -> dict[str, Any]:
    """Regenerate and commit the complete in-silico final-assembly bundle."""

    context = load_execution_context(config)
    project_root = context.project_output_path.resolve()
    project_root.mkdir(parents=True, exist_ok=True)
    target_dir = project_root / FINAL_ASSEMBLY_DIRNAME
    staging_dir = Path(
        tempfile.mkdtemp(
            prefix=f".{FINAL_ASSEMBLY_DIRNAME}.staging-",
            dir=project_root,
        )
    ).resolve()
    transaction = FinalAssemblyOutputTransaction(
        project_root=project_root,
        target_dir=target_dir,
        staging_dir=staging_dir,
    )
    executor = design_executor or _execute_one_design
    generated_at = _utc_now()
    constructs: list[dict[str, Any]] = []
    failures: list[dict[str, Any]] = []
    try:
        for plan in context.plans:
            design_id = int(plan["parts_design_id"])
            try:
                constructs.append(executor(context, plan, staging_dir))
            except Exception as exc:
                failures.append(
                    {
                        "parts_design_id": design_id,
                        "assembly_method": plan.get("assembly_method"),
                        "error_type": type(exc).__name__,
                        "message": str(exc),
                    }
                )
        constructs.sort(key=lambda item: int(item["parts_design_id"]))
        failures.sort(key=lambda item: int(item["parts_design_id"]))
        status = _status(len(context.plans), len(constructs))
        method_counts = dict(
            sorted(Counter(str(item["assembly_method"]) for item in constructs).items())
        )
        warnings_list = [THEORETICAL_ASSEMBLY_WARNING]
        if failures:
            warnings_list.append(
                f"{len(failures)} design(s) failed during in-silico execution; "
                "successful designs were retained."
            )
        bundle_fingerprint = stable_json_hash(
            {
                "algorithm_version": FINAL_ASSEMBLY_EXECUTION_ALGORITHM_VERSION,
                "plan_selection_fingerprint": context.plan_selection_fingerprint,
                "status": status,
                "constructs": constructs,
                "failures": failures,
            }
        )
        run_summary: dict[str, Any] = {
            "schema_version": "final_assembly_run_summary.v1",
            "algorithm_version": FINAL_ASSEMBLY_EXECUTION_ALGORITHM_VERSION,
            "status": status,
            "generated_at": generated_at,
            "target_compound_id": context.source_context.target_compound_id,
            "planned_design_count": len(context.plans),
            "succeeded_count": len(constructs),
            "failed_count": len(failures),
            "method_counts": method_counts,
            "source": {
                "manifest_revision": context.manifest_revision,
                "final_assembly_plan_selection_fingerprint": (
                    context.plan_selection_fingerprint
                ),
            },
            "constructs": constructs,
            "failures": failures,
            "bundle_fingerprint": bundle_fingerprint,
            "warnings": warnings_list,
        }
        final_section: dict[str, Any] = {
            "schema_version": FINAL_ASSEMBLY_SCHEMA_VERSION,
            "status": status,
            "result_kind": "in_silico_theoretical_assembly",
            "source": "execute_final_assembly",
            "algorithm_version": FINAL_ASSEMBLY_EXECUTION_ALGORITHM_VERSION,
            "generated_at": generated_at,
            "target_compound_id": context.source_context.target_compound_id,
            "planned_design_count": len(context.plans),
            "succeeded_count": len(constructs),
            "failed_count": len(failures),
            "method_counts": method_counts,
            "output_dir": FINAL_ASSEMBLY_DIRNAME,
            "source_fingerprints": {
                "final_assembly_plan_selection_fingerprint": (
                    context.plan_selection_fingerprint
                ),
                "plasmid_selection_fingerprint": (
                    context.source_context.plasmid_selection_fingerprint
                ),
                "assembled_constructs_fingerprint": (
                    context.source_context.assembled_constructs_fingerprint
                ),
            },
            "constructs": constructs,
            "failures": failures,
            "bundle_fingerprint": bundle_fingerprint,
            "warnings": warnings_list,
        }

        report_payload = build_report_payload(
            manifest=context.manifest,
            final_assembly=final_section,
            run_summary=run_summary,
        )
        markdown, generated_by, report_warning = generate_final_design_report(
            report_payload,
            reporter=reporter,
        )
        report_path = staging_dir / FINAL_DESIGN_REPORT_FILENAME
        report_path.write_text(markdown, encoding="utf-8", newline="\n")
        if report_warning:
            warnings_list.append(report_warning)
        run_summary["warnings"] = warnings_list
        summary_path = staging_dir / RUN_SUMMARY_FILENAME
        write_json_atomic(summary_path, run_summary)
        run_summary_file = _staged_file_info(
            summary_path,
            RUN_SUMMARY_FILENAME,
            file_format="json",
        )
        report_file = _staged_file_info(
            report_path,
            FINAL_DESIGN_REPORT_FILENAME,
            file_format="markdown",
        )
        final_section["run_summary_file"] = run_summary_file
        final_section["report_file"] = report_file
        final_section["warnings"] = warnings_list
        final_section["selection_fingerprint"] = _selection_fingerprint(
            final_section
        )
        design_report_section: dict[str, Any] = {
            "schema_version": FINAL_DESIGN_REPORT_SCHEMA_VERSION,
            "status": "exported",
            "source": "execute_final_assembly",
            "generated_at": generated_at,
            "language": "zh-CN",
            "generated_by": generated_by,
            "report_file": report_file,
            "source_final_assembly_bundle_fingerprint": bundle_fingerprint,
            "warnings": [report_warning] if report_warning else [],
        }
        design_report_section["selection_fingerprint"] = _selection_fingerprint(
            design_report_section
        )

        transaction.install()
        try:
            updated = update_design_manifest(
                context.manifest_path,
                target_compound_id=context.source_context.target_compound_id,
                sections={
                    "final_assembly": final_section,
                    "final_design_report": design_report_section,
                },
                expected_revision=context.manifest_revision,
            )
        except Exception:
            transaction.rollback()
            raise
        cleanup_warning = transaction.finalize()
        if cleanup_warning:
            warnings_list.append(cleanup_warning)
        return {
            "ok": status == "complete",
            "status": status,
            "target_compound_id": context.source_context.target_compound_id,
            "counts": {
                "planned": len(context.plans),
                "succeeded": len(constructs),
                "failed": len(failures),
            },
            "method_counts": method_counts,
            "output_dir": str(target_dir.resolve()),
            "run_summary": str((target_dir / RUN_SUMMARY_FILENAME).resolve()),
            "report": str((target_dir / FINAL_DESIGN_REPORT_FILENAME).resolve()),
            "report_generated_by": generated_by,
            "failures": failures,
            "warnings": warnings_list,
            "manifest_path": str(context.manifest_path.resolve()),
            "manifest_revision": updated["revision"],
        }
    except Exception:
        transaction.cleanup_staging()
        raise


def run_final_assembly_execute(config: Any) -> dict[str, Any]:
    result = execute_final_assembly(config)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return result


__all__ = [
    "DesignExecutor",
    "FinalAssemblyOutputTransaction",
    "execute_final_assembly",
    "run_final_assembly_execute",
]
