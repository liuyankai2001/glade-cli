"""Project-bound preview, generation transaction and registered downloads."""

from __future__ import annotations

import re
import shutil
import uuid
from pathlib import Path
from typing import Any

from src.final_assemble_execute.common import resolve_project_file, sha256_file
from src.final_assemble_plan.common import gc_percent
from src.plasmid_design.catalog import ModuleCatalog
from src.plasmid_design.errors import DesignError
from src.plasmid_design.components import (
    LEGACY_SELECTION_FIELDS,
    legacy_components,
    normalize_components,
)
from src.plasmid_design.project import authenticate_request, read_project
from src.plasmid_design.registration import committed_files, export_design
from src.plasmid_design.sequence import (
    DEFAULT_COMPONENT_ORDER,
    build_design,
    feature_payloads,
    normalize_component_order,
    resolve_enzymes,
)
from src.plasmid_selection.get_plasmid_context import stable_json_hash
from src.write_manifest.store import (
    manifest_update_lock,
    read_design_manifest,
    update_design_manifest,
)


class DesignService:
    def __init__(self, config: Any):
        self.config = config

    def modules(self) -> dict:
        catalog = ModuleCatalog(self.config.data_dir)
        return {
            "resistance": [m.to_payload() for m in catalog.resistance.values()],
            "replication": [m.to_payload() for m in catalog.replication.values()],
            "terminator": [m.to_payload() for m in catalog.terminator.values()],
        }

    def context(self) -> dict:
        try:
            manifest = read_design_manifest(self.config.manifest_output_path)
        except (OSError, ValueError) as exc:
            manifest = {}
            initial_error = {
                "code": "invalid_manifest",
                "message": f"项目 manifest 无法读取：{exc}",
            }
        else:
            initial_error = None
        cds = manifest.get("cds_selection", {})
        host = cds.get("host", {}) if isinstance(cds, dict) else {}
        host = host if isinstance(host, dict) else {}
        response = {
            "project": {
                "target": self.config.target_name,
                "name": self.config.target_name,
                "host": host.get("name", "Escherichia coli MG1655"),
                "manifest_revision": manifest.get("revision", 0),
                "source_fingerprint": "",
            },
            "ready": False,
            "issues": [],
            "construct": None,
            "restriction_enzymes": [],
            "selection": None,
            "result": None,
            "active_job": None,
        }
        if initial_error:
            response["issues"].append(initial_error)
            return response
        try:
            snapshot = read_project(self.config)
        except (DesignError, OSError, ValueError, TypeError) as exc:
            response["issues"].append(
                {"code": getattr(exc, "code", "project_not_ready"), "message": str(exc)}
            )
            return response
        response["project"].update(
            {
                "manifest_revision": snapshot.source.manifest_revision,
                "source_fingerprint": snapshot.fingerprint,
                "host": snapshot.source.host_name,
            }
        )
        construct = snapshot.source.constructs[0]
        response["construct"] = {
            "id": construct.design_id,
            "name": f"完整表达构建 · 方案 {construct.design_id}",
            "length_bp": construct.length_bp,
            "gc_percent": gc_percent(str(snapshot.insert_record.seq)),
            "sequence_sha256": construct.sequence_sha256,
            "features": feature_payloads(snapshot.insert_record),
        }
        try:
            enzymes = resolve_enzymes(
                snapshot.manifest.get("cds_selection", {}).get(
                    "restriction_enzymes", []
                )
            )
            response["restriction_enzymes"] = [e.to_payload() for e in enzymes]
            response["ready"] = True
        except DesignError as exc:
            response["issues"].append({"code": exc.code, "message": str(exc)})
        selection = snapshot.manifest.get("plasmid_selection", {})
        component = (
            selection.get("component_design", {}) if isinstance(selection, dict) else {}
        )
        component = component if isinstance(component, dict) else {}
        if "components" in component:
            try:
                response["selection"] = {
                    "components": normalize_components(
                        component["components"], snapshot.catalog
                    ),
                    **{key: component.get(key) for key in LEGACY_SELECTION_FIELDS},
                }
            except DesignError as exc:
                response["issues"].append({"code": exc.code, "message": str(exc)})
        elif (
            component.get("resistance_id") in snapshot.catalog.resistance
            and component.get("replication_id") in snapshot.catalog.replication
        ):
            try:
                order = normalize_component_order(
                    component.get("component_order", DEFAULT_COMPONENT_ORDER)
                )
            except DesignError as exc:
                response["issues"].append({"code": exc.code, "message": str(exc)})
            else:
                response["selection"] = {
                    "resistance_id": component["resistance_id"],
                    "replication_id": component["replication_id"],
                    "component_order": order,
                    **{
                        f"{role}_id": (
                            component.get(f"{role}_id")
                            if isinstance(component.get(f"{role}_id"), str)
                            and component.get(f"{role}_id")
                            in snapshot.catalog.terminator
                            and snapshot.catalog.terminator[
                                component[f"{role}_id"]
                            ].role
                            == role
                            else None
                        )
                        for role in ("t0", "t1")
                    },
                }
        try:
            response["result"] = self._result(snapshot)
        except DesignError as exc:
            response["issues"].append({"code": exc.code, "message": str(exc)})
        return response

    def _prepare(self, request: dict):
        if "components" in request and request.keys() & LEGACY_SELECTION_FIELDS:
            raise DesignError(
                "组件列表不能与旧的单选参数混用。", code="invalid_components"
            )
        snapshot = read_project(self.config)
        authenticate_request(snapshot, request)
        design = build_design(
            snapshot.catalog,
            snapshot.insert_record,
            snapshot.manifest.get("cds_selection", {}).get("restriction_enzymes", []),
            request.get("resistance_id"),
            request.get("replication_id"),
            snapshot.source.constructs[0].design_id,
            component_order=request.get("component_order", DEFAULT_COMPONENT_ORDER),
            t0_id=request.get("t0_id"),
            t1_id=request.get("t1_id"),
            components=(
                normalize_components(request["components"], snapshot.catalog)
                if "components" in request
                else None
            ),
        )
        design.preview.update(
            {
                "source_fingerprint": snapshot.fingerprint,
                "manifest_revision": snapshot.source.manifest_revision,
            }
        )
        return snapshot, design

    def preview(self, request: dict) -> dict:
        return self._prepare(request)[1].preview

    def generate(self, request: dict, progress=None) -> dict:
        progress = progress or (lambda stage, message: None)
        progress("checking", "正在核对项目与酶切位点")
        snapshot, design = self._prepare(request)
        if not design.preview["valid"]:
            raise DesignError(
                "设计检查未通过，请按提示调整组件或限制酶。",
                code="design_invalid",
                issues=design.preview["issues"],
            )
        project = snapshot.source.project_output_path.resolve()
        generation_id = uuid.uuid4().hex
        parent = project / "plasmid_designs"
        parent.mkdir(exist_ok=True)
        directory = parent / generation_id
        directory.mkdir()
        committed = False
        try:
            progress("exporting", "正在生成骨架、制备片段和目标质粒文件")
            sections = export_design(snapshot, design, directory, generation_id)
            candidate = dict(snapshot.manifest)
            candidate.update(sections)
            candidate["revision"] = snapshot.source.manifest_revision + 1
            snapshot.manifest = candidate
            result = self._result(snapshot)
            if result is None:
                raise DesignError(
                    "导出结果的登记信息无法核对，请重试。", code="artifact_invalid"
                )
            progress("registering", "正在核对源文件并登记 manifest")
            with manifest_update_lock(self.config.manifest_output_path):
                current = read_project(self.config)
                authenticate_request(current, request)
                updated = update_design_manifest(
                    self.config.manifest_output_path,
                    target_compound_id=snapshot.source.target_compound_id,
                    sections=sections,
                    expected_revision=snapshot.source.manifest_revision,
                )
                committed = True
            result["manifest_revision"] = updated["revision"]
            return result
        finally:
            if not committed:
                self._remove_uncommitted(directory, parent)

    @staticmethod
    def _remove_uncommitted(directory: Path, parent: Path):
        resolved, allowed = directory.resolve(), parent.resolve()
        if resolved.parent != allowed or not re.fullmatch(
            r"[0-9a-f]{32}", resolved.name
        ):
            raise ValueError("refusing to remove an unsafe output path")
        if resolved.exists():
            shutil.rmtree(resolved)

    def _result(self, snapshot) -> dict | None:
        try:
            return self._checked_result(snapshot)
        except DesignError:
            raise
        except (OSError, ValueError, TypeError, AttributeError, KeyError) as exc:
            raise DesignError(
                "已登记文件或结果信息无效，请重新生成。", code="artifact_invalid"
            ) from exc

    def _checked_result(self, snapshot) -> dict | None:
        manifest = snapshot.manifest
        selection, final = (
            manifest.get("plasmid_selection", {}),
            manifest.get("final_assembly", {}),
        )
        component = selection.get("component_design", {})
        if (
            not component
            or final.get("status") != "complete"
            or component.get("source_fingerprint") != snapshot.fingerprint
        ):
            return None
        for section in (selection, final, manifest.get("final_assembly_plan", {})):
            unsigned = dict(section)
            recorded = unsigned.pop("selection_fingerprint", None)
            if recorded != stable_json_hash(unsigned):
                raise DesignError(
                    "已登记结果的校验信息发生变化，请重新生成。",
                    code="artifact_invalid",
                )
        if (
            final.get("source_fingerprints", {}).get("plasmid_selection_fingerprint")
            != selection["selection_fingerprint"]
        ):
            return None
        constructs = final.get("constructs", [])
        if len(constructs) != 1:
            return None
        generation_id = component["generation_id"]
        if "components" in component:
            try:
                components = normalize_components(
                    component["components"], snapshot.catalog
                )
            except DesignError as exc:
                raise DesignError(
                    "已登记组件实例无效，请重新生成。", code="artifact_invalid"
                ) from exc
            order = component.get("component_order")
            if order != [c["instance_id"] for c in components]:
                # Old single-selection generations may retain empty T0/T1 slots.
                try:
                    legacy_order = normalize_component_order(order)
                except DesignError as exc:
                    raise DesignError(
                        "已登记组件顺序无效。", code="artifact_invalid"
                    ) from exc
                if legacy_components(component, legacy_order) != components:
                    raise DesignError("已登记组件顺序不一致。", code="artifact_invalid")
        else:
            order = normalize_component_order(
                component.get("component_order", DEFAULT_COMPONENT_ORDER)
            )
            components = normalize_components(
                legacy_components(component, order), snapshot.catalog
            )
        for role in ("t0", "t1"):
            identifier = component.get(f"{role}_id")
            if identifier is not None and (
                not isinstance(identifier, str)
                or identifier not in snapshot.catalog.terminator
                or snapshot.catalog.terminator[identifier].role != role
            ):
                raise DesignError(
                    "已登记终止子引用无效，请重新生成。", code="artifact_invalid"
                )
        if not re.fullmatch(r"[0-9a-f]{32}", generation_id):
            raise DesignError("已登记结果编号无效。", code="artifact_invalid")
        files = []
        for file_id, info in committed_files(manifest).items():
            path = resolve_project_file(
                snapshot.source.project_output_path, info.get("path")
            )
            if not path.is_file() or sha256_file(path) != info.get("file_sha256"):
                raise DesignError(
                    "已登记导出文件缺失或被修改，请重新生成。", code="artifact_invalid"
                )
            files.append(
                {
                    "id": file_id,
                    "label": info.get("label", path.name),
                    "filename": path.name,
                    "format": info.get("format", ""),
                    "url": f"/api/files/{generation_id}/{file_id}",
                }
            )
        sequence = str(constructs[0].get("sequence_sha256", ""))
        final_record = next((i for i in files if i["id"] == "final_genbank"), None)
        if final_record is None:
            raise DesignError("最终质粒 GenBank 未登记。", code="artifact_invalid")
        from Bio import SeqIO

        record = SeqIO.read(
            resolve_project_file(
                snapshot.source.project_output_path,
                committed_files(manifest)["final_genbank"]["path"],
            ),
            "genbank",
        )
        return {
            "id": generation_id,
            "resistance_id": component["resistance_id"],
            "replication_id": component["replication_id"],
            "t0_id": component.get("t0_id"),
            "t1_id": component.get("t1_id"),
            "component_order": order,
            "components": components,
            "sequence_sha256": sequence,
            "length_bp": constructs[0]["length_bp"],
            "gc_percent": gc_percent(str(record.seq)),
            "manifest_revision": manifest["revision"],
            "source_fingerprint": snapshot.fingerprint,
            "files": files,
            "warnings": final.get("warnings", []),
        }

    def download(self, generation_id: str, file_id: str) -> Path:
        snapshot = read_project(self.config)
        result = self._result(snapshot)
        if result is None or result["id"] != generation_id:
            raise DesignError(
                "该结果已过期或被新的结果替换，请刷新页面。", code="stale_result"
            )
        files = committed_files(snapshot.manifest)
        if file_id not in files:
            raise DesignError("该导出文件未登记。", code="artifact_not_found")
        return resolve_project_file(
            snapshot.source.project_output_path, files[file_id]["path"]
        )
