"""Ordered component instances shared by requests, assembly and saved selections."""

from __future__ import annotations

import re

from src.plasmid_design.errors import DesignError

COMPONENT_TYPES = {"resistance", "replication", "expression", "t0", "t1"}
LEGACY_SELECTION_FIELDS = {
    "resistance_id",
    "replication_id",
    "t0_id",
    "t1_id",
    "component_order",
}


def normalize_components(value, catalog=None) -> list[dict]:
    if not isinstance(value, list) or not value:
        raise DesignError("请提供有序的组件列表。", code="invalid_components")
    result, identifiers = [], set()
    counts = dict.fromkeys(COMPONENT_TYPES, 0)
    for item in value:
        if not isinstance(item, dict) or set(item) - {
            "instance_id",
            "component_type",
            "module_id",
        }:
            raise DesignError("组件实例格式无效。", code="invalid_components")
        identifier, kind = item.get("instance_id"), item.get("component_type")
        if (
            not isinstance(identifier, str)
            or not re.fullmatch(r"[A-Za-z0-9_.:-]{1,120}", identifier)
            or identifier in identifiers
            or not isinstance(kind, str)
            or kind not in COMPONENT_TYPES
        ):
            raise DesignError(
                "组件实例 ID 必须唯一，且组件类型必须有效。", code="invalid_components"
            )
        identifiers.add(identifier)
        counts[kind] += 1
        component = {"instance_id": identifier, "component_type": kind}
        module_id = item.get("module_id")
        if kind == "expression":
            if "module_id" in item:
                raise DesignError(
                    "完整表达构建由当前项目提供。", code="invalid_components"
                )
        else:
            if (
                not isinstance(module_id, str)
                or not module_id.strip()
                or len(module_id) > 120
            ):
                raise DesignError("组件需要有效的模块 ID。", code="invalid_components")
            component["module_id"] = module_id
            if catalog is not None:
                try:
                    module = (
                        catalog.get_terminator(module_id)
                        if kind in ("t0", "t1")
                        else (
                            catalog.get_resistance(module_id)
                            if kind == "resistance"
                            else catalog.get_replication(module_id)
                        )
                    )
                except (KeyError, ValueError) as exc:
                    raise DesignError(
                        "所选模块 ID 不在当前组件库中。", code="unknown_module"
                    ) from exc
                if kind in ("t0", "t1") and module.role != kind:
                    raise DesignError(
                        f"{module.name} 不能作为 {kind.upper()} 使用。",
                        code="invalid_terminator_role",
                    )
        result.append(component)
    if (
        counts["resistance"] < 1
        or counts["replication"] != 1
        or counts["expression"] != 1
    ):
        raise DesignError(
            "需要至少一个抗性组件、一个复制模块和一个完整表达构建。",
            code="invalid_components",
        )
    return result


def legacy_components(selection: dict, order) -> list[dict]:
    result = []
    for kind in order:
        if kind == "expression":
            result.append({"instance_id": kind, "component_type": kind})
        elif selection.get(f"{kind}_id") is not None:
            result.append(
                {
                    "instance_id": kind,
                    "component_type": kind,
                    "module_id": selection[f"{kind}_id"],
                }
            )
    return result
