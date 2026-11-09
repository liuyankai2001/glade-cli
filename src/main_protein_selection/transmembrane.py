"""Structured membrane annotations and configurable expression-score policy."""
from __future__ import annotations

import json
from collections.abc import Mapping
from copy import deepcopy
from typing import Any

POLICY_VERSION = "transmembrane_policy.v1"
RISK_MESSAGE = "跨膜风险：具有跨膜区注释，需关注宿主表达适配。"


def unknown_annotation() -> dict[str, Any]:
    return {"status": "unknown", "regions": [], "source": None}


def validate_options(allow: bool, penalty: int, bonus: int) -> None:
    if type(allow) is not bool:
        raise ValueError("allow_transmembrane must be a boolean")
    for name, value in (("transmembrane_penalty", penalty), ("no_transmembrane_bonus", bonus)):
        if type(value) is not int or value < 0:
            raise ValueError(f"{name} must be a non-negative integer")


def config_options(config: Any) -> dict[str, Any]:
    """Missing fields preserve legacy callers; RunConfig supplies the new defaults."""
    options = {
        "allow_transmembrane": getattr(config, "candidate_protein_allow_transmembrane", False),
        "transmembrane_penalty": getattr(config, "candidate_protein_transmembrane_penalty", 35),
        "no_transmembrane_bonus": getattr(config, "candidate_protein_no_transmembrane_bonus", 10),
    }
    return options


def _position(value: Any) -> int | None:
    if not isinstance(value, Mapping) or value.get("modifier", "EXACT") != "EXACT":
        return None
    number = value.get("value")
    return number if type(number) is int and number > 0 else None


def extract_annotation(entry: Any) -> dict[str, Any]:
    if not isinstance(entry, Mapping) or not entry or not entry.get("sequence"):
        return unknown_annotation()
    features = entry.get("features", [])
    if not isinstance(features, list) or any(not isinstance(item, Mapping) for item in features):
        return unknown_annotation()
    regions = []
    for item in features:
        if str(item.get("type", "")).lower() != "transmembrane":
            continue
        location = item.get("location")
        location = location if isinstance(location, Mapping) else {}
        regions.append({"start": _position(location.get("start")), "end": _position(location.get("end"))})
    return {"status": "annotated" if regions else "not_annotated", "regions": regions, "source": "UniProt"}


def annotation(value: Any = None) -> dict[str, Any]:
    """Validate persisted JSON/CSV annotations; absence means unknown, never negative."""
    if value is None or value == "":
        return unknown_annotation()
    if isinstance(value, str):
        value = json.loads(value)
    if not isinstance(value, Mapping) or set(value) != {"status", "regions", "source"}:
        raise ValueError("invalid transmembrane annotation")
    status, regions = value["status"], value["regions"]
    if status not in ("annotated", "not_annotated", "unknown") or not isinstance(regions, list):
        raise ValueError("invalid transmembrane status or regions")
    if value["source"] not in (None, "UniProt"):
        raise ValueError("invalid transmembrane source")
    if (status == "annotated") != bool(regions) or (status != "unknown" and value["source"] != "UniProt"):
        raise ValueError("inconsistent transmembrane annotation")
    for region in regions:
        if not isinstance(region, Mapping) or set(region) != {"start", "end"}:
            raise ValueError("invalid transmembrane interval")
        for number in region.values():
            if number is not None and (type(number) is not int or number < 1):
                raise ValueError("transmembrane positions must be positive integers or null")
        if region["start"] is not None and region["end"] is not None and region["start"] > region["end"]:
            raise ValueError("transmembrane start exceeds end")
    return deepcopy(dict(value))


def display(value: Any = None) -> dict[str, Any]:
    info = annotation(value)
    message = {"annotated": RISK_MESSAGE, "not_annotated": "未检出跨膜注释，实际跨膜状态未确认。",
               "unknown": "跨膜信息未知。"}[info["status"]]
    return {"状态": info["status"], "区间": info["regions"], "提示": message}
