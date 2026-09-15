"""Exact GC bounds and sliding-window audits for user-selected constraints."""

from __future__ import annotations

from collections.abc import Mapping
from decimal import Decimal, InvalidOperation
from typing import Any


def validate_gc_range(minimum: Any, maximum: Any) -> tuple[Decimal, Decimal]:
    try:
        low, high = Decimal(str(minimum)), Decimal(str(maximum))
    except InvalidOperation as exc:
        raise ValueError("--gc-min 和 --gc-max 必须是百分比数值") from exc
    if not low.is_finite() or not high.is_finite() or not Decimal(0) <= low <= high <= Decimal(100):
        raise ValueError("GC 范围必须满足 0 <= gc-min <= gc-max <= 100")
    return low, high


def validate_window(value: Any, length: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not 1 <= value <= length:
        raise ValueError(f"--window 必须是正整数且不超过 CDS 长度 {length} nt")
    return value


def count_bounds(range_percent: list[str], size: int) -> tuple[int, int]:
    low, high = validate_gc_range(*range_percent)
    lo_num, lo_den = low.as_integer_ratio()
    hi_num, hi_den = high.as_integer_ratio()
    lower = -(-lo_num * size // (lo_den * 100))
    upper = hi_num * size // (hi_den * 100)
    if lower > upper:
        raise ValueError(f"指定 GC 范围对 {size} nt 没有可取的整数 GC 数量")
    return lower, upper


def percent_text(value: Decimal) -> str:
    if value == 0:
        return "0"
    text = format(value, "f")
    return text.rstrip("0").rstrip(".") if "." in text else text


def _range(value: Any) -> list[str]:
    if not isinstance(value, list) or len(value) != 2:
        raise ValueError("已保存的 GC 范围格式无效")
    low, high = validate_gc_range(*value)
    return [percent_text(low), percent_text(high)]


def saved_gc_settings(
    previous: Mapping[str, Any], report: Mapping[str, Any], length: int,
) -> dict[str, Any]:
    """Read explicit saved settings, without inventing a new GC constraint."""
    saved = previous.get("gc_settings", report.get("gc_settings"))
    settings = {}
    if saved is not None:
        if not isinstance(saved, Mapping) or set(saved) - {"global", "local"}:
            raise ValueError("已保存的 GC 设置格式无效")
        for scope, value in saved.items():
            if not isinstance(value, Mapping):
                raise ValueError("已保存的 GC 设置格式无效")
            settings[scope] = {"range_percent": _range(value.get("range_percent"))}
            if scope == "local":
                settings[scope]["window_nt"] = validate_window(value.get("window_nt"), length)
    else:
        legacy = previous.get("gc_range_percent")
        if legacy is None and report.get("schema_version") == "protein_to_cds.gc_optimization.v1":
            request = report.get("request", {})
            legacy = [request.get("gc_min"), request.get("gc_max")]
        if legacy is not None:
            settings["global"] = {"range_percent": _range(legacy)}
    for name, setting in settings.items():
        count_bounds(setting["range_percent"], length if name == "global" else setting["window_nt"])
    return settings


def effective_gc_settings(
    previous: Mapping[str, Any], report: Mapping[str, Any],
    low: Decimal, high: Decimal, window: int | None, length: int,
) -> dict[str, Any]:
    """Replace one kind of constraint while keeping the other kind active."""
    settings = saved_gc_settings(previous, report, length)
    scope = "global" if window is None else "local"
    # Canonical decimal strings make 30 and 30.0 the same request.
    settings[scope] = {"range_percent": [percent_text(low), percent_text(high)]}
    if window is not None:
        settings[scope]["window_nt"] = validate_window(window, length)
    for name, setting in settings.items():
        count_bounds(setting["range_percent"], length if name == "global" else setting["window_nt"])
    return settings


def window_gc_audit(sequence: str, setting: Mapping[str, Any]) -> dict[str, Any]:
    window = validate_window(setting["window_nt"], len(sequence))
    lower, upper = count_bounds(setting["range_percent"], window)
    prefix = [0]
    for base in sequence:
        prefix.append(prefix[-1] + (base in "GC"))
    counts = [prefix[start + window] - prefix[start] for start in range(len(sequence) - window + 1)]
    violations = [
        {"start_1based": start + 1, "end_1based": start + window,
         "gc_count": count, "gc_percent": round(100 * count / window, 8)}
        for start, count in enumerate(counts) if not lower <= count <= upper
    ]
    return {
        "window_nt": window, "range_percent": list(setting["range_percent"]),
        "minimum_count": lower, "maximum_count": upper, "window_count": len(counts),
        "minimum_observed_count": min(counts), "maximum_observed_count": max(counts),
        "min_gc_percent": round(100 * min(counts) / window, 8),
        "max_gc_percent": round(100 * max(counts) / window, 8),
        "violation_count": len(violations), "violations": violations,
    }


def gc_settings_audit(sequence: str, settings: Mapping[str, Any]):
    global_audit = local_audit = None
    passed = True
    if "global" in settings:
        lower, upper = count_bounds(settings["global"]["range_percent"], len(sequence))
        count = sequence.count("G") + sequence.count("C")
        global_audit = {"minimum_count": lower, "maximum_count": upper, "observed_count": count}
        passed = lower <= count <= upper
    if "local" in settings:
        local_audit = window_gc_audit(sequence, settings["local"])
        passed = passed and local_audit["violation_count"] == 0
    return global_audit, local_audit, passed


def local_gc_summary(before: Mapping[str, Any], after: Mapping[str, Any]) -> dict[str, Any]:
    fields = ("window_count", "min_gc_percent", "max_gc_percent", "violation_count")
    return {
        "window_nt": after["window_nt"], "range_percent": after["range_percent"],
        "input": {name: before[name] for name in fields},
        "final": {name: after[name] for name in fields},
    }
