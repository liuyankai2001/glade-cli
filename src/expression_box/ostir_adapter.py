"""Small fail-closed adapter around OSTIR's Python API."""

from __future__ import annotations

import hashlib
import warnings
from importlib.metadata import version
from typing import Any, Callable

from src.expression_box.parts_models import RbsPrediction


_MISSING_CLI_WARNING = "RBS Calculator Vienna is missing dependency ViennaRNA!"


def backend_versions() -> dict[str, str]:
    try:
        import RNA
    except ImportError as exc:
        raise RuntimeError(
            "ViennaRNA Python bindings are unavailable; install viennarna==2.7.2"
        ) from exc
    return {
        "ostir": version("ostir"),
        "viennarna": str(getattr(RNA, "__version__", "unknown")),
    }


def _run_ostir() -> Callable[..., list[dict[str, Any]]]:
    try:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message=_MISSING_CLI_WARNING)
            from ostir import run_ostir
    except ImportError as exc:
        raise RuntimeError("OSTIR is unavailable; install ostir==1.1.3") from exc
    backend_versions()
    return run_ostir


def predict_rbs_context(
    *,
    sequence: str,
    intended_start_position: int,
    accession: str,
    part_id: str,
) -> RbsPrediction:
    if intended_start_position < 1:
        raise ValueError("intended OSTIR start position must be positive")
    runner = _run_ostir()
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message=_MISSING_CLI_WARNING)
        results = runner(
            sequence,
            name=f"{accession}:{part_id}",
            threads=1,
            decimal_places=8,
            verbosity=0,
        )
    intended = [
        item
        for item in results
        if int(item.get("start_position") or 0) == intended_start_position
    ]
    if len(intended) != 1:
        raise ValueError(
            f"OSTIR did not return exactly one intended start for {accession} "
            f"with {part_id}"
        )
    item = intended[0]
    expression = float(item.get("expression") or 0.0)
    if expression <= 0:
        raise ValueError(
            f"OSTIR returned a non-positive expression value for {accession} "
            f"with {part_id}"
        )
    return RbsPrediction(
        part_id=part_id,
        accession=accession,
        expression=expression,
        d_g_total=float(item.get("dG_total")),
        intended_start_position=intended_start_position,
        unintended_start_count=sum(
            int(result.get("start_position") or 0) != intended_start_position
            for result in results
        ),
        context_sha256=hashlib.sha256(sequence.encode("utf-8")).hexdigest(),
    )


__all__ = ["backend_versions", "predict_rbs_context"]
