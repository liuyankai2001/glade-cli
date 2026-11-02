"""CodonTransformer inference and shared CDS artifact helpers."""

from __future__ import annotations

import hashlib
import importlib.metadata
import json
import os
import tempfile
import threading
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from src.protein_to_cds.config import (
    CODON_TRANSFORMER_MODEL_DIR,
    HostProfile,
)
from src.protein_to_cds.search_protein_sequence import ProteinSequenceRecord
from src.protein_to_cds.sequence_constraints import normalize_dna

_MODEL_LOCK = threading.RLock()
_MODEL_CACHE: dict[str, tuple[Any, Any]] = {}


class CdsOptimizationError(RuntimeError):
    """Raised when CDS generation cannot produce a valid sequence."""


@dataclass(frozen=True, slots=True)
class CdsOptimizationResult:
    accession: str
    protein_sequence_sha256: str
    input_fingerprint: str
    host_name: str
    codon_transformer_organism_id: int
    raw_sequence: str
    final_sequence: str
    raw_sequence_sha256: str
    final_sequence_sha256: str
    raw_fasta_path: Path
    optimized_fasta_path: Path
    report_path: Path
    report: dict[str, Any]
    reused_existing: bool


def _package_version(name: str) -> str:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "not-installed"


def _stable_fingerprint(payload: dict[str, Any]) -> str:
    canonical = json.dumps(
        payload,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    )
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def normalize_forbidden_motifs(values: Iterable[str]) -> tuple[str, ...]:
    """Validate and deterministically normalize user-supplied DNA motifs."""

    motifs: set[str] = set()
    for index, value in enumerate(values, start=1):
        motif = normalize_dna(value)
        if len(motif) < 4 or not set(motif).issubset(set("ACGT")):
            raise ValueError(
                f"forbidden motif {index} must contain at least four A/C/G/T bases"
            )
        motifs.add(motif)
    return tuple(sorted(motifs))


def _model_identity() -> dict[str, Any]:
    model_path = CODON_TRANSFORMER_MODEL_DIR / "model.safetensors"
    config_path = CODON_TRANSFORMER_MODEL_DIR / "config.json"
    if not model_path.is_file() or not config_path.is_file():
        raise CdsOptimizationError(
            f"local CodonTransformer model is incomplete: {CODON_TRANSFORMER_MODEL_DIR}"
        )
    stat = model_path.stat()
    return {
        "directory": str(CODON_TRANSFORMER_MODEL_DIR.resolve()),
        "model_size": stat.st_size,
        "model_mtime_ns": stat.st_mtime_ns,
        "config_sha256": hashlib.sha256(config_path.read_bytes()).hexdigest(),
        "package_version": _package_version("CodonTransformer"),
    }


def _resolve_device(requested: str) -> tuple[Any, str]:
    try:
        import torch
    except Exception as exc:
        raise CdsOptimizationError(f"PyTorch is unavailable: {exc}") from exc

    normalized = str(requested or "auto").strip().lower()
    if normalized not in {"auto", "cpu", "cuda"}:
        raise ValueError("device must be one of: auto, cpu, cuda")
    if normalized == "cuda" and not torch.cuda.is_available():
        raise CdsOptimizationError("CUDA was requested but is not available")
    selected = "cuda" if normalized == "cuda" else "cpu"
    if normalized == "auto" and torch.cuda.is_available():
        selected = "cuda"
    return torch.device(selected), selected


def _load_model(device_name: str) -> tuple[Any, Any]:
    with _MODEL_LOCK:
        cached = _MODEL_CACHE.get(device_name)
        if cached is not None:
            return cached
        try:
            import torch
            from transformers import AutoTokenizer, BigBirdForMaskedLM
        except Exception as exc:
            raise CdsOptimizationError(
                "CodonTransformer model dependencies could not be imported. "
                "Run `uv sync` and `uv pip check`; incompatible stray "
                f"torchvision/torchaudio packages are a common cause. Detail: {exc}"
            ) from exc

        try:
            tokenizer = AutoTokenizer.from_pretrained(
                str(CODON_TRANSFORMER_MODEL_DIR),
                local_files_only=True,
            )
            model = BigBirdForMaskedLM.from_pretrained(
                str(CODON_TRANSFORMER_MODEL_DIR),
                local_files_only=True,
                low_cpu_mem_usage=False,
                device_map=None,
            )
            meta_names = [
                name
                for name, parameter in model.named_parameters()
                if getattr(parameter, "is_meta", False)
            ]
            if meta_names:
                raise RuntimeError(
                    "model contains unresolved meta parameters: "
                    + ", ".join(meta_names[:5])
                )
            model.to(torch.device(device_name))
            model.eval()
        except Exception as exc:
            raise CdsOptimizationError(
                "failed to load the local CodonTransformer model from "
                f"{CODON_TRANSFORMER_MODEL_DIR}: {exc}"
            ) from exc
        _MODEL_CACHE[device_name] = (tokenizer, model)
        return tokenizer, model


def _write_atomic(path: Path, value: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            "w",
            encoding="utf-8",
            newline="\n",
            dir=path.parent,
            prefix=f".{path.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            handle.write(value)
            handle.flush()
            os.fsync(handle.fileno())
            temporary = Path(handle.name)
        os.replace(temporary, path)
        temporary = None
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def _write_json_atomic(path: Path, payload: dict[str, Any]) -> None:
    _write_atomic(
        path,
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
    )


def _dna_fasta(accession: str, label: str, sequence: str) -> str:
    return f">{accession} {label}\n{sequence}\n"


def _read_dna_fasta(path: Path) -> str:
    if not path.is_file():
        raise FileNotFoundError(path)
    sequence = "".join(
        line.strip()
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.startswith(">")
    )
    sequence = normalize_dna(sequence)
    if not sequence or not set(sequence).issubset(set("ACGT")):
        raise ValueError(f"invalid DNA FASTA: {path}")
    return sequence


def predict_cds_sequence(
    protein: ProteinSequenceRecord,
    host: HostProfile,
    device: str = "auto",
) -> tuple[str, str]:
    """Run CodonTransformer inference without any constraint repair."""
    import torch
    from CodonTransformer.CodonPrediction import predict_dna_sequence

    torch_device, device_name = _resolve_device(device)
    tokenizer, model = _load_model(device_name)
    with _MODEL_LOCK, torch.inference_mode():
        output = predict_dna_sequence(
            protein=protein.sequence,
            organism=host.codon_transformer_organism_id,
            device=torch_device,
            tokenizer=tokenizer,
            model=model,
            attention_type="original_full",
            deterministic=True,
            match_protein=True,
        )
    return normalize_dna(output.predicted_dna), device_name


def clear_model_cache() -> None:
    """Clear the in-process model cache, primarily for tests."""

    with _MODEL_LOCK:
        _MODEL_CACHE.clear()


__all__ = [
    "CdsOptimizationError",
    "CdsOptimizationResult",
    "clear_model_cache",
    "normalize_forbidden_motifs",
    "predict_cds_sequence",
]
