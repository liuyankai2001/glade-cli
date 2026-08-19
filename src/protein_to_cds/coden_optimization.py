from pathlib import Path
import importlib.metadata
import json
import tempfile
import threading
from typing import Any

from langchain.tools import tool
from CodonTransformer.CodonUtils import ORGANISM2ID
import torch
from transformers import AutoTokenizer, BigBirdForMaskedLM
from CodonTransformer.CodonPrediction import predict_dna_sequence

from src.config.cds_config import CDS_MODEL_DIR
from src.runtime.monitor import monitor
from src.tools.common.session_paths import protein_optimized_sequence_dir, protein_sequence_dir
from src.tools.protein_to_cds_tools.sequence_constraints import (
    CdsConstraintError,
    POLICY_VERSION,
    normalize_dna,
    repair_cds,
    sha256_text,
)

# 兼容原有字段名；模型目录由 global_config/cds_config 统一管理。
CODON_TRANSFORMER_MODEL_DIR = CDS_MODEL_DIR
_CODON_TRANSFORMER_LOCK = threading.RLock()
_CODON_TRANSFORMER_CACHE: tuple[Any, torch.nn.Module] | None = None


def _package_version(name: str) -> str:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "not-installed"


def _write_text_atomic(path: Path, value: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w",
        encoding="utf-8",
        dir=path.parent,
        delete=False,
        suffix=".tmp",
    ) as handle:
        handle.write(value)
        temp_path = Path(handle.name)
    temp_path.replace(path)


def _write_json_atomic(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w",
        encoding="utf-8",
        dir=path.parent,
        delete=False,
        suffix=".tmp",
    ) as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
        handle.write("\n")
        temp_path = Path(handle.name)
    temp_path.replace(path)


def _meta_parameter_names(model: torch.nn.Module) -> list[str]:
    return [
        name
        for name, parameter in model.named_parameters()
        if getattr(parameter, "is_meta", False)
    ]


def _raise_if_meta_parameters(model: torch.nn.Module) -> None:
    meta_names = _meta_parameter_names(model)
    if meta_names:
        preview = ", ".join(meta_names[:5])
        if len(meta_names) > 5:
            preview += f", ... ({len(meta_names)} total)"
        raise RuntimeError(
            "CodonTransformer 模型加载后仍存在 meta 参数，无法执行推理: "
            f"{preview}。这通常由并发加载或 torch/transformers 版本兼容问题触发。"
        )


def _load_codon_transformer():
    global _CODON_TRANSFORMER_CACHE

    with _CODON_TRANSFORMER_LOCK:
        if _CODON_TRANSFORMER_CACHE is not None:
            return _CODON_TRANSFORMER_CACHE

        model_dir = str(CODON_TRANSFORMER_MODEL_DIR)
        tokenizer = AutoTokenizer.from_pretrained(
            model_dir,
            local_files_only=True,
        )
        model = BigBirdForMaskedLM.from_pretrained(
            model_dir,
            local_files_only=True,
            low_cpu_mem_usage=False,
            device_map=None,
        )
        _raise_if_meta_parameters(model)
        model.eval()
        _CODON_TRANSFORMER_CACHE = (tokenizer, model)
        return _CODON_TRANSFORMER_CACHE


def _model_load_error(exc: Exception) -> dict[str, Any]:
    return {
        "ok": False,
        "error": "CodonTransformerModelUnavailable",
        "message": (
            "无法从本地 models/CodonTransformer 加载 CodonTransformer 模型，"
            "请检查该目录是否包含 config.json、tokenizer.json 和 model.safetensors。"
        ),
        "model_dir": str(CODON_TRANSFORMER_MODEL_DIR),
        "detail": str(exc),
    }


def read_fasta_sequence(fasta_path: str) -> str:
    """
    读取 FASTA 文件中的氨基酸序列，返回 str 格式。

    Parameters
    ----------
    fasta_path : str
        FASTA 文件路径，例如 "./protein_sequences/P21685.fasta"

    Returns
    -------
    str
        氨基酸序列字符串
    """
    fasta_file = Path(fasta_path)
    if not fasta_file.exists():
        raise FileNotFoundError(f"FASTA 文件不存在: {fasta_path}")

    sequence_lines = []

    with fasta_file.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            # 跳过空行
            if not line:
                continue

            # 跳过 header 行
            if line.startswith(">"):
                continue

            sequence_lines.append(line)
    sequence = "".join(sequence_lines)
    if not sequence:
        raise ValueError(f"FASTA 文件中没有读取到序列: {fasta_path}")
    return sequence



@tool
def search_standard_host_name():
    """
    列出 CodonTransformer 支持的标准宿主名称和 taxonomy ID。
    
    调用时机：用户不知道 codon_optimization 的 organism 参数时。
    返回：宿主名称、taxonomy ID 和可用于优化的编号。
    限制：只读；不下载序列，不运行密码子优化。
    """

    tool_name = "search_standard_host_name"
    monitor.report_start(tool_name)
    try:
        monitor.report_end(tool_name, {"host_count": len(ORGANISM2ID)})
        return ORGANISM2ID
    except Exception as exc:
        monitor.report_error(tool_name, exc)
        raise

@tool
def codon_optimization(
    Uniprot_id: str,
    organism: int,
    additional_forbidden_motifs: list[str] | None = None,
) -> dict[str, Any]:
    """
    对指定 UniProt 蛋白序列进行宿主密码子优化。
    
    调用时机：用户确认 accession 和目标宿主后，需要生成优化 CDS。
    输入：Uniprot_id、organism，以及可选的额外禁用 DNA motif。
    返回：ok、最终 DNA 序列、初始/最终文件路径、约束报告和模型错误信息。
    限制：不写 manifest；用户确认后调用 write_selected_optimized_cds_to_manifest。
    """

    tool_name = "codon_optimization"
    monitor.report_start(tool_name, {"Uniprot_id": Uniprot_id, "organism": organism})
    try:
        monitor.report_running(tool_name, "正在读取蛋白 FASTA 序列...", progress=0.1)
        protein_path = protein_sequence_dir() / f"{Uniprot_id.strip()}.fasta"
        protein = read_fasta_sequence(protein_path)
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        monitor.report_running(tool_name, "正在加载 CodonTransformer 模型...", progress=0.3)
        try:
            tokenizer, model = _load_codon_transformer()
        except Exception as exc:
            result = _model_load_error(exc)
            monitor.report_error(tool_name, result["message"])
            return result
        monitor.report_running(tool_name, "正在执行密码子优化预测...", progress=0.65)
        with _CODON_TRANSFORMER_LOCK:
            _raise_if_meta_parameters(model)
            output = predict_dna_sequence(
                protein=protein,
                organism=organism,
                device=device,
                tokenizer=tokenizer,
                model=model,
                attention_type="original_full",
                deterministic=True
            )
        ans = normalize_dna(output.predicted_dna)
        protein_id = Path(protein_path).stem
        raw_output_path = protein_optimized_sequence_dir() / f"{protein_id}_codon_transformer_raw.txt"
        output_path = protein_optimized_sequence_dir() / f"{protein_id}_optimized.txt"
        report_path = protein_optimized_sequence_dir() / f"{protein_id}_optimization_report.json"
        _write_text_atomic(raw_output_path, ans)
        try:
            repair_result = repair_cds(
                initial_sequence=ans,
                protein_sequence=protein,
                organism_id=organism,
                accession=protein_id,
                additional_forbidden_motifs=additional_forbidden_motifs,
            )
        except CdsConstraintError as exc:
            failure_report = {
                "schema_version": "protein_to_cds.constraint_repair.v2",
                "policy_version": POLICY_VERSION,
                "status": "FAIL",
                "protein": {
                    "accession": protein_id,
                    "sequence_sha256": sha256_text(protein),
                    "length_aa": len(protein),
                },
                "host": {"organism_id": organism},
                "initial": {
                    "sequence_path": str(raw_output_path.resolve()),
                    "sequence_sha256": sha256_text(ans),
                    "length_nt": len(ans),
                },
                "error": str(exc),
            }
            _write_json_atomic(report_path, failure_report)
            result = {
                "ok": False,
                "error": "CdsConstraintOptimizationFailed",
                "message": str(exc),
                "raw_output_path": str(raw_output_path.resolve()),
                "optimization_report_path": str(report_path.resolve()),
            }
            monitor.report_error(tool_name, result["message"])
            return result

        final_sequence = repair_result.pop("final_sequence")
        repair_result["optimizer"]["version"] = _package_version("dnachisel")
        repair_result["model"] = {
            "name": "CodonTransformer",
            "version": _package_version("CodonTransformer"),
            "deterministic": True,
        }
        repair_result["initial"]["sequence_path"] = str(raw_output_path.resolve())
        repair_result["final"]["sequence_path"] = str(output_path.resolve())
        _write_text_atomic(output_path, final_sequence)
        _write_json_atomic(report_path, repair_result)
        result = {
            "ok": True,
            "sequence": final_sequence[:20] + "......",
            "output_path": str(output_path.resolve()),
            "raw_output_path": str(raw_output_path.resolve()),
            "optimization_report_path": str(report_path.resolve()),
            "initial_sequence_sha256": repair_result["initial"]["sequence_sha256"],
            "final_sequence_sha256": repair_result["final"]["sequence_sha256"],
            "constraint_policy_version": POLICY_VERSION,
            "gate_status": repair_result["status"],
            "metrics": {
                "initial": {
                    "cai": repair_result["initial"]["cai"],
                    "gc_percent": repair_result["initial"]["gc_percent"],
                    "forbidden_site_count": repair_result["initial"]["forbidden_site_count"],
                    "extreme_local_gc_window_count": repair_result["initial"]["extreme_local_gc_window_count"],
                },
                "final": {
                    "cai": repair_result["final"]["cai"],
                    "gc_percent": repair_result["final"]["gc_percent"],
                    "forbidden_site_count": repair_result["final"]["forbidden_site_count"],
                    "extreme_local_gc_window_count": repair_result["final"]["extreme_local_gc_window_count"],
                },
                "changes": repair_result["changes"],
            },
            "message": f"优化成功，约束修复后的序列已保存在：{output_path}",
        }
        monitor.report_end(
            tool_name,
            {
                "output_path": str(output_path.resolve()),
                "optimization_report_path": str(report_path.resolve()),
                "gate_status": repair_result["status"],
            },
        )
        return result
    except Exception as exc:
        monitor.report_error(tool_name, exc)
        raise


if __name__ == "__main__":
    organism = 51
    print(codon_optimization.invoke({"Uniprot_id": "P21685", "organism": organism}))
