"""Resolve a download's selected chassis without registering or writing hosts."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import requests

from src.config.config import GEM_DIR
from src.main_protein_selection.genome_identity import read_genome_identity
from src.main_protein_selection.host_registry import get_host_preset, get_taxonomy_snapshot, load_registry
from src.main_protein_selection.taxonomy_compatibility import (
    ChassisTaxonomyProfile,
    UNIPROT_TAXONOMY_URL,
    _profile_from_record,
    resolve_chassis_taxonomy,
)


def resolve_download_host(config: Any, session: requests.Session) -> ChassisTaxonomyProfile:
    manifest = {}
    manifest_path = getattr(config, "manifest_output_path", None)
    if manifest_path and Path(manifest_path).is_file():
        manifest = json.loads(Path(manifest_path).read_text(encoding="utf-8"))
        if not isinstance(manifest, dict):
            raise ValueError("设计清单必须是 JSON 对象")
    selection = manifest.get("main_enzyme_selection") or {}
    if not isinstance(selection, dict):
        raise ValueError("设计清单宿主信息格式错误")
    keys = list(dict.fromkeys(key for key in (
        getattr(config, "chassis_key", None), selection.get("chassis_key"),
    ) if key))
    presets = [get_host_preset(key) for key in keys]
    if len({preset["taxon_id"] for preset in presets}) > 1:
        raise ValueError("配置与设计清单的底盘宿主不一致")
    key = keys[0] if keys else None

    genome_path = getattr(config, "genome_file", None)
    if genome_path:
        path = Path(genome_path).expanduser()
        if not path.is_absolute():
            path = Path(config.root_dir) / path
        identity = read_genome_identity(path)
        if presets and presets[0]["taxon_id"] != identity.taxon_id:
            raise ValueError("配置基因组与设计清单的底盘宿主不一致")
        taxon_id, scientific_name = identity.taxon_id, identity.name
        if key is None:
            key = next((name for name, value in load_registry()["chassis_taxon_presets"].items()
                        if value["taxon_id"] == taxon_id), None)
    elif key:
        taxon_id, scientific_name = presets[0]["taxon_id"], presets[0]["name"]
    else:
        default_model = Path(getattr(config, "gem_dir", GEM_DIR)) / "iML1515.json"
        model_path = getattr(config, "model_path", None)
        if model_path is not None and Path(model_path).resolve() != default_model.resolve():
            raise ValueError("无法确定自定义模型的底盘宿主，请配置基因组或改用 --uniprot 编号下载")
        key = "ecoli_mg1655"
        preset = get_host_preset(key)
        taxon_id, scientific_name = preset["taxon_id"], preset["name"]

    if key:
        return resolve_chassis_taxonomy(key, session=session, cache_root=None)

    # Newly supplied genomes need not be inserted into the persistent registry.
    try:
        response = session.get(f"{UNIPROT_TAXONOMY_URL}/{taxon_id}", timeout=20)
        response.raise_for_status()
        profile = _profile_from_record(f"taxon_{taxon_id}", response.json(), source="uniprot_live")
        if profile.taxon_id != taxon_id:
            raise ValueError("UniProt taxonomy 返回了不同的底盘编号")
        return profile
    except (requests.RequestException, ValueError, TypeError, AttributeError) as exc:
        snapshot = get_taxonomy_snapshot(taxon_id)
        if snapshot is not None:
            return _profile_from_record(f"taxon_{taxon_id}", snapshot, source="builtin_snapshot")
        return ChassisTaxonomyProfile(
            chassis_key=f"taxon_{taxon_id}", taxon_id=taxon_id,
            scientific_name=scientific_name, rank="unknown", lineage=(),
            source="unavailable", status="unknown", message=str(exc),
        )
