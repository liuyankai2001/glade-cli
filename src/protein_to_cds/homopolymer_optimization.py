"""Batch homopolymer editing through the shared CDS optimization pipeline."""

import sys

from src.protein_to_cds.restriction_optimization import _optimize_batch


def optimize_cds_homopolymer(config):
    return _optimize_batch(config, "homopolymer")


def run_homopolymer_optimization(config):
    from src.info_show.cds_info import format_cds_info, get_cds_info
    try:
        result = optimize_cds_homopolymer(config)
    except (ValueError, OSError, RuntimeError) as exc:
        print(f"同聚物消除失败：{exc}", file=sys.stderr)
        raise SystemExit(2) from exc
    view = get_cds_info(config)
    print(format_cds_info(view))
    return {**view, **result}
