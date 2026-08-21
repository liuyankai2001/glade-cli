"""Configuration constants for in-silico final-assembly execution."""

from __future__ import annotations


FINAL_ASSEMBLY_SCHEMA_VERSION = "final_assembly.v2"
FINAL_DESIGN_REPORT_SCHEMA_VERSION = "final_design_report.v2"
FINAL_ASSEMBLY_EXECUTION_ALGORITHM_VERSION = "final_assembly_execute.v1.0.0"

FINAL_ASSEMBLY_DIRNAME = "final_assembly"
RUN_SUMMARY_FILENAME = "run_summary.json"
FINAL_DESIGN_REPORT_FILENAME = "final_design_report_zh.md"

THEORETICAL_ASSEMBLY_WARNING = (
    "These files are in-silico assembly designs and do not demonstrate "
    "successful wet-lab construction or biological performance."
)


__all__ = [
    "FINAL_ASSEMBLY_DIRNAME",
    "FINAL_ASSEMBLY_EXECUTION_ALGORITHM_VERSION",
    "FINAL_ASSEMBLY_SCHEMA_VERSION",
    "FINAL_DESIGN_REPORT_FILENAME",
    "FINAL_DESIGN_REPORT_SCHEMA_VERSION",
    "RUN_SUMMARY_FILENAME",
    "THEORETICAL_ASSEMBLY_WARNING",
]
