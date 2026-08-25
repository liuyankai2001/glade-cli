from __future__ import annotations

import csv
import io
from dataclasses import dataclass


MAX_SOURCE_BYTES = 64 * 1024
MAX_SINK_BYTES = 50 * 1024 * 1024


class InputValidationError(ValueError):
    pass


@dataclass(frozen=True)
class ValidatedCsv:
    row_count: int


def validate_compound_csv(data: bytes, *, kind: str) -> ValidatedCsv:
    if kind not in {"source", "sink"}:
        raise ValueError(f"unsupported compound CSV kind: {kind}")
    size_limit = MAX_SOURCE_BYTES if kind == "source" else MAX_SINK_BYTES
    if not data:
        raise InputValidationError(f"{kind} file is empty")
    if len(data) > size_limit:
        raise InputValidationError(
            f"{kind} file exceeds the {size_limit}-byte size limit"
        )

    try:
        text = data.decode("utf-8-sig")
    except UnicodeDecodeError as exc:
        raise InputValidationError(f"{kind} file must be UTF-8 CSV") from exc

    try:
        rows = list(csv.reader(io.StringIO(text)))
    except csv.Error as exc:
        raise InputValidationError(f"{kind} file is not valid CSV: {exc}") from exc

    if not rows:
        raise InputValidationError(f"{kind} file is empty")
    header = [value.strip().lower() for value in rows[0][:2]]
    if header != ["name", "inchi"]:
        raise InputValidationError(
            f"{kind} header must begin with Name,InChI"
        )

    compounds = [row for row in rows[1:] if any(value.strip() for value in row)]
    if kind == "source" and len(compounds) != 1:
        raise InputValidationError("source file must contain exactly one compound")
    if kind == "sink" and not compounds:
        raise InputValidationError("sink file must contain at least one compound")

    for index, row in enumerate(compounds, start=2):
        if len(row) < 2:
            raise InputValidationError(
                f"{kind} row {index} must contain Name and InChI"
            )
        name, inchi = row[0].strip(), row[1].strip()
        if not name:
            raise InputValidationError(f"{kind} row {index} has an empty Name")
        if not inchi.startswith(("InChI=1/", "InChI=1S/")):
            raise InputValidationError(
                f"{kind} row {index} does not contain an InChI v1 value"
            )

    return ValidatedCsv(row_count=len(compounds))

