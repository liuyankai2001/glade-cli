"""Invalidate current calculation references without deleting cached source data."""

from pathlib import Path, PurePosixPath


ENZYME_RESULT_PATTERN = "enzyme_selection/candidate_*/retropath_enzyme_selection.json"


def invalidate_result_references(output_dir: Path, *patterns: str) -> None:
    root = Path(output_dir).expanduser().resolve()
    for pattern in patterns:
        relative = PurePosixPath(pattern)
        if relative.is_absolute() or ".." in relative.parts or "\\" in pattern:
            raise ValueError(f"unsafe result-reference pattern: {pattern}")
        for reference in root.glob(pattern):
            # A parent directory symlink must not redirect an unlink elsewhere.
            if not reference.parent.resolve().is_relative_to(root):
                raise ValueError(f"unsafe previous result reference: {reference}")
            if reference.is_symlink() or reference.is_file():
                reference.unlink()
