"""Remove only reviewed, feature-free intervals and retain an exact coordinate map."""

from copy import deepcopy


def extract_reviewed_gaps(sequence, features, spans):
    intervals = []
    for span in spans:
        start, end = span.get("start_1based"), span.get("end_1based")
        if (
            type(start) is not int
            or type(end) is not int
            or not 1 <= start <= end <= len(sequence)
            or not isinstance(span.get("evidence"), str)
            or not span["evidence"].strip()
        ):
            raise ValueError(
                "gap deletion requires valid coordinates and reviewed evidence"
            )
        intervals.append((start - 1, end))
    intervals.sort()
    if any(a[1] > b[0] for a, b in zip(intervals, intervals[1:])):
        raise ValueError("reviewed gap intervals overlap")
    for feature in features:
        for part in feature["parts"]:
            if any(
                start < part["end"] and end > part["start"] for start, end in intervals
            ):
                raise ValueError("gap deletion overlaps a source feature")
    retained, fragments, cursor, length = [], [], 0, 0
    for start, end in [*intervals, (len(sequence), len(sequence))]:
        if start > cursor:
            fragment = sequence[cursor:start]
            fragments.append(fragment)
            retained.append(
                {
                    "original_start_1based": cursor + 1,
                    "original_end_1based": start,
                    "module_start_1based": length + 1,
                    "module_end_1based": length + len(fragment),
                }
            )
            length += len(fragment)
        cursor = end
    if not fragments:
        raise ValueError("gap deletion cannot remove the entire module")
    mapped = deepcopy(features)
    for feature in mapped:
        for part in feature["parts"]:
            part["start"] -= sum(
                end - start for start, end in intervals if end <= part["start"]
            )
            part["end"] -= sum(
                end - start for start, end in intervals if end <= part["end"]
            )
    return "".join(fragments), mapped, retained
