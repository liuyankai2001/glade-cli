"""Check added HTTP dependencies without solving/upgrading the validated core.

The historical Conda environment has inconsistent upstream pip metadata. Global
`pip check` is therefore not an appropriate guard for this unchanged runtime.
"""
import platform
from importlib.metadata import requires, version

from packaging.requirements import Requirement


assert platform.python_version() == "3.9.0", platform.python_version()
expected = {
    "bactabolize": "1.0.6", "cobra": "0.21.0", "biopython": "1.79",
    "fastapi": "0.115.12", "starlette": "0.46.2", "uvicorn": "0.34.3",
    "python-multipart": "0.0.20", "h11": "0.16.0",
}
for name, pinned in expected.items():
    actual = version(name)
    assert actual == pinned, (name, actual, pinned)
for name in ("fastapi", "starlette", "uvicorn", "python-multipart"):
    for text in requires(name) or []:
        requirement = Requirement(text)
        if requirement.marker and not requirement.marker.evaluate({"extra": ""}):
            continue
        actual = version(requirement.name)
        assert requirement.specifier.contains(actual), (name, text, actual)
print("Pinned core versions and HTTP dependency requirements verified")
