import pytest

from app.runner import status_for_return_code, wrapper_environment


@pytest.mark.parametrize(
    "return_code,status",
    [
        (0, "succeeded"),
        (10, "source_in_sink"),
        (11, "no_solution"),
        (1, "failed"),
        (2, "failed"),
        (3, "failed"),
        (4, "failed"),
        (99, "failed"),
    ],
)
def test_return_code_mapping(return_code, status):
    assert status_for_return_code(return_code) == status


def test_wrapper_environment_removes_conda_prefix(monkeypatch):
    monkeypatch.setenv("CONDA_PREFIX", "/opt/conda")
    monkeypatch.setenv("LD_LIBRARY_PATH", "/compat:/opt/conda/lib")
    environment = wrapper_environment()
    assert "CONDA_PREFIX" not in environment
    assert environment["LD_LIBRARY_PATH"] == "/compat:/opt/conda/lib"
