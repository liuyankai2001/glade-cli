import pytest

from app.runner import RetroPathRunner, status_for_return_code, wrapper_environment


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


def test_cgroup_memory_reader_is_fail_closed_for_invalid_samples(tmp_path):
    memory_path = tmp_path / "memory.current"
    memory_path.write_text("6442450945\n", encoding="ascii")
    assert RetroPathRunner._read_memory_bytes(memory_path) == 6442450945

    memory_path.write_text("not-a-number\n", encoding="ascii")
    assert RetroPathRunner._read_memory_bytes(memory_path) is None
    assert RetroPathRunner._read_memory_bytes(tmp_path / "missing") is None


def test_cgroup_working_set_subtracts_reclaimable_inactive_file(tmp_path):
    current_path = tmp_path / "memory.current"
    stat_path = tmp_path / "memory.stat"
    current_path.write_text("7000\n", encoding="ascii")
    stat_path.write_text("anon 4000\ninactive_file 2500\n", encoding="ascii")

    assert RetroPathRunner._read_memory_sample(current_path, stat_path) == (
        7000,
        4500,
    )

    stat_path.write_text("anon 4000\n", encoding="ascii")
    assert RetroPathRunner._read_memory_sample(current_path, stat_path) == (
        7000,
        None,
    )
