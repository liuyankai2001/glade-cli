import pytest

from app.runner import RetroPathRunner, status_for_return_code, wrapper_environment, scheduler_completion


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


def test_wrapper_environment_isolates_scheduler_checkpoints(tmp_path):
    environment = wrapper_environment(scheduler_output_dir=tmp_path / "raw" / "scheduler")
    assert environment["GLADE_SCHEDULER_OUTPUT_DIR"] == str(tmp_path / "raw" / "scheduler")


def test_scheduler_missing_or_unfinished_checkpoint_cannot_be_no_solution(tmp_path):
    import json
    assert scheduler_completion(tmp_path, "no_solution")["error"]
    (tmp_path / "scheduler").mkdir()
    checkpoint=tmp_path / "scheduler" / "checkpoint.json"
    checkpoint.write_text(json.dumps({"stage":"before_attempt","pending_count":2}),encoding="utf-8")
    assert scheduler_completion(tmp_path,"no_solution")["error"]
    checkpoint.write_text(json.dumps({"stage":"finished","pending_count":0}),encoding="utf-8")
    assert scheduler_completion(tmp_path,"no_solution")=={"search_complete":True,"error":None,"stage":"finished","pending_count":0}


def test_found_scope_with_remaining_queue_is_still_partial(tmp_path):
    import json
    (tmp_path / "scheduler").mkdir()
    (tmp_path / "scheduler" / "checkpoint.json").write_text(json.dumps({"stage":"target_closed","pending_count":7}),encoding="utf-8")
    summary=scheduler_completion(tmp_path,"succeeded")
    assert summary["error"] is None and summary["search_complete"] is False


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
