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


def test_wrapper_environment_drops_stale_scheduler_override(monkeypatch):
    monkeypatch.setenv("GLADE_SCHEDULER_OUTPUT_DIR", "/old/job/checkpoint")
    assert "GLADE_SCHEDULER_OUTPUT_DIR" not in wrapper_environment()


@pytest.mark.parametrize("mode,expected,code", [
    ("success", "succeeded", None),
    ("no_solution", "no_solution", None),
    ("timeout", "timed_out", "wall_timeout"),
    ("memory", "failed", "resource_exhausted"),
])
def test_runner_uses_upstream_status_and_resource_monitoring(tmp_path, monkeypatch, mode, expected, code):
    import json
    from dataclasses import replace
    from types import SimpleNamespace
    from app.config import Settings
    import app.runner as runner_module

    inputs = tmp_path / "input"
    inputs.mkdir()
    for name in ("source.csv", "sink.csv"):
        (inputs / name).write_text("Name,InChI\nT,InChI=1S/H2O/h1H2\n")
    # Old scheduler artifacts must not change upstream exit-code interpretation.
    checkpoint = tmp_path / "raw" / "scheduler"
    checkpoint.mkdir(parents=True)
    (checkpoint / "checkpoint.json").write_text("invalid old checkpoint")
    captured = {}

    class Process:
        def __init__(self, command, **kwargs):
            captured.update(command=command, **kwargs)
            self.returncode = 0 if mode == "success" else 11 if mode == "no_solution" else None

        def poll(self):
            return self.returncode

    def terminate(process):
        process.returncode = -15

    monkeypatch.setattr(runner_module.subprocess, "Popen", Process)
    monkeypatch.setattr(RetroPathRunner, "_terminate_process_group", staticmethod(terminate))
    monkeypatch.setattr(RetroPathRunner, "_read_memory_sample", staticmethod(lambda *_: (1000, 1000)))
    monkeypatch.setattr(runner_module.time, "sleep", lambda _: None)
    if mode == "timeout":
        moments = iter([0.0, 2000.0, 2001.0])
        monkeypatch.setattr(runner_module.time, "monotonic", lambda: next(moments))
    monkeypatch.setenv("GLADE_SCHEDULER_OUTPUT_DIR", "/stale")
    settings = replace(Settings.from_env(), memory_limit_bytes=100)
    runtime = SimpleNamespace(wrapper_version="3.9.1", wrapper_reported_version="3.9.0",
        rules_sha256="rules", workflow_sha256="workflow")
    job = {"job_id": "test", "job_dir": str(tmp_path), "created_at": "test",
           "parameters": {"max_steps": 2, "topx": 1, "dmin": 2, "dmax": 16,
                          "mwmax_source": 1000, "msc_timeout": 10}}
    result = RetroPathRunner(settings, runtime).run(job)
    assert result.status == expected and result.failure_code == code
    assert "GLADE_SCHEDULER_OUTPUT_DIR" not in captured["env"]
    manifest = json.loads((tmp_path / "run_manifest.json").read_text())
    assert manifest["scheduler"] is None and manifest["scheduler_patch"] is None
    assert manifest["search_complete"] is (mode in {"success", "no_solution"})
    assert manifest["resource_telemetry"]["resource_exhausted"] is (mode == "memory")
