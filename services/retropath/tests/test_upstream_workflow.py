"""Verify the packaged runtime uses the unmodified upstream search policy."""
import hashlib
import os
from pathlib import Path
import xml.etree.ElementTree as ET
import zipfile

import retropath2_wrapper

from app.config import Settings, inspect_runtime


def test_upstream_workflow_and_native_topx():
    package = Path(retropath2_wrapper.__file__).resolve().parent
    workflow = package / "workflows" / "RetroPath2.0_r20260212.knwf"
    assert hashlib.sha256(workflow.read_bytes()).hexdigest() == (
        "fa07e1c5b66817a72e7a2ea6f743d7924422a0b26f13bbcbc32c14f832090ad1"
    )
    ns = {"k": "http://www.knime.org/2008/09/XMLConfig"}
    with zipfile.ZipFile(workflow) as archive:
        name = next(n for n in archive.namelist() if n.endswith(
            "UPDATE SOURC (#935)/Score Filter (#932)/Top k Selector (#692)/settings.xml"
        ))
        root = ET.fromstring(archive.read(name))
        variable = root.find('k:config[@key="variables"]//k:config[@key="k"]/k:entry[@key="used_variable"]', ns)
        assert variable is not None and variable.get("value") == "topx"
        mode = root.find('k:config[@key="model"]/k:entry[@key="selectionMode"]', ns)
        assert mode is not None and mode.get("value") == "rows"


def test_image_has_no_custom_scheduler():
    package = Path(retropath2_wrapper.__file__).resolve().parent
    assert not list((package / "workflows").glob("*scheduler*.knwf"))
    assert "r20261104_scheduler1" not in (package / "Args.py").read_text()
    assert "r20261109_scheduler1_topx1" not in (package / "Args.py").read_text()
    assert not Path("/opt/retropath/scheduler").exists()
    knime = Path(os.environ.get("RETROPATH_KNIME_DIR", "/opt/knime"))
    assert not list(knime.rglob("org.glade.retropath.scheduler*.jar"))


def test_health_keeps_null_scheduler_compatibility_field():
    settings = Settings.from_env()
    runtime = inspect_runtime(settings)
    health = runtime.public_dict(settings)
    assert health["ready"] is True, health["errors"]
    assert health["workflow_version"] == "r20260212"
    assert health["scheduler_patch"] is None
