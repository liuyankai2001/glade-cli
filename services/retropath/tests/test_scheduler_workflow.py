"""Build-time workflow contract; chemistry nodes must remain the original nodes."""
from pathlib import Path
import importlib.util
import zipfile
import xml.etree.ElementTree as ET


def load_patch():
    path=Path(__file__).resolve().parents[1]/'scheduler'/'patch_workflow.py'
    if not path.exists():
        path=Path('/opt/retropath/scheduler/patch_workflow.py')
    spec=importlib.util.spec_from_file_location('scheduler_workflow_patch',path)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_patch_registers_distinct_workflow_version():
    patch=load_patch()
    assert patch.VERSION=='r20261104_scheduler1'
    assert patch.ORIGINAL_SHA=='fa07e1c5b66817a72e7a2ea6f743d7924422a0b26f13bbcbc32c14f832090ad1'


def test_patch_changes_only_control_and_keeps_original_chemistry(tmp_path):
    patch=load_patch()
    original=Path(__file__).resolve().parents[3]/'tests/C07394_retropath_root_cause_20261104/deployed_workflow.knwf'
    if not original.exists():
        import retropath2_wrapper
        original=Path(retropath2_wrapper.__file__).parent/'workflows/RetroPath2.0_r20260212.knwf'
    output=tmp_path/'patched.knwf'
    record=patch.patch(original,output)
    assert record['workflow_version']==patch.VERSION
    with zipfile.ZipFile(original) as old,zipfile.ZipFile(output) as new:
        for name in old.namelist():
            if '/Reaction Enu (#579)/' in name or '/PARSE (#914)/' in name or '/BUILD (#936)/' in name:
                assert old.read(name)==new.read(name),name
        assert any('Scheduler Select' in name for name in new.namelist())
        namespace={'k':'http://www.knime.org/2008/09/XMLConfig'}
        update=ET.fromstring(new.read(next(n for n in new.namelist() if n.endswith('Core (#1081)/UPDATE SOURC (#935)/workflow.knime'))))
        ids=[e.attrib['value'] for e in update.findall('k:config[@key="nodes"]/k:config/k:entry[@key="id"]',namespace)]
        assert '932' not in ids, 'Bypassed score-filter component must be removed, otherwise KNIME executes its unconnected mandatory inputs'
    assert record['original_sha256']==patch.ORIGINAL_SHA
