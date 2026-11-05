"""Build the control-node plugin and a separately versioned formal workflow."""
from __future__ import annotations
import ast
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile
import zipfile
import retropath2_wrapper
from patch_workflow import patch,VERSION

BASE=Path(__file__).resolve().parent
KNIME=Path('/opt/knime/knime_4.7.0')
JAVA=KNIME/'plugins/org.knime.binary.jre.linux.x86_64_17.0.5.20221116/jre/bin/java'
COMPILER=KNIME/'plugins/org.eclipse.jdt.core_3.30.0.v20220524-1333.jar'
PACKAGE=Path(retropath2_wrapper.__file__).resolve().parent
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def main():
    jars=sorted(set(KNIME.glob('plugins/*.jar'))|set(KNIME.glob('p2/plugins/*.jar'))|set(KNIME.glob('plugins/org.knime.core_*/*.jar')))
    classpath=':'.join(map(str,jars))
    environment=os.environ.copy();environment['_JAVA_OPTIONS']='-Xmx512m';environment['JAVA_TOOL_OPTIONS']='-XX:ActiveProcessorCount=2'
    with tempfile.TemporaryDirectory(prefix='glade-scheduler-') as temporary:
        classes=Path(temporary)
        command=[str(JAVA),'-cp',str(COMPILER),'org.eclipse.jdt.internal.compiler.batch.Main','-proc:none','-17',
                 '-cp',classpath,'-d',str(classes),str(BASE/'SchedulerState.java'),str(BASE/'SchedulerNodes.java'),str(BASE/'tests/SchedulerPolicyTest.java')]
        result=subprocess.run(command,env=environment,capture_output=True,text=True,timeout=120)
        (BASE/'compile.log').write_text(result.stdout+result.stderr,encoding='utf-8')
        if result.returncode:raise RuntimeError((BASE/'compile.log').read_text())
        result=subprocess.run([str(JAVA),'-ea','-cp',str(classes),'org.glade.retropath.scheduler.SchedulerPolicyTest'],env=environment,capture_output=True,text=True,timeout=30)
        (BASE/'policy_tests.log').write_text(result.stdout+result.stderr,encoding='utf-8')
        if result.returncode:raise RuntimeError((BASE/'policy_tests.log').read_text())
        plugin=KNIME/'plugins/org.glade.retropath.scheduler_1.0.0.jar'
        manifest='Manifest-Version: 1.0\r\nBundle-ManifestVersion: 2\r\nBundle-SymbolicName: org.glade.retropath.scheduler;singleton:=true\r\nBundle-Version: 1.0.0\r\nBundle-Name: GLADE RetroPath Scheduler\r\nBundle-Vendor: GLADE\r\nBundle-RequiredExecutionEnvironment: JavaSE-17\r\nRequire-Bundle: org.knime.core,com.google.gson\r\nExport-Package: org.glade.retropath.scheduler\r\n\r\n'
        factories=[('InitFactory','Scheduler Initialize',2,1),('SelectFactory','Scheduler Select',2,2),('AdvanceFactory','Scheduler Advance',4,2)]
        plugin_xml='<?xml version="1.0"?><plugin><extension point="org.knime.workbench.repository.nodes">'+''.join(
            f'<node id="glade_{name}" category-path="/meta" factory-class="org.glade.retropath.scheduler.SchedulerNodes${name}"/>' for name,_,_,_ in factories)+'</extension></plugin>'
        payload={'META-INF/MANIFEST.MF':manifest.encode(),'plugin.xml':plugin_xml.encode()}
        for file in classes.rglob('*.class'):
            if 'SchedulerPolicyTest' not in file.name:payload[file.relative_to(classes).as_posix()]=file.read_bytes()
        for factory,title,inputs,outputs in factories:
            description='<?xml version="1.0" encoding="UTF-8"?><knimeNode xmlns="http://knime.org/node/v2.10" type="Manipulator"><name>'+title+'</name><shortDescription>Persistent per-intermediate RetroPath scheduling.</shortDescription><fullDescription><intro><p>Control only; uses original chemistry subflows.</p></intro></fullDescription><ports>'
            description+=''.join(f'<inPort index="{i}" name="Input {i}">Original workflow table.</inPort>' for i in range(inputs))
            description+=''.join(f'<outPort index="{i}" name="Output {i}">Scheduled workflow table.</outPort>' for i in range(outputs))+'</ports></knimeNode>'
            for name in (factory,'SchedulerNodes$'+factory):payload['org/glade/retropath/scheduler/'+name+'.xml']=description.encode()
        with zipfile.ZipFile(plugin,'w',zipfile.ZIP_DEFLATED) as output:
            for name,data in sorted(payload.items()):
                info=zipfile.ZipInfo(name,date_time=(1980,1,1,0,0,0));info.compress_type=zipfile.ZIP_DEFLATED;info.external_attr=0o100644<<16;output.writestr(info,data)
    bundles=KNIME/'configuration/org.eclipse.equinox.simpleconfigurator/bundles.info'
    text=bundles.read_text();line='org.glade.retropath.scheduler,1.0.0,plugins/org.glade.retropath.scheduler_1.0.0.jar,4,false\n'
    if 'org.glade.retropath.scheduler,' in text:raise RuntimeError('Scheduler bundle already registered')
    bundles.write_text(text.rstrip()+'\n'+line)
    original=PACKAGE/'workflows/RetroPath2.0_r20260212.knwf'
    output=PACKAGE/f'workflows/RetroPath2.0_{VERSION}.knwf'
    record=patch(original,output)
    args=PACKAGE/'Args.py';before=args.read_text();lines=before.splitlines(keepends=True)
    choice=None
    for call in ast.walk(ast.parse(before)):
        if not isinstance(call,ast.Call) or not any(isinstance(a,ast.Constant) and a.value=='--rp2_version' for a in call.args):continue
        choice=next(k.value for k in call.keywords if k.arg=='choices')
    if not isinstance(choice,ast.List):raise RuntimeError('Pinned wrapper version choices not found')
    offset=sum(len(v) for v in lines[:choice.end_lineno-1])+choice.end_col_offset-1
    after=before[:offset]+', '+repr(VERSION)+before[offset:]
    ast.parse(after);args.write_text(after)
    record.update(plugin_sha256=sha(plugin),plugin_path=str(plugin),
        wrapper_args_original_sha256=hashlib.sha256(before.encode()).hexdigest(),wrapper_args_patched_sha256=sha(args),
        source_sha256={p.name:sha(p) for p in BASE.glob('*') if p.is_file() and p.suffix in {'.py','.java'}},
        tests='SchedulerPolicyTest:9 scenarios PASS',
        strategy='Sequential node execution; topx bounds completed nodes per logical batch; all overflow remains pending.')
    (BASE/'patch_manifest.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(record),flush=True)
if __name__=='__main__':main()
