"""Deterministic control-only patch of the pinned deployed KNIME workflow."""
from __future__ import annotations
import copy
import hashlib
from pathlib import Path
import zipfile
import xml.etree.ElementTree as ET

VERSION='r20261104_scheduler1'
ORIGINAL_SHA='fa07e1c5b66817a72e7a2ea6f743d7924422a0b26f13bbcbc32c14f832090ad1'
NS='http://www.knime.org/2008/09/XMLConfig'
ET.register_namespace('',NS)
ET.register_namespace('xsi','http://www.w3.org/2001/XMLSchema-instance')
def tag(name):return '{'+NS+'}'+name
def config(parent,key):return parent.find(f'{{{NS}}}config[@key="{key}"]')
def entry(parent,key):return parent.find(f'{{{NS}}}entry[@key="{key}"]')
def put(parent,key,value,kind='xstring'):
    value_entry=entry(parent,key)
    if value_entry is None:value_entry=ET.SubElement(parent,tag('entry'),key=key,type=kind)
    value_entry.set('value',str(value));value_entry.attrib.pop('isnull',None)
def xml(root):return ET.tostring(root,encoding='utf-8',xml_declaration=True)
def patch(original:Path,output:Path):
    source=original.read_bytes();actual=hashlib.sha256(source).hexdigest()
    if actual!=ORIGINAL_SHA:raise ValueError('Unsupported original workflow SHA256: '+actual)
    with zipfile.ZipFile(original) as archive:
        members={name:archive.read(name) for name in archive.namelist()}
        prefix=next(n.removesuffix('workflow.knime') for n in members if n.endswith('Core (#1081)/workflow.knime'))
        core=ET.fromstring(members[prefix+'workflow.knime']);nodes=config(core,'nodes');connections=config(core,'connections')
        def values(element):return {e.attrib.get('key'):e.attrib.get('value') for e in element if e.tag==tag('entry')}
        def remove_connections(predicate):
            for connection in list(connections):
                if predicate(values(connection)):connections.remove(connection)
        remove_connections(lambda v:v.get('sourceID') in {'780','784'} or v.get('destID') in {'780','784'})
        for node in list(nodes):
            if values(node).get('id') in {'780','784'}:nodes.remove(node)
        replaced={('781','782','2','2'),('781','914','2','2'),('781','909','0','1'),
                  ('940','782','1','4'),('935','783','3','3'),('936','783','1','1'),('936','963','1','2')}
        remove_connections(lambda v:(v.get('sourceID'),v.get('destID'),v.get('sourcePort'),v.get('destPort')) in replaced)
        def connect(source,destination,source_port,destination_port):
            item=ET.SubElement(connections,tag('config'),key='scheduler_connection_'+str(len(connections)))
            for key,value in [('sourceID',source),('destID',destination),('sourcePort',source_port),('destPort',destination_port)]:put(item,key,value,'xint')
        new_nodes=[(1200,'Scheduler Initialize','InitFactory',1),(1201,'Scheduler Select','SelectFactory',2),(1202,'Scheduler Advance','AdvanceFactory',2)]
        template=ET.fromstring(members[prefix+'Java Edit Variable (#909)/settings.xml'])
        for identifier,name,factory,outputs in new_nodes:
            path=f'{name} (#{identifier})/settings.xml'
            node=ET.SubElement(nodes,tag('config'),key='node_'+str(identifier))
            for key,value,kind in [('id',identifier,'xint'),('node_settings_file',path,'xstring'),
                                   ('node_is_meta','false','xboolean'),('node_type','NativeNode','xstring'),
                                   ('ui_classname','org.knime.core.node.workflow.NodeUIInformation','xstring')]:put(node,key,value,kind)
            settings=copy.deepcopy(template)
            for key in ['flow_stack','model','variables','ports','factory_settings']:
                old=config(settings,key)
                if old is not None:settings.remove(old)
                ET.SubElement(settings,tag('config'),key=key)
            put(settings,'factory','org.glade.retropath.scheduler.SchedulerNodes$'+factory)
            for key,value in [('node-name',name),('name',name),('node-bundle-name','GLADE RetroPath Scheduler'),
                              ('node-bundle-symbolic-name','org.glade.retropath.scheduler'),('node-bundle-version','1.0.0'),
                              ('node-bundle-vendor','GLADE')]:put(settings,key,value)
            for index in range(1,outputs+1):
                port=ET.SubElement(config(settings,'ports'),tag('config'),key='port_'+str(index));put(port,'index',index,'xint')
            members[prefix+path]=xml(settings)
        for edge in [(938,1200,3,1),(938,1200,2,2),(1081,1200,1,0),(1200,781,1,2),
                     (781,1201,2,1),(940,1201,1,2),(1201,782,1,2),(1201,782,2,4),
                     (1201,914,1,2),(1201,909,0,1),(1201,1049,0,1),
                     (936,1202,1,1),(914,1202,2,2),(914,1202,3,3),(935,1202,3,4),
                     (1202,783,1,1),(1202,783,2,3),(936,963,1,2),(783,1080,1,2)]:connect(*edge)
        members[prefix+'workflow.knime']=xml(core)
        loop_path=prefix+'Recursive Loop End _2 ports_ (#783)/settings.xml'
        loop=ET.fromstring(members[loop_path]);model=config(loop,'model')
        put(model,'CFG_MaxNrIterations',2147483647,'xint');put(model,'CFG_AddIterationColumn','false','xboolean')
        variables=config(loop,'variables')
        for old in list(variables):variables.remove(old)
        members[loop_path]=xml(loop)
        for identifier,variable in [(909,'scheduler_depth'),(1049,'scheduler_attempt')]:
            path=prefix+f'Java Edit Variable (#{identifier})/settings.xml';root=ET.fromstring(members[path])
            inputs=config(config(root,'model'),'inVars')
            for value in inputs.iter(tag('entry')):
                if value.attrib.get('key')=='Name' and value.attrib.get('value')=='currentIteration':value.set('value',variable)
            members[path]=xml(root)
        # Keep original mass filtering and score aggregation; remove only the beam cut.
        path=prefix+'UPDATE SOURC (#935)/workflow.knime';update=ET.fromstring(members[path]);wires=config(update,'connections')
        update_nodes=config(update,'nodes')
        for node in list(update_nodes):
            if values(node).get('id')=='932':update_nodes.remove(node)
        for wire in list(wires):
            v=values(wire)
            if v.get('sourceID')=='932' or v.get('destID')=='932':wires.remove(wire)
        wire=ET.SubElement(wires,tag('config'),key='scheduler_retain_all_scored_sources')
        for key,value in [('sourceID',924),('destID',942),('sourcePort',1),('destPort',1)]:put(wire,key,value,'xint')
        members[path]=xml(update)
        output.parent.mkdir(parents=True,exist_ok=True)
        with zipfile.ZipFile(output,'w',zipfile.ZIP_DEFLATED) as target:
            for name,data in sorted(members.items()):
                info=zipfile.ZipInfo(name,date_time=(1980,1,1,0,0,0));info.compress_type=zipfile.ZIP_DEFLATED;info.external_attr=0o100644<<16
                target.writestr(info,data)
    with zipfile.ZipFile(original) as old,zipfile.ZipFile(output) as new:
        changed=[n for n in old.namelist() if old.read(n)!=new.read(n)]
        assert all('/Reaction Enu (#579)/' not in n and '/PARSE (#914)/' not in n and '/BUILD (#936)/' not in n for n in changed)
    return {'workflow_version':VERSION,'original_sha256':actual,'patched_sha256':hashlib.sha256(output.read_bytes()).hexdigest(),
            'changed_original_members':changed,'added_members':sorted(set(members)-set(old.namelist())),
            'chemistry':'Original standardization, reaction enumeration, PARSE/MCS, BUILD and output subflows retained.'}
