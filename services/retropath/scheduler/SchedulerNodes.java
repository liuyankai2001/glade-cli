package org.glade.retropath.scheduler;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.*;
import java.util.*;
import org.knime.core.data.*;
import org.knime.core.data.collection.CollectionCellFactory;
import org.knime.core.data.collection.CollectionDataValue;
import org.knime.core.data.collection.SetCell;
import org.knime.core.data.def.*;
import org.knime.core.node.*;
import org.knime.core.node.workflow.WorkflowManager;
import org.knime.core.node.workflow.NodeContainer;
import org.knime.core.node.workflow.SubNodeContainer;
import org.knime.core.node.workflow.NodeContext;

/** Only controls execution. All molecules/reactions are made by the original KNIME nodes. */
public final class SchedulerNodes {
    private static final Gson JSON=new GsonBuilder().disableHtmlEscaping().create();
    private static Session session;
    private enum Mode {INIT, SELECT, ADVANCE}
    public static final class InitFactory extends BaseFactory {public InitFactory(){super(Mode.INIT);}}
    public static final class SelectFactory extends BaseFactory {public SelectFactory(){super(Mode.SELECT);}}
    public static final class AdvanceFactory extends BaseFactory {public AdvanceFactory(){super(Mode.ADVANCE);}}
    public abstract static class BaseFactory extends NodeFactory<Model> {
        private final Mode mode;
        protected BaseFactory(Mode mode){this.mode=mode;}
        @Override public Model createNodeModel(){return new Model(mode);}
        @Override public int getNrNodeViews(){return 0;}
        @Override public NodeView<Model> createNodeView(int index,Model model){throw new IllegalArgumentException();}
        @Override public boolean hasDialog(){return false;}
        @Override public NodeDialogPane createNodeDialogPane(){throw new UnsupportedOperationException();}
    }
    private static String text(DataCell cell) {
        if(cell==null||cell.isMissing())return "";
        if(cell instanceof StringValue value)return value.getStringValue();
        return cell.toString();
    }
    private static List<String> strings(DataCell cell) {
        List<String> result=new ArrayList<>();
        if(cell instanceof CollectionDataValue values){for(DataCell value:values)result.add(text(value));}
        else if(cell!=null&&!cell.isMissing())result.add(text(cell));
        return result;
    }
    private static double number(DataCell cell) {
        if(cell instanceof DoubleValue v)return v.getDoubleValue();
        if(cell instanceof IntValue v)return v.getIntValue();
        if(cell instanceof CollectionDataValue values){double v=Double.NaN;
            for(DataCell c:values){double x=number(c);v=Double.isNaN(v)?x:(session.state.minimize?Math.min(v,x):Math.max(v,x));}return v;}
        return Double.parseDouble(text(cell));
    }
    static String identity(String inchi) {
        // Match the original Lite InChI String Manipulation nodes exactly:
        // truncate at the first stereochemical, charge or protonation layer.
        return inchi.split("/[btqp]",2)[0];
    }
    private static DataCell cell(DataRow row,DataTableSpec spec,String name) {
        int index=spec.findColumnIndex(name);
        if(index<0)throw new IllegalStateException("Required column "+name+" missing from "+spec);
        return row.getCell(index);
    }
    private static String id(DataRow row,DataTableSpec spec) {
        int lite=spec.findColumnIndex("Lite InChI");
        return identity(text(lite>=0?row.getCell(lite):cell(row,spec,"InChI")));
    }
    private static DataTableSpec sourceSpec(DataTableSpec spec) {
        int compound=spec.findColumnIndex("Compound");
        if(compound<0)throw new IllegalArgumentException("Original standardized source table has no Compound column: "+spec);
        return new DataTableSpec(new DataColumnSpecCreator("Name",SetCell.getCollectionType(StringCell.TYPE)).createSpec(),
            new DataColumnSpecCreator("InChI",StringCell.TYPE).createSpec(),
            new DataColumnSpecCreator("SMILES",StringCell.TYPE).createSpec(),spec.getColumnSpec(compound));
    }
    private static DataCell[] sourceCells(DataRow row,DataTableSpec spec) {
        List<DataCell> names=new ArrayList<>();for(String name:strings(cell(row,spec,"Name")))names.add(new StringCell(name));
        return new DataCell[]{CollectionCellFactory.createSetCell(names),cell(row,spec,"InChI"),
            cell(row,spec,"SMILES"),cell(row,spec,"Compound")};
    }
    private static final class Session {
        final SchedulerState state;
        final DataTableSpec sourceSpec;
        final Map<String,DataCell[]> sources=new HashMap<>();
        long provenanceReplays=0;
        final Path directory;
        Session(SchedulerState state,DataTableSpec spec) throws IOException {
            this.state=state;this.sourceSpec=spec;
            String output=System.getenv("GLADE_SCHEDULER_OUTPUT_DIR");
            if(output==null||output.isBlank())throw new IOException("GLADE_SCHEDULER_OUTPUT_DIR is required");
            directory=Paths.get(output).toAbsolutePath();Files.createDirectories(directory);
            if("1".equals(System.getenv("GLADE_SCHEDULER_DIAGNOSTICS"))){
                WorkflowManager context=NodeContext.getContext().getWorkflowManager();
                Thread monitor=new Thread(()->{while(true){try{Thread.sleep(100);List<Map<String,String>> rows=new ArrayList<>();diagnostics(context,rows);if(rows.size()>10)atomic("workflow_diagnostics.json",rows);}catch(InterruptedException stop){return;}catch(Exception ignored){}}},"scheduler-diagnostics");
                monitor.setDaemon(true);monitor.start();
            }
        }
        void diagnostics(NodeContainer node,List<Map<String,String>> rows){
            rows.add(Map.of("id",node.getID().toString(),"name",node.getName(),"state",node.getNodeContainerState().toString(),"message",String.valueOf(node.getNodeMessage())));
            if(node instanceof WorkflowManager manager)for(NodeContainer child:manager.getNodeContainers())diagnostics(child,rows);
            else if(node instanceof SubNodeContainer component)for(NodeContainer child:component.getNodeContainers())diagnostics(child,rows);
        }
        void atomic(String name,Object value) throws IOException {
            Path temp=directory.resolve(name+".tmp");
            Files.writeString(temp,JSON.toJson(value)+"\n",StandardCharsets.UTF_8);
            Files.move(temp,directory.resolve(name),StandardCopyOption.REPLACE_EXISTING,StandardCopyOption.ATOMIC_MOVE);
        }
        void event(Object value) throws IOException {
            Files.writeString(directory.resolve("events.jsonl"),JSON.toJson(value)+"\n",StandardCharsets.UTF_8,
                StandardOpenOption.CREATE,StandardOpenOption.APPEND);
        }
        void checkpoint(String stage) throws IOException {
            Map<String,List<String>> origins=new LinkedHashMap<>();
            for(var value:sources.entrySet())origins.put(value.getKey(),strings(value.getValue()[0]));
            atomic("checkpoint.json",Map.of("stage",stage,"workflow_version","r20261104_scheduler1",
                "at",java.time.Instant.now().toString(),"state",state,"pending_count",state.pendingCount(),
                "retained_edge_count",state.edges.size(),"processing_batch",state.completedNodes/state.topx,
                "source_origins",origins,"targets_closed",state.targetsClosed(),"provenance_replays",provenanceReplays));
            atomic("pending.json",state.pending());
        }
        void mergeOrigins(String key,List<String> incoming){
            DataCell[] cells=sources.get(key);
            if(cells==null)return;
            LinkedHashSet<String> names=new LinkedHashSet<>(strings(cells[0]));
            if(!names.addAll(incoming))return;
            List<DataCell> values=new ArrayList<>();for(String name:names)values.add(new StringCell(name));
            cells[0]=CollectionCellFactory.createSetCell(values);
            var node=state.nodes.get(key);
            // A new root may reach a shared intermediate after it was expanded.
            // Replay original chemistry with merged source provenance so unchanged
            // PARSE/BUILD/SQL scope outputs gain the missing originating root.
            if(node!=null&&node.eligible&&node.depth<state.maxSteps&&node.status.equals("processed")){
                node.status="pending";node.diameter=state.dmax;provenanceReplays++;
            }
        }
        BufferedDataTable nextTable(ExecutionContext exec) throws Exception {
            BufferedDataContainer out=exec.createDataContainer(sourceSpec);
            var node=state.next();
            if(node!=null){
                DataCell[] values=sources.get(node.id);
                if(values==null)throw new IllegalStateException("Lost original source cells for "+node.id);
                out.addRowToTable(new DefaultRow("scheduler_"+node.sequence,values));
            }
            out.close();return out.getTable();
        }
    }
    public static final class Model extends NodeModel {
        private final Mode mode;
        Model(Mode mode){super(mode==Mode.ADVANCE?4:2,mode==Mode.INIT?1:2);this.mode=mode;}
        private int integer(String name,int fallback){var value=getAvailableFlowVariables().get(name);return value==null?fallback:value.getIntValue();}
        private String string(String name,String fallback){var value=getAvailableFlowVariables().get(name);return value==null?fallback:value.getStringValue();}
        private void variables(int depth,int diameter){pushFlowVariableInt("scheduler_depth",depth);pushFlowVariableInt("scheduler_diameter",diameter);pushFlowVariableInt("scheduler_attempt",session==null?0:(int)session.state.completedAttempts);}
        @Override protected DataTableSpec[] configure(DataTableSpec[] input) throws InvalidSettingsException {
            variables(0,16);
            if(mode==Mode.INIT)return new DataTableSpec[]{sourceSpec(input[0])};
            if(mode==Mode.SELECT)return new DataTableSpec[]{input[0],input[1]};
            return new DataTableSpec[]{withIteration(input[0]),sourceSpec(input[3])};
        }
        @Override protected BufferedDataTable[] execute(BufferedDataTable[] input,ExecutionContext exec) throws Exception {
            exec.checkCanceled();
            if(mode==Mode.INIT){
                Set<String> sink=new LinkedHashSet<>();var spec=input[1].getDataTableSpec();
                for(DataRow row:input[1])if(number(cell(row,spec,"Initial Sink"))==1)sink.add(id(row,spec));
                var state=new SchedulerState(integer("max-steps",10),integer("topx",100),integer("dmin",2),integer("dmax",16),
                    !string("score_mode_label","MIN").startsWith("MAX"),sink);
                session=new Session(state,sourceSpec(input[0].getDataTableSpec()));
                spec=input[0].getDataTableSpec();
                for(DataRow row:input[0]){String key=id(row,spec);state.roots.add(key);state.offer(key,0,0,true);session.sources.put(key,sourceCells(row,spec));}
                BufferedDataTable next=session.nextTable(exec);session.checkpoint("initialized");
                return new BufferedDataTable[]{next};
            }
            if(session==null)throw new IllegalStateException("Scheduler initialization did not execute");
            if(mode==Mode.SELECT){
                var node=session.state.active;
                int diameter=node==null?session.state.dmax:node.diameter;
                variables(node==null?0:node.depth,diameter);
                BufferedDataContainer rules=exec.createDataContainer(input[1].getDataTableSpec());
                int index=input[1].getDataTableSpec().findColumnIndex("Diameter");
                for(DataRow row:input[1])if((int)number(row.getCell(index))==diameter)rules.addRowToTable(row);
                rules.close();session.checkpoint("before_attempt");
                return new BufferedDataTable[]{input[0],rules.getTable()};
            }
            var current=session.state.active;
            if(current==null)return new BufferedDataTable[]{iterationTable(input[0],0,exec),session.nextTable(exec)};
            var productSpec=input[1].getDataTableSpec();var transformSpec=input[2].getDataTableSpec();var eligibleSpec=input[3].getDataTableSpec();
            Map<String,LinkedHashSet<String>> groups=new LinkedHashMap<>();
            for(DataRow row:input[1]){
                String key=id(row,productSpec);
                for(String transform:strings(cell(row,productSpec,"Transformation ID")))groups.computeIfAbsent(transform,k->new LinkedHashSet<>()).add(key);
            }
            Set<String> eligible=new HashSet<>();
            for(DataRow row:input[3]){String key=id(row,eligibleSpec);eligible.add(key);session.sources.putIfAbsent(key,sourceCells(row,eligibleSpec));}
            List<SchedulerState.Edge> edges=new ArrayList<>();Map<String,Double> scores=new LinkedHashMap<>();
            for(DataRow row:input[2]){
                double score=number(cell(row,transformSpec,"Score"));
                for(String transform:strings(cell(row,transformSpec,"Transformation ID"))){
                    Set<String> products=groups.get(transform);
                    if(products==null||products.isEmpty())continue;
                    for(String rule:strings(cell(row,transformSpec,"Rule ID")))edges.add(new SchedulerState.Edge(current.id,List.copyOf(products),rule,score,current.diameter));
                    for(String key:products)scores.merge(key,score,session.state.minimize?Math::min:Math::max);
                }
            }
            List<SchedulerState.Product> products=new ArrayList<>();
            List<String> origins=strings(session.sources.get(current.id)[0]);
            for(var entry:scores.entrySet()){
                session.mergeOrigins(entry.getKey(),origins);
                var previous=session.state.nodes.get(entry.getKey());
                products.add(new SchedulerState.Product(entry.getKey(),entry.getValue(),eligible.contains(entry.getKey())||(previous!=null&&previous.eligible)));
            }
            int oldEdges=session.state.edges.size(),depth=current.depth;
            var attempt=session.state.finish(edges,products);
            for(int i=oldEdges;i<session.state.edges.size();i++)session.event(Map.of("type","edge","edge",session.state.edges.get(i)));
            session.event(Map.of("type","attempt","attempt",attempt));
            // Closure must include every product, bounded by chemical depth, never by loop rounds.
            BufferedDataTable next=session.nextTable(exec);
            session.checkpoint(session.state.pendingCount()==0?(session.state.targetsClosed()?"target_closed":"finished"):"after_attempt");
            return new BufferedDataTable[]{iterationTable(input[0],depth,exec),next};
        }
        private static DataTableSpec withIteration(DataTableSpec original){
            List<DataColumnSpec> columns=new ArrayList<>();boolean found=false;
            for(DataColumnSpec column:original){if(column.getName().equals("Iteration")){columns.add(new DataColumnSpecCreator("Iteration",IntCell.TYPE).createSpec());found=true;}else columns.add(column);}
            if(!found)columns.add(new DataColumnSpecCreator("Iteration",IntCell.TYPE).createSpec());
            return new DataTableSpec(columns.toArray(DataColumnSpec[]::new));
        }
        private static BufferedDataTable iterationTable(BufferedDataTable original,int depth,ExecutionContext exec) throws Exception {
            var spec=withIteration(original.getDataTableSpec());int index=spec.findColumnIndex("Iteration");
            BufferedDataContainer out=exec.createDataContainer(spec);
            for(DataRow row:original){DataCell[] cells=new DataCell[spec.getNumColumns()];for(int i=0;i<row.getNumCells();i++)cells[i]=row.getCell(i);cells[index]=new IntCell(depth);out.addRowToTable(new DefaultRow(row.getKey(),cells));}
            out.close();return out.getTable();
        }
        @Override protected void reset(){}
        @Override protected void saveSettingsTo(NodeSettingsWO settings){settings.addString("scheduler_version","r20261104_scheduler1");}
        @Override protected void loadValidatedSettingsFrom(NodeSettingsRO settings){}
        @Override protected void validateSettings(NodeSettingsRO settings){}
        @Override protected void loadInternals(File directory,ExecutionMonitor monitor){}
        @Override protected void saveInternals(File directory,ExecutionMonitor monitor){}
    }
}
