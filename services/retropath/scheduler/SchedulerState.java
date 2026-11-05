package org.glade.retropath.scheduler;

import java.util.*;

/** Chemistry-free policy used by the formal KNIME scheduling control nodes. */
public final class SchedulerState {
    public static final class Node {
        public final String id;
        public int depth, diameter;
        public double score;
        public final long sequence;
        public boolean eligible;
        public String status;
        Node(String id,int depth,double score,long sequence,int diameter,boolean eligible) {
            this.id=id;this.depth=depth;this.score=score;this.sequence=sequence;
            this.diameter=diameter;this.eligible=eligible;this.status="pending";
        }
    }
    public record Product(String id,double score,boolean eligible) {}
    public record Edge(String substrate,List<String> products,String rule,double score,int diameter) {}
    public record Attempt(String node,int depth,int diameter,int edges,int novel_nodes,
                          boolean closes_to_sink,boolean fallback) {}
    public final int maxSteps,topx,dmin,dmax;
    public final boolean minimize;
    public final Set<String> sink;
    public final Map<String,Node> nodes=new LinkedHashMap<>();
    public final List<String> roots=new ArrayList<>();
    public transient final List<Edge> edges=new ArrayList<>();
    private transient final Set<String> signatures=new HashSet<>();
    public Node active;
    public long sequence=0,completedNodes=0,completedAttempts=0;
    public SchedulerState(int steps,int topx,int dmin,int dmax,boolean minimize,Set<String> sink) {
        if(steps<1||topx<1||dmin<0||dmax<dmin)throw new IllegalArgumentException("Invalid scheduler bounds");
        this.maxSteps=steps;this.topx=topx;this.dmin=dmin;this.dmax=dmax;
        this.minimize=minimize;this.sink=new LinkedHashSet<>(sink);
    }
    private boolean better(double score,double previous){return minimize?score<previous:score>previous;}
    public Comparator<Node> order(){return Comparator.comparingInt((Node n)->n.depth)
        .thenComparingDouble(n->minimize?n.score:-n.score).thenComparingLong(n->n.sequence);}
    public boolean offer(String id,int depth,double score,boolean eligible) {
        if(sink.contains(id))return false;
        if(!Double.isFinite(score))throw new IllegalArgumentException("Nonfinite rule score");
        Node n=nodes.get(id);boolean useful=false;
        if(n==null){n=new Node(id,depth,score,sequence++,dmax,eligible);nodes.put(id,n);useful=true;}
        else if(depth<n.depth){n.depth=depth;n.score=score;n.diameter=dmax;n.status="pending";useful=true;}
        else if(depth==n.depth&&better(score,n.score)){n.score=score;}
        n.eligible|=eligible;
        if(!n.eligible)n.status="mass_boundary";
        else if(n.depth>=maxSteps)n.status="depth_boundary";
        return useful&&n.status.equals("pending");
    }
    public int pendingCount(){return (int)nodes.values().stream()
        .filter(n->n.status.equals("pending")||n.status.equals("active")).count();}
    public List<Node> pending(){return nodes.values().stream()
        .filter(n->n.status.equals("pending")||n.status.equals("active")).sorted(order()).toList();}
    public Node next() {
        if(active!=null)return active;
        active=nodes.values().stream().filter(n->n.status.equals("pending")).min(order()).orElse(null);
        if(active!=null)active.status="active";
        return active;
    }
    public Map<String,Integer> closed() {
        Map<String,Integer> height=new HashMap<>();for(String id:sink)height.put(id,0);
        for(int depth=1;depth<=maxSteps;depth++){
            Set<String> add=new HashSet<>();
            for(Edge e:edges)if(!height.containsKey(e.substrate)&&!e.products.isEmpty()
                    &&e.products.stream().allMatch(height::containsKey))add.add(e.substrate);
            if(add.isEmpty())break;
            for(String id:add)height.put(id,depth);
        }
        return height;
    }
    public boolean targetsClosed(){return !roots.isEmpty()&&closed().keySet().containsAll(roots);}
    public Attempt finish(List<Edge> found,List<Product> products) {
        if(active==null)throw new IllegalStateException("No active node");
        Node current=active;int diameter=current.diameter,depth=current.depth,novel=0;
        for(Edge e:found){
            List<String> ids=new ArrayList<>(e.products);Collections.sort(ids);
            String signature=e.substrate+"\t"+e.rule+"\t"+String.join("\t",ids);
            if(signatures.add(signature))edges.add(e);
        }
        for(Product p:products)if(offer(p.id,depth+1,p.score,p.eligible))novel++;
        boolean complete=closed().containsKey(current.id);
        boolean fallback=novel==0&&!complete;
        if(fallback&&diameter-2>=dmin)current.diameter-=2;
        else {current.status="processed";active=null;completedNodes++;}
        completedAttempts++;
        return new Attempt(current.id,depth,diameter,found.size(),novel,complete,fallback);
    }
}
