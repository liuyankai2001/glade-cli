package org.glade.retropath.scheduler;

import java.util.*;

/** Executed with -ea against the same policy class loaded by the KNIME nodes. */
public final class SchedulerPolicyTest {
    static SchedulerState s(int steps, int topx) {
        return new SchedulerState(steps, topx, 2, 16, true, Set.of("S"));
    }
    static SchedulerState.Product p(String id) { return new SchedulerState.Product(id, 0, true); }
    static SchedulerState.Edge e(String parent, int d, String... ids) {
        return new SchedulerState.Edge(parent, List.of(ids), "r"+d, 0, d);
    }
    static void check(boolean value, String text) { if(!value)throw new AssertionError(text); }
    static void mixed() {
        var s=s(10,100);s.offer("A",0,0,true);s.offer("B",0,0,true);
        check(s.next().id.equals("A"),"stable order");
        for(int d=16;d>10;d-=2){check(s.active.diameter==d,"diameter descending");s.finish(List.of(),List.of());}
        s.finish(List.of(e("A",10,"X")),List.of(p("X")));
        check(s.next().id.equals("B"),"B independently starts");
        check(s.active.diameter==16,"B diameter reset");
        s.finish(List.of(),List.of());s.finish(List.of(e("B",14,"Y")),List.of(p("Y")));
        check(s.nodes.containsKey("X")&&s.nodes.containsKey("Y"),"neither intermediate lost");
    }
    static void cycleFallback() {
        var s=s(10,100);s.offer("A",3,0,true);s.next();
        s.finish(List.of(e("A",16,"A")),List.of(p("A")));
        check(s.next().diameter==14,"cycle must lower diameter");
        s.finish(List.of(e("A",14,"B")),List.of(p("B")));
        check(s.nodes.get("B").depth==4,"retries do not spend chemistry depth");
        check(s.edges.size()==2,"keep cyclic edges");
    }
    static void sinkSuccess() {
        var s=s(10,100);s.offer("A",0,0,true);s.next();
        s.finish(List.of(e("A",16,"S")),List.of(p("S")));
        check(s.closed().containsKey("A")&&s.active==null,"sink-only success");
    }
    static void overflow() {
        var s=s(10,100);for(int i=0;i<101;i++)s.offer("A"+i,0,0,true);
        for(int i=0;i<100;i++){var n=s.next();s.finish(List.of(e(n.id,16,"S")),List.of(p("S")));}
        check(s.pendingCount()==1&&s.next().id.equals("A100"),"101st node retained");
    }
    static void falseCycleProof() {
        var s=s(10,100);s.edges.add(e("A",16,"B"));s.edges.add(e("B",16,"A"));
        check(!s.closed().containsKey("A"),"cycle must not prove closure");
    }
    static void depthLimit() {
        var s=s(2,100);s.edges.add(e("A",16,"B"));s.edges.add(e("B",16,"C"));s.edges.add(e("C",16,"S"));
        check(!s.closed().containsKey("A")&&s.closed().containsKey("B"),"proof beyond max steps rejected");
        s.offer("A",1,0,true);s.next();s.finish(List.of(e("A",16,"B")),List.of(p("B")));
        check(s.nodes.get("B").status.equals("depth_boundary"),"last-depth products not expanded");
    }
    static void interruption() {
        var s=s(10,100);s.offer("A",0,0,true);s.offer("B",0,0,true);s.next();
        check(s.pendingCount()==2&&s.nodes.get("A").status.equals("active"),"checkpoint includes active and unprocessed");
    }
    static void priority() {
        var s=s(10,100);s.offer("B",1,9,true);s.offer("C",1,1,true);s.offer("B",1,0,true);
        check(s.next().id.equals("B"),"later lower score improves pending order");
        var max=new SchedulerState(10,100,2,16,false,Set.of());
        max.offer("B",0,1,true);max.offer("C",0,9,true);max.offer("B",0,10,true);
        check(max.next().id.equals("B"),"original score direction respected");
    }
    static void blockedProductsDoNotStopFallback() {
        var s=s(10,100);s.offer("A",0,0,true);s.next();
        s.finish(List.of(e("A",16,"HEAVY")),List.of(new SchedulerState.Product("HEAVY",0,false)));
        check(s.next().diameter==14,"mass-filtered product is not useful progress");
    }
    public static void main(String[] args) {
        mixed();cycleFallback();sinkSuccess();overflow();falseCycleProof();depthLimit();interruption();priority();blockedProductsDoNotStopFallback();
        System.out.println("SchedulerPolicyTest: 9 scenarios PASS");
    }
}
