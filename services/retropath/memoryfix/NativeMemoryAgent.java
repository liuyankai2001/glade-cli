import java.io.IOException;
import java.lang.instrument.Instrumentation;
import java.nio.file.*;
import com.sun.jna.*;

/** Return already-free glibc pages in the actual KNIME JVM; never force GC. */
public class NativeMemoryAgent {
    @Structure.FieldOrder({"arena", "ordblks", "smblks", "hblks", "hblkhd",
            "usmblks", "fsmblks", "uordblks", "fordblks", "keepcost"})
    public static class MallInfo extends Structure implements Structure.ByValue {
        public long arena, ordblks, smblks, hblks, hblkhd;
        public long usmblks, fsmblks, uordblks, fordblks, keepcost;
    }
    public interface LibC extends Library {
        MallInfo mallinfo2();
        int malloc_trim(long padding);
    }
    static long rss() throws IOException {
        for (String line : Files.readAllLines(Path.of("/proc/self/status")))
            if (line.startsWith("VmRSS:")) return Long.parseLong(line.trim().split("\\s+")[1]) * 1024;
        throw new IOException("VmRSS unavailable");
    }
    public static void premain(String options, Instrumentation instrumentation) {
        if (!Boolean.parseBoolean(System.getProperty("glade.nativeTrim.enabled", "true"))) {
            System.err.println("GLADE_NATIVE_TRIM_DISABLED v1");
            return;
        }
        try {
            if (!Platform.isLinux() || Native.SIZE_T_SIZE != 8)
                throw new IllegalStateException("Requires the pinned 64-bit Linux/glibc runtime");
            long interval = Long.getLong("glade.nativeTrim.intervalMs", 30000L);
            long freeThreshold = Long.getLong("glade.nativeTrim.freeBytes", 536870912L);
            long rssThreshold = Long.getLong("glade.nativeTrim.rssBytes", 1610612736L);
            if (interval < 1000 || freeThreshold <= 0 || rssThreshold <= 0)
                throw new IllegalArgumentException("Invalid native trim thresholds");
            LibC libc = Native.load("c", LibC.class);
            libc.mallinfo2();
            rss();
            Thread thread = new Thread(() -> {
                try {
                    while (true) {
                        Thread.sleep(interval);
                        MallInfo value = libc.mallinfo2();
                        long before = rss();
                        if (value.fordblks >= freeThreshold && before >= rssThreshold) {
                            long start = System.nanoTime();
                            int result = libc.malloc_trim(0);
                            long after = rss();
                            System.err.println("GLADE_NATIVE_TRIM v1 before_rss=" + before
                                    + " after_rss=" + after + " free_bytes=" + value.fordblks
                                    + " result=" + result + " duration_ms=" + ((System.nanoTime() - start) / 1000000));
                        }
                    }
                } catch (InterruptedException interrupted) {
                    Thread.currentThread().interrupt();
                } catch (Throwable failure) {
                    // An optional allocator optimization must not abort a chemistry job.
                    System.err.println("GLADE_NATIVE_TRIM_STOPPED v1 " + failure);
                }
            }, "GladeNativeMemoryTrim");
            thread.setDaemon(true);
            thread.setPriority(Thread.MIN_PRIORITY);
            thread.start();
            System.err.println("GLADE_NATIVE_TRIM_STARTED v1 interval_ms=" + interval
                    + " free_threshold=" + freeThreshold + " rss_threshold=" + rssThreshold);
        } catch (Throwable failure) {
            System.err.println("GLADE_NATIVE_TRIM_UNAVAILABLE v1 " + failure);
        }
    }
}
