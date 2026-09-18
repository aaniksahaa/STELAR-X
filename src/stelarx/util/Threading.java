package stelarx.util;

import java.util.List;
import java.util.concurrent.*;
import java.util.function.Consumer;

public class Threading {
    private static ExecutorService executor;
    private static int numThreads;

    public static void start(int threads) {
        numThreads = Math.max(1, threads);
        executor = Executors.newFixedThreadPool(numThreads);
    }

    public static void shutdown() {
        if (executor != null) {
            executor.shutdown();
            try { executor.awaitTermination(60, TimeUnit.SECONDS); }
            catch (InterruptedException e) { executor.shutdownNow(); Thread.currentThread().interrupt(); }
            finally { executor = null; numThreads = 0; }
        }
    }

    public static int getNumThreads() { return numThreads; }
    public static boolean isStarted() { return executor != null; }
    public static Future<?> submit(Runnable r) { return executor.submit(r); }
    public static <T> Future<T> submit(Callable<T> c) { return executor.submit(c); }

    /** Divide list across threads, apply action in parallel, wait for all. */
    public static <T> void processParallel(List<T> items, Consumer<T> action) {
        if (items.isEmpty()) return;
        int chunk = Math.max(1, (items.size() + numThreads - 1) / numThreads);
        int actual = Math.min(numThreads, (items.size() + chunk - 1) / chunk);
        CompletionService<Void> completed = new ExecutorCompletionService<>(executor);
        List<Future<Void>> futures = new java.util.ArrayList<>(actual);
        for (int t = 0; t < actual; t++) {
            int lo = t * chunk, hi = Math.min(lo + chunk, items.size());
            futures.add(completed.submit(() -> {
                for (int i = lo; i < hi; i++) action.accept(items.get(i));
                return null;
            }));
        }
        awaitWorkers(completed, futures, actual);
    }

    /** Divide range [0,count) across threads in parallel, wait for all. */
    public static void processRangeParallel(int count, Consumer<Integer> action) {
        if (count == 0) return;
        int chunk = Math.max(1, (count + numThreads - 1) / numThreads);
        int actual = Math.min(numThreads, (count + chunk - 1) / chunk);
        CompletionService<Void> completed = new ExecutorCompletionService<>(executor);
        List<Future<Void>> futures = new java.util.ArrayList<>(actual);
        for (int t = 0; t < actual; t++) {
            int lo = t * chunk, hi = Math.min(lo + chunk, count);
            futures.add(completed.submit(() -> {
                for (int i = lo; i < hi; i++) action.accept(i);
                return null;
            }));
        }
        awaitWorkers(completed, futures, actual);
    }

    /** Wait for worker completion and surface the original worker failure. */
    private static void awaitWorkers(CompletionService<Void> completed,
                                     List<Future<Void>> futures,
                                     int taskCount) {
        Throwable firstFailure = null;
        try {
            for (int i = 0; i < taskCount; i++) {
                try {
                    completed.take().get();
                } catch (ExecutionException e) {
                    if (firstFailure == null) firstFailure = e.getCause();
                }
            }
        } catch (InterruptedException e) {
            for (Future<Void> future : futures) future.cancel(true);
            Thread.currentThread().interrupt();
            throw new RuntimeException("Interrupted while waiting for parallel workers", e);
        }
        if (firstFailure instanceof RuntimeException runtime) throw runtime;
        if (firstFailure instanceof Error error) throw error;
        if (firstFailure != null) throw new RuntimeException("Parallel worker failed", firstFailure);
    }
}
