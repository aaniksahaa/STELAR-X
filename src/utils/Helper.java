package utils;

import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicLong;

public class Helper {
    public String greet(String name) {
        return "Hello, " + name + "!";
    }
    
    /**
     * Calculates the sum from 1 to n using parallel computation
     * The work is divided among as many threads as there are CPU cores
     */
    public long parallelSum(long n) {
        if (n <= 0) {
            return 0;
        }
        
        int numThreads = Runtime.getRuntime().availableProcessors();
        System.out.println("Using " + numThreads + " threads for parallel computation");
        
        ExecutorService executor = Executors.newFixedThreadPool(numThreads);
        AtomicLong totalSum = new AtomicLong(0);
        
        // Calculate range per thread
        long rangePerThread = n / numThreads;
        long remainder = n % numThreads;
        
        CountDownLatch latch = new CountDownLatch(numThreads);
        
        for (int i = 0; i < numThreads; i++) {
            long start = i * rangePerThread + 1;
            long end = (i + 1) * rangePerThread;
            
            // Add remainder to the last thread
            if (i == numThreads - 1) {
                end += remainder;
            }
            
            final long threadStart = start;
            final long threadEnd = end;
            final int threadId = i;
            
            executor.submit(() -> {
                try {
                    long partialSum = 0;
                    for (long j = threadStart; j <= threadEnd; j++) {
                        partialSum += j;
                    }
                    
                    totalSum.addAndGet(partialSum);
                    System.out.println("Thread " + threadId + " computed sum from " + 
                                     threadStart + " to " + threadEnd + " = " + partialSum);
                } finally {
                    latch.countDown();
                }
            });
        }
        
        try {
            latch.await(); // Wait for all threads to complete
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            throw new RuntimeException("Computation was interrupted", e);
        } finally {
            executor.shutdown();
        }
        
        return totalSum.get();
    }
    
    /**
     * Sequential sum for comparison purposes
     */
    public long sequentialSum(long n) {
        long sum = 0;
        for (long i = 1; i <= n; i++) {
            sum += i;
        }
        return sum;
    }
    
    /**
     * Mathematical formula for sum (for verification)
     */
    public long formulaSum(long n) {
        return n * (n + 1) / 2;
    }
}