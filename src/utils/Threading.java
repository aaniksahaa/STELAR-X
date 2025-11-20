package utils;

import java.util.concurrent.Callable;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.List;
import java.util.ArrayList;
import java.util.function.Consumer;
import java.util.function.Function;

// import org.jocl.cl_context;
// import org.jocl.cl_context_properties;
// import org.jocl.cl_device_id;

public class Threading {

	private static ExecutorService eService;
	// static cl_context_properties contextProperties;
	// static cl_context [] context;
	// static cl_device_id[] usedDevices;
	static String[] deviceVendors;
	private static int numThreads = Runtime.getRuntime().availableProcessors();
	private static int distMatrixChunks = -1;
	
	public static Future submit(Callable c) {
		return Threading.eService.submit(c);
	}
	
	public static void execute(Runnable c) {
		Threading.eService.execute(c);
	}

	public static void shutdown() {
		System.out.println("Shutting down threading");
		if (Threading.eService != null)
			Threading.eService.shutdown();
	}

	public static void startThreading(int t) {
		Threading.numThreads = t;
		if (Threading.numThreads<2) {
			throw new RuntimeException("Sorry, at least two threads are needed."); 
		}
		Threading.eService = Executors.newFixedThreadPool(Threading.numThreads);
		System.out.println("There are " + Threading.getNumThreads() + " threads used to run.");
		
	}

	public static int getNumThreads() {
		return numThreads;
	}
	
	public static int getDistMatrixChunkSize() {
		return distMatrixChunks;
	}

	static void setDistMatrixChunkSize(int chunks) {
		distMatrixChunks  = chunks;
	}
	
	/**
	 * Utility method to process a list in parallel using multiple threads.
	 * Each thread processes a chunk of the list.
	 * 
	 * @param <T> The type of items in the list
	 * @param items The list of items to process
	 * @param processor Function to process each item
	 * @param verbose Whether to print progress messages
	 */
	public static <T> void processListParallel(List<T> items, Consumer<T> processor, boolean verbose) {
		if (items.isEmpty()) {
			return;
		}
		
		int numThreads = getNumThreads();
		if (numThreads < 2) {
			// Fallback to sequential processing
			if (verbose) {
				System.out.println("Processing " + items.size() + " items sequentially...");
			}
			for (T item : items) {
				processor.accept(item);
			}
			return;
		}
		
		if (verbose) {
			System.out.println("Processing " + items.size() + " items in parallel with " + numThreads + " threads...");
		}
		
		startThreading(numThreads);
		
		int chunkSize = (items.size() + numThreads - 1) / numThreads;
		CountDownLatch latch = new CountDownLatch(numThreads);
		
		for (int i = 0; i < numThreads; i++) {
			final int startIdx = i * chunkSize;
			final int endIdx = Math.min(startIdx + chunkSize, items.size());
			final int threadId = i;
			
			if (startIdx >= items.size()) {
				latch.countDown();
				continue;
			}
			
			execute(() -> {
				try {
					if (verbose) {
						System.out.println("Thread " + threadId + " processing items " + startIdx + " to " + (endIdx-1));
					}
					
					for (int j = startIdx; j < endIdx; j++) {
						processor.accept(items.get(j));
					}
					
					if (verbose) {
						System.out.println("Thread " + threadId + " completed processing " + (endIdx - startIdx) + " items");
					}
				} finally {
					latch.countDown();
				}
			});
		}
		
		try {
			latch.await();
		} catch (InterruptedException e) {
			Thread.currentThread().interrupt();
			throw new RuntimeException("Parallel processing was interrupted", e);
		} finally {
			shutdown();
		}
		
		if (verbose) {
			System.out.println("Parallel processing completed.");
		}
	}
	
	/**
	 * Utility method to process a list in parallel and collect results.
	 * Each thread processes a chunk of the list and returns results.
	 * 
	 * @param <T> The type of items in the list
	 * @param <R> The type of results
	 * @param items The list of items to process
	 * @param processor Function to process each item and return a result
	 * @param verbose Whether to print progress messages
	 * @return List of results from processing
	 */
	public static <T, R> List<R> processListParallelWithResults(List<T> items, Function<T, R> processor, boolean verbose) {
		List<R> results = new ArrayList<>();
		
		if (items.isEmpty()) {
			return results;
		}
		
		int numThreads = getNumThreads();
		if (numThreads < 2) {
			// Fallback to sequential processing
			if (verbose) {
				System.out.println("Processing " + items.size() + " items sequentially...");
			}
			for (T item : items) {
				R result = processor.apply(item);
				if (result != null) {
					results.add(result);
				}
			}
			return results;
		}
		
		if (verbose) {
			System.out.println("Processing " + items.size() + " items in parallel with " + numThreads + " threads...");
		}
		
		startThreading(numThreads);
		
		int chunkSize = (items.size() + numThreads - 1) / numThreads;
		CountDownLatch latch = new CountDownLatch(numThreads);
		List<List<R>> threadResults = new ArrayList<>(numThreads);
		
		// Initialize thread results lists
		for (int i = 0; i < numThreads; i++) {
			threadResults.add(new ArrayList<>());
		}
		
		for (int i = 0; i < numThreads; i++) {
			final int startIdx = i * chunkSize;
			final int endIdx = Math.min(startIdx + chunkSize, items.size());
			final int threadId = i;
			final List<R> threadResultList = threadResults.get(i);
			
			if (startIdx >= items.size()) {
				latch.countDown();
				continue;
			}
			
			execute(() -> {
				try {
					if (verbose) {
						System.out.println("Thread " + threadId + " processing items " + startIdx + " to " + (endIdx-1));
					}
					
					for (int j = startIdx; j < endIdx; j++) {
						R result = processor.apply(items.get(j));
						if (result != null) {
							threadResultList.add(result);
						}
					}
					
					if (verbose) {
						System.out.println("Thread " + threadId + " completed: processed " + (endIdx - startIdx) + 
										 " items, generated " + threadResultList.size() + " results");
					}
				} finally {
					latch.countDown();
				}
			});
		}
		
		try {
			latch.await();
		} catch (InterruptedException e) {
			Thread.currentThread().interrupt();
			throw new RuntimeException("Parallel processing was interrupted", e);
		} finally {
			shutdown();
		}
		
		// Combine results from all threads
		for (List<R> threadResultList : threadResults) {
			results.addAll(threadResultList);
		}
		
		if (verbose) {
			System.out.println("Parallel processing completed. Total results: " + results.size());
		}
		
		return results;
	}
}
