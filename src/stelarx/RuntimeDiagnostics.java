package stelarx;

import stelarx.gpu.GPUDistanceMatrix;
import stelarx.gpu.GPUDPBuilder;
import stelarx.gpu.GPUSimilarityMatrix;
import stelarx.gpu.GPUWeightCalculator;

import java.nio.file.Files;
import java.nio.file.Path;

/** Human-readable, input-free installation and hardware self-check. */
public final class RuntimeDiagnostics {
    private RuntimeDiagnostics() {}

    public static void print(Config cfg) {
        Runtime rt = Runtime.getRuntime();
        GPUWeightCalculator.Probe gpu = GPUWeightCalculator.probe();

        System.out.println("STELAR-X DIAGNOSTICS");
        row("Version", Main.VERSION);
        row("OS", prop("os.name") + " " + prop("os.version"));
        row("Architecture", prop("os.arch"));
        row("Java runtime", prop("java.runtime.name") + " " + prop("java.runtime.version"));
        row("Bundled runtime home", prop("java.home"));
        row("Application image", System.getProperty("jpackage.app-path", "not packaged"));
        row("Working directory", prop("user.dir"));
        row("CPU cores", Integer.toString(rt.availableProcessors()));
        row("Maximum JVM memory", formatBytes(rt.maxMemory()));
        row("Requested compute", cfg.getRequestedComputeMode().name());
        row("Selected compute", cfg.getComputeMode().name());
        row("Selection detail", cfg.getComputeModeDetail());
        row("Native library path", prop("java.library.path"));

        System.out.println();
        System.out.println("Native backends");
        backend("weight / CUDA probe", gpu.libraryLoaded(),
            gpu.libraryLoaded() ? gpu.detail() : GPUWeightCalculator.getLoadError());
        backend("cross-tree DP", GPUDPBuilder.tryLoad(), GPUDPBuilder.getLoadError());
        backend("distance matrix", GPUDistanceMatrix.tryLoad(), GPUDistanceMatrix.getLoadError());
        backend("similarity matrix", GPUSimilarityMatrix.tryLoad(), GPUSimilarityMatrix.getLoadError());

        System.out.println();
        System.out.println("CUDA");
        row("Usable", gpu.cudaAvailable() ? "yes" : "no (CPU fallback is available)");
        if (gpu.cudaAvailable()) {
            row("Device", gpu.deviceName());
            row("Compute capability", gpu.computeMajor() + "." + gpu.computeMinor());
            row("CUDA devices", Integer.toString(gpu.deviceCount()));
            row("Driver API version", cudaVersion(gpu.driverVersion()));
            row("Bundled runtime version", cudaVersion(gpu.runtimeVersion()));
            row("VRAM", gpu.freeMiB() + " MiB free / " + gpu.totalMiB() + " MiB total");
        } else {
            row("Reason", gpu.detail());
        }

        System.out.println();
        System.out.println("Filesystem");
        Path cwd = Path.of(prop("user.dir"));
        row("Current directory readable", Boolean.toString(Files.isReadable(cwd)));
        row("Current directory writable", Boolean.toString(Files.isWritable(cwd)));
        System.out.println();
        System.out.println("Result: " + (cfg.getComputeMode() == Config.ComputeMode.GPU
            ? "GPU execution is ready." : "CPU execution is ready."));
    }

    private static void backend(String name, boolean loaded, String detail) {
        row(name, loaded ? "loaded" : "unavailable");
        if (!loaded && detail != null && !detail.isBlank()) row("  reason", detail);
    }

    private static String prop(String key) {
        return System.getProperty(key, "unknown");
    }

    private static String cudaVersion(int v) {
        if (v <= 0) return "unknown";
        return (v / 1000) + "." + ((v % 1000) / 10);
    }

    private static String formatBytes(long bytes) {
        if (bytes >= (1L << 30)) return String.format("%.1f GiB", bytes / (double)(1L << 30));
        return String.format("%.1f MiB", bytes / (double)(1L << 20));
    }

    private static void row(String label, String value) {
        System.out.printf("  %-28s %s%n", label + ":", value == null ? "" : value);
    }
}
