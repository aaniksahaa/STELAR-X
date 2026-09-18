#!/usr/bin/env python3
"""Accuracy and resource-regression benchmark for STELAR-X.

This is deliberately a small, repeatable development-host benchmark rather
than a noisy microbenchmark. It validates every S×I combination, exact scaling
under replicated gene observations, CPU/CUDA parity, and optionally compares
S1/I2 against a separately built reference checkout.
"""

from __future__ import annotations

import argparse
import json
import pathlib
import re
import statistics
import subprocess
import tempfile

from test_stelarx_differential import oracle_score, parse_newick

ROOT = pathlib.Path(__file__).resolve().parents[1]
GENES = ROOT / "example" / "all_gt_37.tre"
SPECIES = ROOT / "example" / "true_37.tre"
METHODS = ("I1", "I2", "I3", "I4")
PRESETS = ("S1", "S2", "S3")
TIME = pathlib.Path("/usr/bin/time")


def parse_decorated_newick(text: str):
    """Remove branch lengths before feeding the independent topology parser."""
    return parse_newick(re.sub(r":[^,();\s]+", "", text))


def score_from_output(output: str) -> int:
    patterns = (
        r"^TRIPLET_SCORE:\s*([0-9]+)\s*$",
        r"(?:Final triplet score\s*=|Triplet score\s+)([0-9]+)",
    )
    for pattern in patterns:
        matches = re.findall(pattern, output, re.MULTILINE | re.IGNORECASE)
        if matches:
            return int(matches[-1])
    raise AssertionError("program output contained no rooted-triplet score")


def timed_run(command: list[str], metric_path: pathlib.Path,
              cwd: pathlib.Path, require_gpu: bool = False) -> dict:
    wrapped = [str(TIME), "-f", "%e\t%M", "-o", str(metric_path), *command]
    run = subprocess.run(wrapped, cwd=cwd, text=True,
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if run.returncode:
        raise AssertionError(f"command failed ({run.returncode}):\n{' '.join(command)}\n{run.stdout}")
    elapsed_text, rss_text = metric_path.read_text().strip().split("\t")
    if require_gpu and "[STELAR-X GPU] weight" not in run.stdout:
        raise AssertionError("strict CUDA run did not execute a GPU weight phase")
    return {
        "seconds": float(elapsed_text),
        "rss_kib": int(rss_text),
        "score": score_from_output(run.stdout),
        "stdout": run.stdout,
    }


def java_command(root: pathlib.Path, compute: str, *args: str) -> list[str]:
    return [
        "java", "-Djava.library.path=" + str(root / "native"),
        "-cp", str(root / "build"), "stelarx.Main", compute, "-q", *args,
    ]


def median_benchmark(root: pathlib.Path, compute: str, repeats: int,
                     work: pathlib.Path, label: str, extra: tuple[str, ...] = ()) -> dict:
    samples = []
    trees = []
    for repeat in range(repeats):
        output_tree = work / f"{label}-{repeat}.tre"
        command = java_command(
            root, compute, "-i", str(GENES), "-o", str(output_tree),
            "--search-space", "S1", "--intersection-method", "I2", *extra)
        samples.append(timed_run(command, work / f"{label}-{repeat}.metrics", root,
                                 require_gpu=compute == "--gpu-strict"))
        trees.append(output_tree.read_text().strip())
    if len(set(trees)) != 1 or len({sample["score"] for sample in samples}) != 1:
        raise AssertionError(f"{label} was nondeterministic across repeats")
    return {
        "seconds_median": statistics.median(sample["seconds"] for sample in samples),
        "rss_kib_median": int(statistics.median(sample["rss_kib"] for sample in samples)),
        "score": samples[0]["score"],
        "tree": trees[0],
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--require-gpu", action="store_true")
    parser.add_argument("--reference-dir", type=pathlib.Path,
                        help="built reference checkout for S1/I2 resource comparison")
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--scales", type=int, nargs="+", default=(1, 4, 16))
    parser.add_argument("--json", type=pathlib.Path, help="optional result file")
    args = parser.parse_args()
    if args.repeats < 1 or any(scale < 1 for scale in args.scales):
        parser.error("repeats and scales must be positive")
    if not TIME.is_file():
        parser.error("/usr/bin/time is required for peak-RSS measurement")

    results: dict = {"all_modes": {}, "scaling": {}, "baseline": {}}
    large_genes = [parse_decorated_newick(line) for line in GENES.read_text().splitlines()
                   if line.strip()]
    oracle = oracle_score(parse_decorated_newick(SPECIES.read_text()), large_genes)
    if oracle != 1_390_544:
        raise AssertionError(f"37-taxon independent oracle changed: {oracle}")
    results["independent_fixed_tree_score"] = oracle

    with tempfile.TemporaryDirectory(prefix="stelarx-scalability-") as directory:
        work = pathlib.Path(directory)
        current = median_benchmark(ROOT, "--cpu", args.repeats, work, "current")
        results["baseline"]["stelarx"] = {k: v for k, v in current.items() if k != "tree"}

        if args.reference_dir is not None:
            reference_root = args.reference_dir.resolve()
            reference = median_benchmark(
                reference_root, "--cpu", args.repeats, work, "reference",
                ("--rooted", "--no-anchor-outgroup"))
            results["baseline"]["reference"] = {
                k: v for k, v in reference.items() if k != "tree"
            }
            time_ratio = current["seconds_median"] / max(reference["seconds_median"], 0.001)
            rss_ratio = current["rss_kib_median"] / max(reference["rss_kib_median"], 1)
            results["baseline"]["time_ratio"] = time_ratio
            results["baseline"]["rss_ratio"] = rss_ratio
            # Wide enough for normal host noise, strict enough to catch a material regression.
            if time_ratio > 1.35:
                raise AssertionError(f"S1/I2 CPU time regression: ratio={time_ratio:.3f}")
            if rss_ratio > 1.25:
                raise AssertionError(f"S1/I2 CPU RSS regression: ratio={rss_ratio:.3f}")

        computes = ["--cpu"] + (["--gpu-strict"] if args.require_gpu else [])
        canonical: dict[tuple[str, str], tuple[int, str]] = {}
        for compute in computes:
            compute_name = "gpu" if compute == "--gpu-strict" else "cpu"
            for preset in PRESETS:
                preset_tree = None
                preset_score = None
                for method in METHODS:
                    label = f"mode-{compute_name}-{preset}-{method}"
                    output_tree = work / f"{label}.tre"
                    command = java_command(
                        ROOT, compute, "-i", str(GENES), "-o", str(output_tree),
                        "--search-space", preset, "--intersection-method", method)
                    sample = timed_run(command, work / f"{label}.metrics", ROOT,
                                       require_gpu=compute == "--gpu-strict")
                    tree = output_tree.read_text().strip()
                    results["all_modes"][label] = {
                        key: sample[key] for key in ("seconds", "rss_kib", "score")
                    }
                    if preset_tree is None:
                        preset_tree, preset_score = tree, sample["score"]
                    elif tree != preset_tree or sample["score"] != preset_score:
                        raise AssertionError(f"{compute_name}/{preset} differs across I1-I4")
                    if compute_name == "cpu":
                        canonical[(preset, method)] = (sample["score"], tree)
                    elif canonical[(preset, method)] != (sample["score"], tree):
                        raise AssertionError(f"CPU/GPU mismatch for {preset}/{method}")
                if compute_name == "cpu":
                    inferred_oracle = oracle_score(parse_decorated_newick(preset_tree), large_genes)
                    if inferred_oracle != preset_score:
                        raise AssertionError(
                            f"independent inferred-tree score mismatch for {preset}: "
                            f"{inferred_oracle} != {preset_score}")
                    results["all_modes"][f"oracle-{preset}"] = {
                        "score": inferred_oracle
                    }

        base_lines = [line for line in GENES.read_text().splitlines() if line.strip()]
        for scale in args.scales:
            scaled_path = work / f"genes-x{scale}.tre"
            scaled_path.write_text("\n".join(base_lines * scale) + "\n")
            for compute in computes:
                compute_name = "gpu" if compute == "--gpu-strict" else "cpu"
                base_rss = None
                base_time = None
                for method in METHODS:
                    label = f"scale-{compute_name}-x{scale}-{method}"
                    command = java_command(
                        ROOT, compute, "-i", str(scaled_path),
                        "--score-species-tree", str(SPECIES),
                        "--intersection-method", method)
                    sample = timed_run(command, work / f"{label}.metrics", ROOT,
                                       require_gpu=compute == "--gpu-strict")
                    if sample["score"] != oracle * scale:
                        raise AssertionError(
                            f"nonlinear score at x{scale}/{compute_name}/{method}: "
                            f"{sample['score']} != {oracle * scale}")
                    results["scaling"][label] = {
                        key: sample[key] for key in ("seconds", "rss_kib", "score")
                    }
                    if scale == min(args.scales):
                        base_rss = sample["rss_kib"] if base_rss is None else max(base_rss, sample["rss_kib"])
                        base_time = sample["seconds"] if base_time is None else max(base_time, sample["seconds"])

        # Aggregate scaling guards across methods. JVM startup dominates x1, so
        # use generous bounds that detect explosions without making timing flaky.
        for compute_name in ("cpu", "gpu") if args.require_gpu else ("cpu",):
            low = min(args.scales)
            high = max(args.scales)
            low_rows = [results["scaling"][f"scale-{compute_name}-x{low}-{m}"] for m in METHODS]
            high_rows = [results["scaling"][f"scale-{compute_name}-x{high}-{m}"] for m in METHODS]
            low_rss = max(row["rss_kib"] for row in low_rows)
            high_rss = max(row["rss_kib"] for row in high_rows)
            low_time = max(row["seconds"] for row in low_rows)
            high_time = max(row["seconds"] for row in high_rows)
            if high_rss > low_rss + 512 * 1024:
                raise AssertionError(f"{compute_name} RSS grew by more than 512 MiB")
            if high_time > low_time * (2.0 * high / low) + 2.0:
                raise AssertionError(f"{compute_name} time scaling is unexpectedly superlinear")

    if args.json:
        args.json.write_text(json.dumps(results, indent=2, sort_keys=True) + "\n")

    print("STELAR-X scalability/accuracy benchmark: PASS")
    print(f"  independent 37-taxon fixed-tree score: {oracle}")
    print(f"  current S1/I2 CPU median: {current['seconds_median']:.3f}s, "
          f"{current['rss_kib_median'] / 1024:.1f} MiB RSS ({args.repeats} runs)")
    if args.reference_dir is not None:
        reference = results["baseline"]["reference"]
        print(f"  reference S1/I2 CPU median: {reference['seconds_median']:.3f}s, "
              f"{reference['rss_kib_median'] / 1024:.1f} MiB RSS")
        print(f"  current/reference: time={results['baseline']['time_ratio']:.3f}x, "
              f"RSS={results['baseline']['rss_ratio']:.3f}x")
    print(f"  all-mode parity: {len(computes) * len(PRESETS) * len(METHODS)} "
          "S×I×compute runs")
    print(f"  scaling parity: {len(results['scaling'])} fixed-tree runs at "
          + ", ".join(f"x{scale}" for scale in args.scales))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
