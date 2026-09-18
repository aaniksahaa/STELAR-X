#!/usr/bin/env python3
"""Deterministic randomized differential tests for STELAR-X fixed-tree scoring.

The expected value is computed by direct rooted-triple enumeration.  This test
does not reuse STELAR-X clusters, partitions, hashes, or weight formulas.
"""

from __future__ import annotations

import argparse
import itertools
import os
import pathlib
import random
import re
import subprocess
import tempfile

ROOT = pathlib.Path(__file__).resolve().parents[1]
METHODS = ("I1", "I2", "I3", "I4")


class Node:
    def __init__(self, name: str | None = None, children: list["Node"] | None = None):
        self.name = name
        self.children = children or []
        self.taxa: frozenset[str] = frozenset()


def leaf(name: str) -> Node:
    return Node(name=name)


def annotate(node: Node) -> frozenset[str]:
    if node.name is not None:
        node.taxa = frozenset((node.name,))
    else:
        node.taxa = frozenset().union(*(annotate(child) for child in node.children))
    return node.taxa


def newick(node: Node, mirror: bool = False) -> str:
    if node.name is not None:
        return node.name
    children = list(reversed(node.children)) if mirror else node.children
    return "(" + ",".join(newick(child, mirror) for child in children) + ")"


def parse_newick(text: str) -> Node:
    """Small independent parser for the generated, undecorated Newick subset."""
    stack: list[Node | str] = []
    token = ""
    for char in text.strip().rstrip(";"):
        if char == "(":
            stack.append(char)
        elif char in ",)":
            if token.strip():
                stack.append(leaf(token.strip()))
                token = ""
            if char == ")":
                children: list[Node] = []
                while stack and stack[-1] != "(":
                    item = stack.pop()
                    if not isinstance(item, Node):
                        raise ValueError("malformed Newick")
                    children.append(item)
                if not stack:
                    raise ValueError("unbalanced Newick")
                stack.pop()
                stack.append(Node(children=list(reversed(children))))
        else:
            token += char
    if token.strip():
        stack.append(leaf(token.strip()))
    if len(stack) != 1 or not isinstance(stack[0], Node):
        raise ValueError("malformed Newick")
    annotate(stack[0])
    return stack[0]


def resolved_pair(node: Node, triple: tuple[str, str, str]) -> tuple[str, str] | None:
    if node.name is not None:
        return None
    locations = [next(i for i, child in enumerate(node.children) if taxon in child.taxa)
                 for taxon in triple]
    if locations[0] == locations[1] != locations[2]:
        return tuple(sorted((triple[0], triple[1])))
    if locations[0] == locations[2] != locations[1]:
        return tuple(sorted((triple[0], triple[2])))
    if locations[1] == locations[2] != locations[0]:
        return tuple(sorted((triple[1], triple[2])))
    if locations[0] == locations[1] == locations[2]:
        return resolved_pair(node.children[locations[0]], triple)
    return None


def oracle_score(species: Node, genes: list[Node]) -> int:
    total = 0
    for gene in genes:
        for triple in itertools.combinations(sorted(species.taxa & gene.taxa), 3):
            total += resolved_pair(species, triple) == resolved_pair(gene, triple) \
                and resolved_pair(gene, triple) is not None
    return total


def random_subtree(names: list[str], rng: random.Random, poly_probability: float) -> Node:
    if len(names) == 1:
        return leaf(names[0])
    shuffled = names[:]
    rng.shuffle(shuffled)
    if len(names) >= 3 and rng.random() < poly_probability:
        arity = rng.randint(3, min(5, len(names)))
    else:
        arity = 2
    cuts = sorted(rng.sample(range(1, len(names)), arity - 1))
    groups: list[list[str]] = []
    start = 0
    for end in cuts + [len(names)]:
        groups.append(shuffled[start:end])
        start = end
    return Node(children=[random_subtree(group, rng, poly_probability) for group in groups])


def random_rooted_tree(names: list[str], rng: random.Random,
                       poly_probability: float = 0.0) -> Node:
    """Build an explicitly rooted tree: the supplied top-level root is binary."""
    shuffled = names[:]
    rng.shuffle(shuffled)
    cut = rng.randint(1, len(shuffled) - 1)
    root = Node(children=[random_subtree(shuffled[:cut], rng, poly_probability),
                          random_subtree(shuffled[cut:], rng, poly_probability)])
    annotate(root)
    return root


class Case:
    def __init__(self, name: str, species: Node, genes: list[Node]):
        self.name = name
        self.species = species
        self.genes = genes
        self.expected = oracle_score(species, genes)


def fixed_cases() -> list[Case]:
    raw = [
        ("minimal-three", "(A,(B,C));", ["(A,(B,C));", "(B,(A,C));"]),
        ("balanced-duplicates", "((A,B),(C,D));",
         ["((A,B),(C,D));", "((A,B),(C,D));", "((A,C),(B,D));"]),
        ("complete-internal-polytomy", "((A,B),((C,D),(E,F)));",
         ["((A,B,C),(D,(E,F)));", "((A,(B,C)),(D,E,F));",
          "(((A,B),C),(D,(E,F)));"]),
        ("incomplete-polytomies", "(((A,B),C),((D,E),(F,G)));",
         ["((A,B,C),(D,E));", "((A,C),(F,D,E));", "((B,C),(D,(F,G)));",
          "((A,B),(E,G));"]),
        ("deep-and-shallow", "(A,(B,(C,(D,(E,(F,G))))));",
         ["(A,(B,(C,(D,(E,(F,G))))));", "(((A,B,C),D),(E,F,G));",
          "((A,(C,E)),(B,D,F,G));"]),
        # These two internal polytomies have the same unordered set of four
        # groups, but a different group is outside the polytomous node.  A
        # rooted partition key must therefore keep the complement distinguished.
        ("polytomy-complement-dedup", "(((a1,b),(a2,c)),d);",
         ["(((a1,a2),b,c),d);", "(((a1,a2),b,d),c);"]),
    ]
    cases = []
    for name, species_text, gene_texts in raw:
        cases.append(Case(name, parse_newick(species_text),
                          [parse_newick(text) for text in gene_texts]))
    return cases


def generated_cases(count: int, seed: int) -> list[Case]:
    rng = random.Random(seed)
    cases: list[Case] = []
    for index in range(count):
        n = rng.randint(4, 10)
        taxa = [f"T{i}" for i in range(n)]
        species = random_rooted_tree(taxa, rng)
        genes = [random_rooted_tree(taxa, rng, rng.uniform(0.0, 0.75))]
        for _ in range(rng.randint(3, 10)):
            kept = [taxon for taxon in taxa if rng.random() >= rng.uniform(0.05, 0.45)]
            if len(kept) < 3:
                kept = rng.sample(taxa, 3)
            genes.append(random_rooted_tree(kept, rng, rng.uniform(0.0, 0.85)))
        if index % 4 == 0:
            genes.extend((genes[0], genes[-1]))  # frequency/deduplication coverage
        cases.append(Case(f"random-{index:02d}-n{n}-g{len(genes)}", species, genes))

    # Explicit metamorphic variants: line order and child order cannot alter a score;
    # duplicating all observations must exactly double it.
    if cases:
        base = cases[0]
        cases.append(Case("metamorphic-reversed-lines", base.species, list(reversed(base.genes))))
        mirrored_species = parse_newick(newick(base.species, mirror=True) + ";")
        mirrored_genes = [parse_newick(newick(gene, mirror=True) + ";") for gene in base.genes]
        cases.append(Case("metamorphic-mirrored-children", mirrored_species, mirrored_genes))
        cases.append(Case("metamorphic-doubled-genes", base.species, base.genes + base.genes))
        assert cases[-1].expected == 2 * base.expected
    return cases


def run_score(case: Case, method: str, work: pathlib.Path, compute: str,
              numeric: str = "long", threads: int | None = None) -> int:
    genes_path = work / f"{case.name}.genes.tre"
    species_path = work / f"{case.name}.species.tre"
    if not genes_path.exists():
        genes_path.write_text("".join(newick(gene) + ";\n" for gene in case.genes))
        species_path.write_text(newick(case.species) + ";\n")

    command = [
        "java", "-Djava.library.path=" + str(ROOT / "native"),
        "-cp", str(ROOT / "build"), "stelarx.Main", compute, "-q",
        "-i", str(genes_path), "--score-species-tree", str(species_path),
        "--intersection-method", method,
    ]
    env = os.environ.copy()
    if numeric != "long":
        env["STELARX_WEIGHT_FORCE_DOUBLE"] = "1"
        command += ["--large-n-score-type", numeric]
    if threads is not None:
        command += ["--threads", str(threads)]
    run = subprocess.run(command, cwd=ROOT, env=env, text=True,
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if run.returncode:
        raise AssertionError(f"{case.name}/{method}/{compute}/{numeric} failed:\n{run.stdout}")
    match = re.search(r"^TRIPLET_SCORE:\s*([0-9]+(?:[.]0+)?)\s*$",
                      run.stdout, re.MULTILINE)
    if not match:
        raise AssertionError(f"{case.name}/{method} emitted no score:\n{run.stdout}")
    if compute == "--gpu-strict" and "[STELAR-X GPU] weight" not in run.stdout:
        raise AssertionError(f"{case.name}/{method} did not execute a GPU weight phase")
    return int(float(match.group(1)))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cases", type=int, default=12,
                        help="number of seeded random cases (plus fixed/metamorphic cases)")
    parser.add_argument("--seed", type=lambda value: int(value, 0), default=0x5E1A7)
    parser.add_argument("--gpu", action="store_true", help="require strict CUDA execution")
    parser.add_argument("--no-build", action="store_true")
    parser.add_argument("--skip-numeric", action="store_true",
                        help="skip forced Int128/Double path checks")
    args = parser.parse_args()
    if args.cases < 0:
        parser.error("--cases must be non-negative")
    if not args.no_build:
        subprocess.run([str(ROOT / "build.sh")], cwd=ROOT, check=True,
                       stdout=subprocess.DEVNULL)

    cases = fixed_cases() + generated_cases(args.cases, args.seed)
    compute = "--gpu-strict" if args.gpu else "--cpu"
    assertions = 0
    with tempfile.TemporaryDirectory(prefix="stelarx-differential-") as directory:
        work = pathlib.Path(directory)
        for case in cases:
            scores = {method: run_score(case, method, work, compute) for method in METHODS}
            if set(scores.values()) != {case.expected}:
                raise AssertionError(f"{case.name}: oracle={case.expected}, actual={scores}")
            assertions += len(scores)

        # Exercise serial and parallel CPU implementations on representative binary,
        # polytomous, and incomplete inputs. GPU kernels are thread-count independent.
        if not args.gpu:
            for case in (cases[0], cases[2], cases[3], cases[-1]):
                for threads in (1, min(4, os.cpu_count() or 1)):
                    actual = run_score(case, "I2", work, compute, threads=threads)
                    if actual != case.expected:
                        raise AssertionError(f"{case.name}: threads={threads}, score={actual}")
                    assertions += 1

        # Force the normally large-N-only arithmetic paths on small exact oracles.
        if not args.skip_numeric:
            for case in (cases[2], cases[3], cases[5]):
                for numeric in ("int128", "double"):
                    for method in METHODS:
                        actual = run_score(case, method, work, compute, numeric=numeric)
                        if actual != case.expected:
                            raise AssertionError(
                                f"{case.name}/{method}/{numeric}: "
                                f"oracle={case.expected}, actual={actual}")
                        assertions += 1

    print(f"STELAR-X differential oracle: PASS ({len(cases)} cases, "
          f"{assertions} scorer executions, seed=0x{args.seed:x}, "
          f"compute={'CUDA-strict' if args.gpu else 'CPU'})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
