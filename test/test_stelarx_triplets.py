#!/usr/bin/env python3
"""Independent rooted-triplet oracle and scoring-path parity test."""

from __future__ import annotations

import itertools
import os
import pathlib
import re
import subprocess
import sys
import tempfile

ROOT = pathlib.Path(__file__).resolve().parents[1]
GENES = ROOT / "test/input/test_5taxa.tre"
CANDIDATE = ROOT / "test/input/stelar_candidate_5taxa.tre"
UNROOTED = ROOT / "test/input/tc9_unrooted_simple.tre"
POLYTOMY_GENES = ROOT / "test/input/stelar_polytomy_5taxa.tre"
INCOMPLETE_POLYTOMY_GENES = ROOT / "test/input/stelar_polytomy_incomplete_6taxa.tre"
SIX_TAXON_CANDIDATE = ROOT / "test/input/stelar_candidate_6taxa.tre"
INCOMPLETE_GENES = ROOT / "test/input/test_incomplete.tre"
CHILD_ONLY_GENES = ROOT / "test/input/child_only_dedup_incomplete.tre"
CHILD_ONLY_CANDIDATE = ROOT / "test/input/child_only_candidate_7taxa.tre"


class Node:
    def __init__(self, children=None, name=None):
        self.children = children or []
        self.name = name


def parse_newick(text: str) -> Node:
    stack: list[Node | str] = []
    token = ""
    for char in text.strip().rstrip(";"):
        if char == "(":
            stack.append(char)
        elif char in ",)":
            if token.strip():
                stack.append(Node(name=token.strip()))
                token = ""
            if char == ")":
                children = []
                while stack and stack[-1] != "(":
                    children.append(stack.pop())
                if not stack:
                    raise ValueError("unbalanced Newick")
                stack.pop()
                stack.append(Node(children=list(reversed(children))))
        else:
            token += char
    if len(stack) != 1 or not isinstance(stack[0], Node):
        raise ValueError("malformed Newick")
    return stack[0]


def leaves(node: Node) -> set[str]:
    if node.name is not None:
        return {node.name}
    return set().union(*(leaves(child) for child in node.children))


def rooted_pair(node: Node, triple: tuple[str, str, str]):
    if node.name is not None:
        return None
    child_sets = [leaves(child) for child in node.children]
    locations = [next(i for i, group in enumerate(child_sets) if taxon in group)
                 for taxon in triple]
    if locations[0] == locations[1] != locations[2]:
        return tuple(sorted(triple[:2]))
    if locations[0] == locations[2] != locations[1]:
        return tuple(sorted((triple[0], triple[2])))
    if locations[1] == locations[2] != locations[0]:
        return tuple(sorted(triple[1:]))
    if locations[0] == locations[1] == locations[2]:
        return rooted_pair(node.children[locations[0]], triple)
    return None


def oracle_score(species: Node, genes: list[Node]) -> int:
    taxa = sorted(leaves(species))
    total = 0
    for gene in genes:
        present = leaves(gene)
        for triple in itertools.combinations(taxa, 3):
            if set(triple) <= present:
                species_pair = rooted_pair(species, triple)
                gene_pair = rooted_pair(gene, triple)
                total += gene_pair is not None and species_pair == gene_pair
    return total


def stelar_score(method: str, genes_path: pathlib.Path = GENES,
                 candidate_path: pathlib.Path = CANDIDATE,
                 numeric: str = "long") -> int:
    command = [
        "java", "-cp", str(ROOT / "build"), "stelarx.Main",
        "-i", str(genes_path), "--score-species-tree", str(candidate_path),
        "--cpu", "--intersection-method", method, "-q",
    ]
    env = os.environ.copy()
    if numeric != "long":
        env["STELARX_WEIGHT_FORCE_DOUBLE"] = "1"
        command += ["--large-n-score-type", numeric]
    run = subprocess.run(command, cwd=ROOT, env=env, text=True,
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if run.returncode:
        raise AssertionError(f"{method} failed:\n{run.stdout}")
    match = re.search(r"^TRIPLET_SCORE:\s*(\d+)(?:[.]0+)?\s*$", run.stdout, re.MULTILINE)
    if not match:
        raise AssertionError(f"{method} emitted no triplet score:\n{run.stdout}")
    return int(match.group(1))


def verify_root_preserving_completion() -> None:
    with tempfile.TemporaryDirectory(prefix="stelarx-completion-") as directory:
        completed_path = pathlib.Path(directory) / "completed.tre"
        species_path = pathlib.Path(directory) / "species.tre"
        command = [
            "java", "-cp", str(ROOT / "build"), "stelarx.Main",
            "-i", str(INCOMPLETE_GENES), "--cpu", "-q", "--search-space", "S2",
            "--dump-completed-gene-trees", str(completed_path), "-o", str(species_path),
        ]
        run = subprocess.run(command, cwd=ROOT, text=True,
                             stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        if run.returncode:
            raise AssertionError(f"root-preserving completion failed:\n{run.stdout}")

        originals = [parse_newick(line) for line in INCOMPLETE_GENES.read_text().splitlines()
                     if line.strip()]
        completed = [parse_newick(line) for line in completed_path.read_text().splitlines()
                     if line.strip()]
        assert len(completed) == len(originals)
        all_taxa = set().union(*(leaves(tree) for tree in originals))
        for index, (before, after) in enumerate(zip(originals, completed)):
            assert leaves(after) == all_taxa, (index, leaves(after), all_taxa)
            for triple in itertools.combinations(sorted(leaves(before)), 3):
                assert rooted_pair(before, triple) == rooted_pair(after, triple), (index, triple)


def main() -> int:
    if os.environ.get("STELARX_SKIP_BUILD") != "1":
        subprocess.run([str(ROOT / "scripts" / "build.sh")], cwd=ROOT, check=True,
                       stdout=subprocess.DEVNULL)
    genes = [parse_newick(line) for line in GENES.read_text().splitlines() if line.strip()]
    candidate = parse_newick(CANDIDATE.read_text())
    expected = oracle_score(candidate, genes)
    assert expected == 21, expected
    scores = {method: stelar_score(method) for method in ("I1", "I2", "I3", "I4")}
    assert set(scores.values()) == {expected}, scores

    polytomy_genes = [parse_newick(line) for line in POLYTOMY_GENES.read_text().splitlines()
                      if line.strip()]
    polytomy_expected = oracle_score(candidate, polytomy_genes)
    polytomy_scores = {
        method: stelar_score(method, POLYTOMY_GENES)
        for method in ("I1", "I2", "I3", "I4")
    }
    assert set(polytomy_scores.values()) == {polytomy_expected}, polytomy_scores

    incomplete_polytomy_genes = [
        parse_newick(line) for line in INCOMPLETE_POLYTOMY_GENES.read_text().splitlines()
        if line.strip()
    ]
    six_taxon_candidate = parse_newick(SIX_TAXON_CANDIDATE.read_text())
    incomplete_polytomy_expected = oracle_score(
        six_taxon_candidate, incomplete_polytomy_genes)
    incomplete_polytomy_scores = {
        method: stelar_score(method, INCOMPLETE_POLYTOMY_GENES, SIX_TAXON_CANDIDATE)
        for method in ("I1", "I2", "I3", "I4")
    }
    assert set(incomplete_polytomy_scores.values()) == {incomplete_polytomy_expected}, \
        incomplete_polytomy_scores

    # Targeted regression: identical binary and polytomous child collections occur
    # under different incomplete-tree outside sets. Child-only dedup must merge them
    # without changing the independently enumerated rooted-triplet score.
    child_only_genes = [
        parse_newick(line) for line in CHILD_ONLY_GENES.read_text().splitlines()
        if line.strip()
    ]
    child_only_candidate = parse_newick(CHILD_ONLY_CANDIDATE.read_text())
    child_only_expected = oracle_score(child_only_candidate, child_only_genes)
    child_only_scores = {
        method: stelar_score(method, CHILD_ONLY_GENES, CHILD_ONLY_CANDIDATE)
        for method in ("I1", "I2", "I3", "I4")
    }
    assert set(child_only_scores.values()) == {child_only_expected}, child_only_scores
    child_only_wide_scores = {
        (method, numeric): stelar_score(
            method, CHILD_ONLY_GENES, CHILD_ONLY_CANDIDATE, numeric)
        for numeric in ("int128", "double")
        for method in ("I1", "I2", "I3", "I4")
    }
    assert set(child_only_wide_scores.values()) == {child_only_expected}, \
        child_only_wide_scores

    verify_root_preserving_completion()

    with tempfile.TemporaryDirectory(prefix="stelarx-expected-crash-") as directory:
        reject_env = os.environ.copy()
        reject_env["STELARX_CRASH_DIR"] = str(pathlib.Path(directory) / "crash_logs")
        rejected = subprocess.run(
            ["java", "-cp", str(ROOT / "build"), "stelarx.Main",
             "-i", str(UNROOTED), "--cpu", "-q"],
            cwd=ROOT, env=reject_env, text=True,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    assert rejected.returncode != 0
    assert "never roots input trees arbitrarily" in rejected.stdout

    print(f"PASS: binary={expected} {scores}; polytomy={polytomy_expected} "
          f"{polytomy_scores}; incomplete-polytomy={incomplete_polytomy_expected} "
          f"{incomplete_polytomy_scores}; child-only={child_only_expected} "
          f"{child_only_scores} (+wide modes); completion root preserved; unrooted rejected")
    return 0


if __name__ == "__main__":
    sys.exit(main())
