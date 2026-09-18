#!/usr/bin/env python3
"""Regression tests for quoted/unquoted Newick taxon identity."""

from __future__ import annotations

import pathlib
import subprocess
import sys
import tempfile

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import rf  # noqa: E402


UNQUOTED = "((Acanthisitta_chloris,B),(C,D));\n"
QUOTED = "(('Acanthisitta_chloris',B),(C,D));\n"
MISMATCHED = "((Acanthisitta_chloris,B),(C,E));\n"


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="stelarx-labels-") as directory:
        work = pathlib.Path(directory)
        plain = work / "plain.tre"
        quoted = work / "quoted.tre"
        mismatched = work / "mismatched.tre"
        rooted = work / "rooted.tre"
        cleaned = work / "cleaned.tre"
        plain.write_text(UNQUOTED)
        quoted.write_text(QUOTED)
        mismatched.write_text(MISMATCHED)

        assert rf.compare_trees(str(plain), str(quoted)) == 0.0

        try:
            rf.compare_trees(str(plain), str(mismatched))
        except ValueError as exc:
            message = str(exc)
            assert "different taxon sets" in message
            assert "D" in message and "E" in message
        else:
            raise AssertionError("RF accepted genuinely different taxon sets")

        subprocess.run(
            [sys.executable, str(ROOT / "scripts" / "root_by_outgroups.py"),
             "-i", str(plain), "-o", str(rooted),
             "-og", "Acanthisitta_chloris", "--num-workers", "1", "-q"],
            cwd=ROOT, check=True,
        )
        rooted_text = rooted.read_text()
        assert "Acanthisitta_chloris" in rooted_text
        assert "'Acanthisitta_chloris'" not in rooted_text

        subprocess.run(
            [sys.executable, str(ROOT / "scripts" / "clean.py"),
             "-i", str(quoted), "-o", str(cleaned),
             "--num-workers", "1", "--deterministic"],
            cwd=ROOT, check=True, stdout=subprocess.DEVNULL,
        )
        cleaned_text = cleaned.read_text()
        assert "Acanthisitta_chloris" in cleaned_text
        assert "'Acanthisitta_chloris'" not in cleaned_text

    print("Quoted/unquoted taxon normalization and RF validation: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
