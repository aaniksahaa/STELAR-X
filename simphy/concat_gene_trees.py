#!/usr/bin/env python3
# concat_gene_trees.py
from pathlib import Path
import re
import sys

def natural_key(p: Path):
    m = re.search(r'(\d+)', p.stem)
    return int(m.group(1)) if m else 0

def collect_into(dirpath: Path):
    files = sorted(dirpath.glob("g_trees*.trees"), key=natural_key)
    if not files:
        return False
    out = dirpath / "all_gt.tre"
    with out.open("w", encoding="utf-8") as w:
        for f in files:
            text = f.read_text(encoding="utf-8")
            text = text.replace("[&R]", "")  # drop any prefix if present
            text = text.replace("\r", "").replace("\n", "").strip()  # single-line
            if not text.endswith(";"):
                text += ";"
            w.write(text + "\n")
    print(f"Wrote {out} ({len(files)} trees)")
    return True

def main():
    root = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("out")
    if root.is_dir() and any(d.is_dir() for d in root.iterdir()):
        # If given 'out', process each numeric subdir
        any_done = False
        for sub in sorted([p for p in root.iterdir() if p.is_dir()], key=natural_key):
            any_done |= collect_into(sub)
        if not any_done:
            print("No g_trees*.trees files found.")
    else:
        # Process a single replicate directory
        if not collect_into(root):
            print(f"No g_trees*.trees files found in {root}")

if __name__ == "__main__":
    main()
