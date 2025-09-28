#!/usr/bin/env python3
"""
run_paup_consensus.py

Usage:
  python3 run_paup_consensus.py -i <gene_trees> -o <output_prefix>
      [--greedy-start N] [--greedy-end N] [--greedy-step N] [--no-clean]
      [--paup-cmd "path/to/paup -n"]

Produces:
  <prefix>.strict.tree
  <prefix>.majority.tree
  <prefix>.greedy.<P>.tree  (for each greedy percent P)
  <prefix>.consensus.trees  (merged file containing all consensus trees)
"""
import argparse
import os
import re
import shlex
import subprocess
from pathlib import Path
from typing import List, Optional

DEFAULT_PAUP_CMD = os.environ.get("PAUP_CMD", "./paup -n")


def read_lines_strip(path: Path) -> List[str]:
    with path.open() as fh:
        return [ln.rstrip() for ln in fh if ln.strip()]


def convert_newick_to_nexus(tree_file: Path, nexus_file: Path) -> None:
    """Create a NEXUS 'Begin trees' file containing all trees from tree_file."""
    lines = read_lines_strip(tree_file)
    # ensure parent dir exists for nexus_file (should already exist but safe)
    if nexus_file.parent and not nexus_file.parent.exists():
        nexus_file.parent.mkdir(parents=True, exist_ok=True)
    with nexus_file.open("w") as out:
        out.write("#NEXUS\n\nBegin trees;\n")
        for i, ln in enumerate(lines):
            tree_ln = ln
            if "[&U]" not in tree_ln:
                tree_ln = "[&U] " + tree_ln
            out.write(f"tree TREE{i} = {tree_ln}\n")
        out.write("End;\n")


def get_taxa_list(tree_file: Path) -> List[str]:
    """Extract taxa names from Newick file (naively remove branch lengths and tokens)."""
    lines = read_lines_strip(tree_file)
    tokens = set()
    for ln in lines:
        # remove branch lengths like :0.123 or :1.2e-03
        ln2 = re.sub(r":-?\d+(?:\.\d+)?(?:[eE]-?\d+)?", "", ln)
        ln2 = re.sub(r"\[&U\]", "", ln2)
        # replace common separators with spaces
        ln2 = re.sub(r"[(),;:\[\]]", " ", ln2)
        for tok in ln2.split():
            if tok:
                tokens.add(tok)
    return sorted(tokens)


def write_paup_script(
    nexus_file: Path,
    paup_file: Path,
    strict_out: Path,
    major_out: Path,
    greedy_base: Path,
    greedy_start: int,
    greedy_end: int,
    greedy_step: int,
    taxa_list: List[str],
) -> None:
    """Write the PAUP script that will be fed to 'paup -n <script>'."""
    num_taxa = len(taxa_list)
    taxlabels = " ".join(taxa_list)
    greedy_percents = list(range(greedy_start, greedy_end + 1, greedy_step))

    # ensure paup_file parent exists
    if paup_file.parent and not paup_file.parent.exists():
        paup_file.parent.mkdir(parents=True, exist_ok=True)

    with paup_file.open("w") as out:
        out.write("#NEXUS\n\n")
        out.write("begin taxa;\n")
        out.write(f"    dimensions ntax={num_taxa};\n")
        out.write(f"    taxlabels {taxlabels} ;\n")
        out.write("end;\n\n")
        out.write("begin paup;\n")
        out.write("    set autoclose = yes warntree = no warnreset = no notifybeep = no monitor = yes taxlabels = full;\n")
        out.write("    set criterion = parsimony;\n")
        out.write("    set increase = auto;\n")
        out.write(f"    gettrees file = {nexus_file.name} allblocks = yes warntree = no unrooted = yes;\n")
        out.write(f"    contree all / strict = yes treefile = {strict_out.name} replace;\n")
        out.write(f"    contree all / strict = no majrule = yes treefile = {major_out.name} replace;\n")
        for p in greedy_percents:
            out.write(f"    contree all / strict = no majrule = yes percent={p} treefile = {greedy_base.name}.{p} replace;\n")
        out.write(f"    contree all / strict = no majrule = yes le50 = yes treefile = {greedy_base.name}.le50 replace;\n")
        out.write("end;\n\nquit warntsave = no;\n")


def run_paup_script(paup_cmd: str, paup_script: Path, stdout_log: Path, stderr_log: Path) -> int:
    """Run PAUP command and write stdout/stderr to provided logs. Return exit code."""
    # Make sure log parent exists
    if stdout_log.parent and not stdout_log.parent.exists():
        stdout_log.parent.mkdir(parents=True, exist_ok=True)
    if stderr_log.parent and not stderr_log.parent.exists():
        stderr_log.parent.mkdir(parents=True, exist_ok=True)

    cmd_list = shlex.split(paup_cmd) + [str(paup_script)]
    with stdout_log.open("wb") as outfh, stderr_log.open("wb") as errfh:
        try:
            rc = subprocess.call(cmd_list, stdout=outfh, stderr=errfh)
        except FileNotFoundError:
            errfh.write(b"PAUP executable not found. Check PAUP_CMD or your PATH.\n")
            rc = 127
    return rc


def find_paup_output_file(base: Path) -> Optional[Path]:
    """Try plausible filenames PAUP might produce for a treefile base."""
    candidates = [
        base,
        Path(str(base) + ".tre"),
        base.with_suffix(".tre"),
        base.with_suffix(".nex"),
        base.with_suffix(".tre.nex"),
        base.with_suffix(".tree"),
    ]
    for p in candidates:
        if p.exists():
            return p
    # fallback: try literal names (base + ".tre", etc.)
    for ext in ("", ".tre", ".nex", ".tre.nex", ".tree"):
        p = Path(str(base) + ext)
        if p.exists():
            return p
    return None


def extract_newick_from_paup_output(paup_file: Path) -> Optional[str]:
    """Attempt to extract the first Newick string from a PAUP output file content."""
    try:
        content = paup_file.read_text()
    except Exception:
        return None
    # look for tree ... = [&U] (something; )
    m = re.search(r"tree\s+[^=]+=.*?\[&U\]\s*(.*?;)", content, flags=re.S | re.I)
    if m:
        return m.group(1).strip()
    # fallback: first parenthesized group until semicolon
    m2 = re.search(r"(\([^;]+;)", content, flags=re.S)
    if m2:
        return m2.group(1).strip()
    return None


def convert_paup_output_to_newick(src_base: Path, out_newick: Path) -> bool:
    """Find PAUP output file for 'src_base' and write its Newick to out_newick.
    Returns True if successful."""
    found = find_paup_output_file(src_base)
    if not found:
        # nothing to convert
        return False
    newick = extract_newick_from_paup_output(found)
    if not newick:
        return False
    # ensure parent dir for out_newick exists
    if out_newick.parent and not out_newick.parent.exists():
        out_newick.parent.mkdir(parents=True, exist_ok=True)
    out_newick.write_text(newick + "\n")
    return True


def merge_consensus_trees(tree_files: List[str], output_file: Path) -> bool:
    """
    Merge all consensus tree files into a single file.
    Splits content by semicolons and newlines, then merges with semicolons and newlines.
    
    Args:
        tree_files (List[str]): List of tree file paths to merge
        output_file (Path): Output file path for merged trees
        
    Returns:
        bool: True if successful, False otherwise
    """
    if not tree_files:
        print("Warning: No tree files to merge")
        return False
    
    print(f"Merging {len(tree_files)} consensus trees into {output_file}")
    
    merged_trees = []
    
    for tree_file in tree_files:
        tree_path = Path(tree_file)
        if not tree_path.exists():
            print(f"Warning: Tree file {tree_file} does not exist, skipping")
            continue
            
        print(f"Processing: {tree_path.name}")
        try:
            with open(tree_path, 'r') as f:
                content = f.read().strip()
                
            # Split by semicolon and newline, filter out empty strings
            trees = [tree.strip() for tree in content.replace('\n', ';').split(';') if tree.strip()]
            merged_trees.extend(trees)
            
        except Exception as e:
            print(f"Error reading {tree_file}: {e}")
            continue
    
    if not merged_trees:
        print("Warning: No trees found to merge")
        return False
    
    # Write merged trees to output file
    try:
        # Ensure parent directory exists
        if output_file.parent and not output_file.parent.exists():
            output_file.parent.mkdir(parents=True, exist_ok=True)
            
        with open(output_file, 'w') as f:
            for tree in merged_trees:
                f.write(tree + ';\n')
        
        print(f"Successfully merged {len(merged_trees)} trees into '{output_file}'")
        return True
        
    except Exception as e:
        print(f"Error writing to output file '{output_file}': {e}")
        return False


def main():
    ap = argparse.ArgumentParser(description="Run PAUP consensus (strict, majority, greedy ranges).")
    ap.add_argument("-i", "--input", required=True, help="Gene trees file (one Newick per line).")
    ap.add_argument("-o", "--output", required=True, help="Output prefix.")
    ap.add_argument("--greedy-start", type=int, default=5, help="Greedy start percent (default 5).")
    ap.add_argument("--greedy-end", type=int, default=80, help="Greedy end percent (default 80).")
    ap.add_argument("--greedy-step", type=int, default=5, help="Greedy step (default 5).")
    ap.add_argument("--no-clean", action="store_true", help="Keep all generated files (by default only .consensus.trees is kept)")
    ap.add_argument("--paup-cmd", default=DEFAULT_PAUP_CMD, help="PAUP command (default from PAUP_CMD env or './paup -n').")
    args = ap.parse_args()

    gt = Path(args.input)
    if not gt.exists():
        raise SystemExit(f"Input file not found: {gt}")

    out_prefix = Path(args.output)

    # ensure parent directories for the prefix exist (so a/b/c.nexus can be created)
    if out_prefix.parent and not out_prefix.parent.exists():
        try:
            out_prefix.parent.mkdir(parents=True, exist_ok=True)
        except Exception as e:
            raise SystemExit(f"Failed to create parent directories for output prefix '{out_prefix}': {e}")

    nexus_file = out_prefix.with_suffix(".nexus")
    paup_script = out_prefix.with_suffix(".paup")
    strict_base = out_prefix.with_suffix(".strict")
    major_base = out_prefix.with_suffix(".majority")
    greedy_base = out_prefix.with_suffix(".greedy")

    # sanity checks for greedy range
    if args.greedy_start < 0 or args.greedy_end < args.greedy_start or args.greedy_step <= 0:
        raise SystemExit("Invalid greedy start/end/step parameters.")

    # create nexus file
    convert_newick_to_nexus(gt, nexus_file)

    # taxa list
    taxa = get_taxa_list(gt)

    # write paup script
    write_paup_script(
        nexus_file=nexus_file,
        paup_file=paup_script,
        strict_out=strict_base,
        major_out=major_base,
        greedy_base=greedy_base,
        greedy_start=args.greedy_start,
        greedy_end=args.greedy_end,
        greedy_step=args.greedy_step,
        taxa_list=taxa,
    )

    stdout_log = paup_script.with_suffix(paup_script.suffix + ".stdout")
    stderr_log = paup_script.with_suffix(paup_script.suffix + ".stderr")

    # run paup
    print("Running PAUP (capturing stdout/stderr)...")
    rc = run_paup_script(args.paup_cmd, paup_script, stdout_log, stderr_log)
    if rc != 0:
        print(f"PAUP returned non-zero exit code {rc}. See {stdout_log} and {stderr_log} for details.")

    # convert strict and majority
    final_trees = []
    if convert_paup_output_to_newick(strict_base, out_prefix.with_suffix(".strict.tree")):
        final_trees.append(str(out_prefix.with_suffix(".strict.tree")))
    else:
        print("Warning: strict consensus not found or not created by PAUP.")

    if convert_paup_output_to_newick(major_base, out_prefix.with_suffix(".majority.tree")):
        final_trees.append(str(out_prefix.with_suffix(".majority.tree")))
    else:
        print("Warning: majority consensus not found or not created by PAUP.")

    # le50 consensus (majority rule with le50 = yes)
    le50_base = Path(str(greedy_base) + ".le50")
    if convert_paup_output_to_newick(le50_base, out_prefix.with_suffix(".le50.tree")):
        final_trees.append(str(out_prefix.with_suffix(".le50.tree")))
    else:
        print("Warning: le50 consensus not found or not created by PAUP.")

    # greedy percent outputs
    for p in range(args.greedy_start, args.greedy_end + 1, args.greedy_step):
        base = Path(str(greedy_base) + f".{p}")
        outtree = Path(str(out_prefix) + f".greedy.{p}.tree")
        ok = convert_paup_output_to_newick(base, outtree)
        if ok:
            final_trees.append(str(outtree))
        else:
            print(f"NOTE: greedy {p}% produced no tree file (or it used an unexpected name).")

    # merge all consensus trees into a single file
    consensus_output = out_prefix.with_suffix(".consensus.trees")
    merge_success = merge_consensus_trees(final_trees, consensus_output)
    if merge_success:
        final_trees.append(str(consensus_output))

    # cleanup files: by default keep only .consensus.trees, unless --no-clean is specified
    if not args.no_clean:
        # Remove all intermediate files
        for p in [nexus_file, paup_script, strict_base, major_base]:
            if p.exists():
                try:
                    p.unlink()
                except Exception:
                    pass
        
        # Remove individual tree files (keep only the merged .consensus.trees)
        files_to_remove = []
        files_to_remove.append(out_prefix.with_suffix(".strict.tree"))
        files_to_remove.append(out_prefix.with_suffix(".majority.tree"))
        files_to_remove.append(out_prefix.with_suffix(".le50.tree"))
        
        # Remove greedy tree files
        for p in range(args.greedy_start, args.greedy_end + 1, args.greedy_step):
            files_to_remove.append(Path(str(out_prefix) + f".greedy.{p}.tree"))
        
        # Remove all the individual tree files
        for tree_file in files_to_remove:
            if tree_file.exists():
                try:
                    tree_file.unlink()
                except Exception:
                    pass
        
        # Remove le50 base file and greedy base files (PAUP intermediate outputs)
        le50_candidate = Path(str(greedy_base) + ".le50")
        for ext in ("", ".tre", ".nex", ".tre.nex"):
            f = Path(str(le50_candidate) + ext)
            if f.exists():
                try:
                    f.unlink()
                except Exception:
                    pass
        
        for p in range(args.greedy_start, args.greedy_end + 1, args.greedy_step):
            candidate = Path(str(greedy_base) + f".{p}")
            for ext in ("", ".tre", ".nex", ".tre.nex"):
                f = Path(str(candidate) + ext)
                if f.exists():
                    try:
                        f.unlink()
                    except Exception:
                        pass

        # Keep stderr log for debugging, but remove stdout
        if stdout_log.exists():
            try:
                stdout_log.unlink()
            except Exception:
                pass

    if args.no_clean:
        print("Wrote trees:")
        for t in final_trees:
            print(t)
    else:
        print("Wrote trees:")
        print(str(consensus_output))
        print("(Individual tree files cleaned up. Use --no-clean to keep all files.)")
    print("Done.")


if __name__ == "__main__":
    main()
