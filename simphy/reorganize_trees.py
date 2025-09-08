#!/usr/bin/env python3
# reorganize_trees.py
from pathlib import Path
import re
import sys
import shutil

def natural_key(p: Path):
    """Extract numeric part from path for sorting"""
    m = re.search(r'(\d+)', p.stem)
    return int(m.group(1)) if m else 0

def reorganize_trees(root_dir: Path):
    """Reorganize tree files from numeric subdirs to R-prefixed subdirs"""
    if not root_dir.is_dir():
        print(f"Error: {root_dir} is not a directory")
        return False
    
    # Find all numeric subdirectories
    numeric_subdirs = [d for d in root_dir.iterdir() if d.is_dir() and d.name.isdigit()]
    
    if not numeric_subdirs:
        print(f"No numeric subdirectories found in {root_dir}")
        return False
    
    # Sort by numeric value
    numeric_subdirs.sort(key=natural_key)
    
    processed_count = 0
    successfully_processed_dirs = []
    
    for subdir in numeric_subdirs:
        # Create corresponding R-prefixed directory
        r_dir_name = f"R{subdir.name}"
        r_dir_path = root_dir / r_dir_name
        
        # Files we need to copy
        gene_tree_file = subdir / "all_gt.tre"
        species_tree_file = subdir / "s_tree.trees"
        
        # Check if required files exist
        missing_files = []
        if not gene_tree_file.exists():
            missing_files.append("all_gt.tre")
        if not species_tree_file.exists():
            missing_files.append("s_tree.trees")
        
        if missing_files:
            print(f"Skipping {subdir.name}: missing files {missing_files}")
            continue
        
        # Create R directory if it doesn't exist
        r_dir_path.mkdir(exist_ok=True)
        
        # Copy and clean the gene tree file (remove _0_0 patterns)
        gene_content = gene_tree_file.read_text(encoding="utf-8")
        cleaned_gene_content = gene_content.replace("_0_0", "")
        (r_dir_path / "all_gt.tre").write_text(cleaned_gene_content, encoding="utf-8")
        
        # Copy the species tree file as-is
        shutil.copy2(species_tree_file, r_dir_path / "s_tree.trees")
        
        print(f"Created {r_dir_name}/ with gene trees and species tree from {subdir.name}/")
        processed_count += 1
        successfully_processed_dirs.append(subdir)
    
    if processed_count > 0:
        print(f"\nSuccessfully processed {processed_count} directories")
        
        # Clean up original numeric directories
        print("Cleaning up original directories...")
        for subdir in successfully_processed_dirs:
            try:
                shutil.rmtree(subdir)
                print(f"Deleted {subdir.name}/")
            except Exception as e:
                print(f"Warning: Could not delete {subdir.name}/: {e}")
        
        print(f"Cleanup complete!")
        return True
    else:
        print("No directories were processed")
        return False

def main():
    root = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("out")
    
    if not root.exists():
        print(f"Error: Directory {root} does not exist")
        return
    
    reorganize_trees(root)

if __name__ == "__main__":
    main()
