#!/usr/bin/env python3
import importlib.util
from pathlib import Path
import sys

SCRIPT_PATH = Path(__file__).with_name("analyze-dataset.py")
spec = importlib.util.spec_from_file_location("analyze_dataset_impl", SCRIPT_PATH)
module = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(module)

if __name__ == "__main__":
    sys.exit(module.main())
