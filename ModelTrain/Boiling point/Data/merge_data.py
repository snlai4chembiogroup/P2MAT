#!/usr/bin/env python
######################################################
# Developer : Methun Kamruzzaman
# Date      : May 15, 2026
# Purpose   : Merges numbered CSV split parts back into original files and
#             removes the split parts.
######################################################
"""
merge_data.py
-------------
Reconstructs original CSV files from their numbered split parts, then deletes
the parts once the merged file is successfully written.

For each base name that has numbered parts (e.g. train_cho_1.csv, train_cho_2.csv),
this script merges them back in order into the original file (train_cho.csv),
skipping the duplicate header rows from parts 2 onward.

Usage:
    python merge_data.py

The merged file overwrites any existing file with the original name.
Split part files are removed after a successful merge.
"""

import os
import re
from collections import defaultdict

DATA_DIR = os.path.dirname(os.path.abspath(__file__))


def merge_parts(base: str, parts: list[str]) -> None:
    out_path = os.path.join(DATA_DIR, f"{base}.csv")
    parts_sorted = sorted(parts, key=lambda f: int(re.search(r"_(\d+)\.csv$", f).group(1)))

    print(f"Merging {len(parts_sorted)} parts → {os.path.basename(out_path)}")
    total_rows = 0

    with open(out_path, "w", encoding="utf-8") as out:
        for i, part_name in enumerate(parts_sorted):
            part_path = os.path.join(DATA_DIR, part_name)
            with open(part_path, "r", encoding="utf-8") as fh:
                header = fh.readline()
                if i == 0:
                    out.write(header)   # write header once from first part
                for line in fh:
                    out.write(line)
                    total_rows += 1
            print(f"  merged {part_name}")

    size_mb = os.path.getsize(out_path) / 1e6
    print(f"  → {os.path.basename(out_path)} ({size_mb:.1f} MB, {total_rows} data rows)")

    for part_name in parts_sorted:
        os.remove(os.path.join(DATA_DIR, part_name))
        print(f"  deleted {part_name}")
    print()


if __name__ == "__main__":
    # Group split files by their base name.
    groups: dict[str, list[str]] = defaultdict(list)
    for fname in os.listdir(DATA_DIR):
        m = re.search(r"^(.+)_(\d+)\.csv$", fname)
        if m:
            groups[m.group(1)].append(fname)

    if not groups:
        print("No split files found. Nothing to merge.")
    else:
        print(f"Found {len(groups)} group(s) to merge:\n")
        for base, parts in sorted(groups.items()):
            merge_parts(base, parts)
        print("Done.")
