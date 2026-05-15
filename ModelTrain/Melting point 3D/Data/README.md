# Data Files

This folder contains the train/test CSV datasets used for melting point (3D descriptors) model training.

| File | Rows | Size | Status |
|---|---|---|---|
| `train_set.csv` | 3 589 | 83 MB | Split → `train_set_1.csv` … `train_set_4.csv` |
| `train_chon.csv` | 1 562 | 37 MB | Split → `train_chon_1.csv`, `train_chon_2.csv` |
| `train_cho.csv` | 1 231 | 28 MB | Split → `train_cho_1.csv`, `train_cho_2.csv` |
| `test_set.csv` | 897 | 21 MB | No split needed |
| `test_chon.csv` | 393 | 9 MB | No split needed |
| `test_cho.csv` | 297 | 7 MB | No split needed |

Large files are split into numbered parts (each < 24 MB) to comply with GitHub's 25 MB file size limit.

## Merging Split Files

From this directory, run:

```bash
python merge_data.py
```

This will automatically detect all numbered parts, merge each group back into the original file, and delete the parts once merging is complete.

**Example output:**

```
Merging 4 parts → train_set.csv
  merged train_set_1.csv
  merged train_set_2.csv
  merged train_set_3.csv
  merged train_set_4.csv
  → train_set.csv (82.8 MB, 3589 data rows)
  deleted train_set_1.csv
  ...
```

After running, only the original CSV files will remain in this folder.
