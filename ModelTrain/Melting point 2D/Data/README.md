# Data Files

This folder contains the train/test CSV datasets used for melting point (2D descriptors) model training.

| File | Rows | Size | Status |
|---|---|---|---|
| `train_set.csv` | 15 428 | 250 MB | Split → `train_set_1.csv` … `train_set_11.csv` |
| `train_chon.csv` | 6 431 | 106 MB | Split → `train_chon_1.csv` … `train_chon_5.csv` |
| `train_cho.csv` | 5 583 | 91 MB | Split → `train_cho_1.csv` … `train_cho_4.csv` |
| `test_chon.csv` | 1 647 | 27 MB | Split → `test_chon_1.csv`, `test_chon_2.csv` |
| `test_set.csv` | 3 860 | 63 MB | Split → `test_set_1.csv`, `test_set_2.csv`, `test_set_3.csv` |
| `test_cho.csv` | 1 374 | 22 MB | No split needed |

Large files are split into numbered parts (each < 24 MB) to comply with GitHub's 25 MB file size limit.

## Merging Split Files

From this directory, run:

```bash
python merge_data.py
```

This will automatically detect all numbered parts, merge each group back into the original file, and delete the parts once merging is complete.

**Example output:**

```
Merging 11 parts → train_set.csv
  merged train_set_1.csv
  ...
  merged train_set_11.csv
  → train_set.csv (249.8 MB, 15428 data rows)
  deleted train_set_1.csv
  ...
```

After running, only the original CSV files will remain in this folder.
