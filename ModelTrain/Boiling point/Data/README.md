# Data Files

This folder contains the train/test CSV datasets used for boiling point model training.

| File | Rows | Size | Status |
|---|---|---|---|
| `train_set.csv` | 3 764 | 59 MB | Split → `train_set_1.csv`, `train_set_2.csv`, `train_set_3.csv` |
| `train_cho.csv` | 2 303 | 37 MB | Split → `train_cho_1.csv`, `train_cho_2.csv` |
| `train_chon.csv` | 576 | 9 MB | No split needed |
| `test_set.csv` | 941 | 15 MB | No split needed |
| `test_cho.csv` | 542 | 9 MB | No split needed |
| `test_chon.csv` | 161 | 3 MB | No split needed |

Large files are split into numbered parts (each < 24 MB) to comply with GitHub's 25 MB file size limit.

## Merging Split Files

From this directory, run:

```bash
python merge_data.py
```

This will automatically detect all numbered parts, merge each group back into the original file, and delete the parts once merging is complete.

**Example output:**

```
Merging 2 parts → train_cho.csv
  merged train_cho_1.csv
  merged train_cho_2.csv
  → train_cho.csv (36.5 MB, 2303 data rows)
  deleted train_cho_1.csv
  deleted train_cho_2.csv
```

After running, only the original CSV files will remain in this folder.
