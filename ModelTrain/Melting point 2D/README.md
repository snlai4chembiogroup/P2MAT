# QSAR Melting Point Prediction Pipeline

A modular machine-learning pipeline for predicting molecular melting points (MP\_K, Kelvin)
from molecular descriptors using 6 regression models, 4 feature-selection strategies,
and 3 data transformations — all optimised with Bayesian hyperparameter tuning
and 10-fold cross-validation.

---

## Project structure

```
Melting point 2D/
├── Data/
│   ├── train_set.csv          original training set (15 428 molecules)
│   ├── test_set.csv           original test set     ( 3 860 molecules)
│   ├── train_cho.csv          CHO subset (no N)     – generated  ( 5 583 train / 1 374 test)
│   ├── test_cho.csv
│   ├── train_chon.csv         CHON subset (has N)   – generated  ( 6 431 train / 1 647 test)
│   ├── test_chon.csv
│   ├── split_data.py          splits large CSVs into <24 MB parts for GitHub
│   └── merge_data.py          merges numbered parts back into original files
├── saved_models/              best model artefacts per experiment
│   └── {dataset}/{fs}/{transform}/{model}/pipeline.joblib
├── results/
│   ├── all_results.csv        full metric table (all 216 experiments)
│   ├── summary_test_mae.csv   pivot table – best test MAE per combo
│   ├── summary_train_mae.csv  pivot table – best train MAE per combo
│   └── stacking_results.csv   stacking ensemble results
├── src/
│   ├── config.py              paths, constants, hyperparameter search spaces
│   ├── data_loader.py         CSV loading, cleaning, feature/target split
│   ├── feature_selection.py   4 feature-selection strategies
│   ├── models.py              model factory (6 estimators)
│   ├── dnn_model.py           PyTorch DNN with MPS + sklearn wrapper
│   ├── trainer.py             BayesSearchCV experiment runner
│   ├── evaluator.py           MAE / R² / RMSE computation + results table
│   └── stacking.py            stacking ensemble over LightGBM + XGBoost + DNN
├── saved_models_cleaned/      retrained model artefacts (outlier-cleaned training sets)
│   └── {dataset}/{fs}/{transform}/{model}/pipeline.joblib
├── results_cleaned/
│   ├── comparison.csv         before vs. after outlier removal metric comparison
│   └── all_results_cleaned.csv
├── boruta.py                  Boruta feature selection implementation
├── split_datasets.py          generates CHO and CHON subsets
├── main.py                    experiment orchestrator (216 base experiments)
├── stack.py                   stacking ensemble entry point
├── outlier_detection.py       Williams-plot outlier detection + AD analysis
├── retrain_cleaned.py         retrain best models on response-outlier-cleaned data
├── shap_analysis.py           SHAP feature importance analysis
├── williams_plot.py           Williams plot visualisation
├── uncertainty_analysis.py    prediction uncertainty estimation
├── similarity_score.py        Tanimoto-based training set similarity scoring
├── pred_plot.py               predicted vs. actual scatter plots
├── dataset_similarity.py      inter-dataset Tanimoto similarity analysis
├── methodology_diagram.py     pipeline methodology figure for the manuscript
├── schematic_diagram.py       schematic figure for the manuscript
└── README.md
```

---

## Datasets

| Name  | Filter                        | Train  | Test  |
|-------|-------------------------------|--------|-------|
| full  | all molecules                 | 15 428 | 3 860 |
| cho   | nN == 0 (C, H, O only)        | 7 074  | ~1 770|
| chon  | nN > 0  (C, H, O + N)         | 8 354  | ~2 090|

Each dataset uses ~1 240 molecular descriptors (PaDEL / CDK) as features
and `MP_K` (melting point in Kelvin) as the target variable.

Dropped columns: `#PubChem`, `key`, `name`, `smiles`, `mpC`, `scaffold`  
(identifiers, duplicate target, or non-numeric strings)

> **Note:** Large CSV files in `Data/` are split into numbered parts for GitHub
> (e.g. `train_set_1.csv` … `train_set_11.csv`).  
> Run `python Data/merge_data.py` to reconstruct the originals before use.

---

## Experiment axes

### Models (6)

| Key         | Algorithm                  | Notes                                              |
|-------------|----------------------------|----------------------------------------------------|
| linear      | Ridge Regression           | L2-regularised linear model                        |
| elasticnet  | ElasticNet                 | L1+L2, replaces SVR — fast, handles multicollinearity |
| extratrees  | Extra Trees Regressor      | Randomised splits, replaces RF — 2-4× faster       |
| lgbm        | LightGBM                   | Histogram + leaf-wise, replaces GBR — 10-20× faster |
| xgb         | XGBoost                    | hist tree method, n\_jobs=-1                       |
| dnn         | Deep Neural Network        | PyTorch + MPS (Apple Silicon)                      |

DNN architecture: `Input → [Linear → BatchNorm1d → Activation → Dropout] × n_layers → Linear(1)`  
Loss: Huber loss · Optimiser: Adam · LR scheduler: ReduceLROnPlateau · Early stopping

**DNN two-stage training strategy** (to accelerate HPO):  
Training the DNN 200 times during a standard 10-fold × 20-iteration Bayesian search would require up to 40 000 epoch-runs. Instead, a two-stage approach is used:

| Stage | CV folds | Iterations | Max epochs | Patience | Purpose |
|-------|----------|-----------|------------|----------|---------|
| HPO search | 5 | 15 | 50 | 10 | Fast candidate evaluation (~3 750 epoch-runs) |
| Final refit | — | — | 200 | 20 | Full training of the best hyperparameters on all training data |

This delivers approximately **10× fewer epoch-runs** during HPO with no loss in final model quality, since the winner is always fully retrained from scratch after the search completes.

### Data transformations (3)

| Key      | Description                                      |
|----------|--------------------------------------------------|
| none     | Raw features, no scaling                         |
| standard | StandardScaler (mean 0, unit variance)           |
| quantile | QuantileTransformer (normal output distribution) |

### Feature selection methods (4)

| Key         | Description                                                                                       |
|-------------|---------------------------------------------------------------------------------------------------|
| whole       | All features (no selection)                                                                       |
| boruta      | Boruta with XGBoost estimator — selects all confirmed + tentative features                       |
| pca         | PCA with adaptive component selection (see detail below)                                          |
| correlation | Pearson correlation filter: select features correlated with target, remove inter-correlated ones  |

**PCA selection detail:**  
1. Fit full PCA on the (scaled) training data.  
2. Find the minimum number of PCs (`k`) that together explain ≥ 90 % of variance.  
3. If `k ≥ 10` → use `k` as-is and report the actual variance (e.g. 90.3 %).  
4. If `k < 10` → relax the variance constraint and force `k = 10`, then report the variance those 10 PCs actually explain (e.g. 93.1 %). A log message flags when relaxation occurs.  
5. If fewer than 10 features exist in total → cap at the number of available features (100 % variance).  

The final component count and explained variance are stored in `n_components_` and `actual_variance_` on the fitted selector.

**Correlation selection detail:**  
1. Keep features with |corr(feature, target)| ≥ 0.05  
2. Sort by descending |corr with target|  
3. Greedily discard features where |corr with any accepted feature| ≥ 0.85  
This avoids multicollinearity while preserving predictive signal.

---

## Hyperparameter optimisation

- **Algorithm**: Bayesian optimisation (`scikit-optimize` `BayesSearchCV`)
- **Scoring**: `neg_mean_absolute_error`
- **Random initial points**: 5 per experiment

| Model | CV folds | Bayesian iterations | Notes |
|-------|----------|-------------------|-------|
| linear, elasticnet, extratrees, lgbm, xgb | 10 | 20 | Standard budget |
| dnn | 5 | 15 | Reduced HPO budget; winner fully retrained (see DNN two-stage strategy above) |

- **Refit**: best hyperparameters are always refitted on the full training set

### Hyperparameter search spaces

| Model       | Parameters                                                                      |
|-------------|---------------------------------------------------------------------------------|
| linear      | alpha ∈ [1e-4, 1e3] (log-uniform)                                               |
| elasticnet  | alpha (log-uniform), l1\_ratio ∈ [0, 1]                                         |
| extratrees  | n\_estimators, max\_depth, min\_samples\_split/leaf, max\_features              |
| lgbm        | n\_estimators, learning\_rate, max\_depth, num\_leaves, subsample, colsample, reg |
| xgb         | n\_estimators, learning\_rate, max\_depth, subsample, colsample, reg            |
| dnn         | n\_layers, hidden\_size, dropout, lr, batch\_size, weight\_decay, act           |

---

## Evaluation metrics

Reported for both training and test phases:

| Metric | Formula                                 |
|--------|-----------------------------------------|
| MAE    | Mean Absolute Error (K)                 |
| R²     | Coefficient of Determination            |
| RMSE   | Root Mean Squared Error (K)             |

---

## Hardware optimisation (MacBook Pro M2 Max)

- All sklearn models use `n_jobs=-1` (all performance cores)
- XGBoost uses `tree_method="hist"` (optimised for modern CPUs)
- PyTorch DNN uses **MPS** (Metal Performance Shaders) when available, falls back to CUDA then CPU
- BayesSearchCV for DNN uses `n_jobs=1` (MPS doesn't support multi-process GPU)
- DNN HPO uses a reduced budget (5-fold, 15 iter, 50 epochs) then a full refit (200 epochs) — ~10× fewer epoch-runs during search

---

## Setup & usage

### 1. Environment

```bash
conda activate qsar
pip install scikit-optimize   # if not already installed
```

Required packages: `pandas`, `numpy`, `scikit-learn`, `xgboost`, `torch`, `scipy`,
`joblib`, `scikit-optimize`

### 2. Generate CHO / CHON subsets

```bash
conda run -n qsar python split_datasets.py
```

### 3. Run the full pipeline (216 experiments)

```bash
conda run -n qsar python main.py
```

### 4. Run a subset (e.g. quick test)

```bash
# Only full dataset, no transformation, whole features, ExtraTrees and XGB
conda run -n qsar python main.py \
    --datasets full \
    --transforms none \
    --fs whole \
    --models extratrees xgb
```

### 5. Run datasets sequentially across separate sessions

Each run automatically loads any existing `all_results.csv` and **appends** new
results — previous datasets are never overwritten. You can safely run one dataset
at a time across separate sessions:

```bash
# Session 1
conda run -n qsar python main.py --datasets full ...

# Session 2 (full results are preserved and cho results are appended)
conda run -n qsar python main.py --datasets cho ...

# Session 3 (full + cho results are preserved and chon results are appended)
conda run -n qsar python main.py --datasets chon ...
```

The summary tables (`summary_test_mae.csv`, `summary_train_mae.csv`) are always
rebuilt from the combined results of all completed sessions.

### 6. Automatic load-or-train logic (default behaviour)

Before running any experiment the pipeline applies the following priority check —
**no flags needed**:

```
results already in all_results.csv  →  SKIP  (nothing to do)
model saved, results entry missing  →  LOAD  → predict on train+test → append results
no saved model and no results       →  TRAIN → HPO → refit → save → append results
```

This means:
- Re-running `main.py` on an already-completed dataset is instant (all SKIPped).
- If a run was interrupted mid-way, only the missing combinations are LOADed or TRAINed; completed ones are untouched.
- Results missing from `all_results.csv` but whose model file exists are recovered automatically via LOAD — no retraining needed.

To force retraining even when a model and results already exist, use `--retrain`:

```bash
conda run -n qsar python main.py --datasets cho --retrain
```

### 7. Run the stacking ensemble

Stacking requires the base models (LightGBM, XGBoost, DNN) to already be trained
by `main.py` for the desired combinations.

```bash
# Run stacking over all 36 combinations (3 datasets × 3 transforms × 4 FS methods)
conda run -n qsar python stack.py

# Specific subset
conda run -n qsar python stack.py --datasets full cho --fs whole boruta

# Force re-run even if stacking results already exist
conda run -n qsar python stack.py --retrain
```

---

## Stacking ensemble

The stacking ensemble combines four model families —
**LightGBM**, **XGBoost**, **Extra Trees**, and **Deep Neural Network** — into a single
level-2 model using a Ridge regression meta-learner.

### Why these four?

LightGBM and XGBoost consistently outperform the linear baselines on molecular
descriptor data. Extra Trees adds diversity through its randomised split strategy,
which produces different error patterns than the boosting models. The DNN, while
weaker alone, captures complementary non-linear patterns that further improve
ensemble diversity.

### Architecture

```
                          ┌─────────────┐
                          │   X_train   │
                          └──────┬──────┘
            ┌─────────────┬──────┴──────┬─────────────┐
            ▼             ▼             ▼             ▼
       LightGBM        XGBoost     ExtraTrees        DNN
  (best_params HPO) (best_params) (best_params) (best_params HPO)
            │             │             │             │
       OOF preds     OOF preds     OOF preds     OOF preds
            └─────────────┴──────┬──────┴─────────────┘
                                 ▼
                      Ridge meta-learner
                      (trained on OOF matrix)
                                 │
                                 ▼
                          Final prediction
```

### Training procedure

1. **Load base-model artifacts** — for each (dataset, transform, fs_method) combo,
   load the three saved `pipeline.joblib` files produced by `main.py`.

2. **Out-of-fold (OOF) predictions** — to train the meta-learner without
   data leakage, fresh model instances are created using each base model's
   saved best hyperparameters and evaluated with 5-fold cross-validation:

   | Base model  | OOF epochs / config            |
   |-------------|-------------------------------|
   | LightGBM    | best_params from HPO (full)   |
   | XGBoost     | best_params from HPO (full)   |
   | Extra Trees | best_params from HPO (full)   |
   | DNN         | best_params, 50 epochs (fast) |

   The DNN uses a reduced epoch budget during OOF folds (same as the HPO
   stage) to keep the stacking run practical.

3. **Meta-learner training** — a Ridge regression (`alpha=1.0`) is fitted
   on the `[n_train × 4]` OOF prediction matrix:

   ```
   meta_train = [oof_lgbm | oof_xgb | oof_extratrees | oof_dnn]  shape: (n_train, 4)
   meta.fit(meta_train, y_train)
   ```

4. **Test prediction** — each base model's fully-trained saved artifact
   is used to generate test-set predictions. These are stacked and passed
   through the meta-learner:

   ```
   meta_test = [pred_lgbm | pred_xgb | pred_extratrees | pred_dnn]  shape: (n_test, 4)
   y_pred = meta.predict(meta_test)
   ```

5. **Evaluation** — train MAE is computed from OOF predictions (unbiased);
   test MAE uses the stacked test predictions.

6. **Artifacts saved** to `saved_models/{dataset}/{fs}/{transform}/stack/pipeline.joblib`:

   ```python
   {
       "base_models":    ["lgbm", "xgb", "extratrees", "dnn"],
       "meta_learner":   Ridge(...),        # fitted meta-learner
       "n_folds":        5,
       "meta_weights":   {"lgbm": ..., "xgb": ..., "extratrees": ..., "dnn": ...},
       "meta_intercept": float,
       "train_metrics":  {"mae": ..., "r2": ..., "rmse": ...},
       "test_metrics":   {"mae": ..., "r2": ..., "rmse": ...},
   }
   ```

### CLI options

| Flag           | Default | Description                                 |
|----------------|---------|---------------------------------------------|
| `--datasets`   | all     | `full`, `cho`, `chon`                       |
| `--transforms` | all     | `none`, `standard`, `quantile`              |
| `--fs`         | all     | `whole`, `boruta`, `pca`, `correlation`     |
| `--retrain`    | off     | Re-run even if stacking result already saved |

### Output

- **`results/stacking_results.csv`** — one row per (dataset, transform, fs_method)
  with columns: `dataset`, `transform`, `fs_method`, `model` (`stack`),
  `n_features`, `meta_weights`, `train_mae`, `train_r2`, `train_rmse`,
  `test_mae`, `test_r2`, `test_rmse`, `elapsed_sec`
- A summary pivot table is printed on completion showing Test MAE by
  dataset × feature-selection × transformation.

---

## Outlier detection and applicability domain

### Overview

After training the three best stacking ensembles, a Williams-plot-based outlier
detection step is applied to each (dataset, fs\_method, transform) combination
via `outlier_detection.py`. The goal is to (1) identify and characterise
problematic training compounds, (2) flag test-set compounds that lie outside the
model's applicability domain (AD), and (3) measure the impact of outlier removal
on test-set performance.

### Workflow

1. **Load predictions** — OOF (out-of-fold) training predictions and test
   predictions are loaded from the stacking pipeline cache.

2. **Compute standardised residuals** — raw residuals for both splits are
   scaled by the **training residual standard deviation** (consistent with
   standard AD practice so that both splits share the same reference scale):

   ```
   std_res_train = (y_true - y_pred_oof) / σ_train
   std_res_test  = (y_true - y_pred_test) / σ_train   # σ_train used for both
   ```

3. **Compute leverage** — the leverage $h_i$ of each compound is derived from
   the meta-feature matrix (the OOF prediction matrix used to train the Ridge
   meta-learner):

   ```
   H = X (XᵀX)⁻¹ Xᵀ    →    hᵢ = diagonal element
   ```

   The warning leverage threshold is:

   $$h^* = \frac{3(k+1)}{n}$$

   where $k$ is the number of base models (4) and $n$ is the number of
   training compounds.

4. **Classify training-set compounds** — each training compound falls into
   one of four categories:

   | Category | Condition | Action |
   |---|---|---|
   | Normal | $\|e_i\| \leq 3$ and $h_i \leq h^*$ | Kept |
   | Response outlier | $\|e_i\| > 3$ only | **Removed** in retrain |
   | Leverage outlier | $h_i > h^*$ only | Kept (no prediction error) |
   | Bad leverage point | $\|e_i\| > 3$ **and** $h_i > h^*$ | **Removed** in retrain |

   Pure leverage outliers are retained because they are structurally unusual
   but not poorly predicted — removing them is not statistically justified.

5. **Flag test-set compounds** — the same thresholds are applied to the test
   set. Flagged test compounds are **not removed**; they are used to report
   split performance metrics (all compounds vs. in-AD compounds only).

6. **Outputs**:
   - `results/outlier_williams_{dataset}.png` — annotated Williams plot
     (training and test panels, colour-coded by outlier category)
   - `results/outliers_smiles_{dataset}.csv` — flagged training compounds
     with SMILES, true/predicted MP, residual, standardised residual,
     leverage, and outlier type

### Running outlier detection

```bash
conda run -n qsar python outlier_detection.py
```

Must be run **after** `stack.py` has completed for the three best combinations,
since it reads the stacking OOF/test prediction cache.

### Retraining on cleaned data

`retrain_cleaned.py` uses the outlier CSV produced above to retrain the three
best stacking models after removing response outliers from the training set:

```bash
conda run -n qsar python retrain_cleaned.py
```

Cleaned model artefacts are written to `saved_models_cleaned/` and comparison
metrics to `results_cleaned/`. The original `saved_models/` directory is
**not modified**.

### Summary statistics (Williams plot)

| Column | Description |
|---|---|
| $n_\text{train}$ | Training set size |
| $n_\text{test}$ | Test set size |
| $h^*$ | Leverage warning threshold $= 3(k+1)/n$ |
| $\|e\|>3$ | Training response outliers ($\|\text{std\_res}\|>3$) |
| $h>h^*$ | Training leverage outliers (pure, no large residual) |
| Both | Training bad leverage points ($\|e\|>3$ **and** $h>h^*$) |
| Ts\_flag | Test compounds outside AD ($\|e\|>3$ or $h>h^*$) |
| MAE$_\text{all}$ | Test MAE on all test compounds (K) |
| MAE$_\text{clean}$ | Test MAE excluding flagged test compounds (K) |
| $R^2_\text{all}$ / $R^2_\text{clean}$ | Corresponding $R^2$ values |

---

## Output files

| File                             | Contents                                          |
|----------------------------------|---------------------------------------------------|
| `results/all_results.csv`        | All completed experiments, all metrics (appended across runs) |
| `results/summary_test_mae.csv`   | Pivot: best test MAE per (dataset, transform, model) |
| `results/summary_train_mae.csv`  | Same for training phase                           |
| `results/stacking_results.csv`   | Stacking ensemble results (one row per combo)     |
| `saved_models/.../pipeline.joblib` | Scaler + selector + trained model per experiment |
| `saved_models/.../dnn_weights.pt`  | DNN weights (PyTorch state\_dict)               |
| `saved_models/.../stack/pipeline.joblib` | Meta-learner + weights per stacking combo |
| `saved_models_cleaned/.../pipeline.joblib` | Same artefacts retrained after outlier removal |
| `results/outlier_williams_{dataset}.png` | Annotated Williams plot (train + test panels) |
| `results/outliers_smiles_{dataset}.csv`  | Flagged training compounds with SMILES + diagnostics |
| `results_cleaned/comparison.csv`   | Before vs. after outlier removal metric comparison |

---

## Results table format

```
Dataset | Data Transformation | Ridge | ElasticNet | ExtraTrees | LightGBM | XGBoost | DNN
--------+---------------------+-------+------------+------------+----------+---------+----
FULL    | No Transformation   |  ...  |    ...     |    ...     |   ...    |   ...   | ...
FULL    | Standard Scaler     |  ...  |    ...     |    ...     |   ...    |   ...   | ...
...
```

Each cell shows the **best test MAE (K)** across all feature-selection methods
for that (dataset, transformation, model) combination.

---

## Module overview

| Module                   | Responsibility                                       |
|--------------------------|------------------------------------------------------|
| `src/config.py`          | All constants, paths, and hyperparameter spaces      |
| `src/data_loader.py`     | Load CSVs, impute, drop non-numeric / constant cols  |
| `src/feature_selection.py` | WholeSelector, BorutaSelector, PCASelector, CorrelationSelector |
| `src/models.py`          | `get_model(name)` factory                            |
| `src/dnn_model.py`       | `_QSARNet` (PyTorch), `DNNRegressor` (sklearn wrap)  |
| `src/trainer.py`         | `run_experiment()` — HPO → save; `load_and_evaluate()` — load saved model → evaluate |
| `src/evaluator.py`       | `compute_metrics()`, `build_results_table()`, `print_results_table()` |
| `src/stacking.py`        | `run_stacking_experiment()` — OOF → Ridge meta-learner → save; `stacking_exists()` (base: lgbm, xgb, extratrees, dnn) |
| `main.py`                | CLI, dataset caching, load-or-train logic, checkpointing (216 base experiments) |
| `stack.py`               | CLI, stacking ensemble over all (dataset, transform, fs) combos |
| `split_datasets.py`      | One-time CHO / CHON file generation                  |
| `outlier_detection.py`   | Williams-plot AD analysis: leverage + standardised residuals, saves plots + SMILES CSVs |
| `retrain_cleaned.py`     | Retrains three best stacking models after removing response outliers; saves to `saved_models_cleaned/` |
| `shap_analysis.py`       | SHAP feature importance plots for the best models    |
| `williams_plot.py`       | Standalone Williams plot visualisation               |
| `uncertainty_analysis.py`| Prediction uncertainty estimation                    |
| `similarity_score.py`    | Tanimoto-based training set similarity scoring       |
| `pred_plot.py`           | Predicted vs. actual scatter plots                   |
| `dataset_similarity.py`  | Inter-dataset Tanimoto similarity analysis           |
| `methodology_diagram.py` | Pipeline methodology figure for the manuscript       |
| `schematic_diagram.py`   | Schematic figure for the manuscript                  |
