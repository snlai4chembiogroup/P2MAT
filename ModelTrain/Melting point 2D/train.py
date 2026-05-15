######################################################
# Developer : Methun Kamruzzaman
# Date      : April 9, 2026
# Purpose   : Orchestrates the complete melting point (2D) QSAR machine-learning pipeline.
######################################################
"""
train.py
-------
Orchestrates the complete QSAR machine-learning pipeline.

Experiment grid
---------------
  Datasets            : full | cho | chon
  Transformations     : none | standard | quantile
  Feature selections  : whole | boruta | pca | correlation
  Models              : linear | elasticnet | extratrees | lgbm | xgb | dnn

  Total               : 3 × 3 × 4 × 6 = 216 experiments

Each experiment uses 10-fold cross-validation with Bayesian optimisation
(neg_mean_absolute_error) to find the best hyperparameters, then evaluates
on the held-out test set.

Usage
-----
    conda run -n qsar python train.py [--datasets full cho chon]
                                      [--transforms none standard quantile]
                                      [--fs whole boruta pca correlation]
                                      [--models linear elasticnet extratrees lgbm xgb dnn]
                                      [--retrain]

    Default : if a saved model exists → load and evaluate it (no HPO).
              if no saved model exists → run full HPO and save.
    --retrain : force retraining even when a saved model already exists.

Outputs
-------
    results/all_results.csv        – full metric table
    results/summary_test_mae.csv   – pivot table (best test MAE per model)
    results/summary_train_mae.csv  – same for training phase
"""

import argparse
import os
import sys
import json
import warnings
import pandas as pd

warnings.filterwarnings("ignore")

# ---------------------------------------------------------------------------
# Path setup so src/ is importable from the project root
# ---------------------------------------------------------------------------
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from src.config import (
    DATASETS,
    TRANSFORMATIONS,
    FEATURE_SELECTIONS,
    MODEL_NAMES,
    RESULTS_DIR,
    MODELS_DIR,
)
from src.data_loader import load_dataset
from src.trainer import run_experiment, load_and_evaluate, _pipeline_path
from src.evaluator import (
    build_results_table,
    build_full_report,
    print_results_table,
)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(
        description="QSAR Melting Point Prediction Pipeline"
    )
    p.add_argument("--datasets",   nargs="+", default=list(DATASETS.keys()),
                   choices=list(DATASETS.keys()),
                   help="Datasets to run (default: all)")
    p.add_argument("--transforms", nargs="+", default=TRANSFORMATIONS,
                   choices=TRANSFORMATIONS,
                   help="Transformations to run (default: all)")
    p.add_argument("--fs",         nargs="+", default=FEATURE_SELECTIONS,
                   choices=FEATURE_SELECTIONS,
                   help="Feature selection methods to run (default: all)")
    p.add_argument("--models",     nargs="+", default=MODEL_NAMES,
                   choices=MODEL_NAMES,
                   help="Models to run (default: all)")
    p.add_argument("--retrain",    action="store_true",
                   help="Force retraining even if a saved model already exists")
    return p.parse_args()


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _model_exists(dataset, fs, transform, model):
    """Return True if a saved pipeline.joblib exists for this combo."""
    return os.path.exists(_pipeline_path(dataset, fs, transform, model))


def _total_experiments(datasets, transforms, fs_list, models):
    return len(datasets) * len(transforms) * len(fs_list) * len(models)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    args = parse_args()

    total = _total_experiments(args.datasets, args.transforms, args.fs, args.models)
    print(f"\n{'='*70}")
    print(f"  QSAR Melting Point Prediction  –  {total} experiments")
    print(f"  Datasets   : {args.datasets}")
    print(f"  Transforms : {args.transforms}")
    print(f"  FS methods : {args.fs}")
    print(f"  Models     : {args.models}")
    print(f"{'='*70}\n")

    # Always load existing results so previous runs are never lost.
    results_path = os.path.join(RESULTS_DIR, "all_results.csv")
    existing_results = []
    completed_keys   = set()          # (dataset, transform, fs, model) already in CSV
    if os.path.exists(results_path):
        old_df = pd.read_csv(results_path)
        existing_results = old_df.to_dict("records")
        completed_keys = {
            (r["dataset"], r["transform"], r["fs_method"], r["model"])
            for r in existing_results
        }
        print(f"Found {len(existing_results)} existing results — "
              f"new results will be appended.\n")

    new_results = []
    done_count  = 0
    skip_count  = 0
    exp_num     = 0

    # Cache loaded datasets to avoid repeated I/O
    _data_cache = {}

    for dataset in args.datasets:
        ds_cfg = DATASETS[dataset]

        if dataset not in _data_cache:
            print(f"\nLoading dataset: {dataset.upper()} "
                  f"({ds_cfg['train']} / {ds_cfg['test']}) ...", flush=True)
            X_train, X_test, y_train, y_test, feat_names = load_dataset(
                ds_cfg["train"], ds_cfg["test"]
            )
            print(f"  → {X_train.shape[0]} train  |  {X_test.shape[0]} test  "
                  f"|  {X_train.shape[1]} features")
            _data_cache[dataset] = (X_train, X_test, y_train, y_test)

        X_train, X_test, y_train, y_test = _data_cache[dataset]

        for transform in args.transforms:
            for fs in args.fs:
                for model in args.models:
                    exp_num += 1
                    tag = f"{dataset}/{fs}/{transform}/{model}"

                    combo_key     = (dataset, transform, fs, model)
                    has_results   = combo_key in completed_keys
                    already_built = _model_exists(dataset, fs, transform, model)

                    # Priority order (unless --retrain is set):
                    #   1. Results already in CSV        → SKIP  (nothing to do)
                    #   2. Model saved, no results entry → LOAD  (evaluate only)
                    #   3. No model, no results          → TRAIN (full HPO)
                    if has_results and not args.retrain:
                        skip_count += 1
                        print(f"[{exp_num:3d}/{total}] SKIP  {tag}  (results already recorded)")
                        continue
                    elif already_built and not args.retrain:
                        print(f"\n[{exp_num:3d}/{total}] LOAD  {tag}  (model saved, results missing)")
                        action = "load"
                    else:
                        reason = "(--retrain)" if args.retrain else ""
                        print(f"\n[{exp_num:3d}/{total}] TRAIN {tag} {reason}")
                        action = "train"

                    try:
                        if action == "load":
                            result = load_and_evaluate(
                                dataset    = dataset,
                                transform  = transform,
                                fs_method  = fs,
                                model_name = model,
                                X_train    = X_train,
                                X_test     = X_test,
                                y_train    = y_train,
                                y_test     = y_test,
                                verbose    = True,
                            )
                        else:
                            result = run_experiment(
                                dataset    = dataset,
                                transform  = transform,
                                fs_method  = fs,
                                model_name = model,
                                X_train    = X_train,
                                X_test     = X_test,
                                y_train    = y_train,
                                y_test     = y_test,
                                verbose    = True,
                            )

                        new_results.append(result)
                        done_count += 1

                        # Checkpoint: save after every experiment
                        _checkpoint(existing_results + new_results, results_path)

                    except Exception as exc:
                        print(f"  [ERROR] {tag}: {exc}")
                        import traceback
                        traceback.print_exc()

    # -----------------------------------------------------------------------
    # Final reporting
    # -----------------------------------------------------------------------
    all_results = existing_results + new_results

    if not all_results:
        print("No results to report.")
        return

    full_df = build_full_report(all_results)
    full_df.to_csv(results_path, index=False)

    test_table  = build_results_table(all_results, phase="test")
    train_table = build_results_table(all_results, phase="train")

    test_table.to_csv(os.path.join(RESULTS_DIR, "summary_test_mae.csv"))
    train_table.to_csv(os.path.join(RESULTS_DIR, "summary_train_mae.csv"))

    print_results_table(test_table,  title="Best Test  MAE (K) by Dataset × Transformation × Model")
    print_results_table(train_table, title="Best Train MAE (K) by Dataset × Transformation × Model")

    print(f"Results saved to: {RESULTS_DIR}")
    print(f"Models  saved to: {MODELS_DIR}")
    print(f"\nCompleted: {done_count} new  |  Skipped: {skip_count}  |  Total: {total}\n")


def _checkpoint(results: list, path: str):
    """Save intermediate results after each experiment."""
    if results:
        pd.DataFrame(results).to_csv(path, index=False)


if __name__ == "__main__":
    main()
