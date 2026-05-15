######################################################
# Developer : Methun Kamruzzaman
# Date      : April 12, 2026
# Purpose   : Entry point for stacking ensemble training over all pipeline
#             combinations (melting point 2D).
######################################################
"""
stack.py
--------
Entry point for running stacking ensembles (LightGBM + XGBoost + DNN)
over all (dataset, transform, fs_method) combinations.

Requires the base models to already be trained by main.py.

Usage
-----
    conda run -n qsar python stack.py
    conda run -n qsar python stack.py --datasets full cho --fs whole boruta
    conda run -n qsar python stack.py --retrain
"""

import argparse
import os
import sys
import warnings
import pandas as pd

warnings.filterwarnings("ignore")

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from src.config import (
    DATASETS,
    TRANSFORMATIONS,
    FEATURE_SELECTIONS,
    RESULTS_DIR,
    TRANSFORM_DISPLAY_NAMES,
)
from src.data_loader import load_dataset
from src.stacking import run_stacking_experiment, stacking_exists, BASE_MODELS


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="QSAR Stacking Ensemble: LightGBM + XGBoost + DNN"
    )
    p.add_argument(
        "--datasets",  nargs="+", default=list(DATASETS.keys()),
        choices=list(DATASETS.keys()),
        help="Datasets to run (default: all)",
    )
    p.add_argument(
        "--transforms", nargs="+", default=TRANSFORMATIONS,
        choices=TRANSFORMATIONS,
        help="Transformations to run (default: all)",
    )
    p.add_argument(
        "--fs", nargs="+", default=FEATURE_SELECTIONS,
        choices=FEATURE_SELECTIONS, dest="fs_methods",
        help="Feature selection methods to run (default: all)",
    )
    p.add_argument(
        "--retrain", action="store_true",
        help="Re-run even if a stacking result already exists",
    )
    p.add_argument(
        "--retrain-oof", action="store_true", dest="retrain_oof",
        help="Recompute OOF predictions for all base models (ignores cache)",
    )
    return p.parse_args()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args  = parse_args()
    total = len(args.datasets) * len(args.transforms) * len(args.fs_methods)

    print(f"\n{'='*70}")
    print(f"  QSAR Stacking Ensemble  –  up to {total} experiments")
    print(f"  Base models : {list(BASE_MODELS)}")
    print(f"  Datasets    : {args.datasets}")
    print(f"  Transforms  : {args.transforms}")
    print(f"  FS methods  : {args.fs_methods}")
    print(f"{'='*70}\n")

    out_path = os.path.join(RESULTS_DIR, "stacking_results.csv")

    # Load previously saved stacking results
    existing_results = []
    completed_keys   = set()
    if os.path.exists(out_path):
        old_df = pd.read_csv(out_path)
        existing_results = old_df.to_dict("records")
        completed_keys = {
            (r["dataset"], r["transform"], r["fs_method"])
            for r in existing_results
        }
        print(f"Found {len(existing_results)} existing stacking results — "
              f"new results will be appended.\n")

    new_results = []
    data_cache  = {}
    exp_num     = 0

    for dataset in args.datasets:
        ds_cfg = DATASETS[dataset]

        if dataset not in data_cache:
            print(f"Loading dataset: {dataset.upper()} ...")
            X_train, X_test, y_train, y_test, _ = load_dataset(
                ds_cfg["train"], ds_cfg["test"]
            )
            print(f"  → {X_train.shape[0]} train | {X_test.shape[0]} test "
                  f"| {X_train.shape[1]} features")
            data_cache[dataset] = (X_train, X_test, y_train, y_test)

        X_train, X_test, y_train, y_test = data_cache[dataset]

        for transform in args.transforms:
            for fs_method in args.fs_methods:
                exp_num += 1
                key = (dataset, transform, fs_method)
                tag = f"{dataset}/{fs_method}/{transform}/stack"

                if key in completed_keys and not args.retrain:
                    print(f"[{exp_num:3d}/{total}] SKIP  {tag}  (already recorded)")
                    continue

                print(f"\n[{exp_num:3d}/{total}] RUN   {tag}")

                try:
                    result = run_stacking_experiment(
                        dataset=dataset,
                        transform=transform,
                        fs_method=fs_method,
                        X_train=X_train,
                        X_test=X_test,
                        y_train=y_train,
                        y_test=y_test,
                        retrain_oof=args.retrain_oof,
                        verbose=True,
                    )

                    if result is not None:
                        new_results.append(result)
                        # Checkpoint after each experiment
                        all_so_far = existing_results + new_results
                        pd.DataFrame(all_so_far).to_csv(out_path, index=False)

                except Exception as exc:
                    print(f"  [ERROR] {tag}: {exc}")
                    import traceback
                    traceback.print_exc()

    # -----------------------------------------------------------------------
    # Final summary
    # -----------------------------------------------------------------------
    all_results = existing_results + new_results

    if not all_results:
        print("\nNo stacking results to report.")
        return

    df = pd.DataFrame(all_results)
    df.to_csv(out_path, index=False)

    print(f"\n{'='*80}")
    print("  STACKING RESULTS — Test MAE by Dataset × Feature Selection")
    print(f"{'='*80}")

    summary = df.pivot_table(
        index   = ["dataset", "fs_method"],
        columns = "transform",
        values  = "test_mae",
        aggfunc = "min",
    )
    summary.columns = [
        TRANSFORM_DISPLAY_NAMES.get(c, c) for c in summary.columns
    ]
    print(summary.to_string(float_format=lambda x: f"{x:.3f}"))
    print(f"{'='*80}")
    print(f"\nResults saved to : {out_path}")
    print(f"Models  saved to : saved_models/{{dataset}}/{{fs}}/{{transform}}/stack/\n")


if __name__ == "__main__":
    main()
