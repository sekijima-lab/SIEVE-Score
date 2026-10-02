"""Cross-release comparison: run with 0.18 first, then the modern environment.

The legacy interpreter is validation-only. Inputs contain only numeric/string
NumPy arrays, not pickled models. This script is compatible with Python 3.5.
"""
import argparse
import json
from pathlib import Path
import sys
import time

import numpy as np
import scipy
from scipy.stats import spearmanr
import sklearn
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

# Fixed before the modern comparison, not relaxed to fit observed results.
LIMITS = {'probability_mae': 0.02, 'probability_max_error': 0.10,
          'rank_spearman_min': 0.98, 'auc_abs_difference': 0.02,
          'top_fraction_active_count_difference': 2}


def metrics(y, scores):
    order = np.argsort(-scores, kind='mergesort')
    result = {'auc': float(roc_auc_score(y, scores))}
    for fraction in (0.01, 0.1):
        n = int(np.ceil(len(y) * fraction))
        hits = int(y[order[:n]].sum())
        key = 'ef1' if fraction == 0.01 else 'ef10'
        result[key] = float((hits / float(n)) / y.mean())
        result[key + '_hits'] = hits
    return result


def run(args):
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    results = []
    failed = []
    for target in args.targets:
        data = np.load(str(Path(args.inputs) / (target + '.npz')), allow_pickle=False)
        X, y = data['X'], data['y']
        if args.baseline:
            folds = np.zeros(len(y), dtype=np.int64)
            splitter = StratifiedKFold(n_splits=5, shuffle=True, random_state=1729)
            for fold, (_, test) in enumerate(splitter.split(X, y)):
                folds[test] = fold
        else:
            reference = np.load(str(Path(args.reference) / (target + '.npz')), allow_pickle=False)
            folds = reference['folds']
            from models import random_forest
        saved = {'folds': folds, 'y': y}
        cases = [(seed, False) for seed in (0, 1, 2)] + [(0, True)]
        for seed, with_docking in cases:
            started = time.time()
            features = X if with_docking else X[:, :-1]
            predictions = np.zeros(len(y), dtype=np.float64)
            importance = np.zeros(features.shape[1], dtype=np.float64)
            for fold in range(5):
                train, test = np.where(folds != fold)[0], np.where(folds == fold)[0]
                if args.baseline:
                    # These are the original scoring.py RF settings, with a fixed seed.
                    clf = RandomForestClassifier(n_estimators=1000, criterion='gini',
                        max_features=6, random_state=seed, n_jobs=args.n_jobs)
                else:
                    clf = random_forest(random_state=seed, n_jobs=args.n_jobs)
                clf.fit(features[train], y[train])
                predictions[test] = clf.predict_proba(features[test])[:, 1]
                importance += clf.feature_importances_ / 5.0
            key = 'seed{0}_dock{1}'.format(seed, int(with_docking))
            saved[key] = predictions
            saved[key + '_importance'] = importance
            record = dict(target=target, seed=seed, with_docking=with_docking,
                n_rows=len(y), n_features=features.shape[1], n_trees=1000, n_jobs=args.n_jobs,
                metrics=metrics(y, predictions), elapsed_seconds=time.time()-started)
            if not args.baseline:
                old = reference[key]
                old_metrics = metrics(y, old)
                checks = {'probability_mae':float(np.mean(np.abs(predictions-old))),
                    'probability_max_error':float(np.max(np.abs(predictions-old))),
                    'rank_spearman':float(spearmanr(predictions, old)[0]),
                    'auc_abs_difference':abs(record['metrics']['auc']-old_metrics['auc']),
                    'ef1_active_count_difference':abs(record['metrics']['ef1_hits']-old_metrics['ef1_hits']),
                    'ef10_active_count_difference':abs(record['metrics']['ef10_hits']-old_metrics['ef10_hits']),
                    'feature_importance_l1_difference':float(np.sum(np.abs(importance-reference[key+'_importance'])))}
                passed = (checks['probability_mae'] <= LIMITS['probability_mae']
                    and checks['probability_max_error'] <= LIMITS['probability_max_error']
                    and checks['rank_spearman'] >= LIMITS['rank_spearman_min']
                    and checks['auc_abs_difference'] <= LIMITS['auc_abs_difference']
                    and checks['ef1_active_count_difference'] <= LIMITS['top_fraction_active_count_difference']
                    and checks['ef10_active_count_difference'] <= LIMITS['top_fraction_active_count_difference'])
                record.update(reference_metrics=old_metrics, comparison=checks, passed=passed)
                if not passed:
                    failed.append(target + '/' + key)
            results.append(record)
            print(json.dumps(record), flush=True)
        np.savez_compressed(str(output / (target + '.npz')), **saved)
    summary = {'python':sys.version, 'numpy':np.__version__, 'scipy':scipy.__version__,
        'sklearn':sklearn.__version__, 'baseline':args.baseline, 'limits':LIMITS,
        'cases':results, 'failed_cases':failed}
    (output / 'results.json').write_text(json.dumps(summary, indent=2))
    return bool(failed)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--reference')
    parser.add_argument('--baseline', action='store_true')
    parser.add_argument('--n-jobs', type=int, default=1)
    parser.add_argument('--targets', nargs='+', choices=['CAG', 'FXI', 'HIV'],
                        default=['CAG', 'FXI', 'HIV'])
    args = parser.parse_args()
    if not args.baseline and not args.reference:
        parser.error('--reference is required for the modern comparison')
    sys.exit(run(args))
