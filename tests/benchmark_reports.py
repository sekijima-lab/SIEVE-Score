"""Compare actual CV report calculations using fixed legacy probabilities/folds.

Run once with the original scoring.py and 0.18, then with modern scoring.py.
This isolates reporting/API changes from retraining and CSV loader corrections.
"""
import argparse
from contextlib import redirect_stdout
import importlib.util
import io
import json
import os
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np
import sklearn
import sklearn.model_selection

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


class FrozenClassifier:
    def __init__(self, predictions, folds, importance):
        self.predictions, self.folds = predictions, folds
        self.feature_importances_ = importance
        self.fold = -1

    def fit(self, X, y):
        self.fold += 1
        return self

    def predict_proba(self, X):
        p = self.predictions[self.folds == self.fold]
        return np.column_stack((1-p, p))


def run(args):
    output = Path(args.output).absolute()
    output.mkdir(parents=True, exist_ok=True)
    os.chdir(str(output))
    sys.path.insert(0, str(Path(args.scoring_path).resolve().parent))
    spec = importlib.util.spec_from_file_location('report_scoring', args.scoring_path)
    scoring = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(scoring)
    results = []
    for target in ('CAG', 'FXI', 'HIV'):
        data = np.load(str(Path(args.inputs)/(target+'.npz')), allow_pickle=False)
        reference = np.load(str(Path(args.reference)/(target+'.npz')), allow_pickle=False)
        folds = reference['folds']

        class FrozenSplitter:
            def __init__(self, *a, **kw):
                pass

            def split(self, X, y):
                for fold in range(5):
                    yield np.where(folds != fold)[0], np.where(folds == fold)[0]

        original = sklearn.model_selection.StratifiedKFold
        sklearn.model_selection.StratifiedKFold = FrozenSplitter
        try:
            for docking in (False, True):
                key = 'seed0_dock{0}'.format(int(docking))
                clf = FrozenClassifier(reference[key], folds, reference[key+'_importance'])
                options = SimpleNamespace(use_docking_score=docking, random_state=0,
                    reverse=False, active=None, decoy=None, title=target,
                    output=str(output/(target+'_'+key+'.csv')))
                with redirect_stdout(io.StringIO()):
                    aucs, efs = scoring.cv_accuracy_plot(clf, data['X'], data['y'],
                        data['names'], 'SIEVE-Score_RF', options)
                results.append(dict(target=target, with_docking=docking,
                    aucs=np.asarray(aucs).tolist(), efs=np.asarray(efs).tolist()))
        finally:
            sklearn.model_selection.StratifiedKFold = original
    (output/'report-metrics.json').write_text(json.dumps(dict(
        python=sys.version, sklearn=sklearn.__version__, cases=results), indent=2))
    print('Saved six CV report comparisons.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs', required=True)
    parser.add_argument('--reference', required=True)
    parser.add_argument('--scoring-path', required=True)
    parser.add_argument('--output', required=True)
    run(parser.parse_args())
