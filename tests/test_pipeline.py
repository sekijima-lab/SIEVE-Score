"""Exercise the supported precomputed-CSV boundary and actual CLI workflows."""
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from interaction_io import load_interaction
from models import random_forest
from calc_ef import calc_ef, get_y_score_from_result


class PipelineTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.work = Path(self.tmp.name)
        rng = np.random.RandomState(11)
        self.names = ['compound %d' % i for i in range(120)]
        self.y = np.repeat([0, 1], 60)
        values = rng.normal(size=(120, 6))
        values[:, 0] += self.y * 2
        self.frame = pd.DataFrame(values, index=self.names,
                                   columns=['R%d_vdw' % i for i in range(6)])
        self.frame.insert(0, 'ishit', self.y)
        self.frame['docking_score'] = -self.y + rng.normal(size=120)
        self.frame.index.name = '# title'
        self.input = self.work / 'train.interaction'
        self.frame.to_csv(self.input)

    def tearDown(self):
        self.tmp.cleanup()

    def cli(self, script='SIEVE-Score.py', *args):
        result = subprocess.run([sys.executable, str(ROOT / script), '-i', str(self.input),
            '-z', '--random_state', '7', *args], cwd=self.work,
            capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        return result

    def test_loader_keeps_first_compound_and_all_features(self):
        names, y, features, X = load_interaction(self.input)
        self.assertEqual(list(names), self.names)
        self.assertEqual(features, list(self.frame.columns[1:]))
        np.testing.assert_array_equal(y, self.y)
        np.testing.assert_allclose(X, self.frame.iloc[:, 1:].to_numpy())
        ignore = self.work / 'ignore.txt'
        ignore.write_text('compound 0\ncompound 4\n')
        names, _, _, _ = load_interaction(self.work / 'train.maegz', ignore=ignore)
        self.assertEqual(len(names), 118)
        self.assertNotIn('compound 0', names)
        self.assertNotIn('schrodinger', sys.modules)

    def test_bad_input_rejected(self):
        frame = self.frame.copy()
        frame.iloc[0, 1] = np.inf
        frame.to_csv(self.input)
        with self.assertRaisesRegex(ValueError, 'finite'):
            load_interaction(self.input)
        self.frame.drop(columns='ishit').to_csv(self.input)
        with self.assertRaisesRegex(ValueError, 'ishit'):
            load_interaction(self.input)
        self.input.write_text('# title,ishit,x,x,docking_score\na,1,2,3,4\n')
        with self.assertRaisesRegex(ValueError, 'unique'):
            load_interaction(self.input)
        self.input.unlink()
        with self.assertRaisesRegex(FileNotFoundError, 'separate'):
            load_interaction(self.input)

    def test_rf_settings_and_parallel_reproducibility(self):
        X = self.frame.iloc[:, 1:-1].to_numpy()
        a, b = random_forest(7, 1), random_forest(7, 2)
        self.assertEqual(a.n_estimators, 1000)
        self.assertEqual(a.max_features, 6)
        self.assertEqual(a.criterion, 'gini')
        np.testing.assert_allclose(a.fit(X, self.y).predict_proba(X),
                                   b.fit(X, self.y).predict_proba(X), atol=1e-14)

    def test_cv_rf_and_svm(self):
        for model in ['RF', 'SVM']:
            self.cli('SIEVE-Score.py', '--model', model)
            name = 'SIEVE-Score_RF' if model == 'RF' else 'SIEVE-SVM'
            scores = pd.read_csv(self.work / ('result_' + name + '_scores.csv'), header=None)
            self.assertEqual(len(scores), 120)
            self.assertEqual(set(scores[0]), set(self.names))
            self.assertTrue(scores[1].between(0, 1).all())
            self.assertTrue(scores[1].is_monotonic_decreasing)
            self.assertEqual(np.loadtxt(self.work / 'SIEVE-Score_auc.csv').shape, (7,))
            self.assertTrue((self.work / 'SIEVE-Score_auc.png').stat().st_size > 1000)
        importance = pd.read_csv(self.work / 'SIEVE-Score.importance', header=None)
        self.assertEqual(list(importance[0]), list(self.frame.columns[1:-1]))
        self.assertAlmostEqual(importance[1].sum(), 1)

    def test_screen_flags_and_feature_order(self):
        test = self.work / 'test.interaction'
        self.frame.to_csv(test)
        self.cli('SIEVE-Score.py', '-m', 'screen', '--testdata', str(test))
        scores = pd.read_csv(self.work / 'SIEVE-Score.csv')
        X = self.frame.iloc[:, 1:-1].to_numpy()
        expected = random_forest(7).fit(X, self.y).predict_proba(X)[:, 1]
        np.testing.assert_allclose(scores.set_index('name').loc[self.names, 'score'], expected)
        self.cli('SIEVE-Score.py', '-m', 'screen', '--testdata', str(test), '--use_docking_score', '--reverse')
        scores = pd.read_csv(self.work / 'SIEVE-Score.csv')
        X = self.frame.iloc[:, 1:].to_numpy()
        expected = random_forest(7).fit(X, self.y).predict_proba(X)[:, 0]
        np.testing.assert_allclose(scores.set_index('name').loc[self.names, 'score'], expected)
        self.assertEqual(len(pd.read_csv(self.work / 'SIEVE-Score.importance', header=None)), 7)
        from screening import screening
        from types import SimpleNamespace
        with self.assertRaisesRegex(ValueError, 'names/order'):
            screening(self.names, self.y, ['a', 'docking_score'], X,
                      self.names, self.y, ['b', 'docking_score'], X, SimpleNamespace())

    def test_unlabelled_screen_has_scores_without_report(self):
        test = self.work / 'test.interaction'
        frame = self.frame.copy()
        frame['ishit'] = 0
        frame.to_csv(test)
        self.cli('SIEVE-Score.py', '-m', 'screen', '--testdata', str(test))
        self.assertEqual(len(pd.read_csv(self.work / 'SIEVE-Score.csv')), 120)
        self.assertFalse((self.work / 'SIEVE-Score_auc.png').exists())

    def test_other_modes(self):
        self.cli('SIEVE-Score.py', '-m', 'datasize', '--n_splits', '2')
        self.assertEqual(pd.read_csv(self.work / 'SIEVE-Score.csv', index_col=0).shape, (3, 2))
        small = self.frame[['ishit', 'R0_vdw', 'R1_vdw', 'docking_score']]
        small.to_csv(self.input)
        self.cli('SIEVE-Score.py', '-m', 'paramsearch')
        grid = pd.read_csv(self.work / 'SIEVE-Score.csv')
        self.assertEqual(list(grid.max_features), ['sqrt'])
        self.cli('SIEVE-Score.py', '-m', 'paramsearch', '--model', 'SVM')
        self.assertEqual(len(pd.read_csv(self.work / 'SIEVE-Score.csv')), 9)
        self.cli('SIEVE-Score.py', '-m', 'comparesvm')
        self.assertTrue((self.work / 'result_SIEVE-SVM_scores.csv').exists())
        self.assertTrue((self.work / 'result_SIEVE-Score_RF_scores.csv').exists())
        self.cli('multi_importance.py', '--use_docking_score', '--n_iter', '2')
        self.assertEqual(pd.read_csv(self.work / 'multiple_importance.csv').shape, (2, 3))
        self.cli('extract_gscore.py')
        self.assertEqual(len(pd.read_csv(self.work / 'glide_score.csv')), 120)

    def test_histogram_csv_tools(self):
        histogram = self.work / 'distances.csv'
        histogram.write_text('1,2,3,4\n2,3,4,5\n3,4,5,6\n')
        for script in ['make_histogram.py', 'make_histogram_2.py']:
            output = self.work / script.removesuffix('.py')
            result = subprocess.run([sys.executable, str(ROOT / script),
                str(histogram), str(output)], cwd=self.work,
                capture_output=True, text=True, env={**__import__('os').environ,
                                                     'MPLBACKEND': 'Agg'})
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertTrue(Path(str(output) + '_a.png').stat().st_size > 1000)

    def test_enrichment_and_numeric_csv(self):
        self.assertEqual(calc_ef([1] * 100 + [0] * 4900, threshold=.1), 10)
        self.assertEqual(calc_ef([1] * 100 + [0] * 4900, threshold=.01), 50)
        with self.assertRaises(ValueError):
            calc_ef([1, 0], threshold=0)
        result, active = self.work / 'scores.csv', self.work / 'active.csv'
        result.write_text('compound 0,1e-8\ncompound 1,0.9\n')
        active.write_text('compound 1,1\ncompound 0,0\n')
        y, score = get_y_score_from_result(result, active)
        np.testing.assert_array_equal(y, [1, 0])
        np.testing.assert_allclose(score, [.9, 1e-8])


if __name__ == '__main__':
    unittest.main()
