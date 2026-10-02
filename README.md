# SIEVE-Score: precomputed interaction workflow

The machine-learning pipeline now runs independently of Schrödinger on Python
3.12. Generate the `.interaction` CSV from Glide output in your separate
Schrödinger environment, then copy that CSV into this environment. The converter
in `read_interaction.py` is unchanged and is not supported by the new environment.

## Isolated installation

Use a new directory/environment; do not upgrade a shared or Schrödinger environment.
The tested interpreter is Python 3.12.15. `.python-version` selects that security release.

```sh
uv venv --python 3.12.15 .venv
uv pip install --python .venv/bin/python -r requirements.txt
uv pip check --python .venv/bin/python
```

If a manager does not yet distribute 3.12.15, install that interpreter separately
and pass its full executable path to `uv venv --python`. Do not silently fall
back to an older patch release.

Alternatively, use a separate micromamba environment containing Python 3.12.15 and
pip, then run `python -m pip install -r requirements.txt` inside it. The pinned
requirements include the transitive dependencies tested on macOS arm64; other
platforms must run the validation suite before production use.

## Input contract

The CSV has one header row and one row per compound, as emitted by the original
converter:

```csv
# title,ishit,RES1_vdw,RES1_coul,RES1_hbond,docking_score
compound A,1,-2.5,-0.1,0,-7.2
compound B,0,-1.0,0,0,-5.1
```

- Compound titles and feature names must be unique. Titles may contain spaces.
- `ishit` and every feature must be finite numeric values. With `-z`, labels
  greater than zero are active; all other labels are inactive.
- Feature columns preserve CSV order; `docking_score` is placed last. Training
  and screening CSVs must contain identical feature names in identical order.
- `--ignore` accepts one exact compound title per line.
- `--hits` does not override labels in a precomputed CSV. Label that CSV correctly
  in the conversion environment.
- `-i docking.maegz` still finds the sibling `docking.interaction`; direct CSV
  input is preferred. A missing sidecar requires the separate conversion environment.

## Evaluate and screen

```sh
.venv/bin/python SIEVE-Score.py -i training.interaction -z --random_state 0
.venv/bin/python SIEVE-Score.py -i training.interaction -z --random_state 0 \
  -m screen --testdata screening.interaction -o screening.csv
```

The Random Forest keeps 1,000 trees, Gini criterion, `max_features=6`, bootstrap
sampling and the original remaining tree settings. `--use_docking_score` includes
the final docking feature; otherwise it is excluded. `--random_state` now controls
both model fitting and shuffled splits; omission preserves nondeterministic runs.
`--nprocs` controls parallel fitting. `--model SVM`, `paramsearch`, `datasize`, and
`comparesvm` are also supported. RF parameter search maps the removed `auto`
setting to its classification equivalent, `sqrt`. SVM parameter search now
searches C={1,10,100}, gamma={0.01,0.1,1}; its previous grid was invalid.

CV writes `result_SIEVE-Score_RF_scores.csv` (title, probability, no header),
`SIEVE-Score.importance` (feature, mean importance), and ROC/AUC/EF artifacts
beside the working directory. Screening writes a headed score CSV to `-o` and
adds evaluation artifacts when both active and inactive test labels exist.
A screening CSV with all labels zero produces scores without an AUC/EF report.
`--active`/`--decoy` are optional numeric counts for enrichment calculations;
when omitted, counts come from the evaluated labels. Passing molecule filenames
for these counts still requires Schrödinger.

`extract_gscore.py` and `multi_importance.py` accept the same precomputed input.
The latter uses seed+i for repetition i when a seed is supplied.

## Validation and behavior changes

```sh
.venv/bin/python -m unittest discover -s tests -v
```

See `validation/REPORT.md` for comparisons against the paper's scikit-learn 0.18
release. These compare newly fitted models on identical rows, features, seeds,
and frozen folds, rather than assuming cross-version pickle compatibility.

The old loader accidentally discarded the first real compound and omitted the
first feature name from importance output. Both are corrected. Numeric score
sorting and exact-title filtering are corrected as well; these changes can alter
results independently of the library upgrade. Modern `StratifiedKFold` may produce
different folds even with the same seed. Exact reproduction requires using the
archived folds in the benchmark, not just specifying a CLI seed. Stable sorting preserves input order for ties. Some fold EF10 values change
when their cutoffs fall inside tied scores; the validation report quantifies
these differences separately from the model comparison. Historical paper
metrics are not claimed to be reproduced.

The standalone research plotting scripts with hard-coded personal paths
(`curation.py`, `plot_importance.py`) and the Schrödinger converter are outside the
supported workflow. Do not load legacy pickle/joblib models into the new runtime;
retrain from interaction data.
