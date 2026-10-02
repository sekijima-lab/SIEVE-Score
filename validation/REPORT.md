# Python / scikit-learn migration validation

Validated 2026-10-03. Original source commit: `ada055bc41758dd05e98951f9d84f221741e1f17`.

## Result and scope

The Random Forest migration meets all predeclared engineering comparison criteria
in all **12 conditions**, using **44,191 supplied compounds** across CAG, FXI and
HIV. This supports retraining the interaction-based workflow on Python 3.12.15
and scikit-learn 1.9.1. It does not establish exact reproduction of the paper's
reported metrics or cross-version compatibility of persisted models.

Across the 12 conditions:

- Largest mean absolute probability difference: 0.0000090779.
- Largest single-compound probability difference: 0.002000.
- Lowest Spearman rank correlation: 0.9991769196.
- Largest absolute pooled out-of-fold ROC-AUC difference: 0.0009655062.
- Active counts in the top 1% and 10%: identical in every condition; hence pooled
  EF1 and EF10 are identical under the specified stable tie policy.
- Largest L1 difference in fold-mean feature importance: 0.0013682830
  (diagnostic, not an acceptance criterion).

Full metrics and thresholds are in `comparison.json`; hashes, row counts and
feature counts are in `input-manifest.json`. Every condition fits five forests
of 1,000 trees in each environment. Seeds 0/1/2 exclude docking score; seed 0
also includes it. Folds are generated once by sklearn 0.18 (split seed 1729)
and reused unchanged. Both environments retain every CSV compound row.
The legacy run uses one fitting worker; the final run uses four. Parallel fitting
is separately tested against one worker on a labeled fixture.

## Actual report calculations and ties

Six additional comparisons run the original and migrated `cv_accuracy_plot`
functions with identical archived probabilities, frozen folds and compound order,
so model fitting and loader corrections cannot confound the reporting comparison.
All five-fold AUCs, interpolated mean-ROC AUCs and Glide AUCs agree to floating-point
precision (maximum difference 2.22e-16). EF1 is identical in all six comparisons.

EF10 can change when a cutoff crosses a group of compounds with identical scores.
The old code used an unstable quicksort and reversed its order. The migration uses
stable descending sorting, retaining input order for ties. With fixed probabilities,
three of the six conditions change some fold EF10 values by one selected active:

- CAG including docking: mean fold EF10 changes from 1.0 to 2.0.
- FXI excluding docking: mean fold EF10 changes from 3.328249 to 2.995932.
- HIV excluding docking: mean fold EF10 changes from 2.065006 to 2.398114.

The maximum single-fold EF10 difference is 1.666667. The affected cutoffs are
inside tied probability groups (61–72 compounds in the affected folds). This is
a tie-order effect, **not** a Random Forest or AUC discrepancy. The unchanged
pooled EF results above must not be interpreted as exact equality of every
historical fold EF. See `report-comparison.json` for all numbers. Users requiring
historical tie ordering should preserve their original ranked outputs.

## Environments and security rationale

The legacy reconstruction uses Python 3.5.5, numpy 1.11.3, scipy 0.18.1,
scikit-learn **0.18 exactly**, pandas 0.19.2 and matplotlib 2.0.2 under macOS x86_64
/Rosetta. Only sklearn 0.18 is specified by the user as the paper's version;
the other versions reconstruct a compatible period environment, not an exact
assertion about the authors' original installation. It is validation-only.

The production candidate uses Python **3.12.15**, numpy **2.5.3**, scipy **1.18.1**,
scikit-learn **1.9.1**, pandas **3.0.6** and matplotlib **3.11.2** on macOS arm64.
All 18 installed dependencies are pinned in `requirements.txt`, pass `uv pip check`,
and have no matching known OSV advisories in the recorded version query snapshot.
This is a point-in-time check, not a guarantee against undisclosed vulnerabilities.

Python 3.12.15 is the current 3.12 security release and incorporates interpreter
security fixes ([official release](https://www.python.org/downloads/release/python-31215/)).
The legacy sklearn release is outside the project's supported security releases
([security policy](https://github.com/scikit-learn/scikit-learn/security/policy)).
OSV matches sklearn 0.18 to CVE-2024-5206 and disputed pickle-deserialization
advisories. CVE-2024-5206 concerns TfidfVectorizer, which this RF workflow does not
use; it is not evidence of an exploitable RF-specific issue
([advisory](https://github.com/advisories/GHSA-jw8x-6495-233v)).
Upgrading does not make untrusted pickle/joblib loading safe. The supported
pipeline retrains from CSV and uses no serialized models. The benchmark archives
numeric-only NPZ files loaded with `allow_pickle=False`.

Because the installed managers did not yet offer Python 3.12.15 binaries, the
final local interpreter was built from the official source archive in a new
workspace prefix, with OpenSSL 3.6.4. Source SHA-256:
`c2c4321961fab0fb999d66e0cecf521c2ab3994c7992873ea99e306c1094fd5a`.
Existing Python, micromamba, uv and Schrödinger environments were not upgraded.

## Supported workflow verification

Nine automated tests cover first-row and feature-name preservation, exact-title
filtering, sidecar lookup without importing Schrödinger, invalid CSV rejection,
RF worker reproducibility, RF/SVM CV score/importance/AUC outputs, screening with
and without docking score and reversed probabilities, feature-order validation,
unlabeled screening, RF and SVM parameter search, datasize analysis, comparesvm,
repeated importance, docking-score extraction, histogram tools and numeric EF CSV
handling. All pass locally on Python 3.12.15. GitHub rejected adding a CI workflow because
the authenticated OAuth token lacks the `workflow` scope. The workflow is omitted
from this PR; no remote CI or Ubuntu validation is claimed.

Deliberate fixes can change results independently of library versions:

- Include the first real compound, previously incorrectly discarded.
- Associate every importance value with the correct feature name.
- Sort scores numerically, with an explicit stable tie policy.
- Honor `--use_docking_score` during independent screening.
- Preserve compound titles containing spaces when applying `--ignore`.
- Replace removed sklearn/scipy/pandas/matplotlib APIs and repair previously
  failing comparison, parameter-search and default datasize paths.

Schrödinger conversion is unchanged and unvalidated. Precomputed CSVs form the
supported boundary. Research scripts with hard-coded personal paths and other
hardware/platform combinations are not certified by this validation.

Reproduction instructions are in `REPRODUCE.md`; numeric reference predictions,
folds and fold-mean importances are committed in `baseline/`.
