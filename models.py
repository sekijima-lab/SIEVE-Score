"""Keep the published estimator settings explicit across sklearn versions."""
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC


def random_forest(random_state=None, n_jobs=1, max_features=6):
    return RandomForestClassifier(
        n_estimators=1000,
        criterion="gini",
        max_features=max_features,
        max_depth=None,
        min_samples_split=2,
        min_samples_leaf=1,
        bootstrap=True,
        class_weight=None,
        random_state=random_state,
        n_jobs=n_jobs,
    )


def support_vector(random_state=None, probability=True):
    return SVC(C=10, kernel="rbf", degree=3, gamma=0.1,
               cache_size=1000, probability=probability, random_state=random_state)
