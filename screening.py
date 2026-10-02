import numpy as np
import pandas as pd
from os.path import splitext
import logging
logger = logging.getLogger(__name__)
from models import random_forest, support_vector


def write_importance(clf, interaction_name, args):
    try:
        importances = clf.feature_importances_
        if not args.use_docking_score:
            interaction_name = interaction_name[:-1]
        importance = pd.DataFrame(np.array([interaction_name, importances]).T)
        outfile = splitext(args.output)[0] + ".importance"
        importance.to_csv(outfile, header=None, index=None)
        return importances
    except AttributeError:
        import traceback
        traceback.print_exc()
        logger.debug("feature importance is not available. Not Forests?")
        return None


def report(y_test, y_score, docking_score, model_name, args):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from sklearn.metrics import roc_curve, auc
    # Compute ROC curve and area the curve
    fpr, tpr, thresholds = roc_curve(y_test, y_score)
    roc_auc = auc(fpr, tpr)
    plt.plot(fpr, tpr, lw=1, label='%s (AUC = %0.3f)'
             % (model_name, roc_auc))

    # glide score, reverse order
    fpr, tpr, thresholds = roc_curve(y_test, docking_score * (-1))
    docking_auc = auc(fpr, tpr)
    plt.plot(fpr, tpr, 'r--', lw=1, label='Glide SP (AUC = %0.3f)' % docking_auc)

    plt.plot([0, 1], [0, 1], '--', lw=1, color=(0.6, 0.6, 0.6), label='Random')

    # formatting
    plt.xlim([-0.05, 1.05])
    plt.ylim([-0.05, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC: ' + args.title)
    plt.legend(loc="lower right", fontsize=9)
    outfile = splitext(args.output)[0] + "_auc.png"
    plt.savefig(outfile)
    plt.clf()

    # EF
    score_order = np.argsort(-y_score, kind="stable")
    sorted_y = y_test[score_order]
    docking_sorted_y = y_test[np.argsort(docking_score, kind="stable")]
    from calc_ef import calc_ef
    ef10 = calc_ef(sorted_y, args.active, args.decoy, threshold=0.1)
    ef1 = calc_ef(sorted_y, args.active, args.decoy, threshold=0.01)
    docking_ef10 = calc_ef(docking_sorted_y, args.active, args.decoy, threshold=0.1)
    docking_ef1 = calc_ef(docking_sorted_y, args.active, args.decoy, threshold=0.01)

    # save
    outputfile = open(splitext(args.output)[0]+"_scores.csv", "w")
    outputfile.write("target,method,auc,ef10,ef1\n")
    outputfile.write(",%s,%.3f,%.3f,%.3f\n" % (model_name, roc_auc, ef10, ef1))
    outputfile.write(",%s,%.3f,%.3f,%.3f\n" % ("Glide",docking_auc, docking_ef10, docking_ef1))
    outputfile.close()
    return {"ef10":ef10, "ef1":ef1, "roc_auc":roc_auc}


def do_screen(clf, X_train, X_test, y_train, y_test, cpd_names_test,
              model_name, interaction_name, docking_score_test, args):
    clf = clf.fit(X_train, y_train)
    probas_ = clf.predict_proba(X_test)

    # Record the requested class probability.
    if args.reverse == True:
        probas = probas_[:, 0]
    else:
        probas = probas_[:, 1]

    if args.model == "RF":
        feature_importance = write_importance(clf, interaction_name, args)

    score = probas.ravel()
    if y_test is None:
        result = pd.DataFrame({"name": cpd_names_test, "score": score})
    else:
        result = pd.DataFrame({"name": cpd_names_test, "score": score, "ishit": y_test})
    result = result.sort_values("score", ascending=False, kind="stable")
    result.to_csv(args.output, sep=",", index=None)

    if y_test is None or np.unique(y_test).size < 2:
        logger.info("No two-class test labels: scores saved without AUC/EF evaluation.")
        return None
    else:
        return report(y_test, score, docking_score_test, model_name, args)


def screening(cpd_names, label_data, interaction_name, features,
              cpd_names_test, label_data_test, interaction_name_test, features_test, args):

    if list(interaction_name) != list(interaction_name_test):
        raise ValueError("Training and test interaction feature names/order must match.")
    if args.model == "RF":
        clf = random_forest(args.random_state, args.nprocs)
        model_name = 'SIEVE-Score_RF'

    elif args.model == "SVM":
        clf = support_vector(args.random_state)
        model_name = 'SIEVE-SVM'

    docking_score_test = features_test[:, -1]
    X_train = features if args.use_docking_score else features[:, :-1]
    X_test = features_test if args.use_docking_score else features_test[:, :-1]

    y_train = np.array([1 if x > 0 else 0 for x in label_data])
    if np.unique(y_train).size != 2:
        raise ValueError("Training labels must include both active and inactive compounds.")
    if label_data_test is None:
        y_test = None
    else:
        y_test = np.array([1 if x > 0 else 0 for x in label_data_test])

    return do_screen(clf, X_train, X_test, y_train, y_test, cpd_names_test,
              model_name, interaction_name, docking_score_test, args)
