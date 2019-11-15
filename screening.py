import numpy as np
import pandas as pd
from os.path import splitext
import logging
logger = logging.getLogger(__name__)


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
    tpr = [0.0] + tpr
    fpr = [0.0] + fpr
    roc_auc = auc(fpr, tpr)
    plt.plot(fpr, tpr, lw=1, label='%s (AUC = %0.3f)'
             % (model_name, roc_auc))

    # glide score, reverse order
    fpr, tpr, thresholds = roc_curve(y_test, docking_score * (-1))
    tpr = [0.0] + tpr
    fpr = [0.0] + fpr
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
    score_order = np.argsort(y_score)[::-1]
    sorted_score = y_score[score_order]
    sorted_y = y_test[score_order]
    docking_sorted_y = y_test[docking_score.argsort()]
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
    from sklearn.metrics import roc_curve, auc
    clf = clf.fit(X_train, y_train)
    probas_ = clf.predict_proba(X_test)
    
    # record scores, having bug for some target?
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
    result = result.sort_values("score", ascending=False)
    result.to_csv(args.output, sep=",")

    if y_test is None:
        return None
    else:
        return report(y_test, score, docking_score_test, model_name, args)
        

def screening(cpd_names, label_data, interaction_name, features, 
              cpd_names_test, label_data_test, interaction_name_test, features_test, args):

    random_state = 0
    if args.model == "RF":
        from sklearn.ensemble import RandomForestClassifier as RFC
        clf = RFC(n_estimators=1000, criterion='gini', max_features=6, 
                         n_jobs=args.nprocs)
        model_name = 'SIEVE-Score_RF'

    elif args.model == "SVM":
        from sklearn.svm import SVC
        clf = SVC(C=10, kernel="rbf", degree=3, gamma=0.1,
                         cache_size=1000, probability=True)
        model_name = 'SIEVE-SVM'

    docking_score = features[:, -1]
    docking_score_test = features_test[:, -1]
    X_train = features[:, :-1]
    X_test = features_test[:, :-1]

    y_train = np.array([1 if x > 0 else 0 for x in label_data])
    if label_data_test is None:
        y_test = None
    else:
        y_test = np.array([1 if x > 0 else 0 for x in label_data_test])

    return do_screen(clf, X_train, X_test, y_train, y_test, cpd_names_test, 
              model_name, interaction_name, docking_score_test, args)
