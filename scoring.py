import numpy as np
import logging
logger = logging.getLogger(__name__)


def training(data, args):

    title = data[:, 0]
    labels = data[:, 1].astype('int')
    features = data[:, 3:].astype('float')

    #parameter search (SVM)
    from sklearn import grid_search   
    clf.fit(features, labels)
    
    return clf


def cv_accuracy(clf, features, labels):
    from sklearn import cross_validation
    scores = cross_validation.cross_val_score(clf, features, labels, cv=5)
    print("Accuracy: %0.2f (+/- %0.2f)" %(scores.mean(), scores.std()*2))

    
def scoring(test_data, args, clf):
    features = test_data[:, 3:].astype('float')
    return clf.predict(features)


def scoring_eval(data, args):

    outputfile = args.output
    score_type = args.score_type

    title = data[:, 0]
    label_data = data[:, 1].astype('int')
    labels = np.array([1 if x>0 else 0 for x in label_data])    
    features = data[:, 3:].astype('float')

    if labels.size < 5:
        logger.warning('Need more data to train.')
    elif labels[labels>0].size == 0:
        logger.warning('Training data have only non-hits.')
    elif (args.zeroneg and labels[labels<=0].size == 0 or
          not args.zeroneg and labels[labels<0].size == 0):
        logger.warning('Training data have only hits.')


    #SVM, Grid search by CV
    from sklearn import svm, grid_search
    from sklearn.cross_validation import StratifiedKFold

    svc = svm.SVC(kernel="rbf", degree=3, probability=True, cache_size=200)
    print(labels)
    skf = StratifiedKFold(labels, 5)
    param_grid = [{'kernel': ['rbf'],
                    'gamma': [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 10, 100, 1000],
                    'C': [1e-2, 1e-1, 1, 10, 100]
                   }]
    clf = grid_search.GridSearchCV(svc, param_grid, cv=skf,
                                   scoring=None, n_jobs=-1)
    clf.fit(features, labels)

    grid_scores = clf.grid_scores_
    
    np.savetxt(outputfile, grid_scores, fmt="%s", delimiter=",")
    
    print(clf.best_params_, clf.best_score_)
    score = clf.predict_proba(features)[:,1]

    rank = np.argsort(score)[::-1][:args.propose]
    cpdname = title[rank]
    score = score[rank]
    label = labels[rank]

    result = np.array(zip(cpdname, score, label))
    #print(result)
    result = np.dstack((cpdname, score, label))

    #test
    #np.savetxt(outputfile, result, fmt="%s", delimiter=",")
    logger.info('Saved SIEVE-Score.')

    return cpdname, score, label


def sort_delete_duplicate(score):
    # sort by title, score
    score = sorted(score, key=lambda x: float(x[2]), reverse=True)
    score = sorted(score, key=lambda x: x[1], reverse=True)

    ret = []
    for d in score:
        cpdname = d[1]
        if len(ret)>0:
            lastcpdname = ret[-1][1]
        else:
            lastcpdname = None

        if cpdname != lastcpdname:
            ret.append(d)

    # re-sort by score
    ret = sorted(ret, key=lambda x: float(x[2]), reverse=True)
    return ret
