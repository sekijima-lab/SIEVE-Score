from sklearn.metrics import roc_curve, roc_auc_score
import csv
import matplotlib.pyplot as plt
from os.path import splitext

def plot_ROC(n_actives,n_decoys, title, outputname,
            results, show):
    colors = ['b','g','r','c','m','y','k']

    results_file=results[0::3]
    results_type=results[1::3]
    results_name=results[2::3]
    
    ys=[]
    scores=[]

    for i in range(len(results_file)):
        y,score = GetYScoreFromResult(results_file[i], results_type[i])
        ys.append(y)
        scores.append(score)

    plt.figure(figsize=(4,4), dpi=150)
    SetupROCCurvePlot(plt,title)

    aucs = []
    for i in range(len(results_file)):
        aucs.append(AddROCCurve(plt, ys[i], scores[i], colors[i], results_name[i], n_actives, n_decoys))
    
    SaveROCCurvePlot(plt, outputname, show, True)
    with open(splitext(outputfile)[0]+'result.txt',"w")\
             as f_result:
        
        auc_sp = 0
        auc_proposed = 0
        for i in range(len(results_file)):
            f_result.write(results_name[i]+
                            ", "+str(aucs[i])+"\n")
            if results_name[i]=='Glide_SP':
                auc_sp = aucs[i]
            elif 'proposed' in results_name[i]:
                auc_proposed = aucs[i]

        if auc_sp == 0 or auc_proposed == 0:
            f_result.write("Cannot compared.")
        elif auc_proposed - auc_sp > 0.001:
            f_result.write("proposed win!")
        elif auc_proposed - auc_sp < -0.001:
            f_result.write("proposed lose...")
        else:
            f_result.write("draw")

def GetYScoreFromResult(filename,datatype):
    y=[]
    score=[]
    data = csv.reader(open(filename, 'rb'), delimiter=',', quotechar='#')

    if datatype == 'dock':
        #print(filename)
        for line in data:
            y.append(1 if "ZINC" not in line[0] else 0)
            score.append(float(line[1]))
    
    else:
    #elif datatype=='stat':
        for line in data:
            y.append(int(line[3]))
            score.append(float(line[2]))

    return y, score

def GetRates(y, scores, n_actives, n_decoys):

    tpr = [0.0]  # true positive rate
    fpr = [0.0]  # false positive rate

    foundactives = 0.0
    founddecoys = 0.0
    for idx, score in enumerate(scores):
        if y[idx]==1:
            foundactives += 1.0
        else:
            founddecoys += 1.0

        tpr.append(foundactives / float(n_actives))
        fpr.append(founddecoys / float(n_decoys))

    return tpr, fpr
    
def SetupROCCurvePlot(plt,title):

    plt.xlabel("False Positive rate", fontsize=14)
    plt.ylabel("True Positive rate", fontsize=14)
    #plt.title("ROC Curve ("+title+")", fontsize=14)
    plt.title(title[2:], fontsize=14)
def AddROCCurve(plt, actives, scores, color, label, n_actives, n_decoys):
    tpr, fpr = GetRates(actives, scores, n_actives, n_decoys)
    #print(actives,label)
    if "Glide" in label:
        scores = [-x for x in scores] #reverse order
    print(tpr)

    auc_tmp = roc_auc_score(actives, scores)
    auc_tmp = (auc_tmp*tpr[-1])#+((1-fpr[-1])*tpr[-1]) #adjust final position
    auc = round(auc_tmp,3)

    label = label+", auc="+str(auc)

    if "Glide" in label:
        plt.plot(fpr, tpr, linestyle='dashed', color=color, linewidth=2, label=label)
    else:
        plt.plot(fpr, tpr, color=color, linewidth=2, label=label)
    return auc

def SaveROCCurvePlot(plt, outputname, show, randomline=True):

    if randomline:
        x = [0.0, 1.0]
        plt.plot(x, x, linestyle='dashed', color='red', linewidth=2, label='random')

    plt.xlim(0.0, 1.0)
    plt.ylim(0.0, 1.0)
    #plt.legend(fontsize=10, loc='best')
    plt.tight_layout()
    plt.savefig(outputname)
    if show:
        plt.show()
    

if __name__ == '__main__':
    import sys
    x = sys.argv
    if len(x)<=5 and len(x)%3!=2:
        print("usage: n_actives, n_decoys, graph_title, graph_filename,"
        " results(file, type, name)...")
    
    n_actives = int(x[1])
    n_decoys = int(x[2])
    title = x[3]
    outputname = x[4]

    show = x[5]
    if show in ["False", "false", "0", "No", "no"]:
        show = False
    else:
        show = True

    results = x[6:]

    plot_ROC(n_actives,n_decoys,title,outputname,
            results, show)
    

