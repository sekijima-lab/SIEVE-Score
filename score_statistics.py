from sklearn.metrics import roc_curve, roc_auc_score
import csv
import matplotlib.pyplot as plt
from os.path import splitext

def plot_ROC(n_actives,n_decoys, title, fname, results,
             legend, show, glidefile):
    colors = ['b','g','b','g','r','c','m','y','k']

    results_file=results[0::3]
    results_type=results[1::3]
    results_name=results[2::3]
    
    ys=[]
    scores=[]
    data_types=[]
    glide_y, glide_score, _ = GetYScoreFromResult(glidefile, 'dock', [],[])

    
    for i in range(len(results_file)):
        y,score,t = GetYScoreFromResult(results_file[i], results_type[i],
                                      glide_y, glide_score)
        ys.append(y)
        scores.append(score)
        data_types.append(t)

    plt.figure(figsize=(8,8), dpi=150)
    SetupROCCurvePlot(plt,title)

    aucs = []
    ef10s = []
    ef1s = []
    for i in range(len(results_file)):
        auc, ef10, ef1 = AddROCCurve(plt, ys[i], scores[i], colors[i],
                         results_name[i], n_actives, n_decoys, data_types[i])
        aucs.append(auc)
        ef10s.append(ef10)
        ef1s.append(ef1)


    SaveROCCurvePlot(plt, fname, show, legend, True)
    SaveAUCEF(fname, results, aucs, ef10s, ef1s)



def SaveAUCEF(fname, results, aucs, ef10s, ef1s):

    results_file=results[0::3]
    results_type=results[1::3]
    results_name=results[2::3]
    
    glide_auc = proposed_auc = 0
    with open(splitext(fname)[0]+"_result.txt","w") as f_result:
        f_result.write(str(title)+"\n")
        for i in range(len(results_file)):
            f_result.write(results_name[i]+"_AUC, "+str(aucs[i])+"\n")
            f_result.write(results_name[i]+"_EF10, "+str(ef10s[i])+"\n")
            f_result.write(results_name[i]+"_EF1, "+str(ef1s[i])+"\n")

            if "Glide" in results_name[i]:
                glide_auc = aucs[i]
                
            elif "proposed" in results_name[i]:
                proposed_auc = aucs[i]

        if proposed_auc==0 or glide_auc==0:
            pass
        elif proposed_auc - glide_auc > 0.001:
            f_result.write("proposed, win\n")
        elif proposed_auc - glide_auc < -0.001:
            f_result.write("proposed, lose\n")
        else:
            f_result.write("proposed, draw\n")
            

def GetYScoreFromResult(filename,datatype,glide_y, glide_score):
    y=[]
    score=[]

    try:
        data = csv.reader(open(filename, 'rb'), delimiter=',', quotechar='#')
    except IOError:
        print("postprocess was cancelled. auc is the same as glide SP.")
        return glide_y, glide_score, "dock"
    
    if 'dock' in datatype:
        #print(filename)
        for line in data:
            if "ZINC" in line[0] or "TEMP" in line[0]:
                y.append(0)
            else:
                y.append(1)
            
            score.append(float(line[1]))
    
    else:
    #datatype=='stat':
        for line in data:
            y.append(int(line[3]))
            score.append(float(line[2]))

    return y, score, datatype

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
    plt.title("ROC Curve ("+title+")", fontsize=14)
    
def AddROCCurve(plt, actives, scores, color, label, n_actives, n_decoys, type):
    tpr, fpr = GetRates(actives, scores, n_actives, n_decoys)
    #print(actives,label)
    if "dock" in type:
        scores = [-x for x in scores] #reverse order
    #print(tpr)

    auc_tmp = roc_auc_score(actives, scores)
    auc_tmp = (auc_tmp*tpr[-1])#+((1-fpr[-1])*tpr[-1]) #adjust final position
    auc = round(auc_tmp,3)

    ef_10 = tpr[len(tpr)//10]*10
    ef_1  = tpr[len(tpr)//100]*100
    label = label+", auc="+str(auc)

    if "Glide" in label:
        plt.plot(fpr, tpr, linestyle='dashed', color=color,
                 linewidth=2, label=label)
    else:
        plt.plot(fpr, tpr, color=color, linewidth=2, label=label)
    return auc, ef_10, ef_1

def SaveROCCurvePlot(plt, fname, show, legend, randomline=True):

    if randomline:
        x = [0.0, 1.0]
        plt.plot(x, x, linestyle='dashed', color='red',
                 linewidth=2, label='random')

    plt.xlim(0.0, 1.0)
    plt.ylim(0.0, 1.0)
    if legend:
        plt.legend(fontsize=10, loc='best')

    plt.tight_layout()
    plt.savefig(fname)
    if show:
        plt.show()
    

if __name__ == '__main__':
    import sys
    x = sys.argv
    if len(x)<=5 and len(x)%3!=2:
        print("usage: n_actives, n_decoys, graph_title, graph_filename,"
        "legend, show, glidefile, results(file, type, name)...")
    
    n_actives = int(x[1])
    n_decoys = int(x[2])
    title = x[3]
    fname = x[4]
    legend = x[5]
    if legend in ["False", "false", "0", "No", "no"]:
        legend = False
    else:
        legend = True

    show = x[6]
    if show in ["False", "false", "0", "No", "no"]:
        show = False
    else:
        show = True

    glidefile = x[7]
    results = x[8:]

    plot_ROC(n_actives,n_decoys,title,fname,results,legend,show,glidefile)
    

