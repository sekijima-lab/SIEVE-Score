from os.path import splitext
import numpy as np
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA

from read_interaction import read_interaction
from save import save_to_csv
from scoring import scoring

def glide_KMeans_parameter(x):
    #read_interaction.py
    inter_array = np.array(read_interaction(x['i'],x['hits']))
    interactions = np.array(inter_array[1:,2:],dtype='float')

    #PCA
    PCA_model = PCA(n_components=3)
    pca = PCA_model.fit(interactions).transform(interactions)

    #k-means
    kmeans_model = KMeans(n_clusters=x['cl'], random_state=1000)\
                   .fit(interactions)
    labels = kmeans_model.labels_

    #save.py
    result = save_to_csv(x['o'], PCA_model, pca, inter_array, labels)

    #analyze, visualize
   
    scoring(result,x['o'],x['p'],x['m'],x['propose'],
                        x['zeroneg'],x['score_correction'])
        #save_to_maegz(x['i'],labels,score)

    #plot(pca, inter_array, labels, x['cl'], x['o'],
            #x['skip'], x['title'], x['show'])

    print('\n*****Process Complete.*****\n')
    with open(splitext(x['o'])[0]+'.log','a') as f_log:
        f_log.write('\n*****Process Complete.*****\n')

def plot(pca,inter_array,labels,n_clus,outputfile, skip, graph_title, show):
    if show == False:
        import matplotlib
        matplotlib.use('Agg') # for remote machine

    import matplotlib.pyplot as plt

    ishit = np.array(inter_array[1:,1],dtype='int')

    def plot_(x_axis, y_axis):
        for c, clus_num in zip(colors, range(n_clus)):
            cond = np.multiply(labels==clus_num,ishit == 0)
            disp = np.multiply(cond, skipping)
            if np.any(cond):
                plt.scatter(pca[disp, x_axis],pca[disp, y_axis],\
                            c=c, label='#'+str(clus_num), s=20, marker='^')
        for c, clus_num in zip(colors, range(n_clus)):
            cond = np.multiply(labels==clus_num,ishit < 0)
            if np.any(cond):
                plt.scatter(pca[cond, x_axis],pca[cond, y_axis],\
                            c=c, s=50, marker='x')
        for c, clus_num in zip(colors, range(n_clus)):
            cond = np.multiply(labels==clus_num,ishit > 0)
            if np.any(cond):
                plt.scatter(pca[cond, x_axis],pca[cond, y_axis],\
                            c=c, s=50, marker='o')

        plt.title(graph_title+", K-means clustering")
        plt.xlabel("PC"+str(x_axis+1),fontsize=16)
        plt.ylabel("PC"+str(y_axis+1),fontsize=16)
        plt.legend(fontsize=10)

    skipping = [True if i%skip==0 else False for i in range(pca[:,0].size)]
    skipping = np.array(skipping)
    colors = [plt.cm.rainbow(i/float(n_clus), 1) for i in range(n_clus)]

    plt.figure(figsize=(12,9))

    plt.subplot(221)
    plot_(0,1)

    plt.subplot(222)
    plot_(0,2)

    plt.subplot(223)
    plot_(1,2)

    outputname = splitext(outputfile)[0]+'.png'
    plt.savefig(outputname,dpi=150)
    print('Saved visualized image to '+outputname)

    if show: plt.show()

if __name__ == '__main__':
    import sys
    from options import Input_func as Input

    #options.py
    default_option = {'i': None, 'o': None, 'hits': 'hits.txt',
                      'cl': 5, 'skip': 1, 'title': None,
                      'p': 1, 'm': -1, 'propose': 1000, 'show': True,
                      'score': True, 'zeroneg':False, 'score_correction':False}
    option = Input(default_option,sys.argv)
    x = option
    with open(splitext(x['o'])[0]+'.log','w') as f_log:
        f_log.write('options:\n'+str(x)+'\n')

    glide_KMeans_parameter(x)