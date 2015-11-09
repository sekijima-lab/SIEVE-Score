import numpy as np
import schrodinger.structure as structure
import csv

def score_func(a, b, coef):
    dist = np.linalg.norm(a-b, 2)
    if dist < 1:
        return coef
    else:
        return coef / dist

def scoring_noPCA_hist(inter_array,outputfile,plus,minus,propose,zeroneg=False,correction=False):
    data  = inter_array[1:,:]
    title = data[:,0]
    ishit = data[:,1].astype('int')
    xyz   = data[:,3:].astype('float')   

    hits  = np.array([True if x>0 else False for x in ishit])
    if zeroneg==False:
        nonhits=np.array([True if x<0 else False for x in ishit])
    elif zeroneg==True:
        nonhits=np.array([True if x<=0 else False for x in ishit])
    #<= for DUD-E set, < for others

    hits_xyz = data[hits,3:].astype('float')
    numhits = ishit[hits].shape[0]
    nonhits_xyz = data[nonhits,3:].astype('float')
    known_xyz = np.r_[hits_xyz, nonhits_xyz]

    print(numhits,plus,minus,known_xyz.shape)

    score = np.zeros(data.shape[0])

    """if known_xyz != []: 
        for i in range(data.shape[0]):
            for j in range(len(known_xyz)):
                if j<numhits:
                    load = plus
                else:
                    load = minus
    """
    pp = [0]*20
    pn = [0]*20
    nn = [0]*20
    mk_index = lambda i,j: min(int(np.linalg.norm(xyz[i]-xyz[j], 2)),19)

    for i in range(data.shape[0]):
        for j in range(i+1,data.shape[0]):
            if ishit[i]>0 and ishit[j]>0:
                pp[mk_index(i,j)]+=1
            elif ((ishit[i]>0 and ishit[j]<=0) or
                  (ishit[i]<=0 and ishit[j]>0)):
                pn[mk_index(i,j)]+=1
            elif ishit[i]<=0 and ishit[j]<=0:
                nn[mk_index(i,j)]+=1

    pp = [float(x)/sum(pp) for x in pp]
    pn = [float(x)/sum(pn) for x in pn]
    nn = [float(x)/sum(nn) for x in nn]

    """
    X = [i+0.5 for i in range(20)]
    plt.bar(x,pp, label="active-active", color = "blue")
    plt.bar(x,pn, label="active-inactive", color = "red")
    plt.xlabel("Distance",fontsize=14)
    plt.ylabel("Percentage",fontsize=14)
    plt.xticks([i for i in range(20)])
    plt.title(outputfile+"_active")
    plt.savefig(outputfile+"_a.eps")

    plt.clf()
    plt.bar(x,nn, label="inactive-inactive", color = "blue")
    plt.bar(x,pn, label="inactive-active", color = "red")
    plt.xlabel("Distance",fontsize=14)
    plt.ylabel("Percentage",fontsize=14)
    plt.xticks([i for i in range(20)])
    plt.title(outputfile+"_inactive")
    plt.savefig(outputfile+"_i.eps")
    """
    with open(outputfile, 'w') as f_out:
        csvWriter = csv.writer(f_out)
        csvWriter.writerow(pp)
        csvWriter.writerow(pn)
        csvWriter.writerow(nn)

