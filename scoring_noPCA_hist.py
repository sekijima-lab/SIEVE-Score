import numpy as np
import schrodinger.structure as structure
import csv
import random

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
    numnonhits = ishit[nonhits].shape[0]
    numdata = ishit.shape[0]
    nonhits_xyz = data[nonhits,3:].astype('float')
    known_xyz = np.r_[hits_xyz, nonhits_xyz]

    print(numhits,plus,minus,known_xyz.shape)

    """
    score = np.zeros(data.shape[0])

    if known_xyz != []: 
        for i in range(data.shape[0]):
            for j in range(len(known_xyz)):
                if j<numhits:
                    load = plus
                else:
                    load = minus
    """
    pp = []
    pn = []
    nn = []
    
    random.seed()
    dist = lambda i,j: np.linalg.norm(known_xyz[i]-known_xyz[j],2)
    ri = random.randint
    for _ in range(10000):
        x=0
        while x<=0:
            x = dist(ri(0,numhits-1),ri(0,numhits-1))
        pp.append(x)
        
        x=0
        while x<=0:
            x = dist(ri(0,numhits-1),ri(numhits,numdata-1))
        pn.append(x)
        
        x=0
        while x<=0:
            x = dist(ri(numhits,numdata-1),ri(numhits,numdata-1))
        nn.append(x)

    with open(outputfile, 'w') as f_out:
        csvWriter = csv.writer(f_out)
        csvWriter.writerow(pp)
        csvWriter.writerow(pn)
        csvWriter.writerow(nn)

