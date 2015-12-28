import numpy as np
import schrodinger.structure as structure
from os.path import splitext
from scipy.spatial.distance import cosine
import random
import csv
import math

def score_func(a, b, coef, cutoff):
    dist = np.linalg.norm(a-b, 2)
    if dist < cutoff:
        return coef / cutoff
    else:
        return coef / dist

def scoring_noPCA(inter_array,outputfile,plus,minus,propose,
                  cutoff,zeroneg,correction,threshold):
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

    #print(numhits,plus,minus,known_xyz.shape)
#calculate positive-negative distance histgram
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

    name,ext = splitext(outputfile)
    with open(name+"_hist"+ext, 'w') as f_out:
        csvWriter = csv.writer(f_out)
        csvWriter.writerow(pp)
        csvWriter.writerow(pn)
        csvWriter.writerow(nn)

#calculate cos similarity of pp-pn
    _range = (0.0,200.0)
    _bins = math.ceil(math.log(len(pp),2))+1
        #bins = 15 when N=10000, Sturges' formula
    pp_hist, bins = np.histogram(pp,bins=_bins,range=_range,density=True)
    pn_hist, bins = np.histogram(pn,bins=_bins,range=_range,density=True)
    nn_hist, bins = np.histogram(nn,bins=_bins,range=_range,density=True)

    cos_sim = 1 - cosine(pp_hist, pn_hist)
    with open(splitext(outputfile)[0]+'.log','a') as f_log:
        f_log.write('Similarity of positive-negative was:'+\
                    str(cos_sim))
    
#if cos_similarity > threshold: return without do anything
    if cos_sim > threshold:
        with open(splitext(outputfile)[0]+'.log','a') as f_log:
            f_log.write('Similarity was over the threshold.'+\
                        'We did not change the result of glide.')
            print('Similarity was over the threshold.'+\
                        'We did not change the result of glide.')
        return cos_sim
    
#else: calculate score
        
    score = np.zeros(data.shape[0])
    if known_xyz != []: 
        for i in range(data.shape[0]):
            for j in range(len(known_xyz)):
                if j<numhits:
                    load = plus
                else:
                    load = minus
                
                score[i] = score[i] +\
                    score_func(xyz[i],known_xyz[j],load,cutoff)

    saved = (np.argsort(score)[::-1])[:propose] #reversed,cut
    saved_comp = title[saved]
    saved_score = score[saved]
    saved_ishit = ishit[saved]

    if correction == True:
        for i in range(len(saved_score)):
            if saved_ishit[i]>0:
                saved_score[i]-plus
            elif saved_ishit[i]<0 or (zeroneg==True and saved_ishit[i]==0):
                saved_score[i]-minus
            #<= for DUD-E set, < for others
    print(saved_score)
    result = np.dstack((saved, saved_comp, saved_score, saved_ishit))[0].astype('str')
 
    #save
    """
    if hit_neighbor != []:
        out_hits = outputfile.rstrip('.csv')+'_hits_neighbor.csv'
        np.savetxt(out_hits, hit_neighbor, fmt="%s", delimiter=",")
    if nonhit_neighbor != []:
        out_nonhits = outputfile.rstrip('.csv')+'_nonhits_neighbor.csv'
        np.savetxt(out_nonhits, nonhit_neighbor, fmt="%s", delimiter=",")
    """

    if known_xyz != []: 
        #print('Saved hits/nonhits neighbor.')
        out_rank = outputfile
        np.savetxt(out_rank, result, fmt="%s", delimiter=",")
        print('Saved high-score compounds.')

    return cos_sim

def save_to_maegz(inputfile,labels,score):
    reader = structure.StructureReader(inputfile)
    write_maegz = inputfile.rstrip(".maegz")+"_cluster.maegz"
    writer = structure.StructureWriter(write_maegz)

    index=0
    for st in reader:
        prop = st.property

        if 'r_i_docking_score' not in prop.keys(): #protein? error?
            writer.append(st)
            continue
        else:
            st.property['i_Clustering_Cluster#'] = labels[index]
            st.property['i_Clustering_PCAScore'] = score[index]
            writer.append(st)
        index=index+1

    reader.close()
    writer.close()
    print('Saved Cluster info to '+write_maegz)
