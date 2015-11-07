import numpy as np
import schrodinger.structure as structure

def score_func(a, b, coef):
    dist=np.linalg.norm(a-b, 2)
    coed = float(coef)
    if dist < 1:
        return coef
    else:
        return coef / dist

def scoring(result,outputfile,plus,minus,propose,zeroneg=False,correction=False):
    data  = result[1:,:]
    title = data[:,0]
    ishit = data[:,1].astype('int')
    xyz   = data[:,3:6].astype('float')   

    hits  = np.array([True if x>0 else False for x in ishit])
    if zeroneg==False:
        nonhits=np.array([True if x<0 else False for x in ishit])
    elif zeroneg==True:
        nonhits=np.array([True if x<=0 else False for x in ishit])
    #<= for DUD-E set, < for others

    hits_xyz = data[hits,3:6].astype('float')
    numhits = ishit[hits].shape[0]
    nonhits_xyz = data[nonhits,3:6].astype('float')
    known_xyz = np.r_[hits_xyz, nonhits_xyz]

    print(numhits,plus,minus,known_xyz.shape)

    score = np.zeros(data.shape[0])
    if known_xyz != []: 
        for i in range(data.shape[0]):
            for j in range(len(known_xyz)):
                if j<numhits:
                    load = plus
                else:
                    load = minus
                
                score[i] = score[i] + score_func(xyz[i],known_xyz[j],load)

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

    return score

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
