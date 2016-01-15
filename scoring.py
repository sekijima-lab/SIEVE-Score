import numpy as np

def score_func(a, b, coef, cutoff):
    dist = np.linalg.norm(a-b, 2)
    return coef / max(dist, cutoff)

def scoring(inter_array,outputfile,plus,minus,propose,
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
    nonhits_xyz = data[nonhits,3:].astype('float')
    known_xyz = np.r_[hits_xyz, nonhits_xyz]

#calculate score
        
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
    result = np.dstack((saved, saved_comp,
                        saved_score, saved_ishit))[0].astype('str')
    print(result)

    if known_xyz != []: 
        out_rank = outputfile
        np.savetxt(out_rank, result, fmt="%s", delimiter=",")
        print('Saved SIEVE-Score.')

    return 1

def save_to_maegz(inputfile,labels,score):
    import schrodinger.structure as structure
    
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
        index+=1

    reader.close()
    writer.close()
    print('Saved Cluster info to '+write_maegz)
