from os.path import splitext
import numpy as np
from schrodinger import structure

def save_to_csv(outputfile, PCA_model, pca, inter_array, labels):
    #save loadings
    loadings = PCA_model.components_.tolist()
    for i in range(3):
        loadings[i].insert(0,'PC'+str(i))
    loadings.insert(0, inter_array[0,2:].tolist())
    
    loadings[0].insert(0, 'residues')
    loadings = np.array(loadings, dtype='str')
    loadings = np.transpose(loadings)
    loadingfile = splitext(outputfile)[0]+'_loading.csv'
    np.savetxt(loadingfile, loadings, fmt="%s", delimiter=",")
    with open(splitext(outputfile)[0]+'.log','a') as f_log:
        f_log.write('\nSaved loadings to '+loadingfile+'\n')
    print('Saved loadings to '+loadingfile)
 
    #save results
    result_title = ['Title']+(inter_array[1:,0].tolist())
    result_ishit = ['ishit']+(inter_array[1:,1].tolist())
    result_label = ['cluster']+(labels.tolist())
    result_x     = ['PC1']+(pca[:,0].tolist())
    result_y     = ['PC2']+(pca[:,1].tolist())
    result_z     = ['PC3']+(pca[:,2].tolist())

    output_pca = splitext(outputfile)[0]+'_cluster.csv'
    result = np.transpose(np.array([result_title,result_ishit,result_label,
                                    result_x, result_y, result_z],dtype='str'))
    np.savetxt(output_pca, result, fmt="%s", delimiter=",")
    print('Saved cluster/PCA info to '+outputfile)
    with open(splitext(outputfile)[0]+'.log','a') as f_log:
        f_log.write('Saved cluster/PCA info to '+outputfile+'\n')

    return result
