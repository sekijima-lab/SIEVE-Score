from os.path import splitext
import numpy as np
from schrodinger import structure

def save_to_csv_noPCA(outputfile, interactions, inter_array, labels):
    #save results
    result_title = inter_array[1:,0].tolist()
    result_ishit = inter_array[1:,1].tolist()
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
