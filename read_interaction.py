from schrodinger import structure
import numpy as np

def read_interaction(inputfile,hitsfile):
    reader = structure.StructureReader(inputfile)

    reader.next() #skip protein
    st = reader.next()

    residues = []
    result = []
    legend = ['# title', 'ishit']


    for x in st.property.keys():
        if x.startswith('r_glide_res:') and\
           x.endswith('_Eint'):
            res = x.replace('r_glide_res:','').replace('_Eint','')
            residues.append(res)

    residues = sorted(residues)
    [legend.append(i) for i in residues]
    result.append(legend)

    reader.close()  

    try:
        hits = np.loadtxt(hitsfile, delimiter=',', comments='#',
               dtype={'names': ('title', 'ishit'), 'formats': ('S16', '<i2')})
        if hits.shape == ():
            hits = np.array([hits])
    except IOError:
        hits = None


    reader = structure.StructureReader(inputfile) #rewind iter
    for st in reader:
        prop = st.property
        inter=[]

        if 'r_i_docking_score' not in prop.keys(): #protein? error?
            continue
        
        inter.append(prop['s_m_title'])

        if hits != None and prop['s_m_title'] in hits['title']:
            place = np.where(hits['title']==prop['s_m_title'])
            i=int(place[0])            
            inter.append(hits[i][1])
        else:
            inter.append(0)

        for x in residues:
            inter.append(float(prop['r_glide_res:'+x+'_Eint']))

        result.append(inter)
        
            
    reader.close()
    return result
