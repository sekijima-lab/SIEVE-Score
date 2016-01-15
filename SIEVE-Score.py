from os.path import splitext
import numpy as np

from read_interaction import read_interaction
from scoring import scoring

def SIEVE(x):
    inter_array = np.array(read_interaction(x['i'],x['hits']))
    interactions = np.array(inter_array[1:,2:],dtype='float')

    #print(inter_array)

    sim = scoring(inter_array,x['o'],x['p'],x['m'],x['propose'],
                        x['cutoff'],x['zeroneg'],
                        x['score_correction'],x['threshold'])

    print('\n*****Process Complete.*****\n')
    with open(splitext(x['o'])[0]+'.log','a') as f_log:
        f_log.write('\n*****Process Complete.*****\n')


if __name__ == '__main__':
    import sys
    from options import Input_func as Input

    #options.py
    default_option = {'i': [], 'o': None, 'hits': 'hits.txt',
                      'cutoff': 1,'threshold': 1.0,
                      'zeroneg':False, 'score_correction':False,
                      'p': 1, 'm': -1,'p_opt': False,
                      'actives':None, 'decoys':None, 'propose':100000
                      }

    if "-file" in sys.argv:
        option_file = sys.argv[sys.argv.index("-file") + 1]
        o = open(option_file, 'r').readlines()
        o = [s.replace("\n","") for s in o]
        o = o + sys.argv
    else:
        o = sys.argv

    option = Input(default_option,o)

    x = option
    with open(splitext(x['o'])[0]+'.log','w') as f_log:
        f_log.write('options:\n'+str(x)+'\n')

    SIEVE(x)
