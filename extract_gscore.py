import numpy as np
import pandas as pd

def extract_gscore(args):
    # read interaction
    import os.path
    input_interaction = os.path.splitext(args.input[0])[0] + ".interaction"
    if len(args.input) == 1 and os.path.exists(input_interaction):
        inter_array = np.genfromtxt(input_interaction, comments=None, delimiter=",", dtype=None)
        inter_array = inter_array.tolist()
        inter_array = np.asarray(inter_array)
    else:
        from read_interaction import read_interaction
        inter_array = np.array(read_interaction(args.input, args.hits)).astype(str)
    
    inter_array = pd.DataFrame(inter_array[1:, :], columns=inter_array[0, :].astype(str), dtype=str)
    
    interactions = inter_array.iloc[:, [0, 1, -1]]
    print(interactions)
    interactions.to_csv("glide_score.csv", sep=",", index=None)
    print('\n*****Process Complete.*****\n')



if __name__ == '__main__':
    # options.py
    from options import Input_func as Input

    args = Input()

    extract_gscore(args)
