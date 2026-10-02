import pandas as pd

def extract_gscore(args):
    from interaction_io import load_interaction
    names, labels, _, features = load_interaction(args.input, args.hits, args.ignore)
    interactions = pd.DataFrame({"# title": names, "ishit": labels,
                                 "docking_score": features[:, -1]})
    print(interactions)
    interactions.to_csv("glide_score.csv", sep=",", index=None)
    print('\n*****Process Complete.*****\n')



if __name__ == '__main__':
    # options.py
    from options import Input_func as Input

    args = Input()

    extract_gscore(args)
