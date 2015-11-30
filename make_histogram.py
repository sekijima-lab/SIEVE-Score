import pylab as plt
import csv
import sys
import math

def make_histogram(inputfile, outputfile):
    with open(inputfile, 'rb') as f:
        reader = csv.reader(f)
        pp, pn, nn = [l for l in reader]
    pp = [float(x) for x in pp]
    pn = [float(x) for x in pn]
    nn = [float(x) for x in nn]

    print(pp)

    _range = (0.0,200.0)
    bins = math.ceil(math.log(len(pp),2))+1
    #bins = 15 when N=10000, Sturges' formula
    normed = True

    plt.figure()
    plt.hist(pp, bins=bins, range=_range, label="active-active", alpha=0.5, normed=True, color="b")
    plt.hist(pn, bins=bins, range=_range, label="active-inactive", alpha=0.5, normed=True, color="r")
    plt.xlabel("Distance",fontsize=14)
    plt.ylabel("Percentage",fontsize=14)
    plt.title(outputfile+"_active")
    plt.savefig(outputfile+"_a.png")
    plt.show()

    plt.clf()
    

    plt.hist(nn, bins=bins, range=_range, label="inactive-active", alpha=0.5, normed=True, color="b")
    plt.hist(pn, bins=bins, range=_range, label="inactive-inactive", alpha=0.5, normed=True, color="r")
    plt.xlabel("Distance",fontsize=14)
    plt.ylabel("Percentage",fontsize=14)
    plt.title(outputfile+"_inactive")
    plt.savefig(outputfile+"_i.png")
    plt.show()

if __name__=='__main__':

    inputfile = sys.argv[1]
    outputfile = sys.argv[2]
    make_histogram(inputfile, outputfile)
