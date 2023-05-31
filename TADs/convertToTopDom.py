#!/home/gdstantonlab/mxw010/.conda/envs/hic/bin/python

import cooler
import numpy as np
import scipy.sparse as sparse
import re
import pandas as pd
from scipy.sparse import save_npz
import argparse
  
parser = argparse.ArgumentParser(description='convert to TopDom format')
parser.add_argument('--chr', type=str, nargs=1,
                    help='chromosome to convert')
parser.add_argument('--cell', type=str, nargs=1,
                    help = "cell line")
parser.add_argument('--res', type= int, nargs =1,
                    help="resolution",
                    default=50000)
args = parser.parse_args()
chr = args.chr[0]
cell = args.cell[0]
resolution = args.res[0]

c = cooler.Cooler('/home/gdstantonlab/mxw010/Data/SPICE-C/' + cell + '/cools/' + cell + '.mcool::/resolutions/' + str(resolution))
#python index [a:b] goes from a to b-1
#mat = c.matrix(sparse=True)[:]
#save_npz("freq.npz", mat)

if resolution < int(1e6): 
    res = str(int(resolution / 1000)) + 'kb'
else: 
    res = str(int(resolution / 1000000)) + 'mb'

print("writing matrix of "  + chr + " at resolution " + res + "...")
mat = c.matrix().fetch(chr)
np.savetxt("./data_" + res + "/freq_" + chr, mat, delimiter="\t")

print("writing chr information of " + chr + "...")
#bed information columns
info = c.bins().fetch(chr)[['chrom', 'start', 'end']][:]
info.to_csv("./data_" + res + "/chrom_" + chr, sep="\t", header=False,index=False)                 
