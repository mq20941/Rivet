#! /usr/bin/env python
"""\
%(prog)s [-o outfile] <yodafile>
  Calculate correlation coefficient for variables in specified 1D Profile Histogram from input YODA file. 
"""

import yoda, argparse, sys, math
from scipy.stats import pearsonr

parser = argparse.ArgumentParser()
parser.add_argument('inputfile', help="input yodafile")
#parser.add_argument("-o", "--output", default="-", dest="OUTPUT_FILE", metavar="PATH",
                    #help="write output to specified path")
parser.add_argument('plot')
args = parser.parse_args()

input = args.inputfile
aos = yoda.read(input)

cor_plot = args.plot

p_cor = {}
for aopath, ao in aos.items(): 
    path = ao.annotation("Path")
    if type(ao) is yoda.core.Profile1D:
        if cor_plot in path and "[" not in path and "RAW" not in path:
                p_cor[aopath] = ao

for k, v in p_cor.items():
    V = v.mkScatter()
    points = V.points()
    cor_x = []
    cor_y = []
    for p in points:
        cor_x.append(p.x())
        cor_y.append(p.y())
    corr, pval = pearsonr(cor_x, cor_y)
    plot = "Plot:" + k
    print plot
    print('Pearsons correlation: %.3f' % corr)
    print('P-value: %.3f' % pval)