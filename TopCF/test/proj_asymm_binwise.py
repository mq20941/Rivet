#! /usr/bin/env python
"""\
%(prog)s [-o outfile] <yodafile>
  Create asymmetry histograms from input YODA file. 
  Takes difference of histograms in each bin and weights by x coord. 
"""

import yoda, argparse, sys, math

parser = argparse.ArgumentParser()
parser.add_argument('inputfile', help="input yodafile")
parser.add_argument("-o", "--output", default="-", dest="OUTPUT_FILE", metavar="PATH",
                    help="write output to specified path")
parser.add_argument("-n", "--normalize", dest="NORMALIZE_HISTOS", action='store_true',
                    help="normalize asymm histograms")
parser.add_argument("-N", "--no-normalize", dest="NORMALIZE_HISTOS", action='store_false',
                    help="do not normalize asymm histograms, default setting")
parser.set_defaults(NORMALIZE_HISTOS=False)
args = parser.parse_args()

## Put the incoming objects into a dict from each path to a list of 1D Histograms 
input = args.inputfile
aos = yoda.read(input)
h_perp_pos = {} ##create empty dict for each histogram type 
h_perp_neg = {}
h_para_pos = {}
h_para_neg = {}
for aopath, ao in aos.items(): ##items=keys:values
    path = ao.annotation("Path")
    if type(ao) is yoda.core.Histo1D:
        if "pullvectorj1_perpproj_all_signed_pos" in path:
                h_perp_pos[aopath] = ao
        if "pullvectorj1_perpproj_all_signed_neg" in path:
                h_perp_neg[aopath] = ao
        if "pullvectorj1_paraproj_all_signed_pos" in path:
                h_para_pos[aopath] = ao
        if "pullvectorj1_paraproj_all_signed_neg" in path:
                h_para_neg[aopath] = ao

h_perp_asymm = {}
for k1, v1 in h_perp_pos.items():
    for k2, v2 in h_perp_neg.items(): ##v1, v2 are members of class of 1D Histograms
        path1 = v1.annotation("Path")
        path2 = v2.annotation("Path")
        path1 = path1.replace("_pos", "")
        path2 = path2.replace("_neg", "")
        if path1 == path2:
            numBins = v1.numBins()
            xMax = v1.xMax()
            xMin = v1.xMin()
            path = v1.annotation("Path")
            path = path.replace("pullvectorj1_perpproj_all_signed_pos", "perpproj_asymm")
            a = yoda.Histo1D(numBins, xMin, xMax, path, "")
            edges = a.xEdges()
            for i in range(0, numBins-1):
                pos_bin = v1.bin(i)
                neg_bin = v2.bin(i)
                binwidth = edges[i+1] - edges[i]
                diff = pos_bin.height() - neg_bin.height()
                a.fillBin(i, diff*binwidth) ##when histogram filled with weight w it gives that bin area w, not height w
                #asymm = diff*(edges[i+1]/2)
                #a.fillBin(i, asymm)
            if args.NORMALIZE_HISTOS == True:
                if a.integral > 0:
                    a.normalize() 
            h_perp_asymm[path] = a

h_para_asymm = {}
for k1, v1 in h_para_pos.items():
    for k2, v2 in h_para_neg.items(): ##v1, v2 are members of class of 1D Histograms
        path1 = v1.annotation("Path")
        path2 = v2.annotation("Path")
        path1 = path1.replace("_pos", "")
        path2 = path2.replace("_neg", "")
        if path1 == path2:
            numBins = v1.numBins()
            xMax = v1.xMax()
            xMin = v1.xMin()
            path = v1.annotation("Path")
            path = path.replace("pullvectorj1_paraproj_all_signed_pos", "paraproj_asymm")
            a = yoda.Histo1D(numBins, xMin, xMax, path, "")
            edges = a.xEdges()
            for i in range(0, numBins-1):
                pos_bin = v1.bin(i)
                neg_bin = v2.bin(i)
                binwidth = edges[i+1] - edges[i]
                diff = pos_bin.height() - neg_bin.height()
                a.fillBin(i, diff*binwidth)
                #asymm = diff*(edges[i+1]/2)
                #a.fillBin(i, asymm)
            if args.NORMALIZE_HISTOS == True:
                if a.integral > 0:
                    a.normalize()  
            h_para_asymm[path] = a

h_perp_asymm.update(h_para_asymm)
h_asymm = h_perp_asymm

## Write output
yoda.write(h_asymm, args.OUTPUT_FILE)
