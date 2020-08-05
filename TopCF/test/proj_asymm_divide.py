#! /usr/bin/env python
"""\
%(prog)s [-o outfile] <yodafile>
  Create asymmetry histograms (normal and Marzani) from input YODA file. 
  Takes difference of two histograms and plots as 2D Scatter. 
"""

import yoda, argparse, sys, math

parser = argparse.ArgumentParser()
parser.add_argument('inputfile', help="input yodafile")
parser.add_argument("-o", "--output", default="-", dest="OUTPUT_FILE", metavar="PATH",
                    help="write output to specified path")
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
        
s_perp_asymm = {}
for k1, v1 in h_perp_pos.items():
    for k2, v2 in h_perp_neg.items(): ##v1, v2 are members of class of 1D Histograms
        path1 = v1.annotation("Path")
        path2 = v2.annotation("Path")
        path1 = path1.replace("_pos", "")
        path2 = path2.replace("_neg", "")
        if path1 == path2:
            #v1.normalize()
            #v2.normalize()
            a = v1 - v2 ##subtracts bin areas (<=> bin heights)
            b = v1 + v2 
            patha = v1.annotation("Path")
            patha = patha.replace("pullvectorj1_perpproj_all_signed_pos", "perpproj_asymm")
            A = a.divideBy(b) ##creates asymmetry 2D scatter
            A.setAnnotation("Path", patha)
            s_perp_asymm[patha] = A

s_para_asymm = {}
for k1, v1 in h_para_pos.items():
    for k2, v2 in h_para_neg.items(): ##v1, v2 are members of class of 1D Histograms
        path1 = v1.annotation("Path")
        path2 = v2.annotation("Path")
        path1 = path1.replace("_pos", "")
        path2 = path2.replace("_neg", "")
        if path1 == path2:
            #v1.normalize()
            #v2.normalize()
            a = v1 - v2
            b = v1 + v2
            patha = v1.annotation("Path")
            patha = patha.replace("pullvectorj1_paraproj_all_signed_pos", "paraproj_asymm")
            A = a.divideBy(b)
            A.setAnnotation("Path", patha)
            s_para_asymm[patha] = A

s_perp_asymm.update(s_para_asymm)
s_asymm = s_perp_asymm
      
## Write output
yoda.write(s_asymm, args.OUTPUT_FILE)