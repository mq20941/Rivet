#! /usr/bin/env python
"""\
%(prog)s [-o outfile] <yodafile1> <yodafile2>...
  e.g. %(prog)s run1.yoda run2.yoda run3.yoda  (unweighted merging of three runs)
  Merge 2D asymmetry scatters from multiple YODA files. 
  The list of points from the various Scatter2D objects are concatenated.
"""

import yoda, argparse, sys, math

parser = argparse.ArgumentParser()
parser.add_argument("INFILES", nargs="+", help="datafile1 datafile2 [...]") #nargs="+", all command line arguments gathered into a list
parser.add_argument("-o", "--output", default="-", dest="OUTPUT_FILE", metavar="PATH",
                    help="write output to specified path")
args = parser.parse_args()

## Put the incoming objects into a dict from each path to a list of 2D scatters
analysisobjects_in = {} #create empty dict
for fa in args.INFILES:
    filename = fa
    aos = yoda.read(filename)
    for aopath, ao in aos.items(): #items = keys:values
        path = ao.annotation("Path")
        if type(ao) is yoda.core.Scatter2D:
            if "asymm" in path:
                analysisobjects_in.setdefault(aopath, []).append(ao)

analysisobjects_out = {}
for p, aos in analysisobjects_in.items():

    aotype = type(aos[0])
    
    ## Make a copy of the first object as the basis for merging (suitable for all Scatter types)
    ao_out = aos[0].clone()

    ## Merge scatters by combining points 
     
    for ao in aos[1:]:
        ao_out.combineWith(ao)

    ## Put the output AO into the output dict
    analysisobjects_out[p] = ao_out

## Write output
yoda.writeYODA(analysisobjects_out, args.OUTPUT_FILE)

