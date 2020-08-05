#!/usr/bin/env python

import sys, os
if len(sys.argv) < 3:
    print "Usage: %s INPUT.yoda OUTPUT.yoda" % sys.argv[0]
    print "\nThe INPUT.yoda must contain CMS_2015_I1370682 and CMS_2015_I1370682_internal"
    print "to extract correction functions and the correction function should not be applied"
    sys.exit()

from math import *
import yoda
src = yoda.readYODA(sys.argv[1])

hists = {}
## Get histogram pairs
for k in src.viewkeys():
    module, hName = k.split('/')[-2:]
    if 'CMS_2015_I1370682' not in module: continue

    if hName not in hists: hists[hName] = []

    hists[hName].append(src[k])

## Get correction factors
corrFactors = {}
for k in hists.keys():
    if len(hists[k]) != 2: continue
    h1, h2 = hists[k]
    # put particle level front, parton level latter.
    if 'internal' in h1.path: h2, h1 = h1, h2

    cfs = []
    for i, b in enumerate(h1.bins):
        #x = (b.xEdges[0]+b.xEdges[1])/2
        #xErr = (b.xEdges[1]-b.xEdges[0])/2
        y1, y1Err = h1.bins[i].height, h1.bins[i].heightErr
        y2, y2Err = h2.bins[i].height, h2.bins[i].heightErr

        r, rErr = 1, 1000
        if y1 > 0: r = y2/y1 ## We are calculating particle -> parton in full PS
        if y1 > 0 and y2 > 0: rErr = r*sqrt((y1Err/y1)**2+(y2Err/y2)**2)

        cfs.append([r, rErr])
    corrFactors[h1.path] = cfs[:]

## Print out correction functions in C++ code format
for path in sorted(corrFactors.keys()):
    hName = path.split('/')[-1]

    print '//', hName,
    print '= {',
    for x in corrFactors[path]: print ("%f," % x[0]),
    print '};'

## Apply correction functions and write output yoda file
out = {}
for k in sorted(src.viewkeys()):
    if k not in corrFactors:
        out[k] = src[k]
    else:
        f = corrFactors[k]
        #out[k] = src[k].clone()
        out[k] = yoda.Scatter2D(src[k].path, src[k].title)
        #out[k].reset()
        for i, b in enumerate(src[k].bins):
            x = (b.xEdges[0]+b.xEdges[1])/2
            xErr = (b.xEdges[1]-b.xEdges[0])/2
            sf = f[i][0]
            #out[k].addPoint(x, b.area*f[i], xErr, b.areaErr*f[i])
            out[k].addPoint(x, b.height*sf, xErr, b.heightErr*sf)
            #out[k].fillBin(i, f[i]*b.sumW)

yoda.writeYODA(out, sys.argv[2])
