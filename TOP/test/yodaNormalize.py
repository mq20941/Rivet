#!/usr/bin/env python

import sys, os
if len(sys.argv) < 3:
    print "Usage: %s INPUT.yoda OUTPUT.yoda" % sys.argv[0]
    print "\nScale all histograms by unit area"
    sys.exit()

from math import *
import yoda
src = yoda.readYODA(sys.argv[1])

out = {}
for k in src.viewkeys():
    out[k] = src[k]
    if type(out[k]) == yoda.core.Scatter2D:
        sumW = 0.
        for p in out[k].points:
            width = abs(p.xErrs[0])+abs(p.xErrs[1])
            area = width*p.y
            sumW += area
        for p in out[k].points:
            p.scaleY(1./sumW)
    else:
        if out[k].integral() == 0: continue
        out[k].normalize()

yoda.writeYODA(out, sys.argv[2])
