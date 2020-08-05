#!/usr/bin/env python

import yoda, argparse

parser = argparse.ArgumentParser(description='Print mean of histograms from input yodafile') #input yodafile from commandline
parser.add_argument('inputfile', help="input yodafile")
parser.add_argument('histogram')

args = parser.parse_args()

input = args.inputfile
aos = yoda.read(input) 
max = len(aos) #number of objects in input yodafile

hist = args.histogram

for i in range(0,max):
    h = aos.values()[i] 
    hpath = h.annotation("Path")
    htype = h.annotation("Type")
    if "Histo1D" in htype: 
        if hist in hpath:
            print h.xMean() 
            print hpath

