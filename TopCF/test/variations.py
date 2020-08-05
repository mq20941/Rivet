#!/usr/bin/env python

"""\
%(prog)s [options] <yodafile> <txt file>
Produce yodafile for all histograms with a given variation within 
an input yodafile, produced by default in wd. 
Variations are specified in txt file in wd.
"""
import yoda, argparse, shutil, os, os.path
from os import path

parser = argparse.ArgumentParser(description='Produce variation histograms in new yodafile from input yodafile') #input yodafile from commandline
parser.add_argument('inputfile', help="input yodafile")
parser.add_argument('variations', help="input txt file containing variations")
parser.add_argument('-o', '--outputdirectory', default="", help="directory for yodafile output")
args = parser.parse_args()

directory = args.outputdirectory
nwd = "./" + directory

if not path.exists(nwd) : 
     os.mkdir(nwd)

input = args.inputfile
aos = yoda.read(input) #reads input yodafile
max = len(aos) #number of objects in input yodafile

variations = args.variations
vars = open(variations, "r") #vars is a file


#obtain variations in list
nvars = 0 
variation = []
for x in vars : 
     nvars += 1
     x = x.strip()
     variation.append(x)

for j in range(0,nvars) : #runs over each variation, loops over input yodafile nvars times
     v = variation[j]
     o = v.replace('=', '_').replace('.', '') #creates output yodafile names, dependent on input data model
     o = "CP5_" + o + ".yoda"
     hvariation = {} #create empty dict
     for i in range(0,max) : #loops over all hists in input yodafile
          h = aos.values()[i] #aos.values is a list 
          hpath = h.annotation("Path")
          if v in hpath : 
               y = "[" + v + "]"
               hpath = hpath.replace(y, '')
               h.setAnnotation("Path", hpath)
               hvariation[hpath] = h
           
     yoda.write(hvariation, o)
     cwd = "./" + o
     fwd = nwd + "/" + o
     shutil.move(cwd, fwd)

vars.close()