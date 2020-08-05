#!/usr/bin/python

import sys, getopt
import os.path

def main(argv):
    inputanalysis = ''
    outputanalysis = ''
    try:
        opts, args = getopt.getopt(argv,"hi:o:",["ianalysis=","oanalysis="])
    except getopt.GetoptError:
        print('test.py -i <inputanalysis> -o <outputanalysis>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('test.py -i <inputanalysis> -o <outputanalysis>')
            sys.exit()
        elif opt in ("-i", "--ianalysis"):
            inputanalysis = arg
        elif opt in ("-o", "--oanalysis"):
            outputanalysis = arg
    if (inputanalysis == '' and outputanalysis == ''):
        print('test.py -i <inputanalysis> -o <outputanalysis>')
        sys.exit(2)
    print('Input analysis  = ' + inputanalysis)
    print('Output analysis = ' + outputanalysis)
    
    if (inputanalysis == outputanalysis):
        print 'Input and output must not be the same!'
        sys.exit(2)
    
    f = open('src/'+inputanalysis+'.cc', 'r')
    w = open('src/'+outputanalysis+'.cc', 'w')
    
    for line in f:
        w.write(line.replace(inputanalysis, outputanalysis))
    
    f.close()
    w.close()
    
    datafiles = ['.info', '.plot', '.yoda']
    
    for d in datafiles:
        if not os.path.isfile('data/'+inputanalysis+d): continue
        f = open('data/'+inputanalysis+d, 'r')
        w = open('data/'+outputanalysis+d, 'w')
        
        for line in f:
            w.write(line.replace(inputanalysis, outputanalysis))
        
        f.close()
        w.close()


if __name__ == "__main__":
   main(sys.argv[1:])
