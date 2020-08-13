#!/usr/bin/env python

import yoda, argparse, math

parser = argparse.ArgumentParser(description='Create ROC scatter from input yodafile') #input yodafile from commandline
parser.add_argument('inputfile', help="input yodafile")

args = parser.parse_args()

input = args.inputfile
aos = yoda.read(input) 
max = len(aos) #number of objects in input yodafile

s1_PT = {}
s2_PT = {}
s1_TP = {}
s2_TP = {}

for i in range(0,max):
    h = aos.values()[i] 
    hpath = h.annotation("Path")
    htype = h.annotation("Type")
    if "Profile1D" in htype: 
        if "PT_s1" in hpath and "[" not in hpath and "RAW" not in hpath:
            s1_PT[hpath] = h
        if "TP_s1" in hpath and "[" not in hpath and "RAW" not in hpath:
            s1_TP[hpath] = h
        if "PT_s2" in hpath and "[" not in hpath and "RAW" not in hpath:
            s2_PT[hpath] = h
        if "TP_s2" in hpath and "[" not in hpath and "RAW" not in hpath:
            s2_TP[hpath] = h


for k1, v1 in s1_TP.items():
    for k2, v2 in s1_PT.items():
        path1 = v1.annotation("Path")
        path2 = v2.annotation("Path")
        path1 = path1.replace("TP_", "")
        path2 = path2.replace("PT_", "")
        if path1 == path2:
            V1 = v1.mkScatter()
            xpoints = V1.points()
            ROCx = []
            for p in xpoints:
                ROCx.append(p.y())
            V2 = v2.mkScatter()
            ypoints = V2.points()
            ROCy = []
            for p in ypoints:
                ROCy.append(p.y())
            npoints = len(ROCx)
            s_ROC_s1 = yoda.core.Scatter2D("/MC_ColourReconstruction/ROC_s1")
            for j in range(0, npoints):
                if not math.isnan(ROCx[j]) and not math.isnan(ROCy[j]):
                    s_ROC_s1.addPoint(ROCx[j], ROCy[j])
            file = open("ROC.txt","w")
            sx = []
            sy = []
            for p in s_ROC_s1.points():
                sx.append(p.x())
                sy.append(p.y())
            X = str(sx)
            X = X + "\n"
            Y = str(sy)
            Y = Y + "\n"
            file.write(X)
            file.write(Y)
            file.close()
            
for k1, v1 in s2_TP.items():
    for k2, v2 in s2_PT.items():
        path1 = v1.annotation("Path")
        path2 = v2.annotation("Path")
        path1 = path1.replace("TP_", "")
        path2 = path2.replace("PT_", "")
        if path1 == path2:
            V1 = v1.mkScatter()
            xpoints = V1.points()
            ROCx = []
            for p in xpoints:
                ROCx.append(p.y())
            V2 = v2.mkScatter()
            ypoints = V2.points()
            ROCy = []
            for p in ypoints:
                ROCy.append(p.y())
            npoints = len(ROCx)
            s_ROC_s2 = yoda.core.Scatter2D("/MC_ColourReconstruction/ROC_s2")
            for j in range(0, npoints):
                if not math.isnan(ROCx[j]) and not math.isnan(ROCy[j]):
                    s_ROC_s2.addPoint(ROCx[j], ROCy[j])
            file = open("ROC.txt","a")
            sx = []
            sy = []
            for p in s_ROC_s2.points():
                sx.append(p.x())
                sy.append(p.y())
            X = str(sx)
            X = X + "\n"
            Y = str(sy)
            Y = Y + "\n"
            file.write(X)
            file.write(Y)
            file.close()

            


            
            

