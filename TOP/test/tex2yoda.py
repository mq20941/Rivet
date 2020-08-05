#!/usr/bin/env python
#TEX input converter for TOP-12-028 paper

BASEDIR='papers/TOP-12-028/arxiv-0/'
INPUT=BASEDIR+'supplemental_material.tex'

files = []
for l in open(INPUT).readlines():
    l = l.strip()
    if not l.startswith('\\input{figures'): continue
    files.append(l[7:-1])

for i, f in enumerate(files):
    name = f.split('/')[-1].split('.')[0]
    lines = open(BASEDIR+f).readlines()
    doBinCenter = len(lines[0][:-2].split('&')) > 5

    print "# BEGIN YODA_SCATTER2D /REF/CMS_TOP_12_028/h%02d_%s" % (i, name)
    print "Path=/REF/CMS_TOP_12_028/h%02d_%s" % (i, name)
    print "Title=%s" % name
    print "Type=Scatter2D"
    print "# xval   xerr-   xerr+   yval  yerr-   yerr+"
    for l in lines:
        l = l.strip().replace('$', '')
        item = [x.strip() for x in l[:-2].split('&')]
        xlow, xhigh = [float(x) for x in item[0][1:-1].split(',')]
        yval = item[-4].replace('\\cdot', '*').replace('{', '').replace('}', '')
        yval = yval.replace('10^', '1E')
        yval = eval(yval)
        stat = float(item[-3])
        syst = float(item[-2])
        total = float(item[-1])

        if doBinCenter: xval = float(item[1].split()[0])
        else: xval = (xlow+xhigh)/2
        err = yval*total/100

        if doBinCenter: print xval, (xlow+xhigh)/2

        print "%4.10e %4.10e %4.10e %4.10e %4.10e %4.10e" % (xval, xval-xlow, xhigh-xval, yval, err, err)
    print "# END YODA_SCATTER2D"
    print "\n"

