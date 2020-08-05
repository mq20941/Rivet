# Setup instructions

```
# Work in /storage/<username>
cmsrel CMSSW_11_1_0_pre8
cd CMSSW_11_1_0_pre8/src
cmsenv

git clone ssh://git@gitlab.cern.ch:7999/ejclemen/Rivet.git -b TopCF

source Rivet/rivetSetup.sh

scram b -j 8

cd Rivet/TopCF/test

../../../condor-comforter/cmsRunCondor/cmsRunCondor.py --filelist=CP5.txt --inputFile=../../rivetSetup.sh --unitsPerJob=1 --eventsPerJob=-1 --totalNumJobs=20 --outputDir=/hdfs/user/mq20941/TopCF/CP5/ --mc config.py

#  To check progress of jobs, do:
condor_q

# If it works, you will see yoda files appearing in /hdfs/user/mq20941/TopCF/CP5/ as they complete
# Combine with:
rivet-merge -o MC.yoda /hdfs/user/mq20941/TopCF/CP5/*yoda flag -o specifies output histo file path (default = Rivet.yoda) 

# Plot combined file with MC.yoda as usual
rivet-mkhtml --mc-errs --no-weights MC.yoda

```

If these ran, you should be in /storage/<username>/CMSSW_11_1_0_pre8/src/Rivet/TopCF/test
and see the files test.yoda and rivet-plots/index.html.


# Official README below - can ignore these

# CMS Rivet repository

This GitLab repository contains Rivet routines for CMS analyses and sample configuration files to run them.

Please consult the README files in the subdirectories for information about the plugins contained there.

Please follow the [contribution guide](CONTRIBUTING.md) for developing your plugins.

## Installation

First, create a personal fork of this repository: https://gitlab.cern.ch/cms-gen/Rivet/forks/new

    cmsrel CMSSW_11_0_1
    cd CMSSW_11_0_1/src
    cmsenv

    git clone ssh://git@gitlab.cern.ch:7999/${USER}/Rivet.git
    cd Rivet
    git remote add cms-gen ssh://git@gitlab.cern.ch:7999/cms-gen/Rivet.git
    git pull cms-gen master

    source rivetSetup.sh
    scram b -j8

