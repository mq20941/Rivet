## The Les Houches plugin
### Description
This plugin has been developped for the Les Houches ttbar working group.

It fullfills the requirements from:

http://phystev.cnrs.fr/wiki/2015:groups:tools:ttjets

More details can be found here:

https://twiki.cern.ch/twiki/bin/viewauth/CMS/TOPRivetForLesHouches

and the result there:

http://ebouvier.web.cern.ch/ebouvier/TOPRivetForLesHouches/plots/CMS\_AN\_PseudoTop/index.html

### How-to-setup (and run)

    cp $CMSSW_BASE/src/TopMonteCarlo/Configuration/analysis/Hadronizer_pythia8_cff.py $CMSSW_BASE/src/TopMonteCarlo/Configuration/python/Hadronizer_pythia8_cff.py
    cd $CMSSW_BASE/src/TopMonteCarlo/Configuration/python
    cmsRun runRivietWithPythia8.py
    cmsRun runRivietWithHerwig.py
    rivet-mkhtml -c $CMSSW_BASE/src/TopMonteCarlo/RivetTop/data/CMS_LesHouches2015.plot Pythia8.yoda:'Powheg+Pythia 8' Herwig.yoda:'Powheg+Herwig++'
    firefox plots/index.html &
