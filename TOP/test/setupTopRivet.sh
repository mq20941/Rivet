#!/bin/bash
source $CMSSW_BASE/src/GeneratorInterface/RivetInterface/test/rivetSetup.sh

export RIVET_REF_PATH=$RIVET_REF_PATH:$CMSSW_BASE/src/Rivet/Top/data
export RIVET_INFO_PATH=$RIVET_INFO_PATH:$CMSSW_BASE/src/Rivet/Top/data
export RIVET_PLOT_PATH=$RIVET_PLOT_PATH:$CMSSW_BASE/src/Rivet/Top/data

which yodamerge &> /dev/null || GETYODA=1
if [ $GETYODA==1 ]; then
  eval `scram tool info yoda | grep YODA_BASE`
  export PATH=$PATH:$YODA_BASE/bin
fi
