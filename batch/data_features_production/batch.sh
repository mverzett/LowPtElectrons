#! /bin/bash

sbox=$PWD
cd /afs/cern.ch/work/m/mverzett/RK94New/src/CMGTools/LowPtElectrons
eval `scramv1 runtime -sh`
cd $sbox

cmsRun /afs/cern.ch/work/m/mverzett/RK94New/src/CMGTools/LowPtElectrons/run/raw2features.py $@
cp *.root /eos/cms/store/cmst3/user/mverzett/low_pt_electrons/.
rm *.root