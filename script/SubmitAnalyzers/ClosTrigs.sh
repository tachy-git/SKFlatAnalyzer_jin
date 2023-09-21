#!/bin/sh
ERA=$1

SKFlat.py -a ClosDiLepTrigs -i DYJets -n 10 -e $ERA --userflags RunDiMu --python &
SKFlat.py -a ClosDiLepTrigs -i TTLL_powheg -n 10 -e $ERA --userflags RunDiMu --python &
SKFlat.py -a ClosDiLepTrigs -i TTLL_powheg -n 10 -e $ERA --userflags RunEMu --python &
SKFlat.py -a ClosTriLepTrigs -i DYJets --skim SkimTree_SS2lOR3l -n 10 -e $ERA --userflags Skim1E2Mu --python &
SKFlat.py -a ClosTriLepTrigs -i DYJets --skim SkimTree_SS2lOR3l -n 10 -e $ERA --userflags Skim3Mu --python &
SKFlat.py -a ClosTriLepTrigs -i TTLL_powheg --skim SkimTree_SS2lOR3l -n 10 -e $ERA --userflags Skim1E2Mu --python &
SKFlat.py -a ClosTriLepTrigs -i TTLL_powheg --skim SkimTree_SS2lOR3l -n 10 -e $ERA --userflags Skim3Mu --python &
SKFlat.py -a ClosTriLepTrigs -l SampleLists/signalSamples.txt -n 10 -e $ERA --userflags Skim1E2Mu --python &
SKFlat.py -a ClosTriLepTrigs -l SampleLists/signalSamples.txt -n 10 -e $ERA --userflags Skim3Mu --python &

