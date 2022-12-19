#!/bin/bash
ANALYZER=$1
ERA=$2
SKFlat.py -a $ANALYZER -i DoubleMuon --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --python &
SKFlat.py -a $ANALYZER -i TTLL_powheg --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --python &
SKFlat.py -a $ANALYZER -i DYJets --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --python &
SKFlat.py -a $ANALYZER -i DYJets10to50_MG --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --python &
SKFlat.py -a $ANALYZER -i WZTo3LNu_amcatnlo --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --python &
SKFlat.py -a $ANALYZER -i ZZTo4L_powheg --skim SkimTree_SS2lOR3l -n 30 -e ${ERA} --python &
SKFlat.py -a $ANALYZER -i ZGToLLG --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --python &
SKFlat.py -a $ANALYZER -i ttWToLNu --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --python &
SKFlat.py -a $ANALYZER -i ttZToLLNuNu --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --python &
SKFlat.py -a $ANALYZER -i ttHToNonbb --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --python &
SKFlat.py -a $ANALYZER -i tZq --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --python &
SKFlat.py -a $ANALYZER -i tHq --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --python &
SKFlat.py -a $ANALYZER -i WWW --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --python &
SKFlat.py -a $ANALYZER -i WWZ --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --python &
SKFlat.py -a $ANALYZER -i WZZ --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --python &
SKFlat.py -a $ANALYZER -i ZZZ --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --python &
SKFlat.py -a $ANALYZER -i WWG --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --python &
SKFlat.py -a $ANALYZER -i TTG --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --python &
SKFlat.py -a $ANALYZER -i TTTT --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --python &
SKFlat.py -a $ANALYZER -i VBF_HToZZTo4L --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --python &
SKFlat.py -a $ANALYZER -i GluGluHToZZTo4L --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --python &
SKFlat.py -a $ANALYZER -i TTToHcToWAToMuMu_MHc-70_MA-15 -n 10 -e ${ERA} --python &
SKFlat.py -a $ANALYZER -i TTToHcToWAToMuMu_MHc-70_MA-40 -n 10 -e ${ERA} --python &
SKFlat.py -a $ANALYZER -i TTToHcToWAToMuMu_MHc-70_MA-65 -n 10 -e ${ERA} --python &
SKFlat.py -a $ANALYZER -i TTToHcToWAToMuMu_MHc-100_MA-15 -n 10 -e ${ERA} --python &
SKFlat.py -a $ANALYZER -i TTToHcToWAToMuMu_MHc-100_MA-60 -n 10 -e ${ERA} --python &
SKFlat.py -a $ANALYZER -i TTToHcToWAToMuMu_MHc-100_MA-95 -n 10 -e ${ERA} --python &
SKFlat.py -a $ANALYZER -i TTToHcToWAToMuMu_MHc-130_MA-15 -n 10 -e ${ERA} --python &
SKFlat.py -a $ANALYZER -i TTToHcToWAToMuMu_MHc-130_MA-55 -n 10 -e ${ERA} --python &
SKFlat.py -a $ANALYZER -i TTToHcToWAToMuMu_MHc-130_MA-90 -n 10 -e ${ERA} --python &
SKFlat.py -a $ANALYZER -i TTToHcToWAToMuMu_MHc-130_MA-125 -n 10 -e ${ERA} --python &
SKFlat.py -a $ANALYZER -i TTToHcToWAToMuMu_MHc-160_MA-15 -n 10 -e ${ERA} --python & 
SKFlat.py -a $ANALYZER -i TTToHcToWAToMuMu_MHc-160_MA-85 -n 10 -e ${ERA} --python &
SKFlat.py -a $ANALYZER -i TTToHcToWAToMuMu_MHc-160_MA-120 -n 10 -e ${ERA} --python &
SKFlat.py -a $ANALYZER -i TTToHcToWAToMuMu_MHc-160_MA-155 -n 10 -e ${ERA} --python &
