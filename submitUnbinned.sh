#!/bin/bash
ERA=$1
CHANNEL=$2
NETWORK=$3

SKFlat.py -a SkimPromptTree -i DoubleMuon --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL,$NETWORK --python &
SKFlat.py -a SkimNonpromptTree -i DoubleMuon --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL,$NETWORK --python &
SKFlat.py -a SkimPromptTree -i TTLL_powheg --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL,$NETWORK --python & 
SKFlat.py -a SkimPromptTree -i DYJets --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL,$NETWORK --python &
SKFlat.py -a SkimPromptTree -i DYJets10to50_MG --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL,$NETWORK --python & 
SKFlat.py -a SkimPromptTree -i WZTo3LNu_amcatnlo --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL,$NETWORK --python &
SKFlat.py -a SkimPromptTree -i ZZTo4L_powheg --skim SkimTree_SS2lOR3l -n 30 -e ${ERA} --userflags $CHANNEL,$NETWORK --python &
SKFlat.py -a SkimPromptTree -i ZGToLLG --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL,$NETWORK --python &
SKFlat.py -a SkimPromptTree -i ttWToLNu --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL,$NETWORK --python &
SKFlat.py -a SkimPromptTree -i ttZToLLNuNu --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL,$NETWORK --python &
SKFlat.py -a SkimPromptTree -i ttHToNonbb --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL,$NETWORK --python &
SKFlat.py -a SkimPromptTree -i tZq --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL,$NETWORK --python &
SKFlat.py -a SkimPromptTree -i tHq -n 10 -e ${ERA} --userflags $CHANNEL,$NETWORK --python &
SKFlat.py -a SkimPromptTree -i WWW --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL,$NETWORK --python &
SKFlat.py -a SkimPromptTree -i WWZ --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL,$NETWORK --python &
SKFlat.py -a SkimPromptTree -i WZZ --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL,$NETWORK --python & 
SKFlat.py -a SkimPromptTree -i ZZZ --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL,$NETWORK --python &
SKFlat.py -a SkimPromptTree -i WWG -n 10 -e ${ERA} --userflags $CHANNEL,$NETWORK --python &
SKFlat.py -a SkimPromptTree -i TTG --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL,$NETWORK --python &
SKFlat.py -a SkimPromptTree -i TTTT -n 10 -e ${ERA} --userflags $CHANNEL,$NETWORK --python &
SKFlat.py -a SkimPromptTree -i VBF_HToZZTo4L --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL,$NETWORK --python &
SKFlat.py -a SkimPromptTree -i GluGluHToZZTo4L --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL,$NETWORK --python &
SKFlat.py -a SkimPromptTree -i TTToHcToWAToMuMu_MHc-70_MA-15 -n 10 -e ${ERA} --userflags $CHANNEL,$NETWORK --python & 
SKFlat.py -a SkimPromptTree -i TTToHcToWAToMuMu_MHc-70_MA-40 -n 10 -e ${ERA} --userflags $CHANNEL,$NETWORK --python &
SKFlat.py -a SkimPromptTree -i TTToHcToWAToMuMu_MHc-70_MA-65 -n 10 -e ${ERA} --userflags $CHANNEL,$NETWORK --python &
SKFlat.py -a SkimPromptTree -i TTToHcToWAToMuMu_MHc-100_MA-15 -n 10 -e ${ERA} --userflags $CHANNEL,$NETWORK --python &
SKFlat.py -a SkimPromptTree -i TTToHcToWAToMuMu_MHc-100_MA-60 -n 10 -e ${ERA} --userflags $CHANNEL,$NETWORK --python &
SKFlat.py -a SkimPromptTree -i TTToHcToWAToMuMu_MHc-100_MA-95 -n 10 -e ${ERA} --userflags $CHANNEL,$NETWORK --python &
SKFlat.py -a SkimPromptTree -i TTToHcToWAToMuMu_MHc-130_MA-15 -n 10 -e ${ERA} --userflags $CHANNEL,$NETWORK --python &
SKFlat.py -a SkimPromptTree -i TTToHcToWAToMuMu_MHc-130_MA-55 -n 10 -e ${ERA} --userflags $CHANNEL,$NETWORK --python &
SKFlat.py -a SkimPromptTree -i TTToHcToWAToMuMu_MHc-130_MA-90 -n 10 -e ${ERA} --userflags $CHANNEL,$NETWORK --python &
SKFlat.py -a SkimPromptTree -i TTToHcToWAToMuMu_MHc-130_MA-125 -n 10 -e ${ERA} --userflags $CHANNEL,$NETWORK --python &
SKFlat.py -a SkimPromptTree -i TTToHcToWAToMuMu_MHc-160_MA-15 -n 10 -e ${ERA} --userflags $CHANNEL,$NETWORK --python & 
SKFlat.py -a SkimPromptTree -i TTToHcToWAToMuMu_MHc-160_MA-85 -n 10 -e ${ERA} --userflags $CHANNEL,$NETWORK --python &
SKFlat.py -a SkimPromptTree -i TTToHcToWAToMuMu_MHc-160_MA-120 -n 10 -e ${ERA} --userflags $CHANNEL,$NETWORK --python &
SKFlat.py -a SkimPromptTree -i TTToHcToWAToMuMu_MHc-160_MA-155 -n 10 -e ${ERA} --userflags $CHANNEL,$NETWORK --python &
