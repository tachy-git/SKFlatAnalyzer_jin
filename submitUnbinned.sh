#!/bin/bash
ERA=$1
CHANNEL=$2
NETWORK=$3
MEMORY=$4

SKFlat.py -a PromptUnbinned -i DoubleMuon --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL,$NETWORK --memory $MEMORY --python &
SKFlat.py -a MatrixUnbinned -i DoubleMuon --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL,$NETWORK --memory $MEMORY --python &
SKFlat.py -a PromptUnbinned -i TTLL_powheg --skim SkimTree_SS2lOR3l -n 20 -e ${ERA} --userflags $CHANNEL,$NETWORK --memory $MEMORY --python & 
SKFlat.py -a PromptUnbinned -i DYJets --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL,$NETWORK --memory $MEMORY --python &
SKFlat.py -a PromptUnbinned -i DYJets10to50_MG --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL,$NETWORK --memory $MEMORY --python & 
SKFlat.py -a PromptUnbinned -i WZTo3LNu_amcatnlo --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL,$NETWORK --memory $MEMORY --python &
SKFlat.py -a PromptUnbinned -i ZZTo4L_powheg --skim SkimTree_SS2lOR3l -n 30 -e ${ERA} --userflags $CHANNEL,$NETWORK --memory $MEMORY --python &
SKFlat.py -a PromptUnbinned -i ZGToLLG --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL,$NETWORK --memory $MEMORY --python &
SKFlat.py -a PromptUnbinned -i ttWToLNu --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL,$NETWORK --memory $MEMORY --python &
SKFlat.py -a PromptUnbinned -i ttZToLLNuNu --skim SkimTree_SS2lOR3l -n 20 -e ${ERA} --userflags $CHANNEL,$NETWORK --memory $MEMORY --python &
SKFlat.py -a PromptUnbinned -i ttHToNonbb --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL,$NETWORK --memory $MEMORY --python &
SKFlat.py -a PromptUnbinned -i tZq --skim SkimTree_SS2lOR3l -n 20 -e ${ERA} --userflags $CHANNEL,$NETWORK --memory $MEMORY --python &
SKFlat.py -a PromptUnbinned -i tHq --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL,$NETWORK --memory $MEMORY --python &
SKFlat.py -a PromptUnbinned -i WWW --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL,$NETWORK --memory $MEMORY --python &
SKFlat.py -a PromptUnbinned -i WWZ --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL,$NETWORK --memory $MEMORY --python &
SKFlat.py -a PromptUnbinned -i WZZ --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL,$NETWORK --memory $MEMORY --python & 
SKFlat.py -a PromptUnbinned -i ZZZ --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL,$NETWORK --memory $MEMORY --python &
SKFlat.py -a PromptUnbinned -i WWG --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL,$NETWORK --memory $MEMORY --python &
SKFlat.py -a PromptUnbinned -i TTG --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL,$NETWORK --memory $MEMORY --python &
SKFlat.py -a PromptUnbinned -i TTTT --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL,$NETWORK --memory $MEMORY --python &
SKFlat.py -a PromptUnbinned -i VBF_HToZZTo4L --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL,$NETWORK --memory $MEMORY --python &
SKFlat.py -a PromptUnbinned -i GluGluHToZZTo4L --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL,$NETWORK --memory $MEMORY --python &
SKFlat.py -a PromptUnbinned -i TTToHcToWAToMuMu_MHc-70_MA-15 -n 10 -e ${ERA} --userflags $CHANNEL,$NETWORK --memory $MEMORY --python & 
SKFlat.py -a PromptUnbinned -i TTToHcToWAToMuMu_MHc-70_MA-40 -n 10 -e ${ERA} --userflags $CHANNEL,$NETWORK --memory $MEMORY --python & 
SKFlat.py -a PromptUnbinned -i TTToHcToWAToMuMu_MHc-70_MA-65 -n 10 -e ${ERA} --userflags $CHANNEL,$NETWORK --memory $MEMORY --python &
SKFlat.py -a PromptUnbinned -i TTToHcToWAToMuMu_MHc-100_MA-15 -n 10 -e ${ERA} --userflags $CHANNEL,$NETWORK --memory $MEMORY --python &
SKFlat.py -a PromptUnbinned -i TTToHcToWAToMuMu_MHc-100_MA-60 -n 10 -e ${ERA} --userflags $CHANNEL,$NETWORK --memory $MEMORY --python &
SKFlat.py -a PromptUnbinned -i TTToHcToWAToMuMu_MHc-100_MA-95 -n 10 -e ${ERA} --userflags $CHANNEL,$NETWORK --memory $MEMORY --python &
SKFlat.py -a PromptUnbinned -i TTToHcToWAToMuMu_MHc-130_MA-15 -n 10 -e ${ERA} --userflags $CHANNEL,$NETWORK --memory $MEMORY --python &
SKFlat.py -a PromptUnbinned -i TTToHcToWAToMuMu_MHc-130_MA-55 -n 10 -e ${ERA} --userflags $CHANNEL,$NETWORK --memory $MEMORY --python &
SKFlat.py -a PromptUnbinned -i TTToHcToWAToMuMu_MHc-130_MA-90 -n 10 -e ${ERA} --userflags $CHANNEL,$NETWORK --memory $MEMORY --python &
SKFlat.py -a PromptUnbinned -i TTToHcToWAToMuMu_MHc-130_MA-125 -n 10 -e ${ERA} --userflags $CHANNEL,$NETWORK --memory $MEMORY --python &
SKFlat.py -a PromptUnbinned -i TTToHcToWAToMuMu_MHc-160_MA-15 -n 10 -e ${ERA} --userflags $CHANNEL,$NETWORK --memory $MEMORY --python & 
SKFlat.py -a PromptUnbinned -i TTToHcToWAToMuMu_MHc-160_MA-85 -n 10 -e ${ERA} --userflags $CHANNEL,$NETWORK --memory $MEMORY --python &
SKFlat.py -a PromptUnbinned -i TTToHcToWAToMuMu_MHc-160_MA-120 -n 10 -e ${ERA} --userflags $CHANNEL,$NETWORK --memory $MEMORY --python &
SKFlat.py -a PromptUnbinned -i TTToHcToWAToMuMu_MHc-160_MA-155 -n 10 -e ${ERA} --userflags $CHANNEL,$NETWORK --memory $MEMORY --python &
