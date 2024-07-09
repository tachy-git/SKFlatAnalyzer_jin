#!/bin/bash
ERA=$1
CHANNEL=$2
MEMORY=$3

if [[ $CHANNEL == "Skim1E2Mu" ]]; then
    DATASTREAM="MuonEG"
elif [[ $CHANNEL == "Skim3Mu" ]]; then
    DATASTREAM="DoubleMuon"
else
    echo "Wrong channel $CHANNEL"
    exit 1
fi

SKFlat.py -a PromptEstimator -i $DATASTREAM --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL --memory $MEMORY --python &
SKFlat.py -a MatrixEstimator -i $DATASTREAM --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL --memory $MEMORY --python &
#SKFlat.py -a MatrixEstimator -i DYJets_MG --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL --memory $MEMORY --python &
#SKFlat.py -a MatrixEstimator -i DYJets10to50_MG --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL --memory $MEMORY --python &
#SKFlat.py -a MatrixEstimator -i TTG --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL --memory $MEMORY --python &
#SKFlat.py -a MatrixEstimator -i WWG --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL --memory $MEMORY --python &
SKFlat.py -a PromptEstimator -i TTLL_powheg --skim SkimTree_SS2lOR3l -n 30 -e ${ERA} --userflags $CHANNEL,RunSyst --memory $MEMORY --python & 
SKFlat.py -a PromptEstimator -i DYJets_MG --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL,RunSyst --memory $MEMORY --python &
SKFlat.py -a PromptEstimator -i DYJets10to50_MG --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL,RunSyst --memory $MEMORY --python & 
SKFlat.py -a PromptEstimator -i WZTo3LNu_amcatnlo --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL,RunSyst --memory $MEMORY --python &
SKFlat.py -a PromptEstimator -i ZZTo4L_powheg --skim SkimTree_SS2lOR3l -n 30 -e ${ERA} --userflags $CHANNEL,RunSyst --memory $MEMORY --python &
SKFlat.py -a PromptEstimator -i ZGToLLG --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL,RunSyst --memory $MEMORY --python &
SKFlat.py -a PromptEstimator -i ttWToLNu --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL,RunSyst --memory $MEMORY --python &
SKFlat.py -a PromptEstimator -i ttZToLLNuNu --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL,RunSyst --memory $MEMORY --python &
SKFlat.py -a PromptEstimator -i ttHToNonbb --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL,RunSyst --memory $MEMORY --python &
SKFlat.py -a PromptEstimator -i tZq --skim SkimTree_SS2lOR3l -n 20 -e ${ERA} --userflags $CHANNEL,RunSyst --memory $MEMORY --python &
SKFlat.py -a PromptEstimator -i tHq --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL,RunSyst --memory $MEMORY --python &
SKFlat.py -a PromptEstimator -i WWW --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL,RunSyst --memory $MEMORY --python &
SKFlat.py -a PromptEstimator -i WWZ --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL,RunSyst --memory $MEMORY --python &
SKFlat.py -a PromptEstimator -i WZZ --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL,RunSyst --memory $MEMORY --python & 
SKFlat.py -a PromptEstimator -i ZZZ --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL,RunSyst --memory $MEMORY --python &
SKFlat.py -a PromptEstimator -i WWG --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL,RunSyst --memory $MEMORY --python &
SKFlat.py -a PromptEstimator -i TTG --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL,RunSyst --memory $MEMORY --python &
SKFlat.py -a PromptEstimator -i TTTT --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL,RunSyst --memory $MEMORY --python &
SKFlat.py -a PromptEstimator -i VBF_HToZZTo4L --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL,RunSyst --memory $MEMORY --python &
SKFlat.py -a PromptEstimator -i GluGluHToZZTo4L --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL,RunSyst --memory $MEMORY --python &
#SKFlat.py -a PromptEstimator -i TTToHcToWAToMuMu_MHc-70_MA-15 -n 10 -e ${ERA} --userflags $CHANNEL,RunSyst --memory $MEMORY --python & 
#SKFlat.py -a PromptEstimator -i TTToHcToWAToMuMu_MHc-70_MA-40 -n 10 -e ${ERA} --userflags $CHANNEL,RunSyst --memory $MEMORY --python &
#SKFlat.py -a PromptEstimator -i TTToHcToWAToMuMu_MHc-70_MA-65 -n 10 -e ${ERA} --userflags $CHANNEL,RunSyst --memory $MEMORY --python &
#SKFlat.py -a PromptEstimator -i TTToHcToWAToMuMu_MHc-100_MA-15 -n 10 -e ${ERA} --userflags $CHANNEL,RunSyst --memory $MEMORY --python &
#SKFlat.py -a PromptEstimator -i TTToHcToWAToMuMu_MHc-100_MA-60 -n 10 -e ${ERA} --userflags $CHANNEL,RunSyst --memory $MEMORY --python &
#SKFlat.py -a PromptEstimator -i TTToHcToWAToMuMu_MHc-100_MA-95 -n 10 -e ${ERA} --userflags $CHANNEL,RunSyst --memory $MEMORY --python &
#SKFlat.py -a PromptEstimator -i TTToHcToWAToMuMu_MHc-130_MA-15 -n 10 -e ${ERA} --userflags $CHANNEL,RunSyst --memory $MEMORY --python &
#SKFlat.py -a PromptEstimator -i TTToHcToWAToMuMu_MHc-130_MA-55 -n 10 -e ${ERA} --userflags $CHANNEL,RunSyst --memory $MEMORY --python &
#SKFlat.py -a PromptEstimator -i TTToHcToWAToMuMu_MHc-130_MA-90 -n 10 -e ${ERA} --userflags $CHANNEL,RunSyst --memory $MEMORY --python &
#SKFlat.py -a PromptEstimator -i TTToHcToWAToMuMu_MHc-130_MA-125 -n 10 -e ${ERA} --userflags $CHANNEL,RunSyst --memory $MEMORY --python &
#SKFlat.py -a PromptEstimator -i TTToHcToWAToMuMu_MHc-160_MA-15 -n 10 -e ${ERA} --userflags $CHANNEL,RunSyst --memory $MEMORY --python & 
#SKFlat.py -a PromptEstimator -i TTToHcToWAToMuMu_MHc-160_MA-85 -n 10 -e ${ERA} --userflags $CHANNEL,RunSyst --memory $MEMORY --python &
#SKFlat.py -a PromptEstimator -i TTToHcToWAToMuMu_MHc-160_MA-120 -n 10 -e ${ERA} --userflags $CHANNEL,RunSyst --memory $MEMORY --python &
#SKFlat.py -a PromptEstimator -i TTToHcToWAToMuMu_MHc-160_MA-155 -n 10 -e ${ERA} --userflags $CHANNEL,RunSyst --memory $MEMORY --python &
