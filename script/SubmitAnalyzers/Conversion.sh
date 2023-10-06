#!/bin/bash
ERA=$1
CHANNEL=$2

if [[ $CHANNEL == "Skim1E2Mu" ]]; then
    DATASTREAM="MuonEG"
elif [[ $CHANNEL == "Skim3Mu" ]]; then
    DATASTREAM="DoubleMuon"
else
    echo "Wrong channel $CHANNEL"
    exit 1
fi
SKFlat.py -a MeasConversion -i $DATASTREAM --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL  --python &
SKFlat.py -a MeasConvMatrix -i $DATASTREAM --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL  --python &
SKFlat.py -a MeasConversion -i DYJets_MG  --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL  --python &
SKFlat.py -a MeasConversion -i DYJets10to50_MG   --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL  --python & 
SKFlat.py -a MeasConversion -i WZTo3LNu_amcatnlo --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL  --python &
SKFlat.py -a MeasConversion -i ZZTo4L_powheg --skim SkimTree_SS2lOR3l -n 30 -e ${ERA} --userflags $CHANNEL  --python &
SKFlat.py -a MeasConversion -i ZGToLLG --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL  --python &
SKFlat.py -a MeasConversion -i ttWToLNu --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL  --python &
SKFlat.py -a MeasConversion -i ttZToLLNuNu --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL  --python &
SKFlat.py -a MeasConversion -i ttHToNonbb --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL  --python &
SKFlat.py -a MeasConversion -i tZq --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL  --python &
SKFlat.py -a MeasConversion -i tHq --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL  --python &
SKFlat.py -a MeasConversion -i WWW --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL  --python &
SKFlat.py -a MeasConversion -i WWZ --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL  --python &
SKFlat.py -a MeasConversion -i WZZ --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL  --python & 
SKFlat.py -a MeasConversion -i ZZZ --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL  --python &
SKFlat.py -a MeasConversion -i WWG --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL  --python &
SKFlat.py -a MeasConversion -i TTG --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL  --python &
SKFlat.py -a MeasConversion -i TTTT --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL  --python &
SKFlat.py -a MeasConversion -i VBF_HToZZTo4L --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL  --python &
SKFlat.py -a MeasConversion -i GluGluHToZZTo4L --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL  --python &
