#!/bin/bash
ERA=$1
CHANNEL=$2

SKFlat.py -a MatrixEstimator -i MuonEG --skim SkimTree_SS2lOR3l -n 10 -e $ERA --userflags $CHANNEL --memory 4000 --python &
SKFlat.py -a MatrixEstimator -i DYJets_MG --skim SkimTree_SS2lOR3l -n 10 -e $ERA --userflags $CHANNEL --memory 4000 --python &
SKFlat.py -a MatrixEstimator -i DYJets10to50_MG --skim SkimTree_SS2lOR3l -n 10 -e $ERA --userflags $CHANNEL --memory 4000 --python &
SKFlat.py -a MatrixEstimator -i TTG --skim SkimTree_SS2lOR3l -n 10 -e $ERA --userflags $CHANNEL --memory 4000 --python &
SKFlat.py -a MatrixEstimator -i WWG --skim SkimTree_SS2lOR3l -n 10 -e $ERA --userflags $CHANNEL --memory 4000 --python &
