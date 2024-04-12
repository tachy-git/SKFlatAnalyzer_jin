#!/bin/sh
echo "Start testing..."

rm -rf $SKFlatOutputDir/$SKFlatV/TutorialRun $SKFlatOutputDir/$SKFlatV/ExampleRun
SKFlat.py -a ExampleRun -i DYJets -e 2016a -n 10 --reduction 100
SKFlat.py -a ExampleRun -l SampleLists/fakeSamples.txt -e 2016b -n 10 --reduction 100
SKFlat.py -a TutorialRun -i DYJets -e 2017 -n 10 --python --reduction 100
SKFlat.py -a TutorialRun -l SampleLists/fakeSamples.txt -e 2017 -n 10 --python --reduction 100
SKFlat.py -a ExampleRun -i DYJets -e 2018 -n 10 --userflags RunSyst,RunNewPDF --reduction 100
SKFlat.py -a ExampleRun -i SingleMuon -e 2016preVFP -n 10 --reduction 100
SKFlat.py -a ExampleRun -i SingleMuon -p C -e 2017 -n 10 --reduction 100
