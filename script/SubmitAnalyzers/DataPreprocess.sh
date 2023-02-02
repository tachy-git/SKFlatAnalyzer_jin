SKFlat.py -a DataPreprocess -l SampleLists/signalSamples.txt -e 2016preVFP -n 10 --userflags Skim1E2Mu &
SKFlat.py -a DataPreprocess -l SampleLists/signalSamples.txt -e 2016postVFP -n 10 --userflags Skim1E2Mu &
SKFlat.py -a DataPreprocess -l SampleLists/signalSamples.txt -e 2017 -n 10 --userflags Skim1E2Mu &
SKFlat.py -a DataPreprocess -l SampleLists/signalSamples.txt -e 2018 -n 10 --userflags Skim1E2Mu &
SKFlat.py -a DataPreprocess -l SampleLists/bkgSamples.txt -e 2016preVFP --skim SkimTree_SS2lOR3l -n 10 --userflags Skim1E2Mu &
SKFlat.py -a DataPreprocess -l SampleLists/bkgSamples.txt -e 2016postVFP --skim SkimTree_SS2lOR3l -n 10 --userflags Skim1E2Mu &
SKFlat.py -a DataPreprocess -l SampleLists/bkgSamples.txt -e 2017 --skim SkimTree_SS2lOR3l -n 10 --userflags Skim1E2Mu &
SKFlat.py -a DataPreprocess -l SampleLists/bkgSamples.txt -e 2018 --skim SkimTree_SS2lOR3l -n 10 --userflags Skim1E2Mu &

SKFlat.py -a DataPreprocess -l SampleLists/signalSamples.txt -e 2016preVFP -n 10 --userflags Skim3Mu &
SKFlat.py -a DataPreprocess -l SampleLists/signalSamples.txt -e 2016postVFP -n 10 --userflags Skim3Mu &
SKFlat.py -a DataPreprocess -l SampleLists/signalSamples.txt -e 2017 -n 10 --userflags Skim3Mu &
SKFlat.py -a DataPreprocess -l SampleLists/signalSamples.txt -e 2018 -n 10 --userflags Skim3Mu &
SKFlat.py -a DataPreprocess -l SampleLists/bkgSamples.txt -e 2016preVFP --skim SkimTree_SS2lOR3l -n 10 --userflags Skim3Mu &
SKFlat.py -a DataPreprocess -l SampleLists/bkgSamples.txt -e 2016postVFP --skim SkimTree_SS2lOR3l -n 10 --userflags Skim3Mu &
SKFlat.py -a DataPreprocess -l SampleLists/bkgSamples.txt -e 2017 --skim SkimTree_SS2lOR3l -n 10 --userflags Skim3Mu &
SKFlat.py -a DataPreprocess -l SampleLists/bkgSamples.txt -e 2018 --skim SkimTree_SS2lOR3l -n 10 --userflags Skim3Mu &
