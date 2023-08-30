#!/bin/bash
ERA=$1
CHANNEL=$2
DATASTREAM=""

if [[ $ERA == "2016"* ]]; then
    if [[ $CHANNEL == "MeasFakeMu"* ]]; then
        DATASTREAM="DoubleMuon"
    elif [[ $CHANNEL == "MeasFakeEl"* ]]; then
        DATASTREAM="DoubleEG"
    else
        echo "Wrong channel $CHANNEL"
        exit 1
    fi
    SKFlat.py -a MeasFakeRateV3 -i $DATASTREAM                     -n 10 -e $ERA --userflags $CHANNEL --python &
    SKFlat.py -a MeasFakeRateV3 -i TTLJ_powheg                     -n 30 -e $ERA --userflags $CHANNEL,RunSyst --python &
    SKFlat.py -a MeasFakeRateV3 -i TTLL_powheg                     -n 20 -e $ERA --userflags $CHANNEL,RunSyst --python &
    SKFlat.py -a MeasFakeRateV3 -i WJets_MG                        -n 20 -e $ERA --userflags $CHANNEL,RunSyst --python &
    SKFlat.py -a MeasFakeRateV3 -i DYJets                          -n 20 -e $ERA --userflags $CHANNEL,RunSyst --python &
    SKFlat.py -a MeasFakeRateV3 -i DYJets10to50_MG                 -n 10 -e $ERA --userflags $CHANNEL,RunSyst --python &
    SKFlat.py -a MeasFakeRateV3 -i WW_pythia                       -n 10 -e $ERA --userflags $CHANNEL,RunSyst --python &
    SKFlat.py -a MeasFakeRateV3 -i WZ_pythia                       -n 10 -e $ERA --userflags $CHANNEL,RunSyst --python &
    SKFlat.py -a MeasFakeRateV3 -i ZZ_pythia                       -n 10 -e $ERA --userflags $CHANNEL,RunSyst --python &
    SKFlat.py -a MeasFakeRateV3 -i SingleTop_sch_Lep               -n 10 -e $ERA --userflags $CHANNEL,RunSyst --python &
    SKFlat.py -a MeasFakeRateV3 -i SingleTop_tch_top_Incl          -n 10 -e $ERA --userflags $CHANNEL,RunSyst --python &
    SKFlat.py -a MeasFakeRateV3 -i SingleTop_tch_antitop_Incl      -n 10 -e $ERA --userflags $CHANNEL,RunSyst --python &
    SKFlat.py -a MeasFakeRateV3 -i SingleTop_tW_top_NoFullyHad     -n 10 -e $ERA --userflags $CHANNEL,RunSyst --python &
    SKFlat.py -a MeasFakeRateV3 -i SingleTop_tW_antitop_NoFullyHad -n 10 -e $ERA --userflags $CHANNEL,RunSyst --python &
elif [[ $ERA == "2017" ]]; then
    if [[ $CHANNEL == "MeasFakeMu"* ]]; then
        DATASTREAM="DoubleMuon"
    elif [[ $CHANNEL == "MeasFakeEl"* ]]; then
        DATASTREAM="SingleElectron"
    else
        echo "Wrong channel $CHANNEL"
        exit 1
    fi
    SKFlat.py -a MeasFakeRateV3 -i $DATASTREAM                     -n 10 -e $ERA --userflags $CHANNEL --python &
    SKFlat.py -a MeasFakeRateV3 -i TTLJ_powheg                     -n 30 -e $ERA --userflags $CHANNEL,RunSyst --python &
    SKFlat.py -a MeasFakeRateV3 -i TTLL_powheg                     -n 20 -e $ERA --userflags $CHANNEL,RunSyst --python &
    SKFlat.py -a MeasFakeRateV3 -i WJets_MG                        -n 20 -e $ERA --userflags $CHANNEL,RunSyst --python &
    SKFlat.py -a MeasFakeRateV3 -i DYJets                          -n 20 -e $ERA --userflags $CHANNEL,RunSyst --python &
    SKFlat.py -a MeasFakeRateV3 -i DYJets10to50_MG                 -n 10 -e $ERA --userflags $CHANNEL,RunSyst --python &
    SKFlat.py -a MeasFakeRateV3 -i WW_pythia                       -n 15 -e $ERA --userflags $CHANNEL,RunSyst --python &
    SKFlat.py -a MeasFakeRateV3 -i WZ_pythia                       -n 15 -e $ERA --userflags $CHANNEL,RunSyst --python &
    SKFlat.py -a MeasFakeRateV3 -i ZZ_pythia                       -n 15 -e $ERA --userflags $CHANNEL,RunSyst --python &
    SKFlat.py -a MeasFakeRateV3 -i SingleTop_sch_Lep               -n 15 -e $ERA --userflags $CHANNEL,RunSyst --python &
    SKFlat.py -a MeasFakeRateV3 -i SingleTop_tch_top_Incl          -n 15 -e $ERA --userflags $CHANNEL,RunSyst --python &
    SKFlat.py -a MeasFakeRateV3 -i SingleTop_tch_antitop_Incl      -n 15 -e $ERA --userflags $CHANNEL,RunSyst --python &
    SKFlat.py -a MeasFakeRateV3 -i SingleTop_tW_top_NoFullyHad     -n 15 -e $ERA --userflags $CHANNEL,RunSyst --python &
    SKFlat.py -a MeasFakeRateV3 -i SingleTop_tW_antitop_NoFullyHad -n 15 -e $ERA --userflags $CHANNEL,RunSyst --python &
elif [[ $ERA == "2018" ]]; then
    if [[ $CHANNEL == "MeasFakeMu"* ]]; then
        DATASTREAM="DoubleMuon"
    elif [[ $CHANNEL == "MeasFakeEl"* ]]; then
        DATASTREAM="EGamma"
    else
        echo "Wrong channel $CHANNEL"
        exit 1
    fi
    SKFlat.py -a MeasFakeRateV3 -i $DATASTREAM                     -n 10 -e $ERA --userflags $CHANNEL --python &
    SKFlat.py -a MeasFakeRateV3 -i TTLJ_powheg                     -n 40 -e $ERA --userflags $CHANNEL,RunSyst --python &
    SKFlat.py -a MeasFakeRateV3 -i TTLL_powheg                     -n 30 -e $ERA --userflags $CHANNEL,RunSyst --python &
    SKFlat.py -a MeasFakeRateV3 -i WJets_MG                        -n 30 -e $ERA --userflags $CHANNEL,RunSyst --python &
    SKFlat.py -a MeasFakeRateV3 -i DYJets                          -n 30 -e $ERA --userflags $CHANNEL,RunSyst --python &
    SKFlat.py -a MeasFakeRateV3 -i DYJets10to50_MG                 -n 10 -e $ERA --userflags $CHANNEL,RunSyst --python &
    SKFlat.py -a MeasFakeRateV3 -i WW_pythia                       -n 20 -e $ERA --userflags $CHANNEL,RunSyst --python &
    SKFlat.py -a MeasFakeRateV3 -i WZ_pythia                       -n 20 -e $ERA --userflags $CHANNEL,RunSyst --python &
    SKFlat.py -a MeasFakeRateV3 -i ZZ_pythia                       -n 20 -e $ERA --userflags $CHANNEL,RunSyst --python &
    SKFlat.py -a MeasFakeRateV3 -i SingleTop_sch_Lep               -n 20 -e $ERA --userflags $CHANNEL,RunSyst --python &
    SKFlat.py -a MeasFakeRateV3 -i SingleTop_tch_top_Incl          -n 20 -e $ERA --userflags $CHANNEL,RunSyst --python &
    SKFlat.py -a MeasFakeRateV3 -i SingleTop_tch_antitop_Incl      -n 20 -e $ERA --userflags $CHANNEL,RunSyst --python &
    SKFlat.py -a MeasFakeRateV3 -i SingleTop_tW_top_NoFullyHad     -n 20 -e $ERA --userflags $CHANNEL,RunSyst --python &
    SKFlat.py -a MeasFakeRateV3 -i SingleTop_tW_antitop_NoFullyHad -n 20 -e $ERA --userflags $CHANNEL,RunSyst --python &
else
    echo "Wrong era $ERA"
    exit 1
fi
