#!/bin/bash
ERA=$1
CHANNEL=$2
FLAG=$3
DATASTREAM=""
QCD=""

if [[ $ERA == "2016"* ]]; then
    if [[ $CHANNEL == "MeasFakeMu"* ]]; then
        DATASTREAM="DoubleMuon"
        QCD="qcdMuEnriched"
    elif [[ $CHANNEL == "MeasFakeEl"* ]]; then
        DATASTREAM="DoubleEG"
        QCD="qcdEMEnriched"
    else
        echo "Wrong channel $CHANNEL"
        exit 1
    fi
    SKFlat.py -a MeasFakeRateV4 -i $DATASTREAM               -n 10 -e $ERA --userflags $CHANNEL,$FLAG &
    SKFlat.py -a MeasFakeRateV4 -l SampleLists/fakeHeavy.txt -n 20 -e $ERA --userflags $CHANNEL,$FLAG &
    SKFlat.py -a MeasFakeRateV4 -l SampleLists/fakeLight.txt -n 10 -e $ERA --userflags $CHANNEL,$FLAG &
    SKFlat.py -a MeasFakeRateV4 -l SampleLists/${QCD}.txt    -n 10 -e $ERA --userflags $CHANNEL,$FLAG &
elif [[ $ERA == "2017" ]]; then
    if [[ $CHANNEL == "MeasFakeMu"* ]]; then
        DATASTREAM="DoubleMuon"
        QCD="qcdMuEnriched"
    elif [[ $CHANNEL == "MeasFakeEl"* ]]; then
        DATASTREAM="SingleElectron"
        QCD="qcdEMEnriched"
    else
        echo "Wrong channel $CHANNEL"
        exit 1
    fi
    SKFlat.py -a MeasFakeRateV4 -i $DATASTREAM               -n 10 -e $ERA --userflags $CHANNEL,$FLAG &
    SKFlat.py -a MeasFakeRateV4 -l SampleLists/fakeHeavy.txt -n 20 -e $ERA --userflags $CHANNEL,$FLAG &
    SKFlat.py -a MeasFakeRateV4 -l SampleLists/fakeLight.txt -n 15 -e $ERA --userflags $CHANNEL,$FLAG &
    SKFlat.py -a MeasFakeRateV4 -l SampleLists/${QCD}.txt    -n 15 -e $ERA --userflags $CHANNEL,$FLAG &
elif [[ $ERA == "2018" ]]; then
    if [[ $CHANNEL == "MeasFakeMu"* ]]; then
        DATASTREAM="DoubleMuon"
        QCD="qcdMuEnriched"
    elif [[ $CHANNEL == "MeasFakeEl"* ]]; then
        DATASTREAM="EGamma"
        QCD="qcdEMEnriched"
    else
        echo "Wrong channel $CHANNEL"
        exit 1
    fi
    SKFlat.py -a MeasFakeRateV4 -i $DATASTREAM               -n 15 -e $ERA --userflags $CHANNEL,$FLAG &
    SKFlat.py -a MeasFakeRateV4 -l SampleLists/fakeHeavy.txt -n 30 -e $ERA --userflags $CHANNEL,$FLAG &
    SKFlat.py -a MeasFakeRateV4 -l SampleLists/fakeLight.txt -n 20 -e $ERA --userflags $CHANNEL,$FLAG &
    SKFlat.py -a MeasFakeRateV4 -l SampleLists/${QCD}.txt    -n 20 -e $ERA --userflags $CHANNEL,$FLAG &
else
    echo "Wrong era $ERA"
    exit 1
fi
