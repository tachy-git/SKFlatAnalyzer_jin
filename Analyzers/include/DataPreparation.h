#ifndef DataPreparation_h
#define DataPreparation_h

#include "AnalyzerCore.h"

class DataPreparation: public AnalyzerCore {
public:
    bool Skim1E2Mu, Skim3Mu;    // flags
    vector<TString> ElectronIDs, MuonIDs;
    vector<TString> DblMuTriggers, EMuTriggers;

    


    DataPreparation();
    ~DataPreparation();
    void initializeAnalyzer();
    void executeEvent();
};

#endif
