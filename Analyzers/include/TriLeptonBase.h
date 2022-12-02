#ifndef TriLeptonBase_h
#define TriLeptonBase_h

#include "AnalyzerCore.h"

class TriLeptonBase: public AnalyzerCore {
public:
    bool Skim1E2Mu, Skim3Mu;    // flags
    vector<TString> ElectronIDs, MuonIDs;
    vector<TString> DblMuTriggers, EMuTriggers;
    TriLeptonBase();
    ~TriLeptonBase();
    void initializeAnalyzer();
    void executeEvent();
};

#endif
