#ifndef TriLeptonBase_h
#define TriLeptonBase_h

#include "AnalyzerCore.h"

class TriLeptonBase: public AnalyzerCore {
public:
    bool Skim1E2Mu, Skim3Mu;    // flags
    vector<TString> ElectronIDs, MuonIDs;
    vector<TString> DblMuTriggers, EMuTriggers;

    TH2D *hMuonIDSF, *hMu17Leg1_Data, *hMu17Leg1_MC, *hMu8Leg2_Data, *hMu8Leg2_MC;
    TriLeptonBase();
    ~TriLeptonBase();
    void initializeAnalyzer();
    void executeEvent();

    double getMuonIDSF(Muon &mu, int sys);
    double getTriggerEff(Muon &mu, TString histkey, bool isDataEff, int sys);
};

#endif
