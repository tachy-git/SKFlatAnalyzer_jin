#ifndef DiLeptonBase_h
#define DiLeptonBase_h

#include "AnalyzerCore.h"

class DiLeptonBase: public AnalyzerCore {
public:
    bool RunDiMu, RunEMu; // channel
    bool RunSyst;         // systematic run
    vector<TString> ElectronIDs, MuonIDs;
    vector<TString> DblMuTriggers, EMuTriggers;

    TH2D *hMuonIDSF;
    TH2D *hMu17Leg1_Data, *hMu17Leg1_MC;
    TH2D *hMu8Leg2_Data, *hMu8Leg2_MC;
    DiLeptonBase();
    ~DiLeptonBase();

    void initializeAnalyzer();
    void executeEvent();

    double getMuonRecoSF(const Muon &mu, int sys);
    double getMuonIDSF(const Muon &mu, int sys);
    double getTriggerEff(const Muon &mu, TString histkey, bool isDATA, int sys);
    double getDblMuTriggerEff(vector<Muon> &muons, bool isDATA, int sys);
    double getDblMuTriggerSF(vector<Muon> &muons, int sys);
};

#endif