#ifndef TriLeptonBase_h
#define TriLeptonBase_h

#include "AnalyzerCore.h"

class TriLeptonBase: public AnalyzerCore {
public:
    bool Skim1E2Mu, Skim3Mu;    // channel flags 
    bool DenseNet, GraphNet;    // network flags
    bool ScaleVar, WeightVar;   // systematics flags
    bool RunSyst;
    bool FakeStudy;             // for AcceptanceStudy, if FakeStudy is true, reverse prompt matching will be performed
    vector<TString> ElectronIDs, MuonIDs;
    vector<TString> DblMuTriggers, EMuTriggers;
    vector<TString> MASSPOINTs;

    TH2D *hMuonIDSF;
    TH2D *hMu17Leg1_Data, *hMu17Leg1_MC;
    TH2D *hMu8Leg2_Data, *hMu8Leg2_MC;
    TH2D *hElIDSF;
    TH2D *hEl23Leg1_Data, *hEl23Leg1_MC;
    TH2D *hEl12Leg2_Data, *hEl12Leg2_MC;
    TH2D *hMuFR, *hElFR;
    TriLeptonBase();
    ~TriLeptonBase();
    void initializeAnalyzer();
    void executeEvent();

    double getMuonRecoSF(const Muon &mu, int sys);
    double getMuonIDSF(const Muon &mu, int sys);
    double getEleIDSF(const Electron &ele, int sys);
    double getTriggerEff(const Muon &mu, TString histkey, bool isDataEff, int sys);
    double getTriggerEff(const Electron &ele, TString histkey, bool isDATA, int sys);
    double getDblMuTriggerEff(vector<Muon> &muons, bool isDATA, int sys);
    double getDblMuTriggerSF(vector<Muon> &muons, int sys);
    double getEMuTriggerEff(vector<Electron> &electrons, vector<Muon> &muons, bool isDATA, int sys);
    double getEMuTriggerSF(vector<Electron> &electrons, vector<Muon> &muons, int sys);
    double getDZEfficiency(TString SFkey, bool isDATA);
    double getMuonFakeProb(const Muon &mu, int sys);
    double getElectronFakeProb(const Electron &ele, int sys);
    double getFakeWeight(const vector<Muon> &muons, const vector<Electron> &electrons, int sys);
};

#endif
