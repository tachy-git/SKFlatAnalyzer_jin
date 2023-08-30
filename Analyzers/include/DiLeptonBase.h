#ifndef DiLeptonBase_h
#define DiLeptonBase_h

#include "AnalyzerCore.h"

class DiLeptonBase: public AnalyzerCore {
public:
    bool RunDiMu, RunEMu; // channel
    //bool MeasFakeMu, MeasFakeEl; // for fake measurement
    // Divide paths...
    bool MeasFakeMu8, MeasFakeMu17;
    bool MeasFakeEl8, MeasFakeEl12, MeasFakeEl23;
    bool RunSyst;         // systematic run
    vector<TString> ElectronIDs, MuonIDs;
    vector<TString> DblMuTriggers, EMuTriggers;

    TH2D *hMuonIDSF;
    TH2D *hMu17Leg1_Data, *hMu17Leg1_MC;
    TH2D *hMu8Leg2_Data, *hMu8Leg2_MC;
    TH2D *hElIDSF;
    TH2D *hEl23Leg1_Data, *hEl23Leg1_MC;
    TH2D *hEl12Leg2_Data, *hEl12Leg2_MC;
    TH1D *hNPVEl8_Data, *hNPVEl12_Data, *hNPVEl23_Data;
    TH1D *hNPVEl8_MC, *hNPVEl12_MC, *hNPVEl23_MC;
    TH1D *hNPVMu8_Data, *hNPVMu17_Data;
    TH1D *hNPVMu8_MC, *hNPVMu17_MC;
    DiLeptonBase();
    ~DiLeptonBase();

    void initializeAnalyzer();
    void executeEvent();

    double getMuonRecoSF(const Muon &mu, int sys);
    double getMuonIDSF(const Muon &mu, int sys);
    double getEleIDSF(const Electron &ele, int sys);
    double getTriggerEff(const Muon &mu, TString histkey, bool isDATA, int sys);
    double getTriggerEff(const Electron &ele, TString histkey, bool isDATA, int sys);
    double getDblMuTriggerEff(vector<Muon> &muons, bool isDATA, int sys);
    double getDblMuTriggerSF(vector<Muon> &muons, int sys);
    double getEMuTriggerEff(vector<Electron> &electrons, vector<Muon> &muons, bool isDATA, int sys);
    double getEMuTriggerSF(vector<Electron> &electrons, vector<Muon> &muons, int sys);
    double getDZEfficiency(TString SFkey, bool isDATA);
    double getNPVReweight(unsigned int NPV, TString &path);
};

#endif
