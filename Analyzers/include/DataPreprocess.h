#ifndef DataPreprocess_h
#define DataPreprocess_h

#include "AnalyzerCore.h"

class DataPreprocess: public AnalyzerCore {
private:
    TTree *Events;
    // METVector
    float METvPt, METvPhi;

    // muons
    unsigned int nMuons;
    float MuonPtColl[3];
    float MuonEtaColl[3];
    float MuonPhiColl[3];
    float MuonMassColl[3];
    int   MuonChargeColl[3];
    bool  MuonLabelColl[3];

    // electrons
    unsigned int nElectrons;
    float ElectronPtColl[1];
    float ElectronEtaColl[1];
    float ElectronPhiColl[1];
    float ElectronMassColl[1];
    int   ElectronChargeColl[1];
    bool  ElectronLabelColl[1];

    // jets
    unsigned int nJets;
    float JetPtColl[20];
    float JetEtaColl[20];
    float JetPhiColl[20];
    float JetMassColl[20];
    float JetChargeColl[20];
    float JetBtagScoreColl[20];
    bool  JetIsBtaggedColl[20];
    bool  JetLabelColl[20];

public:
    bool Skim1E2Mu, Skim3Mu;
    bool MatchChargedHiggs;
    vector<TString> ElectronIDs, MuonIDs;
    vector<TString> DblMuTriggers, EMuTriggers;

    DataPreprocess();
    ~DataPreprocess();
    void initializeAnalyzer();
    void executeEvent();

    bool isWinDecay(vector<Gen> &chargedDecay);
    vector<Gen> replaceWtoDecays(vector<Gen> &chargedDecay, vector<Gen> &truth);
    vector<Gen> findChargedDecays(vector<Gen> &truth);
};

#endif
