#include "TriLeptonBase.h"

TriLeptonBase::TriLeptonBase() {}
TriLeptonBase::~TriLeptonBase() {
    delete hMuonIDSF;
    delete hMu17Leg1_Data;
    delete hMu17Leg1_MC;
    delete hMu8Leg2_Data;
    delete hMu8Leg2_MC;
    delete hMuonFR;
    delete hElectronFR;
}

void TriLeptonBase::initializeAnalyzer() {
    // flags
    Skim1E2Mu = HasFlag("Skim1E2Mu");
    Skim3Mu = HasFlag("Skim3Mu");
    DenseNet = HasFlag("DenseNet");
    GraphNet = HasFlag("GraphNet");
    ScaleVar = HasFlag("ScaleVar");
    WeightVar = HasFlag("WeightVar");

    // triggers & ID settings
    if (DataEra == "2016preVFP") {
        DblMuTriggers = {
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",
            "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",
        };
        EMuTriggers = {
            "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",
            "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",
            "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v",
            "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v"
        };
        MuonIDs = {"HcToWATight", "HcToWALoose", "HcToWAVeto"};
        ElectronIDs = {"HcToWATight16a", "HcToWALoose16a", "HcToWAVeto16a"};
    }
    else if (DataEra == "2016postVFP") {
        DblMuTriggers = {
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",
            "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",
        };
        EMuTriggers = {
            "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",
            "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",
            "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v",
            "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v"
        };
        MuonIDs = {"HcToWATight", "HcToWALoose", "HcToWAVeto"};
        ElectronIDs = {"HcToWATight16b", "HcToWALoose16b", "HcToWAVeto16b"};
    }
    else if (DataEra == "2017") {
        DblMuTriggers = {
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",
            "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v"
        };
        EMuTriggers = {
            "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",
            "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"
        };
        MuonIDs = {"HcToWATight", "HcToWALoose", "HcToWAVeto"};
        ElectronIDs = {"HcToWATight17", "HcToWALoose17", "HcToWAVeto17"};
    }
    else if (DataEra == "2018") {
        DblMuTriggers = {
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",
            "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v"
        };
        EMuTriggers = {
            "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",
            "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",
            "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v"
        };
        MuonIDs = {"HcToWATight", "HcToWALoose", "HcToWAVeto"};
        ElectronIDs = {"HcToWATight18", "HcToWALoose18", "HcToWAVeto18"};
    }
    else {
        cerr << "[TriLeptonBase::initializeAnalyzer] Wrong era " << DataEra << endl;
        exit(EXIT_FAILURE);
    }

    TString datapath = getenv("DATA_DIR");
    // muonID
    TString muonIDpath = datapath + "/" + GetEra() + "/ID/Muon";
    TFile* fMuonID = new TFile(muonIDpath+"/efficiency_TopHN_IDIso.root");
    hMuonIDSF = (TH2D*)fMuonID->Get("SF_fabs(probe_eta)_probe_pt");
    hMuonIDSF->SetDirectory(0);
    fMuonID->Close();

    // doublemuon trigger
    TFile* fMu17Leg1 = new TFile(muonIDpath+"/efficiency_Mu17Leg1_DoubleMuonTriggers.root");
    hMu17Leg1_Data = (TH2D*)fMu17Leg1->Get("muonEffi_data_fabs(probe_eta)_probe_pt");
    hMu17Leg1_MC = (TH2D*)fMu17Leg1->Get("muonEffi_mc_fabs(probe_eta)_probe_pt");
    hMu17Leg1_Data->SetDirectory(0);
    hMu17Leg1_MC->SetDirectory(0);
    fMu17Leg1->Close();

    TFile* fMu8Leg2 = new TFile(muonIDpath+"/efficiency_Mu8Leg2_DoubleMuonTriggers.root");
    hMu8Leg2_Data = (TH2D*)fMu8Leg2->Get("muonEffi_data_fabs(probe_eta)_probe_pt");
    hMu8Leg2_MC = (TH2D*)fMu8Leg2->Get("muonEffi_mc_fabs(probe_eta)_probe_pt");
    hMu8Leg2_Data->SetDirectory(0);
    hMu8Leg2_MC->SetDirectory(0);
    fMu8Leg2->Close();

    // muon fake rate
    TFile* fMuonFR = new TFile(muonIDpath+"/fakerate_TopHN.root");
    hMuonFR = (TH2D*)fMuonFR->Get("FR_cent_TopHNT_TopHNL");
    hMuonFR->SetDirectory(0);
    fMuonFR->Close();

    // electron fake rate
    hElectronFR = nullptr;

    // Jet tagger
    vector<JetTagging::Parameters> jtps;
    jtps.emplace_back(JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Medium, JetTagging::incl, JetTagging::mujets));
    mcCorr->SetJetTaggingParameters(jtps);
}

void TriLeptonBase::executeEvent() {
    double a;
}

double TriLeptonBase::getMuonIDSF(Muon &mu, int sys) {
    double pt = max(mu.Pt(), 10.);
    double eta = min(fabs(mu.Eta()), 2.4);
    int thisBin = hMuonIDSF->FindBin(eta, pt);
    double value = hMuonIDSF->GetBinContent(thisBin);
    double error = hMuonIDSF->GetBinError(thisBin);
    
    return value + int(sys)*error;
}

double TriLeptonBase::getTriggerEff(Muon &mu, TString histkey, bool isDataEff, int sys) {
    TH2D *h = nullptr;
    double pt = mu.Pt();
    double eta = fabs(mu.Eta());
    if (histkey == "Mu17Leg1" && isDataEff) {
        h = hMu17Leg1_Data;
        if (pt < 16.) pt = 16.;
        if (pt > 200.) pt = 199.;
        if (eta > 2.4) eta = 2.39;
    }
    else if (histkey == "Mu17Leg1" && (!isDataEff)) {
        h = hMu17Leg1_MC;
        if (pt < 16.) pt = 16.;
        if (pt > 200.) pt = 199.;
        if (eta > 2.4) eta = 2.39;
    }
    else if (histkey == "Mu8Leg2" && isDataEff) {
        h = hMu8Leg2_Data;
        if (pt < 10.) pt = 10.;
        if (pt > 200.) pt = 199.;
        if (eta > 2.4) eta = 2.39;
    }
    else if (histkey == "Mu8Leg2" && (!isDataEff)) {
        h = hMu8Leg2_MC;
        if (pt < 10.) pt = 10.;
        if (pt > 200.) pt = 199.;
        if (eta > 2.4) eta = 2.39;
    }
    else {
        cerr << "[TriLeptonBase::getTriggerEff] Wrong combination of histkey and isDataEff" << endl;
        cerr << "[TriLeptonBase::getTriggerEff] histkey = " << histkey << endl;
        cerr << "[TriLeptonBase::getTriggerEff] isDataEff = " << isDataEff << endl;
    }

    int thisBin = h->FindBin(eta, pt);
    double value = h->GetBinContent(thisBin);
    double error = h->GetBinError(thisBin);

    return value + int(sys)*error;
}

double TriLeptonBase::getMuonFakeProb(const Muon &mu, int sys) {
    if (mu.PassID(MuonIDs[0]))
        return 1.;
    double ptCorr = mu.Pt()*(1.+max(0., mu.MiniRelIso()-0.1));
    int thisBin = hMuonFR->FindBin(ptCorr, fabs(mu.Eta()));
    double value = hMuonFR->GetBinContent(thisBin);
    double error = hMuonFR->GetBinError(thisBin);
    
    return value + sys*error;
}

double TriLeptonBase::getElectronFakeProb(const Electron &ele, int sys) {
    if (ele.PassID(ElectronIDs[0]))
        return 1.;
    double ptCorr = ele.Pt()*(1.+max(0., ele.MiniRelIso()-0.1));
    int thisBin = hElectronFR->FindBin(ptCorr, fabs(ele.Eta()));
    double value = hElectronFR->GetBinContent(thisBin);
    double error = hElectronFR->GetBinError(thisBin);

    return value + sys*error;
}

double TriLeptonBase::getFakeWeight(const vector<Muon> &muons, const vector<Electron> &electrons, int sys) {
    double weight = -1.;
    for (const auto &mu: muons) {
        double fr = getMuonFakeProb(mu, sys);
        weight *= -1.*(fr / (1.-fr));
    }
    for (const auto &ele: electrons) {
        double fr = getElectronFakeProb(ele, sys);
        weight *= -1.*(fr / (1.-fr));
    }
    return weight;
}
