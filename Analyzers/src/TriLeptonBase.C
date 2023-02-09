#include "TriLeptonBase.h"

TriLeptonBase::TriLeptonBase() {}
TriLeptonBase::~TriLeptonBase() {
    delete hMuonIDSF;
    delete hMu17Leg1_Data;
    delete hMu17Leg1_MC;
    delete hMu8Leg2_Data;
    delete hMu8Leg2_MC;
    delete hMuonFR;
    delete hMuonFRUp;
    delete hMuonFRDown;
    delete hElectronFR;
    delete hElectronFRUp;
    delete hElectronFRDown;
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
    hMuonFRUp = (TH2D*)fMuonFR->Get("FRErrUp_Tot_TopHNT_TopHNL");
    hMuonFRDown = (TH2D*)fMuonFR->Get("FRErrDown_Tot_TopHNT_TopHNL");
    hMuonFR->SetDirectory(0);
    hMuonFRUp->SetDirectory(0);
    hMuonFRDown->SetDirectory(0);
    fMuonFR->Close();

    // electron fake rate
    hElectronFR = nullptr;
    hElectronFRUp = nullptr;
    hElectronFRDown = nullptr;

    // Jet tagger
    vector<JetTagging::Parameters> jtps;
    jtps.emplace_back(JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Medium, JetTagging::incl, JetTagging::mujets));
    mcCorr->SetJetTaggingParameters(jtps);
}

void TriLeptonBase::executeEvent() {

    if (! PassMETFilter()) return;
    
    // object definition
    Event ev = GetEvent();
    vector<Muon>     rawMuons = GetAllMuons();
    vector<Electron> rawElectrons = GetAllElectrons();
    vector<Jet>      rawJets = GetAllJets();
    Particle         METv = ev.GetMETVector();
    vector<Gen>      truth = GetGens();

    vector<Muon> vetoMuons = SelectMuons(rawMuons, MuonIDs.at(2), 10., 2.4);
    vector<Muon> tightMuons = SelectMuons(vetoMuons, MuonIDs.at(0), 10., 2.4);
    vector<Electron> vetoElectrons = SelectElectrons(rawElectrons, ElectronIDs.at(2), 10., 2.5);
    vector<Electron> tightElectrons = SelectElectrons(vetoElectrons, ElectronIDs.at(0), 10., 2.5);
    vector<Jet> jets = SelectJets(rawJets, "tight", 20., 2.4);
    jets = JetsVetoLeptonInside(jets, vetoElectrons, vetoMuons, 0.4);
    vector<Jet> bjets;
    for (const auto &j: jets) {
        const double btagScore = j.GetTaggerResult(JetTagging::DeepJet);
        const double wp = mcCorr->GetJetTaggingCutValue(JetTagging::DeepJet, JetTagging::Medium);
        if (btagScore > wp) bjets.emplace_back(j);
    }
    
    std::sort(vetoMuons.begin(), vetoMuons.end(), PtComparing);
    std::sort(tightMuons.begin(), tightMuons.end(), PtComparing);
    std::sort(vetoElectrons.begin(), vetoElectrons.end(), PtComparing);
    std::sort(tightElectrons.begin(), tightElectrons.end(), PtComparing);
    std::sort(jets.begin(), jets.end(), PtComparing);
    std::sort(bjets.begin(), bjets.end(), PtComparing);

    // baseline event selection
    const bool is3Mu = (tightMuons.size() == 3 && vetoMuons.size() == 3 && tightElectrons.size() == 0 && vetoElectrons.size() == 0);
    const bool is1E2Mu = (tightMuons.size() == 2 && vetoMuons.size() == 2 && tightElectrons.size() == 1 && vetoElectrons.size() == 1);

    if (Skim1E2Mu && (! is1E2Mu)) return;
    if (Skim3Mu && (! is3Mu)) return;

    // prompt matching
    vector<Muon> promptMuons;
    vector<Electron> promptElectrons;
    for (const auto &mu: tightMuons) 
        if (GetLeptonType(mu, truth) > 0) promptMuons.emplace_back(mu);
    for (const auto &ele: tightElectrons)
        if (GetLeptonType(ele, truth) > 0)  promptElectrons.emplace_back(ele);

    if (! (promptMuons.size() == tightMuons.size())) return;
    if (! (promptElectrons.size() == tightElectrons.size())) return;

    TString channel = "";
    // 1E2Mu Baseline
    // 1. pass EMuTriggers
    // 2. Exact 2 tight muons and 1 tight electrons, no additional leptons
    // 3. Exists OS muon pair wiht mass > 12 GeV
    // 4. At least two jets, one bjet
    if (Skim1E2Mu) {
        if (! ev.PassTrigger(EMuTriggers)) return;
        const bool passLeadMu = tightMuons.at(0).Pt() > 25. && tightElectrons.at(0).Pt() > 15.;
        const bool passLeadEle = tightMuons.at(0).Pt() > 10. && tightElectrons.at(0).Pt() > 25.;
        if (! (passLeadMu || passLeadEle)) return;
        if (! (tightMuons.at(0).Charge() + tightMuons.at(1).Charge() == 0)) return;
        Particle pair = tightMuons.at(0) + tightMuons.at(1);
        if (! (pair.M() > 12.)) return;
        if (! (jets.size() >= 2)) return;
        if (! (bjets.size() >= 1)) return;
        channel = "SR1E2Mu";
    }
    // 3Mu baseline
    // 1. pass DblMuTriggers
    // 2. Exact 3 tight muons, no additional leptons
    // 3. Exist OS muon pair,
    // 4. All OS muon pair mass > 12 GeV
    // 5. At least two jets and one bjet
    else if (Skim3Mu) {
        if (! ev.PassTrigger(DblMuTriggers)) return;
        const Muon &mu1 = tightMuons.at(0);
        const Muon &mu2 = tightMuons.at(1);
        const Muon &mu3 = tightMuons.at(2);
        if (! (mu1.Pt() > 20.)) return;
        if (! (fabs(mu1.Charge()+mu2.Charge()+mu3.Charge()) == 1)) return;
        
        Particle pair1, pair2;
        if (mu1.Charge() == mu2.Charge()) {
            pair1 = mu1 + mu3;
            pair2 = mu2 + mu3;
        }
        else if (mu1.Charge() == mu3.Charge()) {
            pair1 = mu1 + mu2;
            pair2 = mu2 + mu3;
        }
        else {  // mu2.Charge() == mu3.Charge()
            pair1 = mu1 + mu2;
            pair2 = mu1 + mu3;
        }
        if (! (pair1.M() > 12.)) return;
        if (! (pair2.M() > 12.)) return;
        if (! (jets.size() >= 2)) return;
        if (! (bjets.size() >= 1)) return;
        channel = "SR3Mu";
    }
    else {
        cerr << "[TriLeptonBase::ExecuteEvent] Please set flag for the channel" << endl;
        exit(EXIT_FAILURE);
    }


    // Now signal region distributions
    double weight = 1.;
    weight *= MCweight();
    weight *= ev.GetTriggerLumi("Full");
    weight *= GetPrefireWeight(0);
    weight *= GetPileUpWeight(nPileUp, 0);

    // Fill objects
    for (unsigned int i = 0; i < tightMuons.size(); i++) {
        TString histkey = channel+"/muons/"+TString::Itoa(i+1, 10);
        Muon &mu = tightMuons.at(i);
        FillHist(histkey+"/pt", mu.Pt(), weight, 300, 0., 300.);
        FillHist(histkey+"/eta", mu.Eta(), weight, 48, -2.4, 2.4);
        FillHist(histkey+"/phi", mu.Phi(), weight, 64, -3.2, 3.2);
    }
    for (unsigned int i = 0; i < tightElectrons.size(); i++) {
        TString histkey = channel+"/electrons/"+TString::Itoa(i+1, 10);
        Electron &ele = tightElectrons.at(i);
        FillHist(histkey+"/pt", ele.Pt(), weight, 300, 0., 300.);
        FillHist(histkey+"/eta", ele.Eta(), weight, 48, -2.4, 2.4);
        FillHist(histkey+"/phi", ele.Phi(), weight, 64, -3.2, 3.2);
    }
    for (unsigned int i = 0; i < jets.size(); i++) {
        TString histkey = channel+"/jets/"+TString::Itoa(i+1, 10);
        Jet &jet = jets.at(i);
        FillHist(histkey+"/pt", jet.Pt(), weight, 300, 0., 300.);
        FillHist(histkey+"/eta", jet.Eta(), weight, 48, -2.4, 2.4);
        FillHist(histkey+"/phi", jet.Phi(), weight, 64, -3.2, 3.2);
    }
    for (unsigned int i = 0; i < bjets.size(); i++) {
        TString histkey = channel+"/bjets/"+TString::Itoa(i+1, 10);
        Jet &bjet = bjets.at(i);
        FillHist(histkey+"/pt", bjet.Pt(), weight, 300, 0., 300.);
        FillHist(histkey+"/eta", bjet.Eta(), weight, 48, -2.4, 2.4);
        FillHist(histkey+"/phi", bjet.Phi(), weight, 64, -3.2, 3.2);
    }
    FillHist(channel+"/jets/size", jets.size(), weight, 30, 0., 30.);
    FillHist(channel+"/bjets/size", bjets.size(), weight, 20, 0., 20.);
    FillHist(channel+"/MissingPT", METv.Pt(), weight, 300, 0., 300.);
    FillHist(channel+"/MissingPhi", METv.Phi(), weight, 64, -3.2, 3.2);


    // Now signal study
    // First find two signal muons from A
    if (Skim3Mu) {
        vector<Muon> signalMuons, nonSignalMuons;

        for (const auto &mu: tightMuons) {
            if (GetLeptonType(mu, truth) == 2) signalMuons.emplace_back(mu);
            else                               nonSignalMuons.emplace_back(mu);
        }

        if (! (signalMuons.size() == 2 && nonSignalMuons.size() == 1)) return;
        Muon promptMu = nonSignalMuons.at(0);
        Muon signalMuSS, signalMuOS;
        if (promptMu.Charge() == signalMuons.at(0).Charge()) {
            signalMuSS = signalMuons.at(0);
            signalMuOS = signalMuons.at(1);
        }
        else {
            signalMuSS = signalMuons.at(1);
            signalMuOS = signalMuons.at(0);
        }
    
        FillHist(channel+"/signalMuSS/pt", signalMuSS.Pt(), weight, 300, 0., 300.);
        FillHist(channel+"/signalMuSS/eta", signalMuSS.Eta(), weight, 48, -2.4, 2.4);
        FillHist(channel+"/signalMuSS/phi", signalMuSS.Phi(), weight, 64, -3.2, 3.2);
        FillHist(channel+"/signalMuSS/MT", (signalMuSS+METv).Mt(), weight, 300, 0., 300.);
        FillHist(channel+"/signalMuOS/pt", signalMuOS.Pt(), weight, 300, 0., 300.);
        FillHist(channel+"/signalMuOS/eta", signalMuOS.Eta(), weight, 48, -2.4, 2.4);
        FillHist(channel+"/signalMuOS/phi", signalMuOS.Phi(), weight, 64, -3.2, 3.2);
        FillHist(channel+"/promptMu/pt", promptMu.Pt(), weight, 300, 0., 300.);
        FillHist(channel+"/promptMu/eta", promptMu.Eta(), weight, 48, -2.4, 2.4);
        FillHist(channel+"/promptMu/phi", promptMu.Phi(), weight, 64, -3.2, 3.2); 
        FillHist(channel+"/promptMu/MT", (promptMu+METv).Mt(), weight, 300, 0., 300.);
        FillHist(channel+"/deltaR/promptPair", signalMuOS.DeltaR(promptMu), weight, 100, 0., 5.);
        FillHist(channel+"/deltaR/signalPair", signalMuOS.DeltaR(signalMuSS), weight, 100, 0., 5.);
    }

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

double TriLeptonBase::getDblMuTriggerEff(vector<Muon> &muons, bool isDATA, int sys) {
    // check no. of muons
    if (! (muons.size() == 3)) {
        cerr << "[TriLeptonBase::getDblMuTriggerEff] Wrong no. of muons " << muons.size() << endl;
        exit(EXIT_FAILURE);
    }
    
    Muon &mu1 = muons.at(0);
    Muon &mu2 = muons.at(1);
    Muon &mu3 = muons.at(2);


    double case1 = getTriggerEff(mu1, "Mu17Leg1", isDATA, sys);
           case1 *= getTriggerEff(mu2, "Mu8Leg2", isDATA, sys);
    double case2 = 1. - getTriggerEff(mu1, "Mu17Leg1", isDATA, sys);
           case2 *= getTriggerEff(mu2, "Mu17Leg1", isDATA, sys);
           case2 *= getTriggerEff(mu3, "Mu8Leg2", isDATA, sys);
    double case3 = getTriggerEff(mu1, "Mu17Leg1", isDATA, sys);
           case3 *= 1. - getTriggerEff(mu2, "Mu8Leg2", isDATA, sys);
           case3 *= getTriggerEff(mu3, "Mu8Leg2", isDATA, sys);
    return case1 + case2 + case3;
}

double TriLeptonBase::getDblMuTriggerSF(vector<Muon> &muons, int sys) {
   double effData = getDblMuTriggerEff(muons, true, sys); 
   double effMC   = getDblMuTriggerEff(muons, false, sys);
   if (effMC == 0 || effData == 0)
       return 1.;

   return effData / effMC;
}

double TriLeptonBase::getMuonFakeProb(const Muon &mu, int sys) {
    double ptCorr = mu.Pt()*(1.+max(0., mu.MiniRelIso()-0.1));
    double absEta = fabs(mu.Eta());
    if (ptCorr < 10.) ptCorr = 10.;
    if (ptCorr > 50.) ptCorr = 49.9;
    if (absEta > 2.4) absEta = 2.399;

    int thisBin = hMuonFR->FindBin(ptCorr, absEta);
    double value = hMuonFR->GetBinContent(thisBin);
    double error = 0.;
    if (sys == 1) error = hMuonFRUp->GetBinContent(thisBin);
    if (sys == -1) error = hMuonFRDown->GetBinContent(thisBin);
    
    return value + error;
}

double TriLeptonBase::getElectronFakeProb(const Electron &ele, int sys) {
    double ptCorr = ele.Pt()*(1.+max(0., ele.MiniRelIso()-0.1));
    int thisBin = hElectronFR->FindBin(ptCorr, fabs(ele.Eta()));
    double value = hElectronFR->GetBinContent(thisBin);
    double error = 0.;
    if (sys == 1) error = hElectronFRUp->GetBinError(thisBin);
    if (sys == -1) error = hElectronFRDown->GetBinError(thisBin);

    return value + error;
}

double TriLeptonBase::getFakeWeight(const vector<Muon> &muons, const vector<Electron> &electrons, int sys) {
    double weight = -1.;
    for (const auto &mu: muons) {
        if (mu.PassID(MuonIDs.at(0))) continue;
        double fr = getMuonFakeProb(mu, sys);
        weight *= -1.*(fr / (1.-fr));
    }
    for (const auto &ele: electrons) {
        if (ele.PassID(ElectronIDs.at(0))) continue;
        double fr = getElectronFakeProb(ele, sys);
        weight *= -1.*(fr / (1.-fr));
    }
    return weight;
}
