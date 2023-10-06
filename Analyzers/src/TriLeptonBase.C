#include "TriLeptonBase.h"

TriLeptonBase::TriLeptonBase() {}
TriLeptonBase::~TriLeptonBase() {
    delete hMuonIDSF;
    delete hMu17Leg1_Data;
    delete hMu17Leg1_MC;
    delete hMu8Leg2_Data;
    delete hMu8Leg2_MC;
    delete hMuFR;
    delete hElFR;
}

void TriLeptonBase::initializeAnalyzer() {
    // flags
    Skim1E2Mu = HasFlag("Skim1E2Mu");
    Skim3Mu = HasFlag("Skim3Mu");
    DenseNet = HasFlag("DenseNet");
    GraphNet = HasFlag("GraphNet");
    ScaleVar = HasFlag("ScaleVar");
    WeightVar = HasFlag("WeightVar");
    FakeStudy = HasFlag("FakeStudy");

    // triggers & ID settings
    if (DataEra == "2016preVFP") {
        DblMuTriggers = {
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v",
            "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v",
        };
        EMuTriggers = {
            "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v",
            "HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v",
        };
        MuonIDs = {"HcToWATight", "HcToWALoose", "HcToWAVeto"};
        ElectronIDs = {"HcToWATight16a", "HcToWALoose16a", "HcToWAVeto16a"};
    }
    else if (DataEra == "2016postVFP") {
        DblMuTriggers = {
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",
            "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",
            "HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v"
        };
        EMuTriggers = {
            "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",
            "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v"
        };
        MuonIDs = {"HcToWATight", "HcToWALoose", "HcToWAVeto"};
        ElectronIDs = {"HcToWATight16b", "HcToWALoose16b", "HcToWAVeto16b"};
    }
    else if (DataEra == "2017") {
        DblMuTriggers = {
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",
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
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v"
        };
        EMuTriggers = {
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

    MASSPOINTs = {"MHc-70_MA-15", "MHc-70_MA-40", "MHc-70_MA-65",
                  "MHc-100_MA-15", "MHc-100_MA-60", "MHc-100_MA-95",
                  "MHc-130_MA-15", "MHc-130_MA-55", "MHc-130_MA-90", "MHc-130_MA-125",
                  "MHc-160_MA-15", "MHc-160_MA-85", "MHc-160_MA-120", "MHc-160_MA-155"};

    TString datapath = getenv("DATA_DIR");
    // muon ID
    TString muonIDpath = datapath + "/" + GetEra() + "/ID/Muon";
    TFile *fMuonID = new TFile(muonIDpath+"/efficiency_TopHNT.root");
    hMuonIDSF = (TH2D*)fMuonID->Get("sf");  hMuonIDSF->SetDirectory(0);
    fMuonID->Close();

    // muon trigger legs
    TFile *fMu17Leg1 = new TFile(muonIDpath+"/efficiency_Mu17Leg1.root");
    hMu17Leg1_Data = (TH2D*)fMu17Leg1->Get("data"); hMu17Leg1_Data->SetDirectory(0);
    hMu17Leg1_MC = (TH2D*)fMu17Leg1->Get("sim");    hMu17Leg1_MC->SetDirectory(0);
    fMu17Leg1->Close();

    TFile* fMu8Leg2 = new TFile(muonIDpath+"/efficiency_Mu8Leg2.root");
    hMu8Leg2_Data = (TH2D*)fMu8Leg2->Get("data"); hMu8Leg2_Data->SetDirectory(0);
    hMu8Leg2_MC = (TH2D*)fMu8Leg2->Get("sim");    hMu8Leg2_MC->SetDirectory(0);
    fMu8Leg2->Close();

    // ele ID
    TString eleIDPath = datapath + "/" + GetEra() + "/ID/Electron";
    TFile *fEleID = new TFile(eleIDPath+"/efficiency_TopHNT.root");
    hElIDSF = (TH2D*)fEleID->Get("sf"); hElIDSF->SetDirectory(0);
    fEleID->Close();

    TFile *fEl23Leg1 = new TFile(eleIDPath+"/efficiency_El23Leg1.root");
    hEl23Leg1_Data = (TH2D*)fEl23Leg1->Get("data"); hEl23Leg1_Data->SetDirectory(0);
    hEl23Leg1_MC = (TH2D*)fEl23Leg1->Get("sim");    hEl23Leg1_MC->SetDirectory(0);
    fEl23Leg1->Close();

    TFile *fEl12Leg2 = new TFile(eleIDPath+"/efficiency_El12Leg2.root");
    hEl12Leg2_Data = (TH2D*)fEl12Leg2->Get("data"); hEl12Leg2_Data->SetDirectory(0);
    hEl12Leg2_MC = (TH2D*)fEl12Leg2->Get("sim");    hEl12Leg2_MC->SetDirectory(0);
    fEl12Leg2->Close();

    // muon fake rate
    if (FakeStudy) {
        TFile* fMuFR = new TFile(muonIDpath+"/fakerate_TopHNT_TopHNL_qcd.root");
        hMuFR = (TH2D*)fMuFR->Get("fakerate");
        hMuFR->SetDirectory(0);
        fMuFR->Close();

        // electron fake rate
        TFile* fElFR = new TFile(eleIDPath+"/fakerate_TopHNT_TopHNL_qcd.root");
        hElFR = (TH2D*)fElFR->Get("fakerate");
        hElFR->SetDirectory(0);
        fElFR->Close();
    }
    else {
        TFile* fMuFR = new TFile(muonIDpath+"/fakerate_TopHNT_TopHNL.root");
        hMuFR = (TH2D*)fMuFR->Get("fakerate");
        hMuFR->SetDirectory(0);
        fMuFR->Close();

        // electron fake rate
        TFile* fElFR = new TFile(eleIDPath+"/fakerate_TopHNT_TopHNL.root");
        hElFR = (TH2D*)fElFR->Get("fakerate");
        hElFR->SetDirectory(0);
        fElFR->Close();
    }

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
        Particle sigPair = signalMuOS + signalMuSS;
        Particle bkgPair = signalMuOS + promptMu;
        
        Particle lowPair = sigPair.M() < bkgPair.M() ? sigPair : bkgPair;
        Particle highPair = sigPair.M() > bkgPair.M() ? sigPair : bkgPair;
        bool isSignalMassLow = sigPair.M() < bkgPair.M();
        FillHist(channel+"/lowPair/pt", lowPair.Pt(), weight, 300, 0., 300.);
        FillHist(channel+"/lowPair/eta", lowPair.Eta(), weight, 100, -5., 5.);
        FillHist(channel+"/lowPair/phi", lowPair.Phi(), weight, 64, -3.2, 3.2);
        FillHist(channel+"/lowPair/mass", lowPair.M(), weight, 300, 0., 300.);
        FillHist(channel+"/highPair/pt", highPair.Pt(), weight, 300, 0., 300.);
        FillHist(channel+"/highPair/eta", highPair.Eta(), weight, 100, -5., 5.);
        FillHist(channel+"/highPair/phi", highPair.Phi(), weight, 64, -3.2, 3.2);
        FillHist(channel+"/highPair/mass", highPair.M(), weight, 300, 0., 300.);
        FillHist(channel+"/isSignalMassLow", isSignalMassLow, weight, 2, 0., 2.);
    }

}

double TriLeptonBase::getMuonRecoSF(const Muon &mu, int sys) {
    const double abseta = fabs(mu.Eta());

    if (DataEra == "2016preVFP") {
        if (abseta < 0.9) {
            const double value = 0.9998229551300333;
            const double stat = 0.0001538802103231026;
            const double syst = 0.0003540897399334497;
            return value + sys*sqrt(pow(stat, 2)+pow(syst, 2));
        }
        else if (abseta < 1.2) {
            const double value = 1.0001593416915515;
            const double stat = 0.00019861903120026457;
            const double syst = 0.00031024592139106355;
            return value + sys*sqrt(pow(stat, 2)+pow(syst, 2)); 
        }
        else if (abseta < 2.1) {
            const double value = 0.9998936144006075;
            const double stat = 0.00012188589514012365;
            const double syst = 0.00021277119878493345;
            return value + sys*sqrt(pow(stat, 2)+pow(syst, 2)); 
        }
        else if (abseta < 2.4) {
            const double value = 0.9990268820042745;
            const double stat = 0.00027638902644996395;
            const double syst = 0.0019462359914510508;
            return value + sys*sqrt(pow(stat, 2)+pow(syst, 2)); 
        }
        else {
            cerr << "[TriLeptonBase::getMuonRecoSF] wrong muon eta value " << abseta << endl;
            exit(EXIT_FAILURE); 
        }
    }
    else if (DataEra == "2016postVFP") {
        if (abseta < 0.9) {
            const double value = 1.0000406419782646;
            const double stat = 0.00010260291858070426;
            const double syst = 0.0014366927652431664;
            return value + sys*sqrt(pow(stat, 2)+pow(syst, 2));
        }
        else if (abseta < 1.2) {
            const double value = 0.9997959311146515;
            const double stat = 0.00019912837537507789;
            const double syst = 0.0010917857343065423;
            return value + sys*sqrt(pow(stat, 2)+pow(syst, 2));
        }
        else if (abseta < 2.1) {
            const double value = 0.9994928400570587;
            const double stat = 0.00012513847429973846;
            const double syst = 0.0014814654032937547;
            return value + sys*sqrt(pow(stat, 2)+pow(syst, 2));
        }
        else if (abseta < 2.4) {
            const double value = 0.9990728619505579;
            const double stat = 0.0002754474704705526;
            const double syst = 0.0017364778744567663;
            return value + sys*sqrt(pow(stat, 2)+pow(syst, 2));
        }
        else {
            cerr << "[TriLeptonBase::getMuonRecoSF] wrong muon eta value " << abseta << endl;
            exit(EXIT_FAILURE);
        }
    }
    else if (DataEra == "2017") {
        if (abseta < 0.9) {
            const double value = 0.9996742562806361;
            const double stat = 7.650191371261136e-05;
            const double syst = 0.0006514874387277825;
            return value + sys*sqrt(pow(stat, 2)+pow(syst, 2));
        }
        else if (abseta < 1.2) {
            const double value = 0.9997813602035737;
            const double stat = 0.00014496238686164667;
            const double syst = 0.0004372795928526685;
            return value + sys*sqrt(pow(stat, 2)+pow(syst, 2));
        }
        else if (abseta < 2.1) {
            const double value = 0.9994674742459532;
            const double stat = 7.739510750489317e-05;
            const double syst = 0.0010650515080936618;
            return value + sys*sqrt(pow(stat, 2)+pow(syst, 2));
        }
        else if (abseta < 2.4) {
            const double value = 0.9993566412630517;
            const double stat = 0.00022835790507860388;
            const double syst = 0.0011810962222705494;
            return value + sys*sqrt(pow(stat, 2)+pow(syst, 2));
        }
        else {
            cerr << "[TriLeptonBase::getMuonRecoSF] wrong muon eta value " << abseta << endl;
            exit(EXIT_FAILURE);
        }
    }
    else if (DataEra == "2018") {
        if (abseta < 0.9) {
            const double value = 0.9998088006315689;
            const double stat = 6.498845788247257e-05;
            const double syst = 0.0003823987368622994;
            return value + sys*sqrt(pow(stat, 2)+pow(syst, 2));
        }
        else if (abseta < 1.2) {
            const double value = 0.999754701980269;
            const double stat = 0.00011054079511271507;
            const double syst = 0.0005221124230931915;
            return value + sys*sqrt(pow(stat, 2)+pow(syst, 2));
        }
        else if (abseta < 2.1) {
            const double value = 0.9995842791862117;
            const double stat = 7.574443994874554e-05;
            const double syst = 0.0008314416275765346;
            return value + sys*sqrt(pow(stat, 2)+pow(syst, 2));
        }
        else if (abseta < 2.4) {
            const double value = 0.9990341741614288;
            const double stat = 0.00019911479235592246;
            const double syst = 0.0017237408292350668;
            return value + sys*sqrt(pow(stat, 2)+pow(syst, 2));
        }
        else {
            cerr << "[TriLeptonBase::getMuonRecoSF] wrong muon eta value " << abseta << endl;
            exit(EXIT_FAILURE);
        }
    }
    else {
        cerr << "[TriLeptonBase::getMuonRecoSF] not implemented era" << DataEra << endl;
        exit(EXIT_FAILURE);
    }
}

double TriLeptonBase::getMuonIDSF(const Muon &mu, int sys) {
    double pt = mu.Pt();
    double eta = fabs(mu.Eta());
    if (pt < 10.) pt = 10;
    if (pt >= 200) pt = 199.;
    if (eta > 2.4) eta = 2.39;
    int thisBin = hMuonIDSF->FindBin(eta, pt);
    double value = hMuonIDSF->GetBinContent(thisBin);
    double error = hMuonIDSF->GetBinError(thisBin);
    
    return value + int(sys)*error;
}

double TriLeptonBase::getEleIDSF(const Electron &ele, int sys) {
    double pt = ele.Pt();
    double eta = ele.scEta();
    if (pt < 10.) pt = 10.;
    if (pt >= 500.) pt = 499.;
    if (eta < -2.5) eta = -2.499;
    if (eta > 2.5) eta = 2.499;

    int thisBin = hElIDSF->FindBin(eta, pt);
    double value = hElIDSF->GetBinContent(thisBin);
    double error = hElIDSF->GetBinError(thisBin);

    return value + float(sys)*error;
}


double TriLeptonBase::getTriggerEff(const Muon &mu, TString histkey, bool isDataEff, int sys) {
    TH2D *h = nullptr;
    double pt = mu.Pt();
    double eta = fabs(mu.Eta());
    if (histkey == "Mu17Leg1" && isDataEff) {
        h = hMu17Leg1_Data;
        if (pt < 20.) pt = 20.;
        if (pt > 200.) pt = 199.;
        if (eta > 2.4) eta = 2.39;
    }
    else if (histkey == "Mu17Leg1" && (!isDataEff)) {
        h = hMu17Leg1_MC;
        if (pt < 20.) pt = 20.;
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

double TriLeptonBase::getTriggerEff(const Electron &ele, TString histkey, bool isDATA, int sys) {
    TH2D *h = nullptr;
    double pt = ele.Pt();
    double eta = ele.scEta();
    if (eta < -2.5) eta = -2.499;
    if (eta >= 2.5) eta = 2.499;
    if (histkey == "El23Leg1" && isDATA) {
        h = hEl23Leg1_Data;
        if (pt < 25.) pt = 25.;
        if (pt >= 200.) pt = 199.;
    }
    else if (histkey == "El23Leg1" && (!isDATA)) {
        h = hEl23Leg1_MC;
        if (pt < 25.) pt = 25.;
        if (pt > 200.) pt = 199.;
    }
    else if (histkey == "El12Leg2" && isDATA) {
        h = hEl12Leg2_Data;
        if (pt < 15.) pt = 15.;
        if (pt > 200.) pt = 199.;
    }
    else if (histkey == "El12Leg2" && (!isDATA)) {
        h = hEl12Leg2_MC;
        if (pt < 15.) pt = 15.;
        if (pt > 200.) pt = 199.;
    }
    else {
        cerr << "[TriLeptonBase::getTriggerEff] Wrong combination of histkey and isDataEff" << endl;
        cerr << "[TriLeptonBase::getTriggerEff] histkey = " << histkey << endl;
        cerr << "[TriLeptonBase::getTriggerEff] isDATA = " << isDATA << endl;
    }
    int thisBin = h->FindBin(eta, pt);
    double value = h->GetBinContent(thisBin);
    double error = h->GetBinError(thisBin);

    return value + float(sys)*error;
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

// WARNING: Mu23 leg is not measured! Temporarily use Mu17 leg
double TriLeptonBase::getEMuTriggerEff(vector<Electron> &electrons, vector<Muon> &muons, bool isDATA, int sys) {
    // check no. of electrons and muons
    if (! (electrons.size() == 1 && muons.size() == 2)) {
        cerr << "[TriLeptonBase::getEMuTriggerEff] Wrong no. of leptons " << electrons.size() << " " << muons.size() << endl;
        exit(EXIT_FAILURE);
    }

    Electron &ele = electrons.at(0);
    Muon     &mu1 = muons.at(0);
    Muon     &mu2 = muons.at(1);

    double case1 = getTriggerEff(mu1, "Mu8Leg2", isDATA, sys) + (1.-getTriggerEff(mu1, "Mu8Leg2", isDATA, sys)*getDZEfficiency("EMu", isDATA))*getTriggerEff(mu2, "Mu8Leg2", isDATA, sys);
    double case2 = mu2.Pt() > 25. ? getTriggerEff(mu1, "Mu17Leg1", isDATA, sys) + (1.-getTriggerEff(mu1, "Mu17Leg1", isDATA, sys)*getDZEfficiency("EMu", isDATA))*getTriggerEff(mu2, "Mu17Leg1", isDATA, sys) 
                                  : getTriggerEff(mu1, "Mu17Leg1", isDATA, sys); 

    double eff_el = (mu1.Pt() > 25. || mu2.Pt() > 25.) ? getTriggerEff(ele, "El12Leg2", isDATA, sys) : getTriggerEff(ele, "El23Leg1", isDATA, sys);
    double eff_mu = (ele.Pt() > 25.) ? case1 : case2;
    return eff_el * eff_mu;
}

double TriLeptonBase::getDZEfficiency(TString SFKey, bool isDATA) {
    double eff = 0.;
    if (SFKey.Contains("DiMu")) {
        if (DataEra == "2016postVFP") eff = isDATA ? 0.9798 : 0.9969;
        else if (DataEra == "2017")   eff = 0.9958;
        else                          eff = 1.;
    }
    else if (SFKey.Contains("DiElIso")) {
        if (DataEra=="2016preVFP")       eff = 0.986;
        else if (DataEra=="2016postVFP") eff = 0.980;
        else                             eff = 1.;
    }
    else if (SFKey.Contains("EMu")){
        if(DataEra=="2016postVFP") eff = isDATA ? 0.9648:0.9882;
        //else if(DataEra=="2017"  ) Eff = 0.9951; //for now included in muleg
        else                       eff = 1.;
    }
    else {
        eff = 1.;
    }

    return eff;
}

double TriLeptonBase::getDblMuTriggerSF(vector<Muon> &muons, int sys) {
   double effData = getDblMuTriggerEff(muons, true, sys); 
   double effMC   = getDblMuTriggerEff(muons, false, -1*sys);
   if (effMC == 0 || effData == 0)
       return 1.;

   return effData / effMC;
}

double TriLeptonBase::getEMuTriggerSF(vector<Electron> &electrons, vector<Muon> &muons, int sys) {
    double effData = getEMuTriggerEff(electrons, muons, true, sys);
    double effMC = getEMuTriggerEff(electrons, muons, false, -1*sys);
    if (effMC == 0 || effData == 0) return 1.;
    return effData / effMC;
}

double TriLeptonBase::getMuonFakeProb(const Muon &mu, int sys) {
    double ptCorr = mu.Pt()*(1.+max(0., mu.MiniRelIso()-0.1));
    double absEta = fabs(mu.Eta());
    if (ptCorr < 10.) ptCorr = 10.;
    if (ptCorr > 50.) ptCorr = 49.9;
    if (absEta > 2.4) absEta = 2.399;

    int thisBin = hMuFR->FindBin(absEta, ptCorr);
    double value = hMuFR->GetBinContent(thisBin);
    double error = hMuFR->GetBinError(thisBin);
    
    return value + double(sys)*error;
}

double TriLeptonBase::getElectronFakeProb(const Electron &ele, int sys) {
    double ptCorr = ele.Pt()*(1.+max(0., ele.MiniRelIso()-0.1));
    double absEta = fabs(ele.scEta());
    if (ptCorr < 10.) ptCorr = 10.;
    if (ptCorr > 50.) ptCorr = 49.9;
    if (absEta > 2.5) absEta = 2.499;
    int thisBin = hElFR->FindBin(absEta, ptCorr);
    double value = hElFR->GetBinContent(thisBin);
    double error = hElFR->GetBinError(thisBin);

    return value + double(sys)*error;
}

double TriLeptonBase::getFakeWeight(const vector<Muon> &muons, const vector<Electron> &electrons, int sys) {
    double weight = -1.;
    for (const auto &mu: muons) {
        if (mu.PassID(MuonIDs.at(0))) continue;
        const double fr = getMuonFakeProb(mu, sys);
        weight *= -1.*(fr / (1.-fr));
    }
    for (const auto &ele: electrons) {
        if (ele.PassID(ElectronIDs.at(0))) continue;
        const double fr = getElectronFakeProb(ele, sys);
        weight *= -1.*(fr / (1.-fr));
    }
    return weight;
}
