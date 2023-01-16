#include "DataPreprocess.h"

DataPreprocess::DataPreprocess() {
    // Link Tree Contents
    Events = new TTree("Events", "Events");
    
    // METv
    Events->Branch("METvPt", &METvPt);
    Events->Branch("METvPhi", &METvPhi);

    // muons
    Events->Branch("nMuons", &nMuons);
    Events->Branch("MuonPtColl", MuonPtColl, "MuonPtColl[nMuons]/F");
    Events->Branch("MuonEtaColl", MuonEtaColl, "MuonEtaColl[nMuons]/F");
    Events->Branch("MuonPhiColl", MuonPhiColl, "MuonPhiColl[nMuons]/F");
    Events->Branch("MuonMassColl", MuonMassColl, "MuonMassColl[nMuons]/F");
    Events->Branch("MuonChargeColl", MuonChargeColl, "MuonChargeColl[nMuons]/I");
    Events->Branch("MuonLabelColl", MuonLabelColl, "MuonLabelColl[nMuons]/O");

    // electrons
    Events->Branch("nElectrons", &nElectrons);
    Events->Branch("ElectronPtColl", ElectronPtColl, "ElectronPtColl[nElectrons]/F");
    Events->Branch("ElectronEtaColl", ElectronEtaColl, "ElectronEtaColl[nElectrons]/F");
    Events->Branch("ElectronPhiColl", ElectronPhiColl, "ElectronPhiColl[nElectrons]/F");
    Events->Branch("ElectronMassColl", ElectronMassColl, "ElectronMassColl[nElectrons]/F");
    Events->Branch("ElectronChargeColl", ElectronChargeColl, "ElectronChargeColl[nElectrons]/I");
    Events->Branch("ElectronLabelColl", ElectronLabelColl, "ElectronLabelColl[nElectrons]/O");

    // jets
    Events->Branch("nJets", &nJets);
    Events->Branch("JetPtColl", JetPtColl, "JetPtColl[nJets]/F");
    Events->Branch("JetEtaColl", JetEtaColl, "JetEtaColl[nJets]/F");
    Events->Branch("JetPhiColl", JetPhiColl, "JetPhiColl[nJets]/F");
    Events->Branch("JetMassColl", JetMassColl, "JetMassColl[nJets]/F");
    Events->Branch("JetChargeColl", JetChargeColl, "JetChargeColl[nJets]/F");
    Events->Branch("JetBtagScoreColl", JetBtagScoreColl, "JetBtagScoreColl[nJets]/F");
    Events->Branch("JetLabelColl", JetLabelColl, "JetLabelColl[nJets]/O");
}

DataPreprocess::~DataPreprocess() {
    outfile->cd();
    Events->Write();
}

void DataPreprocess::initializeAnalyzer() {
    // flags
    Skim1E2Mu = HasFlag("Skim1E2Mu");
    Skim3Mu = HasFlag("Skim3Mu");

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
        cerr << "[diLepControlRegion::initializeAnalyzer] Wrong era " << DataEra << endl;
        exit(EXIT_FAILURE);
    }

    // Jet tagger
    vector<JetTagging::Parameters> jtps;
    jtps.emplace_back(JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Medium, JetTagging::incl, JetTagging::mujets));
    mcCorr->SetJetTaggingParameters(jtps);
}

void DataPreprocess::executeEvent() {

    if (! PassMETFilter()) return;

    Event ev = GetEvent();
    vector<Muon> rawMuons = GetAllMuons();
    vector<Electron> rawElectrons = GetAllElectrons();
    vector<Jet> rawJets = GetAllJets();
    Particle METv = ev.GetMETVector();

    // object definition
    vector<Muon> vetoMuons = SelectMuons(rawMuons, MuonIDs.at(2), 10., 2.4);
    vector<Muon> looseMuons = SelectMuons(vetoMuons, MuonIDs.at(1), 10., 2.4);
    vector<Electron> vetoElectrons = SelectElectrons(rawElectrons, ElectronIDs.at(2), 10., 2.5);
    vector<Electron> looseElectrons = SelectElectrons(vetoElectrons, ElectronIDs.at(1), 10., 2.5);
    vector<Jet> jets = SelectJets(rawJets, "tight", 15., 2.4);
    // jets = JetsVetoLeptonInside(jets, vetoElectrons, vetoMuons, 0.4);
    vector<Jet> bjets;
    for (const auto &jet: jets) {
        const double this_discr = jet.GetTaggerResult(JetTagging::DeepJet);
        if (this_discr > mcCorr->GetJetTaggingCutValue(JetTagging::DeepJet, JetTagging::Medium))
            bjets.emplace_back(jet);
    }

    // sort objects
    std::sort(vetoMuons.begin(), vetoMuons.end(), PtComparing);
    std::sort(looseMuons.begin(), looseMuons.end(), PtComparing);
    std::sort(vetoElectrons.begin(), vetoElectrons.end(), PtComparing);
    std::sort(looseElectrons.begin(), looseElectrons.end(), PtComparing);
    std::sort(jets.begin(), jets.end(), PtComparing);
    std::sort(bjets.begin(), bjets.end(), PtComparing);

    // event selection
    const bool is3Mu = looseMuons.size() == 3 && looseElectrons.size() == 0 && vetoMuons.size() == 3 && vetoElectrons.size() == 0;
    const bool is1E2Mu = looseMuons.size() == 2 && looseElectrons.size() == 1 && vetoMuons.size() == 2 && vetoElectrons.size() == 1;
    if (! (is3Mu || is1E2Mu)) return;

    TString channel = "";
    // Reduced baseline for 1E2Mu channel
    // only require OS muon pair mass > 12 GeV
    // at least two jets, one bjet
    if (is1E2Mu) {
        Muon &mu1 = looseMuons.at(0);
        Muon &mu2 = looseMuons.at(1);
        if (! (mu1.Charge() + mu2.Charge() == 0)) return;
        Particle pair = mu1+mu2;
        if (! (pair.M() > 12.)) return;
        if (! (jets.size() >= 2)) return;
        if (! (bjets.size() >= 1)) return;
        channel = "SR1E2Mu";
    }
    // Reduced baseline for 3Mu channel
    // Exist OS pair and all OS pair mass > 12 GeV
    // at least two jets, one bjet
    else {      // is3Mu
        Muon &mu1 = looseMuons.at(0);
        Muon &mu2 = looseMuons.at(1);
        Muon &mu3 = looseMuons.at(2);
        if (! (fabs(mu1.Charge()+mu2.Charge()+mu3.Charge()) == 1)) return;

        // charge combination
        Particle pair1, pair2;
        if (mu1.Charge() == mu2.Charge()) {
            pair1 = mu1 + mu3;
            pair2 = mu2 + mu3;
        }
        else if (mu1.Charge() == mu3.Charge()) {
            pair1 = mu1 + mu2;
            pair2 = mu2 + mu3;
        }
        else { // mu2.Charge() == mu3.Charge()
            pair1 = mu1 + mu2;
            pair2 = mu1 + mu3;
        }
        if (! (pair1.M() > 12.)) return;
        if (! (pair2.M() > 12.)) return;
        if (! (jets.size() >= 2)) return;
        if (! (bjets.size() >= 1)) return;
        channel = "SR3Mu";
    }
    // end event selection
    
    // store informations
    METvPt = METv.Pt();
    METvPhi = METv.Phi();
    nMuons = looseMuons.size();
    for (unsigned int i = 0; i < nMuons; i++) {
        const Muon &mu = looseMuons.at(i);
        MuonPtColl[i] = mu.Pt();
        MuonEtaColl[i] = mu.Eta();
        MuonPhiColl[i] = mu.Phi();
        MuonMassColl[i] = mu.M();
        MuonChargeColl[i] = mu.Charge();
        MuonLabelColl[i] = false;
    }
    nElectrons = looseElectrons.size();
    for (unsigned int i = 0; i < nElectrons; i++) {
        const Electron &ele = looseElectrons.at(i);
        ElectronPtColl[i] = ele.Pt();
        ElectronEtaColl[i] = ele.Eta();
        ElectronPhiColl[i] = ele.Phi();
        ElectronMassColl[i] = ele.M();
        ElectronChargeColl[i] = ele.Charge();
        ElectronLabelColl[i] = false;
    }
    nJets = jets.size();
    for (unsigned int i = 0; i < nJets; i++) {
        const Jet &jet = jets.at(i);
        JetPtColl[i] = jet.Pt();
        JetEtaColl[i] = jet.Eta();
        JetPhiColl[i] = jet.Phi();
        JetMassColl[i] = jet.M();
        JetChargeColl[i] = jet.Charge();
        JetBtagScoreColl[i] = jet.GetTaggerResult(JetTagging::DeepJet);
        JetLabelColl[i] = false;
    }

    // start matching for signal sample
    if (MCSample.Contains("MHc")) {
        // Get charged Higgs mass point
        // e.g. TTToHcToWAToMuMu_MHc-130_MA-90 -> MHc-130 -> 130
        TObjArray *tokens;
        TString token;
        tokens = MCSample.Tokenize("_");
        token = ((TObjString*) tokens->At(1))->GetString();
        tokens = token.Tokenize("-");
        token = ((TObjString*) tokens->At(1))->GetString();
        int mHc = token.Atoi();

        // Get Charged Higgs decays
        vector<Gen> truth = GetGens();
        vector<Gen> chargedDecays = findChargedDecays(truth);
        if (! (chargedDecays.size() == 3)) return;

        // find A candidates
        vector<Muon> signalColl, promptColl;
        for (const auto &mu: looseMuons) {
            if (GetLeptonType(mu, truth) == 2)
                signalColl.emplace_back(mu);
            else if (GetLeptonType(mu, truth) > 0)
                promptColl.emplace_back(mu);
            else
                continue;
        }
        if (! (signalColl.size() == 2)) return;
        Particle ACand = signalColl.at(0) + signalColl.at(1);
    
        // Ajj
        if (0 < fabs(chargedDecays.at(1).PID()) && (fabs(chargedDecays.at(1).PID()) < 6)) {
            // partons should reside in acceptance
            const Gen &parton1 = chargedDecays.at(1);
            const Gen &parton2 = chargedDecays.at(2);
            if (! (parton1.Pt() > 15.)) return;
            if (! (parton2.Pt() > 15.)) return;
            if (! (fabs(parton1.Eta()) < 2.5)) return;
            if (! (fabs(parton2.Eta()) < 2.5)) return;

            // Find nearest jets
            Jet *j1 = nullptr; double dR1 = 0.3;
            Jet *j2 = nullptr; double dR2 = 0.3;
            for (auto &jet: jets) {
                if (parton1.DeltaR(jet) < dR1) {
                    j1 = &jet;
                    dR1 = parton1.DeltaR(jet);
                }
                if (parton2.DeltaR(jet) < dR2) {
                    j2 = &jet;
                    dR2 = parton2.DeltaR(jet);
                }
            }
            if (!j1 || !j2 ) return;
            if (j1 == j2) return;
            Particle WCand = *j1 + *j2;
            Particle ChargedHiggs = ACand + WCand;
            if (! (fabs(ChargedHiggs.M() - mHc) < 20.)) return;
            
            for (unsigned int i = 0; i < nMuons; i++) {
                const Muon &mu = looseMuons.at(i);
                if (GetLeptonType(mu, truth) == 2)
                    MuonLabelColl[i] = true;
            }
            for (unsigned int i = 0; i < nJets; i++) {
                Jet *j = &(jets.at(i));
                if (j == j1) JetLabelColl[i] = true;
                if (j == j2) JetLabelColl[i] = true;
            }
        }
        else if (fabs(chargedDecays.at(1).PID()) == 11 || fabs(chargedDecays.at(2).PID()) == 11) {
            if (! channel.Contains("1E2Mu")) return;
            Gen *eleGen = nullptr;
            for (auto &gen: chargedDecays) {
                if (fabs(gen.PID()) == 11) eleGen = &gen;
            }

            const Electron &ele = looseElectrons.at(0);
            if (! (GetLeptonType(ele, truth) > 0)) return;
            if (! (ele.DeltaR(*eleGen) > 0.1)) return;

            Particle WCand = ele + METv;
            Particle ChargedHiggs = WCand + ACand;
            for (unsigned int i = 0; i < nMuons; i++)
                MuonLabelColl[i] = true;
            for (unsigned int i = 0; i < nElectrons; i++)
                ElectronLabelColl[i] = true;
        }
        else if (fabs(chargedDecays.at(1).PID()) == 13 || fabs(chargedDecays.at(2).PID()) == 13) {
            if (! channel.Contains("3Mu")) return;
            if (! (promptColl.size() == 1)) return;

            Gen *muGen = nullptr;
            for (auto &gen: chargedDecays) {
                if (fabs(gen.PID()) == 13) muGen = &gen;
            }

            const Muon &mu = promptColl.at(0);
            if (! (mu.DeltaR(*muGen) > 0.1)) return;
            Particle WCand = mu + METv;
            Particle ChargedHiggs = WCand + ACand;
            for (unsigned int i = 0; i < nMuons; i++)
                MuonLabelColl[i] = true;
        }
        else {  // A tau nu
            return;
        }
    }
    // backgrounds
    Events->Fill();
    return;
}

bool DataPreprocess::isWinDecay(vector<Gen> &chargedDecay) {
    for (const auto &gen: chargedDecay) {
        if (fabs(gen.PID() == 24))
            return true;
        else 
            continue;
    }
    return false;
}

vector<Gen> DataPreprocess::replaceWtoDecays(vector<Gen> &chargedDecay, vector<Gen> &truth) {
    while (isWinDecay(chargedDecay)) {
        // find W
        Gen Wgen;
        vector<Gen> chargedDecayNonW;
        for (auto &gen: chargedDecay) {
            if (fabs(gen.PID() == 24))
                Wgen = gen;
            else
                chargedDecayNonW.emplace_back(gen);
        }

        // find W decays
        vector<Gen> Wdecays;
        for (auto &gen: truth) {
            if (gen.MotherIndex() == Wgen.Index())
                Wdecays.emplace_back(gen);
        }

        // update charged decay
        chargedDecay.clear();
        for (auto &gen: chargedDecayNonW)
            chargedDecay.emplace_back(gen);
        for (auto &gen: Wdecays)
            chargedDecay.emplace_back(gen);
    }
    
    // remove photons
    vector<Gen> out;
    for (auto &gen: chargedDecay) {
        if (gen.PID() == 22)
            continue;
        else
            out.emplace_back(gen);
    }

    return out;
}

vector<Gen> DataPreprocess::findChargedDecays(vector<Gen> &truth) {
    // 1. Find the charged Higgs and its decays
    vector<int> chargedIdx;
    for (auto &gen: truth)
        if (fabs(gen.PID()) == 37) 
            chargedIdx.emplace_back(gen.Index());

    // find decays
    vector<Gen> chargedDecay;
    for (auto &gen: truth) {
        // H+ itself
        if (fabs(gen.PID()) == 37) continue;

        // check isDecay from H+
        int motherIdx = gen.MotherIndex();
        bool isChargedDecay = false;
        for (auto cIdx: chargedIdx)
            if (motherIdx == cIdx) isChargedDecay = true;
        if (! isChargedDecay) continue;

        chargedDecay.emplace_back(gen);
    }

    // replace W with W decays
    chargedDecay = replaceWtoDecays(chargedDecay, truth);

    return chargedDecay;
}
