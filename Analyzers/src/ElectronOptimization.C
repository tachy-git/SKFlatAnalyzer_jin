#include "ElectronOptimization.h"

ElectronOptimization::ElectronOptimization() {}
ElectronOptimization::~ElectronOptimization() {
    outfile->cd();
    Events->Write();
}

void ElectronOptimization::initializeAnalyzer() {
    // Link Tree Contents
    Events = new TTree("Events", "Events");
    Events->Branch("nElectrons", &nElectrons);
    Events->Branch("ptCorr", ptCorr, "ptCorr[nElectrons]/F");
    Events->Branch("scEta", scEta, "scEta[nElectrons]/F");
    Events->Branch("MVANoIso", MVANoIso, "MVANoIso[nElectrons]/F");
    Events->Branch("MiniRelIso", MiniRelIso, "MiniRelIso[nElectrons]/F");
    Events->Branch("DeltaR", DeltaR, "DeltaR[nElectrons]/F");
    Events->Branch("SIP3D", SIP3D, "SIP3D[nElectrons]/F");
    Events->Branch("PassMVANoIsoWP90", PassMVANoIsoWP90, "PassMVANoIsoWP90[nElectrons]/O");
    Events->Branch("PassMVANoIsoWPLoose", PassMVANoIsoWPLoose, "PassMVANoIsoWPLoose[nElectrons]/O");
    Events->Branch("NearestJetFlavour", NearestJetFlavour, "NearestJetFlavour[nElectrons]/I");
    Events->Branch("genWeight", &genWeight);

    // ID settings
    ElectronVetoID = "HcToWAVeto";

    // Jet tagger
    vector<JetTagging::Parameters> jtps;
    jtps.emplace_back(JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Medium, JetTagging::incl, JetTagging::mujets));
    mcCorr->SetJetTaggingParameters(jtps);
}

void ElectronOptimization::executeEvent() {
    if (!PassMETFilter()) return;

    // object definition
    Event ev = GetEvent();
    //vector<Muon>     rawMuons = GetAllMuons();
    vector<Electron> rawElectrons = GetAllElectrons();
    vector<Jet>      rawJets = GetAllJets();
    Particle         METv = ev.GetMETVector();
    vector<Gen>      truth = GetGens();

    //vector<Muon> vetoMuons = SelectMuons(rawMuons, "HcToWAVeto", 10., 2.4);
    //vector<Muon> tightMuons = SelectMuons(vetoMuons, "HcToWATight", 10., 2.4);
    vector<Electron> vetoElectrons = SelectElectrons(rawElectrons, ElectronVetoID, 10., 2.5);
    vector<Jet> jets = SelectJets(rawJets, "tight", 20., 2.4);
    vector<Jet> bjets;
    for (const auto &j: jets) {
        const double btagScore = j.GetTaggerResult(JetTagging::DeepJet);
        const double wp = mcCorr->GetJetTaggingCutValue(JetTagging::DeepJet, JetTagging::Medium);
        if (btagScore > wp)
            bjets.emplace_back(j);
    }

    //std::sort(vetoMuons.begin(), vetoMuons.end(), PtComparing);
    //std::sort(tightMuons.begin(), tightMuons.end(), PtComparing);
    std::sort(vetoElectrons.begin(), vetoElectrons.end(), PtComparing);
    std::sort(jets.begin(), jets.end(), PtComparing);
    std::sort(bjets.begin(), bjets.end(), PtComparing);

    // baselind event selection
	// check if nonprompt electron exists
	vector<Electron> fakeElectrons;
	for (const auto &el: vetoElectrons) {
	  if (GetLeptonType(el, truth) > 0) continue;
	  fakeElectrons.emplace_back(el);
	}
	const bool existElectron = fakeElectrons.size() > 0;
    
	// Event selection
    if (! existElectron) return;
    //if (! (tightMuons.at(0).Charge() + tightMuons.at(1).Charge() == 0)) return;
    //Particle pair = tightMuons.at(0) + tightMuons.at(1);
    //if (! (pair.M() > 12.)) return;
    if (! (jets.size() > 0)) return;
    //if (! (bjets.size() >= 1)) return;


    // Update branches
    nElectrons = fakeElectrons.size();
	genWeight = MCweight() * ev.GetTriggerLumi("Full") * GetPrefireWeight(0) * GetPileUpWeight(nPileUp, 0);
	for (unsigned int i = 0; i < nElectrons; i++) {
        const auto &el = fakeElectrons.at(i);
        // Find the nearest jet
        Jet nearest_jet = jets.at(0);
        for (const auto &j: jets) {
            if (j.DeltaR(el) > nearest_jet.DeltaR(el))
                continue;
            nearest_jet = j;
        }
        ptCorr[i]     = el.Pt() * (1.+max(0., el.MiniRelIso()-0.1));
        scEta[i]      = el.scEta();
        MVANoIso[i]   = el.MVANoIso();
        MiniRelIso[i] = el.MiniRelIso();
        DeltaR[i]     = el.DeltaR(nearest_jet);
        SIP3D[i]      = (el.IP3D() / el.IP3Derr());
        PassMVANoIsoWP90[i] = el.PassID("passMVAID_noIso_WP90");
        PassMVANoIsoWPLoose[i] = el.PassID("passMVAID_noIso_WPLoose");
        NearestJetFlavour[i] = nearest_jet.GenHFHadronMatcherFlavour();
    }
    Events->Fill();
    return;
}
