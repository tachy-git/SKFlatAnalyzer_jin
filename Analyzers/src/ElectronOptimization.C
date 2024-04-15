#include "ElectronOptimization.h"

ElectronOptimization::ElectronOptimization() {}
ElectronOptimization::~ElectronOptimization() {}

void ElectronOptimization::initializeAnalyzer() {
    // ID settings
    if (DataEra == "2016preVFP") { ElectronVetoID = "HcToWAVeto16a"; }
    else if (DataEra == "2016postVFP") { ElectronVetoID = "HcToWAVeto16b"; }
    else if (DataEra == "2017") { ElectronVetoID = "HcToWAVeto17"; }
    else if (DataEra == "2018") { ElectronVetoID = "HcToWAVeto18"; }
    else {
        cerr << "Wrong era: " << DataEra << endl;
        exit(EXIT_FAILURE);
    
    }

    // Jet tagger
    vector<JetTagging::Parameters> jtps;
    jtps.emplace_back(JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Medium, JetTagging::incl, JetTagging::mujets));
    mcCorr->SetJetTaggingParameters(jtps);
}

void ElectronOptimization::executeEvent() {
    if (!PassMETFilter()) return;

    // object definition
    Event ev = GetEvent();
    vector<Muon>     rawMuons = GetAllMuons();
    vector<Electron> rawElectrons = GetAllElectrons();
    vector<Jet>      rawJets = GetAllJets();
    Particle         METv = ev.GetMETVector();
    vector<Gen>      truth = GetGens();

    vector<Muon> vetoMuons = SelectMuons(rawMuons, "HcToWAVeto", 10., 2.4);
    vector<Muon> tightMuons = SelectMuons(vetoMuons, "HcToWATight", 10., 2.4);
    vector<Electron> vetoElectrons = SelectElectrons(rawElectrons, ElectronVetoID, 10., 2.5);
    vector<Jet> jets = SelectJets(rawJets, "tight", 20., 2.4);
    vector<Jet> bjets;
    for (const auto &j: jets) {
        const double btagScore = j.GetTaggerResult(JetTagging::DeepJet);
        const double wp = mcCorr->GetJetTaggingCutValue(JetTagging::DeepJet, JetTagging::Medium);
        if (btagScore > wp)
            bjets.emplace_back(j);
    }

    std::sort(vetoMuons.begin(), vetoMuons.end(), PtComparing);
    std::sort(tightMuons.begin(), tightMuons.end(), PtComparing);
    std::sort(vetoElectrons.begin(), vetoElectrons.end(), PtComparing);
    std::sort(jets.begin(), jets.end(), PtComparing);
    std::sort(bjets.begin(), bjets.end(), PtComparing);

    // baselind event selection
    const bool is1E2Mu = (tightMuons.size() == 2 && vetoMuons.size() == 2 && vetoElectrons.size() == 1);

    // Event selection
    if (!is1E2Mu) return;
    if (! (tightMuons.at(0).Charge() + tightMuons.at(1).Charge() == 0)) return;
    Particle pair = tightMuons.at(0) + tightMuons.at(1);
    if (! (pair.M() > 12.)) return;
    if (! (jets.size() >= 2)) return;
    if (! (bjets.size() >= 1)) return;


    // Check nearby jet flavour
    const auto ele = vetoElectrons.at(0);
    const bool isPrompt = (GetLeptonType(ele, truth) > 0);
    if (isPrompt) return;
    
    // Define eta region
    TString region = "";
    if (ele.scEta() < 0.8) region = "Barrel";
    else if (ele.scEta() < 1.479) region = "Transition";
    else if (ele.scEta() < 2.5) region = "Endcap";
    else return;

    // Define category
    TString category = "";
    // find nearest jet
    Jet nearestJet = jets.at(0);
    for (const auto &j: jets) {
        if (j.DeltaR(ele) < nearestJet.DeltaR(ele))
        nearestJet = j;
    }      
    if (ele.DeltaR(nearestJet) > 0.3) return;
    // check jet flavour
    const int jetFlavour = nearestJet.GenHFHadronMatcherFlavour();
    if (jetFlavour == 5) {
        category = "bjet";
    } else if (jetFlavour == 4) {
        category = "cjet";
    } else {
        category = "ljet";
    }

    const double weight = MCweight() * ev.GetTriggerLumi("Full") * GetPrefireWeight(0) * GetPileUpWeight(nPileUp, 0);
    const double ptCorr = ele.Pt() * (1.+max(0., ele.MiniRelIso()-0.1));
    FillHist(region+"/"+category, ptCorr, ele.MVANoIso(), ele.MiniRelIso(), weight,
                                 100, 0., 100., 
                                 100, -1., 1.,
                                 100, 0., 1.0);

}
