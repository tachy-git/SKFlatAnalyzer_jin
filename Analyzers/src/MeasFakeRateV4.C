#include "MeasFakeRateV4.h"

MeasFakeRateV4::MeasFakeRateV4() {}
MeasFakeRateV4::~MeasFakeRateV4() {
    if (RunSyst) delete hNPV_SF;
}

void MeasFakeRateV4::initializeAnalyzer() {
    // Userflags
    MeasFakeMu8 = HasFlag("MeasFakeMu8");
    MeasFakeMu17 = HasFlag("MeasFakeMu17");
    MeasFakeEl8 = HasFlag("MeasFakeEl8");
    MeasFakeEl12 = HasFlag("MeasFakeEl12");
    MeasFakeEl23 = HasFlag("MeasFakeEl23");
    MeasFakeMu = MeasFakeMu8 || MeasFakeMu17;
    MeasFakeEl = MeasFakeEl8 ||  MeasFakeEl12 || MeasFakeEl23;
    RunSyst = HasFlag("RunSyst");
    RunSystSimple = HasFlag("RunSystSimple"); // Only check the effect of the selection variations

    // binnings
    cout <<  "Setting binning" << endl;
    if (MeasFakeMu) {
        ptcorr_bins = {10., 15., 20., 30., 50., 100., 200.};
        abseta_bins = {0., 0.9, 1.6, 2.4};
    } else if (MeasFakeEl) {
        ptcorr_bins = {10., 15., 20., 25., 35., 50., 100., 200.};
        abseta_bins = {0., 0.8, 1.479, 2.5};
    } else {
        cerr << "[MeasFakeRateV4::initializeAnalyzer] No path specified by userflags" << endl;
        exit(EXIT_FAILURE);
    }

    // ID settings
    cout << "Setting IDs" << endl;
    if (DataEra == "2016preVFP") {
        MuID.SetIDs("HcToWATight", "HcToWALoose", "HcToWAVeto");
        ElID.SetIDs("HcToWATight16a", "HcToWALoose16a", "HcToWAVeto16a");
    } else if (DataEra == "2016postVFP") {
        MuID.SetIDs("HcToWATight", "HcToWALoose", "HcToWAVeto");
        ElID.SetIDs("HcToWATight16b", "HcToWALoose16b", "HcToWAVeto16b");
    } else if (DataEra == "2017") {
        MuID.SetIDs("HcToWATight", "HcToWALoose", "HcToWAVeto");
        ElID.SetIDs("HcToWATight17", "HcToWALoose17", "HcToWAVeto17");
    } else if (DataEra == "2018") {
        MuID.SetIDs("HcToWATight", "HcToWALoose", "HcToWAVeto");
        ElID.SetIDs("HcToWATight18", "HcToWALoose18", "HcToWAVeto18");
    } else {
        cerr << "[MeasFakeRateV4::initializeAnalzyer] Wrong era " << DataEra << endl;
        exit(EXIT_FAILURE);
    }

    // Trigger Settings
    cout << "Setting triggers" << endl;
    if (MeasFakeMu8) {
       isoSglLepTrig = "HLT_Mu8_TrkIsoVVL_v";
       trigSafePtCut = 10.;
    } else if (MeasFakeMu17) {
       isoSglLepTrig = "HLT_Mu17_TrkIsoVVL_v";
       trigSafePtCut = 20.;
    } else if (MeasFakeEl8) {
       isoSglLepTrig = "HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v";
       trigSafePtCut = 10.;
    } else if (MeasFakeEl12) {
       isoSglLepTrig = "HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v";
       trigSafePtCut = 15.;
    } else if (MeasFakeEl23) {
       isoSglLepTrig = "HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v";
       trigSafePtCut = 25.;
    } else {
       cerr << "[MeasFakeRateV4::initializeAnalyzer] No trigger specified by userflags" << endl;
       exit(EXIT_FAILURE);
    }

    // Jet Tagger
    vector<JetTagging::Parameters> jtps;
    jtps.emplace_back(JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Medium, JetTagging::incl, JetTagging::mujets));
    mcCorr->SetJetTaggingParameters(jtps);

    // Systematics
    weightVariations = {"Central"};
    scaleVariations = {};
    selectionVariations = {};
    
    if (RunSystSimple) {
        selectionVariations = {"MotherJetPtUp", "MotherJetPtDown", "RequireHeavyTag"};
    } else if (RunSyst) {
        if (MeasFakeMu && !IsDATA) {
            weightVariations = {"Central", 
                                "PileupReweight",
                                "L1PrefireUp", "L1PrefireDown",
                                "MuonRecoSFUp", "MuonRecoSFDown"};
        }
        if (MeasFakeEl && !IsDATA) {
            weightVariations = {"Central", 
                                "PileupReweight",
                                "L1PrefireUp", "L1PrefireDown",
                                "ElectronRecoSFUp", "ElectronRecoSFDown"};
        }
        scaleVariations = {"JetResUp", "JetResDown",
                           "JetEnUp", "JetEnDown",
                           "MuonEnUp", "MuonEnDown",
                           "ElectronEnUp", "ElectronEnDown",
                           "ElectronResUp", "ElectronResDown"};
        selectionVariations = {"MotherJetPtUp", "MotherJetPtDown", "RequireHeavyTag"};
    }

    if (IsDATA) {
        systematics = {"Central"};
        systematics.insert(systematics.end(), scaleVariations.begin(), scaleVariations.end());
        systematics.insert(systematics.end(), selectionVariations.begin(), selectionVariations.end());
    } else {
        systematics = weightVariations;
        systematics.insert(systematics.end(), scaleVariations.begin(), scaleVariations.end());
        systematics.insert(systematics.end(), selectionVariations.begin(), selectionVariations.end());
    }

    // link histograms
    if (RunSyst) {
        TString PUPath = TString(getenv("DATA_DIR")) + "/" + GetEra() + "/PileUp";
        cout << "Linking NPV scale factor histograms from " << PUPath << endl;
        TFile *fNPV = nullptr;
        if (MeasFakeMu8) {
            fNPV = TFile::Open(PUPath+"/NPV_MeasFakeMu8.root");
            hNPV_SF = (TH1D*)fNPV->Get("NPV_SF"); hNPV_SF->SetDirectory(0);
        } else if (MeasFakeMu17) {
            fNPV = TFile::Open(PUPath+"/NPV_MeasFakeMu17.root");
            hNPV_SF = (TH1D*)fNPV->Get("NPV_SF"); hNPV_SF->SetDirectory(0);
        } else if (MeasFakeEl8) {
            fNPV = TFile::Open(PUPath+"/NPV_MeasFakeEl8.root");
            hNPV_SF = (TH1D*)fNPV->Get("NPV_SF"); hNPV_SF->SetDirectory(0);
        } else if (MeasFakeEl12) {
            fNPV = TFile::Open(PUPath+"/NPV_MeasFakeEl12.root");
            hNPV_SF = (TH1D*)fNPV->Get("NPV_SF"); hNPV_SF->SetDirectory(0);
        } else if (MeasFakeEl23) {
            fNPV = TFile::Open(PUPath+"/NPV_MeasFakeEl23.root");
            hNPV_SF = (TH1D*)fNPV->Get("NPV_SF"); hNPV_SF->SetDirectory(0);
        } else {
            cerr << "[MeasFakeRateV4::initializeAnalyzer] No path specified by userflags" << endl;
            exit(EXIT_FAILURE);
        }
        fNPV->Close();

        if (!hNPV_SF) {
            cerr << "[MeasFakeRateV4::initializeAnalyzer] Failed to retrieve histogrm from files" << endl;
            exit(EXIT_FAILURE);
        }
    }
}

void MeasFakeRateV4::executeEvent() {
    rawMuons = GetAllMuons();
    rawElectrons = GetAllElectrons();
    rawJets = GetAllJets();
    truth = GetGens();
    
    NonpromptParameter param;
    // loop over IDs
    vector<TString> IDs = {"loose", "tight"};
    for (const auto &ID: IDs) {
        param.SetDefault(); // Set all the parameters as "Central", with no ID
        param.SetID(ID);
        executeEventFrom(param); // It will include weightVariations

        // Scale Variations
        for (const auto &scale: scaleVariations) {
            param.SetDefault();
            param.SetID(ID);
            param.SetScale(scale);
            executeEventFrom(param);
        }
        // Selection Variations
        for (const auto &selection: selectionVariations) {
            param.SetDefault();
            param.SetID(ID);
            param.SetSelection(selection);
            executeEventFrom(param);
        }
    }
}

void MeasFakeRateV4::executeEventFrom(NonpromptParameter &param) {
    if (!PassMETFilter()) return;
    Event ev = GetEvent();
    Particle METv = ev.GetMETVector();
    if (! ev.PassTrigger(isoSglLepTrig)) return;

    // start of object definition
    vector<Muon> allMuons = rawMuons;
    vector<Electron> allElectrons = rawElectrons;
    vector<Jet> allJets = rawJets;
    ApplyScaleVariation(allMuons, allElectrons, allJets, param.GetScale());
    
    // select objects
    TString baseMuID, baseElID;
    if (param.GetID() == "loose") {
        baseMuID = MuID.GetLooseID();
        baseElID = ElID.GetLooseID();
    } else if (param.GetID() == "tight") {
        baseMuID = MuID.GetTightID();
        baseElID = ElID.GetTightID();
    } else {
        cerr << "[MeasFakeRateV4::executeEventWith] Wrong ID " << param.GetID() << endl;
        exit(EXIT_FAILURE);
    }
    
    vector<Muon> vetoMuons = SelectMuons(allMuons, MuID.GetVetoID(), 10., 2.4);
    vector<Muon> muons = SelectMuons(allMuons, baseMuID, 10., 2.4);
    vector<Electron> vetoElectrons = SelectElectrons(allElectrons, ElID.GetVetoID(), 10., 2.5);
    vector<Electron> electrons = SelectElectrons(allElectrons, baseElID, 10., 2.5);

    // apply selection variations to jet pt cut
    const double jetPtCut = GetJetPtCut(param.GetSelection());
    vector<Jet> jets = SelectJets(allJets, "tight", jetPtCut, 2.4);
    jets = JetsVetoLeptonInside(jets, vetoElectrons, vetoMuons, 0.4);
    vector<Jet> bjets = SelectBJets(jets, JetTagging::DeepJet, JetTagging::Medium);
    // end of object definitions, scale and selection variations
    // HeavyTag will be set on SelectEvent function

    // event selection
    TString channel = "";
    if (param.GetSelection() == "RequireHeavyTag") {
        channel = SelectEvent(muons, vetoMuons, electrons, vetoElectrons, bjets);
    } else {
        channel = SelectEvent(muons, vetoMuons, electrons, vetoElectrons, jets);
    }
    if (channel == "") return;
    
    // Fill Histograms
    // Not the prefix of the histograms will be channel+"/"+param.GetName()+"/"+...
    if (param.GetScale() != "Central" || param.GetSelection() != "Central") {
        const double weight = GetEventWeight(param, ev, muons, electrons, jets);
        FillObjects(channel, param.GetName(), muons, electrons, jets, bjets, METv, weight);
    } else {
        for (const auto &syst: weightVariations) {
            param.SetSyst(syst);
            const double weight = GetEventWeight(param, ev, muons, electrons, jets);
            FillObjects(channel, param.GetName(), muons, electrons, jets, bjets, METv, weight);
        }
    }
}

void MeasFakeRateV4::ApplyScaleVariation(vector<Muon> &muons, vector<Electron> &electrons, vector<Jet> &jets, const TString &scale) {
    if (scale == "Central") {
        // do nothing
    } else if (scale == "MuonEnUp") {
        muons = ScaleMuons(muons, 1);
    } else if (scale == "MuonEnDown") {
        muons = ScaleMuons(muons, -1);
    } else if (scale == "ElectronEnUp") {
        electrons = ScaleElectrons(electrons, 1);
    } else if (scale == "ElectronEnDown") {
        electrons = ScaleElectrons(electrons, -1);
    } else if (scale == "ElectronResUp") {
        electrons = SmearElectrons(electrons, 1);
    } else if (scale == "ElectronResDown") {
        electrons = SmearElectrons(electrons, -1);
    } else if (scale == "JetEnUp") {
        jets = ScaleJets(jets, 1);
    } else if (scale == "JetEnDown") {
        jets = ScaleJets(jets, -1);
    } else if (scale == "JetResUp") {
        jets = SmearJets(jets, 1);
    } else if (scale == "JetResDown") {
        jets = SmearJets(jets, -1);
    } else {
        cerr << "[MeasFakeRateV4::executeEventWith] Wrong scale variation " << scale << endl;
        exit(EXIT_FAILURE);
    }
}

double MeasFakeRateV4::GetJetPtCut(const TString &selection) {
    if (selection == "MotherJetPtUp") {
        return 60.;
    } else if (selection == "MotherJetPtDown") {
        return 30.;
    } else {
        return 40.;
    }
}

TString MeasFakeRateV4::SelectEvent(const vector<Muon> &muons, const vector<Muon> &vetoMuons, const vector<Electron> &electrons, const vector<Electron> &vetoElectrons, const vector<Jet> &jets) {
    TString channel = "";
    // Lepton Multiplicities
    const bool SglMu = (muons.size() == 1 && vetoMuons.size() == 1 && electrons.size() == 0 && vetoElectrons.size() == 0);
    const bool SglEl = (electrons.size() == 1 && vetoElectrons.size() == 1 && muons.size() == 0 && vetoMuons.size() == 0);
    const bool DblMu = (muons.size() == 2 && vetoMuons.size() == 2 && electrons.size() == 0 && vetoElectrons.size() == 0);
    const bool DblEl = (electrons.size() == 2 && vetoElectrons.size() == 2 && muons.size() == 0 && vetoMuons.size() == 0);

    // Muon channels
    if (MeasFakeMu) {
        if (! (SglMu || DblMu)) return channel;
        if (! (muons.at(0).Pt() > trigSafePtCut)) return channel;
        if (! (jets.size() > 0)) return channel;

        // Divide channel by lepton multiplicity
        if (SglMu) {
            bool existAwayJet = false;
            for (const auto &j: jets) {
                if (j.DeltaR(muons.at(0)) > 0.7) {
                    existAwayJet = true;
                    break;
                }
            }
            if (! existAwayJet) return channel;
            channel = "Inclusive";
        } else { // DblMu
            const Particle ZCand = muons.at(0) + muons.at(1);
            const bool isOnZ = (fabs(ZCand.M() - 91.2) < 15.);
            if (! isOnZ) return channel;
            channel = "ZEnriched";
        }
    } else if (MeasFakeEl) {
        if (! (SglEl || DblEl)) return channel;
        if (! (electrons.at(0).Pt() > trigSafePtCut)) return channel;
        if (! (jets.size() > 0)) return channel;

        // Divide channel by lepton multiplicity
        if (SglEl) {
            bool existAwayJet = false;
            for (const auto &j: jets) {
                if (j.DeltaR(electrons.at(0)) > 0.7) {
                    existAwayJet = true;
                    break;
                }
            }
            if (! existAwayJet) return channel;
            channel = "Inclusive";
        } else { // DblEl
            const Particle ZCand = electrons.at(0) + electrons.at(1);
            const bool isOnZ = (fabs(ZCand.M() - 91.2) < 15.);
            if (! isOnZ) return channel;
            channel = "ZEnriched";
        }
    }
    return channel;
}

double MeasFakeRateV4::GetEventWeight(const NonpromptParameter &param, Event &ev, vector<Muon> &muons, vector<Electron> &electrons, vector<Jet> &jets) {
    double weight = 1.;
    if (!IsDATA) {
        weight = MCweight()*ev.GetTriggerLumi("Full");
        // L1 Prefire
        if (param.GetSyst() == "L1PrefireUp") {
            weight *= GetPrefireWeight(1);
        } else if (param.GetSyst() == "L1PrefireDown") {
            weight *= GetPrefireWeight(-1);
        } else {
            weight *= GetPrefireWeight(0);
        }

        // Top PT reweight
        if (MCSample.Contains("TTLL") || MCSample.Contains("TTLJ")) {
            weight *= mcCorr->GetTopPtReweight(truth);
        }

        // Only apply lepton reco scale factor
        for (const auto &mu: muons) {
            if (param.GetSyst() == "MuonRecoSFUp") {
                weight *= mcCorr->MuonReco_SF(mu.Eta(), 1);
            } else if (param.GetSyst() == "MuonRecoSFDown") {
                weight *= mcCorr->MuonReco_SF(mu.Eta(), -1);
            } else {
                weight *= mcCorr->MuonReco_SF(mu.Eta(), 0);
            }
        }
        for (const auto &el: electrons) {
            if (param.GetSyst() == "ElectronRecoSFUp") {
                weight *= mcCorr->ElectronReco_SF(el.scEta(), el.Pt(), 1);
            } else if (param.GetSyst() == "ElectronRecoSFDown") {
                weight *= mcCorr->ElectronReco_SF(el.scEta(), el.Pt(),-1);
            } else {
                weight *= mcCorr->ElectronReco_SF(el.scEta(), el.Pt(), 0);
            }
        }

        // nPV reweight
        if (param.GetSyst() == "PileupReweight") {
            weight *= GetPileUpWeight(nPileUp, 0);
        } else if (RunSyst) {
            weight *= GetNPVReweight(nPV);
        } else {
        }
        if (param.GetSelection() == "RequireHeavyTag") {
            JetTagging::Parameters jtp = JetTagging::Parameters(JetTagging::DeepJet,
                                                                JetTagging::Medium,
                                                                JetTagging::incl,
                                                                JetTagging::mujets);
            weight *= mcCorr->GetBTaggingReweight_1a(jets, jtp);
        }
    }
    return weight;
}

double MeasFakeRateV4::GetNPVReweight(unsigned int nPV) {
    nPV = max(1, min(int(nPV), 69));
    const unsigned int bin = hNPV_SF->FindBin(nPV);
    double this_sf = hNPV_SF->GetBinContent(bin);
    this_sf = fabs(this_sf-1.) < 0.5 ? this_sf : 1.;
    // for 2016a and 2016b, some bins (nPV > 60) are empty...
    // almost no event would be affected by this bins, just apply 1
    return this_sf;
}

TString MeasFakeRateV4::FindBin(const double ptcorr, const double abseta) {
    int ptcorr_idx = -1;
    int abseta_idx = -1;
    for (int i = 0; i < ptcorr_bins.size()-1; i++) {
        if (ptcorr_bins.at(i) <= ptcorr && ptcorr < ptcorr_bins.at(i+1)) {
            ptcorr_idx = i;
            break;
        }
    }
    for (int i = 0; i < abseta_bins.size()-1; i++) {
        if (abseta_bins.at(i) <= abseta && abseta < abseta_bins.at(i+1)) {
            abseta_idx = i;
            break;
        }
    }
    if (ptcorr_idx == -1)
        ptcorr_idx = ptcorr_bins.size()-2;

    // for eta bins, we will use EB1, EB2, EE
    TString etaBin = "";
    if (abseta_idx == 0) {
        etaBin = "EB1";
    } else if (abseta_idx == 1) {
        etaBin = "EB2";
    } else if (abseta_idx == 2){
        etaBin = "EE";
    } else {
        cerr << "[MeasFakeRateV4::FindBin] Wrong abseta index " << abseta_idx << endl;
        exit(EXIT_FAILURE);
    }

    TString formattedString = TString::Format("ptcorr_%dto%d_%s",
                                              static_cast<int>(ptcorr_bins.at(ptcorr_idx)),
                                              static_cast<int>(ptcorr_bins.at(ptcorr_idx+1)),
                                              etaBin.Data());
    formattedString.ReplaceAll(".", "p");
    // for debug
    // cout << "ptcorr: " << ptcorr << " abseta: " << abseta << " bin: " << formattedString << endl;

    return formattedString;
}

void MeasFakeRateV4::FillObjects(const TString &channel,
                                 const TString &prefix,
                                 const vector<Muon> &muons,
                                 const vector<Electron> &electrons,
                                 const vector<Jet> &jets,
                                 const vector<Jet> &bjets,
                                 const Particle &METv,
                                 const double &weight) {
    // Divide into MeasFakeMu and MeasFakeEl first
    if (MeasFakeMu && channel == "Inclusive") {
        const Muon &mu = muons.at(0);
        const double mT = MT(mu, METv);
        const double mTfix = TMath::Sqrt(2.*35.*METv.Pt()*(1.-TMath::Cos(mu.DeltaPhi(METv))));
        const double ptcorr = mu.Pt()*(1.+max(0., mu.MiniRelIso()-0.1));
        const double abseta = fabs(mu.Eta());
        const TString thisbin = FindBin(ptcorr, abseta);

        FillHist(channel+"/"+prefix+"/muon/pt", mu.Pt(), weight, 300, 0., 300.);
        FillHist(channel+"/"+prefix+"/muon/eta", mu.Eta(), weight, 48, -2.4, 2.4);
        FillHist(channel+"/"+prefix+"/muon/phi", mu.Phi(), weight, 64, -3.2, 3.2);
        FillHist(channel+"/"+prefix+"/muon/ptcorr", ptcorr, weight, ptcorr_bins);
        FillHist(channel+"/"+prefix+"/muon/abseta", abseta, weight, abseta_bins);
        FillHist(channel+"/"+prefix+"/MT", mT, weight, 600, 0., 300.);
        FillHist(channel+"/"+prefix+"/MTfix", mTfix, weight, 600, 0., 300.);
        FillHist(channel+"/"+prefix+"/MET", METv.Pt(), weight, 500, 0., 500.);
        FillHist(channel+"/"+prefix+"/nJets", jets.size(), weight, 10, 0., 10.);
        FillHist(channel+"/"+prefix+"/nBJets", bjets.size(), weight, 5, 0., 5.);
        FillHist(channel+"/"+prefix+"/nPV", nPV, weight, 70, 0., 70.);
        FillHist(channel+"/"+prefix+"/nPU", nPileUp, weight, 70, 0., 70.);
        FillHist(channel+"/"+prefix+"/MTvsMET", mT, METv.Pt(), weight, 500, 0., 500., 500, 0., 500.);
        FillHist(channel+"/"+prefix+"/muon/pt", mu.Pt(), weight, 300, 0., 300.);
        FillHist(thisbin+"/"+channel+"/"+prefix+"/muon/ptcorr", ptcorr, weight, 200, 0., 200.);
        FillHist(thisbin+"/"+channel+"/"+prefix+"/muon/abseta", abseta, weight, 24, 0., 2.4);
        FillHist(thisbin+"/"+channel+"/"+prefix+"/MT", mT, weight, 500, 0., 500.);
        FillHist(thisbin+"/"+channel+"/"+prefix+"/MTfix", mTfix, weight, 600, 0., 300.);
        FillHist(thisbin+"/"+channel+"/"+prefix+"/MET", METv.Pt(), weight, 600, 0., 300.);
        
        // Fill subchannel
        TString subchannel = "";
        if (mT < 25. && METv.Pt() < 25.) {
            subchannel = "QCDEnriched";
        } else if (mT > 60.) {
            subchannel = "WEnriched";
        } else {
            return;
        }
        FillHist(thisbin+"/"+subchannel+"/"+prefix+"/muon/pt", mu.Pt(), weight, 200, 0., 200.);
        FillHist(thisbin+"/"+subchannel+"/"+prefix+"/muon/eta", mu.Eta(), weight, 48, -2.4, 2.4);
        FillHist(thisbin+"/"+subchannel+"/"+prefix+"/muon/ptcorr", ptcorr, weight, 200, 0., 200.);
        FillHist(thisbin+"/"+subchannel+"/"+prefix+"/muon/abseta", abseta, weight, 24, 0., 2.4);
        FillHist(thisbin+"/"+subchannel+"/"+prefix+"/MT", mT, weight, 300, 0., 300.);
        FillHist(thisbin+"/"+subchannel+"/"+prefix+"/MET", METv.Pt(), weight, 300, 0., 300.);
    } else if (MeasFakeMu && channel == "ZEnriched") {
        const Particle ZCand = muons.at(0) + muons.at(1);
        FillHist(channel+"/"+prefix+"/ZCand/mass", ZCand.M(), weight, 40, 75., 115.);
        FillHist(channel+"/"+prefix+"/ZCand/pt", ZCand.Pt(), weight, 300, 0., 300.);
        FillHist(channel+"/"+prefix+"/ZCand/eta", ZCand.Eta(), weight, 100, -5., 5.);
        FillHist(channel+"/"+prefix+"/ZCand/phi", ZCand.Phi(), weight, 64, -3.2, 3.2);
        FillHist(channel+"/"+prefix+"/muons/1/pt", muons.at(0).Pt(), weight, 300, 0., 300.);
        FillHist(channel+"/"+prefix+"/muons/1/eta", muons.at(0).Eta(), weight, 48, -2.4, 2.4);
        FillHist(channel+"/"+prefix+"/muons/1/phi", muons.at(0).Phi(), weight, 64, -3.2, 3.2);
        FillHist(channel+"/"+prefix+"/muons/2/pt", muons.at(1).Pt(), weight, 300, 0., 300.);
        FillHist(channel+"/"+prefix+"/muons/2/eta", muons.at(1).Eta(), weight, 48, -2.4, 2.4);
        FillHist(channel+"/"+prefix+"/muons/2/phi", muons.at(1).Phi(), weight, 64, -3.2, 3.2);
        FillHist(channel+"/"+prefix+"/nJets", jets.size(), weight, 10, 0., 10.);
        FillHist(channel+"/"+prefix+"/nBJets", bjets.size(), weight, 5, 0., 5.);
        FillHist(channel+"/"+prefix+"/nPV", nPV, weight, 70, 0., 70.);
        FillHist(channel+"/"+prefix+"/nPU", nPileUp, weight, 70, 0., 70.);
    } else if (MeasFakeEl && channel == "Inclusive") {
        const Electron &el = electrons.at(0);
        const double mT = MT(el, METv);
        const double mTfix = TMath::Sqrt(2.*35.*METv.Pt()*(1.-TMath::Cos(el.DeltaPhi(METv))));
        const double ptcorr = el.Pt()*(1.+max(0., el.MiniRelIso()-0.1));
        const double abseta = fabs(el.scEta());
        const TString thisbin = FindBin(ptcorr, abseta);
        
        FillHist(channel+"/"+prefix+"/electron/pt", el.Pt(), weight, 300, 0., 300.);
        FillHist(channel+"/"+prefix+"/electron/eta", el.Eta(), weight, 50, -254, 2.5);
        FillHist(channel+"/"+prefix+"/electron/phi", el.Phi(), weight, 64, -3.2, 3.2);
        FillHist(channel+"/"+prefix+"/electron/ptcorr", ptcorr, weight, ptcorr_bins);
        FillHist(channel+"/"+prefix+"/electron/abseta", abseta, weight, abseta_bins);
        FillHist(channel+"/"+prefix+"/MT", mT, weight, 600, 0., 300.);
        FillHist(channel+"/"+prefix+"/MTfix", mTfix, weight, 600, 0., 300.);
        FillHist(channel+"/"+prefix+"/MET", METv.Pt(), weight, 500, 0., 500.);
        FillHist(channel+"/"+prefix+"/nJets", jets.size(), weight, 10, 0., 10.);
        FillHist(channel+"/"+prefix+"/nBJets", bjets.size(), weight, 5, 0., 5.);
        FillHist(channel+"/"+prefix+"/nPV", nPV, weight, 70, 0., 70.);
        FillHist(channel+"/"+prefix+"/nPU", nPileUp, weight, 70, 0., 70.);
        FillHist(channel+"/"+prefix+"/MTvsMET", mT, METv.Pt(), weight, 500, 0., 500., 500, 0., 500.);
        FillHist(thisbin+"/"+channel+"/"+prefix+"/electron/ptcorr", ptcorr, weight, ptcorr_bins);
        FillHist(thisbin+"/"+channel+"/"+prefix+"/electron/abseta", abseta, weight, abseta_bins);
        FillHist(thisbin+"/"+channel+"/"+prefix+"/MT", mT, weight, 600, 0., 300.);
        FillHist(thisbin+"/"+channel+"/"+prefix+"/MTfix", mTfix, weight, 600, 0., 300.);
        FillHist(thisbin+"/"+channel+"/"+prefix+"/MET", METv.Pt(), weight, 500, 0., 500.);

        // Fill subchannel
        TString subchannel = "";
        if (mT < 25. && METv.Pt() < 25.) {
            subchannel = "QCDEnriched";
        } else if (mT > 60.) {
            subchannel = "WEnriched";
        } else {
            return;
        }
        FillHist(thisbin+"/"+subchannel+"/"+prefix+"/electron/pt", el.Pt(), weight, 300, 0., 300.);
        FillHist(thisbin+"/"+subchannel+"/"+prefix+"/electron/eta", el.Eta(), weight, 50, -2.5, 2.5);
        FillHist(thisbin+"/"+subchannel+"/"+prefix+"/electron/ptcorr", ptcorr, weight, 300, 0., 300.);
        FillHist(thisbin+"/"+subchannel+"/"+prefix+"/electron/abseta", abseta, weight, 25, 0., 2.5);
        FillHist(thisbin+"/"+subchannel+"/"+prefix+"/MT", mT, weight, 300, 0., 300.);
        FillHist(thisbin+"/"+subchannel+"/"+prefix+"/MET", METv.Pt(), weight, 300, 0., 300.);
    } else if (MeasFakeEl && channel == "ZEnriched") {
        const Particle ZCand = electrons.at(0) + electrons.at(1);    
        FillHist(channel+"/"+prefix+"/ZCand/mass", ZCand.M(), weight, 40, 75., 115.);
        FillHist(channel+"/"+prefix+"/ZCand/pt", ZCand.Pt(), weight, 300, 0., 300.);
        FillHist(channel+"/"+prefix+"/ZCand/eta", ZCand.Eta(), weight, 100, -5., 5.);
        FillHist(channel+"/"+prefix+"/ZCand/phi", ZCand.Phi(), weight, 64, -3.2, 3.2);
        FillHist(channel+"/"+prefix+"/electrons/1/pt", electrons.at(0).Pt(), weight, 300, 0., 300.);
        FillHist(channel+"/"+prefix+"/electrons/1/eta", electrons.at(0).Eta(), weight, 50, -2.5, 2.5);
        FillHist(channel+"/"+prefix+"/electrons/1/phi", electrons.at(0).Phi(), weight, 64, -3.2, 3.2);
        FillHist(channel+"/"+prefix+"/electrons/2/pt", electrons.at(1).Pt(), weight, 300, 0., 300.);
        FillHist(channel+"/"+prefix+"/electrons/2/eta", electrons.at(1).Eta(), weight, 50, -2.5, 2.5);
        FillHist(channel+"/"+prefix+"/electrons/2/phi", electrons.at(1).Phi(), weight, 64, -3.2, 3.2);
        FillHist(channel+"/"+prefix+"/nJets", jets.size(), weight, 10, 0., 10.);
        FillHist(channel+"/"+prefix+"/nBJets", bjets.size(), weight, 5, 0., 5.);
        FillHist(channel+"/"+prefix+"/nPV", nPV, weight, 70, 0., 70.);
        FillHist(channel+"/"+prefix+"/nPU", nPileUp, weight, 70, 0., 70.);
    } else {
        cerr << "[MeasFakeRateV4::FillObjects] Wrong channel " << channel << endl;
        exit(EXIT_FAILURE);
    }
}
