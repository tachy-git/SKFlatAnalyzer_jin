from ROOT import gSystem
from ROOT import DiLeptonBase
from ROOT import TString
from ROOT.std import vector
from ROOT.JetTagging import Parameters as jParameters
from ROOT import Muon, Electron, Jet
from array import array
gSystem.Load("/cvmfs/cms.cern.ch/slc7_amd64_gcc900/external/lhapdf/6.2.3/lib/libLHAPDF.so")


class MeasFakeRate(DiLeptonBase):
    def __init__(self):
        super().__init__()
        # at this point, DiLeptonBase::initializeAnalyzer has not been called
        
    def initializePyAnalyzer(self):
        super().initializeAnalyzer()
        
        # link flags
        self.channel = None
        self.run_syst = False
        
        if super().MeasNormEl: self.channel = "MeasNormEl"
        if super().MeasNPVEl:  self.channel = "MeasNPVEl"
        if super().MeasFakeEl: self.channel = "MeasFakeEl"
        if super().MeasNormMu: self.channel = "MeasNormMu"
        if super().MeasNPVMu:  self.channel = "MeasNPVMu"
        if super().MeasFakeMu: self.channel = "MeasFakeMu" 
        if super().RunSyst:  self.run_syst = True
        
        # triggers and binning definition
        self.triggerList = vector[TString]()
        if "El" in self.channel:
            # ptcorr < 20: Ele8, ptcorr < 35: Ele12, ptcorr > 35: Ele23
            self.triggerList.emplace_back("HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v")
            self.triggerList.emplace_back("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v")
            self.triggerList.emplace_back("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v")
            self.ptcorr_bins = array('d', [10., 15., 20., 25., 35., 50., 70., 100.])
            self.eta_bins = array('d', [0., 0.8, 1.579, 2.5])
        if "Mu" in self.channel:
            # ptcorr < 30: Mu8, ptcoor > 30: Mu8
            self.triggerList.emplace_back("HLT_Mu8_TrkIsoVVL_v")
            self.triggerList.emplace_back("HLT_Mu17_TrkIsoVVL_v")
            self.ptcorr_bins = array('d', [10., 15., 20., 30., 50., 70., 100.])
            self.eta_bins = array('d', [0., 0.9, 1.6, 2.4])
       
        # systematics
        self.weightVariations = ["Central"]
        self.scaleVariations = []
        self.selectionVariations = []
        if self.run_syst:
            if self.channel == "MeasNormEl":
                self.weightVariations += ["PileupReweight",
                                          "L1PrefireUp", "L1PrefireDown",
                                          "ElectronRecoSFUp", "ElectronRecoSFDown",
                                          "HeavyTagUpUnCorr", "HeavyTagDownUnCorr",
                                          "LightTagUpUnCorr", "LightTagDownUnCorr"]
            if self.channel == "MeasNormMu": 
                self.weightVariations += ["PileupReweight",
                                          "L1PrefireUp", "L1PrefireDown",
                                          "MuonRecoSFUp", "MuonRecoSFDown",
                                          "HeavyTagUpUnCorr", "HeavyTagDownUnCorr",
                                          "LightTagUpUnCorr", "LightTagDownUnCorr"]
            self.scaleVariations += ["JetResUp", "JetResDown",
                                     "JetEnUp", "JetEnDown",
                                     "ElectronResUp", "ElectronResDown",
                                     "ElectronEnUp", "ElectronEnDown",
                                     "MuonEnUp", "MuonEnDown"]
            if self.channel in ["MeasFakeMu", "MeasFakeEl"]:
                self.selectionVariations += ["RequireHeavyTag", "MotherJetPtUp", "MotherJetPtDown"]
        self.systematics = self.weightVariations + self.scaleVariations + self.selectionVariations
    
    def executeEvent(self):
        if not super().PassMETFilter(): return None
        ev = super().GetEvent()
        rawMuons = super().GetAllMuons()
        rawElectrons = super().GetAllElectrons()
        rawJets = super().GetAllJets()
        METv = ev.GetMETVector()
        truth = super().GetGens() if not super().IsDATA else None
        
        # Central scale
        vetoMuons, looseMuons, tightMuons, vetoElectrons, looseElectrons, tightElectrons, jets, bjets = self.defineObjects(rawMuons, rawElectrons, rawJets)
        channel = self.selectEvent(ev, vetoMuons, looseMuons, vetoElectrons, looseElectrons, jets, bjets, METv)
        if not channel is None:
            objects = {"looseMuons": looseMuons,
                       "looseElectrons": looseElectrons,
                       "tightMuons": tightMuons,
                       "tightElectrons": tightElectrons,
                       "jets": jets,
                       "bjets": bjets,
                       "METv": METv
                       }
            for syst in self.weightVariations:
                weight = self.getWeight(ev, looseMuons, looseElectrons, jets, truth, syst)
                self.fillObjects(ev, objects, weight, syst)
                
        # Scale Variations
        for syst in self.scaleVariations:
            vetoMuons, looseMuons, tightMuons, vetoElectrons, looseElectrons, tightElectrons, jets, bjets = self.defineObjects(rawMuons, rawElectrons, rawJets, syst)
            channel = self.selectEvent(ev, vetoMuons, looseMuons, vetoElectrons, looseElectrons, jets, bjets, METv)
            if channel is None: continue
            objects = {"looseMuons": looseMuons,
                       "looseElectrons": looseElectrons,
                       "tightMuons": tightMuons,
                       "tightElectrons": tightElectrons,
                       "jets": jets,
                       "bjets": bjets,
                       "METv": METv
                       }
            weight = self.getWeight(ev, looseMuons, looseElectrons, jets, truth, syst)
            self.fillObjects(ev, objects, weight, syst)
        
        # Selection Variations
        for syst in self.selectionVariations:
            vetoMuons, looseMuons, tightMuons, vetoElectrons, looseElectrons, tightElectrons, jets, bjets = self.defineObjects(rawMuons, rawElectrons, rawJets, syst)
            channel = self.selectEvent(ev, vetoMuons, looseMuons, vetoElectrons, looseElectrons, jets, bjets, METv, syst)
            if channel is None: continue
            objects = {"looseMuons": looseMuons,
                       "looseElectrons": looseElectrons,
                       "tightMuons": tightMuons,
                       "tightElectrons": tightElectrons,
                       "jets": jets,
                       "bjets": bjets,
                       "METv": METv
                       }
            weight = self.getWeight(ev, looseMuons, looseElectrons, jets, truth, syst)
            self.fillObjects(ev, objects, weight, syst)
    
    def defineObjects(self, rawMuons, rawElectrons, rawJets, syst="Central"):
        # first copy objects
        allMuons = rawMuons
        allElectrons = rawElectrons
        allJets = rawJets
        
        # check the syst argument
        if not syst in self.systematics:
            print(f"[PromptEstimator::defineObjects] Wrong systematics {syst}")
        if syst == "MuonEnUp":         allMuons = super().ScaleMuons(allMuons, +1)
        if syst == "MuonEnDown":       allMuons = super().ScaleMuons(allMuons, -1)
        if syst == "ElectronResUp":    allElectrons = super().SmearElectrons(allElectrons, +1)
        if syst == "ElectronsResDown": allElectrons = super().SmearElectrons(allElectrons, -1)
        if syst == "ElectronEnUp":     allElectrons = super().ScaleElectrons(allElectrons, +1)
        if syst == "ElectronEnDown":   allElectrons = super().ScaleElectrons(allElectrons, -1)
        if syst == "JetResUp":         allJets = super().SmearJets(allJets, +1)
        if syst == "JetResDown":       allJets = super().SmearJets(allJets, -1)
        if syst == "JetEnUp":          allJets = super().ScaleJets(allJets, +1)
        if syst == "JetEnDown":        allJets = super().ScaleJets(allJets, -1)
        
        vetoMuons = super().SelectMuons(allMuons, super().MuonIDs[2], 10., 2.4)
        looseMuons = super().SelectMuons(vetoMuons, super().MuonIDs[1], 10., 2.4)
        tightMuons = super().SelectMuons(looseMuons, super().MuonIDs[0], 10., 2.4)
        vetoElectrons = super().SelectElectrons(allElectrons, super().ElectronIDs[2], 10., 2.5)
        looseElectrons = super().SelectElectrons(vetoElectrons, super().ElectronIDs[1], 10., 2.5)
        tightElectrons = super().SelectElectrons(looseElectrons, super().ElectronIDs[0], 10., 2.5)

        if syst == "MotherJetPtUp":     jetPtCut = 60.
        elif syst == "MotherJetPtDown": jetPtCut = 30.
        else:                           jetPtCut = 40.
        jets = super().SelectJets(allJets, "tight", jetPtCut, 2.4)
        jets = super().JetsVetoLeptonInside(jets, vetoElectrons, vetoMuons, 0.4)
        bjets = vector[Jet]()
        for jet in jets:
            btagScore = jet.GetTaggerResult(3)                  # DeepJet score
            wp = super().mcCorr.GetJetTaggingCutValue(3, 1)     # DeepJet Medium
            if btagScore > wp: bjets.emplace_back(jet)
            
        vetoMuons = vector[Muon](sorted(vetoMuons, key=lambda x: x.Pt(), reverse=True))
        tightMuons = vector[Muon](sorted(tightMuons, key=lambda x: x.Pt(), reverse=True))
        vetoElectrons = vector[Electron](sorted(vetoElectrons, key=lambda x: x.Pt(), reverse=True))
        tightElectrons = vector[Electron](sorted(tightElectrons, key=lambda x: x.Pt(), reverse=True))
        jets = vector[Jet](sorted(jets, key=lambda x: x.Pt(), reverse=True))
        bjets = vector[Jet](sorted(bjets, key=lambda x: x.Pt(), reverse=True))
        
        return (vetoMuons, looseMuons, tightMuons, vetoElectrons, looseElectrons, tightElectrons, jets, bjets)
    
    def selectEvent(self, event, vetoMuons, looseMuons, vetoElectrons, looseElectrons, jets, bjets, METv, syst="Central"):
        NML, NMV = len(looseMuons), len(vetoMuons)
        NEL, NEV = len(looseElectrons), len(vetoElectrons)
        if self.channel == "MeasNormEl":
            if not (NEL == 2 and NEV == 2 and NML == 0 and NMV == 0): return None
            if not event.PassTrigger(self.triggerList): return None
            ZCand = looseElectrons.at(0) + looseElectrons.at(1)
            if not (60. < ZCand.M() and ZCand.M() < 120.): return None
            return self.channel
        if self.channel == "MeasFakeEl":
            if not (NEL == 1 and NEV == 1 and NML == 0 and NMV == 0): return None
            if not event.PassTrigger(self.triggerList): return None
            if syst == "RequireHeavyTag":
                if not bjets.size() >= 1: return None
                if not bjets.at(0).DeltaR(looseElectrons.at(0)) > 1.: return None
            else:
                if not jets.size() >= 1: return None
                if not jets.at(0).DeltaR(looseElectrons.at(0)) > 1.: return None
            if not METv.Pt() < 25.: return None
            if not (METv+looseElectrons.at(0)).Mt() < 25.: return None
            return self.channel
        if self.channel == "MeasNPVEl":
            if not (NEL == 1 and NEV == 1 and NML == 0 and NMV == 0): return None
            if not event.PassTrigger(self.triggerList): return None
            if not jets.size() >= 1: return None
            if not jets.at(0).DeltaR(looseElectrons.at(0)) > 1.: return None
            #if not METv.Pt() < 25.: return None
            #if not (METv+looseElectrons.at(0)).Mt() < 25.: return None
            return self.channel
        if self.channel == "MeasNormMu":
            if not (NEL == 0 and NEV == 0 and NML == 2 and NMV == 2): return None
            if not event.PassTrigger(self.triggerList): return None
            ZCand = looseMuons.at(0) + looseMuons.at(1)
            if not (abs(ZCand.M()-91.2) < 15.): return None 
            return self.channel
        if self.channel == "MeasFakeMu":
            if not (NEL == 0 and NEV == 0 and NML == 1 and NMV == 1): return None
            if not event.PassTrigger(self.triggerList): return None
            if syst == "RequireHeavyTag":
                if not bjets.size() >= 1: return None
                if not bjets.at(0).DeltaR(looseMuons.at(0)) > 1.: return None
            else:
                if not jets.size() >= 1: return None
                if not jets.at(0).DeltaR(looseMuons.at(0)) > 1.: return None
            if not METv.Pt() < 25.: return None
            if not (METv+looseMuons.at(0)).Mt() < 25.: return None 
            return self.channel        
        if self.channel == "MeasNPVMu":
            if not (NEL == 0 and NEV == 0 and NML == 1 and NMV == 1): return None
            if not event.PassTrigger(self.triggerList): return None
            if not jets.size() >= 1: return None
            if not jets.at(0).DeltaR(looseMuons.at(0)) > 1.: return None
            #if not METv.Pt() < 25.: return None
            #if not (METv+looseMuons.at(0)).Mt() < 25.: return None
            return self.channel
        return None
    
    def getWeight(self, event, muons, electrons, jets, truth, syst="Central"):
        weight = 1.
        if not syst in self.systematics:
            print(f"[FakeMeasurement::getWeight] Wrong systematic {syst}")
            exit(1)

        if not super().IsDATA:
            weight *= super().MCweight()
            weight *= event.GetTriggerLumi("Full")
            if syst == "L1PrefireUp":     w_prefire = super().GetPrefireWeight(1)
            elif syst == "L1PrefireDown": w_prefire = super().GetPrefireWeight(-1)
            else:                         w_prefire = super().GetPrefireWeight(0)
            
            w_topptweight = 1.
            if "TTLL" in super().MCSample or "TTLJ" in super().MCSample:
                w_topptweight = super().mcCorr.GetTopPtReweight(truth)
            weight *= w_prefire
            weight *= w_topptweight
            if "NPV" in self.channel:
                return weight

            # Do not apply lepton scale factors
            # it can be absorbed in the prompt normalization scale
            w_muonRecoSF = 1.
            #w_muonIDSF = 1.
            w_eleRecoSF = 1.
            #w_eleIDSF = 1.
            for mu in muons:
                if syst == "MuonRecoSFUp":     w_muonRecoSF *= super().getMuonRecoSF(mu, 1)
                elif syst == "MuonRecoSFDown": w_muonRecoSF *= super().getMuonRecoSF(mu, -1)
                else:                          w_muonRecoSF *= super().getMuonRecoSF(mu, 0)
                #if syst == "MuonIDSFUp":     w_muonIDSF *= self.getMuonIDSF(mu, 1)
                #elif syst == "MuonIDSFDown": w_muonIDSF *= self.getMuonIDSF(mu, -1)
                #else:                        w_muonIDSF *= self.getMuonIDSF(mu, 0)

            for ele in electrons:
                if syst == "ElectronRecoSFUp":   w_eleRecoSF *= super().mcCorr.ElectronReco_SF(ele.scEta(), ele.Pt(), 1)
                elif syst == "ElectronRecoSFUp": w_eleRecoSF *= super().mcCorr.ElectronReco_SF(ele.scEta(), ele.Pt(), -1)
                else:                            w_eleRecoSF *= super().mcCorr.ElectronReco_SF(ele.scEta(), ele.Pt(), 0)
                #if syst == "EleIDSFUp":       w_eleIDSF *= self.getEleIDSF(ele, 1)
                #elif syst == "EleIDSFDown":   w_eleIDSF *= self.getEleIDSF(ele, -1)
                #else:                         w_eleIDSF *= self.getEleIDSF(ele, 0)
            
            weight *= w_muonRecoSF
            #weight *= w_muonIDSF
            weight *= w_eleRecoSF
            #weight *= w_eleIDSF
            
            jtp = jParameters(3, 1, 0, 1)    # DeepJet, Medium, incl, mujets
            vjets = vector[Jet]()
            for j in jets: vjets.emplace_back(j)
            if syst == "HeavyTagUpUnCorr":     w_btag = super().mcCorr.GetBTaggingReweight_1a(vjets, jtp, "SystUpHTag")
            elif syst == "HeavyTagDownUnCorr": w_btag = super().mcCorr.GetBTaggingReweight_1a(vjets, jtp, "SystDownHTag")
            elif syst == "HeavyTagUpCorr":     w_btag = super().mcCorr.GetBTaggingReweight_1a(vjets, jtp, "SystUpHTagCorr")
            elif syst == "HeavyTagDownCorr":   w_btag = super().mcCorr.GetBTaggingReweight_1a(vjets, jtp, "SystDownHTagCorr")
            elif syst == "LightTagUpUnCorr":   w_btag = super().mcCorr.GetBTaggingReweight_1a(vjets, jtp, "SystUpLTag")
            elif syst == "LightTagDownUnCorr": w_btag = super().mcCorr.GetBTaggingReweight_1a(vjets, jtp, "SystDownLTag")
            elif syst == "LightTagUpCorr":     w_btag = super().mcCorr.GetBTaggingReweight_1a(vjets, jtp, "SystUpLTagCorr")
            elif syst == "LightTagDownCorr":   w_btag = super().mcCorr.GetBTaggingReweight_1a(vjets, jtp, "SystDownLTagCorr")
            else:                              w_btag = super().mcCorr.GetBTaggingReweight_1a(vjets, jtp)
            weight *= w_btag
        
        return weight
    
    def fillObjects(self, event, objects, weight, syst):
        looseMuons = objects["looseMuons"]
        tightMuons = objects["tightMuons"]
        looseElectrons = objects["looseElectrons"]
        tightElectrons = objects["tightElectrons"]
        jets = objects["jets"]
        bjets = objects["bjets"]
        METv = objects["METv"]
        
        if self.channel == "MeasNormMu":
            pair = looseMuons.at(0) + looseMuons.at(1)
            # loose
            if event.PassTrigger("HLT_Mu8_TrkIsoVVL_v"):
                trigSuffix = "Mu8"
                if syst == "PileupReweight": npvweight = super().GetPileUpWeight(super().nPileUp, 0)
                else:                        npvweight = self.getNPVReweight(super().nPV, trigSuffix)
                for idx, mu in enumerate(looseMuons, start=1):
                    super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/muons/{idx}/pt", mu.Pt(), weight*npvweight, 300, 0., 300.)
                    super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/muons/{idx}/eta", mu.Eta(), weight*npvweight, 48, -2.4, 2.4)
                    super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/muons/{idx}/phi", mu.Phi(), weight*npvweight, 64, -3.2, 3.2)
                    super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/muons/{idx}/mass", mu.M(), weight*npvweight, 10, 0., 1.)
                for idx, jet in enumerate(jets, start=1):
                    super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/jets/{idx}/pt", jet.Pt(), weight*npvweight, 300, 0., 300.)
                    super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/jets/{idx}/eta", jet.Eta(), weight*npvweight, 48, -2.4, 2.4)
                    super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/jets/{idx}/phi", jet.Phi(), weight*npvweight, 64, -3.2, 3.2)
                    super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/jets/{idx}/mass", jet.M(), weight*npvweight, 100, 0., 100.)
                for idx, bjet in enumerate(bjets, start=1):
                    super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/bjets/{idx}/pt", bjet.Pt(), weight*npvweight, 300, 0., 300.)
                    super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/bjets/{idx}/eta", bjet.Eta(), weight*npvweight, 48, -2.4, 2.4)
                    super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/bjets/{idx}/phi", bjet.Phi(), weight*npvweight, 64, -3.2, 3.2)
                    super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/bjets/{idx}/mass", bjet.M(), weight*npvweight, 100, 0., 100.)
                super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/jets/size", jets.size(), weight*npvweight, 20, 0., 20.)
                super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/bjets/size", bjets.size(), weight*npvweight, 15, 0., 15.)
                super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/METv/pt", METv.Pt(), weight*npvweight, 300, 0., 300.)
                super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/pair/pt", pair.Pt(), weight*npvweight, 300, 0., 300.)
                super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/pair/eta", pair.Eta(), weight*npvweight, 100, -5., 5.)
                super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/pair/phi", pair.Phi(), weight*npvweight, 64, -3.2, 3.2)
                super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/pair/mass", pair.M(), weight*npvweight, 80, 50., 130.)
                super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/nPV", super().nPV, weight*npvweight, 100, 0., 100.)
                super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/nPileUp", super().nPileUp, weight*npvweight, 100, 0., 100.) 
            if event.PassTrigger("HLT_Mu17_TrkIsoVVL_v") and looseMuons.at(0).Pt() > 20.: 
                trigSuffix = "Mu17"
                if syst == "PileupReweight": npvweight = super().GetPileUpWeight(super().nPileUp, 0)
                else:                        npvweight = self.getNPVReweight(super().nPV, trigSuffix)
                for idx, mu in enumerate(looseMuons, start=1):
                    super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/muons/{idx}/pt", mu.Pt(), weight*npvweight, 300, 0., 300.)
                    super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/muons/{idx}/eta", mu.Eta(), weight*npvweight, 48, -2.4, 2.4)
                    super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/muons/{idx}/phi", mu.Phi(), weight*npvweight, 64, -3.2, 3.2)
                    super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/muons/{idx}/mass", mu.M(), weight*npvweight, 10, 0., 1.)
                for idx, jet in enumerate(jets, start=1):
                    super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/jets/{idx}/pt", jet.Pt(), weight*npvweight, 300, 0., 300.)
                    super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/jets/{idx}/eta", jet.Eta(), weight*npvweight, 48, -2.4, 2.4)
                    super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/jets/{idx}/phi", jet.Phi(), weight*npvweight, 64, -3.2, 3.2)
                    super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/jets/{idx}/mass", jet.M(), weight*npvweight, 100, 0., 100.)
                for idx, bjet in enumerate(bjets, start=1):
                    super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/bjets/{idx}/pt", bjet.Pt(), weight*npvweight, 300, 0., 300.)
                    super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/bjets/{idx}/eta", bjet.Eta(), weight*npvweight, 48, -2.4, 2.4)
                    super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/bjets/{idx}/phi", bjet.Phi(), weight*npvweight, 64, -3.2, 3.2)
                    super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/bjets/{idx}/mass", bjet.M(), weight*npvweight, 100, 0., 100.)
                super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/jets/size", jets.size(), weight*npvweight, 20, 0., 20.)
                super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/bjets/size", bjets.size(), weight*npvweight, 15, 0., 15.)
                super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/METv/pt", METv.Pt(), weight*npvweight, 300, 0., 300.)
                super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/pair/pt", pair.Pt(), weight*npvweight, 300, 0., 300.)
                super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/pair/eta", pair.Eta(), weight*npvweight, 100, -5., 5.)
                super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/pair/phi", pair.Phi(), weight*npvweight, 64, -3.2, 3.2)
                super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/pair/mass", pair.M(), weight*npvweight, 80, 50., 130.)
                super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/nPV", super().nPV, weight*npvweight, 100, 0., 100.)
            # tight
            if tightMuons.size() != 2: return None
            if event.PassTrigger("HLT_Mu8_TrkIsoVVL_v"):
                trigSuffix = "Mu8"
                if syst == "PileupReweight": npvweight = super().GetPileUpWeight(super().nPileUp, 0)
                else:                        npvweight = self.getNPVReweight(super().nPV, trigSuffix)
                for idx, mu in enumerate(looseMuons, start=1):
                    super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/muons/{idx}/pt", mu.Pt(), weight*npvweight, 300, 0., 300.)
                    super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/muons/{idx}/eta", mu.Eta(), weight*npvweight, 48, -2.4, 2.4)
                    super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/muons/{idx}/phi", mu.Phi(), weight*npvweight, 64, -3.2, 3.2)
                    super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/muons/{idx}/mass", mu.M(), weight*npvweight, 10, 0., 1.)
                for idx, jet in enumerate(jets, start=1):
                    super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/jets/{idx}/pt", jet.Pt(), weight*npvweight, 300, 0., 300.)
                    super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/jets/{idx}/eta", jet.Eta(), weight*npvweight, 48, -2.4, 2.4)
                    super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/jets/{idx}/phi", jet.Phi(), weight*npvweight, 64, -3.2, 3.2)
                    super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/jets/{idx}/mass", jet.M(), weight*npvweight, 100, 0., 100.)
                for idx, bjet in enumerate(bjets, start=1):
                    super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/bjets/{idx}/pt", bjet.Pt(), weight*npvweight, 300, 0., 300.)
                    super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/bjets/{idx}/eta", bjet.Eta(), weight*npvweight, 48, -2.4, 2.4)
                    super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/bjets/{idx}/phi", bjet.Phi(), weight*npvweight, 64, -3.2, 3.2)
                    super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/bjets/{idx}/mass", bjet.M(), weight*npvweight, 100, 0., 100.)
                super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/jets/size", jets.size(), weight*npvweight, 20, 0., 20.)
                super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/bjets/size", bjets.size(), weight*npvweight, 15, 0., 15.)
                super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/METv/pt", METv.Pt(), weight*npvweight, 300, 0., 300.)
                super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/pair/pt", pair.Pt(), weight*npvweight, 300, 0., 300.)
                super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/pair/eta", pair.Eta(), weight*npvweight, 100, -5., 5.)
                super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/pair/phi", pair.Phi(), weight*npvweight, 64, -3.2, 3.2)
                super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/pair/mass", pair.M(), weight*npvweight, 80, 50., 130.)
                super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/nPV", super().nPV, weight*npvweight, 100, 0., 100.)
            if event.PassTrigger("HLT_Mu17_TrkIsoVVL_v") and tightMuons.at(0).Pt() > 20.: 
                trigSuffix = "Mu17"
                if syst == "PileupReweight": npvweight = super().GetPileUpWeight(super().nPileUp, 0)
                else:                        npvweight = self.getNPVReweight(super().nPV, trigSuffix)
                for idx, mu in enumerate(tightMuons, start=1):
                    super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/muons/{idx}/pt", mu.Pt(), weight*npvweight, 300, 0., 300.)
    }
                    super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/muons/{idx}/eta", mu.Eta(), weight*npvweight, 48, -2.4, 2.4)
                    super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/muons/{idx}/phi", mu.Phi(), weight*npvweight, 64, -3.2, 3.2)
                    super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/muons/{idx}/mass", mu.M(), weight*npvweight, 10, 0., 1.)
                for idx, jet in enumerate(jets, start=1):
                    super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/jets/{idx}/pt", jet.Pt(), weight*npvweight, 300, 0., 300.)
                    super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/jets/{idx}/eta", jet.Eta(), weight*npvweight, 48, -2.4, 2.4)
                    super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/jets/{idx}/phi", jet.Phi(), weight*npvweight, 64, -3.2, 3.2)
                    super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/jets/{idx}/mass", jet.M(), weight*npvweight, 100, 0., 100.)
                for idx, bjet in enumerate(bjets, start=1):
                    super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/bjets/{idx}/pt", bjet.Pt(), weight*npvweight, 300, 0., 300.)
                    super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/bjets/{idx}/eta", bjet.Eta(), weight*npvweight, 48, -2.4, 2.4)
                    super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/bjets/{idx}/phi", bjet.Phi(), weight*npvweight, 64, -3.2, 3.2)
                    super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/bjets/{idx}/mass", bjet.M(), weight*npvweight, 100, 0., 100.)
                super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/jets/size", jets.size(), weight*npvweight, 20, 0., 20.)
                super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/bjets/size", bjets.size(), weight*npvweight, 15, 0., 15.)
                super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/METv/pt", METv.Pt(), weight*npvweight, 300, 0., 300.)
                super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/pair/pt", pair.Pt(), weight*npvweight, 300, 0., 300.)
                super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/pair/eta", pair.Eta(), weight*npvweight, 100, -5., 5.)
                super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/pair/phi", pair.Phi(), weight*npvweight, 64, -3.2, 3.2)
                super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/pair/mass", pair.M(), weight*npvweight, 80, 50., 130.)
                super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/nPV", super().nPV, weight*npvweight, 100, 0., 100.)
            return None
        if self.channel == "MeasNPVMu":
            if event.PassTrigger("HLT_Mu8_TrkIsoVVL_v"):
                trigSuffix = "Mu8"
                super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/nPV", super().nPV, weight, 100, 0., 100.)
                super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/Rho", super().Rho, weight, 100, 0., 100.)
            if event.PassTrigger("HLT_Mu17_TrkIsoVVL_v") and looseMuons.at(0).Pt() > 20.:
                trigSuffix = "Mu17"
                super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/nPV", super().nPV, weight, 100, 0., 100.)
                super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/Rho", super().Rho, weight, 100, 0., 100.)
            
            if tightMuons.size() != 1:    return None
            if event.PassTrigger("HLT_Mu8_TrkIsoVVL_v"):
                trigSuffix = "Mu8"
                super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/nPV", super().nPV, weight, 100, 0., 100.)
                super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/Rho", super().Rho, weight, 100, 0., 100.)
            if event.PassTrigger("HLT_Mu17_TrkIsoVVL_v") and looseMuons.at(0).Pt() > 20.:
                trigSuffix = "Mu17"
                super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/nPV", super().nPV, weight, 100, 0., 100.)
                super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/Rho", super().Rho, weight, 100, 0., 100.) 
            return None
        # later
        if self.channel == "MeasFakeMu":
            mu = looseMuons.at(0)
            ptCorr = mu.Pt()*(1.+max(0., mu.MiniRelIso()-0.1))
            abseta = abs(mu.Eta())
            
            # trigger selection
            measLowPtCorr = ptCorr < 30. and event.PassTrigger("HLT_Mu8_TrkIsoVVL_v")
            measHighPtCorr = ptCorr > 30. and event.PassTrigger("HLT_Mu17_TrkIsoVVL_v")
            if not (measLowPtCorr or measHighPtCorr): return None
        
        
        if self.channel == "MeasNormEl":
            pair = looseElectrons.at(0) + looseElectrons.at(1)
            # loose
            if event.PassTrigger("HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v"):
                trigSuffix = "Ele8"
                if syst == "PileupReweight": npvweight = super().GetPileUpWeight(super().nPileUp, 0)
                else:                        npvweight = self.getNPVReweight(super().nPV, trigSuffix)
                for idx, mu in enumerate(looseElectrons, start=1):
                    super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/electrons/{idx}/pt", mu.Pt(), weight*npvweight, 300, 0., 300.)
                    super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/electrons/{idx}/eta", mu.Eta(), weight*npvweight, 48, -2.4, 2.4)
                    super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/electrons/{idx}/phi", mu.Phi(), weight*npvweight, 64, -3.2, 3.2)
                    super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/electrons/{idx}/mass", mu.M(), weight*npvweight, 10, 0., 1.)
                for idx, jet in enumerate(jets, start=1):
                    super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/jets/{idx}/pt", jet.Pt(), weight*npvweight, 300, 0., 300.)
                    super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/jets/{idx}/eta", jet.Eta(), weight*npvweight, 48, -2.4, 2.4)
                    super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/jets/{idx}/phi", jet.Phi(), weight*npvweight, 64, -3.2, 3.2)
                    super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/jets/{idx}/mass", jet.M(), weight*npvweight, 100, 0., 100.)
                for idx, bjet in enumerate(bjets, start=1):
                    super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/bjets/{idx}/pt", bjet.Pt(), weight*npvweight, 300, 0., 300.)
                    super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/bjets/{idx}/eta", bjet.Eta(), weight*npvweight, 48, -2.4, 2.4)
                    super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/bjets/{idx}/phi", bjet.Phi(), weight*npvweight, 64, -3.2, 3.2)
                    super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/bjets/{idx}/mass", bjet.M(), weight*npvweight, 100, 0., 100.)
                super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/jets/size", jets.size(), weight*npvweight, 20, 0., 20.)
                super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/bjets/size", bjets.size(), weight*npvweight, 15, 0., 15.)
                super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/METv/pt", METv.Pt(), weight*npvweight, 300, 0., 300.)
                super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/pair/pt", pair.Pt(), weight*npvweight, 300, 0., 300.)
                super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/pair/eta", pair.Eta(), weight*npvweight, 100, -5., 5.)
                super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/pair/phi", pair.Phi(), weight*npvweight, 64, -3.2, 3.2)
                super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/pair/mass", pair.M(), weight*npvweight, 80, 50., 130.)
                super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/nPV", super().nPV, weight*npvweight, 100, 0., 100.)
            
            if event.PassTrigger("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v") and looseElectrons.at(0).Pt() > 15.:
                trigSuffix = "Ele12"
                if syst == "PileupReweight": npvweight = super().GetPileUpWeight(super().nPileUp, 0)
                else:                        npvweight = self.getNPVReweight(super().nPV, trigSuffix)
                for idx, mu in enumerate(looseElectrons, start=1):
                    super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/electrons/{idx}/pt", mu.Pt(), weight*npvweight, 300, 0., 300.)
                    super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/electrons/{idx}/eta", mu.Eta(), weight*npvweight, 48, -2.4, 2.4)
                    super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/electrons/{idx}/phi", mu.Phi(), weight*npvweight, 64, -3.2, 3.2)
                    super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/electrons/{idx}/mass", mu.M(), weight*npvweight, 10, 0., 1.)
                for idx, jet in enumerate(jets, start=1):
                    super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/jets/{idx}/pt", jet.Pt(), weight*npvweight, 300, 0., 300.)
                    super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/jets/{idx}/eta", jet.Eta(), weight*npvweight, 48, -2.4, 2.4)
                    super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/jets/{idx}/phi", jet.Phi(), weight*npvweight, 64, -3.2, 3.2)
                    super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/jets/{idx}/mass", jet.M(), weight*npvweight, 100, 0., 100.)
                for idx, bjet in enumerate(bjets, start=1):
                    super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/bjets/{idx}/pt", bjet.Pt(), weight*npvweight, 300, 0., 300.)
                    super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/bjets/{idx}/eta", bjet.Eta(), weight*npvweight, 48, -2.4, 2.4)
                    super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/bjets/{idx}/phi", bjet.Phi(), weight*npvweight, 64, -3.2, 3.2)
                    super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/bjets/{idx}/mass", bjet.M(), weight*npvweight, 100, 0., 100.)
                super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/jets/size", jets.size(), weight*npvweight, 20, 0., 20.)
                super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/bjets/size", bjets.size(), weight*npvweight, 15, 0., 15.)
                super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/METv/pt", METv.Pt(), weight*npvweight, 300, 0., 300.)
                super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/pair/pt", pair.Pt(), weight*npvweight, 300, 0., 300.)
                super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/pair/eta", pair.Eta(), weight*npvweight, 100, -5., 5.)
                super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/pair/phi", pair.Phi(), weight*npvweight, 64, -3.2, 3.2)
                super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/pair/mass", pair.M(), weight*npvweight, 80, 50., 130.) 
                super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/nPV", super().nPV, weight*npvweight, 100, 0., 100.)

            if event.PassTrigger("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v") and looseElectrons.at(0).Pt() > 25.:
                trigSuffix = "Ele23"
                if syst == "PileupReweight": npvweight = super().GetPileUpWeight(super().nPileUp, 0)
                else:                        npvweight = self.getNPVReweight(super().nPV, trigSuffix)
                for idx, mu in enumerate(looseElectrons, start=1):
                    super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/electrons/{idx}/pt", mu.Pt(), weight*npvweight, 300, 0., 300.)
                    super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/electrons/{idx}/eta", mu.Eta(), weight*npvweight, 48, -2.4, 2.4)
                    super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/electrons/{idx}/phi", mu.Phi(), weight*npvweight, 64, -3.2, 3.2)
                    super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/electrons/{idx}/mass", mu.M(), weight*npvweight, 10, 0., 1.)
                for idx, jet in enumerate(jets, start=1):
                    super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/jets/{idx}/pt", jet.Pt(), weight*npvweight, 300, 0., 300.)
                    super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/jets/{idx}/eta", jet.Eta(), weight*npvweight, 48, -2.4, 2.4)
                    super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/jets/{idx}/phi", jet.Phi(), weight*npvweight, 64, -3.2, 3.2)
                    super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/jets/{idx}/mass", jet.M(), weight*npvweight, 100, 0., 100.)
                for idx, bjet in enumerate(bjets, start=1):
                    super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/bjets/{idx}/pt", bjet.Pt(), weight*npvweight, 300, 0., 300.)
                    super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/bjets/{idx}/eta", bjet.Eta(), weight*npvweight, 48, -2.4, 2.4)
                    super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/bjets/{idx}/phi", bjet.Phi(), weight*npvweight, 64, -3.2, 3.2)
                    super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/bjets/{idx}/mass", bjet.M(), weight*npvweight, 100, 0., 100.)
                super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/jets/size", jets.size(), weight*npvweight, 20, 0., 20.)
                super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/bjets/size", bjets.size(), weight*npvweight, 15, 0., 15.)
                super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/METv/pt", METv.Pt(), weight*npvweight, 300, 0., 300.)
                super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/pair/pt", pair.Pt(), weight*npvweight, 300, 0., 300.)
                super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/pair/eta", pair.Eta(), weight*npvweight, 100, -5., 5.)
                super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/pair/phi", pair.Phi(), weight*npvweight, 64, -3.2, 3.2)
                super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/pair/mass", pair.M(), weight*npvweight, 80, 50., 130.)
                super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/nPV", super().nPV, weight*npvweight, 100, 0., 100.)

            if tightElectrons.size() != 2: return None
            if event.PassTrigger("HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v"):
                trigSuffix = "Ele8"
                if syst == "PileupReweight": npvweight = super().GetPileUpWeight(super().nPileUp, 0)
                else:                        npvweight = self.getNPVReweight(super().nPV, trigSuffix)
                for idx, mu in enumerate(looseElectrons, start=1):
                    super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/electrons/{idx}/pt", mu.Pt(), weight*npvweight, 300, 0., 300.)
                    super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/electrons/{idx}/eta", mu.Eta(), weight*npvweight, 48, -2.4, 2.4)
                    super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/electrons/{idx}/phi", mu.Phi(), weight*npvweight, 64, -3.2, 3.2)
                    super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/electrons/{idx}/mass", mu.M(), weight*npvweight, 10, 0., 1.)
                for idx, jet in enumerate(jets, start=1):
                    super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/jets/{idx}/pt", jet.Pt(), weight*npvweight, 300, 0., 300.)
                    super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/jets/{idx}/eta", jet.Eta(), weight*npvweight, 48, -2.4, 2.4)
                    super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/jets/{idx}/phi", jet.Phi(), weight*npvweight, 64, -3.2, 3.2)
                    super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/jets/{idx}/mass", jet.M(), weight*npvweight, 100, 0., 100.)
                for idx, bjet in enumerate(bjets, start=1):
                    super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/bjets/{idx}/pt", bjet.Pt(), weight*npvweight, 300, 0., 300.)
                    super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/bjets/{idx}/eta", bjet.Eta(), weight*npvweight, 48, -2.4, 2.4)
                    super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/bjets/{idx}/phi", bjet.Phi(), weight*npvweight, 64, -3.2, 3.2)
                    super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/bjets/{idx}/mass", bjet.M(), weight*npvweight, 100, 0., 100.)
                super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/jets/size", jets.size(), weight*npvweight, 20, 0., 20.)
                super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/bjets/size", bjets.size(), weight*npvweight, 15, 0., 15.)
                super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/METv/pt", METv.Pt(), weight*npvweight, 300, 0., 300.)
                super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/pair/pt", pair.Pt(), weight*npvweight, 300, 0., 300.)
                super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/pair/eta", pair.Eta(), weight*npvweight, 100, -5., 5.)
                super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/pair/phi", pair.Phi(), weight*npvweight, 64, -3.2, 3.2)
                super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/pair/mass", pair.M(), weight*npvweight, 80, 50., 130.)
                super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/nPV", super().nPV, weight*npvweight, 100, 0., 100.)

            if event.PassTrigger("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v") and tightElectrons.at(0).Pt() > 15.:
                trigSuffix = "Ele12"
                if syst == "PileupReweight": npvweight = super().GetPileUpWeight(super().nPileUp, 0)
                else:                        npvweight = self.getNPVReweight(super().nPV, trigSuffix)
                for idx, mu in enumerate(tightElectrons, start=1):
                    super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/electrons/{idx}/pt", mu.Pt(), weight*npvweight, 300, 0., 300.)
                    super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/electrons/{idx}/eta", mu.Eta(), weight*npvweight, 48, -2.4, 2.4)
                    super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/electrons/{idx}/phi", mu.Phi(), weight*npvweight, 64, -3.2, 3.2)
                    super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/electrons/{idx}/mass", mu.M(), weight*npvweight, 10, 0., 1.)
                for idx, jet in enumerate(jets, start=1):
                    super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/jets/{idx}/pt", jet.Pt(), weight*npvweight, 300, 0., 300.)
                    super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/jets/{idx}/eta", jet.Eta(), weight*npvweight, 48, -2.4, 2.4)
                    super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/jets/{idx}/phi", jet.Phi(), weight*npvweight, 64, -3.2, 3.2)
                    super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/jets/{idx}/mass", jet.M(), weight*npvweight, 100, 0., 100.)
                for idx, bjet in enumerate(bjets, start=1):
                    super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/bjets/{idx}/pt", bjet.Pt(), weight*npvweight, 300, 0., 300.)
                    super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/bjets/{idx}/eta", bjet.Eta(), weight*npvweight, 48, -2.4, 2.4)
                    super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/bjets/{idx}/phi", bjet.Phi(), weight*npvweight, 64, -3.2, 3.2)
                    super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/bjets/{idx}/mass", bjet.M(), weight*npvweight, 100, 0., 100.)
                super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/jets/size", jets.size(), weight*npvweight, 20, 0., 20.)
                super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/bjets/size", bjets.size(), weight*npvweight, 15, 0., 15.)
                super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/METv/pt", METv.Pt(), weight*npvweight, 300, 0., 300.)
                super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/pair/pt", pair.Pt(), weight*npvweight, 300, 0., 300.)
                super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/pair/eta", pair.Eta(), weight*npvweight, 100, -5., 5.)
                super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/pair/phi", pair.Phi(), weight*npvweight, 64, -3.2, 3.2)
                super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/pair/mass", pair.M(), weight*npvweight, 80, 50., 130.) 
                super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/nPV", super().nPV, weight*npvweight, 100, 0., 100.)

            if event.PassTrigger("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v") and tightElectrons.at(0).Pt() > 25.:
                trigSuffix = "Ele23"
                if syst == "PileupReweight": npvweight = super().GetPileUpWeight(super().nPileUp, 0)
                else:                        npvweight = self.getNPVReweight(super().nPV, trigSuffix)
                for idx, mu in enumerate(tightElectrons, start=1):
                    super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/electrons/{idx}/pt", mu.Pt(), weight*npvweight, 300, 0., 300.)
                    super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/electrons/{idx}/eta", mu.Eta(), weight*npvweight, 48, -2.4, 2.4)
                    super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/electrons/{idx}/phi", mu.Phi(), weight*npvweight, 64, -3.2, 3.2)
                    super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/electrons/{idx}/mass", mu.M(), weight*npvweight, 10, 0., 1.)
                for idx, jet in enumerate(jets, start=1):
                    super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/jets/{idx}/pt", jet.Pt(), weight*npvweight, 300, 0., 300.)
                    super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/jets/{idx}/eta", jet.Eta(), weight*npvweight, 48, -2.4, 2.4)
                    super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/jets/{idx}/phi", jet.Phi(), weight*npvweight, 64, -3.2, 3.2)
                    super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/jets/{idx}/mass", jet.M(), weight*npvweight, 100, 0., 100.)
                for idx, bjet in enumerate(bjets, start=1):
                    super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/bjets/{idx}/pt", bjet.Pt(), weight*npvweight, 300, 0., 300.)
                    super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/bjets/{idx}/eta", bjet.Eta(), weight*npvweight, 48, -2.4, 2.4)
                    super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/bjets/{idx}/phi", bjet.Phi(), weight*npvweight, 64, -3.2, 3.2)
                    super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/bjets/{idx}/mass", bjet.M(), weight*npvweight, 100, 0., 100.)
                super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/jets/size", jets.size(), weight*npvweight, 20, 0., 20.)
                super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/bjets/size", bjets.size(), weight*npvweight, 15, 0., 15.)
                super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/METv/pt", METv.Pt(), weight*npvweight, 300, 0., 300.)
                super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/pair/pt", pair.Pt(), weight*npvweight, 300, 0., 300.)
                super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/pair/eta", pair.Eta(), weight*npvweight, 100, -5., 5.)
                super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/pair/phi", pair.Phi(), weight*npvweight, 64, -3.2, 3.2)
                super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/pair/mass", pair.M(), weight*npvweight, 80, 50., 130.) 
                super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/nPV", super().nPV, weight*npvweight, 100, 0., 100.)
            return None
        if self.channel == "MeasNPVEl":
            if event.PassTrigger("HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v"):
                trigSuffix = "Ele8"
                super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/nPV", super().nPV, weight, 100, 0., 100.)
                super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/Rho", super().Rho, weight, 100, 0., 100.)
            if event.PassTrigger("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v") and looseElectrons.at(0).Pt() > 15.:
                trigSuffix = "Ele12"
                super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/nPV", super().nPV, weight, 100, 0., 100.)
                super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/Rho", super().Rho, weight, 100, 0., 100.)
            if event.PassTrigger("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v") and looseElectrons.at(0).Pt() > 25.:
                trigSuffix = "Ele23"
                super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/nPV", super().nPV, weight, 100, 0., 100.)
                super().FillHist(f"{self.channel}_{trigSuffix}/loose/{syst}/Rho", super().Rho, weight, 100, 0., 100.) 
            
            if tightElectrons.size() != 1:    return None
            if event.PassTrigger("HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v"):
                trigSuffix = "Ele8"
                super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/nPV", super().nPV, weight, 100, 0., 100.)
                super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/Rho", super().Rho, weight, 100, 0., 100.)
            if event.PassTrigger("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v") and looseElectrons.at(0).Pt() > 15.:
                trigSuffix = "Ele12"
                super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/nPV", super().nPV, weight, 100, 0., 100.)
                super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/Rho", super().Rho, weight, 100, 0., 100.)
            if event.PassTrigger("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v") and looseElectrons.at(0).Pt() > 25.:
                trigSuffix = "Ele23"
                super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/nPV", super().nPV, weight, 100, 0., 100.)
                super().FillHist(f"{self.channel}_{trigSuffix}/tight/{syst}/Rho", super().Rho, weight, 100, 0., 100.)
            return None
        if self.channel == "MeasFakeEl":
            return None
