from ROOT import gSystem
from ROOT import TriLeptonBase
from ROOT import TString
from ROOT.std import vector
from ROOT.JetTagging import Parameters as jParameters
from ROOT import Lepton, Muon, Electron, Jet
gSystem.Load("/cvmfs/cms.cern.ch/slc7_amd64_gcc900/external/lhapdf/6.2.3/lib/libLHAPDF.so")

class MeasConversion(TriLeptonBase):
    def __init__(self):
        super().__init__()
        # at this point, TriLeptonBase::initializeAnalyzer has not been activate
        
    def initializePyAnalyzer(self):
        super().initializeAnalyzer()
        
        ## channel assertion
        try: assert super().Skim1E2Mu or super().Skim3Mu
        except: raise AssertionError(f"Wrong channel flag")
        
        if super().Skim1E2Mu: self.channel = "Skim1E2Mu"
        if super().Skim3Mu: self.channel = "Skim3Mu"
        
        ## systematic sources
        self.weightVariations = []
        self.scaleVariations = []
        if not super().IsDATA:
            if self.channel == "Skim1E2Mu":
                self.weightVariations = ["L1PrefireUp", "L1PrefireDown",
                                         "PileupReweightUp", "PileupReweightDown",
                                         "MuonIDSFUp", "MuonIDSFDown",
                                         "ElectronIDSFUp", "ElectronIDSFDown",
                                         "EMuTrigSFUp", "EMuTrigSFDown",
                                         "HeavyTagUpUnCorr", "HeavyTagDownUnCorr",
                                         "HeavyTagUpCorr", "HeavyTagDownCorr",
                                         "LightTagUpUnCorr", "LightTagDownUnCorr",
                                         "LightTagUpCorr", "LightTagDownCorr"]
        
            if self.channel == "Skim3Mu":
                self.weightVariations = ["L1PrefireUp", "L1PrefireDown",
                                         "PileupReweightUp", "PileupReweightDown",
                                         "MuonIDSFUp", "MuonIDSFDown",
                                         "ElectronIDSFUp", "ElectronIDSFDown",
                                         "DblMuTrigSFUp", "DblMuTrigSFDown",
                                         "HeavyTagUpUnCorr", "HeavyTagDownUnCorr",
                                         "HeavyTagUpCorr", "HeavyTagDownCorr",
                                         "LightTagUpUnCorr", "LightTagDownUnCorr",
                                         "LightTagUpCorr", "LightTagDownCorr"]
            self.scaleVariations = ["JetResUp", "JetResDown",
                                    "JetEnUp", "JetEnDown",
                                    "ElectronResUp", "ElectronResDown",
                                    "ElectronEnUp", "ElectronEnDown",
                                    "MuonEnUp", "MuonEnDown"]
        
        # for data, only central scale with no weight variation
        self.systematics = ["Central"]
        if not super().IsDATA:
            self.systematics += self.weightVariations + self.scaleVariations
            
    def executeEvent(self):
        if not super().PassMETFilter(): return None
        ev = super().GetEvent()
        rawMuons = super().GetAllMuons()
        rawElectrons = super().GetAllElectrons()
        rawJets = super().GetAllJets()
        METv = ev.GetMETVector()
        truth = super().GetGens() if not super().IsDATA else None
        
        # Central scale
        vetoMuons, tightMuons, vetoElectrons, tightElectrons, jets, bjets = self.defineObjects(rawMuons, rawElectrons, rawJets)
        thisChannel = self.selectEvent(ev, truth, vetoMuons, tightMuons, vetoElectrons, tightElectrons, jets, bjets, METv)
        
        if not thisChannel is None:
            pairs = self.makePair(tightMuons)
            # make objects as a dictionary
            objects = {"muons": tightMuons,
                       "electrons": tightElectrons,
                       "jets": jets,
                       "bjets": bjets,
                       "METv": METv,
                       "pairs": pairs,
                       }
            weight = self.getWeight(thisChannel, ev, tightMuons, tightElectrons, jets, syst="Central")
            self.FillObjects(thisChannel, objects, weight, syst="Central")
            
            for syst in self.weightVariations:
                weight = self.getWeight(thisChannel, ev, tightMuons, tightElectrons, jets, syst)
                self.FillObjects(thisChannel, objects, weight, syst)
        
        for syst in self.scaleVariations:
            vetoMuons, tightMuons, vetoElectrons, tightElectrons, jets, bjets = self.defineObjects(rawMuons, rawElectrons, rawJets, syst)
            thisChannel = self.selectEvent(ev, truth, vetoMuons, tightMuons, vetoElectrons, tightElectrons, jets, bjets, METv)
            if not thisChannel is None:
                pairs = self.makePair(tightMuons)
                objects = {"muons": tightMuons,
                           "electrons": tightElectrons,
                           "jets": jets,
                           "bjets": bjets,
                           "METv": METv,
                           "pairs": pairs,
                           }
                weight = self.getWeight(thisChannel, ev, tightMuons, tightElectrons, jets, syst)
                self.FillObjects(thisChannel, objects, weight, syst)
                
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
        tightMuons = super().SelectMuons(vetoMuons, super().MuonIDs[0], 10., 2.4)
        vetoElectrons = super().SelectElectrons(allElectrons, super().ElectronIDs[2], 10., 2.5)
        tightElectrons = super().SelectElectrons(vetoElectrons, super().ElectronIDs[0], 10., 2.5)
        jets = super().SelectJets(allJets, "tight", 20., 2.4)
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
        
        return (vetoMuons, tightMuons, vetoElectrons, tightElectrons, jets, bjets)
    
    def selectEvent(self, event, truth, vetoMuons, tightMuons, vetoElectrons, tightElectrons, jets, bjets, METv):
        is3Mu = (len(tightMuons) == 3 and len(vetoMuons) == 3 and \
                len(tightElectrons) == 0 and len(vetoElectrons) == 0)
        is1E2Mu = len(tightMuons) == 2 and len(vetoMuons) == 2 and \
                  len(tightElectrons) == 1 and len(vetoElectrons) == 1
        
        if self.channel == "Skim1E2Mu" and not is1E2Mu: return None
        if self.channel == "Skim3Mu" and not is3Mu: return None
        
        # prompt matching
        # for patching samples
        if "DYJets" in super().MCSample or "ZGToLLG" in super().MCSample:
            # first do matching
            promptMuons = vector[Muon]()
            promptElectrons = vector[Electron]()
            for mu in tightMuons:
                if super().GetLeptonType(mu, truth) > 0: promptMuons.emplace_back(mu)
            for ele in tightElectrons:
                if super().GetLeptonType(ele, truth) > 0: promptElectrons.emplace_back(ele)

            if promptMuons.size() != tightMuons.size(): return None
            if promptElectrons.size() != tightElectrons.size(): return None

        # for conversion measurement
        leptons = vector[Lepton]()
        for mu in tightMuons: leptons.emplace_back(mu)
        for ele in tightElectrons: leptons.emplace_back(ele)
        if leptons.at(0).Pt() > 20. and leptons.at(1).Pt() > 20. and leptons.at(2).Pt() > 20.:
            self.measure = "HighPT"
        else:
            self.measure = "LowPT"
            
        ##### event selection
        ## 1E2Mu ZGamma
        ## 1. pass EMuTriggers
        ## 2. Exact 2 tight muons and 1 tight electron, no additional lepton
        ## 3. Exists OS muon pair with mass > 12 GeV
        ## 4. |M(mumue) - 91.2| < 10
        ## 5. No b-jet
        if self.channel == "Skim1E2Mu":
            if not event.PassTrigger(super().EMuTriggers): return None
            leptons = vector[Lepton]()
            for mu in tightMuons: leptons.emplace_back(mu)
            for ele in tightElectrons: leptons.emplace_back(ele)
            mu1, mu2, ele = tuple(leptons)
            passLeadMu = mu1.Pt() > 25. and ele.Pt() > 15.
            passLeadEle = mu1.Pt() > 10. and ele.Pt() > 25.
            passSafeCut = passLeadMu or passLeadEle
            if not passSafeCut: return None
            if not mu1.Charge()+mu2.Charge() == 0: return None
            pair = self.makePair(tightMuons)
            isOnZ = abs((mu1+mu2+ele).M() - 91.2) < 10. 
            if not pair.M() > 12.: return None
            if not isOnZ: return None
            if not bjets.size() == 0: return None

            return "ZGamma1E2Mu"
        ## 3Mu ZGamma
        ## 1. pass DblMuTriggers
        ## 2. Exact 3 tight muons, no additional leptons
        ## 3. Exist OS muon pair,
        ## 4. All OS muon pair mass > 12 GeV
        ## 5. |M(mumumu) - 91.2| < 10
        ## 6. No b-jet
        if self.channel == "Skim3Mu":
            if not event.PassTrigger(super().DblMuTriggers): return None
            mu1, mu2, mu3  = tuple(tightMuons)
            if not mu1.Pt() > 20.: return None
            if not mu2.Pt() > 10.: return None
            if not mu3.Pt() > 10.: return None
            if not abs(mu1.Charge()+mu2.Charge()+mu3.Charge()) == 1: return None
            pair1, pair2 = self.makePair(tightMuons)
            isOnZ = abs((mu1+mu2+mu3).M() - 91.2) < 10.
            if not pair1.M() > 12.: return None
            if not pair2.M() > 12.: return None
            if not isOnZ: return None
            if not bjets.size() == 0: return None
            
            return "ZGamma3Mu"
        
        # should not reach here
        raise EOFError("Analyzer should not be reached here. Check your channel condition")
    
    def makePair(self, muons):
        if muons.size() == 2:
            return (muons[0] + muons[1])
        elif muons.size() == 3:
            mu1, mu2, mu3 = tuple(muons)
            if mu1.Charge() == mu2.Charge():
                pair1 = mu1 + mu3
                pair2 = mu2 + mu3
            elif mu1.Charge() == mu3.Charge():
                pair1 = mu1 + mu2
                pair2 = mu2 + mu3
            else:   # mu2.Charge() == mu3.Charge()
                pair1 = mu1 + mu2
                pair2 = mu1 + mu3
            return (pair1, pair2)
        else:
            raise NotImplementedError(f"Wrong number of muons (muons.size())")
        
    # set weight
    def getWeight(self, channel, event, muons, electrons, jets, syst="Central"):
        weight = 1.
        if not syst in self.systematics:
            print(f"[PromptEstimator::getWeight] Wrong systematic {syst}")
            exit(1)

        if not super().IsDATA:
            weight *= super().MCweight() * super().GetKFactor()
            weight *= event.GetTriggerLumi("Full")
            if syst == "L1PrefireUp":     w_prefire = super().GetPrefireWeight(1)
            elif syst == "L1PrefireDown": w_prefire = super().GetPrefireWeight(-1)
            else:                         w_prefire = super().GetPrefireWeight(0)
            
            if syst == "PileupReweightUp":     w_pileup = super().GetPileUpWeight(super().nPileUp, 1)
            elif syst == "PileupReweightDown": w_pileup = super().GetPileUpWeight(super().nPileUp, -1)
            else:                              w_pileup = super().GetPileUpWeight(super().nPileUp, 0)
            
            w_muonIDSF = 1.
            w_eleIDSF = 1.
            w_trigSF = 1.
            if "1E2Mu" in channel:
                for mu in muons:
                    if syst == "MuonIDSFUp":      w_muonIDSF *= self.getMuonRecoSF(mu, 1) * self.getMuonIDSF(mu, 1)
                    elif syst == "MuonIDSFDown":  w_muonIDSF *= self.getMuonRecoSF(mu, -1) * self.getMuonIDSF(mu, -1)
                    else:                         w_muonIDSF *= self.getMuonRecoSF(mu, 0) * self.getMuonIDSF(mu, 0)
                for el in electrons:
                    if syst == "ElectronIDSFUp":     w_eleIDSF *= self.mcCorr.ElectronReco_SF(el.scEta(), el.Pt(), 1) * self.getEleIDSF(el, 1)
                    elif syst == "ElectronIDSFDown": w_eleIDSF *= self.mcCorr.ElectronReco_SF(el.scEta(), el.Pt(), -1) * self.getEleIDSF(el, -1)
                    else:                            w_eleIDSF *= self.mcCorr.ElectronReco_SF(el.scEta(), el.Pt(), 0) * self.getEleIDSF(el, 0) 

                if syst == "EMuTriggerSFUp":     w_trigSF = self.getEMuTriggerSF(electrons, muons, 1)
                elif syst == "EMuTriggerSFDown": w_trigSF = self.getEMuTriggerSF(electrons, muons, -1)
                else:                            w_trigSF = self.getEMuTriggerSF(electrons, muons, 0)
            
            if "3Mu" in channel:
                for mu in muons:
                    if syst == "MuonIDSFUp":      w_muonIDSF *= self.getMuonRecoSF(mu, 1) * self.getMuonIDSF(mu, 1)
                    elif syst == "MuonIDSFDown":  w_muonIDSF *= self.getMuonRecoSF(mu, -1) * self.getMuonIDSF(mu, -1)
                    else:                         w_muonIDSF *= self.getMuonRecoSF(mu, 0) * self.getMuonIDSF(mu, 0)
                if syst == "DblMuTrigSFUp":     w_trigSF = self.getDblMuTriggerSF(muons, 1)
                elif syst == "DblMuTrigSFDown": w_trigSF = self.getDblMuTriggerSF(muons, -1)
                else:                           w_trigSF = self.getDblMuTriggerSF(muons, 0) 
            
            weight *= w_prefire            # print(f"w_prefire: {w_prefire}")
            weight *= w_pileup             # print(f"w_pileup: {w_pileup}")
            weight *= w_muonIDSF           # print(f"muonID: {w_muonIDSF}")
            weight *= w_eleIDSF;           # print(syst, w_eleIDSF)
            weight *= w_trigSF        # print(f"muontrig: {w_dblMuTrigSF}")
            
            # b-tagging
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
        
    def FillObjects(self, channel, objects, weight, syst):
        muons = objects["muons"]
        electrons = objects["electrons"]
        jets = objects["jets"]
        bjets = objects["bjets"]
        METv = objects["METv"]
        ## fill input observables
        for idx, mu in enumerate(muons, start=1):
            super().FillHist(f"{channel}/{self.measure}/{syst}/muons/{idx}/pt", mu.Pt(), weight, 300, 0., 300.)
            super().FillHist(f"{channel}/{self.measure}/{syst}/muons/{idx}/eta", mu.Eta(), weight, 48, -2.4, 2.4)
            super().FillHist(f"{channel}/{self.measure}/{syst}/muons/{idx}/phi", mu.Phi(), weight, 64, -3.2, 3.2)
            super().FillHist(f"{channel}/{self.measure}/{syst}/muons/{idx}/mass", mu.M(), weight, 10, 0., 1.)
        for idx, ele in enumerate(electrons, start=1):
            super().FillHist(f"{channel}/{self.measure}/{syst}/electrons/{idx}/pt", ele.Pt(), weight, 300, 0., 300.)
            super().FillHist(f"{channel}/{self.measure}/{syst}/electrons/{idx}/eta", ele.Eta(), weight, 50, -2.5, 2.5)
            super().FillHist(f"{channel}/{self.measure}/{syst}/electrons/{idx}/Phi", ele.Phi(), weight, 64, -3.2, 3.2)
            super().FillHist(f"{channel}/{self.measure}/{syst}/electrons/{idx}/mass", ele.M(), weight, 100, 0., 1.)
        for idx, jet in enumerate(jets, start=1):
            super().FillHist(f"{channel}/{self.measure}/{syst}/jets/{idx}/pt", jet.Pt(), weight, 300, 0., 300.)
            super().FillHist(f"{channel}/{self.measure}/{syst}/jets/{idx}/eta", jet.Eta(), weight, 48, -2.4, 2.4)
            super().FillHist(f"{channel}/{self.measure}/{syst}/jets/{idx}/phi", jet.Phi(), weight, 64, -3.2, 3.2)
            super().FillHist(f"{channel}/{self.measure}/{syst}/jets/{idx}/mass", jet.M(), weight, 100, 0., 100.)
        for idx, bjet in enumerate(bjets, start=1):
            super().FillHist(f"{channel}/{self.measure}/{syst}/bjets/{idx}/pt", bjet.Pt(), weight, 300, 0., 300.)
            super().FillHist(f"{channel}/{self.measure}/{syst}/bjets/{idx}/eta", bjet.Eta(), weight, 48, -2.4, 2.4)
            super().FillHist(f"{channel}/{self.measure}/{syst}/bjets/{idx}/phi", bjet.Phi(), weight, 64, -3.2, 3.2)
            super().FillHist(f"{channel}/{self.measure}/{syst}/bjets/{idx}/mass", bjet.M(), weight, 100, 0., 100.)
        super().FillHist(f"{channel}/{self.measure}/{syst}/jets/size", jets.size(), weight, 20, 0., 20.)
        super().FillHist(f"{channel}/{self.measure}/{syst}/bjets/size", bjets.size(), weight, 15, 0., 15.)
        super().FillHist(f"{channel}/{self.measure}/{syst}/METv/pt", METv.Pt(), weight, 300, 0., 300.)
        super().FillHist(f"{channel}/{self.measure}/{syst}/METv/phi", METv.Phi(), weight, 64, -3.2, 3.2)
        
        ## fill ZCand
        ZCand = None
        if "1E2Mu" in channel: ZCand = electrons.at(0) + muons.at(0) + muons.at(1)
        if "3Mu" in channel:   ZCand = muons.at(0) + muons.at(1) + muons.at(2)
        
        super().FillHist(f"{channel}/{self.measure}/{syst}/ZCand/pt", ZCand.Pt(), weight, 300, 0., 300.)
        super().FillHist(f"{channel}/{self.measure}/{syst}/ZCand/eta", ZCand.Eta(), weight, 100, -5., 5.)
        super().FillHist(f"{channel}/{self.measure}/{syst}/ZCand/phi", ZCand.Phi(), weight, 64, -3.2, 3.2)
        super().FillHist(f"{channel}/{self.measure}/{syst}/ZCand/mass", ZCand.M(), weight, 200, 0., 200.)
