from ROOT import gSystem
from ROOT import TriLeptonBase
from ROOT import TTree, TString
from ROOT.std import vector
from ROOT.JetTagging import Parameters as jParameters
from ROOT import Lepton, Muon, Electron, Jet
gSystem.Load("/cvmfs/cms.cern.ch/slc7_amd64_gcc900/external/lhapdf/6.2.3/lib/libLHAPDF.so")

from array import array
from itertools import product
from MLTools.helpers import loadModels
from MLTools.helpers import getGraphInput, getGraphScore

class PromptUnbinned(TriLeptonBase):
    def __init__(self):
        # at this point, TriLeptonBase::initializeAnalyzer has not been activated
        super().__init__()
        
    def initializePyAnalyzer(self):
        super().initializeAnalyzer()
        
        ## channel assertion
        try: assert super().Skim1E2Mu or super().Skim3Mu
        except: raise AssertionError(f"Wrong channel flag")
        
        if super().Skim1E2Mu: self.channel = "Skim1E2Mu"
        if super().Skim3Mu: self.channel = "Skim3Mu"
        
        # Not implementing neural networks now
        # self.network = "GraphNeuralNet"
        
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
        
        self.systematics = ["Central"]
        if not super().IsDATA:
            self.systematics += self.weightVariations + self.scaleVariations
        
        self.signalStrings = ["MHc-160_MA-85", "MHc-130_MA-90", "MHc-100_MA-95"]
        self.backgroundStrings = ["nonprompt", "diboson", "ttZ"]

        self.models = loadModels("GraphNeuralNet", self.channel, self.signalStrings, self.backgroundStrings)
        self.__prepareTTree()
        
    def executeEvent(self):
        if not super().PassMETFilter(): return None
        ev = super().GetEvent()
        rawMuons = super().GetAllMuons()
        rawElectrons = super().GetAllElectrons()
        rawJets = super().GetAllJets()
        METv = ev.GetMETVector()
        truth = super().GetGens() if not super().IsDATA else None
        
        ## initialize tree contents
        self.__initTreeContents()
        
        ## fill contents
        vetoMuons, tightMuons, vetoElectrons, tightElectrons, jets, bjets = self.defineObjects(rawMuons, rawElectrons, rawJets)
        thisChannel = self.selectEvent(ev, truth, vetoMuons, tightMuons, vetoElectrons, tightElectrons, jets, bjets, METv)
        if not thisChannel is None:
            vjets = vector[Jet]()
            for j in jets: vjets.emplace_back(j)
            jtp = jParameters(3, 1, 0, 1)   # DeepJet, Medium, incl, mujets
            pairs = self.makePair(tightMuons)
            _, scores = self.evalScore(tightMuons, tightElectrons, jets, bjets, METv)
            
            if thisChannel == "SR1E2Mu":
                self.mass1["Central"][0] = pairs.M()
                self.mass2["Central"][0] = -999.
            else:   # SR3Mu
                self.mass1["Central"][0] = pairs[0].M()
                self.mass2["Central"][0] = pairs[1].M()
                
            for SIG in self.signalStrings:
                self.scoreX[f"{SIG}_Central"][0] = scores[f"{SIG}_vs_nonprompt"]
                self.scoreY[f"{SIG}_Central"][0] = scores[f"{SIG}_vs_diboson"]
                self.scoreZ[f"{SIG}_Central"][0] = scores[f"{SIG}_vs_ttZ"]
            
            self.weight["Central"][0] = 1.
            if super().IsDATA:
                self.tree["Central"].Fill()
                return None
            
            # weight / scale variations for MC
            w_norm = super().MCweight() * super().GetKFactor() * ev.GetTriggerLumi("Full")
            w_l1prefire = super().GetPrefireWeight(0)
            w_l1prefire_up = super().GetPrefireWeight(1)
            w_l1prefire_down = super().GetPrefireWeight(-1)
            w_pileup = super().GetPileUpWeight(super().nPileUp, 0)
            w_pileup_up = super().GetPileUpWeight(super().nPileUp, 1)
            w_pileup_down = super().GetPileUpWeight(super().nPileUp, -1)
            
            sf_muonid = 1.
            sf_muonid_up = 1.
            sf_muonid_down = 1.
            for mu in tightMuons:
                sf_muonid *= self.getMuonRecoSF(mu, 0) * self.getMuonIDSF(mu, 0)
                sf_muonid_up *= self.getMuonRecoSF(mu, 1) * self.getMuonIDSF(mu, 1)
                sf_muonid_down *= self.getMuonRecoSF(mu, -1) * self.getMuonIDSF(mu, -1)
            sf_eleid = 1.
            sf_eleid_up = 1.
            sf_eleid_down = 1.
            for ele in tightElectrons:
                sf_eleid *= self.mcCorr.ElectronReco_SF(ele.scEta(), ele.Pt(), 0) * self.getEleIDSF(ele, 0)
                sf_eleid_up *= self.mcCorr.ElectronReco_SF(ele.scEta(), ele.Pt(), 1) * self.getEleIDSF(ele, 1)
                sf_eleid_down *= self.mcCorr.ElectronReco_SF(ele.scEta(), ele.Pt(), -1) * self.getEleIDSF(ele, -1) 
            
            sf_trig  = 1.
            if thisChannel == "SR1E2Mu":
                sf_trig = self.getEMuTriggerSF(tightElectrons, tightMuons, 0)
                sf_trig_up = self.getEMuTriggerSF(tightElectrons, tightMuons, 1)
                sf_trig_down = self.getEMuTriggerSF(tightElectrons, tightMuons, -1)
            else:
                sf_trig = self.getDblMuTriggerSF(tightMuons, 0)
                sf_trig_up = self.getDblMuTriggerSF(tightMuons, 1)
                sf_trig_down = self.getDblMuTriggerSF(tightMuons, -1)
            sf_btag = self.mcCorr.GetBTaggingReweight_1a(vjets, jtp)
            
            self.weight["Central"][0] = w_norm * w_l1prefire * w_pileup * sf_muonid * sf_eleid * sf_trig * sf_btag
            self.tree["Central"].Fill()
            
            self.weight["L1PrefireUp"][0] = w_norm * w_l1prefire_up * w_pileup * sf_muonid * sf_eleid * sf_trig * sf_btag
            self.weight["L1PrefireDown"][0] = w_norm * w_l1prefire_down * w_pileup * sf_muonid * sf_eleid * sf_trig * sf_btag
            self.weight["PileupReweightUp"][0] = w_norm * w_l1prefire * w_pileup_up * sf_muonid * sf_eleid * sf_trig * sf_btag
            self.weight["PileupReweightDown"][0] = w_norm * w_l1prefire * w_pileup_down * sf_muonid * sf_eleid * sf_trig * sf_btag
            if thisChannel == "SR1E2Mu":
                self.weight["MuonIDSFUp"][0] = w_norm * w_l1prefire * w_pileup * sf_muonid_up * sf_eleid * sf_trig * sf_btag
                self.weight["MuonIDSFDown"][0] =  w_norm * w_l1prefire * w_pileup * sf_muonid_down * sf_eleid * sf_trig * sf_btag
                self.weight["ElectronIDSFUp"][0] =  w_norm * w_l1prefire * w_pileup * sf_muonid * sf_eleid_up * sf_trig * sf_btag
                self.weight["ElectronIDSFDown"][0] =  w_norm * w_l1prefire * w_pileup * sf_muonid * sf_eleid_down * sf_trig * sf_btag
                self.weight["EMuTrigSFUp"][0] =  w_norm * w_l1prefire * w_pileup * sf_muonid * sf_eleid * sf_trig_up * sf_btag
                self.weight["EMuTrigSFDown"][0] =  w_norm * w_l1prefire * w_pileup * sf_muonid * sf_eleid * sf_trig_down * sf_btag
            if thisChannel == "SR3Mu":
                self.weight["MuonIDSFUp"][0] = w_norm * w_l1prefire * w_pileup * sf_muonid_up * sf_eleid * sf_trig * sf_btag
                self.weight["MuonIDSFDown"][0] =  w_norm * w_l1prefire * w_pileup * sf_muonid_down * sf_eleid * sf_trig * sf_btag
                self.weight["DblMuTrigSFUp"][0] =  w_norm * w_l1prefire * w_pileup * sf_muonid * sf_eleid * sf_trig_up * sf_btag
                self.weight["DblMuTrigSFDown"][0] =  w_norm * w_l1prefire * w_pileup * sf_muonid * sf_eleid * sf_trig_down * sf_btag
            for syst in self.weightVariations:
                self.mass1[syst][0] = self.mass1["Central"][0]
                self.mass2[syst][0] = self.mass2["Central"][0]
                for SIG in self.signalStrings:
                    self.scoreX[f"{SIG}_{syst}"][0] = self.scoreX[f"{SIG}_Central"][0]
                    self.scoreY[f"{SIG}_{syst}"][0] = self.scoreY[f"{SIG}_Central"][0]
                    self.scoreZ[f"{SIG}_{syst}"][0] = self.scoreZ[f"{SIG}_Central"][0]
                self.tree[syst].Fill()
            # end of central scale
        
        # No need scale variation for data
        if super().IsDATA:
            return None
        
        # scale varaition
        for syst in self.scaleVariations:
            vetoMuons, tightMuons, vetoElectrons, tightElectrons, jets, bjets = self.defineObjects(rawMuons, rawElectrons, rawJets, syst)
            thisChannel = self.selectEvent(ev, truth, vetoMuons, tightMuons, vetoElectrons, tightElectrons, jets, bjets, METv)
            if thisChannel is None: continue
            vjets = vector[Jet]()
            for j in jets: vjets.emplace_back(j)
            jtp = jParameters(3, 1, 0, 1) # DeepJet, Medium, incl, mujets
            
            pairs = self.makePair(tightMuons)
            data, scores = self.evalScore(tightMuons, tightElectrons, jets, bjets, METv)
            
            if thisChannel == "SR1E2Mu":
                self.mass1[syst][0] = pairs.M()
                self.mass2[syst][0] = -999.
            else:
                self.mass1[syst][0] = pairs[0].M()
                self.mass2[syst][0] = pairs[1].M()
                
            for SIG in self.signalStrings:
                self.scoreX[f"{SIG}_{syst}"][0] = scores[f"{SIG}_vs_nonprompt"]
                self.scoreY[f"{SIG}_{syst}"][0] = scores[f"{SIG}_vs_diboson"]
                self.scoreZ[f"{SIG}_{syst}"][0] = scores[f"{SIG}_vs_ttZ"]
            
            w_norm = super().MCweight() * ev.GetTriggerLumi("Full")
            w_l1prefire = super().GetPrefireWeight(0)
            w_pileup = super().GetPileUpWeight(super().nPileUp, 0)
            sf_muonid = 1.
            for mu in tightMuons:
                sf_muonid *= super().getMuonRecoSF(mu, 0) * self.getMuonIDSF(mu, 0)
            sf_eleid = 1.
            for ele in tightElectrons:
                sf_eleid *= self.mcCorr.ElectronReco_SF(ele.scEta(), ele.Pt(), 0) * self.getEleIDSF(ele, 0)
            sf_trig = 1.
            if thisChannel == "SR1E2Mu":
                sf_trig = self.getEMuTriggerSF(tightElectrons, tightMuons, 0)
            else:
                sf_trig = self.getDblMuTriggerSF(tightMuons, 0)
            sf_btag = self.mcCorr.GetBTaggingReweight_1a(vjets, jtp)
            self.weight[syst][0] = w_norm * w_l1prefire * w_pileup * sf_muonid * sf_eleid * sf_trig * sf_btag
            self.tree[syst].Fill()
    
    def __prepareTTree(self):
        self.tree = {}
        self.mass1 = {}
        self.mass2 = {}
        self.scoreX = {}
        self.scoreY = {}
        self.scoreZ = {}
        self.weight = {}
        
        for syst in self.systematics:
            thisTree = TTree(f"Events_{syst}", "")
            self.mass1[syst] = array("d", [0.]); thisTree.Branch("mass1", self.mass1[syst], "mass1/D")
            self.mass2[syst] = array("d", [0.]); thisTree.Branch("mass2", self.mass2[syst], "mass2/D")
            for SIG in self.signalStrings:
                # vs nonprompt
                self.scoreX[f"{SIG}_{syst}"] = array("d", [0.])
                thisTree.Branch(f"score_{SIG}_vs_nonprompt", self.scoreX[f"{SIG}_{syst}"], f"score_{SIG}_vs_nonprompt/D")
                # vs diboson
                self.scoreY[f"{SIG}_{syst}"] = array("d", [0.])
                thisTree.Branch(f"score_{SIG}_vs_diboson", self.scoreY[f"{SIG}_{syst}"], f"score_{SIG}_vs_diboson/D")
                # vs ttZ
                self.scoreZ[f"{SIG}_{syst}"] = array("d", [0.])
                thisTree.Branch(f"score_{SIG}_vs_ttZ", self.scoreZ[f"{SIG}_{syst}"], f"score_{SIG}_vs_ttZ/D")
            self.weight[syst] = array("d", [0.]); thisTree.Branch("weight", self.weight[syst], "weight/D")
            thisTree.SetDirectory(0)
            self.tree[syst] = thisTree
        
    def __initTreeContents(self):
        for syst in self.systematics:
            self.mass1[syst][0] = -999.
            self.mass2[syst][0] = -999.
            for SIG in self.signalStrings:
                self.scoreX[f"{SIG}_{syst}"][0] = -999.
                self.scoreY[f"{SIG}_{syst}"][0] = -999.
            self.weight[syst][0] = -999.
            
    def defineObjects(self, rawMuons, rawElectrons, rawJets, syst="Central"):
        # first copy objects
        allMuons = rawMuons
        allElectrons = rawElectrons
        allJets = rawJets

        # systematics assertion
        try: assert syst in ["Central"] + self.scaleVariations
        except: raise NameError(f"Wrong systematic {syst}")
        
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
        # check lepton multiplicity first
        is3Mu = (len(tightMuons) == 3 and len(vetoMuons) == 3 and \
                len(tightElectrons) == 0 and len(vetoElectrons) == 0)
        is1E2Mu = len(tightMuons) == 2 and len(vetoMuons) == 2 and \
                  len(tightElectrons) == 1 and len(vetoElectrons) == 1
                  
                  
        if self.channel == "Skim1E2Mu" and not is1E2Mu: return None
        if self.channel == "Skim3Mu" and not is3Mu: return None
        
        # for conversion samples
        if super().MCSample in ["DYJets_MG", "DYJets10to50_MG", "TTG", "WWG"]:
            # at least one conversion lepton should exist
            # internal conversion: 4, 5
            # external conversion: -5, -6
            convMuons = vector[Muon]()
            convElectrons = vector[Electron]()
            for mu in tightMuons:
                if super().GetLeptonType(mu, truth) in [4, 5, -5, -6]: convMuons.emplace_back(mu)
            for ele in tightElectrons:
                if super().GetLeptonType(ele, truth) in [4, 5, -5, -6]: convElectrons.emplace_back(ele)
            if self.channel == "Skim1E2Mu":
                if convElectrons.size() == 0: return None
            if self.channel == "Skim3Mu":
                if convMuons.size() == 0: return None

        ## 1E2Mu baseline
        ## 1. pass EMuTriggers
        ## 2. Exact 2 tight muons and 1 tight electron, no additional lepton
        ## 3. Exists OS muon pair with mass > 12 GeV
        ## 4. At least two jets, at least one b-jet
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
            if not pair.M() > 12.: return None
            if not jets.size() >= 2: return None
            if not bjets.size() >= 1: return None
            return "SR1E2Mu"
        
        ## 3Mu baseline
        ## 1. pass DblMuTriggers
        ## 2. Exact 3 tight muons, no additional leptons
        ## 3. Exist OS muon pair,
        ## 4. All OS muon pair mass > 12 GeV
        ## 5. At least two jets, at least one b-jet
        if self.channel == "Skim3Mu":
            if not event.PassTrigger(super().DblMuTriggers): return None
            mu1, mu2, mu3  = tuple(tightMuons)
            if not mu1.Pt() > 20.: return None
            if not mu2.Pt() > 10.: return None
            if not mu3.Pt() > 10.: return None
            if not abs(mu1.Charge()+mu2.Charge()+mu3.Charge()) == 1: return None
            pair1, pair2 = self.makePair(tightMuons)
            if not pair1.M() > 12.: return None
            if not pair2.M() > 12.: return None
            if not jets.size() >= 2: return None
            if not bjets.size() >= 1: return None
            return "SR3Mu"
        
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
    
    #### Get scores for each event
    def evalScore(self, muons, electrons, jets, bjets, METv):
        scores = {}
        data = getGraphInput(muons, electrons, jets, bjets, METv)
        for sig, bkg in product(self.signalStrings, self.backgroundStrings):
            scores[f"{sig}_vs_{bkg}"] = getGraphScore(self.models[f"{sig}_vs_{bkg}"], data)
        
        return data, scores

    def WriteHist(self):
        super().outfile.cd()
        for syst in self.systematics:
            self.tree[syst].Write()
