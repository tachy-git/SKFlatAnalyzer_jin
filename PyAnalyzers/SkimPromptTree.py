from ROOT import gSystem
from ROOT import TriLeptonBase
from ROOT import std, TString
from ROOT import TTree
from ROOT.JetTagging import Parameters as jParameters
from ROOT import Lepton, Muon, Electron, Jet
gSystem.Load("/cvmfs/cms.cern.ch/slc7_amd64_gcc900/external/lhapdf/6.2.3/lib/libLHAPDF.so")

import os
import numpy as np
import pandas as pd
import torch

from array import array
from itertools import product, combinations
from MLTools.models import SNN, SNNLite, ParticleNet, ParticleNetLite
from MLTools.helpers import evtToGraph
from MLTools.formats import NodeParticle


class SkimPromptTree(TriLeptonBase):
    def __init__(self):
        super().__init__()
        # at this point, TriLeptonBase::initializeAnalyzer has not been activate
                
    def initializePyAnalyzer(self):
        super().initializeAnalyzer()
        
        if super().Skim1E2Mu: 
            self.channel = "Skim1E2Mu"
        elif super().Skim3Mu: 
            self.channel = "Skim3Mu"
        else:
            print("Wrong channel")
            exit(1)

        if super().DenseNet:   self.network = "DenseNeuralNet"
        elif super().GraphNet: self.network = "GraphNeuralNet"
        else:
            print("Wrong network")
            exit(1)
        
        self.scaleVariations = ["Central",
                                "JetResUp", "JetResDown",
                                "JetEnUp", "JetEnDown",
                                "ElectronResUp", "ElectronResDown",
                                "ElectronEnUp", "ElectronEnDown",
                                "MuonEnUp", "MuonEnDown"]
        
        self.signalStrings = [
                "MHc-70_MA-15", "MHc-70_MA-40", "MHc-70_MA-65",
                "MHc-100_MA-15", "MHc-100_MA-60", "MHc-100_MA-95",
                "MHc-130_MA-15", "MHc-130_MA-55", "MHc-130_MA-90", "MHc-130_MA-125",
                "MHc-160_MA-15", "MHc-160_MA-85", "MHc-160_MA-120", "MHc-160_MA-155"]
        self.backgroundStrings = ["TTLL_powheg", "ttX"]
        
        self.__loadModels()
        self.__prepareTTree()
        
    def executeEvent(self):
        if not super().PassMETFilter(): return None
        ev = super().GetEvent()
        rawMuons = super().GetAllMuons()
        rawElectrons = super().GetAllElectrons()
        rawJets = super().GetAllJets()
        METv = ev.GetMETVector()
        truth = super().GetGens() if not super().IsDATA else None
        
        # initialize contents
        self.__initTreeContents(self.n_entry)
        
        # fill contents   
        if not super().IsDATA:
            self.w_norm[0] = super().MCweight()*ev.GetTriggerLumi("Full")
            self.w_l1prefire[0] = super().GetPrefireWeight(0)
            self.w_l1prefire_up[0] = super().GetPrefireWeight(1)
            self.w_l1prefire_down[0] = super().GetPrefireWeight(-1)
            self.w_pileup[0] = super().GetPileUpWeight(super().nPileUp, 0)
            self.w_pileup_up[0] = super().GetPileUpWeight(super().nPileUp, 1)
            self.w_pileup_down[0] = super().GetPileUpWeight(super().nPileUp, -1) 
            
        for syst in self.scaleVariations:
            vetoMuons, tightMuons, vetoElectrons, tightElectrons, jets, bjets = self.defineObjects(rawMuons, rawElectrons, rawJets, syst)
            channel = self.selectEvent(ev, truth, vetoMuons, tightMuons, vetoElectrons, tightElectrons, jets, bjets, METv)
            if channel is None: continue
            vjets = std.vector[Jet]()
            for j in jets: vjets.emplace_back(j)
            jtp = jParameters(3, 1, 0, 1)    # DeepJet, Medium, incl, mujets
            
            pairs = self.makePair(tightMuons)
            data, scores = self.evalScore(tightMuons, tightElectrons, jets, bjets, METv)
            
            if channel == "SR1E2Mu": 
                self.dict_mass1[syst][0] = pairs.M()
                self.dict_mass2[syst][0] = 0
            else:
                self.dict_mass1[syst][0] = pairs[0].M()
                self.dict_mass2[syst][0] = pairs[1].M()
            for SIG, BKG in product(self.signalStrings, self.backgroundStrings):
                self.dict_score[f"{syst}/score_{SIG}_vs_{BKG}"][0] = scores[f"{SIG}_vs_{BKG}"]
            sf_muonid = 1.
            for mu in tightMuons:
                sf_muonid *= self.getMuonIDSF(mu, 0)
            self.sf_muonID[syst][0] = sf_muonid
            self.sf_dblmutrig[syst][0] = self.getDblMuTriggerSF(tightMuons, 0)
            self.sf_btag_central[syst][0] = self.mcCorr.GetBTaggingReweight_1a(vjets, jtp)

            if not syst == "Central": continue
            sf_muonid_up = 1.
            sf_muonid_down = 1.
            for mu in tightMuons: 
                sf_muonid_up *= self.getMuonIDSF(mu, 1)
                sf_muonid_down *= self.getMuonIDSF(mu, -1)
            self.sf_muonID_up[syst][0] = sf_muonid_up
            self.sf_muonID_down[syst][0] = sf_muonid_down
            self.sf_dblmutrig_up[syst][0] = self.getDblMuTriggerSF(tightMuons, 1)
            self.sf_dblmutrig_down[syst][0] = self.getDblMuTriggerSF(tightMuons, -1)
            self.sf_btag_htag_up_uncorr[syst][0] = self.mcCorr.GetBTaggingReweight_1a(vjets, jtp, "SystUpHTag")
            self.sf_btag_htag_down_uncorr[syst][0] = self.mcCorr.GetBTaggingReweight_1a(vjets, jtp, "SystDownHTag")
            self.sf_btag_htag_up_corr[syst][0] = self.mcCorr.GetBTaggingReweight_1a(vjets, jtp, "SystUpHTagCorr")
            self.sf_btag_htag_down_corr[syst][0] = self.mcCorr.GetBTaggingReweight_1a(vjets, jtp, "SystDownHTagCorr")
            self.sf_btag_ltag_up_uncorr[syst][0] = self.mcCorr.GetBTaggingReweight_1a(vjets, jtp, "SystUpLTag")
            self.sf_btag_ltag_down_uncorr[syst][0] = self.mcCorr.GetBTaggingReweight_1a(vjets, jtp, "SystDownLTag")
            self.sf_btag_ltag_up_corr[syst][0] = self.mcCorr.GetBTaggingReweight_1a(vjets, jtp, "SystUpLTagCorr")
            self.sf_btag_ltag_down_corr[syst][0] = self.mcCorr.GetBTaggingReweight_1a(vjets, jtp, "SystDownLTagCorr")
        self.tree.Fill()

    def __loadModels(self):
        self.models = {}

        for sig, bkg in product(self.signalStrings, self.backgroundStrings):
            csv = pd.read_csv(f"{os.environ['DATA_DIR']}/FullRun2/{self.network}/{self.channel}/summary/summary_{sig}_vs_{bkg}.txt", 
                              sep=",\s", 
                              engine="python", 
                              header=None).transpose()
            modelPath = f"{os.environ['DATA_DIR']}/FullRun2/{self.network}/{self.channel}/models/{sig}_vs_{bkg}.pt"
            if self.network == "DenseNeuralNet":
                modelArch = csv[0][4]
                if self.channel == "Skim1E2Mu":
                    if modelArch == "SNN": model = SNN(41, 2)
                    else:                  model = SNNLite(41, 2)
                if self.channel == "Skim3Mu" :
                    if modelArch == "SNN": model = SNN(47, 2)
                    else:                  model = SNNLite(47, 2)
            else:               # GraphNeuralNet
                modelArch, dropout_p, readout = csv[0][3:6]
                if modelArch == "ParticleNet": model = ParticleNet(9, 2, dropout_p, readout)
                else:                          model = ParticleNetLite(9, 2, dropout_p, readout)
            model.load_state_dict(torch.load(modelPath, map_location=torch.device("cpu")))
            model.eval()
            self.models[f"{sig}_vs_{bkg}"] = model
            
    def __prepareTTree(self):
        self.tree = TTree("Events", "")
        self.dict_mass1 = {}
        self.dict_mass2 = {}
        self.dict_score = {}
        self.w_norm = array("f", [0.])    # mcweight * lumi
        self.w_l1prefire = array("f", [0.])
        self.w_l1prefire_up = array("f", [0.])
        self.w_l1prefire_down = array("f", [0.])
        self.w_pileup = array("f", [0.])
        self.w_pileup_up = array("f", [0.])
        self.w_pileup_down = array("f", [0.])
        self.sf_muonID = {}
        self.sf_muonID_up = {}
        self.sf_muonID_down = {}
        self.sf_dblmutrig = {}
        self.sf_dblmutrig_up = {}
        self.sf_dblmutrig_down = {}
        self.sf_btag_central = {}
        self.sf_btag_htag_up_uncorr = {}
        self.sf_btag_htag_down_uncorr = {}
        self.sf_btag_htag_up_corr = {}
        self.sf_btag_htag_down_corr = {}
        self.sf_btag_ltag_up_uncorr = {}
        self.sf_btag_ltag_down_uncorr = {}
        self.sf_btag_ltag_up_corr = {}
        self.sf_btag_ltag_down_corr = {}
        self.n_entry = 0
        
        for syst in self.scaleVariations:
            self.dict_mass1[syst] = array("f", [0.])
            self.dict_mass2[syst] = array("f", [0.])
            for SIG, BKG in product(self.signalStrings, self.backgroundStrings):
                self.dict_score[f"{syst}/score_{SIG}_vs_{BKG}"] = array("f", [0.])
            self.sf_muonID[syst] = array("f", [0.])
            self.sf_muonID_up[syst] = array("f", [0.])
            self.sf_muonID_down[syst] = array("f", [0.])
            self.sf_dblmutrig[syst] = array("f", [0.])
            self.sf_dblmutrig_up[syst] = array("f", [0.])
            self.sf_dblmutrig_down[syst] = array("f", [0.])
            self.sf_btag_central[syst] = array("f", [0.])
            self.sf_btag_htag_up_uncorr[syst] = array("f", [0.])
            self.sf_btag_htag_down_uncorr[syst] = array("f", [0.])
            self.sf_btag_htag_up_corr[syst] = array("f", [0.])
            self.sf_btag_htag_down_corr[syst] = array("f", [0.])
            self.sf_btag_ltag_up_uncorr[syst] = array("f", [0.])
            self.sf_btag_ltag_down_uncorr[syst] = array("f", [0.])
            self.sf_btag_ltag_up_corr[syst] = array("f", [0.])
            self.sf_btag_ltag_down_corr[syst] = array("f", [0.])
        
        self.tree.Branch("w_norm", self.w_norm, "w_norm/F")
        self.tree.Branch("w_l1prefire", self.w_l1prefire, "w_l1prefire/F")
        self.tree.Branch("w_l1prefire_up", self.w_l1prefire_up, "w_l1prefire_up/F")
        self.tree.Branch("w_l1prefire_down", self.w_l1prefire_down, "w_l1prefire_down/F")
        self.tree.Branch("w_pileup", self.w_pileup, "w_pileup/F")
        self.tree.Branch("w_pileup_up", self.w_pileup_up, "w_pileup_up/F")
        self.tree.Branch("w_pileup_down", self.w_pileup_down, "w_pileup_down/F") 
        for syst in self.scaleVariations:    
            self.tree.Branch(f"{syst}_mass1", self.dict_mass1[syst], f"{syst}_mass1/F")
            self.tree.Branch(f"{syst}_mass2", self.dict_mass1[syst], f"{syst}_mass1/F")
            for SIG, BKG in product(self.signalStrings, self.backgroundStrings): 
                self.tree.Branch(f"{syst}_score_{SIG}_vs_{BKG}",
                                 self.dict_score[f"{syst}/score_{SIG}_vs_{BKG}"],
                                 f"{syst}_score_{SIG}_vs_{BKG}/F") 
            self.tree.Branch(f"{syst}_sf_muonID", self.sf_muonID[syst], f"{syst}_sf_muonID/F")
            self.tree.Branch(f"{syst}_sf_muonID_up", self.sf_muonID_up[syst], f"{syst}_sf_muonID_up/F")
            self.tree.Branch(f"{syst}_sf_muonID_down", self.sf_muonID_down[syst], f"{syst}_sf_muonID_down/F") 
            self.tree.Branch(f"{syst}_sf_dblmutrig", self.sf_dblmutrig[syst], f"{syst}_sf_dblmutrig/F") 
            self.tree.Branch(f"{syst}_sf_dblmutrig_up", self.sf_dblmutrig_up[syst], f"{syst}_sf_dblmutrig_up/F")
            self.tree.Branch(f"{syst}_sf_dblmutrig_down", self.sf_dblmutrig_down[syst], f"{syst}_sf_dblmutrig_down/F")
            self.tree.Branch(f"{syst}_sf_btag_central", self.sf_btag_central[syst], f"{syst}_sf_btag_central/F")
            self.tree.Branch(f"{syst}_sf_btag_htag_up_uncorr", self.sf_btag_htag_up_uncorr[syst], f"{syst}_sf_btag_htag_up_uncorr/F")
            self.tree.Branch(f"{syst}_sf_btag_htag_down_uncorr", self.sf_btag_htag_down_uncorr[syst], f"{syst}_sf_btag_htag_down_uncorr/F") 
            self.tree.Branch(f"{syst}_sf_btag_htag_up_corr", self.sf_btag_htag_up_corr[syst], f"{syst}_sf_btag_htag_up_corr/F")
            self.tree.Branch(f"{syst}_sf_btag_htag_down_corr", self.sf_btag_htag_down_corr[syst], f"{syst}_sf_btag_htag_down_corr/F")
            self.tree.Branch(f"{syst}_sf_btag_ltag_up_uncorr", self.sf_btag_ltag_up_uncorr[syst], f"{syst}_sf_btag_ltag_up_uncorr/F")
            self.tree.Branch(f"{syst}_sf_btag_ltag_down_uncorr", self.sf_btag_ltag_down_uncorr[syst], f"{syst}_sf_btag_ltag_down_uncorr/F") 
            self.tree.Branch(f"{syst}_sf_btag_ltag_up_corr", self.sf_btag_ltag_up_corr[syst], f"{syst}_sf_btag_ltag_up_corr/F")
            self.tree.Branch(f"{syst}_sf_btag_ltag_down_corr", self.sf_btag_ltag_down_corr[syst], f"{syst}_sf_btag_ltag_down_corr/F") 

    def __initTreeContents(self, n):
        self.w_norm[0] = -999.
        self.w_l1prefire[0] = -999.
        self.w_l1prefire_up[0] = -999.
        self.w_l1prefire_down[0] = -999.
        self.w_pileup[0] = -999.
        self.w_pileup_up[0] = -999.
        self.w_pileup_down[0] = -999.
        for syst in self.scaleVariations:
            self.dict_mass1[syst][0] = -999.
            self.dict_mass2[syst][0] = -999.
            for SIG, BKG in product(self.signalStrings, self.backgroundStrings):
                self.dict_score[f"{syst}/score_{SIG}_vs_{BKG}"][0] = -999.
            self.sf_muonID[syst][0] = -999.
            self.sf_muonID_up[syst][0] = -999.
            self.sf_muonID_down[syst][0] = -999.
            self.sf_dblmutrig[syst][0] = -999.
            self.sf_dblmutrig_up[syst][0] = -999.
            self.sf_dblmutrig_down[syst][0] = -999.
            self.sf_btag_central[syst][0] = -999.
            self.sf_btag_htag_up_uncorr[syst][0] = -999.
            self.sf_btag_htag_down_uncorr[syst][0] = -999.
            self.sf_btag_htag_up_corr[syst][0] = -999.
            self.sf_btag_htag_down_corr[syst][0] = -999.
            self.sf_btag_ltag_up_uncorr[syst][0] = -999.
            self.sf_btag_ltag_down_uncorr[syst][0] = -999.
            self.sf_btag_ltag_up_corr[syst][0] = -999.
            self.sf_btag_ltag_down_corr[syst][0] = -999.
        
    def defineObjects(self, rawMuons, rawElectrons, rawJets, syst="Central"):
        # first copy objects
        allMuons = rawMuons
        allElectrons = rawElectrons
        allJets = rawJets
        
        # check the syst argument
        if not syst in self.scaleVariations:
            print(f"[SkimPromptTree::defineObjects] Wrong scale {syst}")
            exit(1)
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
        bjets = std.vector[Jet]()
        for jet in jets:
            btagScore = jet.GetTaggerResult(3)                  # DeepJet score
            wp = super().mcCorr.GetJetTaggingCutValue(3, 1)     # DeepJet Medium
            if btagScore > wp: bjets.emplace_back(jet)
            
        vetoMuons = std.vector[Muon](sorted(vetoMuons, key=lambda x: x.Pt(), reverse=True))
        tightMuons = std.vector[Muon](sorted(tightMuons, key=lambda x: x.Pt(), reverse=True))
        vetoElectrons = std.vector[Electron](sorted(vetoElectrons, key=lambda x: x.Pt(), reverse=True))
        tightElectrons = std.vector[Electron](sorted(tightElectrons, key=lambda x: x.Pt(), reverse=True))
        jets = std.vector[Jet](sorted(jets, key=lambda x: x.Pt(), reverse=True))
        bjets = std.vector[Jet](sorted(bjets, key=lambda x: x.Pt(), reverse=True))
        
        return (vetoMuons, tightMuons, vetoElectrons, tightElectrons, jets, bjets)
    
    def selectEvent(self, event, truth, vetoMuons, tightMuons, vetoElectrons, tightElectrons, jets, bjets, METv):
        # check lepton multiplicity first
        is3Mu = (len(tightMuons) == 3 and len(vetoMuons) == 3 and \
                len(tightElectrons) == 0 and len(vetoElectrons) == 0)
        is1E2Mu = len(tightMuons) == 2 and len(vetoMuons) == 2 and \
                  len(tightElectrons) == 1 and len(vetoElectrons) == 1
        
        if self.channel == "Skim1E2Mu":
            if not is1E2Mu: return None
        if self.channel == "Skim3Mu":
            if not is3Mu: return None

        # prompt matching
        if not super().IsDATA:
            promptMuons = std.vector[Muon]()
            promptElectrons = std.vector[Electron]()
            for mu in tightMuons:
                if super().GetLeptonType(mu, truth) > 0: promptMuons.emplace_back(mu)
            for ele in tightElectrons:
                if super().GetLeptonType(ele, truth) > 0: promptElectrons.emplace_back(ele)

            if promptMuons.size() != tightMuons.size(): return None
            if promptElectrons.size() != tightElectrons.size(): return None

        # for patching samples
        if "DYJets" in super().MCSample:
            leptons = std.vector[Lepton]()
            for mu in tightMuons: leptons.emplace_back(mu)
            for ele in tightElectrons: leptons.emplace_back(ele)
            if leptons.at(0).Pt() > 20. and leptons.at(1).Pt() > 20. and leptons.at(2).Pt() > 20.:
                return None
        if "ZGToLLG" in super().MCSample:
            leptons = std.vector[Lepton]()
            for mu in tightMuons: leptons.emplace_back(mu)
            for ele in tightElectrons: leptons.emplace_back(ele)
            if leptons.at(0).Pt() < 20. or leptons.at(1).Pt() < 20. or leptons.at(2).Pt() < 20.:
                return None

        ## 1E2Mu baseline
        ## 1. pass EMuTriggers
        ## 2. Exact 2 tight muons and 1 tight electron, no additional lepton
        ## 3. Exists OS muon pair with mass > 12 GeV
        ## 4. At least two jets
        if self.channel == "Skim1E2Mu":
            if not event.PassTrigger(super().EMuTriggers): return None
            leptons = std.vector[Lepton]()
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
        ## 5. At least two jets
        elif self.channel == "Skim3Mu":
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
        else:
            return None
    
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
            print(f"[PromptEstimator::makePair] wrong no. of muons {muons.size}")
    
    #### Get scores for each event
    def evalScore(self, muons, electrons, jets, bjets, METv):
        scores = {}
        if self.network == "DenseNeuralNet": 
            data = self.getDenseInput(muons, electrons, jets, bjets, METv)
            for sig, bkg in product(self.signalStrings, self.backgroundStrings):
                scores[f"{sig}_vs_{bkg}"] = self.getDenseScore(f"{sig}_vs_{bkg}", data)
        else:                                
            data = self.getGraphInput(muons, electrons, jets, METv)
            for sig, bkg in product(self.signalStrings, self.backgroundStrings):
                scores[f"{sig}_vs_{bkg}"] = self.getGraphScore(f"{sig}_vs_{bkg}", data)
        return data, scores
    
    def getDenseInput(self, tightMuons, tightElectrons, jets, bjets, METv):
        inputs = []
        dRjets = []
        for idx1, idx2 in combinations(range(jets.size()), 2):
            dRjets.append(jets[idx1].DeltaR(jets[idx2]))
        dRjets = sum(dRjets) / len(dRjets)

        if self.channel == "Skim1E2Mu":
            mu1 = tightMuons.at(0)
            mu2 = tightMuons.at(1)
            ele = tightElectrons.at(0)
            j1 = jets.at(0)
            j2 = jets.at(1)
            MT = (ele+METv).Mt()
            
            inputs.append(mu1.E())
            inputs.append(mu1.Px())
            inputs.append(mu1.Py())
            inputs.append(mu1.Pz())
            inputs.append(mu1.Charge())
            inputs.append(mu2.E())
            inputs.append(mu2.Px())
            inputs.append(mu2.Py())
            inputs.append(mu2.Pz())
            inputs.append(mu2.Charge())
            inputs.append(ele.E())
            inputs.append(ele.Px())
            inputs.append(ele.Py())
            inputs.append(ele.Pz())
            inputs.append(ele.Charge())
            inputs.append(j1.E())
            inputs.append(j1.Px())
            inputs.append(j1.Py())
            inputs.append(j1.Pz())
            inputs.append(j1.Charge())
            inputs.append(j1.GetTaggerResult(3))
            inputs.append(j2.E())
            inputs.append(j2.Px())
            inputs.append(j2.Py())
            inputs.append(j2.Pz())
            inputs.append(j2.Charge())
            inputs.append(j2.GetTaggerResult(3))
            inputs.append(mu1.DeltaR(mu2))
            inputs.append(mu1.DeltaR(ele))
            inputs.append(mu2.DeltaR(ele))
            inputs.append(j1.DeltaR(j2))
            inputs.append(sum([j.Pt() for j in jets]))
            inputs.append(sum([l.Pt() for l in tightElectrons+tightMuons]))
            inputs.append(MT)
            inputs.append(METv.Pt())
            inputs.append(jets.size())
            inputs.append(bjets.size())
            inputs.append(dRjets)
            inputs.append(sum([j.GetTaggerResult(3) for j in jets]) / jets.size())
        elif self.channel == "Skim3Mu":
            mu1 = tightMuons.at(0)
            mu2 = tightMuons.at(1)
            mu3 = tightMuons.at(2)
            j1 = jets.at(0)
            j2 = jets.at(1)
            MT1 = (mu1+METv).Mt()
            MT2 = (mu2+METv).Mt()
            MT3 = (mu3+METv).Mt()

            inputs.append(mu1.E())
            inputs.append(mu1.Px())
            inputs.append(mu1.Py())
            inputs.append(mu1.Pz())
            inputs.append(mu1.Charge())
            inputs.append(mu2.E())
            inputs.append(mu2.Px())
            inputs.append(mu2.Py())
            inputs.append(mu2.Pz())
            inputs.append(mu2.Charge())
            inputs.append(mu3.E())
            inputs.append(mu3.Px())
            inputs.append(mu3.Py())
            inputs.append(mu3.Pz())
            inputs.append(mu3.Charge())
            inputs.append(j1.E())
            inputs.append(j1.Px())
            inputs.append(j1.Py())
            inputs.append(j1.Pz())
            inputs.append(j1.Charge())
            inputs.append(j1.GetTaggerResult(3))
            inputs.append(j2.E())
            inputs.append(j2.Px())
            inputs.append(j2.Py())
            inputs.append(j2.Pz())
            inputs.append(j2.Charge())
            inputs.append(j2.GetTaggerResult(3))
            inputs.append(mu1.DeltaR(mu2))
            inputs.append(mu1.DeltaR(mu3))
            inputs.append(mu2.DeltaR(mu3))
            inputs.append(j1.DeltaR(mu1))
            inputs.append(j1.DeltaR(mu2))
            inputs.append(j1.DeltaR(mu3))
            inputs.append(j2.DeltaR(mu1))
            inputs.append(j2.DeltaR(mu2))
            inputs.append(j2.DeltaR(mu3))
            inputs.append(j1.DeltaR(j2))
            inputs.append(sum([j.Pt() for j in jets]))
            inputs.append(sum([l.Pt() for l in tightMuons]))
            inputs.append(MT1)
            inputs.append(MT2)
            inputs.append(MT3)
            inputs.append(METv.Pt())
            inputs.append(jets.size())
            inputs.append(bjets.size())
            inputs.append(dRjets)
            inputs.append(sum([j.GetTaggerResult(3) for j in jets])/jets.size())
        return torch.FloatTensor([inputs])
        
    def getDenseScore(self, modelKey, data):
        with torch.no_grad():
            out = self.models[modelKey](data)
        return out.numpy()[0][1]

    def getGraphInput(self, tightMuons, tightElectrons, jets, METv):
        particles = []
        for muon in tightMuons:
            node = NodeParticle()
            node.isMuon = True
            node.SetPtEtaPhiM(muon.Pt(), muon.Eta(), muon.Phi(), muon.M())
            node.charge = muon.Charge()
            particles.append(node)
        for ele in tightElectrons:
            node = NodeParticle()
            node.isElectron = True
            node.SetPtEtaPhiM(ele.Pt(), ele.Eta(), ele.Phi(), ele.M())
            node.charge = ele.Charge()
            particles.append(node)
        for jet in jets:
            node = NodeParticle()
            node.isJet = True
            node.SetPtEtaPhiM(jet.Pt(), jet.Eta(), jet.Phi(), jet.M())
            node.charge = jet.Charge()
            node.btagScore = jet.GetTaggerResult(3)
            particles.append(node)
        missing = NodeParticle()
        missing.SetPtEtaPhiM(METv.Pt(), 0., METv.Phi(), 0.)
        particles.append(missing)
        nodeList = []
        for particle in particles:
            nodeList.append([particle.E(), particle.Px(), particle.Py(), particle.Pz(),
                             particle.Charge(), particle.BtagScore(), 
                             particle.IsMuon(), particle.IsElectron(), particle.IsJet()])
        data = evtToGraph(nodeList, y=None, k=4)
        return data

    def getGraphScore(self, modelKey, data):
        with torch.no_grad():
            out = self.models[modelKey](data.x, data.edge_index)
        return out.numpy()[0][1]
    
    def WriteHist(self):
        super().outfile.cd()
        self.tree.Write()
