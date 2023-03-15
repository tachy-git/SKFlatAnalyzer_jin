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


class MatrixUnbinned(TriLeptonBase):
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
        
        self.systematics = ["Central",
                            "NonpromptUp",
                            "NonpromptDown"]
        
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
        
        # initialize contents
        self.__initTreeContents()
        
        # fill contents
        vetoMuons, looseMuons, tightMuons, vetoElectrons, looseElectrons, tightElectrons, jets, bjets = self.defineObjects(rawMuons, rawElectrons, rawJets)
        channel = self.selectEvent(ev, vetoMuons, looseMuons, tightMuons, vetoElectrons, looseElectrons, tightElectrons, jets, bjets, METv)

        if channel is None: return None

        pairs = self.makePair(looseMuons)
        data, scores = self.evalScore(looseMuons, looseElectrons, jets, bjets, METv)
            
        if channel == "SR1E2Mu": 
            self.mass1[0] = pairs.M()
            self.mass2[0] = 0
        else:
            self.mass1[0] = pairs[0].M()
            self.mass2[0] = pairs[1].M()
        
        for SIG, BKG in product(self.signalStrings, self.backgroundStrings):
            self.dict_score[f"score_{SIG}_vs_{BKG}"][0] = scores[f"{SIG}_vs_{BKG}"]
        
        self.w_nonprompt[0] = super().getFakeWeight(looseMuons, looseElectrons, 0)
        self.w_nonprompt_up[0] = super().getFakeWeight(looseMuons, looseElectrons, 1)
        self.w_nonprompt_down[0] = super().getFakeWeight(looseMuons, looseElectrons, -1)
        self.tree.Fill()

    def __loadModels(self):
        self.models = {}

        for sig, bkg in product(self.signalStrings, self.backgroundStrings):
            csv = pd.read_csv(f"{os.environ['DATA_DIR']}/FullRun2/{self.network}/{self.channel}/results/summary_{sig}_vs_{bkg}.txt", 
                              sep=",\s", 
                              engine="python", 
                              header=None).transpose()
            modelPath = f"{os.environ['DATA_DIR']}/FullRun2/{self.network}/{self.channel}/models/{sig}_vs_{bkg}.pt"
            if self.network == "DenseNeuralNet":
                modelArch = csv[0][3]
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
        self.mass1 = array("f", [0.])
        self.mass2 = array("f", [0.])
        self.dict_score = {}
        self.w_nonprompt = array("f", [0.])
        self.w_nonprompt_up = array("f", [0.])
        self.w_nonprompt_down = array("f", [0.])
 
        for SIG, BKG in product(self.signalStrings, self.backgroundStrings):
            self.dict_score[f"score_{SIG}_vs_{BKG}"] = array("f", [0.])
        
        self.tree.Branch("Central_mass1", self.mass1, "mass1/F")
        self.tree.Branch("Central_mass2", self.mass2, "mass2/F")
        for SIG, BKG in product(self.signalStrings, self.backgroundStrings): 
            self.tree.Branch(f"Central_score_{SIG}_vs_{BKG}",
                             self.dict_score[f"score_{SIG}_vs_{BKG}"],
                             f"Central_score_{SIG}_vs_{BKG}/F") 
        self.tree.Branch("w_nonprompt", self.w_nonprompt, "w_nonprompt/F")
        self.tree.Branch("w_nonprompt_up", self.w_nonprompt_up, "w_nonprompt_up/F")
        self.tree.Branch("w_nonprompt_down", self.w_nonprompt_down, "w_nonprompt_down/F")

    def __initTreeContents(self):
        self.mass1[0] = -999.
        self.mass2[0] = -999.
        self.w_nonprompt[0] = -999.
        self.w_nonprompt_up[0] = -999.
        self.w_nonprompt_down[0] = -999.
        for SIG, BKG in product(self.signalStrings, self.backgroundStrings):
            self.dict_score[f"score_{SIG}_vs_{BKG}"][0] = -999.

    def defineObjects(self, rawMuons, rawElectrons, rawJets, syst="Central"):
        # first copy objects
        allMuons = rawMuons
        allElectrons = rawElectrons
        allJets = rawJets

        vetoMuons = super().SelectMuons(allMuons, super().MuonIDs[2], 10., 2.4)
        looseMuons = super().SelectMuons(vetoMuons, super().MuonIDs[1], 10., 2.4)
        tightMuons = super().SelectMuons(looseMuons, super().MuonIDs[0], 10., 2.4)
        vetoElectrons = super().SelectElectrons(allElectrons, super().ElectronIDs[2], 10., 2.5)
        looseElectrons = super().SelectElectrons(vetoElectrons, super().ElectronIDs[1], 10., 2.5)
        tightElectrons = super().SelectElectrons(looseElectrons, super().ElectronIDs[0], 10., 2.5)
        jets = super().SelectJets(allJets, "tight", 20., 2.4)
        jets = super().JetsVetoLeptonInside(jets, vetoElectrons, vetoMuons, 0.4)
        bjets = std.vector[Jet]()
        for jet in jets:
            btagScore = jet.GetTaggerResult(3)                  # DeepJet score
            wp = super().mcCorr.GetJetTaggingCutValue(3, 1)     # DeepJet Medium
            if btagScore > wp: bjets.emplace_back(jet)

        vetoMuons = std.vector[Muon](sorted(vetoMuons, key=lambda x: x.Pt(), reverse=True))
        looseMuons = std.vector[Muon](sorted(looseMuons, key=lambda x: x.Pt(), reverse=True))
        tightMuons = std.vector[Muon](sorted(tightMuons, key=lambda x: x.Pt(), reverse=True))
        vetoElectrons = std.vector[Electron](sorted(vetoElectrons, key=lambda x: x.Pt(), reverse=True))
        looseElectrons = std.vector[Electron](sorted(looseElectrons, key=lambda x: x.Pt(), reverse=True))
        tightElectrons = std.vector[Electron](sorted(tightElectrons, key=lambda x: x.Pt(), reverse=True))
        jets = std.vector[Jet](sorted(jets, key=lambda x: x.Pt(), reverse=True))
        bjets = std.vector[Jet](sorted(bjets, key=lambda x: x.Pt(), reverse=True))

        return (vetoMuons, looseMuons, tightMuons, vetoElectrons, looseElectrons, tightElectrons, jets, bjets)

    def selectEvent(self, event, vetoMuons, looseMuons, tightMuons, vetoElectrons, looseElectrons, tightElectrons, jets, bjets, METv):
        is3Mu = (looseMuons.size() == 3 and vetoMuons.size() == 3 and \
                 looseElectrons.size() == 0 and vetoElectrons.size() == 0)
        is1E2Mu = (looseMuons.size() == 2 and vetoMuons.size() == 2 and \
                   looseElectrons.size() == 1 and vetoElectrons.size() == 1)

        #### not all leptons tight
        if self.channel == "Skim1E2Mu":
            if not is1E2Mu: return None
            if tightMuons.size() == looseMuons.size(): return None
            if tightElectrons.size() == looseElectrons.size(): return None

        if self.channel == "Skim3Mu":
            if not is3Mu: return None
            if tightMuons.size() == looseMuons.size(): return None

        channel = ""
        ## 1E2Mu baseline
        ## 1. pass EMuTriggers
        ## 2. Exact 2 tight muons and 1 tight electron, no additional lepton
        ## 3. Exists OS muon pair with mass > 12 GeV
        ## 4. At least two jets
        if self.channel == "Skim1E2Mu":
            if not event.PassTrigger(super().EMuTriggers): return None
            leptons = std.vector[Lepton]()
            for mu in looseMuons: leptons.emplace_back(mu)
            for ele in looseElectrons: leptons.emplace_back(ele)
            mu1, mu2, ele = tuple(leptons)
            passLeadMu = mu1.Pt() > 25. and ele.Pt() > 15.
            passLeadEle = mu1.Pt() > 10. and ele.Pt() > 25.
            passSafeCut = passLeadMu or passLeadEle
            if not passSafeCut: return None
            if not mu1.Charge()+mu2.Charge() == 0: return None
            pair = self.makePair(looseMuons)
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
            mu1, mu2, mu3  = tuple(looseMuons)
            if not mu1.Pt() > 20.: return None
            if not mu2.Pt() > 10.: return None
            if not mu3.Pt() > 10.: return None
            if not abs(mu1.Charge()+mu2.Charge()+mu3.Charge()) == 1: return None
            pair1, pair2 = self.makePair(looseMuons)
            if not pair1.M() > 12.: return None
            if not pair2.M() > 12.: return None
            if not jets.size() >= 2: return None
            # orthogonality of SR and CR done by bjet multiplicity
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
