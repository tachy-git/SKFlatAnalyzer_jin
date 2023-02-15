from ROOT import gSystem
from ROOT import TriLeptonBase
from ROOT import std, TString
from ROOT.JetTagging import Parameters as jParameters
from ROOT import Lepton, Muon, Electron, Jet
gSystem.Load("/cvmfs/cms.cern.ch/slc7_amd64_gcc900/external/lhapdf/6.2.3/lib/libLHAPDF.so")

import os
import numpy as np
import pandas as pd
import torch

from itertools import product, combinations
from MLTools.models import SNN, SNNLite, ParticleNet, ParticleNetLite
from MLTools.helpers import evtToGraph
from MLTools.formats import NodeParticle

class NonpromptEstimator(TriLeptonBase):
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
            
        self.systematics = ["Central", "NonpromptUp", "NonpromptDown"]     
        self.signalStrings = [
                "MHc-70_MA-15", "MHc-70_MA-40", "MHc-70_MA-65",
                "MHc-100_MA-15", "MHc-100_MA-60", "MHc-100_MA-95",
                "MHc-130_MA-15", "MHc-130_MA-55", "MHc-130_MA-90", "MHc-130_MA-125",
                "MHc-160_MA-15", "MHc-160_MA-85", "MHc-160_MA-120", "MHc-160_MA-155"]
        self.backgroundStrings = ["TTLL_powheg", "ttX"]
    
        self.__loadModels()
        
    def __loadModels(self):
        self.models = {}

        for sig, bkg in product(self.signalStrings, self.backgroundStrings):
            csv = pd.read_csv(f"{os.environ['DATA_DIR']}/FullRun2/{self.network}/{self.channel}/summary/summary_{sig}_vs_{bkg}.txt",
                              sep=",\s",
                              engine="python",
                              header=None).transpose()
            modelArch, dropout_p, readout = csv[0][3:6]
            modelPath = f"{os.environ['DATA_DIR']}/FullRun2/{self.network}/{self.channel}/models/{sig}_vs_{bkg}.pt"
            if self.network == "DenseNeuralNet":
                if self.channel == "Skim1E2Mu":
                    if modelArch == "SNN": model = SNN(41, 2)
                    else:                  model = SNNLite(41, 2)
                if self.channel == "Skim3Mu" :
                    if modelArch == "SNN": model = SNN(47, 2)
                    else:                  model = SNNLite(47, 2)
            else:               # GraphNeuralNet
                if modelArch == "ParticleNet": model = ParticleNet(9, 2, dropout_p, readout)
                else:                          model = ParticleNetLite(9, 2, dropout_p, readout)
            model.load_state_dict(torch.load(modelPath, map_location=torch.device("cpu")))
            model.eval()
            self.models[f"{sig}_vs_{bkg}"] = model

    def executeEvent(self):
        if not super().PassMETFilter(): return None
        ev = super().GetEvent()
        rawMuons = super().GetAllMuons()
        rawElectrons = super().GetAllElectrons()
        rawJets = super().GetAllJets()
        METv = ev.GetMETVector()
        
        vetoMuons, looseMuons, tightMuons, vetoElectrons, looseElectrons, tightElectrons, jets, bjets = self.defineObjects(rawMuons, rawElectrons, rawJets)
        channel = self.selectEvent(ev, vetoMuons, looseMuons, tightMuons, vetoElectrons, looseElectrons, tightElectrons, jets, bjets, METv)
        
        if channel is None: return None
        
        pairs = self.makePair(looseMuons)
        data, scores = self.evalScore(looseMuons, looseElectrons, jets, bjets, METv)
        objects = {"muons": looseMuons,
                   "electrons": looseElectrons,
                   "jets": jets,
                   "bjets": bjets,
                   "METv": METv,
                   "pairs": pairs,
                   "data": data,
                   "scores": scores
                    }
        for syst in self.systematics:
            if syst == "Central":
                weight = super().getFakeWeight(looseMuons, looseElectrons, 0)
            elif syst == "NonpromptUp":
                weight = super().getFakeWeight(looseMuons, looseElectrons, 1)
            elif syst == "NonpromptDown":
                weight = super().getFakeWeight(looseMuons, looseElectrons, -1)
            else:
                print(f"[NonpromptEstimator::executeEvent] Wrong systematics {syst}")
                exit(1)
            self.FillObjects(channel, objects, weight, syst)
        
        
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

            # orthogonality of SR and CR done by bjet multiplicity
            if bjets.size() >= 1:
                channel == "SR1E2Mu"
            else:
                mZ = 91.2
                isOnZ = abs(pair.M() - mZ) < 10.
                if isOnZ: channel = "ZFake1E2Mu"
                else:
                    if abs((mu1+mu2+ele).M() - mZ) < 10.: channel = "ZGamma1E2Mu"
                    else: return None

        ## 3Mu baseline
        ## 1. pass DblMuTriggers
        ## 2. Exact 3 tight muons, no additional leptons
        ## 3. Exist OS muon pair,
        ## 4. All OS muon pair mass > 12 GeV
        ## 5. At least two jets
        else:
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
            if bjets.size() >= 1:
                channel = "SR3Mu"
            else:
                mZ = 91.2
                isOnZ = abs(pair1.M() - mZ) < 10. or abs(pair2.M() - mZ) < 10.
                if isOnZ: channel = "ZFake3Mu"
                else:
                    if abs((mu1+mu2+mu3).M() - mZ) < 10.: channel = "ZGamma3Mu"
                    else: return None
        return channel
        
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
    
    def FillObjects(self, channel, objects, weight, syst):
        muons = objects["muons"]
        electrons = objects["electrons"]
        jets = objects["jets"]
        bjets = objects["bjets"]
        METv = objects["METv"]
        pairs = objects["pairs"]
        data = objects["data"]
        scores = objects["scores"]
        
        ## fill input observables
        for idx, mu in enumerate(muons, start=1):
            super().FillHist(f"{channel}/{syst}/muons/{idx}/pt", mu.Pt(), weight, 300, 0., 300.)
            super().FillHist(f"{channel}/{syst}/muons/{idx}/eta", mu.Eta(), weight, 48, -2.4, 2.4)
            super().FillHist(f"{channel}/{syst}/muons/{idx}/phi", mu.Phi(), weight, 64, -3.2, 3.2)
            super().FillHist(f"{channel}/{syst}/muons/{idx}/mass", mu.M(), weight, 10, 0., 1.)
        for idx, ele in enumerate(electrons, start=1):
            super().FillHist(f"{channel}/{syst}/electrons/{idx}/pt", ele.Pt(), weight, 300, 0., 300.)
            super().FillHist(f"{channel}/{syst}/electrons/{idx}/eta", ele.Eta(), weight, 50, -2.5, 2.5)
            super().FillHist(f"{channel}/{syst}/electrons/{idx}/Phi", ele.Phi(), weight, 64, -3.2, 3.2)
            super().FillHist(f"{channel}/{syst}/electrons/{idx}/mass", ele.M(), weight, 100, 0., 1.)
        for idx, jet in enumerate(jets, start=1):
            super().FillHist(f"{channel}/{syst}/jets/{idx}/pt", jet.Pt(), weight, 300, 0., 300.)
            super().FillHist(f"{channel}/{syst}/jets/{idx}/eta", jet.Eta(), weight, 48, -2.4, 2.4)
            super().FillHist(f"{channel}/{syst}/jets/{idx}/phi", jet.Phi(), weight, 64, -3.2, 3.2)
            super().FillHist(f"{channel}/{syst}/jets/{idx}/mass", jet.M(), weight, 100, 0., 100.)
        for idx, bjet in enumerate(bjets, start=1):
            super().FillHist(f"{channel}/{syst}/bjets/{idx}/pt", bjet.Pt(), weight, 300, 0., 300.)
            super().FillHist(f"{channel}/{syst}/bjets/{idx}/eta", bjet.Eta(), weight, 48, -2.4, 2.4)
            super().FillHist(f"{channel}/{syst}/bjets/{idx}/phi", bjet.Phi(), weight, 64, -3.2, 3.2)
            super().FillHist(f"{channel}/{syst}/bjets/{idx}/mass", bjet.M(), weight, 100, 0., 100.)
        super().FillHist(f"{channel}/{syst}/jets/size", jets.size(), weight, 20, 0., 20.)
        super().FillHist(f"{channel}/{syst}/bjets/size", bjets.size(), weight, 15, 0., 15.)
        super().FillHist(f"{channel}/{syst}/METv/pt", METv.Pt(), weight, 300, 0., 300.)
        super().FillHist(f"{channel}/{syst}/METv/phi", METv.Phi(), weight, 64, -3.2, 3.2)

        # Fill ZCands
        if "1E2Mu" in channel:
            if not "ZGamma" in channel: ZCand = pairs       # pairs == pair
            else:                       ZCand = muons.at(0) + muons.at(1) + electrons.at(0) 
            super().FillHist(f"{channel}/{syst}/ZCand/pt", ZCand.Pt(), weight, 300, 0., 300.)
            super().FillHist(f"{channel}/{syst}/ZCand/eta", ZCand.Eta(), weight, 100, -5., 5.)
            super().FillHist(f"{channel}/{syst}/ZCand/phi", ZCand.Phi(), weight, 64, -3.2, 3.2)
            super().FillHist(f"{channel}/{syst}/ZCand/mass", ZCand.M(), weight, 200, 0., 200.)
        else:
            if not "ZGamma" in channel:
                pair1, pair2 = pairs
                mZ = 91.2
                if abs(pair1.M() - mZ) < abs(pair2.M() - mZ): ZCand, nZCand = pair1, pair2
                else:                                         ZCand, nZCand = pair2, pair1
            else:
                ZCand = muons.at(0) + muons.at(1) + muons.at(2)
                nZCand = NodeParticle()
            super().FillHist(f"{channel}/{syst}/ZCand/pt", ZCand.Pt(), weight, 300, 0., 300.)
            super().FillHist(f"{channel}/{syst}/ZCand/eta", ZCand.Eta(), weight, 100, -5., 5.)
            super().FillHist(f"{channel}/{syst}/ZCand/phi", ZCand.Phi(), weight, 64, -3.2, 3.2)
            super().FillHist(f"{channel}/{syst}/ZCand/mass", ZCand.M(), weight, 200, 0., 200.)
            super().FillHist(f"{channel}/{syst}/nZCand/pt", nZCand.Pt(), weight, 300, 0., 300.)
            super().FillHist(f"{channel}/{syst}/nZCand/eta", nZCand.Eta(), weight, 100, -5., 5.)
            super().FillHist(f"{channel}/{syst}/nZCand/phi", nZCand.Phi(), weight, 64, -3.2, 3.2)
            super().FillHist(f"{channel}/{syst}/nZCand/mass", nZCand.M(), weight, 200, 0., 200.)
         
        # Fill inputs for the network
        if self.network == "DenseNeuralNet":
            data = data[0].numpy()
            if self.channel == "Skim1E2Mu":
                super().FillHist(f"{channel}/{syst}/inputs/mu1_energy", data[0], weight, 300, 0., 300.)
                super().FillHist(f"{channel}/{syst}/inputs/mu1_px", data[1], weight, 500, -250., 250.)
                super().FillHist(f"{channel}/{syst}/inputs/mu1_py", data[2], weight, 500, -250., 250.)
                super().FillHist(f"{channel}/{syst}/inputs/mu1_pz", data[3], weight, 500, -250., 250.)
                super().FillHist(f"{channel}/{syst}/inputs/mu1_charge", data[4], weight, 3, -1, 2)
                super().FillHist(f"{channel}/{syst}/inputs/mu2_energy", data[5], weight, 300, 0., 300.)
                super().FillHist(f"{channel}/{syst}/inputs/mu2_px", data[6], weight, 500, -250., 250.)
                super().FillHist(f"{channel}/{syst}/inputs/mu2_py", data[7], weight, 500, -250., 250.)
                super().FillHist(f"{channel}/{syst}/inputs/mu2_pz", data[8], weight, 500, -250., 250.)
                super().FillHist(f"{channel}/{syst}/inputs/mu2_charge", data[9], weight, 3, -1, 2) 
                super().FillHist(f"{channel}/{syst}/inputs/ele_energy", data[10], weight, 300, 0., 300.)
                super().FillHist(f"{channel}/{syst}/inputs/ele_px", data[11], weight, 500, -250., 250.)
                super().FillHist(f"{channel}/{syst}/inputs/ele_py", data[12], weight, 500, -250., 250.)
                super().FillHist(f"{channel}/{syst}/inputs/ele_pz", data[13], weight, 500, -250., 250.)
                super().FillHist(f"{channel}/{syst}/inputs/ele_charge", data[14], weight, 3, -1, 2)
                super().FillHist(f"{channel}/{syst}/inputs/j1_energy", data[15], weight, 300, 0., 300.)
                super().FillHist(f"{channel}/{syst}/inputs/j1_px", data[16], weight, 500, -250., 250.)
                super().FillHist(f"{channel}/{syst}/inputs/j1_py", data[17], weight, 500, -250., 250.)
                super().FillHist(f"{channel}/{syst}/inputs/j1_pz", data[18], weight, 500, -250., 250.)
                super().FillHist(f"{channel}/{syst}/inputs/j1_charge", data[19], weight, 200, -1, 1)
                super().FillHist(f"{channel}/{syst}/inputs/j1_btagScore", data[20], weight, 100, 0., 1)
                super().FillHist(f"{channel}/{syst}/inputs/j2_energy", data[21], weight, 300, 0., 300.)
                super().FillHist(f"{channel}/{syst}/inputs/j2_px", data[22], weight, 500, -250., 250.)
                super().FillHist(f"{channel}/{syst}/inputs/j2_py", data[23], weight, 500, -250., 250.)
                super().FillHist(f"{channel}/{syst}/inputs/j2_pz", data[24], weight, 500, -250., 250.)
                super().FillHist(f"{channel}/{syst}/inputs/j2_charge", data[25], weight, 200, -1, 1)
                super().FillHist(f"{channel}/{syst}/inputs/j2_btagScore", data[26], weight, 100, 0., 1)
                super().FillHist(f"{channel}/{syst}/inputs/dR_mu1mu2", data[27], weight, 100, 0., 5.)
                super().FillHist(f"{channel}/{syst}/inputs/dR_mu1ele", data[28], weight, 100, 0., 5.)
                super().FillHist(f"{channel}/{syst}/inputs/dR_mu2ele", data[29], weight, 100, 0., 5.)
                super().FillHist(f"{channel}/{syst}/inputs/dR_j1ele", data[30], weight, 100, 0., 5.)
                super().FillHist(f"{channel}/{syst}/inputs/dR_j2ele", data[31], weight, 100, 0., 5.)
                super().FillHist(f"{channel}/{syst}/inputs/dR_j1j2", data[32], weight, 100, 0., 5.) 
                super().FillHist(f"{channel}/{syst}/inputs/HT", data[33], weight, 1000, 0., 1000.)
                super().FillHist(f"{channel}/{syst}/inputs/LT", data[34], weight, 800, 0., 800.)
                super().FillHist(f"{channel}/{syst}/inputs/MT", data[35], weight, 300, 0., 300.)
                super().FillHist(f"{channel}/{syst}/inputs/MET", data[36], weight, 300, 0., 300.)
                super().FillHist(f"{channel}/{syst}/inputs/Nj", data[37], weight, 20, 0., 20.)
                super().FillHist(f"{channel}/{syst}/inputs/Nb", data[38], weight, 20, 0., 20.) 
                super().FillHist(f"{channel}/{syst}/inputs/avg_dRjets", data[39], weight, 100, 0., 5.)
                super().FillHist(f"{channel}/{syst}/inputs/avg_btagScore", data[40], weight, 100, 0., 1.)
            else:       # Skim3Mu
                super().FillHist(f"{channel}/{syst}/inputs/mu1_energy", data[0], weight, 300, 0., 300.)
                super().FillHist(f"{channel}/{syst}/inputs/mu1_px", data[1], weight, 500, -250., 250.)
                super().FillHist(f"{channel}/{syst}/inputs/mu1_py", data[2], weight, 500, -250., 250.)
                super().FillHist(f"{channel}/{syst}/inputs/mu1_pz", data[3], weight, 500, -250., 250.)
                super().FillHist(f"{channel}/{syst}/inputs/mu1_charge", data[4], weight, 3, -1, 2)
                super().FillHist(f"{channel}/{syst}/inputs/mu2_energy", data[5], weight, 300, 0., 300.)
                super().FillHist(f"{channel}/{syst}/inputs/mu2_px", data[6], weight, 500, -250., 250.)
                super().FillHist(f"{channel}/{syst}/inputs/mu2_py", data[7], weight, 500, -250., 250.)
                super().FillHist(f"{channel}/{syst}/inputs/mu2_pz", data[8], weight, 500, -250., 250.)
                super().FillHist(f"{channel}/{syst}/inputs/mu2_charge", data[9], weight, 3, -1, 2)
                super().FillHist(f"{channel}/{syst}/inputs/mu3_energy", data[10], weight, 300, 0., 300.)
                super().FillHist(f"{channel}/{syst}/inputs/mu3_px", data[11], weight, 500, -250., 250.)
                super().FillHist(f"{channel}/{syst}/inputs/mu3_py", data[12], weight, 500, -250., 250.)
                super().FillHist(f"{channel}/{syst}/inputs/mu3_pz", data[13], weight, 500, -250., 250.)
                super().FillHist(f"{channel}/{syst}/inputs/mu3_charge", data[14], weight, 3, -1, 2)
                super().FillHist(f"{channel}/{syst}/inputs/j1_energy", data[15], weight, 300, 0., 300.)
                super().FillHist(f"{channel}/{syst}/inputs/j1_px", data[16], weight, 500, -250., 250.)
                super().FillHist(f"{channel}/{syst}/inputs/j1_py", data[17], weight, 500, -250., 250.)
                super().FillHist(f"{channel}/{syst}/inputs/j1_pz", data[18], weight, 500, -250., 250.)
                super().FillHist(f"{channel}/{syst}/inputs/j1_charge", data[19], weight, 200, -1, 1)
                super().FillHist(f"{channel}/{syst}/inputs/j1_btagScore", data[20], weight, 100, 0., 1)
                super().FillHist(f"{channel}/{syst}/inputs/j2_energy", data[21], weight, 300, 0., 300.)
                super().FillHist(f"{channel}/{syst}/inputs/j2_px", data[22], weight, 500, -250., 250.)
                super().FillHist(f"{channel}/{syst}/inputs/j2_py", data[23], weight, 500, -250., 250.)
                super().FillHist(f"{channel}/{syst}/inputs/j2_pz", data[24], weight, 500, -250., 250.)
                super().FillHist(f"{channel}/{syst}/inputs/j2_charge", data[25], weight, 200, -1, 1)
                super().FillHist(f"{channel}/{syst}/inputs/j2_btagScore", data[26], weight, 100, 0., 1)
                super().FillHist(f"{channel}/{syst}/inputs/dR_mu1mu2", data[27], weight, 100, 0., 5.)
                super().FillHist(f"{channel}/{syst}/inputs/dR_mu1mu3", data[28], weight, 100, 0., 5.)
                super().FillHist(f"{channel}/{syst}/inputs/dR_mu2mu3", data[29], weight, 100, 0., 5.)
                super().FillHist(f"{channel}/{syst}/inputs/dR_j1mu1", data[30], weight, 100, 0., 5.)
                super().FillHist(f"{channel}/{syst}/inputs/dR_j1mu2", data[31], weight, 100, 0., 5.)
                super().FillHist(f"{channel}/{syst}/inputs/dR_j1mu3", data[32], weight, 100, 0., 5.)  
                super().FillHist(f"{channel}/{syst}/inputs/dR_j2mu1", data[33], weight, 100, 0., 5.)
                super().FillHist(f"{channel}/{syst}/inputs/dR_j2mu2", data[34], weight, 100, 0., 5.)
                super().FillHist(f"{channel}/{syst}/inputs/dR_j2mu3", data[35], weight, 100, 0., 5.)
                super().FillHist(f"{channel}/{syst}/inputs/dR_j1j2", data[36], weight, 100, 0., 5.)
                super().FillHist(f"{channel}/{syst}/inputs/HT", data[37], weight, 1000, 0., 1000.)
                super().FillHist(f"{channel}/{syst}/inputs/LT", data[38], weight, 800, 0., 800.)
                super().FillHist(f"{channel}/{syst}/inputs/MT1", data[39], weight, 300, 0., 300.)
                super().FillHist(f"{channel}/{syst}/inputs/MT2", data[40], weight, 300, 0., 300.)
                super().FillHist(f"{channel}/{syst}/inputs/MT3", data[41], weight, 300, 0., 300.) 
                super().FillHist(f"{channel}/{syst}/inputs/MET", data[42], weight, 300, 0., 300.)
                super().FillHist(f"{channel}/{syst}/inputs/Nj", data[43], weight, 20, 0., 20.)
                super().FillHist(f"{channel}/{syst}/inputs/Nb", data[44], weight, 20, 0., 20.) 
                super().FillHist(f"{channel}/{syst}/inputs/avg_dRjets", data[45], weight, 100, 0., 5.)
                super().FillHist(f"{channel}/{syst}/inputs/avg_btagScore", data[46], weight, 100, 0., 1.)
        else:       # network == GraphNeuralNet
            for idx, mu in enumerate(muons, start=1):
                super().FillHist(f"{channel}/{syst}/inputs/muons/{idx}/energy", mu.E(), weight, 300, 0., 300.)
                super().FillHist(f"{channel}/{syst}/inputs/muons/{idx}/px", mu.Px(), weight, 500, -250., 250.)
                super().FillHist(f"{channel}/{syst}/inputs/muons/{idx}/py", mu.Py(), weight, 500, -250., 250.)
                super().FillHist(f"{channel}/{syst}/inputs/muons/{idx}/pz", mu.Pz(), weight, 500, -250., 250.)
                super().FillHist(f"{channel}/{syst}/inputs/muons/{idx}/charge", mu.Charge(), weight, 3, -1, 2) 
            for idx, ele in enumerate(electrons, start=1):
                super().FillHist(f"{channel}/{syst}/inputs/electrons/{idx}/energy", ele.E(), weight, 300, 0., 300.)
                super().FillHist(f"{channel}/{syst}/inputs/electrons/{idx}/px", ele.Px(), weight, 500, -250., 250.)
                super().FillHist(f"{channel}/{syst}/inputs/electrons/{idx}/py", ele.Py(), weight, 500, -250., 250.)
                super().FillHist(f"{channel}/{syst}/inputs/electrons/{idx}/pz", ele.Pz(), weight, 500, -250., 250.)
                super().FillHist(f"{channel}/{syst}/inputs/electrons/{idx}/charge", ele.Charge(), weight, 3, -1, 2) 
            for idx, jet in enumerate(jets, start=1):
                super().FillHist(f"{channel}/{syst}/inputs/jets/{idx}/energy", jet.E(), weight, 300, 0., 300.)
                super().FillHist(f"{channel}/{syst}/inputs/jets/{idx}/px", jet.Px(), weight, 500, -250., 250.)
                super().FillHist(f"{channel}/{syst}/inputs/jets/{idx}/py", jet.Py(), weight, 500, -250., 250.)
                super().FillHist(f"{channel}/{syst}/inputs/jets/{idx}/pz", jet.Pz(), weight, 500, -250., 250.)
                super().FillHist(f"{channel}/{syst}/inputs/jets/{idx}/charge", jet.Charge(), weight, 200, -1, 1)
                super().FillHist(f"{channel}/{syst}/inputs/jets/{idx}/btagScore", jet.GetTaggerResult(3), weight, 100, 0., 1.)
            super().FillHist(f"{channel}/{syst}/inputs/METv/energay", METv.E(), weight, 300, 0., 300.)
            super().FillHist(f"{channel}/{syst}/inputs/METv/px", METv.Px(), weight, 500, -250, 250)
            super().FillHist(f"{channel}/{syst}/inputs/METv/py", METv.Py(), weight, 500, -250, 250)
            super().FillHist(f"{channel}/{syst}/inputs/METv/pz", METv.Pz(), weight, 500, -250, 250) 
        
        # Fill signal dependent distributions 
        for signal in self.signalStrings:
            if "1E2Mu" in channel:
                ACand = pairs
                super().FillHist(f"{channel}/{syst}/{signal}/ACand/pt", ACand.Pt(), weight, 300, 0., 300.)
                super().FillHist(f"{channel}/{syst}/{signal}/ACand/eta", ACand.Eta(), weight, 100, -5., 5.)
                super().FillHist(f"{channel}/{syst}/{signal}/ACand/phi", ACand.Phi(), weight, 64, -3.2, 3.2)
                super().FillHist(f"{channel}/{syst}/{signal}/ACand/mass", ACand.M(), weight, 200, 0., 200.)
            else:
                pair1, pair2 = pairs
                mA = int(signal.split("_")[1].split("-")[1]) 
                if abs(pair1.M() - mA) < abs(pair2.M() - mA): ACand, nACand = pair1, pair2
                else:                                         ACand, nACand = pair2, pair2
                super().FillHist(f"{channel}/{syst}/{signal}/ACand/pt", ACand.Pt(), weight, 300, 0., 300.)
                super().FillHist(f"{channel}/{syst}/{signal}/ACand/eta", ACand.Eta(), weight, 100, -5., 5.)
                super().FillHist(f"{channel}/{syst}/{signal}/ACand/phi", ACand.Phi(), weight, 64, -3.2, 3.2)
                super().FillHist(f"{channel}/{syst}/{signal}/ACand/mass", ACand.M(), weight, 200, 0., 200.)
                super().FillHist(f"{channel}/{syst}/{signal}/nACand/pt", nACand.Pt(), weight, 300, 0., 300.)
                super().FillHist(f"{channel}/{syst}/{signal}/nACand/eta", nACand.Eta(), weight, 100, -5., 5.)
                super().FillHist(f"{channel}/{syst}/{signal}/nACand/phi", nACand.Phi(), weight, 64, -3.2, 3.2)
                super().FillHist(f"{channel}/{syst}/{signal}/nACand/mass", nACand.M(), weight, 300, 0., 300.)
                
            score_TTFake = scores[f"{signal}_vs_TTLL_powheg"]
            score_TTX    = scores[f"{signal}_vs_ttX"]
            super().FillHist(f"{channel}/{syst}/{signal}/score_TTFake", score_TTFake, weight, 100, 0., 1.)
            super().FillHist(f"{channel}/{syst}/{signal}/score_TTX", score_TTX, weight, 100, 0., 1.)
            super().FillHist(f"{channel}/{syst}/{signal}/3D", ACand.M(), 
                             score_TTFake, score_TTX, weight,
                             100, mA-5., mA+5.,
                             100, 0., 1.,
                             100, 0., 1.)
    
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
