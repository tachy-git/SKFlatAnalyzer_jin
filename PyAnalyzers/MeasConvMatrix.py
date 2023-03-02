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

class MeasConvMatrix(TriLeptonBase):
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

        self.systematics = ["Central", "NonpromptUp", "NonpromptDown"]     

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
        objects = {"muons": looseMuons,
                   "electrons": looseElectrons,
                   "jets": jets,
                   "bjets": bjets,
                   "METv": METv,
                   "pairs": pairs,
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
        
        # for conversion measurement
        leptons = std.vector[Lepton]()
        for mu in looseMuons: leptons.emplace_back(mu)
        for ele in looseElectrons: leptons.emplace_back(ele)
        if leptons.at(0).Pt() > 20. and leptons.at(1).Pt() > 20. and leptons.at(2).Pt() > 20.:
            self.measure = "HighPT"
        else:
            self.measure = "LowPT"

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
            
    def FillObjects(self, channel, objects, weight, syst):
        muons = objects["muons"]
        electrons = objects["electrons"]
        jets = objects["jets"]
        bjets = objects["bjets"]
        METv = objects["METv"]
        pairs = objects["pairs"]
        
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

        # Fill ZCands
        if "1E2Mu" in channel:
            if not "ZGamma" in channel: ZCand = pairs       # pairs == pair
            else:                       ZCand = muons.at(0) + muons.at(1) + electrons.at(0) 
            super().FillHist(f"{channel}/{self.measure}/{syst}/ZCand/pt", ZCand.Pt(), weight, 300, 0., 300.)
            super().FillHist(f"{channel}/{self.measure}/{syst}/ZCand/eta", ZCand.Eta(), weight, 100, -5., 5.)
            super().FillHist(f"{channel}/{self.measure}/{syst}/ZCand/phi", ZCand.Phi(), weight, 64, -3.2, 3.2)
            super().FillHist(f"{channel}/{self.measure}/{syst}/ZCand/mass", ZCand.M(), weight, 200, 0., 200.)
        else:
            if not "ZGamma" in channel:
                pair1, pair2 = pairs
                mZ = 91.2
                if abs(pair1.M() - mZ) < abs(pair2.M() - mZ): ZCand, nZCand = pair1, pair2
                else:                                         ZCand, nZCand = pair2, pair1
            else:
                ZCand = muons.at(0) + muons.at(1) + muons.at(2)
                nZCand = NodeParticle()
            super().FillHist(f"{channel}/{self.measure}/{syst}/ZCand/pt", ZCand.Pt(), weight, 300, 0., 300.)
            super().FillHist(f"{channel}/{self.measure}/{syst}/ZCand/eta", ZCand.Eta(), weight, 100, -5., 5.)
            super().FillHist(f"{channel}/{self.measure}/{syst}/ZCand/phi", ZCand.Phi(), weight, 64, -3.2, 3.2)
            super().FillHist(f"{channel}/{self.measure}/{syst}/ZCand/mass", ZCand.M(), weight, 200, 0., 200.)
            super().FillHist(f"{channel}/{self.measure}/{syst}/nZCand/pt", nZCand.Pt(), weight, 300, 0., 300.)
            super().FillHist(f"{channel}/{self.measure}/{syst}/nZCand/eta", nZCand.Eta(), weight, 100, -5., 5.)
            super().FillHist(f"{channel}/{self.measure}/{syst}/nZCand/phi", nZCand.Phi(), weight, 64, -3.2, 3.2)
            super().FillHist(f"{channel}/{self.measure}/{syst}/nZCand/mass", nZCand.M(), weight, 200, 0., 200.)

