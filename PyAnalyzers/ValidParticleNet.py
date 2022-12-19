from ROOT import gSystem
from ROOT import TriLeptonBase
gSystem.Load("/cvmfs/cms.cern.ch/slc7_amd64_gcc900/external/lhapdf/6.2.3/lib/libLHAPDF.so")

import sys; sys.path.insert(0, "/data6/Users/choij/SKFlatAnalyzer/python")
import numpy as np
import pandas as pd
import torch
from MLTools.models import ParticleNet, ParticleNetLite
from MLTools.helpers import evtToGraph, predictProba
from MLTools.formats import NodeParticle


class ValidParticleNet(TriLeptonBase):
    def __init__(self):
        super().__init__()
        self.__loadModels()

    def __loadModels(self):
        self.models = {}
        
        # read csv file
        csv = pd.read_csv("/data6/Users/choij/SKFlatAnalyzer/external/ParticleNet/modelInfo.csv",
                          sep=",\s",
                          engine="python")
        for idx in csv.index:
            sig, bkg = csv.loc[idx, 'signal'], csv.loc[idx, 'background']
            model = csv.loc[idx, 'model']
            if model == "ParticleNet":
                thisModel = ParticleNet(num_features=9, num_classes=2)
            elif model == "ParticleNetLite":
                thisModel = ParticleNetLite(num_features=9, num_classes=2)
            else:
                print(f"[ValidParticleNet] Wrong model {model}")
                exit(1)
            modelPath = f"/data6/Users/choij/SKFlatAnalyzer/external/ParticleNet/models/{sig}_vs_{bkg}.pt"
            thisModel.load_state_dict(torch.load(modelPath, map_location=torch.device("cpu")))
            self.models[f"{sig}_vs_{bkg}"] = thisModel


    def executeEvent(self):
        if not super().PassMETFilter(): return None
        ev = super().GetEvent()
        truthColl = super().GetGens()
        rawMuons = super().GetAllMuons()
        rawElectrons = super().GetAllElectrons()
        rawJets = super().GetAllJets()
        METv = ev.GetMETVector()

        #### Object definition
        vetoMuons = list(super().SelectMuons(rawMuons, super().MuonIDs[2], 10., 2.4))
        looseMuons = list(super().SelectMuons(vetoMuons, super().MuonIDs[1], 10., 2.4))
        tightMuons = list(super().SelectMuons(looseMuons, super().MuonIDs[0], 10., 2.4))
        vetoElectrons = list(super().SelectElectrons(rawElectrons, super().ElectronIDs[2], 10., 2.5))
        looseElectrons = list(super().SelectElectrons(vetoElectrons, super().ElectronIDs[1], 10., 2.5))
        tightElectrons = list(super().SelectElectrons(looseElectrons, super().ElectronIDs[0], 10., 2.5))
        jets = super().SelectJets(rawJets, "tight", 20., 2.4)
        jets = list(super().JetsVetoLeptonInside(jets, vetoElectrons, vetoMuons, 0.4))
        bjets = []
        for jet in jets:
            score = jet.GetTaggerResult(3)                      # DeepJet
            wp = super().mcCorr.GetJetTaggingCutValue(3, 1)     # DeepJet Medium
            if score > wp:
                bjets.append(jet)

        # sort objects
        # some contents are modified after sorting?
        vetoMuons.sort(key=lambda x: x.Pt(), reverse=True)
        looseMuons.sort(key=lambda x: x.Pt(), reverse=True)
        tightMuons.sort(key=lambda x: x.Pt(), reverse=True)
        vetoElectrons.sort(key=lambda x: x.Pt(), reverse=True)
        looseElectrons.sort(key=lambda x: x.Pt(), reverse=True)
        tightElectrons.sort(key=lambda x: x.Pt(), reverse=True)
        jets.sort(key=lambda x: x.Pt(), reverse=True)
        bjets.sort(key=lambda x: x.Pt(), reverse=True)
        
        #### event selection
        is3Mu = (len(looseMuons) == 3 and len(vetoMuons) == 3 and \
                len(looseElectrons) == 0 and len(vetoElectrons) == 0)
        is1E2Mu = len(looseMuons) == 2 and len(vetoMuons) == 2 and \
                  len(looseElectrons) == 1 and len(vetoElectrons) == 1
        if not (is3Mu or is1E2Mu): return None
        ## for fake estimation in data
        if not super().IsDATA:
            if not len(tightMuons) == len(looseMuons): return None
            if not len(tightElectrons) == len(looseElectrons): return None

        channel = ""
        ## 1E2Mu baseline
        ## 1. pass EMuTriggers
        ## 2. Exact 2 tight muons and 1 tight electron, no additional lepton
        ## 3. Exists OS muon pair with mass > 12 GeV
        if is1E2Mu:
            if not ev.PassTrigger(super().EMuTriggers): return None

        ## 3Mu baseline
        ## 1. pass DblMuTriggers
        ## 2. Exact 3 tight muons, no additional leptons
        ## 3. Exist OS muon pair,
        ## 4. All OS muon pair mass > 12 GeV
        else:
            if not ev.PassTrigger(super().DblMuTriggers): return None
            mu1, mu2, mu3  = tuple(looseMuons)
            if not mu1.Pt() > 20.: return None
            if not mu2.Pt() > 10.: return None
            if not mu3.Pt() > 10.: return None
            if not abs(mu1.Charge()+mu2.Charge()+mu3.Charge()) == 1: return None
            if mu1.Charge() == mu2.Charge():
                pair1 = mu1 + mu3
                pair2 = mu2 + mu3
            elif mu1.Charge() == mu3.Charge():
                pair1 = mu1 + mu2
                pair2 = mu2 + mu3
            else:   # mu2.Charge() == mu3.Charge()
                pair1 = mu1 + mu2
                pair2 = mu1 + mu3
            if not pair1.M() > 12.: return None
            if not pair2.M() > 12.: return None

            # orthogonality of SR and CR done by bjet multiplicity
            if len(bjets) >= 1:
                if len(jets) >= 2: channel = "SR3Mu"
                else: return None
            else:
                mZ = 91.2
                isOnZ = abs(pair1.M() - mZ) < 10. or abs(pair2.M() - mZ) < 10.
                if isOnZ: channel = "ZFake3Mu"
                else:
                    if abs((mu1+mu2+mu3).M() - mZ) < 10.: channel = "ZGamma3Mu"
                    else: return None

        if "3Mu" not in channel: return None
        ## event selection done
        ## make a graph
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
                             particle.Charge(), particle.IsMuon(), particle.IsElectron(),
                             particle.IsJet(), particle.BtagScore()])
        data = evtToGraph(nodeList, y=None, k=4)
        for key, model in self.models.items():
            score = predictProba(model, data.x, data.edge_index)
            super().FillHist(f"{channel}/{key}/score", score, 1., 100, 0., 1.) 
