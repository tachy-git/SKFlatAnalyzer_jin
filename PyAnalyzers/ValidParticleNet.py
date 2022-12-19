from ROOT import gSystem
from ROOT import TriLeptonBase
gSystem.Load("/cvmfs/cms.cern.ch/slc7_amd64_gcc900/external/lhapdf/6.2.3/lib/libLHAPDF.so")

import os
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
        csv = pd.read_csv(f"{os.environ['SKFlat_WD']}/external/ParticleNet/modelInfo.csv",
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
            modelPath = f"{os.environ['SKFlat_WD']}/external/ParticleNet/models/{sig}_vs_{bkg}.pt"
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
        ## temporarily no data driven fake estimation
        #if not super().IsDATA:
        if not len(tightMuons) == len(looseMuons): return None
        if not len(tightElectrons) == len(looseElectrons): return None

        channel = ""
        ## 1E2Mu baseline
        ## 1. pass EMuTriggers
        ## 2. Exact 2 tight muons and 1 tight electron, no additional lepton
        ## 3. Exists OS muon pair with mass > 12 GeV
        if is1E2Mu:
            if not ev.PassTrigger(super().EMuTriggers): return None
            mu1, mu2, ele = tuple(looseMuons+looseElectrons)
            passLeadMu = mu1.Pt() > 25. and ele.Pt() > 15.
            passLeadEle = mu1.Pt() > 10. and ele.Pt() > 25.
            passSafeCut = passLeadMu or passLeadEle
            if not passSafeCut: return None
            if not mu1.Charge()+mu2.Charge() == 0: return None
            pair = mu1 + mu2
            if not pair.M() > 12.: return None

            # orthogonality of SR and CR done by bjet multiplicity
            if len(bjets) >= 1:
                if len(jets) >= 2: channel = "SR1E2Mu"
                else: return None
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

        if not ("1E2Mu" in channel or "3Mu" in channel): return None
        ## event selection done
        
        ## set weight
        weight = 1.
        if not super().IsDATA:
            weight *= super().MCweight()
            weight *= ev.GetTriggerLumi("Full") 

        ## fill input observables
        for idx, mu in enumerate(looseMuons):
            super().FillHist(f"{channel}/muons/{idx}/pt", mu.Pt(), weight, 300, 0., 300.)
            super().FillHist(f"{channel}/muons/{idx}/eta", mu.Eta(), weight, 48, -2.4, 2.4)
            super().FillHist(f"{channel}/muons/{idx}/phi", mu.Phi(), weight, 64, -3.2, 3.2)
            super().FillHist(f"{channel}/muons/{idx}/mass", mu.M(), weight, 10, 0., 1.)
        for idx, ele in enumerate(looseElectrons):
            super().FillHist(f"{channel}/electrons/{idx}/pt", ele.Pt(), weight, 300, 0., 300.)
            super().FillHist(f"{channel}/electrons/{idx}/eta", ele.Eta(), weight, 50, -2.5, 2.5)
            super().FillHist(f"{channel}/electrons/{idx}/Phi", ele.Phi(), weight, 64, -3.2, 3.2)
            super().FillHist(f"{channel}/electrons/{idx}/mass", ele.M(), weight, 100, 0., 1.)
        for idx, jet in enumerate(jets):
            super().FillHist(f"{channel}/jets/{idx}/pt", jet.Pt(), weight, 300, 0., 300.)
            super().FillHist(f"{channel}/jets/{idx}/eta", jet.Eta(), weight, 48, -2.4, 2.4)
            super().FillHist(f"{channel}/jets/{idx}/phi", jet.Phi(), weight, 64, -3.2, 3.2)
            super().FillHist(f"{channel}/jets/{idx}/mass", jet.M(), weight, 100, 0., 100.)
            super().FillHist(f"{channel}/jets/{idx}/charge", jet.Charge(), weight, 200, -5., 5.)
            super().FillHist(f"{channel}/jets/{idx}/btagScore", jet.GetTaggerResult(3), weight, 100, 0., 1.)
        for idx, bjet in enumerate(bjets):
            super().FillHist(f"{channel}/bjets/{idx}/pt", bjet.Pt(), weight, 300, 0., 300.)
            super().FillHist(f"{channel}/bjets/{idx}/eta", bjet.Eta(), weight, 48, -2.4, 2.4)
            super().FillHist(f"{channel}/bjets/{idx}/phi", bjet.Phi(), weight, 64, -3.2, 3.2)
            super().FillHist(f"{channel}/bjets/{idx}/mass", bjet.M(), weight, 100, 0., 100.) 
            super().FillHist(f"{channel}/bjets/{idx}/charge", bjet.Charge(), weight, 200, -5., 5.)
            super().FillHist(f"{channel}/bjets/{idx}/btagScore", bjet.GetTaggerResult(3), weight, 100, 0., 1.)
        super().FillHist(f"{channel}/METv/pt", METv.Pt(), weight, 300, 0., 300.)
        super().FillHist(f"{channel}/METv/phi", METv.Phi(), weight, 64, -3.2, 3.2)
    
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
            super().FillHist(f"{channel}/{key}/score", score, weight, 100, 0., 1.)

            # pair
            if "1E2Mu" in channel:
                ACand = pair
                super().FillHist(f"{channel}/{key}/ACand/pt", ACand.Pt(), weight, 300, 0., 300.)
                super().FillHist(f"{channel}/{key}/ACand/eta", ACand.Eta(), weight, 100, -5., 5.)
                super().FillHist(f"{channel}/{key}/ACand/phi", ACand.Phi(), weight, 64, -3.2, 3.2)
                super().FillHist(f"{channel}/{key}/ACand/mass", ACand.M(), weight, 200, 0., 200.)
            else: # 3Mu
                mA = int(string.split("_")[1].split("-")[1])
                if abs(pair1.M() - mA) < abs(pair2.M() - mA):
                    ACand, nACand = pair1, pair2
                else:
                    ACand, nACand = pair2, pair2
                    
                super().FillHist(f"{channel}/{key}/ACand/pt", ACand.Pt(), weight, 300, 0., 300.)                
                super().FillHist(f"{channel}/{key}/ACand/eta", ACand.Eta(), weight, 100, -5., 5.)
                super().FillHist(f"{channel}/{key}/ACand/phi", ACand.Phi(), weight, 64, -3.2, 3.2)
                super().FillHist(f"{channel}/{key}/ACand/mass", ACand.M(), weight, 200, 0., 200.) 
                super().FillHist(f"{channel}/{key}/nACand/pt", nACand.Pt(), weight, 300, 0., 300.)
                super().FillHist(f"{channel}/{key}/nACand/eta", nACand.Eta(), weight, 100, -5., 5.)
                super().FillHist(f"{channel}/{key}/nACand/phi", nACand.Phi(), weight, 64, -3.2, 3.2)
                super().FillHist(f"{channel}/{key}/nACand/mass", nACand.M(), weight, 300, 0., 300.)

#if __name__ == "__main__":
#    m = ValidParticleNet()
#    m.SetTreeName("recoTree/SKFlat")
#    m.IsDATA = False
#    m.MCSample = "TTToHcToWAToMuMu_MHc-130_MA-90"
#    m.xsec = 0.015
#    m.sumSign = 599702.0
#    m.sumW = 3270.46
#    m.IsFastSim = False
#    m.SetEra("2017")
#    if not m.AddFile("/home/choij/workspace/DATA/SKFlat/Run2UltraLegacy_v3/2017/TTToHcToWAToMuMu_MHc-130_MA-90_MultiLepFilter_TuneCP5_13TeV-madgraph-pythia8/SKFlat_Run2UltraLegacy_v3/220714_084244/0000/SKFlatNtuple_2017_MC_14.root"): exit(1)
#    if not m.AddFile("/home/choij/workspace/DATA/SKFlat/Run2UltraLegacy_v3/2017/TTToHcToWAToMuMu_MHc-130_MA-90_MultiLepFilter_TuneCP5_13TeV-madgraph-pythia8/SKFlat_Run2UltraLegacy_v3/220714_084244/0000/SKFlatNtuple_2017_MC_5.root"): exit(1)
#    m.SetOutfilePath("hists.root")
#    m.Init()
#    m.initializeAnalyzer()
#    m.initializeAnalyzerTools()
#    m.SwitchToTempDir()
#    m.Loop()
#    m.WriteHist()
