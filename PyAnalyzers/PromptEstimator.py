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
from MLTools.helpers import evtToGraph, predictProba
from MLTools.formats import NodeParticle


class PromptEstimator(TriLeptonBase):
    def __init__(self):
        super().__init__()
        # at this point, TriLeptonBase.InitializeAnalyzer has not been activate

    def initializePyAnalyzer(self):
        super().initializeAnalyzer()
        if super().Skim1E2Mu: 
            self.channel = "Skim1E2Mu"
            self.features = ["mu1_energy", "mu1_px", "mu1_py", "mu1_pz", "mu1_charge",
                             "mu2_energy", "mu2_px", "mu2_py", "mu2_pz", "mu2_charge",
                             "ele_energy", "ele_px", "ele_py", "ele_pz", "ele_charge",
                             "j1_energy", "j1_px", "j1_py", "j1_pz", "j1_charge", "j1_btagScore",
                             "j2_energy", "j2_px", "j2_py", "j2_pz", "j2_charge", "j2_btagScore",
                             'dR_mu1mu2', 'dR_mu1ele', 'dR_mu2ele', 'dR_j1ele', 'dR_j2ele', "dR_j1j2",
                             "HT", "LT", "MT", "MET", "Nj", "Nb",
                             "avg_dRjets", "avg_btagScore"]
        elif super().Skim3Mu: 
            self.channel = "Skim3Mu"
            self.features = ["mu1_energy", "mu1_px", "mu1_py", "mu1_pz", "mu1_charge",
                             "mu2_energy", "mu2_px", "mu2_py", "mu2_pz", "mu2_charge",
                             "mu3_energy", "mu3_px", "mu3_py", "mu3_pz", "mu3_charge",
                             "j1_energy", "j1_px", "j1_py", "j1_pz", "j1_charge", "j1_btagScore",
                             "j2_energy", "j2_px", "j2_py", "j2_pz", "j2_charge", "j2_btagScore",
                             'dR_mu1mu2', 'dR_mu1mu3', 'dR_mu2mu3',
                             'dR_j1mu1', 'dR_j1mu2', "dR_j1mu3",
                             "dR_j2mu1", "dR_j2mu2", "dR_j2mu3",
                             "dR_j1j2",
                             "HT", "LT", "MT1", "MT2", "MT3", "MET", "Nj", "Nb",
                             "avg_dRjets", "avg_btagScore"]
        else:
            print("Wrong channel")
            exit(1)

        if super().DenseNet:   self.network = "DenseNeuralNet"
        elif super().GraphNet: self.network = "GraphNeuralNet"
        else:
            print("Wrong network")
            exit(1)
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
            modelArch = csv[0][4]
            modelPath = f"{os.environ['DATA_DIR']}/FullRun2/{self.network}/{self.channel}/models/{sig}_vs_{bkg}.pt"
            if self.network == "DenseNeuralNet":
                if self.channel == "Skim1E2Mu": 
                    if modelArch == "SNN": model = SNN(41, 2)
                    else:                  model = SNNLite(41, 2)
                if self.channel == "Skim3Mu" :
                    if modelArch == "SNN": model = SNN(47, 2)
                    else:                  model = SNNLite(47, 2)
            else:               # GraphNeuralNet
                if modelArch == "ParticleNet": model = ParticleNet(9, 2)
                else:                          model = ParticleNetLite(9, 2)
            model.load_state_dict(torch.load(modelPath, map_location=torch.device("cpu")))
            model.eval()
            self.models[f"{sig}_vs_{bkg}"] = model

    def __getDblMuTriggerEff(self, muons, isData, sys):
        assert len(muons) == 3
        mu1, mu2, mu3 = tuple(muons)

        # data
        case1 = super().getTriggerEff(mu1, "Mu17Leg1", isData, sys)
        case1 *= super().getTriggerEff(mu2, "Mu8Leg2", isData, sys)
        case2 = 1.-super().getTriggerEff(mu1, "Mu17Leg1", isData, sys)
        case2 *= super().getTriggerEff(mu2, "Mu17Leg1", isData, sys)
        case2 *= super().getTriggerEff(mu3, "Mu8Leg2", isData, sys)
        case3 = super().getTriggerEff(mu1, "Mu17Leg1", isData, sys)
        case3 *= 1.-super().getTriggerEff(mu2, "Mu8Leg2", isData, sys)
        case3 *= super().getTriggerEff(mu3, "Mu8Leg2", isData, sys)

        eff = case1+case2+case3
        return eff

    def getDblMuTriggerSF(self, muons, sys):
        effData = self.__getDblMuTriggerEff(muons, True, sys)
        effMC = self.__getDblMuTriggerEff(muons, False, sys)
        if effMC == 0 or effData == 0:
            return 1.

        return effData / effMC

    def executeEvent(self):
        if not super().PassMETFilter(): return None
        ev = super().GetEvent()
        rawMuons = super().GetAllMuons()
        rawElectrons = super().GetAllElectrons()
        rawJets = super().GetAllJets()
        METv = ev.GetMETVector()

        #### Object definition
        vetoMuons = super().SelectMuons(rawMuons, super().MuonIDs[2], 10., 2.4)
        tightMuons = super().SelectMuons(vetoMuons, super().MuonIDs[0], 10., 2.4)
        vetoElectrons = super().SelectElectrons(rawElectrons, super().ElectronIDs[2], 10., 2.5)
        tightElectrons = super().SelectElectrons(vetoElectrons, super().ElectronIDs[0], 10., 2.5)
        jets = super().SelectJets(rawJets, "tight", 20., 2.4)
        jets = super().JetsVetoLeptonInside(jets, vetoElectrons, vetoMuons, 0.4)
        bjets = std.vector[Jet]()
        for jet in jets:
            score = jet.GetTaggerResult(3)                      # DeepJet
            wp = super().mcCorr.GetJetTaggingCutValue(3, 1)     # DeepJet Medium
            if score > wp: bjets.emplace_back(jet)

        # sort objects
        sorted(vetoMuons, key=lambda x: x.Pt(), reverse=True)
        sorted(tightMuons, key=lambda x: x.Pt(), reverse=True)
        sorted(vetoElectrons, key=lambda x: x.Pt(), reverse=True)
        sorted(tightElectrons, key=lambda x: x.Pt(), reverse=True)
        sorted(jets, key=lambda x: x.Pt(), reverse=True)
        sorted(bjets, key=lambda x: x.Pt(), reverse=True)

        #### event selection
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
            truth = super().GetGens()
            promptMuons = std.vector[Muon]()
            promptElectrons = std.vector[Electron]()
            for mu in tightMuons:
                if super().GetLeptonType(mu, truth) > 0: promptMuons.emplace_back(mu)
            for ele in tightElectrons:
                if super().GetLeptonType(ele, truth) > 0: promptElectrons.emplace_back(ele)

            if len(promptMuons) != len(tightMuons): return None
            if len(promptElectrons) != len(tightElectrons): return None

        channel = ""
        ## 1E2Mu baseline
        ## 1. pass EMuTriggers
        ## 2. Exact 2 tight muons and 1 tight electron, no additional lepton
        ## 3. Exists OS muon pair with mass > 12 GeV
        ## 4. At least two jets
        if self.channel == "Skim1E2Mu":
            if not ev.PassTrigger(super().EMuTriggers): return None
            leptons = std.vector[Lepton]()
            for mu in tightMuons: leptons.emplace_back(mu)
            for ele in tightElectrons: leptons.emplace_back(ele)
            mu1, mu2, ele = tuple(leptons)
            passLeadMu = mu1.Pt() > 25. and ele.Pt() > 15.
            passLeadEle = mu1.Pt() > 10. and ele.Pt() > 25.
            passSafeCut = passLeadMu or passLeadEle
            if not passSafeCut: return None
            if not mu1.Charge()+mu2.Charge() == 0: return None
            pair = mu1 + mu2
            if not pair.M() > 12.: return None
            if not jets.size() >= 2: return None

            # orthogonality of SR and CR done by bjet multiplicity
            if bjets.size() >= 1:
                channel = "SR1E2Mu"
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
        elif self.channel == "Skim3Mu":
            if not ev.PassTrigger(super().DblMuTriggers): return None
            mu1, mu2, mu3  = tuple(tightMuons)
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

        # if not ("1E2Mu" in channel or "3Mu" in channel): return None

        ## set weight
        weight = 1.
        if not super().IsDATA:
            weight *= super().MCweight()
            weight *= ev.GetTriggerLumi("Full")
            w_prefire = super().GetPrefireWeight(0)
            w_pileup = super().GetPileUpWeight(super().nPileUp, 0)
            w_muonIDSF = 1.
            w_dblMuTrigSF = 1.
            if "3Mu" in channel:
                w_muonIDSF = 1.
                for mu in tightMuons:
                    w_muonIDSF *= super().getMuonIDSF(mu, 0)

                w_dblMuTrigSF = self.getDblMuTriggerSF(tightMuons, 0)
            weight *= w_prefire            # print(f"w_prefire: {w_prefire}")
            weight *= w_pileup             # print(f"w_pileup: {w_pileup}")
            weight *= w_muonIDSF           # print(f"muonID: {w_muonIDSF}")
            weight *= w_dblMuTrigSF        # print(f"muontrig: {w_dblMuTrigSF}")

            # b-tagging
            jtp_DeepJet_Medium = jParameters(3, 1, 0, 1)    # DeepJet, Medium, incl, mujets
            vjets = std.vector[Jet]()
            for j in jets: vjets.emplace_back(j)
            w_btag = super().mcCorr.GetBTaggingReweight_1a(vjets, jtp_DeepJet_Medium)
            weight *= w_btag

        ## fill input observables
        for idx, mu in enumerate(tightMuons, start=1):
            super().FillHist(f"{channel}/muons/{idx}/pt", mu.Pt(), weight, 300, 0., 300.)
            super().FillHist(f"{channel}/muons/{idx}/eta", mu.Eta(), weight, 48, -2.4, 2.4)
            super().FillHist(f"{channel}/muons/{idx}/phi", mu.Phi(), weight, 64, -3.2, 3.2)
            super().FillHist(f"{channel}/muons/{idx}/mass", mu.M(), weight, 10, 0., 1.)
            super().FillHist(f"{channel}/muons/{idx}/energy", mu.E(), weight, 300, 0., 300.)
            super().FillHist(f"{channel}/muons/{idx}/px", mu.Px(), weight, 500, -250., 250.)
            super().FillHist(f"{channel}/muons/{idx}/py", mu.Py(), weight, 500, -250., 250.)
            super().FillHist(f"{channel}/muons/{idx}/pz", mu.Pz(), weight, 500, -250., 250.)
            super().FillHist(f"{channel}/muons/{idx}/charge", mu.Charge(), weight, 3, -1., 1.)
        for idx, ele in enumerate(tightElectrons, start=1):
            super().FillHist(f"{channel}/electrons/{idx}/pt", ele.Pt(), weight, 300, 0., 300.)
            super().FillHist(f"{channel}/electrons/{idx}/eta", ele.Eta(), weight, 50, -2.5, 2.5)
            super().FillHist(f"{channel}/electrons/{idx}/Phi", ele.Phi(), weight, 64, -3.2, 3.2)
            super().FillHist(f"{channel}/electrons/{idx}/mass", ele.M(), weight, 100, 0., 1.)
            super().FillHist(f"{channel}/electrons/{idx}/energy", ele.E(), weight, 300, 0., 300.)
            super().FillHist(f"{channel}/electrons/{idx}/px", ele.Px(), weight, 500, -250., 250.)
            super().FillHist(f"{channel}/electrons/{idx}/py", ele.Py(), weight, 500, -250., 250.)
            super().FillHist(f"{channel}/electrons/{idx}/pz", ele.Pz(), weight, 500, -250., 250.)
            super().FillHist(f"{channel}/electrons/{idx}/charge", ele.Charge(), weight, 3, -1., 1.)
        for idx, jet in enumerate(jets, start=1):
            super().FillHist(f"{channel}/jets/{idx}/pt", jet.Pt(), weight, 300, 0., 300.)
            super().FillHist(f"{channel}/jets/{idx}/eta", jet.Eta(), weight, 48, -2.4, 2.4)
            super().FillHist(f"{channel}/jets/{idx}/phi", jet.Phi(), weight, 64, -3.2, 3.2)
            super().FillHist(f"{channel}/jets/{idx}/mass", jet.M(), weight, 100, 0., 100.)
            super().FillHist(f"{channel}/jets/{idx}/btagScore", jet.GetTaggerResult(3), weight, 100, 0., 1.)
            super().FillHist(f"{channel}/jets/{idx}/energy", jet.E(), weight, 300, 0., 300.)
            super().FillHist(f"{channel}/jets/{idx}/px", jet.Px(), weight, 500, -250., 250.)
            super().FillHist(f"{channel}/jets/{idx}/py", jet.Py(), weight, 500, -250., 250.)
            super().FillHist(f"{channel}/jets/{idx}/pz", jet.Pz(), weight, 500, -250., 250.)
            super().FillHist(f"{channel}/jets/{idx}/charge", jet.Charge(), weight, 200, -1., 1.)
        for idx, bjet in enumerate(bjets, start=1):
            super().FillHist(f"{channel}/bjets/{idx}/pt", bjet.Pt(), weight, 300, 0., 300.)
            super().FillHist(f"{channel}/bjets/{idx}/eta", bjet.Eta(), weight, 48, -2.4, 2.4)
            super().FillHist(f"{channel}/bjets/{idx}/phi", bjet.Phi(), weight, 64, -3.2, 3.2)
            super().FillHist(f"{channel}/bjets/{idx}/mass", bjet.M(), weight, 100, 0., 100.)
            super().FillHist(f"{channel}/bjets/{idx}/btagScore", bjet.GetTaggerResult(3), weight, 100, 0., 1.)
            super().FillHist(f"{channel}/bjets/{idx}/energy", bjet.E(), weight, 300, 0., 300.)
            super().FillHist(f"{channel}/bjets/{idx}/px", bjet.Px(), weight, 500, -250., 250.)
            super().FillHist(f"{channel}/bjets/{idx}/py", bjet.Py(), weight, 500, -250., 250.)
            super().FillHist(f"{channel}/bjets/{idx}/pz", bjet.Pz(), weight, 500, -250., 250.)
            super().FillHist(f"{channel}/bjets/{idx}/charge", mu.Charge(), weight, 200, -1., 1.)
        super().FillHist(f"{channel}/jets/size", jets.size(), weight, 20, 0., 20.)
        super().FillHist(f"{channel}/bjets/size", bjets.size(), weight, 15, 0., 15.)
        super().FillHist(f"{channel}/METv/pt", METv.Pt(), weight, 300, 0., 300.)
        super().FillHist(f"{channel}/METv/phi", METv.Phi(), weight, 64, -3.2, 3.2)
        super().FillHist(f"{channel}/METv/energy", METv.E(), weight, 300, 0., 300.)
        super().FillHist(f"{channel}/METv/px", METv.Px(), weight, 500, -250., 250.)
        super().FillHist(f"{channel}/METv/py", METv.Py(), weight, 500, -250., 250.)
        super().FillHist(f"{channel}/METv/pz", METv.Pz(), weight, 500, -250., 250.)


        ## ZCand and additional inputs
        if self.network == "DenseNeuralNet":
            HT = sum([j.Pt() for j in jets])
            dRjets = []
            for idx1, idx2 in combinations(range(jets.size()), 2):
                dRjets.append(jets[idx1].DeltaR(jets[idx2]))
            evt_dRjets = sum(dRjets) / len(dRjets)
            evt_btagScore = sum([j.GetTaggerResult(3) for j in jets]) / jets.size()
            
            super().FillHist(f"{channel}/HT", HT, weight, 600, 0., 600.)
            super().FillHist(f"{channel}/evt_dRjets", evt_dRjets, weight, 50, 0., 5.)
            super().FillHist(f"{channel}/evt_btagScore", evt_btagScore, weight, 100, 0., 1.)
        
        mZ = 91.2
        if "1E2Mu" in channel:
            ZCand = pair
            super().FillHist(f"{channel}/ZCand/pt", ZCand.Pt(), weight, 300, 0., 300.)
            super().FillHist(f"{channel}/ZCand/eta", ZCand.Eta(), weight, 100, -5., 5.)
            super().FillHist(f"{channel}/ZCand/phi", ZCand.Phi(), weight, 64, -3.2, 3.2)
            super().FillHist(f"{channel}/ZCand/mass", ZCand.M(), weight, 200, 0., 200.)
            
            if self.network == "DenseNeuralNet":
                mu1 = tightMuons.at(0)
                mu2 = tightMuons.at(1)
                ele = tightElectrons.at(0)
                j1 = jets.at(0)
                j2 = jets.at(1)
                MT = (ele+METv).Mt()
                LT = sum([l.Pt() for l in tightElectrons+tightMuons])
                dRmu1mu2 = mu1.DeltaR(mu2)
                dRmu1ele = mu1.DeltaR(ele)
                dRmu2ele = mu2.DeltaR(ele)
                dRj1j2 = j1.DeltaR(j2)
                
                super().FillHist(f"{channel}/MT", MT, weight, 300, 0., 300.)
                super().FillHist(f"{channel}/LT", LT, weight, 600, 0., 600.)
                super().FillHist(f"{channel}/dRmu1mu2", dRmu1mu2, weight, 50, 0., 50.)
                super().FillHist(f"{channel}/dRmu1ele", dRmu1ele, weight, 50, 0., 50.)
                super().FillHist(f"{channel}/dRmu2ele", dRmu2ele, weight, 50, 0., 50.)
                super().FillHist(f"{channel}/dRj1j2", dRj1j2, weight, 50, 0., 5.)  
        else:
            if abs(pair1.M() - mZ) < abs(pair2.M() - mZ):
                ZCand, nZCand = pair1, pair2
            else:
                ZCand, nZCand = pair2, pair1

            super().FillHist(f"{channel}/ZCand/pt", ZCand.Pt(), weight, 300, 0., 300.)
            super().FillHist(f"{channel}/ZCand/eta", ZCand.Eta(), weight, 100, -5., 5.)
            super().FillHist(f"{channel}/ZCand/phi", ZCand.Phi(), weight, 64, -3.2, 3.2)
            super().FillHist(f"{channel}/ZCand/mass", ZCand.M(), weight, 200, 0., 200.)
            super().FillHist(f"{channel}/nZCand/pt", nZCand.Pt(), weight, 300, 0., 300.)
            super().FillHist(f"{channel}/nZCand/eta", nZCand.Eta(), weight, 100, -5., 5.)
            super().FillHist(f"{channel}/nZCand/phi", nZCand.Phi(), weight, 64, -3.2, 3.2)
            super().FillHist(f"{channel}/nZCand/mass", nZCand.M(), weight, 200, 0., 200.)

            if self.network == "DenseNeuralNet":
                mu1 = tightMuons.at(0)
                mu2 = tightMuons.at(1)
                mu3 = tightMuons.at(2)
                j1 = jets.at(0)
                j2 = jets.at(1)
                MT1 = (mu1+METv).Mt()
                MT2 = (mu2+METv).Mt()
                MT3 = (mu3+METv).Mt()
                LT = sum([l.Pt() for l in tightMuons])
                dRmu1mu2 = mu1.DeltaR(mu2)
                dRmu1mu3 = mu1.DeltaR(mu3)
                dRmu2mu3 = mu2.DeltaR(mu3)
                dRj1mu1 = j1.DeltaR(mu1)
                dRj1mu2 = j1.DeltaR(mu2)
                dRj1mu3 = j1.DeltaR(mu3)
                dRj2mu1 = j2.DeltaR(mu1)
                dRj2mu2 = j2.DeltaR(mu2)
                dRj2mu3 = j2.DeltaR(mu3)
                dRj1j2 = j1.DeltaR(j2)
                
                super().FillHist(f"{channel}/MT1", MT1, weight, 300, 0., 300.)
                super().FillHist(f"{channel}/MT2", MT2, weight, 300, 0., 300.)
                super().FillHist(f"{channel}/MT3", MT3, weight, 300, 0., 300.) 
                super().FillHist(f"{channel}/LT", LT, weight, 600, 0., 600.)
                super().FillHist(f"{channel}/dRmu1mu2", dRmu1mu2, weight, 50, 0., 5.)
                super().FillHist(f"{channel}/dRmu1mu3", dRmu1mu3, weight, 50, 0., 5.)
                super().FillHist(f"{channel}/dRmu2mu3", dRmu2mu3, weight, 50, 0., 5.)
                super().FillHist(f"{channel}/dRj1mu1", dRj1mu1, weight, 50, 0., 5.)
                super().FillHist(f"{channel}/dRj1mu2", dRj1mu2, weight, 50, 0., 5.)
                super().FillHist(f"{channel}/dRj1mu3", dRj1mu3, weight, 50, 0., 5.)
                super().FillHist(f"{channel}/dRj2mu1", dRj2mu1, weight, 50, 0., 5.)
                super().FillHist(f"{channel}/dRj2mu2", dRj2mu2, weight, 50, 0., 5.)
                super().FillHist(f"{channel}/dRj2mu3", dRj2mu3, weight, 50, 0., 5.)
                super().FillHist(f"{channel}/dRj1j2", dRj1j2, weight, 50, 0., 5.) 
       
        if self.network == "DenseNeuralNet": data = self.getDenseInput(tightMuons, tightElectrons, jets, bjets, METv) 
        else:                                data = self.getGraphInput(tightMuons, tightElectrons, jets, METv)
        for signal in self.signalStrings:
            mA = int(signal.split("_")[1].split("-")[1])
            if self.channel == "Skim1E2Mu":
                ACand = pair
                super().FillHist(f"{channel}/{signal}/ACand/pt", ACand.Pt(), weight, 300, 0., 300.)
                super().FillHist(f"{channel}/{signal}/ACand/eta", ACand.Eta(), weight, 100, -5., 5.)
                super().FillHist(f"{channel}/{signal}/ACand/phi", ACand.Phi(), weight, 64, -3.2, 3.2)
                super().FillHist(f"{channel}/{signal}/ACand/mass", ACand.M(), weight, 200, 0., 200.)
            else:
                if abs(pair1.M() - mA) < abs(pair2.M() - mA):
                    ACand, nACand = pair1, pair2
                else:
                    ACand, nACand = pair2, pair2

                super().FillHist(f"{channel}/{signal}/ACand/pt", ACand.Pt(), weight, 300, 0., 300.)
                super().FillHist(f"{channel}/{signal}/ACand/eta", ACand.Eta(), weight, 100, -5., 5.)
                super().FillHist(f"{channel}/{signal}/ACand/phi", ACand.Phi(), weight, 64, -3.2, 3.2)
                super().FillHist(f"{channel}/{signal}/ACand/mass", ACand.M(), weight, 200, 0., 200.)
                super().FillHist(f"{channel}/{signal}/nACand/pt", nACand.Pt(), weight, 300, 0., 300.)
                super().FillHist(f"{channel}/{signal}/nACand/eta", nACand.Eta(), weight, 100, -5., 5.)
                super().FillHist(f"{channel}/{signal}/nACand/phi", nACand.Phi(), weight, 64, -3.2, 3.2)
                super().FillHist(f"{channel}/{signal}/nACand/mass", nACand.M(), weight, 300, 0., 300.)
            
            if self.network == "DenseNeuralNet":
                score_TTFake = self.getDenseScore(f"{signal}_vs_TTLL_powheg", data)
                score_TTX = self.getDenseScore(f"{signal}_vs_ttX", data)
            else:
                score_TTFake = self.getGraphScore(f"{signal}_vs_TTLL_powheg", data)
                score_TTX = self.getGraphScore(f"{signal}_vs_ttX", data)
            super().FillHist(f"{channel}/{signal}/score_TTFake", score_TTFake, weight, 100, 0., 1.)
            super().FillHist(f"{channel}/{signal}/score_TTX", score_TTX, weight, 100, 0., 1.)
            super().FillHist(f"{channel}/{signal}/3D",
                             ACand.M(), score_TTFake, score_TTX, weight,
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
                             particle.Charge(), particle.IsMuon(), particle.IsElectron(),
                             particle.IsJet(), particle.BtagScore()])
        data = evtToGraph(nodeList, y=None, k=4)
        return data

    def getGraphScore(self, modelKey, data):
        with torch.no_grad():
            out = self.models[modelKey](data.x, data.edge_index)
        return out.numpy()[0][1]

if __name__ == "__main__":
    m = PromptEstimator()
    m.SetTreeName("recoTree/SKFlat")
    m.IsDATA = False
    m.MCSample = "TTToHcToWAToMuMu_MHc-130_MA-90"
    m.xsec = 0.015
    m.sumSign = 599702.0
    m.sumW = 3270.46
    m.IsFastSim = False
    m.SetEra("2017")
    m.Userflags = std.vector[std.string]()
    m.Userflags.emplace_back("Skim3Mu")
    m.Userflags.emplace_back("DenseNet")
    if not m.AddFile("/home/choij/workspace/DATA/SKFlat/Run2UltraLegacy_v3/2017/TTToHcToWAToMuMu_MHc-130_MA-90_MultiLepFilter_TuneCP5_13TeV-madgraph-pythia8/SKFlat_Run2UltraLegacy_v3/220714_084244/0000/SKFlatNtuple_2017_MC_1.root"): exit(1)
    #if not m.AddFile("/home/choij/workspace/DATA/SKFlat/Run2UltraLegacy_v3/2017/TTToHcToWAToMuMu_MHc-130_MA-90_MultiLepFilter_TuneCP5_13TeV-madgraph-pythia8/SKFlat_Run2UltraLegacy_v3/220714_084244/0000/SKFlatNtuple_2017_MC_8.root"): exit(1)
    m.SetOutfilePath("histsDense.root")
    m.Init()
    m.initializePyAnalyzer()
    m.initializeAnalyzerTools()
    m.SwitchToTempDir()
    m.Loop()
    m.WriteHist()
