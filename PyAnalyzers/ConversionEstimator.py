from ROOT import gSystem
from ROOT import TriLeptonBase
from ROOT import std
from ROOT.JetTagging import Parameters as jParameters
from ROOT import Lepton, Muon, Electron, Jet
gSystem.Load("/cvmfs/cms.cern.ch/slc7_amd64_gcc900/external/lhapdf/6.2.3/lib/libLHAPDF.so")

import os
import pandas as pd
import torch 
from MLTools.models import ParticleNet, ParticleNetLite
from MLTools.helpers import evtToGraph, predictProba
from MLTools.formats import NodeParticle


class ConversionEstimator(TriLeptonBase):
    def __init__(self):
        super().__init__()
        self.__loadModels()
        self.signalStrings = [
                "MHc-70_MA-15", "MHc-70_MA-40", "MHc-70_MA-65",
                "MHc-100_MA-15", "MHc-100_MA-60", "MHc-100_MA-95",
                "MHc-130_MA-15", "MHc-130_MA-55", "MHc-130_MA-90", "MHc-130_MA-125",
                "MHc-160_MA-15", "MHc-160_MA-85", "MHc-160_MA-120", "MHc-160_MA-155"]

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
        if eff == 0.:
            print(mu1.Pt(), mu1.Eta())
            print(mu2.Pt(), mu2.Eta())
            print(mu3.Pt(), mu3.Eta())
        return eff

    def getDblMuTriggerSF(self, muons, sys):
        effData = self.__getDblMuTriggerEff(muons, True, sys)
        effMC = self.__getDblMuTriggerEff(muons, False, sys)

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
        if not (is3Mu or is1E2Mu): return None
        
        # prompt matching
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
        if is1E2Mu:
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
        
        ## for patching sample
        leptons = std.vector[Lepton]()
        for mu in tightMuons: leptons.emplace_back(mu)
        for ele in tightElectrons: leptons.emplace_back(ele)
        if leptons[0].Pt() < 20. or leptons[1].Pt() < 20. or leptons[2].Pt() < 20.:
            measure = "DYJets"
        else:
            measure = "ZGToLLG"
        if not measure in super().MCSample: return None
        ## event selection done
        
        ## set weight
        weight = 1.
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
        for idx, ele in enumerate(tightElectrons, start=1):
            super().FillHist(f"{channel}/electrons/{idx}/pt", ele.Pt(), weight, 300, 0., 300.)
            super().FillHist(f"{channel}/electrons/{idx}/eta", ele.Eta(), weight, 50, -2.5, 2.5)
            super().FillHist(f"{channel}/electrons/{idx}/Phi", ele.Phi(), weight, 64, -3.2, 3.2)
            super().FillHist(f"{channel}/electrons/{idx}/mass", ele.M(), weight, 100, 0., 1.)
        for idx, jet in enumerate(jets, start=1):
            super().FillHist(f"{channel}/jets/{idx}/pt", jet.Pt(), weight, 300, 0., 300.)
            super().FillHist(f"{channel}/jets/{idx}/eta", jet.Eta(), weight, 48, -2.4, 2.4)
            super().FillHist(f"{channel}/jets/{idx}/phi", jet.Phi(), weight, 64, -3.2, 3.2)
            super().FillHist(f"{channel}/jets/{idx}/mass", jet.M(), weight, 100, 0., 100.)
            super().FillHist(f"{channel}/jets/{idx}/charge", jet.Charge(), weight, 200, -5., 5.)
            super().FillHist(f"{channel}/jets/{idx}/btagScore", jet.GetTaggerResult(3), weight, 100, 0., 1.)
        for idx, bjet in enumerate(bjets, start=1):
            super().FillHist(f"{channel}/bjets/{idx}/pt", bjet.Pt(), weight, 300, 0., 300.)
            super().FillHist(f"{channel}/bjets/{idx}/eta", bjet.Eta(), weight, 48, -2.4, 2.4)
            super().FillHist(f"{channel}/bjets/{idx}/phi", bjet.Phi(), weight, 64, -3.2, 3.2)
            super().FillHist(f"{channel}/bjets/{idx}/mass", bjet.M(), weight, 100, 0., 100.) 
            super().FillHist(f"{channel}/bjets/{idx}/charge", bjet.Charge(), weight, 200, -5., 5.)
            super().FillHist(f"{channel}/bjets/{idx}/btagScore", bjet.GetTaggerResult(3), weight, 100, 0., 1.)
        super().FillHist(f"{channel}/METv/pt", METv.Pt(), weight, 300, 0., 300.)
        super().FillHist(f"{channel}/METv/phi", METv.Phi(), weight, 64, -3.2, 3.2)
        
        ## ZCand
        mZ = 91.2
        if "1E2Mu" in channel:
            ZCand = pair
            super().FillHist(f"{channel}/ZCand/pt", ZCand.Pt(), weight, 300, 0., 300.)
            super().FillHist(f"{channel}/ZCand/eta", ZCand.Eta(), weight, 100, -5., 5.)
            super().FillHist(f"{channel}/ZCand/phi", ZCand.Phi(), weight, 64, -3.2, 3.2)
            super().FillHist(f"{channel}/ZCand/mass", ZCand.M(), weight, 200, 0., 200.)
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

        for signal in self.signalStrings:
            mA = int(signal.split("_")[1].split("-")[1])
            if "1E2Mu" in channel:
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
            
            score_TTFake = predictProba(self.models[f"{signal}_vs_TTLL_powheg"], data.x, data.edge_index)
            score_TTX = predictProba(self.models[f"{signal}_vs_ttX"], data.x, data.edge_index)
            super().FillHist(f"{channel}/{signal}/score_TTFake", score_TTFake, weight, 100, 0., 1.)
            super().FillHist(f"{channel}/{signal}/score_TTX", score_TTX, weight, 100, 0., 1.)
            super().FillHist(f"{channel}/{signal}/3D",
                             ACand.M(), score_TTFake, score_TTX, weight,
                             100, mA-5., mA+5.,
                             100, 0., 1.,
                             100, 0., 1.)
        

if __name__ == "__main__":
    m = ConversionEstimator()
    m.SetTreeName("recoTree/SKFlat")
    m.IsDATA = False
    m.MCSample = "TTToHcToWAToMuMu_MHc-130_MA-90"
    m.xsec = 0.015
    m.sumSign = 599702.0
    m.sumW = 3270.46
    m.IsFastSim = False
    m.SetEra("2017")
    if not m.AddFile("/home/choij/workspace/DATA/SKFlat/Run2UltraLegacy_v3/2017/TTToHcToWAToMuMu_MHc-130_MA-90_MultiLepFilter_TuneCP5_13TeV-madgraph-pythia8/SKFlat_Run2UltraLegacy_v3/220714_084244/0000/SKFlatNtuple_2017_MC_14.root"): exit(1)
    if not m.AddFile("/home/choij/workspace/DATA/SKFlat/Run2UltraLegacy_v3/2017/TTToHcToWAToMuMu_MHc-130_MA-90_MultiLepFilter_TuneCP5_13TeV-madgraph-pythia8/SKFlat_Run2UltraLegacy_v3/220714_084244/0000/SKFlatNtuple_2017_MC_5.root"): exit(1)
    m.SetOutfilePath("hists.root")
    m.Init()
    m.initializeAnalyzer()
    m.initializeAnalyzerTools()
    m.SwitchToTempDir()
    m.Loop()
    m.WriteHist()