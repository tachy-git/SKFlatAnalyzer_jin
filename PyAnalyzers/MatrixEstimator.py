from ROOT import gSystem
from ROOT import TriLeptonBase
from ROOT import TString
from ROOT.std import vector
from ROOT.JetTagging import Parameters as jParameters
from ROOT import Lepton, Muon, Electron, Jet, Particle
gSystem.Load("/cvmfs/cms.cern.ch/slc7_amd64_gcc900/external/lhapdf/6.2.3/lib/libLHAPDF.so")

from itertools import product
from MLTools.helpers import loadModels
from MLTools.helpers import getGraphInput, getGraphScore

class MatrixEstimator(TriLeptonBase):
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

        self.systematics = ["Central", "NonpromptUp", "NonpromptDown"]
        self.signalStrings = ["MHc-160_MA-85", "MHc-130_MA-90", "MHc-100_MA-95"]
        self.backgroundStrings = ["nonprompt", "diboson", "ttZ"]
    
        self.models = loadModels("GraphNeuralNet", self.channel, self.signalStrings, self.backgroundStrings)
        
    def executeEvent(self):
        if not super().PassMETFilter(): return None
        ev = super().GetEvent()
        rawMuons = super().GetAllMuons()
        rawElectrons = super().GetAllElectrons()
        rawJets = super().GetAllJets()
        METv = ev.GetMETVector()
        truth = super().GetGens() if not super().IsDATA else None

        vetoMuons, looseMuons, tightMuons, vetoElectrons, looseElectrons, tightElectrons, jets, bjets = self.defineObjects(rawMuons, rawElectrons, rawJets)
        thisChannel = self.selectEvent(ev, truth, vetoMuons, looseMuons, tightMuons, vetoElectrons, looseElectrons, tightElectrons, jets, bjets, METv)
        
        if thisChannel is None: return None
        
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
                raise NotImplementedError(f"[MatrixEstimator::executeEvent] Wrong systematics {syst}")
            self.FillObjects(thisChannel, objects, weight, syst)
        
        
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
        bjets = vector[Jet]()
        for jet in jets:
            btagScore = jet.GetTaggerResult(3)                  # DeepJet score
            wp = super().mcCorr.GetJetTaggingCutValue(3, 1)     # DeepJet Medium
            if btagScore > wp: bjets.emplace_back(jet)
            
        vetoMuons = vector[Muon](sorted(vetoMuons, key=lambda x: x.Pt(), reverse=True))
        looseMuons = vector[Muon](sorted(looseMuons, key=lambda x: x.Pt(), reverse=True))
        tightMuons = vector[Muon](sorted(tightMuons, key=lambda x: x.Pt(), reverse=True))
        vetoElectrons = vector[Electron](sorted(vetoElectrons, key=lambda x: x.Pt(), reverse=True))
        looseElectrons = vector[Electron](sorted(looseElectrons, key=lambda x: x.Pt(), reverse=True))
        tightElectrons = vector[Electron](sorted(tightElectrons, key=lambda x: x.Pt(), reverse=True))
        jets = vector[Jet](sorted(jets, key=lambda x: x.Pt(), reverse=True))
        bjets = vector[Jet](sorted(bjets, key=lambda x: x.Pt(), reverse=True))
        
        return (vetoMuons, looseMuons, tightMuons, vetoElectrons, looseElectrons, tightElectrons, jets, bjets)
    
    def selectEvent(self, event, truth, vetoMuons, looseMuons, tightMuons, vetoElectrons, looseElectrons, tightElectrons, jets, bjets, METv):
        is3Mu = (looseMuons.size() == 3 and vetoMuons.size() == 3 and \
                 looseElectrons.size() == 0 and vetoElectrons.size() == 0)
        is1E2Mu = (looseMuons.size() == 2 and vetoMuons.size() == 2 and \
                   looseElectrons.size() == 1 and vetoElectrons.size() == 1)
        
        #### not all leptons tight
        if self.channel == "Skim1E2Mu":
            if not is1E2Mu: return None
            if (tightMuons.size() == looseMuons.size()) and (tightElectrons.size() == looseElectrons.size()): return None

        if self.channel == "Skim3Mu":
            if not is3Mu: return None
            if tightMuons.size() == looseMuons.size(): return None
       
        # for conversion samples
        if super().MCSample in ["DYJets_MG", "DYJets10to50_MG", "TTG", "WWG"]:
            # at least one conversion lepton should exist
            # internal conversion: 4, 5
            # external conversion: -5, -6
            convMuons = vector[Muon]()
            fakeMuons = vector[Muon]()
            convElectrons = vector[Electron]()
            for mu in looseMuons:
                if super().GetLeptonType(mu, truth) in [4, 5, -5, -6]: convMuons.emplace_back(mu)
                if super().GetLeptonType(mu, truth) in [-1, -2, -3, -4]: fakeMuons.emplace_back(mu)
            for ele in looseElectrons:
                if super().GetLeptonType(ele, truth) in [4, 5, -5, -6]: convElectrons.emplace_back(ele)
            if self.channel == "Skim1E2Mu":
                # remove hadronic contribution
                if not fakeMuons.size() == 0: return None
                if not convElectrons.size() == 1: return None
            if self.channel == "Skim3Mu":
                if not fakeMuons.size() == 0: return None
                if not convMuons.size() == 1: return None

        ## 1E2Mu baseline
        ## 1. pass EMuTriggers
        ## 2. Exact 2 tight muons and 1 tight electron, no additional lepton
        ## 3. Exists OS muon pair with mass > 12 GeV
        ## 4. At least two jets
        if self.channel == "Skim1E2Mu":
            if not event.PassTrigger(super().EMuTriggers): return None
            leptons = vector[Lepton]()
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
            if bjets.size() == 0:
                isOnZ = abs(pair.M() - 91.2) < 10.
                if isOnZ: return "ZFake1E2Mu"
                else:     return None
            else:
                return "SR1E2Mu"

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
            if bjets.size() == 0:
                isOnZ = abs(pair1.M() - 91.2) < 10. or abs(pair2.M() - 91.2) < 10.
                if isOnZ: return "ZFake3Mu"
                else:     return None
            else:
                return "SR3Mu"
        
        raise EOFError

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
            raise NotImplementedError(f"[PromptEstimator::makePair] wrong no. of muons {muons.size}")
        
    #### Get scores for each event
    def evalScore(self, muons, electrons, jets, bjets, METv):
        scores = {}
        data = getGraphInput(muons, electrons, jets, bjets, METv)
        for sig, bkg in product(self.signalStrings, self.backgroundStrings):
            scores[f"{sig}_vs_{bkg}"] = getGraphScore(self.models[f"{sig}_vs_{bkg}"], data)
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
        
        ## fill base observables
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
        
        # Fill discrimination variable
        if "1E2Mu" in channel:
            pair = pairs
            super().FillHist(f"{channel}/{syst}/pair/pt", pair.Pt(), weight, 300, 0., 300.)
            super().FillHist(f"{channel}/{syst}/pair/eta", pair.Eta(), weight, 100, -5., 5.)
            super().FillHist(f"{channel}/{syst}/pair/phi", pair.Phi(), weight, 64, -3.2, 3.2)
            super().FillHist(f"{channel}/{syst}/pair/mass", pair.M(), weight, 200, 0., 200.)
        else:
            pair1, pair2 = pairs
            super().FillHist(f"{channel}/{syst}/stack/pt", pair1.Pt(), weight, 300, 0., 300.)
            super().FillHist(f"{channel}/{syst}/stack/eta", pair1.Eta(), weight, 100, -5., 5.)
            super().FillHist(f"{channel}/{syst}/stack/phi", pair1.Phi(), weight, 64, -3.2, 3.2)
            super().FillHist(f"{channel}/{syst}/stack/mass", pair1.M(), weight, 200, 0., 200.)
            super().FillHist(f"{channel}/{syst}/stack/pt", pair2.Pt(), weight, 300, 0., 300.)
            super().FillHist(f"{channel}/{syst}/stack/eta", pair2.Eta(), weight, 100, -5., 5.)
            super().FillHist(f"{channel}/{syst}/stack/phi", pair2.Phi(), weight, 64, -3.2, 3.2)
            super().FillHist(f"{channel}/{syst}/stack/mass", pair2.M(), weight, 200, 0., 200.)
        
        # Fill ZCands
        if "1E2Mu" in channel:
            ZCand = pairs       # pairs == pair
            super().FillHist(f"{channel}/{syst}/ZCand/pt", ZCand.Pt(), weight, 300, 0., 300.)
            super().FillHist(f"{channel}/{syst}/ZCand/eta", ZCand.Eta(), weight, 100, -5., 5.)
            super().FillHist(f"{channel}/{syst}/ZCand/phi", ZCand.Phi(), weight, 64, -3.2, 3.2)
            super().FillHist(f"{channel}/{syst}/ZCand/mass", ZCand.M(), weight, 200, 0., 200.)
        else:
            pair1, pair2 = pairs
            mZ = 91.2
            if abs(pair1.M() - mZ) < abs(pair2.M() - mZ): ZCand, nZCand = pair1, pair2
            else:                                         ZCand, nZCand = pair2, pair1
            super().FillHist(f"{channel}/{syst}/ZCand/pt", ZCand.Pt(), weight, 300, 0., 300.)
            super().FillHist(f"{channel}/{syst}/ZCand/eta", ZCand.Eta(), weight, 100, -5., 5.)
            super().FillHist(f"{channel}/{syst}/ZCand/phi", ZCand.Phi(), weight, 64, -3.2, 3.2)
            super().FillHist(f"{channel}/{syst}/ZCand/mass", ZCand.M(), weight, 200, 0., 200.)
            super().FillHist(f"{channel}/{syst}/nZCand/pt", nZCand.Pt(), weight, 300, 0., 300.)
            super().FillHist(f"{channel}/{syst}/nZCand/eta", nZCand.Eta(), weight, 100, -5., 5.)
            super().FillHist(f"{channel}/{syst}/nZCand/phi", nZCand.Phi(), weight, 64, -3.2, 3.2)
            super().FillHist(f"{channel}/{syst}/nZCand/mass", nZCand.M(), weight, 200, 0., 200.)
         
        # Fill graph inputs
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
        super().FillHist(f"{channel}/{syst}/inputs/METv/energy", METv.E(), weight, 300, 0., 300.)
        super().FillHist(f"{channel}/{syst}/inputs/METv/px", METv.Px(), weight, 500, -250, 250)
        super().FillHist(f"{channel}/{syst}/inputs/METv/py", METv.Py(), weight, 500, -250, 250)
        super().FillHist(f"{channel}/{syst}/inputs/METv/pz", METv.Pz(), weight, 500, -250, 250) 
        super().FillHist(f"{channel}/{syst}/inputs/Nj", jets.size(), weight, 20, 0., 20.)
        super().FillHist(f"{channel}/{syst}/inputs/Nb", bjets.size(), weight, 20, 0., 20.)
        super().FillHist(f"{channel}/{syst}/inputs/MET", METv.Pt(), weight, 300, 0., 300.)
        for idx, lep in enumerate(list(electrons)+list(muons), start=1):
            super().FillHist(f"{channel}/{syst}/inputs/MT{idx}", (lep+METv).Mt(), weight, 300, 0., 300.)
        
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
            
            score_nonprompt = scores[f"{signal}_vs_nonprompt"]
            score_diboson = scores[f"{signal}_vs_diboson"]
            score_ttZ    = scores[f"{signal}_vs_ttZ"]
            super().FillHist(f"{channel}/{syst}/{signal}/score_nonprompt", score_nonprompt, weight, 100, 0., 1.)
            super().FillHist(f"{channel}/{syst}/{signal}/score_diboson", score_diboson, weight, 100, 0., 1.)
            super().FillHist(f"{channel}/{syst}/{signal}/score_ttZ", score_ttZ, weight, 100, 0., 1.)
