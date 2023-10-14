from ROOT import gSystem
from ROOT import TriLeptonBase
from ROOT import TTree, TString
from ROOT.std import vector
from ROOT import Lepton, Muon, Electron, Jet
gSystem.Load("/cvmfs/cms.cern.ch/slc7_amd64_gcc900/external/lhapdf/6.2.3/lib/libLHAPDF.so")

from array import array
from itertools import product
from MLTools.helpers import loadModels
from MLTools.helpers import getGraphInput, getGraphScore

class MatrixUnbinned(TriLeptonBase):
    def __init__(self):
        super().__init__()
        
    def initializePyAnalyzer(self):
        super().initializeAnalyzer()
        
        ## channel assertion
        try: assert super().Skim1E2Mu or super().Skim3Mu
        except: raise AssertionError(f"Wrong channel flag")
        
        if super().Skim1E2Mu: self.channel = "Skim1E2Mu"
        if super().Skim3Mu: self.channel = "Skim3Mu"
        self.network = "GraphNeuralNet"
        
        self.systematics = ["Central", "NonpromptUp", "NonpromptDown"]

        self.signalStrings = ["MHc-160_MA-85", "MHc-130_MA-90", "MHc-100_MA-95"]
        self.backgroundStrings = ["nonprompt", "diboson", "ttZ"]

        self.models = loadModels(self.network, self.channel, self.signalStrings, self.backgroundStrings)
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
        thisChannel = self.selectEvent(ev, vetoMuons, looseMuons, tightMuons, vetoElectrons, looseElectrons, tightElectrons, jets, bjets, METv)

        if thisChannel is None: return None
        pairs = self.makePair(looseMuons)
        _, scores = self.evalScore(looseMuons, looseElectrons, jets, bjets, METv)
        
        if thisChannel == "SR1E2Mu":
            for syst in self.systematics:
                self.mass1[syst][0] = pairs.M()
                self.mass2[syst][0] = -999.
        else:
            for syst in self.systematics:
                self.mass1[syst][0] = pairs[0].M()
                self.mass2[syst][0] = pairs[1].M()
                
        for SIG in self.signalStrings:
            for syst in self.systematics:
                self.scoreX[f"{SIG}_{syst}"][0] = scores[f"{SIG}_vs_nonprompt"]
                self.scoreY[f"{SIG}_{syst}"][0] = scores[f"{SIG}_vs_diboson"]
                self.scoreZ[f"{SIG}_{syst}"][0] = scores[f"{SIG}_vs_ttZ"]
        
        self.weight["Central"][0] = super().getFakeWeight(looseMuons, looseElectrons, 0)
        self.weight["NonpromptUp"][0] = super().getFakeWeight(looseMuons, looseElectrons, 1)
        self.weight["NonpromptDown"][0] = super().getFakeWeight(looseMuons, looseElectrons, -1)
        
        for syst in self.systematics:
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
                self.scoreZ[f"{SIG}_{syst}"][0] = -999.
            self.weight[syst][0] = -999.
            
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

    def selectEvent(self, event, vetoMuons, looseMuons, tightMuons, vetoElectrons, looseElectrons, tightElectrons, jets, bjets, METv):
        is3Mu = (looseMuons.size() == 3 and vetoMuons.size() == 3 and \
                 looseElectrons.size() == 0 and vetoElectrons.size() == 0)
        is1E2Mu = (looseMuons.size() == 2 and vetoMuons.size() == 2 and \
                   looseElectrons.size() == 1 and vetoElectrons.size() == 1)
        
        if self.channel == "Skim1E2Mu":
            if not is1E2Mu: return None
            if (tightMuons.size() == looseMuons.size()) and (tightElectrons.size() == looseElectrons.size()): return None
        if self.channel == "Skim3Mu":
            if not is3Mu: return None
            if tightMuons.size() == looseMuons.size(): return None
            
        ## 1E2Mu baseline
        ## 1. pass EMuTriggers
        ## 2. Exact 2 tight muons and 1 tight electron, no additional lepton
        ## 3. Exists OS muon pair with mass > 12 GeV
        ## 4. At least two jets, at least one b-jet
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
