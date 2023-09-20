from ROOT import gSystem
from ROOT import TriLeptonBase
from ROOT import TString
from ROOT.std import vector
from ROOT.JetTagging import Parameters as jParameters
from ROOT import Lepton, Muon, Electron, Jet
gSystem.Load("/cvmfs/cms.cern.ch/slc7_amd64_gcc900/external/lhapdf/6.2.3/lib/libLHAPDF.so")

class ClosTriLepTrigs(TriLeptonBase):
    def __init__(self):
        super().__init__() # at this point, DiLeptonBase::initializeAnalyzer has not been called
        
    def initializePyAnalyzer(self):
        super().initializeAnalyzer()
        # link flags
        self.channel = None
        if (super().Skim1E2Mu and super().Skim3Mu):
            print("Wrong channel flag")
            exit(1)
        if super().Skim3Mu:     self.channel = "Skim3Mu"
        if super().Skim1E2Mu:   self.channel = "Skim1E2Mu"
        
    def executeEvent(self):
        # only central with trigger eff up / down performed
        if not super().PassMETFilter(): return None
        ev = super().GetEvent()
        rawMuons = super().GetAllMuons()
        rawElectrons = super().GetAllElectrons()
        rawJets = super().GetAllJets()
        METv = ev.GetMETVector()
        truth = super().GetGens() if not super().IsDATA else None
        
        # Central scale
        vetoMuons, tightMuons, vetoElectrons, tightElectrons, jets, bjets = self.defineObjects(rawMuons, rawElectrons, rawJets)
        channel = self.selectEvent(ev, vetoMuons, tightMuons, vetoElectrons, tightElectrons, jets, bjets, METv)
        
        if not channel is None:
            objects = {"muons": tightMuons,
                       "electrons": tightElectrons,
                       "jets": jets,
                       "bjets": bjets,
                       "METv": METv
                       }
            weight = self.getWeight(channel, ev, tightMuons, tightElectrons, jets, truth)
            self.FillObjects(channel, ev, objects, weight)
        
    def defineObjects(self, rawMuons, rawElectrons, rawJets):
        # first copy objects
        allMuons = rawMuons
        allElectrons = rawElectrons
        allJets = rawJets
        
        vetoMuons = super().SelectMuons(allMuons, super().MuonIDs[2], 10., 2.4)
        tightMuons = super().SelectMuons(vetoMuons, super().MuonIDs[0], 10., 2.4)
        vetoElectrons = super().SelectElectrons(allElectrons, super().ElectronIDs[2], 10., 2.5)
        tightElectrons = super().SelectElectrons(vetoElectrons, super().ElectronIDs[0], 10., 2.5)
        jets = super().SelectJets(allJets, "tight", 20., 2.4)
        jets = super().JetsVetoLeptonInside(jets, vetoElectrons, vetoMuons, 0.4)
        bjets = vector[Jet]()
        for jet in jets:
            btagScore = jet.GetTaggerResult(3)                  # DeepJet score
            wp = super().mcCorr.GetJetTaggingCutValue(3, 1)     # DeepJet Medium
            if btagScore > wp: bjets.emplace_back(jet)
            
        vetoMuons = vector[Muon](sorted(vetoMuons, key=lambda x: x.Pt(), reverse=True))
        tightMuons = vector[Muon](sorted(tightMuons, key=lambda x: x.Pt(), reverse=True))
        vetoElectrons = vector[Electron](sorted(vetoElectrons, key=lambda x: x.Pt(), reverse=True))
        tightElectrons = vector[Electron](sorted(tightElectrons, key=lambda x: x.Pt(), reverse=True))
        jets = vector[Jet](sorted(jets, key=lambda x: x.Pt(), reverse=True))
        bjets = vector[Jet](sorted(bjets, key=lambda x: x.Pt(), reverse=True))
        
        return (vetoMuons, tightMuons, vetoElectrons, tightElectrons, jets, bjets)

    def selectEvent(self, event, vetoMuons, tightMuons, vetoElectrons, tightElectrons, jets, bjets, METv):
        is1E2Mu = (len(tightMuons) == 2 and len(vetoMuons) == 2 and \
                   len(tightElectrons) == 1 and len(vetoElectrons) == 1)
        is3Mu   = (len(tightMuons) == 3 and len(vetoMuons) == 3 and \
                   len(tightElectrons) == 0 and len(vetoElectrons) == 0)

        ## channel assertion
        if self.channel == "Skim1E2Mu": 
            if not is1E2Mu: return None
        if self.channel == "Skim3Mu":   
            if not is3Mu: return None
        
        ### syncronized baseline selection with the signal region
        channel = ""
        ## 1E2Mu baseline
        ## 1. pass EMuTriggers
        ## 2. Exact 2 tight muons and 1 tight electron, no additional lepton
        ## 3. Exists OS muon pair with mass > 12 GeV
        ## 4. At least two jets
        if self.channel == "Skim1E2Mu":
            #if not event.PassTrigger(super().EMuTriggers): return None
            leptons = vector[Lepton]()
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
        if self.channel == "Skim3Mu":
            #if not event.PassTrigger(super().DblMuTriggers): return None
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
        
        # should not be here
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
            
    def getWeight(self, channel, event, muons, electrons, jets, truth):
        weight = 1.

        if not super().IsDATA:
            weight *= super().MCweight()
            weight *= event.GetTriggerLumi("Full")
            w_prefire = super().GetPrefireWeight(0)
            w_pileup = super().GetPileUpWeight(super().nPileUp, 0)
            weight *= w_prefire
            weight *= w_pileup

            w_topptweight = 1.
            if "TTLL" in super().MCSample or "TTLJ" in super().MCSample:
                w_topptweight = super().mcCorr.GetTopPtReweight(truth)
            weight *= w_topptweight

        return weight
    
    def FillObjects(self, channel, evt, objects, weight):
        muons = objects["muons"]
        electrons = objects["electrons"]
        jets = objects["jets"]
        bjets = objects["bjets"]
        METv = objects["METv"]
        
        trigWeight = 1.
        trigWeightUp = 1.
        trigWeightDown = 1.
        if channel == "SR1E2Mu":
            dzEff = self.getDZEfficiency("EMu", isDATA=False) 
            trigWeight = self.getEMuTriggerEff(electrons, muons, False, 0) * dzEff
            trigWeightUp = self.getEMuTriggerEff(electrons, muons, False, 1) * dzEff
            trigWeightDown = self.getEMuTriggerEff(electrons, muons, False, -1) * dzEff
        elif channel == "SR3Mu":
            dzEff = self.getDZEfficiency("DiMu", isDATA=False) 
            trigWeight = self.getDblMuTriggerEff(muons, False, 0) * dzEff
            trigWeightUp = self.getDblMuTriggerEff(muons, False, 1) * dzEff
            trigWeightDown = self.getDblMuTriggerEff(muons, False, -1) * dzEff
        else:
            raise NameError(f"Wrong channel {channel}") 
        
        # get weight
        super().FillHist("sumweight", 0., weight, 5, 0., 5.)
        super().FillHist("sumweight", 1., weight*trigWeight, 5, 0., 5.)
        super().FillHist("sumweight", 2., weight*trigWeightUp, 5, 0., 5.)
        super().FillHist("sumweight", 3., weight*trigWeightDown, 5, 0., 5.)
        
        if (channel == "SR1E2Mu" and not evt.PassTrigger(super().EMuTriggers)): return None
        if (channel == "SR3Mu" and not evt.PassTrigger(super().DblMuTriggers)): return None
        super().FillHist("sumweight", 4., weight, 5, 0., 5.)