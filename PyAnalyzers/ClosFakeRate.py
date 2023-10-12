## Fake rate closure
## Using TTLL sample, compare observed and expected M(mu+mu-) distribution in the signal region
## and the difference will be the systematic variation
from ROOT import gSystem
from ROOT import TriLeptonBase
from ROOT import TString
from ROOT.std import vector
from ROOT.JetTagging import Parameters as jParameters
from ROOT import Lepton, Muon , Electron, Jet
gSystem.Load("/cvmfs/cms.cern.ch/slc7_amd64_gcc900/external/lhapdf/6.2.3/lib/libLHAPDF.so")


class ClosFakeRate(TriLeptonBase):
    def __init__(self):
        super().__init__()  # at this point, TriLeptonBase::initializeAnalyzer has not been called
        
    def initializePyAnalyzer(self):
        super().initializeAnalyzer()
        self.channel = ""
        if super().Skim1E2Mu:   self.channel = "Skim1E2Mu"
        if super().Skim3Mu:     self.channel = "Skim3Mu"
        if not self.channel in ["Skim1E2Mu", "Skim3Mu"]:
            raise KeyError(f"Wrong channel {self.channel}")
        
    def executeEvent(self):
        if not super().PassMETFilter(): return None
        ev = super().GetEvent()
        rawMuons = super().GetAllMuons()
        rawElectrons = super().GetAllElectrons()
        rawJets = super().GetAllJets()
        METv = ev.GetMETVector()
        truth = super().GetGens()
        
        # Get central scale objects & event selection
        vetoMuons, looseMuons, tightMuons, vetoElectrons, looseElectrons, tightElectrons, jets, bjets = self.defineObjects(rawMuons, rawElectrons, rawJets)
        thisChannel = self.selectEvent(ev, truth, vetoMuons, looseMuons, tightMuons, vetoElectrons, looseElectrons, tightElectrons, jets, bjets, METv)
        
        if thisChannel is None: return None
        pairs = self.makePair(looseMuons)
        objects = {"muons": looseMuons,
                   "electrons": looseElectrons,
                   "jets": jets,
                   "bjets": bjets,
                   "METv": METv,
                   "pairs": pairs}
        
        weight = super().MCweight() * ev.GetTriggerLumi("Full") * super().GetPrefireWeight(0) * super().GetPileUpWeight(super().nPileUp, 0)

        if "SR" in thisChannel:
            # all leptons are tight
            self.fillObjects(thisChannel, objects, weight, syst="Central")
        else:
            # not all leptons are tight
            fakeWeight = super().getFakeWeight(looseMuons, looseElectrons, 0)
            fakeWeightUp = super().getFakeWeight(looseMuons, looseElectrons, 1)
            fakeWeightDown = super().getFakeWeight(looseMuons, looseElectrons, -1)
            
            self.fillObjects(thisChannel, objects, weight*fakeWeight, syst="Central")
            self.fillObjects(thisChannel, objects, weight*fakeWeightUp, syst="NonpromptUp")
            self.fillObjects(thisChannel, objects, weight*fakeWeightDown, syst="NonpromptDown") 

    def fillObjects(self, channel, objects, weight, syst):
        muons = objects["muons"]
        electrons = objects["electrons"]
        jets = objects["jets"]
        bjets = objects["bjets"]
        METv = objects["METv"]
        pairs = objects["pairs"]

        for idx, mu in enumerate(muons, start=1):
            super().FillHist(f"{channel}/{syst}/muons/{idx}/pt", mu.Pt(), weight, 300, 0., 300.)
            super().FillHist(f"{channel}/{syst}/muons/{idx}/eta", mu.Eta(), weight, 48, -2.4, 2.4)
            super().FillHist(f"{channel}/{syst}/muons/{idx}/phi", mu.Phi(), weight, 64, -3.2, 3.2)
            super().FillHist(f"{channel}/{syst}/muons/{idx}/mass", mu.M(), weight, 10, 0., 1.)
        for idx, ele in enumerate(electrons, start=1):
            super().FillHist(f"{channel}/{syst}/electrons/{idx}/pt", ele.Pt(), weight, 300, 0., 300.)
            super().FillHist(f"{channel}/{syst}/electrons/{idx}/scEta", ele.scEta(), weight, 50, -2.5, 2.5)
            super().FillHist(f"{channel}/{syst}/electrons/{idx}/phi", ele.Phi(), weight, 64, -3.2, 3.2)
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
        
        if "1E2Mu" in channel:
            ZCand = pairs
            super().FillHist(f"{channel}/{syst}/ZCand/pt", ZCand.Pt(), weight, 300, 0., 300.)
            super().FillHist(f"{channel}/{syst}/ZCand/eta", ZCand.Eta(), weight, 100, -5., 5.)
            super().FillHist(f"{channel}/{syst}/ZCand/phi", ZCand.Phi(), weight, 64, -3.2, 3.2)
            super().FillHist(f"{channel}/{syst}/ZCand/mass", ZCand.M(), weight, 200, 0., 200.)
        else:
            pair1, pair2 = pairs
            super().FillHist(f"{channel}/{syst}/ZCand1/pt", pair1.Pt(), weight, 300, 0., 300.)
            super().FillHist(f"{channel}/{syst}/ZCand1/eta", pair1.Eta(), weight, 100, -5., 5.)
            super().FillHist(f"{channel}/{syst}/ZCand1/phi", pair1.Phi(), weight, 64, -3.2, 3.2)
            super().FillHist(f"{channel}/{syst}/ZCand1/mass", pair1.M(), weight, 200, 0., 200.)
            super().FillHist(f"{channel}/{syst}/ZCand2/pt", pair2.Pt(), weight, 300, 0., 300.)
            super().FillHist(f"{channel}/{syst}/ZCand2/eta", pair2.Eta(), weight, 100, -5., 5.)
            super().FillHist(f"{channel}/{syst}/ZCand2/phi", pair2.Phi(), weight, 64, -3.2, 3.2)
            super().FillHist(f"{channel}/{syst}/ZCand2/mass", pair2.M(), weight, 200, 0., 200.)
            
    def defineObjects(self, rawMuons, rawElectrons, rawJets):
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
            btagScore = jet.GetTaggerResult(0)                  # DeepJet score 3
            wp = super().mcCorr.GetJetTaggingCutValue(0, 2)     # DeepJet Medium 3 1
            if btagScore > wp: bjets.emplace_back(jet)
            
        vetoMuons = vector[Muon](sorted(vetoMuons, key=lambda x: x.Pt(), reverse=True))
        tightMuons = vector[Muon](sorted(tightMuons, key=lambda x: x.Pt(), reverse=True))
        vetoElectrons = vector[Electron](sorted(vetoElectrons, key=lambda x: x.Pt(), reverse=True))
        tightElectrons = vector[Electron](sorted(tightElectrons, key=lambda x: x.Pt(), reverse=True))
        jets = vector[Jet](sorted(jets, key=lambda x: x.Pt(), reverse=True))
        bjets = vector[Jet](sorted(bjets, key=lambda x: x.Pt(), reverse=True))
        
        return (vetoMuons, looseMuons, tightMuons, vetoElectrons, looseElectrons, tightElectrons, jets, bjets)
                
    def selectEvent(self, event, truth, vetoMuons, looseMuons, tightMuons, vetoElectrons, looseElectrons, tightElectrons, jets, bjets, METv):
        is3Mu = (len(looseMuons) == 3 and len(vetoMuons) == 3 and len(looseElectrons) == 0 and len(vetoElectrons) == 0)
        is1E2Mu = (len(looseMuons) == 2 and len(vetoMuons) == 2 and len(looseElectrons) == 1 and len(vetoElectrons) == 1)
        
        if self.channel == "Skim1E2Mu":
            if not is1E2Mu: return None
        if self.channel == "Skim3Mu":
            if not is3Mu: return None
            
        # prompt matching
        promptMuons = vector[Muon]()
        promptElectrons = vector[Electron]()
        for mu in tightMuons:
            if super().GetLeptonType(mu, truth) > 0: promptMuons.emplace_back(mu)
        for ele in tightElectrons:
            if super().GetLeptonType(ele, truth) > 0: promptElectrons.emplace_back(ele)
        
        # not all leptons are prompt
        # if len(promptMuons) == len(looseMuons) and len(promptElectrons) == len(looseElectrons): return None
        
        ## 1E2Mu baseline
        ## 1. pass EMuTriggers
        ## 2. Exact 2 tight muons and 1 tight electron, no additional lepton
        ## 3. Exists OS muon pair with mass > 12 GeV
        ## 4. At least two jets and one b-jet
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
            
            if len(looseMuons) == len(tightMuons) and len(looseElectrons) == len(tightElectrons):
                return "SR1E2Mu"
            else:
                return "SB1E2Mu"
        ## 3Mu baseline
        ## 1. pass DblMuTriggers
        ## 2. Exact 3 tight muons, no additional leptons
        ## 3. Exist OS muon pair,
        ## 4. All OS muon pair mass > 12 GeV
        ## 5. At least two jets and one b-jet
        if self.channel == "Skim3Mu":
            if not event.PassTrigger(super().DblMuTriggers): return None
            mu1, mu2, mu3  = tuple(looseMuons)
            if not mu1.Pt() > 20.: return None
            if not mu2.Pt() > 10.: return None
            
            if not abs(mu1.Charge()+mu2.Charge()+mu3.Charge()) == 1: return None
            pair1, pair2 = self.makePair(looseMuons)
            if not pair1.M() > 12.: return None
            if not pair2.M() > 12.: return None
            if not jets.size() >= 2: return None
            if not bjets.size() >= 1: return None
            
            if len(looseMuons) == len(tightMuons) and len(looseElectrons) == len(tightElectrons):
                return "SR3Mu"
            else:
                return "SB3Mu"
        
        # shoule not be here
        raise EOFError()
    
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
            raise NotImplementedError(f"wrong no. of muons {muons.size()}")
