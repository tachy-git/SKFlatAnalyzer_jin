from ROOT import gSystem
from ROOT import TriLeptonBase
from ROOT import TString
from ROOT.std import vector
from ROOT.JetTagging import Parameters as jParameters
from ROOT import Lepton, Muon, Electron, Jet
gSystem.Load("/cvmfs/cms.cern.ch/slc7_amd64_gcc900/external/lhapdf/6.2.3/lib/libLHAPDF.so")

class MeasConvMatrix(TriLeptonBase):
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
        objects = {"muons": looseMuons,
                   "electrons": looseElectrons,
                   "jets": jets,
                   "bjets": bjets,
                   "METv": METv,
                   "pairs": pairs,
                    }
        for syst in self.systematics:
            weight = self.getWeight(thisChannel, ev, looseMuons, looseElectrons, jets)
            if syst == "Central":
                weight *= super().getFakeWeight(looseMuons, looseElectrons, 0)
            elif syst == "NonpromptUp":
                weight *= super().getFakeWeight(looseMuons, looseElectrons, 1)
            elif syst == "NonpromptDown":
                weight *= super().getFakeWeight(looseMuons, looseElectrons, -1)
            else:
                print(f"[NonpromptEstimator::executeEvent] Wrong systematics {syst}")
                exit(1)
            self.FillObjects(thisChannel, objects, weight, syst)
            
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
            
        # for conversion measurement
        if "DYJets" in super().MCSample or "ZGToLLG" in super().MCSample:
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
 
        ##### event selection
        ## 1E2Mu ZGamma
        ## 1. pass EMuTriggers
        ## 2. Exact 2 loose muons and 1 loose electron, no additional lepton
        ## 3. Exists OS muon pair with mass > 12 GeV
        ## 4. |M(mumue) - 91.2| < 10
        ## 5. No b-jet
        if self.channel == "Skim1E2Mu":
            if not event.PassTrigger(super().EMuTriggers): return None
            mu1, mu2 = tuple(looseMuons)
            ele = looseElectrons.at(0)
            passLeadMu = mu1.Pt() > 25. and ele.Pt() > 15.
            passLeadEle = mu1.Pt() > 10. and ele.Pt() > 25.
            passSafeCut = passLeadMu or passLeadEle
            if not passSafeCut: return None
            if not mu1.Charge()+mu2.Charge() == 0: return None
            pair = self.makePair(looseMuons)
            isOnZ = abs((mu1+mu2+ele).M() - 91.2) < 10. 
            if not pair.M() > 12.: return None
            if not isOnZ: return None
            if not bjets.size() == 0: return None
            
            return "ZGamma1E2Mu"
        
        ## 3Mu ZGamma
        ## 1. pass DblMuTriggers
        ## 2. Exact 3 tight muons, no additional leptons
        ## 3. Exist OS muon pair,
        ## 4. All OS muon pair mass > 12 GeV
        ## 5. |M(mumumu) - 91.2| < 10
        ## 6. No b-jet
        if self.channel == "Skim3Mu":
            if not event.PassTrigger(super().DblMuTriggers): return None
            mu1, mu2, mu3  = tuple(looseMuons)
            if not mu1.Pt() > 20.: return None
            if not mu2.Pt() > 10.: return None
            if not mu3.Pt() > 10.: return None
            if not abs(mu1.Charge()+mu2.Charge()+mu3.Charge()) == 1: return None
            pair1, pair2 = self.makePair(looseMuons)
            isOnZ = abs((mu1+mu2+mu3).M() - 91.2) < 10.
            if not pair1.M() > 12.: return None
            if not pair2.M() > 12.: return None
            if not isOnZ: return None
            if not bjets.size() == 0: return None
            
            return "ZGamma3Mu"
        
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
        
    def getWeight(self, channel, event, muons, electrons, jets, syst="Central"):
        # No ID and trigger SF applied since leptons reside in ID sideband
        weight = 1.
        if not super().IsDATA:
            weight *= super().MCweight() * super().GetKFactor()
            weight *= event.GetTriggerLumi("Full")
            weight *= super().GetPrefireWeight(0)
            weight *= super().GetPileUpWeight(super().nPileUp, 0)

            # b-tagging
            jtp = jParameters(3, 1, 0, 1)    # DeepJet, Medium, incl, mujets
            vjets = vector[Jet]()
            for j in jets: vjets.emplace_back(j)
            weight *= super().mcCorr.GetBTaggingReweight_1a(vjets, jtp)
        return weight

    def getConvLepton(self, electrons, muons):
        if electrons.size() == 1 and muons.size() == 2:
            return electrons.at(0);
        elif muons.size() == 3:
            mu1, mu2, mu3 = tuple(muons)
            if mu1.Charge() == mu2.Charge():
                pair1, pair2 = mu1+mu3, mu2+mu3
                if abs(pair1.M() - 91.2) > abs(pair2.M() - 91.2):
                    return mu1
                else:
                    return mu2
            elif mu1.Charge() == mu3.Charge():
                pair1, pair2 = mu1+mu2, mu2+mu3
                if abs(pair1.M() - 91.2) > abs(pair2.M() - 91.2):
                    return mu1
                else:
                    return mu3
            else:   # mu2.Charge() == mu3.Charge()
                pair1, pair2 = mu1+mu2, mu1+mu3
                if abs(pair1.M() - 91.2) > abs(pair2.M() - 91.2):
                    return mu2
                else:
                    return mu3
        else:
            raise NotImplementedError(f"Wrong number of muons (muons.size())")
        
    def FillObjects(self, channel, objects, weight, syst):
        muons = objects["muons"]
        electrons = objects["electrons"]
        jets = objects["jets"]
        bjets = objects["bjets"]
        METv = objects["METv"]
        
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
        
        ## fill ZCand
        ZCand = None
        if "1E2Mu" in channel: ZCand = electrons.at(0) + muons.at(0) + muons.at(1)
        if "3Mu" in channel:   ZCand = muons.at(0) + muons.at(1) + muons.at(2)
        convLep = self.getConvLepton(electrons, muons)
        
        super().FillHist(f"{channel}/{syst}/ZCand/pt", ZCand.Pt(), weight, 300, 0., 300.)
        super().FillHist(f"{channel}/{syst}/ZCand/eta", ZCand.Eta(), weight, 100, -5., 5.)
        super().FillHist(f"{channel}/{syst}/ZCand/phi", ZCand.Phi(), weight, 64, -3.2, 3.2)
        super().FillHist(f"{channel}/{syst}/ZCand/mass", ZCand.M(), weight, 200, 0., 200.)
        super().FillHist(f"{channel}/{syst}/convLep/pt", convLep.Pt(), weight, 300, 0., 300.)
        super().FillHist(f"{channel}/{syst}/convLep/eta", convLep.Eta(), weight, 50, -2.5, 2.5)
        super().FillHist(f"{channel}/{syst}/convLep/phi", convLep.Phi(), weight, 64, -3.2, 3.2)
