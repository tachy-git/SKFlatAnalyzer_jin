from ROOT import gSystem
from ROOT import TriLeptonBase
from ROOT import TTree, TString
from ROOT.std import vector
from ROOT.JetTagging import Parameters as jParameters
from ROOT import Lepton, Muon, Electron, Jet
gSystem.Load("/cvmfs/cms.cern.ch/slc7_amd64_gcc900/external/lhapdf/6.2.3/lib/libLHAPDF.so")

## Analyzer class for acceptance study
class AcceptanceStudy(TriLeptonBase):
    def __init__(self):
        ## call the constructors for TriLeptonBase -> AnalyzerCore -> SKFlatNtuple
        super().__init__()
        # at this point, TriLeptonBase::initializeAnalyzer has not been activated

    def initializePyAnalyzer(self):
        ## python mode initializeAnalyzer
        super().initializeAnalyzer()

        # channel
        if super().Skim1E2Mu: self.channel = "Skim1E2Mu"
        elif super().Skim3Mu: self.channel = "Skim3Mu"
        else:
            print("Wrong channel")
            exit(1)

    def executeEvent(self):
        # Define objects
        ev = super().GetEvent()
        rawMuons = super().GetAllMuons()
        rawElectrons = super().GetAllElectrons()
        rawJets = super().GetAllJets()
        METv = ev.GetMETVector()
        truth = super().GetGens() if not super().IsDATA else None
        vetoMuons, tightMuons, vetoElectrons, tightElectrons, jets, bjets = self.defineObjects(rawMuons, rawElectrons, rawJets)
        
        # set weights
        weight = super().MCweight() * ev.GetTriggerLumi("Full") 
        weight *= super().GetPrefireWeight(0)
        weight *= super().GetPileUpWeight(super().nPileUp, 0)
        jtp = jParameters(3, 1, 0, 1)    # DeepJet, Medium, incl, mujets
        vjets = vector[Jet]()
        for j in jets: vjets.emplace_back(j)
        weight *= super().mcCorr.GetBTaggingReweight_1a(vjets, jtp)

        ## Event Selection
        super().FillHist("cutflow", 0., weight, 10, 0., 10.)
        
        # MET filter
        if not super().PassMETFilter(): return None
        super().FillHist("cutflow", 1., weight, 10, 0., 10.)

        # check lepton multiplicity
        is3Mu = (len(tightMuons) == 3 and len(vetoMuons) == 3 and \
                len(tightElectrons) == 0 and len(vetoElectrons) == 0)
        is1E2Mu = len(tightMuons) == 2 and len(vetoMuons) == 2 and \
                  len(tightElectrons) == 1 and len(vetoElectrons) == 1
        if self.channel == "Skim1E2Mu":
            if not is1E2Mu: return None
        if self.channel == "Skim3Mu":
            if not is3Mu: return None
        super().FillHist("cutflow", 2., weight, 10, 0., 10.)

        # prompt matching
        promptMuons = vector[Muon]()
        promptElectrons = vector[Electron]()
        for mu in tightMuons:
            if super().GetLeptonType(mu, truth) > 0: promptMuons.emplace_back(mu)
        for ele in tightElectrons:
            if super().GetLeptonType(ele, truth) > 0: promptElectrons.emplace_back(ele)

        if promptMuons.size() != tightMuons.size(): return None
        if promptElectrons.size() != tightElectrons.size(): return None
        super().FillHist("cutflow", 3., weight, 10, 0., 10.)

        # set lepton multiplicity based weights
        # electron ID SF & EMu Trigger SF not applied yet
        for mu in tightMuons:
            weight *= self.getMuonIDSF(mu, 0)
        if self.channel == "Skim3Mu":
            weight *= self.getDblMuTriggerSF(tightMuons, 0)

        ## 1E2Mu Baseline
        ## 1. pass EMuTriggers
        ## 2. Exact 2 tight muons and 1 tight electron, no additional lepton
        ## 3. Exists OS muon pair with mass > 12 GeV
        ## 4. At least two jets
        if self.channel == "Skim1E2Mu":
            # Trigger
            if not ev.PassTrigger(super().EMuTriggers): return None
            leptons = vector[Lepton]()
            for mu in tightMuons: leptons.emplace_back(mu)
            for ele in tightElectrons: leptons.emplace_back(ele)
            mu1, mu2, ele = tuple(leptons)
            passLeadMu = mu1.Pt() > 25. and ele.Pt() > 15.
            passLeadEle = mu1.Pt() > 10. and ele.Pt() > 25.
            passSafeCut = passLeadMu or passLeadEle
            if not passSafeCut: return None
            super().FillHist("cutflow", 4., weight, 10, 0., 10.)
            
            # charge condition
            if not mu1.Charge()+mu2.Charge() == 0: return None
            super().FillHist("cutflow", 5., weight, 10, 0., 10.)

            # mass condition
            pair = self.makePair(tightMuons)
            if not pair.M() > 12.: return None
            super().FillHist("cutflow", 6., weight, 10, 0., 10.)

            # jet multiplicity
            if not jets.size() >= 2: return None
            super().FillHist("cutflow", 7., weight, 10, 0., 10.)

            if not bjets.size() >= 1: return None
            super().FillHist("cutflow", 8., weight, 10, 0., 10.)

        ## 3Mu baseline
        ## 1. pass DblMuTriggers
        ## 2. Exact 3 tight muons, no additional leptons
        ## 3. Exist OS muon pair,
        ## 4. All OS muon pair mass > 12 GeV
        ## 5. At least two jets
        if self.channel == "Skim3Mu":
            # trigger
            if not ev.PassTrigger(super().DblMuTriggers): return None
            mu1, mu2, mu3  = tuple(tightMuons)
            if not mu1.Pt() > 20.: return None
            if not mu2.Pt() > 10.: return None
            if not mu3.Pt() > 10.: return None
            super().FillHist("cutflow", 4., weight, 10, 0., 10.)

            # charge condition
            if not abs(mu1.Charge()+mu2.Charge()+mu3.Charge()) == 1: return None
            super().FillHist("cutflow", 5., weight, 10, 0., 10)

            # mass condition
            pair1, pair2 = self.makePair(tightMuons)
            if not pair1.M() > 12.: return None
            if not pair2.M() > 12.: return None
            super().FillHist("cutflow", 6., weight, 10, 0., 10)

            # jet multiplicity
            if not jets.size() >= 2: return None
            super().FillHist("cutflow", 7., weight, 10, 0., 10)
            
            if not bjets.size() >= 1: return None
            super().FillHist("cutflow", 8., weight, 10, 0., 10)
        
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
            print(f"[AcceptanceStudy::makePair] wrong no. of muons {muons.size}")
