from ROOT import gSystem
from ROOT import DiLeptonBase
from ROOT import TString
from ROOT.std import vector
from ROOT.JetTagging import Parameters as jParameters
from ROOT import Muon, Electron, Jet
gSystem.Load("/cvmfs/cms.cern.ch/slc7_amd64_gcc900/external/lhapdf/6.2.3/lib/libLHAPDF.so")

class ClosDiLepTrigs(DiLeptonBase):
    def __init__(self):
        super().__init__()
        # at this point, DiLeptonBase::initializeAnalyzer has not been called
        
    def initializePyAnalyzer(self):
        super().initializeAnalyzer()
        # link flags
        self.channel = None
        self.run_syst = False
        if (super().RunDiMu and super().RunEMu):
            print("Wrong channel flag")
            exit(1)
        if super().RunDiMu:     self.channel = "RunDiMu"
        if super().RunEMu:      self.channel = "RunEMu"
    
    def executeEvent(self):
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
        isDiMu = (len(tightMuons) == 2 and len(vetoMuons) == 2 and \
                  len(tightElectrons) == 0 and len(vetoElectrons) == 0)
        isEMu = (len(tightMuons) == 1 and len(vetoMuons) == 1 and \
                 len(tightElectrons) == 1 and len(vetoElectrons) == 1)
        
        if self.channel == "RunDiMu":
            if not isDiMu: return None
        if self.channel == "RunEMu":
            if not isEMu: return None
        
        if self.channel == "RunDiMu":
            mu1, mu2 = tightMuons.at(0), tightMuons.at(1)
            if not mu1.Pt() > 20.: return None
            if not mu2.Pt() > 10.: return None
        elif self.channel == "RunEMu":
            ele = tightElectrons.at(0)
            mu = tightMuons.at(0)
            leadMu = mu.Pt() > 25. and ele.Pt() > 15.
            leadEle = ele.Pt() > 25. and mu.Pt() > 10.
            if not (leadMu or leadEle): return None
        else:
            print(f"Wrong channel {self.channel}")
            exit(1)
        return self.channel

    def getWeight(self, channel, event, muons, electrons, jets, truth):
        weight = 1.
        
        if not super().IsDATA:
            weight *= super().MCweight()
            weight *= event.GetTriggerLumi("Full")
            w_prefire = super().GetPrefireWeight(0)
            w_pileup = super().GetPileUpWeight(super().nPileUp, 0)
            weight *= w_prefire            # print(f"w_prefire: {w_prefire}")
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
        if channel == "RunDiMu":
            effDZ = self.getDZEfficiency(channel, isDATA=False)
            trigWeight = self.getDblMuTriggerEff(muons, False, 0) * effDZ
            trigWeightUp = self.getDblMuTriggerEff(muons, False, 1) * effDZ
            trigWeightDown = self.getDblMuTriggerEff(muons, False, -1) * effDZ
        elif channel == "RunEMu":
            effDZ = self.getDZEfficiency(channel, isDATA=False)
            trigWeight = self.getEMuTriggerEff(electrons, muons, False, 0) * effDZ
            trigWeightUp = self.getEMuTriggerEff(electrons, muons, False, 1) * effDZ
            trigWeightDown = self.getEMuTriggerEff(electrons, muons, False, -1) * effDZ 
        else:
            print(f"Wrong channel {channel}")
            exit(1)

        # get weight
        super().FillHist("sumweight", 0., weight, 5, 0., 5.)
        super().FillHist("sumweight", 1., weight*trigWeight, 5, 0., 5.)
        super().FillHist("sumweight", 2., weight*trigWeightUp, 5, 0., 5.)
        super().FillHist("sumweight", 3., weight*trigWeightDown, 5, 0., 5.)
        
        if channel == "RunDiMu" and not evt.PassTrigger(super().DblMuTriggers): return None
        if channel == "RunEMu" and not evt.PassTrigger(super().EMuTriggers): return None
        super().FillHist("sumweight", 4., weight, 5, 0., 5.)
