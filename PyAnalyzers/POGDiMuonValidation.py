import ROOT
from ROOT import gSystem
from ROOT.std import vector
from ROOT import TString
from ROOT import DiLeptonBase
from ROOT import Event, Muon, Jet
gSystem.Load("/cvmfs/cms.cern.ch/slc7_amd64_gcc900/external/lhapdf/6.2.3/lib/libLHAPDF.so")

class POGDiMuonValidation(DiLeptonBase):
    def __init__(self):
        super().__init__()
        # at this point, POGDiMuonBase::initializeAnalyzer has not been called
        
    def initializePyAnalyzer(self):
        super().initializeAnalyzer()
        
        # link flags
        self.run_syst = True if super().RunSyst else False
        
        self.weightVariations = ["Central"]
        self.scaleVariations = []
        if self.run_syst:
            self.weightVariations += ["L1PrefireUp", "L1PrefireDown",
                                      "PileupReweightUp", "PileupReweightDown",
                                      "MuonIDSFUp", "MuonIDSFDown",
                                      "MuonIsoSFUp", "MuonIsoSFDown"]
            self.scaleVariations += ["JetResUp", "JetResDown", 
                                     "JetEnUp", "JetEnDown",
                                     "MuonEnUp", "MuonEnDown"]
        self.systematics = self.weightVariations + self.scaleVariations
    
    def executeEvent(self):
        if not super().PassMETFilter(): return None
        ev = super().GetEvent()
        rawMuons = super().GetAllMuons()
        rawJets = super().GetAllJets()
        METv = ev.GetMETVector()
        truth = super().GetGens()
        
        # Central scale
        muons, jets = self.defineObjects(rawMuons, rawJets)
        objects = {"muons": muons, "jets": jets, "MET": METv, "truth": truth}
        channel = self.selectEvent(ev, objects)
        
        if channel is not None:
            weight = self.getWeight(ev, objects, condition="NoIDSF")
            self.fillObjects(channel, objects, weight, condition="NoIDSF")
            weight = self.getWeight(ev, objects, condition="NoIsoSF")
            self.fillObjects(channel, objects, weight, condition="NoIsoSF")
            for syst in self.weightVariations:
                weight = self.getWeight(ev, objects, syst=syst)
                self.fillObjects(channel, objects, weight, syst=syst)
        
        # Scale Variations
        for scale in self.scaleVariations:
            muons, jets = self.defineObjects(rawMuons, rawJets, scale=scale)
            objects = {"muons": muons, "jets": jets, "MET": METv, "truth": truth}
            channel = self.selectEvent(ev, objects)
            
            if channel is None: continue
            weight = self.getWeight(ev, objects)
            self.fillObjects(channel, objects, weight, syst=scale)
            
    def defineObjects(self, rawMuons: vector[Muon],
                            rawJets: vector[Jet],
                            scale: str="Central"):
        # first copy objects
        allMuons = vector[Muon](rawMuons);  assert allMuons is not rawMuons 
        allJets = vector[Jet](rawJets);     assert allJets is not rawJets
        
        # check the scale argument
        if scale == "MuonEnUp":         allMuons = super().ScaleMuons(allMuons, +1)
        if scale == "MuonEnDown":       allMuons = super().ScaleMuons(allMuons, -1)
        if scale == "JetResUp":         allJets = super().SmearJets(allJets, +1)
        if scale == "JetResDown":       allJets = super().SmearJets(allJets, -1)
        if scale == "JetEnUp":          allJets = super().ScaleJets(allJets, +1)
        if scale == "JetEnDown":        allJets = super().ScaleJets(allJets, -1)
        
        muons = super().SelectMuons(allMuons, "POGMediumWithTightIso", 15., 2.4)
        jets = super().SelectJets(allJets, "tight", 20., 2.4)
        
        muons = vector[Muon](sorted(muons, key=lambda x: x.Pt(), reverse=True))
        jets = vector[Jet](sorted(jets, key=lambda x: x.Pt(), reverse=True))
        
        return (muons, jets)

    def selectEvent(self, ev: Event, objects: dict) -> str:
        muons, jets = objects["muons"], objects["jets"]
        
        if not len(muons) == 2: return None
        if not ev.PassTrigger(super().DblMuTriggers): return None
        mu1, mu2 = tuple(muons)
        if not mu1.Pt() > 20.: return None
        if not mu1.Charge()+mu2.Charge() == 0: return None
        pair = mu1 + mu2
        if not pair.M() > 50.: return None
        return "DiMu"
    
    def getWeight(self, ev: Event, objects: dict, condition: str="fullWeight", syst: str="Central") -> float:
        # scale should be defined in getObject level
        weight = 1.
        if not syst in self.weightVariations:
            raise NotImplementedError(f"Wrong weight variation for syst {syst}")
        
        if super().IsDATA: return weight
        
        muons, jets = objects["muons"], objects["jets"]
        truth = objects["truth"]
        
        weight *= super().MCweight()
        weight *= ev.GetTriggerLumi("Full")
        if syst == "L1PrefireUp":     w_prefire = super().GetPrefireWeight(1)
        elif syst == "L1PrefireDown": w_prefire = super().GetPrefireWeight(-1)
        else:                         w_prefire = super().GetPrefireWeight(0)
            
        if syst == "PileupReweightUp":     w_pileup = super().GetPileUpWeight(super().nPileUp, 1)
        elif syst == "PileupReweightDown": w_pileup = super().GetPileUpWeight(super().nPileUp, -1)
        else:                              w_pileup = super().GetPileUpWeight(super().nPileUp, 0)
        weight *= (w_prefire*w_pileup)
        
        w_zptweight = 1.
        w_topptweight = 1.
        if "TTLL" in super().MCSample or "TTLJ" in super().MCSample:
            w_topptweight = super().mcCorr.GetTopPtReweight(truth)
        weight *= (w_zptweight*w_topptweight)
        
        if condition == "NoIDSF": return weight
        w_muonRecoSF = 1.
        w_muonIDSF = 1.
        w_muonIsoSF = 1.
        
        # apply weights
        for mu in muons:
            w_muonRecoSF *= super().getMuonRecoSF(mu, 0)
            if syst == "MuonIDSFUp":   w_muonIDSF *= super().mcCorr.MuonID_SF("NUM_MediumID_DEN_TrackerMuons", mu.Eta(), mu.Pt(), 1)
            if syst == "MuonIDSFDown": w_muonIDSF *= super().mcCorr.MuonID_SF("NUM_MediumID_DEN_TrackerMuons", mu.Eta(), mu.Pt(), -1)
            else:                      w_muonIDSF *= super().mcCorr.MuonID_SF("NUM_MediumID_DEN_TrackerMuons", mu.Eta(), mu.Pt(), 0)
        weight *= w_muonRecoSF*w_muonIDSF
        if condition == "NoIsoSF": return weight
        
        for mu in muons:
            if syst == "MuonIsoSFUp":   w_muonIsoSF *= super().mcCorr.MuonISO_SF("NUM_TightRelIso_DEN_MediumID", mu.Eta(), mu.Pt(), 1)
            if syst == "MuonIsoSFDown": w_muonIsoSF *= super().mcCorr.MuonISO_SF("NUM_TightRelIso_DEN_MediumID", mu.Eta(), mu.Pt(), -1)
            else:                       w_muonIsoSF *= super().mcCorr.MuonISO_SF("NUM_TightRelIso_DEN_MediumID", mu.Eta(), mu.Pt(), 0)
        weight *= w_muonIsoSF
        return weight
    
    def fillObjects(self, channel: str, objects: dict, weight: float, condition: str="fullWeight", syst: str="Central"):
        muons, jets = objects["muons"], objects["jets"]
        METv = objects["MET"] 
        pair = muons.at(0) + muons.at(1)
        
        ## fill objects
        for idx, mu in enumerate(muons, start=1):
            super().FillHist(f"{channel}/{condition}/{syst}/muons/{idx}/pt", mu.Pt(), weight, 300, 0., 300.)
            super().FillHist(f"{channel}/{condition}/{syst}/muons/{idx}/eta", mu.Eta(), weight, 48, -2.4, 2.4)
            super().FillHist(f"{channel}/{condition}/{syst}/muons/{idx}/phi", mu.Phi(), weight, 64, -3.2, 3.2)
            super().FillHist(f"{channel}/{condition}/{syst}/muons/{idx}/mass", mu.M(), weight, 10, 0., 1.)
        for idx, jet in enumerate(jets, start=1):
            super().FillHist(f"{channel}/{condition}/{syst}/jets/{idx}/pt", jet.Pt(), weight, 300, 0., 300.)
            super().FillHist(f"{channel}/{condition}/{syst}/jets/{idx}/eta", jet.Eta(), weight, 48, -2.4, 2.4)
            super().FillHist(f"{channel}/{condition}/{syst}/jets/{idx}/phi", jet.Phi(), weight, 64, -3.2, 3.2)
            super().FillHist(f"{channel}/{condition}/{syst}/jets/{idx}/mass", jet.M(), weight, 100, 0., 100.)
        super().FillHist(f"{channel}/{condition}/{syst}/jets/size", jets.size(), weight, 20, 0., 20.)
        super().FillHist(f"{channel}/{condition}/{syst}/METv/pt", METv.Pt(), weight, 300, 0., 300.)
        super().FillHist(f"{channel}/{condition}/{syst}/METv/phi", METv.Phi(), weight, 64, -3.2, 3.2)
        super().FillHist(f"{channel}/{condition}/{syst}/pair/pt", pair.Pt(), weight, 300, 0., 300.)
        super().FillHist(f"{channel}/{condition}/{syst}/pair/eta", pair.Eta(), weight, 100, -5., 5.)
        super().FillHist(f"{channel}/{condition}/{syst}/pair/phi", pair.Phi(), weight, 64, -3.2, 3.2)
        super().FillHist(f"{channel}/{condition}/{syst}/pair/mass", pair.M(), weight, 300, 0., 300.)
