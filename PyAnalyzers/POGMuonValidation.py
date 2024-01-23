import ROOT
from ROOT import gSystem, std
from ROOT import TString
from ROOT import DiLeptonBase
from ROOT import Event, Muon, Electron, Jet
gSystem.Load("/cvmfs/cms.cern.ch/slc7_amd64_gcc900/external/lhapdf/6.2.3/lib/libLHAPDF.so")


class POGMuonValidation(DiLeptonBase):
    def __init__(self):
        super().__init__()
        # at this point, DiLeptonBase::initializeAnalyzer has not been called yet

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

        # only fill events if they are selected
        # do not return current event, we still have to check for the scale varaitions
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

            # no need to fill events that are not selected for current scale variation
            if channel is None: continue
            weight = self.getWeight(ev, objects)
            self.fillObjects(channel, objects, weight, syst=scale)

    def defineObjects(self, rawMuons: std.vector[Muon], rawJets: std.vector[Jet], scale: str = "Central"):
        # first deep copy objects
        allMuons = std.vector[Muon](rawMuons)
        allJets = std.vector[Jet](rawJets)
        assert allMuons is not rawMuons
        assert allJets is not rawJets

        # apply scale variation
        if scale == "MuonEnUp":         allMuons = super().ScaleMuons(allMuons, +1)
        if scale == "MuonEnDown":       allMuons = super().ScaleMuons(allMuons, -1)
        if scale == "JetResUp":         allJets = super().SmearJets(allJets, +1)
        if scale == "JetResDown":       allJets = super().SmearJets(allJets, -1)
        if scale == "JetEnUp":          allJets = super().ScaleJets(allJets, +1)
        if scale == "JetEnDown":        allJets = super().ScaleJets(allJets, -1)
    
        # select objects
        # for muons, in the analysis we use muon starting from 10 GeV
        # here, since POGMedium ID has been optimized from 15 GeV, we use 15 GeV
        muons = super().SelectMuons(allMuons, "POGMediumWithTightIso", 15., 2.4)
        jets = super().SelectJets(allJets, "tight", 20., 2.4)

        # sort objects
        muons = std.vector[Muon](sorted(muons, key=lambda x: x.Pt(), reverse=True))
        jets = std.vector[Jet](sorted(jets, key=lambda x: x.Pt(), reverse=True))

        return (muons, jets)

    def selectEvent(self, ev: Event, objects: dict) -> str:
        muons, jets = objects["muons"], objects["jets"]
        
        if not len(muons) == 2: return None
        if not ev.PassTrigger(super().DblMuTriggers): return None
        mu1, mu2 = muons
        if not mu1.Pt() > 20.: return None
        if not mu1.Charge()+mu2.Charge() == 0: return None
        pair = mu1 + mu2
        if not pair.M() > 50.: return None
        return "DiMu"

    def getWeight(self, ev: Event, objects: dict, syst: str = "Central", condition: str = "") -> float:
        weight = 1.
        if not syst in self.systematics:
            raise NotImplementedError(f"Undefined weight variation for syst {syst}")

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

        # apply weights
        w_muonRecoSF = 1.
        w_muonIDSF = 1.
        w_muonIsoSF = 1.
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

    def fillObjects(self, channel: str, objects: dict, weight: float, syst: str = "Central", condition: str=""):
        muons, jets = objects["muons"], objects["jets"]
        METv = objects["MET"]
        pair = muons.at(0) + muons.at(1)

        # make prefix
        prefix_list = [channel, condition, syst]
        try: prefix_list.remove("")
        except: pass
        prefix = "/".join(prefix_list)
        # fill objects
        for idx, mu in enumerate(muons, start=1):
            super().FillHist(f"{prefix}/muons/{idx}/pt", mu.Pt(), weight, 300, 0., 300.)
            super().FillHist(f"{prefix}/muons/{idx}/eta", mu.Eta(), weight, 48, -2.4, 2.4)
            super().FillHist(f"{prefix}/muons/{idx}/phi", mu.Phi(), weight, 64, -3.2, 3.2)
            super().FillHist(f"{prefix}/muons/{idx}/mass", mu.M(), weight, 10, 0., 1.)
        for idx, jet in enumerate(jets, start=1):
            super().FillHist(f"{prefix}/jets/{idx}/pt", jet.Pt(), weight, 300, 0., 300.)
            super().FillHist(f"{prefix}/jets/{idx}/eta", jet.Eta(), weight, 48, -2.4, 2.4)
            super().FillHist(f"{prefix}/jets/{idx}/phi", jet.Phi(), weight, 64, -3.2, 3.2)
            super().FillHist(f"{prefix}/jets/{idx}/mass", jet.M(), weight, 100, 0., 100.)
        super().FillHist(f"{prefix}/jets/size", jets.size(), weight, 20, 0., 20.)
        super().FillHist(f"{prefix}/METv/pt", METv.Pt(), weight, 300, 0., 300.)
        super().FillHist(f"{prefix}/METv/phi", METv.Phi(), weight, 64, -3.2, 3.2)
        super().FillHist(f"{prefix}/pair/pt", pair.Pt(), weight, 300, 0., 300.)
        super().FillHist(f"{prefix}/pair/eta", pair.Eta(), weight, 100, -5., 5.)
        super().FillHist(f"{prefix}/pair/phi", pair.Phi(), weight, 64, -3.2, 3.2)
        super().FillHist(f"{prefix}/pair/mass", pair.M(), weight, 300, 0., 300.)
