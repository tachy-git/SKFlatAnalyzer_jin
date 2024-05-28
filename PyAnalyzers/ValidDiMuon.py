import ROOT
from ROOT import gSystem, std
from ROOT import DiLeptonBase
from ROOT.JetTagging import Parameters as jParameters
from ROOT import Event, Muon, Electron, Jet
gSystem.Load("/cvmfs/cms.cern.ch/slc7_amd64_gcc900/external/lhapdf/6.2.3/lib/libLHAPDF.so")

class Valid_POGMuon(DiLeptonBase):
    def __init__(self):
        super().__init__()
        # at this point, DiLeptonBase::initializeAnalyzer has not been called
    
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
                                      "DblMuTrigSFUp", "DblMuTrigSFDown"]
            self.scaleVariations += ["JetResUp", "JetResDown", 
                                     "JetEnUp", "JetEnDown",
                                     "ElectronResUp", "ElectronResDown", 
                                     "ElectronEnUp", "ElectronEnDown",
                                     "MuonEnUp", "MuonEnDown"]
        self.systematics = self.weightVariations + self.scaleVariations
        
    def executeEvent(self):
        if not super().PassMETFilter(): return None
        ev = super().GetEvent()
        rawMuons = super().GetAllMuons()
        rawElectrons = super().GetAllElectrons()
        rawJets = super().GetAllJets()
        METv = ev.GetMETVector()
        truth = super().GetGens()
        
        # Central Scale
        vetoMuons, tightMuons, vetoElectrons, tightElectrons, jets, bjets = self.defineObjects(rawMuons, rawElectrons, rawJets)
        objects = {"tightMuons": tightMuons,
                   "vetoMuons": vetoMuons,
                   "tightElectrons": tightElectrons,
                   "vetoElectrons": vetoElectrons,
                   "jets": jets,
                   "bjets": bjets,
                   "METv": METv,
                   "truth": truth}
        channel = self.selectEvent(ev, objects)
        if channel is not None:
            weight = self.getWeight(ev, objects, condition="beforeIDSF")
            self.fillObjects(channel, objects, weight, condition="beforeIDSF")
            weight = self.getWeight(ev, objects, condition="beforeIsoSF")
            self.fillObjects(channel, objects, weight, condition="beforeIsoSF")
            weight = self.getWeight(ev, objects, condition="beforeTrigSF")
            self.fillObjects(channel, objects, weight, condition="beforeTrigSF")
            for syst in self.weightVariations:
                weight = self.getWeight(ev, objects, syst=syst)
                self.fillObjects(channel, objects, weight, syst=syst)
        
        # Scale Variations
        for scale in self.scaleVariations:
            vetoMuons, tightMuons, vetoElectrons, tightElectrons, jets, bjets = self.defineObjects(rawMuons, rawElectrons, rawJets, scale=scale)
            objects = {"tightMuons": tightMuons,
                       "vetoMuons": vetoMuons,
                       "tightElectrons": tightElectrons,
                       "vetoElectrons": vetoElectrons,
                       "jets": jets,
                       "bjets": bjets,
                       "METv": METv,
                       "truth": truth} 
            channel = self.selectEvent(ev, objects)
            if channel is None: continue
            weight = self.getWeight(ev, objects)
            self.fillObjects(channel, objects, weight, syst=scale)
            
    def defineObjects(self, rawMuons: std.vector[Muon], 
                            rawElectrons: std.vector[Electron],
                            rawJets: std.vector[Jet],
                            scale: str="Central"):
        # first copy objects
        allMuons = std.vector[Muon](rawMuons)
        allElectrons = std.vector[Electron](rawElectrons)
        allJets = std.vector[Jet](rawJets)
        # make sure all objects are deep-copied
        assert allMuons is not rawMuons
        assert allElectrons is not rawElectrons
        assert allJets is not rawJets
        
        # check the scale argument
        if scale == "MuonEnUp":         allMuons = super().ScaleMuons(allMuons, +1)
        if scale == "MuonEnDown":       allMuons = super().ScaleMuons(allMuons, -1)
        if scale == "ElectronResUp":    allElectrons = super().SmearElectrons(allElectrons, +1)
        if scale == "ElectronsResDown": allElectrons = super().SmearElectrons(allElectrons, -1)
        if scale == "ElectronEnUp":     allElectrons = super().ScaleElectrons(allElectrons, +1)
        if scale == "ElectronEnDown":   allElectrons = super().ScaleElectrons(allElectrons, -1)
        if scale == "JetResUp":         allJets = super().SmearJets(allJets, +1)
        if scale == "JetResDown":       allJets = super().SmearJets(allJets, -1)
        if scale == "JetEnUp":          allJets = super().ScaleJets(allJets, +1)
        if scale == "JetEnDown":        allJets = super().ScaleJets(allJets, -1)
        
        vetoMuons = super().SelectMuons(allMuons, "POGMediumWithLooseIso", 15., 2.4)
        tightMuons = super().SelectMuons(allMuons, "POGMediumWithTightIso", 15., 2.4)
        vetoElectrons = super().SelectElectrons(allElectrons, self.ElectronIDs[2], 10., 2.5)
        tightElectrons = super().SelectElectrons(allElectrons, self.ElectronIDs[0], 10., 2.5)
        jets = super().SelectJets(allJets, "tight", 20., 2.4)
        jets = super().JetsVetoLeptonInside(jets, vetoElectrons, vetoMuons, 0.4)
        
        vetoMuons = std.vector[Muon](sorted(vetoMuons, key=lambda x: x.Pt(), reverse=True))
        tightMuons = std.vector[Muon](sorted(tightMuons, key=lambda x: x.Pt(), reverse=True))
        vetoElectrons = std.vector[Electron](sorted(vetoElectrons, key=lambda x: x.Pt(), reverse=True))
        tightElectrons = std.vector[Electron](sorted(tightElectrons, key=lambda x: x.Pt(), reverse=True))
        jets = std.vector[Jet](sorted(jets, key=lambda x: x.Pt(), reverse=True))
        # b-tagging - DeepJet Medium WP (3, 1)
        wp = super().mcCorr.GetJetTaggingCutValue(3, 1)
        tagged_result = [j for j in jets if j.GetTaggerResult(3) > wp]
        bjets = std.vector[Jet](tagged_result)
        
        return (vetoMuons, tightMuons, vetoElectrons, tightElectrons, jets, bjets)
    
    def selectEvent(self, ev: Event, objects: dict) -> str:
        vetoMuons, tightMuons = objects["vetoMuons"], objects["tightMuons"]
        vetoElectrons, tightElectrons = objects["vetoElectrons"], objects["tightElectrons"]
        #jets, bjets = objects["jets"], objects["bjets"]
        #METv = objects["METv"]
        #truth = objects["truth"]
        
        isDiMu = (len(tightMuons) == 2 and len(vetoMuons) == 2 and len(tightElectrons) == 0 and len(vetoElectrons) == 0)
        if not isDiMu: return None
        
        # event selection
        if not ev.PassTrigger(super().isoMuTriggerName): return None
        mu1, mu2 = tuple(tightMuons)
        if not mu1.Pt() > super().triggerSafePtCut: return None
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
        
        tightMuons, tightElectrons = objects["tightMuons"], objects["tightElectrons"]
        #jets, bjets = objects["jets"], objects["bjets"]
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
        
        if condition == "beforeIDSF": return weight
        w_muonRecoSF = 1.
        w_muonIDSF = 1.
        w_muonIsoSF = 1.
        w_trigSF = 1.
        
        # apply weights
        for mu in tightMuons:
            w_muonRecoSF *= super().getMuonRecoSF(mu, 0)
            
        weight *= (w_muonRecoSF*w_muonIDSF)
        if condition == "beforeIsoSF": return weight
        
        weight *= w_muonIsoSF
        if condition == "beforeTrigSF": return weight
        
        weight *= w_trigSF
        return weight
    
    def fillObjects(self, channel: str, objects: dict, weight: float, condition: str="fullWeight", syst: str="Central"):
        muons, electrons = objects["tightMuons"], objects["tightElectrons"]
        jets, bjets = objects["jets"], objects["bjets"]
        METv = objects["METv"] 
        pair = muons.at(0) + muons.at(1)
        
        ## fill objects
        for idx, mu in enumerate(muons, start=1):
            super().FillHist(f"{channel}/{condition}/{syst}/muons/{idx}/pt", mu.Pt(), weight, 300, 0., 300.)
            super().FillHist(f"{channel}/{condition}/{syst}/muons/{idx}/eta", mu.Eta(), weight, 48, -2.4, 2.4)
            super().FillHist(f"{channel}/{condition}/{syst}/muons/{idx}/phi", mu.Phi(), weight, 64, -3.2, 3.2)
            super().FillHist(f"{channel}/{condition}/{syst}/muons/{idx}/mass", mu.M(), weight, 10, 0., 1.)
        for idx, ele in enumerate(electrons, start=1):
            super().FillHist(f"{channel}/{condition}/{syst}/electrons/{idx}/pt", ele.Pt(), weight, 300, 0., 300.)
            super().FillHist(f"{channel}/{condition}/{syst}/electrons/{idx}/eta", ele.Eta(), weight, 50, -2.5, 2.5)
            super().FillHist(f"{channel}/{condition}/{syst}/electrons/{idx}/phi", ele.Phi(), weight, 64, -3.2, 3.2)
            super().FillHist(f"{channel}/{condition}/{syst}/electrons/{idx}/mass", ele.M(), weight, 100, 0., 1.)
        for idx, jet in enumerate(jets, start=1):
            super().FillHist(f"{channel}/{condition}/{syst}/jets/{idx}/pt", jet.Pt(), weight, 300, 0., 300.)
            super().FillHist(f"{channel}/{condition}/{syst}/jets/{idx}/eta", jet.Eta(), weight, 48, -2.4, 2.4)
            super().FillHist(f"{channel}/{condition}/{syst}/jets/{idx}/phi", jet.Phi(), weight, 64, -3.2, 3.2)
            super().FillHist(f"{channel}/{condition}/{syst}/jets/{idx}/mass", jet.M(), weight, 100, 0., 100.)
        super().FillHist(f"{channel}/{condition}/{syst}/jets/size", jets.size(), weight, 20, 0., 20.)
        for idx, bjet in enumerate(bjets, start=1):
            super().FillHist(f"{channel}/{condition}/{syst}/bjets/{idx}/pt", bjet.Pt(), weight, 300, 0., 300.)
            super().FillHist(f"{channel}/{condition}/{syst}/bjets/{idx}/eta", bjet.Eta(), weight, 48, -2.4, 2.4)
            super().FillHist(f"{channel}/{condition}/{syst}/bjets/{idx}/phi", bjet.Phi(), weight, 64, -3.2, 3.2)
            super().FillHist(f"{channel}/{condition}/{syst}/bjets/{idx}/mass", bjet.M(), weight, 100, 0., 100.)
        super().FillHist(f"{channel}/{condition}/{syst}/bjets/size", bjets.size(), weight, 15, 0., 15.)
        super().FillHist(f"{channel}/{condition}/{syst}/METv/pt", METv.Pt(), weight, 300, 0., 300.)
        super().FillHist(f"{channel}/{condition}/{syst}/METv/phi", METv.Phi(), weight, 64, -3.2, 3.2)
        super().FillHist(f"{channel}/{condition}/{syst}/pair/pt", pair.Pt(), weight, 300, 0., 300.)
        super().FillHist(f"{channel}/{condition}/{syst}/pair/eta", pair.Eta(), weight, 100, -5., 5.)
        super().FillHist(f"{channel}/{condition}/{syst}/pair/phi", pair.Phi(), weight, 64, -3.2, 3.2)
        super().FillHist(f"{channel}/{condition}/{syst}/pair/mass", pair.M(), weight, 300, 0., 300.)
        
        
if __name__ == "__main__":
    m = Valid_POGMuon()
    m.SetTreeName("recoTree/SKFlat")
    m.IsDATA = False
    m.MCSample = "DYJets"
    m.xsec = 6077.22
    m.sumSign = 61192713.0
    m.sumW = 1.545707971038e+12
    m.IsFastSim = False
    m.SetEra("2016preVFP")
    m.Userflags = std.vector[ROOT.TString]()
    m.Userflags.emplace_back("RunSyst")
    if not m.AddFile("/DATA/Users/choij/SKFlat/Run2UltraLegacy_v3/2016preVFP/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220706_092734/0000/SKFlatNtuple_2016preVFP_MC_679.root"): exit(1)
    m.SetOutfilePath("hists.root")
    m.Init()
    m.initializePyAnalyzer()
    m.initializeAnalyzerTools()
    m.SwitchToTempDir()
    m.Loop()
    m.WriteHist()