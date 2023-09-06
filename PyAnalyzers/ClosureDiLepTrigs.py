from ROOT import gSystem
from ROOT import DiLeptonBase
from ROOT import TString
from ROOT.std import vector
from ROOT.JetTagging import Parameters as jParameters
from ROOT import Muon, Electron, Jet
gSystem.Load("/cvmfs/cms.cern.ch/slc7_amd64_gcc900/external/lhapdf/6.2.3/lib/libLHAPDF.so")

class ClosureDiLepTrigs(DiLeptonBase):
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
        if super().RunSyst:     self.run_syst = True
        
        self.weightVariations = ["Central"]
        self.scaleVariations = []
        if self.run_syst:
            self.weightVariations += ["L1PrefireUp", "L1PrefireDown",
                                      "PileupReweightUp", "PileupReweightDown",
                                      "MuonIDSFUp", "MuonIDSFDown",
                                      #"ElectronIDSFUp", "ElectronIDSFDown",
                                      "DblMuTrigSFUp", "DblMuTrigSFDown",
                                      #"EMuTrigSFUp", "EMuTrigSFDown",
                                      #"DYReweightUp", "DYReweightDown",
                                      "HeavyTagUpUnCorr", "HeavyTagDownUnCorr",
                                      "HeavyTagUpCorr", "HeavyTagDownCorr",
                                      "LightTagUpUnCorr", "LightTagDownUnCorr",
                                      "LightTagUpCorr", "LightTagDownCorr"
                                      ]
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
            for syst in self.weightVariations:
                weight = self.getWeight(channel, ev, tightMuons, tightElectrons, jets, truth, syst)
                self.FillObjects(channel, ev, objects, weight, syst)
        
        # Scale variations
        for syst in self.scaleVariations:
            vetoMuons, tightMuons, vetoElectrons, tightElectrons, jets, bjets = self.defineObjects(rawMuons, rawElectrons, rawJets, syst)
            channel = self.selectEvent(ev, vetoMuons, tightMuons, vetoElectrons, tightElectrons, jets, bjets, METv)
            if channel is None: continue
            objects = {"muons": tightMuons,
                       "electrons": tightElectrons,
                       "jets": jets,
                       "bjets": bjets,
                       "METv": METv,
                       }
            weight = self.getWeight(channel, ev, tightMuons, tightElectrons, jets, truth, syst)
            self.FillObjects(channel, ev, objects, weight, syst)
        
    def defineObjects(self, rawMuons, rawElectrons, rawJets, syst="Central"):
        # first copy objects
        allMuons = rawMuons
        allElectrons = rawElectrons
        allJets = rawJets
        
        # check the syst argument
        if not syst in self.systematics:
            print(f"[PromptEstimator::defineObjects] Wrong systematics {syst}")
        if syst == "MuonEnUp":         allMuons = super().ScaleMuons(allMuons, +1)
        if syst == "MuonEnDown":       allMuons = super().ScaleMuons(allMuons, -1)
        if syst == "ElectronResUp":    allElectrons = super().SmearElectrons(allElectrons, +1)
        if syst == "ElectronsResDown": allElectrons = super().SmearElectrons(allElectrons, -1)
        if syst == "ElectronEnUp":     allElectrons = super().ScaleElectrons(allElectrons, +1)
        if syst == "ElectronEnDown":   allElectrons = super().ScaleElectrons(allElectrons, -1)
        if syst == "JetResUp":         allJets = super().SmearJets(allJets, +1)
        if syst == "JetResDown":       allJets = super().SmearJets(allJets, -1)
        if syst == "JetEnUp":          allJets = super().ScaleJets(allJets, +1)
        if syst == "JetEnDown":        allJets = super().ScaleJets(allJets, -1)
        
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
            print("Not implemented yet")
            exit(1)
        else:
            print(f"Wrong channel {self.channel}")
            exit(1)
        return self.channel

    def getWeight(self, channel, event, muons, electrons, jets, truth, syst="Central"):
        weight = 1.
        if not syst in self.systematics:
            print(f"[PromptEstimator::getWeight] Wrong systematic {syst}")
            exit(1)

        if not super().IsDATA:
            weight *= super().MCweight()
            weight *= event.GetTriggerLumi("Full")
            if syst == "L1PrefireUp":     w_prefire = super().GetPrefireWeight(1)
            elif syst == "L1PrefireDown": w_prefire = super().GetPrefireWeight(-1)
            else:                         w_prefire = super().GetPrefireWeight(0)
            
            if syst == "PileupReweightUp":     w_pileup = super().GetPileUpWeight(super().nPileUp, 1)
            elif syst == "PileupReweightDown": w_pileup = super().GetPileUpWeight(super().nPileUp, -1)
            else:                              w_pileup = super().GetPileUpWeight(super().nPileUp, 0)
            
            w_topptweight = 1.
            if "TTLL" in super().MCSample or "TTLJ" in super().MCSample:
                w_topptweight = super().mcCorr.GetTopPtReweight(truth)
            weight *= w_topptweight

            weight *= w_prefire            # print(f"w_prefire: {w_prefire}")
            weight *= w_pileup             # print(f"w_pileup: {w_pileup}")
        
        return weight
    
    def FillObjects(self, channel, evt, objects, weight, syst):
        muons = objects["muons"]
        electrons = objects["electrons"]
        jets = objects["jets"]
        bjets = objects["bjets"]
        METv = objects["METv"]
        
        trigWeight = 1.
        trigWeightUp = 1.
        trigWeightDown = 1.
        if channel == "RunDiMu":
            trigWeight = self.getDblMuTriggerEff(muons, False, 0)
            trigWeightUp = self.getDblMuTriggerEff(muons, False, 1)
            trigWeightDown = self.getDblMuTriggerEff(muons, False, -1)
            trigWeight *= super().getDZEfficiency(channel, isDATA=False)
        elif channel == "RunEMu":
            print("Not implemented yet")
            exit(1)
        else:
            print(f"Wrong channel {channel}")
            exit(1)

        # get weight
        super().FillHist("sumweight", 0., weight, 5, 0., 5.)
        super().FillHist("sumweight", 1., weight*trigWeight, 5, 0., 5.)
        super().FillHist("sumweight", 2., weight*trigWeightUp, 5, 0., 5.)
        super().FillHist("sumweight", 3., weight*trigWeightDown, 5, 0., 5.)
        
        if not evt.PassTrigger(super().DblMuTriggers): return None
        super().FillHist("sumweight", 4., weight, 5, 0., 5.)


if __name__ == "__main__":
    m = ClosureDiLepTrigs()
    m.SetTreeName("recoTree/SKFlat")
    m.IsDATA = False
    m.MCSample = "DYJets"
    m.xsec = 6077.22
    m.sumSign = 61192713.0
    m.sumW = 1.545707971038e+12
    m.IsFastSim = False
    m.SetEra("2016preVFP")
    m.Userflags = vector[TString]()
    m.Userflags.emplace_back("RunDiMu")
    if not m.AddFile("/home/choij/workspace/DATA/SKFlat/Run2UltraLegacy_v3/2016preVFP/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220706_092734/0000/SKFlatNtuple_2016preVFP_MC_964.root"): exit(1)
    m.SetOutfilePath("hists.root")
    m.Init()
    m.initializePyAnalyzer()
    m.initializeAnalyzerTools()
    m.SwitchToTempDir()
    m.Loop()
    m.WriteHist()
