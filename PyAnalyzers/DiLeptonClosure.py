from ROOT import gSystem
from ROOT import DiLeptonBase
from ROOT import TString
from ROOT.std import vector
from ROOT.JetTagging import Parameters as jParameters
from ROOT import Muon, Electron, Jet
gSystem.Load("/cvmfs/cms.cern.ch/slc7_amd64_gcc900/external/lhapdf/6.2.3/lib/libLHAPDF.so")

class DiLeptonClosure(DiLeptonBase):
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

        ## DiMu selection
        #if self.channel == "RunDiMu":
        #    if not event.PassTrigger(super().DblMuTriggers): return None
        #    mu1, mu2 = tuple(tightMuons)
        #    if not mu1.Pt() > 20.: return None
        #    if not mu2.Pt() > 10.: return None
        #    if not mu1.Charge()+mu2.Charge() == 0: return None
        #    pair = mu1 + mu2
        #    if abs(pair.M() - 91.2) < 15.:
        #        if not bjets.size() == 0: return None
        #        if not METv.Pt() < 30.:   return None
        #        return "DYDiMu"
        #    else:
        #        if not jets.size() >= 2:  return None
        #        if not bjets.size() >= 1: return None
        #        if not pair.M() > 12.:    return None
        #        if not METv.Pt() > 40.:   return None
        #        return "TTDiMu"
        ## EMu selection
        #elif self.channel == "RunEMu":
        #    print("Not implemented yet")
        #    exit(1)
        #else:
        #    print(f"Wrong channel {self.channel}")
    
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
            
            w_zptweight = 1.
            w_topptweight = 1.
            #if "DYJets" in super().MCSample:
            #    if syst == "DYReweightUp":     w_zptweight = super().mcCorr.GetOfficialDYReweight(truth, 1)
            #    elif syst == "DYReweightDown": w_zptweight = super().mcCorr.GetOfficialDYReweight(truth, -1)
            #    else:                          w_zptweight = super().mcCorr.GetOfficialDYReweight(truth, 0)
            if "TTLL" in super().MCSample or "TTLJ" in super().MCSample:
                w_topptweight = super().mcCorr.GetTopPtReweight(truth)
            weight *= (w_zptweight * w_topptweight)

            #w_muonRecoSF = 1.
            #w_muonIDSF = 1.
            #w_dblMuTrigSF = 1.
            #for mu in muons:
            #    w_muonRecoSF *= super().getMuonRecoSF(mu, 0)
            #    if syst == "MuonIDSFUp":     w_muonIDSF *= super().getMuonIDSF(mu, 1)
            #    elif syst == "MuonIDSFDown": w_muonIDSF *= super().getMuonIDSF(mu, -1)
            #    else:                        w_muonIDSF *= super().getMuonIDSF(mu, 0)

            #if "DiMu" in channel:
            #    # trigger efficiency
            #    if syst == "DblMuTrigSFUp":     w_dblMuTrigSF = self.getDblMuTriggerSF(muons, 1)
            #    elif syst == "DblMuTrigSFDown": w_dblMuTrigSF = self.getDblMuTriggerSF(muons, -1)
            #    else:                           w_dblMuTrigSF = self.getDblMuTriggerSF(muons, 0)
                # DZ efficiency
            #    w_dblMuTrigSF *= super().getDZEfficiency(channel, isDATA=True)/super().getDZEfficiency(channel, isDATA=False)

            weight *= w_prefire            # print(f"w_prefire: {w_prefire}")
            weight *= w_pileup             # print(f"w_pileup: {w_pileup}")
            #weight *= w_muonRecoSF         # print(f"w_muonRecoSF: {w_muonRecoSF}")
            #weight *= w_muonIDSF           # print(f"muonID: {w_muonIDSF}")
            #weight *= w_dblMuTrigSF        # print(f"muontrig: {w_dblMuTrigSF}")

            # b-tagging
            jtp = jParameters(3, 1, 0, 1)    # DeepJet, Medium, incl, mujets
            vjets = vector[Jet]()
            for j in jets: vjets.emplace_back(j)
            if syst == "HeavyTagUpUnCorr":     w_btag = super().mcCorr.GetBTaggingReweight_1a(vjets, jtp, "SystUpHTag")
            elif syst == "HeavyTagDownUnCorr": w_btag = super().mcCorr.GetBTaggingReweight_1a(vjets, jtp, "SystDownHTag")
            elif syst == "HeavyTagUpCorr":     w_btag = super().mcCorr.GetBTaggingReweight_1a(vjets, jtp, "SystUpHTagCorr")
            elif syst == "HeavyTagDownCorr":   w_btag = super().mcCorr.GetBTaggingReweight_1a(vjets, jtp, "SystDownHTagCorr")
            elif syst == "LightTagUpUnCorr":   w_btag = super().mcCorr.GetBTaggingReweight_1a(vjets, jtp, "SystUpLTag")
            elif syst == "LightTagDownUnCorr": w_btag = super().mcCorr.GetBTaggingReweight_1a(vjets, jtp, "SystDownLTag")
            elif syst == "LightTagUpCorr":     w_btag = super().mcCorr.GetBTaggingReweight_1a(vjets, jtp, "SystUpLTagCorr")
            elif syst == "LightTagDownCorr":   w_btag = super().mcCorr.GetBTaggingReweight_1a(vjets, jtp, "SystDownLTagCorr")
            else:                              w_btag = super().mcCorr.GetBTaggingReweight_1a(vjets, jtp)
            weight *= w_btag
        
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
        super().FillHist("sumweight", 2., wieght*trigWeightUp, 5, 0., 5.)
        super().FillHist("sumweight", 3., weight*trigWeightDown, 5, 0., 5.)
        
        if not evt.PassTrigger(super().DblMuTriggers): return None
        super().FillHist("sumweight", 4., weight, 5, 0., 5.)

