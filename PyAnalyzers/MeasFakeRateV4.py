from ROOT import gSystem
from ROOT import DiLeptonBase
from ROOT import TString
from ROOT.std import vector
from ROOT.JetTagging import Parameters as jParameters
from ROOT import Muon, Electron, Jet
from ROOT import TMath
from itertools import product
gSystem.Load("/cvmfs/cms.cern.ch/slc7_amd64_gcc900/external/lhapdf/6.2.3/lib/libLHAPDF.so")

# from V3 to V4
# Run Analyzer for measure fake rate using two method
# Method1: Fitting for each pt-corr bin, validation in Z-enriched region
# Method2: Scale factor calculation done in Z-enriched region
# Changes: Remove WEnriched region and unnecessary histograms

class MeasFakeRateV4(DiLeptonBase):
    def __init__(self):
        super().__init__() # at this point, DiLeptonBase::initializeAnalyzer has not been called

    def initializePyAnalyzer(self):
        super().initializeAnalyzer()
        
        # link flags
        # minimun weight consider L1prefire and MC weight only
        # full weight additionally considier toppt weight, Lepton Reco SF, NPVweight, b-tagging weight
        self.channel   = None
        self.run_syst  = False
        
        if super().MeasFakeEl8:  self.channel = "MeasFakeEl8"
        if super().MeasFakeEl12: self.channel = "MeasFakeEl12"
        if super().MeasFakeEl23: self.channel = "MeasFakeEl23"
        if super().MeasFakeMu8:  self.channel = "MeasFakeMu8"
        if super().MeasFakeMu17: self.channel = "MeasFakeMu17"
        if super().RunSyst:    self.run_syst = True
        
        # triggers and binning definitions
        # for electron, ptcorr < 20: Ele8, ptcorr < 35: Ele12, ptcorr > 35: Ele23
        # for muon, ptcorr < 30: Mu8, ptcoor > 30: Mu8
        self.triggerList = vector[TString]()
        if self.channel == "MeasFakeEl8":  self.triggerList.emplace_back("HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v")
        if self.channel == "MeasFakeEl12": self.triggerList.emplace_back("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v")
        if self.channel == "MeasFakeEl23": self.triggerList.emplace_back("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v")
        if self.channel == "MeasFakeMu8":  self.triggerList.emplace_back("HLT_Mu8_TrkIsoVVL_v")
        if self.channel == "MeasFakeMu17": self.triggerList.emplace_back("HLT_Mu17_TrkIsoVVL_v") 
        
        self.ptCorr_bins = None
        self.eta_bins = None
        if "MeasFakeEl" in self.channel:
            self.ptCorr_bins = [10., 15., 20., 25., 35., 50., 100.]
            self.eta_bins = [0., 0.8, 1.479, 2.5]
        if "MeasFakeMu" in self.channel:
            self.ptCorr_bins = [10., 15., 20., 30., 50., 100.]
            self.eta_bins = [0., 0.9, 1.6, 2.4]
        
        # systematics
        self.weightVariations = ["Central"]
        self.scaleVariations = []
        self.selectionVariations = []
        if self.run_syst:
            if "MeasFakeEl" in self.channel and not super().IsDATA:
                self.weightVariations += ["PileupReweight",
                                          "L1PrefireUp", "L1PrefireDown",
                                          "ElectronRecoSFUp", "ElectronRecoSFDown"]
            if "MeasFakeMu" in self.channel and not super().IsDATA: 
                self.weightVariations += ["PileupReweight",
                                          "L1PrefireUp", "L1PrefireDown",
                                          "MuonRecoSFUp", "MuonRecoSFDown"]
            self.scaleVariations += ["JetResUp", "JetResDown",
                                     "JetEnUp", "JetEnDown",
                                     "ElectronResUp", "ElectronResDown",
                                     "ElectronEnUp", "ElectronEnDown",
                                     "MuonEnUp", "MuonEnDown"]
            self.selectionVariations += ["RequireHeavyTag", "MotherJetPtUp", "MotherJetPtDown"]
        if super().IsDATA:
            self.systematics = ["Central"] + self.scaleVariations + self.selectionVariations
        else:
            self.systematics = self.weightVariations + self.scaleVariations + self.selectionVariations
        
    def executeEvent(self):
        if not super().PassMETFilter(): return None
        ev = super().GetEvent()
        rawMuons = super().GetAllMuons()
        rawElectrons = super().GetAllElectrons()
        rawJets = super().GetAllJets()
        METv = ev.GetMETVector()
        truth = super().GetGens() if not super().IsDATA else None
        
        # Central scale
        vetoMuons, looseMuons, tightMuons, vetoElectrons, looseElectrons, tightElectrons, jets, bjets = self.defineObjects(rawMuons, rawElectrons, rawJets)
        channel = self.preselectEvent(ev, vetoMuons, looseMuons, vetoElectrons, looseElectrons, jets, bjets, METv)
        if not channel is None:
            objects = {"looseMuons": looseMuons,
                       "looseElectrons": looseElectrons,
                       "tightMuons": tightMuons,
                       "tightElectrons": tightElectrons,
                       "jets": jets,
                       "bjets": bjets,
                       "METv": METv
                       }
            for syst in self.weightVariations:
                weight = self.getWeight(ev, looseMuons, looseElectrons, jets, truth, syst)
                self.fillObjects(channel, objects, weight, syst)
        
        # Scale Variations
        for syst in self.scaleVariations:
            vetoMuons, looseMuons, tightMuons, vetoElectrons, looseElectrons, tightElectrons, jets, bjets = self.defineObjects(rawMuons, rawElectrons, rawJets, syst)
            channel = self.preselectEvent(ev, vetoMuons, looseMuons, vetoElectrons, looseElectrons, jets, bjets, METv)
            if channel is None: continue
            objects = {"looseMuons": looseMuons,
                       "looseElectrons": looseElectrons,
                       "tightMuons": tightMuons,
                       "tightElectrons": tightElectrons,
                       "jets": jets,
                       "bjets": bjets,
                       "METv": METv
                       }
            weight = self.getWeight(ev, looseMuons, looseElectrons, jets, truth, syst)
            self.fillObjects(channel, objects, weight, syst)
            
        # Selection Variations
        for syst in self.selectionVariations:
            vetoMuons, looseMuons, tightMuons, vetoElectrons, looseElectrons, tightElectrons, jets, bjets = self.defineObjects(rawMuons, rawElectrons, rawJets, syst)
            channel = self.preselectEvent(ev, vetoMuons, looseMuons, vetoElectrons, looseElectrons, jets, bjets, METv, syst)
            if channel is None: continue
            objects = {"looseMuons": looseMuons,
                       "looseElectrons": looseElectrons,
                       "tightMuons": tightMuons,
                       "tightElectrons": tightElectrons,
                       "jets": jets,
                       "bjets": bjets,
                       "METv": METv
                       }
            weight = self.getWeight(ev, looseMuons, looseElectrons, jets, truth, syst)
            self.fillObjects(channel, objects, weight, syst)
            
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
        looseMuons = super().SelectMuons(vetoMuons, super().MuonIDs[1], 10., 2.4)
        tightMuons = super().SelectMuons(looseMuons, super().MuonIDs[0], 10., 2.4)
        vetoElectrons = super().SelectElectrons(allElectrons, super().ElectronIDs[2], 10., 2.5)
        looseElectrons = super().SelectElectrons(vetoElectrons, super().ElectronIDs[1], 10., 2.5)
        tightElectrons = super().SelectElectrons(looseElectrons, super().ElectronIDs[0], 10., 2.5)

        if syst == "MotherJetPtUp":     jetPtCut = 50.
        elif syst == "MotherJetPtDown": jetPtCut = 30.
        else:                           jetPtCut = 40.
        jets = super().SelectJets(allJets, "tight", jetPtCut, 2.4)
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
        
        return (vetoMuons, looseMuons, tightMuons, vetoElectrons, looseElectrons, tightElectrons, jets, bjets)
    
    def preselectEvent(self, event, vetoMuons, looseMuons, vetoElectrons, looseElectrons, jets, bjets, METv, syst="Central"):
        NML, NMV = looseMuons.size(), vetoMuons.size()
        NEL, NEV = looseElectrons.size(), vetoElectrons.size()
        if "MeasFakeEl" in self.channel:
            if not event.PassTrigger(self.triggerList): return None
            if not (NEL > 0 and NEV > 0): return None
            if self.channel == "MeasFakeEl12" and not looseElectrons.at(0).Pt() > 15.: return None
            if self.channel == "MeasFakeEl23" and not looseElectrons.at(0).Pt() > 25.: return None
            # single electron
            if (NEL == 1 and NEV == 1 and NML == 0 and NMV == 0):
                if syst == "RequireHeavyTag":
                    if not bjets.size() >= 1: return None
                    if not bjets.at(0).DeltaR(looseElectrons.at(0)) > 1.: return None
                else:
                    if not jets.size() >= 1: return None
                    if not jets.at(0).DeltaR(looseElectrons.at(0)) > 1.: return None
                return "preselSglEl"
            # double electron
            elif (NEL == 2 and NEV == 2 and NML == 0 and NMV == 0):
                ZCand = looseElectrons.at(0)+looseElectrons.at(1)
                if not (60. < ZCand.M() and ZCand.M() < 120.): return None
                return "ZEnriched"
            else:
                return None
        elif "MeasFakeMu" in self.channel:
            if not event.PassTrigger(self.triggerList): return None
            if not (NML > 0 and NMV > 0): return None
            if self.channel == "MeasFakeMu17" and not looseMuons.at(0).Pt() > 20.: return None
            # single muon
            if (NEL == 0 and NEV == 0 and NML == 1 and NMV == 1):
                if syst == "RequireHeavyTag":
                    if not bjets.size() >= 1: return None
                    if not bjets.at(0).DeltaR(looseMuons.at(0)) > 1.: return None
                else:
                    if not jets.size() >= 1: return None
                    if not jets.at(0).DeltaR(looseMuons.at(0)) > 1.: return None
                return "preselSglMu"
            elif (NEL == 0 and NEV == 0 and NML == 2 and NMV == 2):
                ZCand = looseMuons.at(0) + looseMuons.at(1)
                if not (abs(ZCand.M()-91.2) < 15.): return None 
                return "ZEnriched"
            else:
                return None
        else:
            print(f"Wrong channel {self.channel}")
            raise(KeyError)
        
    def getWeight(self, event, muons, electrons, jets, truth, syst="Central"):
        weight = 1.
        if not syst in self.systematics:
            raise KeyError(f"[FakeMeasurement::getWeight] Wrong systematic {syst}")   

        if not super().IsDATA:
            weight *= super().MCweight()
            weight *= event.GetTriggerLumi("Full")
            if syst == "L1PrefireUp":     w_prefire = super().GetPrefireWeight(1)
            elif syst == "L1PrefireDown": w_prefire = super().GetPrefireWeight(-1)
            else:                         w_prefire = super().GetPrefireWeight(0)
            weight *= w_prefire
            
            if not self.run_syst:
                return weight 
                        
            # top pt reweight
            w_topptweight = 1.
            if "TTLL" in super().MCSample or "TTLJ" in super().MCSample:
                w_topptweight = super().mcCorr.GetTopPtReweight(truth)
            weight *= w_topptweight
            
            # Do not apply lepton scale factors
            # it can be absorbed in the prompt normalization scale
            w_muonRecoSF = 1.
            w_eleRecoSF = 1.
            for mu in muons:
                if syst == "MuonRecoSFUp":     w_muonRecoSF *= super().getMuonRecoSF(mu, 1)
                elif syst == "MuonRecoSFDown": w_muonRecoSF *= super().getMuonRecoSF(mu, -1)
                else:                          w_muonRecoSF *= super().getMuonRecoSF(mu, 0)

            for ele in electrons:
                if syst == "ElectronRecoSFUp":   w_eleRecoSF *= super().mcCorr.ElectronReco_SF(ele.scEta(), ele.Pt(), 1)
                elif syst == "ElectronRecoSFUp": w_eleRecoSF *= super().mcCorr.ElectronReco_SF(ele.scEta(), ele.Pt(), -1)
                else:                            w_eleRecoSF *= super().mcCorr.ElectronReco_SF(ele.scEta(), ele.Pt(), 0)
            
            weight *= w_muonRecoSF
            weight *= w_eleRecoSF
            
            if syst == "RequireHeavyTag":
                jtp = jParameters(3, 1, 0, 1)    # DeepJet, medium, incl, mujets
                vjets = vector[Jet]()
                for j in jets: vjets.emplace_back(j)
                w_btag = super().mcCorr.GetBTaggingReweight_1a(vjets, jtp)
                weight *= w_btag
        
        return weight
    
    def fillObjects(self, channel, objects, weight, syst):
        looseMuons = objects["looseMuons"]
        tightMuons = objects["tightMuons"]
        looseElectrons = objects["looseElectrons"]
        tightElectrons = objects["tightElectrons"]
        jets = objects["jets"]
        bjets = objects["bjets"]
        METv = objects["METv"]
        
        if channel == "preselSglMu":
            mu = looseMuons.at(0)
            dPhi = mu.DeltaPhi(METv)
            MT = TMath.Sqrt(2.*mu.Pt()*METv.Pt()*(1.-TMath.Cos(dPhi)))
            MTfix = TMath.Sqrt(2.*35.*METv.Pt()*(1.-TMath.Cos(dPhi)))
            MET = METv.Pt()
            ptCorr = mu.Pt()*(1. + max(0, mu.MiniRelIso()-0.1))
            abseta = abs(mu.Eta())
            prefix = self.findbin(ptCorr, abseta)
            trigSuffix = ""
            if self.channel == "MeasFakeMu8": trigSuffix = "Mu8"
            if self.channel == "MeasFakeMu17": trigSuffix = "Mu17"
            npvweight = 1.
            if self.run_syst:
                if syst == "PileupReweight": npvweight = super().GetPileUpWeight(super().nPileUp, 0)
                else:                        npvweight = self.getNPVReweight(super().nPV, trigSuffix)
                
            # loose
            ID = "loose"
            selection = "Inclusive"
            super().FillHist(f"{selection}/{ID}/{syst}/MT", MT, weight*npvweight, 300, 0., 300.)
            super().FillHist(f"{selection}/{ID}/{syst}/MTfix", MTfix, weight*npvweight, 300, 0., 300.)
            super().FillHist(f"{selection}/{ID}/{syst}/ptCorr", ptCorr, weight*npvweight, 300, 0., 300.)
            super().FillHist(f"{selection}/{ID}/{syst}/abseta", abseta, weight*npvweight, 48, 0., 2.4)
            super().FillHist(f"{selection}/{ID}/{syst}/MET", MET, weight*npvweight, 300, 0., 300.)
            super().FillHist(f"{selection}/{ID}/{syst}/nPV", super().nPV, weight*npvweight, 100, 0., 100.)
            super().FillHist(f"{selection}/{ID}/{syst}/nPileUp", super().nPileUp, weight*npvweight, 100, 0., 100.)
            super().FillHist(f"{selection}/{ID}/{syst}/MT_MET", MT, MET, weight*npvweight, 300, 0., 300., 300, 0., 300.)
            super().FillHist(f"{prefix}/{selection}/{ID}/{syst}/MT", MT, weight*npvweight, 300, 0., 300.)
            super().FillHist(f"{prefix}/{selection}/{ID}/{syst}/MTfix", MTfix, weight*npvweight, 300, 0., 300.)
            super().FillHist(f"{prefix}/{selection}/{ID}/{syst}/ptCorr", ptCorr, weight*npvweight, 300, 0., 300.)
            super().FillHist(f"{prefix}/{selection}/{ID}/{syst}/abseta", abseta, weight*npvweight, 48, 0., 2.4)
            super().FillHist(f"{prefix}/{selection}/{ID}/{syst}/MET", MET, weight*npvweight, 300, 0., 300.)
            super().FillHist(f"{prefix}/{selection}/{ID}/{syst}/nPV", super().nPV, weight*npvweight, 100, 0., 100.)
            super().FillHist(f"{prefix}/{selection}/{ID}/{syst}/nPileUp", super().nPileUp, weight*npvweight, 100, 0., 100.)
            super().FillHist(f"{prefix}/{selection}/{ID}/{syst}/MT_MET", MT, MET, weight*npvweight, 300, 0., 300., 300, 0., 300.)

            if ("QCD" in super().MCSample) or (MT < 25. and MET < 25.):
                selection = "QCDEnriched"
                super().FillHist(f"{selection}/{ID}/{syst}/MT", MT, weight*npvweight, 300, 0., 300.)
                super().FillHist(f"{selection}/{ID}/{syst}/MTfix", MTfix, weight*npvweight, 300, 0., 300.)
                super().FillHist(f"{selection}/{ID}/{syst}/ptCorr", ptCorr, weight*npvweight, 300, 0., 300.)
                super().FillHist(f"{selection}/{ID}/{syst}/abseta", abseta, weight*npvweight, 48, 0., 2.4)
                super().FillHist(f"{selection}/{ID}/{syst}/MET", MET, weight*npvweight, 300, 0., 300.)
                super().FillHist(f"{selection}/{ID}/{syst}/nPV", super().nPV, weight*npvweight, 100, 0., 100.)
                super().FillHist(f"{selection}/{ID}/{syst}/nPileUp", super().nPileUp, weight*npvweight, 100, 0., 100.)
                super().FillHist(f"{selection}/{ID}/{syst}/MT_MET", MT, MET, weight*npvweight, 300, 0., 300., 300, 0., 300.)
                super().FillHist(f"{prefix}/{selection}/{ID}/{syst}/MT", MT, weight*npvweight, 300, 0., 300.)
                super().FillHist(f"{prefix}/{selection}/{ID}/{syst}/MTfix", MTfix, weight*npvweight, 300, 0., 300.)
                super().FillHist(f"{prefix}/{selection}/{ID}/{syst}/ptCorr", ptCorr, weight*npvweight, 300, 0., 300.)
                super().FillHist(f"{prefix}/{selection}/{ID}/{syst}/abseta", abseta, weight*npvweight, 48, 0., 2.4)
                super().FillHist(f"{prefix}/{selection}/{ID}/{syst}/MET", MET, weight*npvweight, 300, 0., 300.)
                super().FillHist(f"{prefix}/{selection}/{ID}/{syst}/nPV", super().nPV, weight*npvweight, 100, 0., 100.)
                super().FillHist(f"{prefix}/{selection}/{ID}/{syst}/nPileUp", super().nPileUp, weight*npvweight, 100, 0., 100.)
                super().FillHist(f"{prefix}/{selection}/{ID}/{syst}/MT_MET", MT, MET, weight*npvweight, 300, 0., 300., 300, 0., 300.)

            # tight
            if tightMuons.size() != 1: return None
            ID = "tight"
            selection = "Inclusive"
            super().FillHist(f"{selection}/{ID}/{syst}/MT", MT, weight*npvweight, 300, 0., 300.)
            super().FillHist(f"{selection}/{ID}/{syst}/MTfix", MTfix, weight*npvweight, 300, 0., 300.)
            super().FillHist(f"{selection}/{ID}/{syst}/ptCorr", ptCorr, weight*npvweight, 300, 0., 300.)
            super().FillHist(f"{selection}/{ID}/{syst}/abseta", abseta, weight*npvweight, 48, 0., 2.4)
            super().FillHist(f"{selection}/{ID}/{syst}/MET", MET, weight*npvweight, 300, 0., 300.)
            super().FillHist(f"{selection}/{ID}/{syst}/nPV", super().nPV, weight*npvweight, 100, 0., 100.)
            super().FillHist(f"{selection}/{ID}/{syst}/nPileUp", super().nPileUp, weight*npvweight, 100, 0., 100.)
            super().FillHist(f"{selection}/{ID}/{syst}/MT_MET", MT, MET, weight*npvweight, 300, 0., 300., 300, 0., 300.)
            super().FillHist(f"{prefix}/{selection}/{ID}/{syst}/MT", MT, weight*npvweight, 300, 0., 300.)
            super().FillHist(f"{prefix}/{selection}/{ID}/{syst}/MTfix", MTfix, weight*npvweight, 300, 0., 300.)
            super().FillHist(f"{prefix}/{selection}/{ID}/{syst}/ptCorr", ptCorr, weight*npvweight, 300, 0., 300.)
            super().FillHist(f"{prefix}/{selection}/{ID}/{syst}/abseta", abseta, weight*npvweight, 48, 0., 2.4)
            super().FillHist(f"{prefix}/{selection}/{ID}/{syst}/MET", MET, weight*npvweight, 300, 0., 300.)
            super().FillHist(f"{prefix}/{selection}/{ID}/{syst}/nPV", super().nPV, weight*npvweight, 100, 0., 100.)
            super().FillHist(f"{prefix}/{selection}/{ID}/{syst}/nPileUp", super().nPileUp, weight*npvweight, 100, 0., 100.)
            super().FillHist(f"{prefix}/{selection}/{ID}/{syst}/MT_MET", MT, MET, weight*npvweight, 300, 0., 300., 300, 0., 300.)

            if ("QCD" in super().MCSample) or (MT < 25. and MET < 25.):
                selection = "QCDEnriched"
                super().FillHist(f"{selection}/{ID}/{syst}/MT", MT, weight*npvweight, 300, 0., 300.)
                super().FillHist(f"{selection}/{ID}/{syst}/MTfix", MTfix, weight*npvweight, 300, 0., 300.)
                super().FillHist(f"{selection}/{ID}/{syst}/ptCorr", ptCorr, weight*npvweight, 300, 0., 300.)
                super().FillHist(f"{selection}/{ID}/{syst}/abseta", abseta, weight*npvweight, 48, 0., 2.4)
                super().FillHist(f"{selection}/{ID}/{syst}/MET", MET, weight*npvweight, 300, 0., 300.)
                super().FillHist(f"{selection}/{ID}/{syst}/nPV", super().nPV, weight*npvweight, 100, 0., 100.)
                super().FillHist(f"{selection}/{ID}/{syst}/nPileUp", super().nPileUp, weight*npvweight, 100, 0., 100.)
                super().FillHist(f"{selection}/{ID}/{syst}/MT_MET", MT, MET, weight*npvweight, 300, 0., 300., 300, 0., 300.)
                super().FillHist(f"{prefix}/{selection}/{ID}/{syst}/MT", MT, weight*npvweight, 300, 0., 300.)
                super().FillHist(f"{prefix}/{selection}/{ID}/{syst}/MTfix", MTfix, weight*npvweight, 300, 0., 300.)
                super().FillHist(f"{prefix}/{selection}/{ID}/{syst}/ptCorr", ptCorr, weight*npvweight, 300, 0., 300.)
                super().FillHist(f"{prefix}/{selection}/{ID}/{syst}/abseta", abseta, weight*npvweight, 48, 0., 2.4)
                super().FillHist(f"{prefix}/{selection}/{ID}/{syst}/MET", MET, weight*npvweight, 300, 0., 300.)
                super().FillHist(f"{prefix}/{selection}/{ID}/{syst}/nPV", super().nPV, weight*npvweight, 100, 0., 100.)
                super().FillHist(f"{prefix}/{selection}/{ID}/{syst}/nPileUp", super().nPileUp, weight*npvweight, 100, 0., 100.)
                super().FillHist(f"{prefix}/{selection}/{ID}/{syst}/MT_MET", MT, MET, weight*npvweight, 300, 0., 300., 300, 0., 300.)
            # end of preselSglMu
        if channel == "ZEnriched" and "MeasFakeMu" in self.channel:
            pair = looseMuons.at(0) + looseMuons.at(1)
            trigSuffix = ""
            if self.channel == "MeasFakeMu8": trigSuffix = "Mu8"
            if self.channel == "MeasFakeMu17": trigSuffix = "Mu17"
            npvweight = 1.
            if self.run_syst:
                if syst == "PileupReweight": npvweight = super().GetPileUpWeight(super().nPileUp, 0)
                else:                        npvweight = self.getNPVReweight(super().nPV, trigSuffix)
            
            # loose
            ID = "loose"
            super().FillHist(f"{channel}/{ID}/{syst}/pair/pt", pair.Pt(), weight*npvweight, 300, 0., 300.)
            super().FillHist(f"{channel}/{ID}/{syst}/pair/eta", pair.Eta(), weight*npvweight, 100, -5., 5.)
            super().FillHist(f"{channel}/{ID}/{syst}/pair/phi", pair.Phi(), weight*npvweight, 64, -3.2, 3.2)
            super().FillHist(f"{channel}/{ID}/{syst}/pair/mass", pair.M(), weight*npvweight, 80, 50., 130.)
            super().FillHist(f"{channel}/{ID}/{syst}/nPV", super().nPV, weight*npvweight, 100, 0., 100.)
            super().FillHist(f"{channel}/{ID}/{syst}/nPileUp", super().nPileUp, weight*npvweight, 100, 0., 100.) 
            
            # tight
            if not tightMuons.size() == 2: return None
            ID = "tight"
            super().FillHist(f"{channel}/{ID}/{syst}/pair/pt", pair.Pt(), weight*npvweight, 300, 0., 300.)
            super().FillHist(f"{channel}/{ID}/{syst}/pair/eta", pair.Eta(), weight*npvweight, 100, -5., 5.)
            super().FillHist(f"{channel}/{ID}/{syst}/pair/phi", pair.Phi(), weight*npvweight, 64, -3.2, 3.2)
            super().FillHist(f"{channel}/{ID}/{syst}/pair/mass", pair.M(), weight*npvweight, 80, 50., 130.)
            super().FillHist(f"{channel}/{ID}/{syst}/nPV", super().nPV, weight*npvweight, 100, 0., 100.)
            super().FillHist(f"{channel}/{ID}/{syst}/nPileUp", super().nPileUp, weight*npvweight, 100, 0., 100.) 
            # end of ZEnriched / MeasFakeMu
        if channel == "preselSglEl":
            ele = looseElectrons.at(0)
            dPhi = ele.DeltaPhi(METv)
            MT = TMath.Sqrt(2.*ele.Pt()*METv.Pt()*(1.-TMath.Cos(dPhi)))
            MTfix = TMath.Sqrt(2.*35.*METv.Pt()*(1.-TMath.Cos(dPhi)))
            MET = METv.Pt()
            ptCorr = ele.Pt()*(1. + max(0, ele.MiniRelIso()-0.1))
            abseta = abs(ele.scEta())
            prefix = self.findbin(ptCorr, abseta)
            trigSuffix = ""
            if self.channel == "MeasFakeEl8":  trigSuffix = "Ele8"
            if self.channel == "MeasFakeEl12": trigSuffix = "Ele12"
            if self.channel == "MeasFakeEl23": trigSuffix = "Ele23"
            npvweight = 1.
            if self.run_syst:
                if syst == "PileupReweight": npvweight = super().GetPileUpWeight(super().nPileUp, 0)
                else:                        npvweight = self.getNPVReweight(super().nPV, trigSuffix)
                
            # loose
            ID = "loose"
            selection = "Inclusive"
            super().FillHist(f"{selection}/{ID}/{syst}/MT", MT, weight*npvweight, 300, 0., 300.)
            super().FillHist(f"{selection}/{ID}/{syst}/MTfix", MTfix, weight*npvweight, 300, 0., 300.)
            super().FillHist(f"{selection}/{ID}/{syst}/ptCorr", ptCorr, weight*npvweight, 300, 0., 300.)
            super().FillHist(f"{selection}/{ID}/{syst}/abseta", abseta, weight*npvweight, 50, 0., 2.5)
            super().FillHist(f"{selection}/{ID}/{syst}/MET", MET, weight*npvweight, 300, 0., 300.)
            super().FillHist(f"{selection}/{ID}/{syst}/nPV", super().nPV, weight*npvweight, 100, 0., 100.)
            super().FillHist(f"{selection}/{ID}/{syst}/nPileUp", super().nPileUp, weight*npvweight, 100, 0., 100.)
            super().FillHist(f"{selection}/{ID}/{syst}/MT_MET", MT, MET, weight*npvweight, 300, 0., 300., 300, 0., 300.)
            super().FillHist(f"{prefix}/{selection}/{ID}/{syst}/MT", MT, weight*npvweight, 300, 0., 300.)
            super().FillHist(f"{prefix}/{selection}/{ID}/{syst}/MTfix", MTfix, weight*npvweight, 300, 0., 300.)
            super().FillHist(f"{prefix}/{selection}/{ID}/{syst}/ptCorr", ptCorr, weight*npvweight, 300, 0., 300.)
            super().FillHist(f"{prefix}/{selection}/{ID}/{syst}/abseta", abseta, weight*npvweight, 50, 0., 2.5)
            super().FillHist(f"{prefix}/{selection}/{ID}/{syst}/MET", MET, weight*npvweight, 300, 0., 300.)
            super().FillHist(f"{prefix}/{selection}/{ID}/{syst}/nPV", super().nPV, weight*npvweight, 100, 0., 100.)
            super().FillHist(f"{prefix}/{selection}/{ID}/{syst}/nPileUp", super().nPileUp, weight*npvweight, 100, 0., 100.)
            super().FillHist(f"{prefix}/{selection}/{ID}/{syst}/MT_MET", MT, MET, weight*npvweight, 300, 0., 300., 300, 0., 300.)

            if ("QCD" in super().MCSample) or (MT < 25. and MET < 25.):
                selection = "QCDEnriched"
                super().FillHist(f"{selection}/{ID}/{syst}/MT", MT, weight*npvweight, 300, 0., 300.)
                super().FillHist(f"{selection}/{ID}/{syst}/MTfix", MTfix, weight*npvweight, 300, 0., 300.)
                super().FillHist(f"{selection}/{ID}/{syst}/ptCorr", ptCorr, weight*npvweight, 300, 0., 300.)
                super().FillHist(f"{selection}/{ID}/{syst}/abseta", abseta, weight*npvweight, 50, 0., 2.5)
                super().FillHist(f"{selection}/{ID}/{syst}/MET", MET, weight*npvweight, 300, 0., 300.)
                super().FillHist(f"{selection}/{ID}/{syst}/nPV", super().nPV, weight*npvweight, 100, 0., 100.)
                super().FillHist(f"{selection}/{ID}/{syst}/nPileUp", super().nPileUp, weight*npvweight, 100, 0., 100.)
                super().FillHist(f"{selection}/{ID}/{syst}/MT_MET", MT, MET, weight*npvweight, 300, 0., 300., 300, 0., 300.)
                super().FillHist(f"{prefix}/{selection}/{ID}/{syst}/MT", MT, weight*npvweight, 300, 0., 300.)
                super().FillHist(f"{prefix}/{selection}/{ID}/{syst}/MTfix", MTfix, weight*npvweight, 300, 0., 300.)
                super().FillHist(f"{prefix}/{selection}/{ID}/{syst}/ptCorr", ptCorr, weight*npvweight, 300, 0., 300.)
                super().FillHist(f"{prefix}/{selection}/{ID}/{syst}/abseta", abseta, weight*npvweight, 50, 0., 2.5)
                super().FillHist(f"{prefix}/{selection}/{ID}/{syst}/MET", MET, weight*npvweight, 300, 0., 300.)
                super().FillHist(f"{prefix}/{selection}/{ID}/{syst}/nPV", super().nPV, weight*npvweight, 100, 0., 100.)
                super().FillHist(f"{prefix}/{selection}/{ID}/{syst}/nPileUp", super().nPileUp, weight*npvweight, 100, 0., 100.)
                super().FillHist(f"{prefix}/{selection}/{ID}/{syst}/MT_MET", MT, MET, weight*npvweight, 300, 0., 300., 300, 0., 300.)

            # tight
            if tightElectrons.size() != 1: return None
            ID = "tight"
            selection = "Inclusive"
            super().FillHist(f"{selection}/{ID}/{syst}/MT", MT, weight*npvweight, 300, 0., 300.)
            super().FillHist(f"{selection}/{ID}/{syst}/MTfix", MTfix, weight*npvweight, 300, 0., 300.)
            super().FillHist(f"{selection}/{ID}/{syst}/ptCorr", ptCorr, weight*npvweight, 300, 0., 300.)
            super().FillHist(f"{selection}/{ID}/{syst}/abseta", abseta, weight*npvweight, 50, 0., 2.5)
            super().FillHist(f"{selection}/{ID}/{syst}/MET", MET, weight*npvweight, 300, 0., 300.)
            super().FillHist(f"{selection}/{ID}/{syst}/nPV", super().nPV, weight*npvweight, 100, 0., 100.)
            super().FillHist(f"{selection}/{ID}/{syst}/nPileUp", super().nPileUp, weight*npvweight, 100, 0., 100.)
            super().FillHist(f"{selection}/{ID}/{syst}/MT_MET", MT, MET, weight*npvweight, 300, 0., 300., 300, 0., 300.)
            super().FillHist(f"{prefix}/{selection}/{ID}/{syst}/MT", MT, weight*npvweight, 300, 0., 300.)
            super().FillHist(f"{prefix}/{selection}/{ID}/{syst}/MTfix", MTfix, weight*npvweight, 300, 0., 300.)
            super().FillHist(f"{prefix}/{selection}/{ID}/{syst}/ptCorr", ptCorr, weight*npvweight, 300, 0., 300.)
            super().FillHist(f"{prefix}/{selection}/{ID}/{syst}/abseta", abseta, weight*npvweight, 50, 0., 2.5)
            super().FillHist(f"{prefix}/{selection}/{ID}/{syst}/MET", MET, weight*npvweight, 300, 0., 300.)
            super().FillHist(f"{prefix}/{selection}/{ID}/{syst}/nPV", super().nPV, weight*npvweight, 100, 0., 100.)
            super().FillHist(f"{prefix}/{selection}/{ID}/{syst}/nPileUp", super().nPileUp, weight*npvweight, 100, 0., 100.)
            super().FillHist(f"{prefix}/{selection}/{ID}/{syst}/MT_MET", MT, MET, weight*npvweight, 300, 0., 300., 300, 0., 300.)

            if ("QCD" in super().MCSample) or (MT < 25. and MET < 25.):
                selection = "QCDEnriched"
                super().FillHist(f"{selection}/{ID}/{syst}/MT", MT, weight*npvweight, 300, 0., 300.)
                super().FillHist(f"{selection}/{ID}/{syst}/MTfix", MTfix, weight*npvweight, 300, 0., 300.)
                super().FillHist(f"{selection}/{ID}/{syst}/ptCorr", ptCorr, weight*npvweight, 300, 0., 300.)
                super().FillHist(f"{selection}/{ID}/{syst}/abseta", abseta, weight*npvweight, 50, 0., 2.5)
                super().FillHist(f"{selection}/{ID}/{syst}/MET", MET, weight*npvweight, 300, 0., 300.)
                super().FillHist(f"{selection}/{ID}/{syst}/nPV", super().nPV, weight*npvweight, 100, 0., 100.)
                super().FillHist(f"{selection}/{ID}/{syst}/nPileUp", super().nPileUp, weight*npvweight, 100, 0., 100.)
                super().FillHist(f"{selection}/{ID}/{syst}/MT_MET", MT, MET, weight*npvweight, 300, 0., 300., 300, 0., 300.)
                super().FillHist(f"{prefix}/{selection}/{ID}/{syst}/MT", MT, weight*npvweight, 300, 0., 300.)
                super().FillHist(f"{prefix}/{selection}/{ID}/{syst}/MTfix", MTfix, weight*npvweight, 300, 0., 300.)
                super().FillHist(f"{prefix}/{selection}/{ID}/{syst}/ptCorr", ptCorr, weight*npvweight, 300, 0., 300.)
                super().FillHist(f"{prefix}/{selection}/{ID}/{syst}/abseta", abseta, weight*npvweight, 50, 0., 2.5)
                super().FillHist(f"{prefix}/{selection}/{ID}/{syst}/MET", MET, weight*npvweight, 300, 0., 300.)
                super().FillHist(f"{prefix}/{selection}/{ID}/{syst}/nPV", super().nPV, weight*npvweight, 100, 0., 100.)
                super().FillHist(f"{prefix}/{selection}/{ID}/{syst}/nPileUp", super().nPileUp, weight*npvweight, 100, 0., 100.)
                super().FillHist(f"{prefix}/{selection}/{ID}/{syst}/MT_MET", MT, MET, weight*npvweight, 300, 0., 300., 300, 0., 300.)
            # end of preselSglEl
        if channel == "ZEnriched" and "MeasFakeEl" in self.channel:
            pair = looseElectrons.at(0) + looseElectrons.at(1)
            trigSuffix = ""
            if self.channel == "MeasFakeEl8":  trigSuffix = "Ele8"
            if self.channel == "MeasFakeEl12": trigSuffix = "Ele12"
            if self.channel == "MeasFakeEl23": trigSuffix = "Ele23"
            npvweight = 1.
            if self.run_syst:
                if syst == "PileupReweight": npvweight = super().GetPileUpWeight(super().nPileUp, 0)
                else:                        npvweight = self.getNPVReweight(super().nPV, trigSuffix)
            
            # loose
            ID = "loose"
            super().FillHist(f"{channel}/{ID}/{syst}/pair/pt", pair.Pt(), weight*npvweight, 300, 0., 300.)
            super().FillHist(f"{channel}/{ID}/{syst}/pair/eta", pair.Eta(), weight*npvweight, 100, -5., 5.)
            super().FillHist(f"{channel}/{ID}/{syst}/pair/phi", pair.Phi(), weight*npvweight, 64, -3.2, 3.2)
            super().FillHist(f"{channel}/{ID}/{syst}/pair/mass", pair.M(), weight*npvweight, 80, 50., 130.)
            super().FillHist(f"{channel}/{ID}/{syst}/nPV", super().nPV, weight*npvweight, 100, 0., 100.)
            super().FillHist(f"{channel}/{ID}/{syst}/nPileUp", super().nPileUp, weight*npvweight, 100, 0., 100.) 
            
            # tight
            if not tightElectrons.size() == 2: return None
            ID = "tight"
            super().FillHist(f"{channel}/{ID}/{syst}/pair/pt", pair.Pt(), weight*npvweight, 300, 0., 300.)
            super().FillHist(f"{channel}/{ID}/{syst}/pair/eta", pair.Eta(), weight*npvweight, 100, -5., 5.)
            super().FillHist(f"{channel}/{ID}/{syst}/pair/phi", pair.Phi(), weight*npvweight, 64, -3.2, 3.2)
            super().FillHist(f"{channel}/{ID}/{syst}/pair/mass", pair.M(), weight*npvweight, 80, 50., 130.)
            super().FillHist(f"{channel}/{ID}/{syst}/nPV", super().nPV, weight*npvweight, 100, 0., 100.)
            super().FillHist(f"{channel}/{ID}/{syst}/nPileUp", super().nPileUp, weight*npvweight, 100, 0., 100.) 
            # end of ZEnriched / MeasFakeEl
            
    def findbin(self, ptCorr, abseta):
        if ptCorr > 100.: 
            #print(f"ptCorr = {ptCorr}")
            ptCorr = 99.
        
        prefix = ""
        # find bin index for ptcorr
        for i in range(len(self.ptCorr_bins)-1):
            if self.ptCorr_bins[i] < ptCorr < self.ptCorr_bins[i+1]:
                prefix += f"ptCorr_{int(self.ptCorr_bins[i])}to{int(self.ptCorr_bins[i+1])}"
                break
        # find bin index for abseta
        for i in range(len(self.eta_bins)-1):
            if self.eta_bins[i] < abseta < self.eta_bins[i+1]:
                prefix += f"_abseta_{str(self.eta_bins[i]).replace('.', 'p')}to{str(self.eta_bins[i+1]).replace('.', 'p')}"
        return prefix
