from ROOT import gSystem
from ROOT import TLorentzVector
from ROOT import TutorialBase
gSystem.Load("/cvmfs/cms.cern.ch/slc7_amd64_gcc900/external/lhapdf/6.2.3/lib/libLHAPDF.so")

class TutorialRun(TutorialBase):
    def __init__(self):
        super().__init__()

    # Override executeEvent function
    def executeEvent(self):
        allMuons = super().GetAllMuons()
        allJets = super().GetAllJets()

        for muonID in self.MuonIDs:
            super().FillHist(f"{muonID}/Cutflow", 0., 1., 10, 0., 10.)
            
            # MET Filter
            if not super().PassMETFilter(): return None
            ev = super().GetEvent()
            METv = ev.GetMETVector()

            # Trigger
            if not ev.PassTrigger(super().IsoMuTriggerName): return None
            super().FillHist(f"{muonID}/Cutflow", 2., 1., 10, 0., 10.)
            muons = list(super().SelectMuons(allMuons, muonID, 20., 2.4))
            jets = list(super().SelectJets(allJets, "tight", 30., 2.4))
            # Sort in pt order
            muons.sort(key=lambda x: x.Pt(), reverse=True)
            jets.sort(key=lambda x: x.Pt(), reverse=True)

            # b-tagging
            bjets = []
            for jet in jets:
                btagScore = jet.GetTaggerResult(3) # which is DeepJet
                wp = super().mcCorr.GetJetTaggingCutValue(3, 1) # DeepJet Medium
                if (btagScore > wp): bjets.append(jet)
            
            # Event selection
            if not (len(muons) == 2): return None
            if not (muons[0].Charge() + muons[1].Charge() == 0): return None
            if not (muons[0].Pt() > super().TriggerSafePtCut): return None
            ZCand = muons[0]+muons[1]
            MZ = 91.2
            if not (abs(ZCand.M() - MZ) < 15.): return None

            weight = 1.
            if not super().IsDATA:
                genWeight = super().MCweight()
                lumi = ev.GetTriggerLumi("Full")
                weight *= (genWeight * lumi)

                # Muon ID efficiency
                #for mu in muons:
                #    print(super().mcCorr.MuonID_SF("NUM_TopHN_DEN_TrackerMuons", mu.Eta(), mu.Pt()))



            super().FillHist(f"{muonID}/ZCand/mass", ZCand.M(), weight, 40, 70., 110.);
            super().FillHist(f"{muonID}/ZCand/pt", ZCand.Pt(), weight, 300, 0., 300.);
            super().FillHist(f"{muonID}/ZCand/eta", ZCand.Eta(), weight, 100, -5., 5.);
            super().FillHist(f"{muonID}/ZCand/phi", ZCand.Phi(), weight, 64, -3.2, 3.2);
            super().FillHist(f"{muonID}/MET", METv.Pt(), weight, 300, 0., 300.);
            super().FillHist(f"{muonID}/bjets/size", len(bjets), weight, 20, 0., 20.);

if __name__ == "__main__":
    m = TutorialRun()
    m.SetTreeName("recoTree/SKFlat")
    m.IsDATA = False
    m.MCSample = "TTToHcToWAToMuMu_MHc-130_MA-90"
    m.xsec = 0.015
    m.sumSign = 599702.0
    m.sumW = 3270.46
    m.IsFastSim = False
    m.SetEra("2017")
    if not m.AddFile("/home/choij/workspace/DATA/SKFlat/Run2UltraLegacy_v3/2017/TTToHcToWAToMuMu_MHc-130_MA-90_MultiLepFilter_TuneCP5_13TeV-madgraph-pythia8/SKFlat_Run2UltraLegacy_v3/220714_084244/0000/SKFlatNtuple_2017_MC_14.root"): exit(1)
    if not m.AddFile("/home/choij/workspace/DATA/SKFlat/Run2UltraLegacy_v3/2017/TTToHcToWAToMuMu_MHc-130_MA-90_MultiLepFilter_TuneCP5_13TeV-madgraph-pythia8/SKFlat_Run2UltraLegacy_v3/220714_084244/0000/SKFlatNtuple_2017_MC_5.root"): exit(1)
    m.SetOutfilePath("hists.root")
    m.Init()
    m.initializeAnalyzer()
    m.initializeAnalyzerTools()
    m.SwitchToTempDir()
    m.Loop()
    m.WriteHist()
