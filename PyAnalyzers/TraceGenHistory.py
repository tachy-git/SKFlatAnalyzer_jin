from ROOT import gSystem
from ROOT import TriLeptonBase
from ROOT import std
from ROOT.JetTagging import Parameters as jParameters
from ROOT import Lepton, Muon, Electron, Jet
gSystem.Load("/cvmfs/cms.cern.ch/slc7_amd64_gcc900/external/lhapdf/6.2.3/lib/libLHAPDF.so")

class TraceGenHistory(TriLeptonBase):
    def __init__(self):
        super().__init__()

    def executeEvent(self):
        if not super().PassMETFilter(): return None
        ev = super().GetEvent()
        rawMuons = super().GetAllMuons()
        rawElectrons = super().GetAllElectrons()
        rawJets = super().GetAllJets()
        METv = ev.GetMETVector()

        #### Object definition
        vetoMuons = super().SelectMuons(rawMuons, super().MuonIDs[2], 10., 2.4)
        tightMuons = super().SelectMuons(vetoMuons, super().MuonIDs[0], 10., 2.4)
        vetoElectrons = super().SelectElectrons(rawElectrons, super().ElectronIDs[2], 10., 2.5)
        tightElectrons = super().SelectElectrons(vetoElectrons, super().ElectronIDs[0], 10., 2.5)
        jets = super().SelectJets(rawJets, "tight", 15., 2.4)
        # jets = super().JetsVetoLeptonInside(jets, vetoElectrons, vetoMuons, 0.4)
        bjets = std.vector[Jet]()
        for jet in jets:
            score = jet.GetTaggerResult(3)                      # DeepJet
            wp = super().mcCorr.GetJetTaggingCutValue(3, 1)     # DeepJet Medium
            if score > wp: bjets.emplace_back(jet)

        # sort objects
        sorted(vetoMuons, key=lambda x: x.Pt(), reverse=True)
        sorted(tightMuons, key=lambda x: x.Pt(), reverse=True)
        sorted(vetoElectrons, key=lambda x: x.Pt(), reverse=True)
        sorted(tightElectrons, key=lambda x: x.Pt(), reverse=True)
        sorted(jets, key=lambda x: x.Pt(), reverse=True)
        sorted(bjets, key=lambda x: x.Pt(), reverse=True)

         #### event selection
        is3Mu = (len(tightMuons) == 3 and len(vetoMuons) == 3 and \
                len(tightElectrons) == 0 and len(vetoElectrons) == 0)
        is1E2Mu = len(tightMuons) == 2 and len(vetoMuons) == 2 and \
                  len(tightElectrons) == 1 and len(vetoElectrons) == 1
        if not (is3Mu or is1E2Mu): return None

        # prompt matching
        if not super().IsDATA:
            truth = super().GetGens()
            promptMuons = std.vector[Muon]()
            promptElectrons = std.vector[Electron]()
            for mu in tightMuons:
                if super().GetLeptonType(mu, truth) > 0: promptMuons.emplace_back(mu)
            for ele in tightElectrons:
                if super().GetLeptonType(ele, truth) > 0: promptElectrons.emplace_back(ele)

            if len(promptMuons) != len(tightMuons): return None
            if len(promptElectrons) != len(tightElectrons): return None

        channel = ""
        ## 1E2Mu baseline
        ## 1. pass EMuTriggers
        ## 2. Exact 2 tight muons and 1 tight electron, no additional lepton
        ## 3. Exists OS muon pair with mass > 12 GeV
        if is1E2Mu:
            if not ev.PassTrigger(super().EMuTriggers): return None
            leptons = std.vector[Lepton]()
            for mu in tightMuons: leptons.emplace_back(mu)
            for ele in tightElectrons: leptons.emplace_back(ele)
            mu1, mu2, ele = tuple(leptons)
            passLeadMu = mu1.Pt() > 25. and ele.Pt() > 15.
            passLeadEle = mu1.Pt() > 10. and ele.Pt() > 25.
            passSafeCut = passLeadMu or passLeadEle
            if not passSafeCut: return None
            if not mu1.Charge()+mu2.Charge() == 0: return None
            pair = mu1 + mu2
            if not pair.M() > 12.: return None

            # orthogonality of SR and CR done by bjet multiplicity
            if len(bjets) >= 1:
                if len(jets) >= 2: channel = "SR1E2Mu"
                else: return None
            else:
                mZ = 91.2
                isOnZ = abs(pair.M() - mZ) < 10.
                if isOnZ: channel = "ZFake1E2Mu"
                else:
                    if abs((mu1+mu2+ele).M() - mZ) < 10.: channel = "ZGamma1E2Mu"
                    else: return None

        ## 3Mu baseline
        ## 1. pass DblMuTriggers
        ## 2. Exact 3 tight muons, no additional leptons
        ## 3. Exist OS muon pair,
        ## 4. All OS muon pair mass > 12 GeV
        else:
            if not ev.PassTrigger(super().DblMuTriggers): return None
            mu1, mu2, mu3  = tuple(tightMuons)
            if not mu1.Pt() > 20.: return None
            if not mu2.Pt() > 10.: return None
            if not mu3.Pt() > 10.: return None
            if not abs(mu1.Charge()+mu2.Charge()+mu3.Charge()) == 1: return None
            if mu1.Charge() == mu2.Charge():
                pair1 = mu1 + mu3
                pair2 = mu2 + mu3
            elif mu1.Charge() == mu3.Charge():
                pair1 = mu1 + mu2
                pair2 = mu2 + mu3
            else:   # mu2.Charge() == mu3.Charge()
                pair1 = mu1 + mu2
                pair2 = mu1 + mu3
            if not pair1.M() > 12.: return None
            if not pair2.M() > 12.: return None

            # orthogonality of SR and CR done by bjet multiplicity
            if len(bjets) >= 1:
                if len(jets) >= 2: channel = "SR3Mu"
                else: return None
            else:
                mZ = 91.2
                isOnZ = abs(pair1.M() - mZ) < 10. or abs(pair2.M() - mZ) < 10.
                if isOnZ: channel = "ZFake3Mu"
                else:
                    if abs((mu1+mu2+mu3).M() - mZ) < 10.: channel = "ZGamma3Mu"
                    else: return None

        if not ("1E2Mu" in channel or "3Mu" in channel): return None
        
        # find Charged Higgs
        chargedIndices = []
        for gen in truth:
            if abs(gen.PID()) == 37:
                chargedIndices.append(gen.Index())

        # find daughters
        chargedDecays = []
        for gen in truth:
            if gen.MotherIndex() in chargedIndices:
                if abs(gen.PID()) == 37: continue
                chargedDecays.append(gen)
        
        chargedDecayIndices = [gen.Index() for gen in chargedDecays]
        chargedDecayPIDs    = [abs(gen.PID()) for gen in chargedDecays]
        
        if channel == "SR1E2Mu":
            if 0 < chargedDecayPIDs[1] and chargedDecayPIDs[1] < 6:   # Ajj
                q1 = truth.at(chargedDecayIndices[1])
                q2 = truth.at(chargedDecayIndices[2])
                
                # partons should be inside acceptance
                if not q1.Pt() > 15.: return None
                if not q2.Pt() > 15.: return None
                if not abs(q1.Eta()) < 2.4: return None
                if not abs(q2.Eta()) < 2.4: return None
                
                # find nearest jets
                j1 = None
                j2 = None
                dR1, dR2 = 0.3, 0.3
                for jet in jets:
                    if q1.DeltaR(jet) < dR1:
                        dR1 = q1.DeltaR(jet)
                        j1 = jet
                    if q2.DeltaR(jet) < dR2:
                        dR2 = q2.DeltaR(jet)
                        j2 = jet
                if not j1: return None
                if not j2: return None
                if j1 is j2: return None
                WCand = j1 + j2

                # muons
                mu1, mu2 = tuple(tightMuons)
                ACand = mu1 + mu2
                ChargedHiggs = ACand + WCand
                if abs(ChargedHiggs.M() - 130.) > 20.: return None
                super().FillHist(f"{channel}/Ajj/mA", ACand.M(), 1., 200, 0., 200.)
                super().FillHist(f"{channel}/Ajj/mW", WCand.M(), 1., 200, 0., 200.)
                super().FillHist(f"{channel}/Ajj/mHc", ChargedHiggs.M(), 1., 200, 0., 200.)
                super().FillHist(f"{channel}/Ajj/dRj1", j1.DeltaR(q1), 1., 500, 0., 5.)
                super().FillHist(f"{channel}/Ajj/dRj2", j2.DeltaR(q2), 1., 500, 0., 5.)

            else:       # Aenu, Ataunu -> Ae3nu
                ele = tightElectrons.at(0)
                mu1, mu2 = tuple(tightMuons)
                ACand = mu1 + mu2
                WCand = ele + METv
                ChargedHiggs = ACand + WCand
                if abs(ChargedHiggs.M() - 130.) > 40.: return None
                super().FillHist(f"{channel}/Aenu/mA", ACand.M(), 1., 200, 0., 200.)
                super().FillHist(f"{channel}/Aenu/mW", WCand.M(), 1., 200, 0., 200.)
                super().FillHist(f"{channel}/Aenu/mHc", ChargedHiggs.M(), 1., 200, 0., 200.)
        elif channel == "SR3Mu":
            if 0 < chargedDecayPIDs[1] and chargedDecayPIDs[1] < 6:   # Ajj
                q1 = truth.at(chargedDecayIndices[1])
                q2 = truth.at(chargedDecayIndices[2])
                
                # partons should be inside acceptance
                if not q1.Pt() > 15.: return None
                if not q2.Pt() > 15.: return None
                if not abs(q1.Eta()) < 2.4: return None
                if not abs(q2.Eta()) < 2.4: return None

                # find nearest jets
                j1 = None
                j2 = None
                dR1, dR2 = 0.3, 0.3
                for jet in jets:
                    if q1.DeltaR(jet) < dR1:
                        dR1 = q1.DeltaR(jet)
                        j1 = jet
                    if q2.DeltaR(jet) < dR2:
                        dR2 = q2.DeltaR(jet)
                        j2 = jet
                if not j1: return None
                if not j2: return None
                if j1 is j2: return None
                WCand = j1 + j2

                # ACand
                signalColl = []
                for mu in tightMuons:
                    if super().GetLeptonType(mu, truth) == 2:
                        signalColl.append(mu)
                if not len(signalColl) == 2: return None
                mu1, mu2 = tuple(signalColl)
                ACand = mu1 + mu2
                ChargedHiggs = ACand + WCand
                if abs(ChargedHiggs.M() - 130.) > 20.: return None
                super().FillHist(f"{channel}/Ajj/mA", ACand.M(), 1., 200, 0., 200.)
                super().FillHist(f"{channel}/Ajj/mW", WCand.M(), 1., 200, 0., 200.)
                super().FillHist(f"{channel}/Ajj/mHc", ChargedHiggs.M(), 1., 200, 0., 200.)
                super().FillHist(f"{channel}/Ajj/dRj1", j1.DeltaR(q1), 1., 500, 0., 5.)
                super().FillHist(f"{channel}/Ajj/dRj2", j2.DeltaR(q2), 1., 500, 0., 5.)
            else:
                # ACand
                signalColl, promptColl = [], []
                for mu in tightMuons:
                    if super().GetLeptonType(mu, truth) == 2:
                        signalColl.append(mu)
                    elif super().GetLeptonType(mu, truth) > 0:
                        promptColl.append(mu)
                    else:
                        pass
                if not len(signalColl) == 2: return None
                if not len(signalColl) == 1: return None
                mu1, mu2 = tuple(signalColl)
                m3 = promptColl[0]
                ACand = mu1 + mu2
                WCand = mu3 + METv
                ChargedHiggs = ACand + WCand
                if abs(ChargedHiggs.M() - 130.) > 40.: return None
                super().FillHist(f"{channel}/Aenu/mA", ACand.M(), 1., 200, 0., 200.)
                super().FillHist(f"{channel}/Aenu/mW", WCand.M(), 1., 200, 0., 200.)
                super().FillHist(f"{channel}/Aenu/mHc", ChargedHiggs.M(), 1., 200, 0., 200.)
        else:       # Control regions
            pass

        
#if __name__ == "__main__":
#    m = TraceGenHistory()
#    m.SetTreeName("recoTree/SKFlat")
#    m.IsDATA = False
#    m.MCSample = "TTToHcToWAToMuMu_MHc-130_MA-90"
#    m.xsec = 0.015
#    m.sumSign = 599702.0
#    m.sumW = 3270.46
#    m.IsFastSim = False
#    m.SetEra("2017")
#    if not m.AddFile("/home/choij/workspace/DATA/SKFlat/Run2UltraLegacy_v3/2017/TTToHcToWAToMuMu_MHc-130_MA-90_MultiLepFilter_TuneCP5_13TeV-madgraph-pythia8/SKFlat_Run2UltraLegacy_v3/220714_084244/0000/SKFlatNtuple_2017_MC_14.root"): exit(1)
#    if not m.AddFile("/home/choij/workspace/DATA/SKFlat/Run2UltraLegacy_v3/2017/TTToHcToWAToMuMu_MHc-130_MA-90_MultiLepFilter_TuneCP5_13TeV-madgraph-pythia8/SKFlat_Run2UltraLegacy_v3/220714_084244/0000/SKFlatNtuple_2017_MC_5.root"): exit(1)
#    m.SetOutfilePath("hists.root")
#    m.Init()
#    m.initializeAnalyzer()
#    m.initializeAnalyzerTools()
#    m.SwitchToTempDir()
#    m.Loop()
#    m.WriteHist()
