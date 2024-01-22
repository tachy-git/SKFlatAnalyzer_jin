if __name__ == "__main__":
    m = MatrixEstimator()
    m.SetTreeName("recoTree/SKFlat")
    m.IsDATA = True
    m.DataStream = "DoubleMuon"
    m.SetEra("2018")
    m.Userflags = vector[TString]()
    m.Userflags.emplace_back("Skim3Mu")
    m.Userflags.emplace_back("DenseNet")
    if not m.AddFile("/home/choij/workspace/DATA/SKFlat/Run2UltraLegacy_v3/2018/SKFlatNtuple_2018_DATA_4.root"): exit(1)
    m.SetOutfilePath("hists.root")
    m.MaxEvent = m.fChain.GetEntries()
    m.Init()
    m.initializePyAnalyzer()
    m.initializeAnalyzerTools()
    m.SwitchToTempDir()
    m.Loop()
    m.WriteHist()


if __name__ == "__main__":
    m = PromptEstimator()
    m.SetTreeName("recoTree/SKFlat")
    m.IsDATA = False
    m.MCSample = "TTToHcToWAToMuMu_MHc-130_MA-90"
    m.xsec = 0.015
    m.sumSign = 599702.0
    m.sumW = 3270.46
    m.IsFastSim = False
    m.SetEra("2017")
    m.Userflags = vector[TString]()
    m.Userflags.emplace_back("Skim3Mu")
    m.Userflags.emplace_back("GraphNet")
    m.Userflags.emplace_back("WeightVar")
    if not m.AddFile("/home/choij/workspace/DATA/SKFlat/Run2UltraLegacy_v3/2017/TTToHcToWAToMuMu_MHc-130_MA-90_MultiLepFilter_TuneCP5_13TeV-madgraph-pythia8/SKFlat_Run2UltraLegacy_v3/220714_084244/0000/SKFlatNtuple_2017_MC_10.root"): exit(1)
    m.SetOutfilePath("hists.root")
    m.Init()
    m.initializePyAnalyzer()
    m.initializeAnalyzerTools()
    m.SwitchToTempDir()
    m.Loop()
    m.WriteHist()

if __name__ == "__main__":
    m = POGDiMuonValidation()
    m.SetTreeName("recoTree/SKFlat")
    m.IsDATA = False
    m.MCSample = "DYJets"
    m.xsec = 6077.22
    m.sumSign = 61192713.0
    m.sumW = 1.545707971038e+12
    m.IsFastSim = False
    m.SetEra("2016preVFP")
    m.Userflags = std.vector[ROOT.TString]()
    m.Userflags.emplace_back("RunDiMu")
    if not m.AddFile("/DATA/Users/choij/SKFlat/Run2UltraLegacy_v3/2016preVFP/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220706_092734/0000/SKFlatNtuple_2016preVFP_MC_679.root"): exit(1)
    m.SetOutfilePath("hists.root")
    m.Init()
    m.initializePyAnalyzer()
    m.initializeAnalyzerTools()
    m.SwitchToTempDir()
    m.Loop()
    m.WriteHist()
