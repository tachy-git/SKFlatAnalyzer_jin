R__LOAD_LIBRARY(/home/choij/workspace/SKFlatAnalyzer/external/lhapdf/lib/libLHAPDF.so)
void runMeasFakeRateV2(){
    MeasFakeRateV2 m;
    m.SetTreeName("recoTree/SKFlat");
    m.MCSample = "DYJets";
    m.IsDATA = false;
    m.xsec = 6077.22;
    m.sumSign = 196329377;
    m.sumW = 132089877.0;
    m.IsFastSim = false;
    m.SetEra("2017");
    m.Userflags.emplace_back("MeasFakeEl12");
    if(!m.AddFile("/DATA/SKFlat/Run2UltraLegacy_v3/2017/MC/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/220614_125548/0000/SKFlatNtuple_2017_MC_1.root")) exit(EIO);
    m.SetOutfilePath("hists.root");
    m.Init();
    m.initializeAnalyzer();
    m.initializeAnalyzerTools();
    m.SwitchToTempDir();
    m.Loop();
    m.WriteHist();
}
