#include "DiLeptonBase.h"

DiLeptonBase::DiLeptonBase() {}

void DiLeptonBase::initializeAnalyzer() {
    // flags
    RunDiMu = HasFlag("RunDiMu");
    RunEMu = HasFlag("RunEMu");
    RunSyst = HasFlag("RunSyst");
    //MeasFakeMu = HasFlag("MeasFakeMu");
    //MeasFakeEl = HasFlag("MeasFakeEl");
    MeasFakeMu8 = HasFlag("MeasFakeMu8");
    MeasFakeMu17 = HasFlag("MeasFakeMu17");
    MeasFakeEl8 = HasFlag("MeasFakeEl8");
    MeasFakeEl12 = HasFlag("MeasFakeEl12");
    MeasFakeEl23 = HasFlag("MeasFakeEl23");

    // triggers & ID settings
    if (DataEra == "2016preVFP") {
        DblMuTriggers = {
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v",
            "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v",
        };
        EMuTriggers = {
            "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v",
            "HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v",
        };
        MuonIDs = {"HcToWATight", "HcToWALoose", "HcToWAVeto"};
        ElectronIDs = {"HcToWATight16a", "HcToWALoose16a", "HcToWAVeto16a"};
    }
    else if (DataEra == "2016postVFP") {
        DblMuTriggers = {
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",
            "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",
            "HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v"
        };
        EMuTriggers = {
            "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",
            "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v"
        };
        MuonIDs = {"HcToWATight", "HcToWALoose", "HcToWAVeto"};
        ElectronIDs = {"HcToWATight16b", "HcToWALoose16b", "HcToWAVeto16b"};
    }
    else if (DataEra == "2017") {
        DblMuTriggers = {
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v"
        };
        EMuTriggers = {
            "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",
            "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"
        };
        MuonIDs = {"HcToWATight", "HcToWALoose", "HcToWAVeto"};
        ElectronIDs = {"HcToWATight17", "HcToWALoose17", "HcToWAVeto17"};
    }
    else if (DataEra == "2018") {
        DblMuTriggers = {
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v"
        };
        EMuTriggers = {
            "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",
            "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v"
        };
        MuonIDs = {"HcToWATight", "HcToWALoose", "HcToWAVeto"};
        ElectronIDs = {"HcToWATight18", "HcToWALoose18", "HcToWAVeto18"};
    }
    else {
        cerr << "[DiLeptonBase::initializeAnalyzer] Wrong era " << DataEra << endl;
        exit(EXIT_FAILURE);
    }

    TString datapath = getenv("DATA_DIR");
    // muon ID
    TString muonIDpath = datapath + "/" + GetEra() + "/ID/Muon";
    TFile *fMuonID = new TFile(muonIDpath+"/efficiency_TopHNT.root");
    hMuonIDSF = (TH2D*)fMuonID->Get("sf");  hMuonIDSF->SetDirectory(0);
    fMuonID->Close();

    // muon trigger legs
    TFile *fMu17Leg1 = new TFile(muonIDpath+"/efficiency_Mu17Leg1.root");
    hMu17Leg1_Data = (TH2D*)fMu17Leg1->Get("data"); hMu17Leg1_Data->SetDirectory(0);
    hMu17Leg1_MC = (TH2D*)fMu17Leg1->Get("sim");    hMu17Leg1_MC->SetDirectory(0);
    fMu17Leg1->Close();

    TFile* fMu8Leg2 = new TFile(muonIDpath+"/efficiency_Mu8Leg2.root");
    hMu8Leg2_Data = (TH2D*)fMu8Leg2->Get("data"); hMu8Leg2_Data->SetDirectory(0);
    hMu8Leg2_MC = (TH2D*)fMu8Leg2->Get("sim");    hMu8Leg2_MC->SetDirectory(0);
    fMu8Leg2->Close();

    // ele ID
    TString eleIDPath = datapath + "/" + GetEra() + "/ID/Electron";
    TFile *fEleID = new TFile(eleIDPath+"/efficiency_TopHNT.root");
    hElIDSF = (TH2D*)fEleID->Get("sf"); hElIDSF->SetDirectory(0);
    fEleID->Close();

    TFile *fEl23Leg1 = new TFile(eleIDPath+"/efficiency_El23Leg1.root");
    hEl23Leg1_Data = (TH2D*)fEl23Leg1->Get("data"); hEl23Leg1_Data->SetDirectory(0);
    hEl23Leg1_MC = (TH2D*)fEl23Leg1->Get("sim");    hEl23Leg1_MC->SetDirectory(0);
    fEl23Leg1->Close();

    TFile *fEl12Leg2 = new TFile(eleIDPath+"/efficiency_El12Leg2.root");
    hEl12Leg2_Data = (TH2D*)fEl12Leg2->Get("data"); hEl12Leg2_Data->SetDirectory(0);
    hEl12Leg2_MC = (TH2D*)fEl12Leg2->Get("sim");    hEl12Leg2_MC->SetDirectory(0);
    fEl12Leg2->Close();

    // NPV distributions
    TString PUPath = datapath + "/" + GetEra() + "/PileUp";
    if (MeasFakeMu8 || MeasFakeMu17) {
        TFile *fNPV_Data = new TFile(PUPath+"/NPVMuon_DATA.root");
        TFile *fNPV_MC   = new TFile(PUPath+"/NPVMuon_MC.root");
        hNPVMu8_Data = (TH1D*)fNPV_Data->Get("Inclusive_Mu8/loose/Central/nPV");  hNPVMu8_Data->SetDirectory(0);
        hNPVMu17_Data = (TH1D*)fNPV_Data->Get("Inclusive_Mu8/loose/Central/nPV"); hNPVMu17_Data->SetDirectory(0);
        hNPVMu8_MC = (TH1D*)fNPV_MC->Get("Inclusive_Mu8/loose/Central/nPV");  hNPVMu8_MC->SetDirectory(0);
        hNPVMu17_MC = (TH1D*)fNPV_MC->Get("Inclusive_Mu8/loose/Central/nPV"); hNPVMu17_MC->SetDirectory(0);
        fNPV_Data->Close();
        fNPV_MC->Close();
        // scale
        hNPVMu8_Data->Scale(1./hNPVMu8_Data->Integral());
        hNPVMu17_Data->Scale(1./hNPVMu17_Data->Integral());
        hNPVMu8_MC->Scale(1./hNPVMu8_MC->Integral());
        hNPVMu17_MC->Scale(1./hNPVMu17_MC->Integral());
    }
    if (MeasFakeEl8 || MeasFakeEl12 || MeasFakeEl23) {
        TFile *fNPV_Data = new TFile(PUPath+"/NPVElectron_DATA.root");
        TFile *fNPV_MC   = new TFile(PUPath+"/NPVElectron_MC.root");
        hNPVEl8_Data = (TH1D*)fNPV_Data->Get("Inclusive_Ele8/loose/Central/nPV");   hNPVEl8_Data->SetDirectory(0);
        hNPVEl12_Data = (TH1D*)fNPV_Data->Get("Inclusive_Ele12/loose/Central/nPV"); hNPVEl12_Data->SetDirectory(0);
        hNPVEl23_Data = (TH1D*)fNPV_Data->Get("Inclusive_Ele23/loose/Central/nPV"); hNPVEl23_Data->SetDirectory(0);
        hNPVEl8_MC = (TH1D*)fNPV_MC->Get("Inclusive_Ele8/loose/Central/nPV");   hNPVEl8_MC->SetDirectory(0);
        hNPVEl12_MC = (TH1D*)fNPV_MC->Get("Inclusive_Ele12/loose/Central/nPV"); hNPVEl12_MC->SetDirectory(0);
        hNPVEl23_MC = (TH1D*)fNPV_MC->Get("Inclusive_Ele23/loose/Central/nPV"); hNPVEl23_MC->SetDirectory(0);
        fNPV_Data->Close();
        fNPV_MC->Close();
        // scale
        hNPVEl8_Data->Scale(1./hNPVEl8_Data->Integral());
        hNPVEl12_Data->Scale(1./hNPVEl12_Data->Integral());
        hNPVEl23_Data->Scale(1./hNPVEl23_Data->Integral());
        hNPVEl8_MC->Scale(1./hNPVEl8_MC->Integral());
        hNPVEl12_MC->Scale(1./hNPVEl12_MC->Integral());
        hNPVEl23_MC->Scale(1./hNPVEl23_MC->Integral());
    }
    
    // Jet tagger
    vector<JetTagging::Parameters> jtps;
    jtps.emplace_back(JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Medium, JetTagging::incl, JetTagging::mujets));
    mcCorr->SetJetTaggingParameters(jtps);
}



double DiLeptonBase::getMuonRecoSF(const Muon &mu, int sys) {
    const double abseta = fabs(mu.Eta());

    if (DataEra == "2016preVFP") {
        if (abseta < 0.9) {
            const double value = 0.9998229551300333;
            const double stat = 0.0001538802103231026;
            const double syst = 0.0003540897399334497;
            return value + sys*sqrt(pow(stat, 2)+pow(syst, 2));
        }
        else if (abseta < 1.2) {
            const double value = 1.0001593416915515;
            const double stat = 0.00019861903120026457;
            const double syst = 0.00031024592139106355;
            return value + sys*sqrt(pow(stat, 2)+pow(syst, 2)); 
        }
        else if (abseta < 2.1) {
            const double value = 0.9998936144006075;
            const double stat = 0.00012188589514012365;
            const double syst = 0.00021277119878493345;
            return value + sys*sqrt(pow(stat, 2)+pow(syst, 2)); 
        }
        else if (abseta < 2.4) {
            const double value = 0.9990268820042745;
            const double stat = 0.00027638902644996395;
            const double syst = 0.0019462359914510508;
            return value + sys*sqrt(pow(stat, 2)+pow(syst, 2)); 
        }
        else {
            cerr << "[DiLeptonBase::getMuonRecoSF] wrong muon eta value " << abseta << endl;
            exit(EXIT_FAILURE); 
        }
    }
    else if (DataEra == "2016postVFP") {
        if (abseta < 0.9) {
            const double value = 1.0000406419782646;
            const double stat = 0.00010260291858070426;
            const double syst = 0.0014366927652431664;
            return value + sys*sqrt(pow(stat, 2)+pow(syst, 2));
        }
        else if (abseta < 1.2) {
            const double value = 0.9997959311146515;
            const double stat = 0.00019912837537507789;
            const double syst = 0.0010917857343065423;
            return value + sys*sqrt(pow(stat, 2)+pow(syst, 2));
        }
        else if (abseta < 2.1) {
            const double value = 0.9994928400570587;
            const double stat = 0.00012513847429973846;
            const double syst = 0.0014814654032937547;
            return value + sys*sqrt(pow(stat, 2)+pow(syst, 2));
        }
        else if (abseta < 2.4) {
            const double value = 0.9990728619505579;
            const double stat = 0.0002754474704705526;
            const double syst = 0.0017364778744567663;
            return value + sys*sqrt(pow(stat, 2)+pow(syst, 2));
        }
        else {
            cerr << "[DiLeptonBase::getMuonRecoSF] wrong muon eta value " << abseta << endl;
            exit(EXIT_FAILURE);
        }
    }
    else if (DataEra == "2017") {
        if (abseta < 0.9) {
            const double value = 0.9996742562806361;
            const double stat = 7.650191371261136e-05;
            const double syst = 0.0006514874387277825;
            return value + sys*sqrt(pow(stat, 2)+pow(syst, 2));
        }
        else if (abseta < 1.2) {
            const double value = 0.9997813602035737;
            const double stat = 0.00014496238686164667;
            const double syst = 0.0004372795928526685;
            return value + sys*sqrt(pow(stat, 2)+pow(syst, 2));
        }
        else if (abseta < 2.1) {
            const double value = 0.9994674742459532;
            const double stat = 7.739510750489317e-05;
            const double syst = 0.0010650515080936618;
            return value + sys*sqrt(pow(stat, 2)+pow(syst, 2));
        }
        else if (abseta < 2.4) {
            const double value = 0.9993566412630517;
            const double stat = 0.00022835790507860388;
            const double syst = 0.0011810962222705494;
            return value + sys*sqrt(pow(stat, 2)+pow(syst, 2));
        }
        else {
            cerr << "[DiLeptonBase::getMuonRecoSF] wrong muon eta value " << abseta << endl;
            exit(EXIT_FAILURE);
        }
    }
    else if (DataEra == "2018") {
        if (abseta < 0.9) {
            const double value = 0.9998088006315689;
            const double stat = 6.498845788247257e-05;
            const double syst = 0.0003823987368622994;
            return value + sys*sqrt(pow(stat, 2)+pow(syst, 2));
        }
        else if (abseta < 1.2) {
            const double value = 0.999754701980269;
            const double stat = 0.00011054079511271507;
            const double syst = 0.0005221124230931915;
            return value + sys*sqrt(pow(stat, 2)+pow(syst, 2));
        }
        else if (abseta < 2.1) {
            const double value = 0.9995842791862117;
            const double stat = 7.574443994874554e-05;
            const double syst = 0.0008314416275765346;
            return value + sys*sqrt(pow(stat, 2)+pow(syst, 2));
        }
        else if (abseta < 2.4) {
            const double value = 0.9990341741614288;
            const double stat = 0.00019911479235592246;
            const double syst = 0.0017237408292350668;
            return value + sys*sqrt(pow(stat, 2)+pow(syst, 2));
        }
        else {
            cerr << "[DiLeptonBase::getMuonRecoSF] wrong muon eta value " << abseta << endl;
            exit(EXIT_FAILURE);
        }
    }
    else {
        cerr << "[DiLeptonBase::getMuonRecoSF] not implemented era" << DataEra << endl;
        exit(EXIT_FAILURE);
    }
}

void DiLeptonBase::executeEvent() {
    return;
}

double DiLeptonBase::getMuonIDSF(const Muon &mu, int sys) {
    double pt = mu.Pt();
    double eta = fabs(mu.Eta());
    if (pt < 10.) pt = 10;
    if (pt >= 200) pt = 199.;
    if (eta > 2.4) eta = 2.39;
    int thisBin = hMuonIDSF->FindBin(eta, pt);
    double value = hMuonIDSF->GetBinContent(thisBin);
    double error = hMuonIDSF->GetBinError(thisBin);
    
    return value + float(sys)*error;
}

double DiLeptonBase::getEleIDSF(const Electron &ele, int sys) {
    double pt = ele.Pt();
    double eta = ele.scEta();
    if (pt < 10.) pt = 10.;
    if (pt >= 500) pt = 499.;
    if (eta < -2.5) eta = -2.499;
    if (eta >= 2.5) eta = 2.499;

    int thisBin = hElIDSF->FindBin(eta, pt);
    double value = hElIDSF->GetBinContent(thisBin);
    double error = hElIDSF->GetBinContent(thisBin);

    return value + float(sys)*error;
}

double DiLeptonBase::getTriggerEff(const Muon &mu, TString histkey, bool isDATA, int sys) {
    TH2D *h = nullptr;
    double pt = mu.Pt();
    double eta = fabs(mu.Eta());
    if (histkey == "Mu17Leg1" && isDATA) {
        h = hMu17Leg1_Data;
        if (pt < 20.) pt = 20.;
        if (pt >= 200.) pt = 199.;
        if (eta > 2.4) eta = 2.39;
    }
    else if (histkey == "Mu17Leg1" && (!isDATA)) {
        h = hMu17Leg1_MC;
        if (pt < 20.) pt = 20.;
        if (pt >= 200.) pt = 199.;
        if (eta > 2.4) eta = 2.39;
    }
    else if (histkey == "Mu8Leg2" && isDATA) {
        h = hMu8Leg2_Data;
        if (pt < 10.) pt = 10.;
        if (pt >= 200.) pt = 199.;
        if (eta > 2.4) eta = 2.39;
    }
    else if (histkey == "Mu8Leg2" && (!isDATA)) {
        h = hMu8Leg2_MC;
        if (pt < 10.) pt = 10.;
        if (pt >= 200.) pt = 199.;
        if (eta > 2.4) eta = 2.39;
    }
    else {
        cerr << "[DiLeptonBase::getTriggerEff] Wrong combination of histkey and isDataEff" << endl;
        cerr << "[DiLeptonBase::getTriggerEff] histkey = " << histkey << endl;
        cerr << "[DiLeptonBase::getTriggerEff] isDATA = " << isDATA << endl;
    }

    int thisBin = h->FindBin(eta, pt);
    double value = h->GetBinContent(thisBin);
    double error = h->GetBinError(thisBin);

    return value + float(sys)*error;
}

double DiLeptonBase::getTriggerEff(const Electron &ele, TString histkey, bool isDATA, int sys) {
    TH2D *h = nullptr;
    double pt = ele.Pt();
    double eta = ele.scEta();
    if (eta < -2.5) eta = -2.499;
    if (eta >= 2.5) eta = 2.499;
    if (histkey == "El23Leg1" && isDATA) {
        h = hEl23Leg1_Data;
        if (pt < 25.) pt = 25.;
        if (pt >= 200.) pt = 199.;
    }
    else if (histkey == "El23Leg1" && (!isDATA)) {
        h = hEl23Leg1_MC;
        if (pt < 25.) pt = 25.;
        if (pt > 200.) pt = 199.;
    }
    else if (histkey == "El12Leg2" && isDATA) {
        h = hEl12Leg2_Data;
        if (pt < 15.) pt = 15.;
        if (pt > 200.) pt = 199.;
    }
    else if (histkey == "El12Leg2" && (!isDATA)) {
        h = hEl12Leg2_MC;
        if (pt < 15.) pt = 15.;
        if (pt > 200.) pt = 199.;
    }
    else {
        cerr << "[DiLeptonBase::getTriggerEff] Wrong combination of histkey and isDataEff" << endl;
        cerr << "[DiLeptonBase::getTriggerEff] histkey = " << histkey << endl;
        cerr << "[DiLeptonBase::getTriggerEff] isDATA = " << isDATA << endl;
    }
    int thisBin = h->FindBin(eta, pt);
    double value = h->GetBinContent(thisBin);
    double error = h->GetBinError(thisBin);

    return value + float(sys)*error;
}

double DiLeptonBase::getDblMuTriggerEff(vector<Muon> &muons, bool isDATA, int sys) {
    // check no. of muons
    if (! (muons.size() == 2)) {
        cerr << "[DiLeptonBase::getDblMuTriggerEff] Wrong no. of muons " << muons.size() << endl;
        exit(EXIT_FAILURE);
    }

    Muon &mu1 = muons.at(0);
    Muon &mu2 = muons.at(1);

    return getTriggerEff(mu1, "Mu17Leg1", isDATA, sys) * getTriggerEff(mu2, "Mu8Leg2", isDATA, sys);
}

double DiLeptonBase::getEMuTriggerEff(vector<Electron> &electrons, vector<Muon> &muons, bool isDATA, int sys) {
    // check no. of leptons
    if (! (electrons.size() == 1 && muons.size() == 1)) {
        cerr << "[DiLeptonBase::getEMuTriggerEff] Wrong no. of leptons " << electrons.size() << " " << muons.size() << endl;
        exit(EXIT_FAILURE);
    }
    Electron &ele = electrons.at(0);
    Muon &mu = muons.at(0);

    double eff_el, eff_mu;
    eff_el = mu.Pt() > 20. ? getTriggerEff(ele, "El12Leg2", isDATA, sys) : getTriggerEff(ele, "El23Leg1", isDATA, sys);
    eff_mu = ele.Pt() > 25. ? getTriggerEff(mu, "Mu8Leg2", isDATA, sys) : getTriggerEff(mu, "Mu17Leg1", isDATA, sys);
    return eff_el * eff_mu;
}

double DiLeptonBase::getDZEfficiency(TString SFKey, bool isDATA) {
    double eff = 0.;
    if (SFKey.Contains("DiMu")) {
        if (DataEra == "2016postVFP") eff = isDATA ? 0.9798 : 0.9969;
        else if (DataEra == "2017")   eff = 0.9958;
        else                          eff = 1.;
    }
    else if (SFKey.Contains("DiElIso")) {
        if (DataEra=="2016preVFP")       eff = 0.986;
        else if (DataEra=="2016postVFP") eff = 0.980;
        else                             eff = 1.;
    }
    else if (SFKey.Contains("EMu")){
        if(DataEra=="2016postVFP") eff = isDATA ? 0.9648:0.9882;
        //else if(DataEra=="2017"  ) Eff = 0.9951; //for now included in muleg
        else                       eff = 1.;
    }
    else {
        eff = 1.;
    }

    return eff;
}

double DiLeptonBase::getDblMuTriggerSF(vector<Muon> &muons, int sys) {
   double effData = getDblMuTriggerEff(muons, true, sys); 
   double effMC   = getDblMuTriggerEff(muons, false, sys);
   if (effMC == 0 || effData == 0)
       return 1.;

   return effData / effMC;
}

double DiLeptonBase::getEMuTriggerSF(vector<Electron> &electrons, vector<Muon> &muons, int sys) {
    double effData = getEMuTriggerEff(electrons, muons, true, sys);
    double effMC = getEMuTriggerEff(electrons, muons, false, sys);
    if (effMC == 0 || effData == 0) return 1.;
    return effData / effMC;
}

double DiLeptonBase::getNPVReweight(unsigned int NPV, TString &path) {
    if (NPV < 1) NPV = 1;
    if (NPV > 70) NPV = 69;

    TH1D *h_data, *h_mc;
    if (path == "Ele8") {
        h_data = hNPVEl8_Data; h_mc = hNPVEl8_MC;
    }
    else if (path == "Ele12") {
        h_data = hNPVEl12_Data; h_mc = hNPVEl12_MC;
    }
    else if (path == "Ele23") {
        h_data = hNPVEl23_Data; h_mc = hNPVEl23_MC;
    }
    else if (path == "Mu8") {
        h_data = hNPVMu8_Data; h_mc = hNPVMu8_MC;
    }
    else if (path == "Mu17") {
        h_data = hNPVMu17_Data; h_mc = hNPVMu17_MC;
    }
    else {
        cerr << "[DiLeptonBase::getNPVReweight] Wrong path " << path << endl;
        exit(EXIT_FAILURE);
    }
    const unsigned int thisBin = h_data->FindBin(NPV);
    return h_data->GetBinContent(thisBin) / h_mc->GetBinContent(thisBin);
}


DiLeptonBase::~DiLeptonBase() {
    delete hMuonIDSF;
    delete hMu17Leg1_Data;
    delete hMu17Leg1_MC;
    delete hMu8Leg2_Data;
    delete hMu8Leg2_MC;
}
