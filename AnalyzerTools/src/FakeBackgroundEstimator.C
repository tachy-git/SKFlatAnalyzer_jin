#include "FakeBackgroundEstimator.h"

FakeBackgroundEstimator::FakeBackgroundEstimator() : 
IgnoreNoHist(false),
HasLooseLepton(false)
{

  histDir = TDirectoryHelper::GetTempDirectory("FakeBackgroundEstimator");

}

void FakeBackgroundEstimator::ReadHistograms(){

  TString datapath = getenv("DATA_DIR");
  datapath = datapath+"/"+GetEra()+"/FakeRate/";

  TDirectory* origDir = gDirectory;

  string elline;
  ifstream in(datapath+"/histmap_Electron.txt");
  while(getline(in,elline)){
    std::istringstream is( elline );
    TString a,b,c,d,e;
    is >> a; // <ID>
    is >> b; // <rootfilename>
    TFile *file = new TFile(datapath+"/"+b);
    TList *histlist = file->GetListOfKeys();
    for(int i=0;i<histlist->Capacity();i++){
      TString this_frname = histlist->At(i)->GetName();
      histDir->cd();
      map_hist_Electron[a+"_"+this_frname] = (TH2D *)file->Get(this_frname)->Clone();
      file->Close();
      delete file;
      origDir->cd();
      //cout << "[FakeBackgroundEstimator::FakeBackgroundEstimator] map_hist_Electron : " << a+"_"+this_frname << endl;
    }
  }

  string elline2;
  ifstream in2(datapath+"/histmap_Muon.txt");
  while(getline(in2,elline2)){
    std::istringstream is( elline2 );
    TString a,b,c,d,e;
    is >> a; // <ID>
    is >> b; // <rootfilename>
    TFile *file = new TFile(datapath+"/"+b);
    TList *histlist = file->GetListOfKeys();
    for(int i=0;i<histlist->Capacity();i++){
      TString this_frname = histlist->At(i)->GetName();
      histDir->cd();
      map_hist_Muon[a+"_"+this_frname] = (TH2D *)file->Get(this_frname)->Clone();
      file->Close();
      delete file;
      origDir->cd();
      //cout << "[FakeBackgroundEstimator::FakeBackgroundEstimator] map_hist_Muon : " << a+"_"+this_frname << endl;
    }
  }

}

FakeBackgroundEstimator::~FakeBackgroundEstimator(){

}

double FakeBackgroundEstimator::GetElectronFakeRate(TString ID, TString key, double eta, double pt, int sys){

  //cout << "[FakeBackgroundEstimator::GetElectronFakeRate] ID = " << ID << ", key = " << key << endl;
  //cout << "[FakeBackgroundEstimator::GetElectronFakeRate] eta = " << eta << ", pt = " << pt << endl;

  double value = 1.;
  double error = 0.;

  eta = fabs(eta);

  if(pt>=200) pt = 199;
  if(eta>=2.5) eta = 2.49;
  //==== HOTFIX FIXME
  if(eta>=1.479 && pt>=150) pt = 149.;

  std::map< TString, TH2D* >::const_iterator mapit;
  mapit = map_hist_Electron.find(ID+"_"+key);

  if(mapit==map_hist_Electron.end()){
    if(IgnoreNoHist) return 1.;
    else{
      cout << "[FakeBackgroundEstimator::GetElectronFakeRate] No"<< ID+"_"+key <<endl;
      exit(ENODATA);
    }
  }

  int this_bin = (mapit->second)->FindBin(pt,eta);
  value = (mapit->second)->GetBinContent(this_bin);
  error = (mapit->second)->GetBinError(this_bin);

  //cout << "[FakeBackgroundEstimator::FakeBackgroundEstimator] value = " << value << endl;

  return value+double(sys)*error;

}

double FakeBackgroundEstimator::GetMuonFakeRate(TString ID, TString key, double eta, double pt, int sys){

  //cout << "[FakeBackgroundEstimator::GetMuonFakeRate] ID = " << ID << ", key = " << key << endl;
  //cout << "[FakeBackgroundEstimator::GetMuonFakeRate] eta = " << eta << ", pt = " << pt << endl;

  double value = 1.;
  double error = 0.;

  eta = fabs(eta);

  if(pt>=200) pt = 199;
  if(eta>=2.5) eta = 2.49;
  //==== HOTFIX FIXME
  if(eta<0.8 && pt>=150) pt = 149.;

  std::map< TString, TH2D* >::const_iterator mapit;
  mapit = map_hist_Muon.find(ID+"_"+key);

  if(mapit==map_hist_Muon.end()){
    if(IgnoreNoHist) return 1.;
    else{
      cout << "[FakeBackgroundEstimator::GetMuonFakeRate] No"<< ID+"_"+key <<endl;
      exit(ENODATA);
    }
  }

  int this_bin = (mapit->second)->FindBin(pt,eta);
  value = (mapit->second)->GetBinContent(this_bin);
  error = (mapit->second)->GetBinError(this_bin);

  //cout << "[FakeBackgroundEstimator::FakeBackgroundEstimator] value = " << value << endl;

  return value+double(sys)*error;

}

double FakeBackgroundEstimator::GetWeight(vector<Lepton *> lepptrs, AnalyzerParameter param, int sys){

  double this_weight = -1.;
  vector<double> FRs;
  for(unsigned int i=0; i<lepptrs.size(); i++){

    double this_fr = -999.;

    if(lepptrs.at(i)->LeptonFlavour()==Lepton::ELECTRON){

      Electron *el = (Electron *)( lepptrs.at(i) );
      if( el->PassID(param.Electron_Tight_ID) ) continue;

      double this_pt = el->Pt();
      if(param.Electron_UsePtCone) this_pt = el->PtCone();
      this_fr = GetElectronFakeRate(param.Electron_FR_ID, param.Electron_FR_Key, fabs(el->scEta()), this_pt, sys);

      this_weight *= -1.*this_fr/(1.-this_fr);

      FRs.push_back(this_fr);
    }
    else{

      Muon *mu = (Muon *)( lepptrs.at(i) );
      if( mu->PassID(param.Muon_Tight_ID) ) continue;

      double this_pt = mu->Pt();
      if(param.Muon_UsePtCone) this_pt = mu->PtCone();
      this_fr = GetMuonFakeRate(param.Muon_FR_ID, param.Muon_FR_Key, fabs(mu->Eta()), this_pt, sys);

      this_weight *= -1.*this_fr/(1.-this_fr);

      FRs.push_back(this_fr);
    }


  }

/*
  for(unsigned int i=0; i<FRs.size(); i++){
    cout << FRs.at(i) << endl;
  }
*/

  if(FRs.size()==0){
    HasLooseLepton = false;
    return 0;
  }
  else{
    HasLooseLepton = true;
    return this_weight;
  }

}




