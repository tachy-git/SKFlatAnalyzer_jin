#include "SkimTree_Dilepton.h"

void SkimTree_Dilepton::initializeAnalyzer(){

  outfile->cd();
  cout << "[SkimTree_Dilepton::initializeAnalyzer()] gDirectory = " << gDirectory->GetName() << endl;
  newtree = fChain->CloneTree(0);

  double_triggers.clear();
  single_electron_triggers.clear();
  single_muon_triggers.clear();
  if(DataYear==2016){
    double_triggers = {
      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v",
      "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v",
      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",
      "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",
      "HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v",
      "HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",              // H
      "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",
      "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v",       // B-G
      "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",      // B-G
      "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",    // H
      "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"    // H
    };
    single_muon_triggers = {
      "HLT_IsoMu24_v",
      "HLT_IsoTkMu24_v",
    };
    single_electron_triggers = {
      "HLT_Ele27_WPTight_Gsf_v",
    };
  }else if(DataYear==2017){
    double_triggers = {
      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",                  // B
      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v",          // C-F
      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v",
      "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v",
      "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",      // C-F
      "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",
      "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"
    };
    single_muon_triggers = {
      "HLT_IsoMu24_v",
      "HLT_IsoMu27_v",
    };
    single_electron_triggers = {
      "HLT_Ele27_WPTight_Gsf_v",
      "HLT_Ele32_WPTight_Gsf_v",
      "HLT_Ele32_WPTight_Gsf_L1DoubleEG_v",
    };
  }else if(DataYear==2018){
    double_triggers = {
      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v",
      "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v",
      "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",
      "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",
      "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"
    };
    single_muon_triggers = {
      "HLT_IsoMu24_v",
    };
    single_electron_triggers = {
      "HLT_Ele27_WPTight_Gsf_v",
      "HLT_Ele28_WPTight_Gsf_v",
      "HLT_Ele32_WPTight_Gsf_v",
    };
  }else{
    cout<<"[SkimTree_Dilepton::initializeAnalyzer] DataYear is wrong : " << DataYear << endl;
  }

  cout << "[SkimTree_Dilepton::initializeAnalyzer] triggers to skim = " << endl;
  for(unsigned int i=0; i<double_triggers.size(); i++){
    cout << "[SkimTree_Dilepton::initializeAnalyzer]   " << double_triggers.at(i) << endl;
  }
  for(unsigned int i=0; i<single_muon_triggers.size(); i++){
    cout << "[SkimTree_Dilepton::initializeAnalyzer]   " << single_muon_triggers.at(i) << endl;
  }
  for(unsigned int i=0; i<single_electron_triggers.size(); i++){
    cout << "[SkimTree_Dilepton::initializeAnalyzer]   " << single_electron_triggers.at(i) << endl;
  }

}

void SkimTree_Dilepton::executeEvent(){

  Event ev;
  ev.SetTrigger(*HLT_TriggerName);

  if( ev.PassTrigger(double_triggers) ){
    newtree->Fill();
    return;
  }

  vector<Muon> muons=GetAllMuons();
  std::sort(muons.begin(),muons.end(),PtComparing);
  vector<Electron> electrons=GetAllElectrons();
  std::sort(electrons.begin(),electrons.end(),PtComparing);
  int nmuon=muons.size();
  int nelectron=electrons.size();

  if(ev.PassTrigger(single_muon_triggers)){
    for(int i=0;i<nmuon;i++){
      for(int j=i+1;j<nmuon;j++){
	if(muons.at(i).Pt()<20) continue;
	if(muons.at(j).Pt()<7) continue;
	if((muons.at(i)+muons.at(j)).M()<48) continue;
	if(newtree->Fill()<0) exit(EIO);
	return;
      }
    }
    for(int i=0;i<nmuon;i++){
      for(int j=0;j<nelectron;j++){
	if(muons.at(i).Pt()<20) continue;
	if(electrons.at(j).Pt()<7) continue;
	if((muons.at(i)+electrons.at(j)).M()<48) continue;
	if(newtree->Fill()<0) exit(EIO);
	return;
      }
    }
  }

  if(ev.PassTrigger(single_electron_triggers)){
    for(int i=0;i<nelectron;i++){
      for(int j=i+1;j<nelectron;j++){
	if(electrons.at(i).Pt()<25) continue;
	if(electrons.at(j).Pt()<7) continue;
	if((electrons.at(i)+electrons.at(j)).M()<48) continue;
	if(newtree->Fill()<0) exit(EIO);
	return;
      }
    }
    for(int i=0;i<nelectron;i++){
      for(int j=0;j<nmuon;j++){
	if(electrons.at(i).Pt()<25) continue;
	if(muons.at(j).Pt()<7) continue;
	if((electrons.at(i)+muons.at(j)).M()<48) continue;
	if(newtree->Fill()<0) exit(EIO);
	return;
      }
    }
  }
}

void SkimTree_Dilepton::executeEventFromParameter(AnalyzerParameter param){

}

SkimTree_Dilepton::SkimTree_Dilepton(){
  newtree=NULL;
}

SkimTree_Dilepton::~SkimTree_Dilepton(){

}

void SkimTree_Dilepton::WriteHist(){

  outfile->mkdir("recoTree");
  outfile->cd("recoTree");
  newtree->Write();
  outfile->cd();

}
