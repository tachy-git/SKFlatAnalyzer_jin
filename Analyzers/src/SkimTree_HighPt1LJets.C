#include "SkimTree_HighPt1LJets.h"

void SkimTree_HighPt1LJets::initializeAnalyzer(){

  outfile->cd();
  cout << "[SkimTree_HighPt1LJets::initializeAnalyzer()] gDirectory = " << gDirectory->GetName() << endl;
  newtree = fChain->CloneTree(0);

  /*
  if(!IsDATA){
    newtree->SetBranchStatus("gen_*",0);
    newtree->SetBranchStatus("LHE_*",0);
    newtree->SetBranchStatus("gen_weight",1); // for MCweight()
  }*/

  triggers.clear();

  if(DataYear==2016){
    triggers = {
      "HLT_IsoMu24_v",
      "HLT_Mu50_v",
      "HLT_TkMu50_v",
      "HLT_Ele27_WPTight_Gsf_v",
      "HLT_Photon175_v",
      "HLT_Ele115_CaloIdVT_GsfTrkIdT_v",
    };
  }
  else if(DataYear==2017){
    triggers = {
      "HLT_IsoMu27_v",
      "HLT_Mu50_v",
      "HLT_OldMu100_v",
      "HLT_TkMu100_v",
      "HLT_Ele35_WPTight_Gsf_v",
      "HLT_Photon200_v",
      "HLT_Ele115_CaloIdVT_GsfTrkIdT_v",
    };
  }
  else if(DataYear==2018){
    triggers = {
      "HLT_IsoMu24_v",
      "HLT_Mu50_v",
      "HLT_OldMu100_v",
      "HLT_TkMu100_v",
      "HLT_Ele32_WPTight_Gsf_v",
      "HLT_Photon200_v",
      "HLT_Ele115_CaloIdVT_GsfTrkIdT_v",
    };
  }
  else{
    cout << "[SkimTree_HighPt1LJets::initializeAnalyzer] DataYear is wrong : " << DataYear << endl;
  }

  cout << "[SkimTree_HighPt1LJets::initializeAnalyzer] triggers to skim = " << endl;
  for(unsigned int i=0; i<triggers.size(); i++){
    cout << "[SkimTree_HighPt1LJets::initializeAnalyzer]   " << triggers.at(i) << endl;
  }

  LeptonPtCut = 18.;
  AK4JetPtCut = 20.;
  AK8JetPtCut = 170.;

}

void SkimTree_HighPt1LJets::executeEvent(){

    //==== Skim 1 ) Triggers
    Event ev;
    ev.SetTrigger(*HLT_TriggerName);
    if(! (ev.PassTrigger(triggers)) ) return;

    //==== Skim 2) Leptons
    vector<Muon> muons = GetMuons("NOCUT", LeptonPtCut, 2.4);
    vector<Electron> electrons = GetElectrons("NOCUT", LeptonPtCut, 2.5);
    int n_muons = muons.size();
    int n_electrons = electrons.size();
    if (n_muons + n_electrons == 0) return;

    //==== Skim 3) Jets
    vector<FatJet> fatjets = puppiCorr->Correct( GetFatJets("tight", AK8JetPtCut, 2.4) ); //==== corret SDMass
    vector<Jet> jets = GetJets("tightLepVeto", AK4JetPtCut, 2.4);
    int n_fatjets = fatjets.size();
    int n_jets = jets.size();
    if ((n_fatjets == 0) && (n_jets <= 1)) return;

    //==== Skim 4) Only fill in the survived events
    newtree->Fill();

}

void SkimTree_HighPt1LJets::executeEventFromParameter(AnalyzerParameter param){

}

SkimTree_HighPt1LJets::SkimTree_HighPt1LJets(){

  newtree = NULL;

}

SkimTree_HighPt1LJets::~SkimTree_HighPt1LJets(){

}

void SkimTree_HighPt1LJets::WriteHist(){

  outfile->mkdir("recoTree");
  outfile->cd("recoTree");
  newtree->Write();
  outfile->cd();

}


