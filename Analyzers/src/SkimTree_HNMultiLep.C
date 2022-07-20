#include "SkimTree_HNMultiLep.h"

void SkimTree_HNMultiLep::initializeAnalyzer(){

  outfile->cd();
  cout << "[SkimTree_HNMultiLep::initializeAnalyzer()] gDirectory = " << gDirectory->GetName() << endl;
  newtree = fChain->CloneTree(0);
  

  triggers_dimu.clear();

  triggers.clear();
  if(DataYear==2016){
    triggers_dimu = {
      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v",
      "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v",
      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",
      "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",
      "HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v",
      "HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",

    };
    triggers = {
      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v",
      "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v",
      "HLT_Mu50_v",
      "HLT_TkMu50_v",
      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",
      "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",
      "HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v",
      "HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",              // H
      "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",
      "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v",     // B-G
      "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",    // B-G
      "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",  // H
      "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"  // H
    };
  }
  else if(DataYear==2017){
    triggers_dimu = {
      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",
      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v",
      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v"

    };

    triggers = {
      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",                  // B
      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v",          // C-F
      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v",
      "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v",
      "HLT_Mu50_v",
      "HLT_TkMu50_v",
      "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",      // C-F
      "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",
      "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"

    };
  }
  else if(DataYear==2018){
    triggers_dimu = {
      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v",

    };
    triggers = {
      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v",
      "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v",
      "HLT_Mu50_v",
      "HLT_TkMu50_v",
      "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",
      "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",
      "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"

    };
  }
  else{
    cout << "[SkimTree_HNMultiLep::initializeAnalyzer] DataYear is wrong : " << DataYear << endl;
  }

  cout << "[SkimTree_HNMultiLep::initializeAnalyzer] triggers to skim = " << endl;
  for(unsigned int i=0; i<triggers.size(); i++){
    cout << "[SkimTree_HNMultiLep::initializeAnalyzer]   " << triggers.at(i) << endl;
  }

}

void SkimTree_HNMultiLep::executeEvent(){

  Event ev;
  ev.SetTrigger(*HLT_TriggerName);

  //==== Skim 1 ) trigger
  if(! (ev.PassTrigger(triggers)) ) return;

  if(this->DataStream == "SingleMuon" && (ev.PassTrigger(triggers_dimu))) return;
  
  std::vector<Muon>     muonPreColl     = GetMuons("HNLoosest", 7., 2.4);
  std::vector<Electron> electronPreColl = GetElectrons("HNLoosest", 8., 2.5);

  std::sort(muonPreColl.begin(), muonPreColl.end(), PtComparing);
  std::sort(electronPreColl.begin(), electronPreColl.end(), PtComparing);
  
  int NEl  = electronPreColl.size();
  int NMu  = muonPreColl.size();
  int NLep = NEl+NMu;
  bool HasSS2lOR3l = false;
  bool LeadLepPt = false;

  bool HasFatJet(false);
  
  vector<FatJet> allfatjets = puppiCorr->Correct( GetFatJets("tight", 200., 2.7) ); //==== corret SDMass                                                     

  HasFatJet = (allfatjets.size() > 0);

  std::vector<Lepton *> leps;
  for(unsigned int i=0; i<electronPreColl.size(); i++){
    Electron& el = electronPreColl.at(i);
    leps.push_back( &el );
  }
  for(unsigned int i=0; i<muonPreColl.size(); i++){
    Muon& mu = muonPreColl.at(i);
    leps.push_back( &mu );
  }

  if      ( NLep >= 3 ){ 
    HasSS2lOR3l = true; 
    if(NMu==0 ){
      if(electronPreColl.at(0).Pt() > 23.) LeadLepPt = true;
    }  
    else if(NEl==0){
      if(muonPreColl.at(0).Pt() > 17.) LeadLepPt = true;
    }
    else LeadLepPt = true;
    
  }
  else if ( NLep == 2 ){
    if(muonPreColl.size()==2){
      if (muonPreColl[0].Charge() == muonPreColl[1].Charge())   HasSS2lOR3l = true;
      else if (HasFatJet) HasSS2lOR3l = true;
    }
    else if(electronPreColl.size()==2){
      
      if( electronPreColl[0].Charge() == electronPreColl[1].Charge() ) HasSS2lOR3l = true;
      else if (HasFatJet) HasSS2lOR3l = true;
    }

    else if(electronPreColl.size() == 1 &&  muonPreColl.size() == 1){
      if( electronPreColl[0].Charge() ==  muonPreColl[0].Charge()) HasSS2lOR3l = true;
      else if (HasFatJet) HasSS2lOR3l = true;
    }
  }
  
  if(NMu==2 && muonPreColl.at(0).Pt()>  17.    ) LeadLepPt = true;
  if(NEl==2 && electronPreColl.at(0).Pt()>23 ) LeadLepPt = true;
  if(NMu==1 && NEl==1 && electronPreColl.at(0).Pt()>23 ) LeadLepPt = true;
  if(NMu==1 && NEl==1 && muonPreColl.at(0).Pt()>23 ) LeadLepPt = true;

  if( !(HasSS2lOR3l && LeadLepPt) ) return;

  //=============================
  //==== If survived, fill tree
  //=============================

  newtree->Fill();

}

void SkimTree_HNMultiLep::executeEventFromParameter(AnalyzerParameter param){

}

SkimTree_HNMultiLep::SkimTree_HNMultiLep(){

  newtree = NULL;

}

SkimTree_HNMultiLep::~SkimTree_HNMultiLep(){

}

void SkimTree_HNMultiLep::WriteHist(){

  outfile->mkdir("recoTree");
  outfile->cd("recoTree");
  newtree->Write();
  outfile->cd();

}
