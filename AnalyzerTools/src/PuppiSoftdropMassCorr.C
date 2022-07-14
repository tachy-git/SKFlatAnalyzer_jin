#include "PuppiSoftdropMassCorr.h"

PuppiSoftdropMassCorr::PuppiSoftdropMassCorr(){

  histDir = TDirectoryHelper::GetTempDirectory("PuppiSoftdropMassCorr");

}

void PuppiSoftdropMassCorr::ReadHistograms(){

  TString datapath = getenv("DATA_DIR");

  TDirectory* origDir = gDirectory;

  TString rootfilepath = datapath+"/"+GetEra()+"/PuppiSoftdropMassCorr/puppiCorr.root";

  TFile *rootfile = new TFile(rootfilepath);

  histDir->cd();
  func_gen = (TF1 *)rootfile->Get("puppiJECcorr_gen")->Clone();
  func_reco_0eta1v3 = (TF1 *)rootfile->Get("puppiJECcorr_reco_0eta1v3")->Clone();
  func_reco_1v3eta2v5 = (TF1 *)rootfile->Get("puppiJECcorr_reco_1v3eta2v5")->Clone();
  rootfile->Close();
  delete rootfile;
  origDir->cd();

}

PuppiSoftdropMassCorr::~PuppiSoftdropMassCorr(){

}

std::vector<FatJet> PuppiSoftdropMassCorr::Correct(std::vector<FatJet> jets){

  std::vector<FatJet> out;
  for(unsigned int i=0; i<jets.size(); i++){

    FatJet this_jet = jets.at(i);

    float genCorr  = 1.;
    float recoCorr = 1.;
    float totalWeight = 1.;

    genCorr = func_gen->Eval( this_jet.Pt() );
    if( fabs(this_jet.Eta())  <= 1.3 ){
      recoCorr = func_reco_0eta1v3->Eval( this_jet.Pt() );
    }
    else{
      recoCorr = func_reco_1v3eta2v5->Eval( this_jet.Pt() );
    }
    
    totalWeight = genCorr * recoCorr;

    this_jet.SetSDMass( this_jet.SDMass() * totalWeight );

    out.push_back( this_jet );

  }

  return out;

}

