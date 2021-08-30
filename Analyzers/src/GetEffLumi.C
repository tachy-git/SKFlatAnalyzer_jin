#include "GetEffLumi.h"

void GetEffLumi::initializeAnalyzer(){
  fChain->SetBranchStatus("*",0);
  fChain->SetBranchStatus("IsData",1);
  if(!IsDATA){
    fChain->SetBranchStatus("gen_weight",1);
    fChain->SetBranchStatus("weight_Scale",1);
  }
}

void GetEffLumi::executeEvent(){

  double MCweight = 1.;
  if(!IsDATA) MCweight = gen_weight>0 ? 1. : -1. ;

  FillHist("sumW", 0, MCweight, 1, 0., 1.);
  if(!IsDATA && weight_Scale->size()>0){
    FillHist("sumW_Reweight", 0, MCweight*weight_Scale->at(0), 1, 0., 1.);
  }

}
