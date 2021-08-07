#include "Jet.h"

ClassImp(Jet)

Jet::Jet() : Particle() {
  j_area=-999.;
  j_partonFlavour=-999;
  j_hadronFlavour=-999;
  j_CSVv2=-999.;
  j_DeepCSV=-999.;
  j_DeepFlavour_b=-999;
  j_DeepFlavour_bb=-999;
  j_DeepFlavour_lepb=-999;
  j_DeepFlavour_c=-999;
  j_DeepFlavour_uds=-999;
  j_DeepFlavour_g=-999;
  j_CvsL=-999.;
  j_CvsB=-999.;
  j_DeepCvsL=-999.;
  j_DeepCvsB=-999.;
  j_chargedHadronEnergyFraction=-999.;
  j_neutralHadronEnergyFraction=-999.;
  j_neutralEmEnergyFraction=-999.;
  j_chargedEmEnergyFraction=-999.;
  j_muonEnergyFraction=-999.;
  j_chargedMultiplicity=-999;
  j_neutralMultiplicity=-999;
  j_PileupJetId=-999.;
  j_En_up=1.;
  j_En_down=1.;;
  j_Res_up = 1.;
  j_Res_down = 1.;
  j_tightJetID=false;
  j_tightLepVetoJetID=false;
}

Jet::~Jet(){

}

void Jet::SetArea(double area){
  j_area = area;
}
void Jet::SetGenFlavours(double pf, double hf){
  j_partonFlavour = pf;
  j_hadronFlavour = hf;
}
void Jet::SetTaggerResults(std::vector<double> ds){
  j_CSVv2             = ds.at(0);
  j_DeepCSV           = ds.at(1);
  j_DeepCvsL          = ds.at(2);
  j_DeepCvsB          = ds.at(3);
  j_DeepFlavour_b     = ds.at(4);
  j_DeepFlavour_bb    = ds.at(5);
  j_DeepFlavour_lepb  = ds.at(6);
  j_DeepFlavour_c     = ds.at(7);
  j_DeepFlavour_uds   = ds.at(8);
  j_DeepFlavour_g     = ds.at(9);
  j_CvsL              = ds.at(10);
  j_CvsB              = ds.at(11);
}
void Jet::SetEnergyFractions(double cH, double nH, double nEM, double cEM, double muE){
  j_chargedHadronEnergyFraction = cH;
  j_neutralHadronEnergyFraction = nH;
  j_neutralEmEnergyFraction = nEM;
  j_chargedEmEnergyFraction = cEM;
  j_muonEnergyFraction = muE;
}
void Jet::SetMultiplicities(double cM, double nM){
  j_chargedMultiplicity = cM;
  j_neutralMultiplicity = nM;
}
void Jet::SetPileupJetId(double v){
  j_PileupJetId = v;
}

void Jet::SetEnShift(double en_up, double en_down){
  j_En_up = en_up;
  j_En_down = en_down;
}

void Jet::SetResShift(double res_up, double res_down){
  j_Res_up = res_up;
  j_Res_down = res_down;
}

void Jet::SetTightJetID(double b){
  j_tightJetID = b;
}
void Jet::SetTightLepVetoJetID(double b){
  j_tightLepVetoJetID = b;
}

bool Jet::PassID(TString ID) const {

  if(ID=="tight") return Pass_tightJetID();
  if(ID=="tightLepVeto") return Pass_tightLepVetoJetID();

  cout << "[Jet::PassID] No id : " << ID << endl;
  exit(ENODATA);

  return false;

}

double Jet::GetTaggerResult(JetTagging::Tagger tg) const {

  if(tg==JetTagging::CSVv2) return j_CSVv2;
  else if(tg==JetTagging::DeepCSV) return j_DeepCSV;
  else if(tg==JetTagging::DeepJet) return j_DeepFlavour_b+j_DeepFlavour_bb+j_DeepFlavour_lepb;
  else if(tg==JetTagging::DeepFlavour_b) return j_DeepFlavour_b;
  else if(tg==JetTagging::DeepFlavour_bb) return j_DeepFlavour_bb;
  else if(tg==JetTagging::DeepFlavour_lepb) return j_DeepFlavour_lepb;
  else if(tg==JetTagging::DeepFlavour_c) return j_DeepFlavour_c;
  else if(tg==JetTagging::DeepFlavour_uds) return j_DeepFlavour_uds;
  else if(tg==JetTagging::DeepFlavour_g) return j_DeepFlavour_g;
  else if(tg==JetTagging::CvsL) return j_CvsL;
  else if(tg==JetTagging::CvsB) return j_CvsB;
  else if(tg==JetTagging::DeepCvsL) return j_DeepCvsL;
  else if(tg==JetTagging::DeepCvsB) return j_DeepCvsB;
  else{
    cout << "[Jet::GetTaggerResult] ERROR; Wrong tagger : " << tg << endl;
    return -999;
  }
}

