#include "FatJet.h"

ClassImp(FatJet)

FatJet::FatJet() : Particle() {
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
  j_lsf=-999.;
  j_lsf_pid=-999;
  j_En_up=1.;
  j_En_down=1.;;
  j_Res_up = 1.;
  j_Res_down = 1.;
  j_tightJetID=false;
  j_tightLepVetoJetID=false;
  j_puppi_tau1 = -999.;
  j_puppi_tau2 = -999.;
  j_puppi_tau3 = -999.;
  j_puppi_tau4 = -999.;
  j_SDMass = -999.;
}

FatJet::~FatJet(){

}

void FatJet::SetArea(double area){
  j_area = area;
}
void FatJet::SetGenFlavours(double pf, double hf){
  j_partonFlavour = pf;
  j_hadronFlavour = hf;
}
void FatJet::SetTaggerResults(std::vector<double> ds){
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
void FatJet::SetEnergyFractions(double cH, double nH, double nEM, double cEM, double muE){
  j_chargedHadronEnergyFraction = cH;
  j_neutralHadronEnergyFraction = nH;
  j_neutralEmEnergyFraction = nEM;
  j_chargedEmEnergyFraction = cEM;
  j_muonEnergyFraction = muE;
}
void FatJet::SetMultiplicities(double cM, double nM){
  j_chargedMultiplicity = cM;
  j_neutralMultiplicity = nM;
}

void FatJet::SetLSF(double lsf, int lsf_PID){
  j_lsf = lsf;
  j_lsf_pid = lsf_PID;
}
double FatJet::LSF() const{
  return j_lsf;
}
double FatJet::LSF_PID() const{
  return j_lsf_pid;
}

void FatJet::SetEnShift(double en_up, double en_down){
  j_En_up = en_up;
  j_En_down = en_down;
}

void FatJet::SetResShift(double res_up, double res_down){
  j_Res_up = res_up;
  j_Res_down = res_down;
}

void FatJet::SetTightJetID(double b){
  j_tightJetID = b;
}
void FatJet::SetTightLepVetoJetID(double b){
  j_tightLepVetoJetID = b;
}

bool FatJet::PassID(TString ID) const {

  if(ID=="tight") return Pass_tightJetID();
  if(ID=="tightLepVeto") return Pass_tightLepVetoJetID();

  cout << "[FatJet::PassID] No id : " << ID << endl;
  exit(ENODATA);

  return false;

}

double FatJet::GetTaggerResult(JetTagging::Tagger tg) const{
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
    cout << "[FatJet::GetTaggerResult] ERROR; Wrong tagger" << endl;
    return -999;
  }
}

void FatJet::SetPuppiTaus(double t1, double t2, double t3, double t4){
  j_puppi_tau1 = t1;
  j_puppi_tau2 = t2;
  j_puppi_tau3 = t3;
  j_puppi_tau4 = t4;
}

void FatJet::SetSDMass(double m){
  j_SDMass = m;
}

