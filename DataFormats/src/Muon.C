#include "Muon.h"

ClassImp(Muon)

Muon::Muon() : Lepton() {
  j_chi2 = 999.;
  j_PFCH04 = -999.;
  j_PFNH04 = -999.;
  j_PFPH04 = -999.;
  j_PU04 = -999.;
  j_trkiso = -999.;
  this->SetLeptonFlavour(MUON);
  j_MiniAODPt = -999.;
  j_MomentumScaleUp = -999.;
  j_MomentumScaleDown = -999.;
  j_TunePPtError = -999.;
  j_MVA = -999.;
  j_lowptMVA = -999.;
  j_softMVA = -999.;
  j_trackerLayers = 0;
}

Muon::~Muon(){
}

void Muon::SetTypeBit(unsigned int typebit){
  j_TypeBit = typebit;
}

void Muon::SetIDBit(unsigned int idbit){
  j_IDBit = idbit;
}

void Muon::SetisPOGHighPt(bool b){
  j_isPOGHighPt = b;
}

void Muon::SetIso(double ch04, double nh04, double ph04, double pu04, double trkiso){
  j_PFCH04 = ch04;
  j_PFNH04 = nh04;
  j_PFPH04 = ph04;
  j_PU04 = pu04;
  j_trkiso = trkiso;
  CalcPFRelIso();
}

void Muon::SetChi2(double chi2){
  j_chi2 = chi2;
}

void Muon::CalcPFRelIso(){
  double absiso = j_PFCH04+std::max( 0., j_PFNH04 + j_PFPH04 - 0.5*j_PU04 );
  //cout << "[Muon::CalcPFRelIso] j_PFCH04 = " << j_PFCH04 << endl;
  //cout << "[Muon::CalcPFRelIso] j_PFNH04 = " << j_PFNH04 << endl;
  //cout << "[Muon::CalcPFRelIso] j_PFPH04 = " << j_PFPH04 << endl;
  //cout << "[Muon::CalcPFRelIso] j_PU04 = " << j_PU04 << endl;
  //cout << "[Muon::CalcPFRelIso] --> absiso = " << absiso << endl;
  this->SetRelIso(absiso/this->Pt());
  //this->SetRelIso(absiso/this->MiniAODPt()); //TODO This is same as IDBit
}

double Muon::EA(){

  double eta = fabs(this->Eta());

  if     (eta<0.8000) return 0.0566;
  else if(eta<1.3000) return 0.0562;
  else if(eta<2.0000) return 0.0363;
  else if(eta<2.2000) return 0.0119;
  else if(eta<2.4000) return 0.0064;
  else return 0.0064;

}

void Muon::SetMiniAODPt(double d){
  j_MiniAODPt = d;
}
void Muon::SetMiniAODTunePPt(double d){
  j_MiniAODTunePPt = d;
}

void Muon::SetMomentumScaleUpDown(double pt_up, double pt_down){
  j_MomentumScaleUp = pt_up;
  j_MomentumScaleDown = pt_down;
}

void Muon::SetTuneP4(double pt, double pt_err, double eta, double phi, double q){
  j_TuneP4.SetPtEtaPhiM(pt,eta,phi,M());
  j_TuneP4.SetCharge(q);
  j_TunePPtError = pt_err;
}

void Muon::SetMVA(double MVA){
  j_MVA = MVA;
  //j_lowptMVA = lowptMVA;
  //j_softMVA = softMVA;
}

bool Muon::PassID(TString ID) const {
  //==== POG
  if(ID=="POGTight") return isPOGTight();
  if(ID=="POGHighPt") return isPOGHighPt();
  if(ID=="POGMedium") return isPOGMedium();
  if(ID=="POGLoose") return isPOGLoose();
  if(ID=="POGTightWithTightIso") return Pass_POGTightWithTightIso();
  if(ID=="POGHighPtWithLooseTrkIso") return Pass_POGHighPtWithLooseTrkIso();
  //==== Customized
  if(ID=="TEST") return Pass_TESTID();

  //==== No cut
  if(ID=="NOCUT") return true;

  cout << "[Electron::PassID] No id : " << ID << endl;
  exit(ENODATA);

  return false;

}
bool Muon::Pass_POGTightWithTightIso() const {
  if(!( isPOGTight() )) return false;
  if(!( RelIso()<0.15 ))  return false;
  return true;
}
bool Muon::Pass_POGHighPtWithLooseTrkIso() const {
  if(!( isPOGHighPt() )) return false;
  if(!( TrkIso()/TuneP4().Pt()<0.1 )) return false;
  return true;
}

//==== TEST ID

bool Muon::Pass_TESTID() const {
  return true;
}

void Muon::SetTrackerLayers(int n){
  j_trackerLayers = n;
}

bool Muon::PassFilter(TString filter) const{
  if( filter=="hltDiMu9Ele9CaloIdLTrackIdLMuonlegL3Filtered9" && j_filterbits&(ULong64_t(1)<<0) ) return true;
  else if( filter=="hltDiMuon178Mass3p8Filtered" && j_filterbits&(ULong64_t(1)<<1) ) return true;
  else if( filter=="hltDiMuon178Mass8Filtered" && j_filterbits&(ULong64_t(1)<<2) ) return true;
  else if( filter=="hltDiMuon178RelTrkIsoFiltered0p4" && j_filterbits&(ULong64_t(1)<<3) ) return true;
  else if( filter=="hltDiMuon178RelTrkIsoFiltered0p4DzFiltered0p2" && j_filterbits&(ULong64_t(1)<<4) ) return true;
  else if( filter=="hltDiMuonGlb17Glb8DzFiltered0p2" && j_filterbits&(ULong64_t(1)<<5) ) return true;
  else if( filter=="hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4" && j_filterbits&(ULong64_t(1)<<6) ) return true;
  else if( filter=="hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4DzFiltered0p2" && j_filterbits&(ULong64_t(1)<<7) ) return true;
  else if( filter=="hltDiMuonGlb17Trk8RelTrkIsoFiltered0p4" && j_filterbits&(ULong64_t(1)<<8) ) return true;
  else if( filter=="hltDiMuonGlb17Trk8RelTrkIsoFiltered0p4DzFiltered0p2" && j_filterbits&(ULong64_t(1)<<9) ) return true;
  else if( filter=="hltDiMuonGlbFiltered17TrkFiltered8" && j_filterbits&(ULong64_t(1)<<10) ) return true;
  else if( filter=="hltDiMuonTrk17Trk8RelTrkIsoFiltered0p4" && j_filterbits&(ULong64_t(1)<<11) ) return true;
  else if( filter=="hltDiMuonTrk17Trk8RelTrkIsoFiltered0p4DzFiltered0p2" && j_filterbits&(ULong64_t(1)<<12) ) return true;
  else if( filter=="hltDiTkMuonTkFiltered17TkFiltered8" && j_filterbits&(ULong64_t(1)<<13) ) return true;
  else if( filter=="hltL1TripleMu553L2TriMuFiltered3L3TriMuFiltered10105" && j_filterbits&(ULong64_t(1)<<14) ) return true;
  else if( filter=="hltL1TripleMu553L2TriMuFiltered3L3TriMuFiltered5" && j_filterbits&(ULong64_t(1)<<15) ) return true;
  else if( filter=="hltL3crIsoL1sMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p09" && j_filterbits&(ULong64_t(1)<<16) ) return true;
  else if( filter=="hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07" && j_filterbits&(ULong64_t(1)<<17) ) return true;
  else if( filter=="hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p07" && j_filterbits&(ULong64_t(1)<<18) ) return true;
  else if( filter=="hltL3fL1DoubleMu155fFiltered17" && j_filterbits&(ULong64_t(1)<<19) ) return true;
  else if( filter=="hltL3fL1DoubleMu155fPreFiltered8" && j_filterbits&(ULong64_t(1)<<20) ) return true;
  else if( filter=="hltL3fL1DoubleMu157fFiltered18" && j_filterbits&(ULong64_t(1)<<21) ) return true;
  else if( filter=="hltL3fL1DoubleMu157fFiltered9" && j_filterbits&(ULong64_t(1)<<22) ) return true;
  else if( filter=="hltL3fL1DoubleMu7EG7f0Filtered9" && j_filterbits&(ULong64_t(1)<<23) ) return true;
  else if( filter=="hltL3fL1Mu6DoubleEG10f0Filtered8" && j_filterbits&(ULong64_t(1)<<24) ) return true;
  else if( filter=="hltL3fL1TripleMu553f0Filtered10105" && j_filterbits&(ULong64_t(1)<<25) ) return true;
  else if( filter=="hltL3fL1TripleMu553f0Filtered1055" && j_filterbits&(ULong64_t(1)<<26) ) return true;
  else if( filter=="hltL3fL1TripleMu553f0PreFiltered555" && j_filterbits&(ULong64_t(1)<<27) ) return true;
  else if( filter=="hltL3fL1sDoubleMu114L1f0L2f10L3Filtered17" && j_filterbits&(ULong64_t(1)<<28) ) return true;
  else if( filter=="hltL3fL1sDoubleMu114L1f0L2f10OneMuL3Filtered17" && j_filterbits&(ULong64_t(1)<<29) ) return true;
  else if( filter=="hltL3fL1sDoubleMu114TkFiltered17Q" && j_filterbits&(ULong64_t(1)<<30) ) return true;
  else if( filter=="hltL3fL1sMu10lqL1f0L2f10L3Filtered17" && j_filterbits&(ULong64_t(1)<<31) ) return true;
  else if( filter=="hltL3fL1sMu15DQlqL1f0L2f10L3Filtered17" && j_filterbits&(ULong64_t(1)<<32) ) return true;
  else if( filter=="hltL3fL1sMu1lqL1f0L2f10L3Filtered17TkIsoFiltered0p4" && j_filterbits&(ULong64_t(1)<<33) ) return true;
  else if( filter=="hltL3fL1sMu22L1f0Tkf24QL3trkIsoFiltered0p09" && j_filterbits&(ULong64_t(1)<<34) ) return true;
  else if( filter=="hltL3fL1sMu5L1f0L2f5L3Filtered8" && j_filterbits&(ULong64_t(1)<<35) ) return true;
  else if( filter=="hltL3fL1sMu5L1f0L2f5L3Filtered8TkIsoFiltered0p4" && j_filterbits&(ULong64_t(1)<<36) ) return true;
  else if( filter=="hltL3pfL1sDoubleMu114L1f0L2pf0L3PreFiltered8" && j_filterbits&(ULong64_t(1)<<37) ) return true;
  else if( filter=="hltL3pfL1sDoubleMu114ORDoubleMu125L1f0L2pf0L3PreFiltered8" && j_filterbits&(ULong64_t(1)<<38) ) return true;
  else if( filter=="hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLDZFilter" && j_filterbits&(ULong64_t(1)<<39) ) return true;
  else if( filter=="hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered23" && j_filterbits&(ULong64_t(1)<<40) ) return true;
  else if( filter=="hltMu23TrkIsoVVLEle8CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered23" && j_filterbits&(ULong64_t(1)<<41) ) return true;
  else if( filter=="hltMu8DiEle12CaloIdLTrackIdLMuonlegL3Filtered8" && j_filterbits&(ULong64_t(1)<<42) ) return true;
  else if( filter=="hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLDZFilter" && j_filterbits&(ULong64_t(1)<<43) ) return true;
  else if( filter=="hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered8" && j_filterbits&(ULong64_t(1)<<44) ) return true;
  else if( filter=="hltMu9Ele9DZFilter" && j_filterbits&(ULong64_t(1)<<45) ) return true;
  else if( filter=="hltTripleTrkMuFiltered5NoVtx" && j_filterbits&(ULong64_t(1)<<46) ) return true;
  else{
    cout<<"[Muon::PassFilter] unknown filter "<<filter<<endl;
    exit(ENODATA);
  }
  return false;
}

bool Muon::PassPath(TString path) const{
  if( path=="HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ_v" && j_pathbits&(ULong64_t(1)<<0) ) return true;
  else if( path=="HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v" && j_pathbits&(ULong64_t(1)<<1) ) return true;
  else if( path=="HLT_IsoMu24_v" && j_pathbits&(ULong64_t(1)<<2) ) return true;
  else if( path=="HLT_IsoMu27_v" && j_pathbits&(ULong64_t(1)<<3) ) return true;
  else if( path=="HLT_IsoTkMu24_v" && j_pathbits&(ULong64_t(1)<<4) ) return true;
  else if( path=="HLT_Mu17_Mu8_SameSign_DZ_v" && j_pathbits&(ULong64_t(1)<<5) ) return true;
  else if( path=="HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v" && j_pathbits&(ULong64_t(1)<<6) ) return true;
  else if( path=="HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v" && j_pathbits&(ULong64_t(1)<<7) ) return true;
  else if( path=="HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v" && j_pathbits&(ULong64_t(1)<<8) ) return true;
  else if( path=="HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v" && j_pathbits&(ULong64_t(1)<<9) ) return true;
  else if( path=="HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v" && j_pathbits&(ULong64_t(1)<<10) ) return true;
  else if( path=="HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v" && j_pathbits&(ULong64_t(1)<<11) ) return true;
  else if( path=="HLT_Mu17_TrkIsoVVL_v" && j_pathbits&(ULong64_t(1)<<12) ) return true;
  else if( path=="HLT_Mu17_v" && j_pathbits&(ULong64_t(1)<<13) ) return true;
  else if( path=="HLT_Mu18_Mu9_SameSign_v" && j_pathbits&(ULong64_t(1)<<14) ) return true;
  else if( path=="HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v" && j_pathbits&(ULong64_t(1)<<15) ) return true;
  else if( path=="HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v" && j_pathbits&(ULong64_t(1)<<16) ) return true;
  else if( path=="HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v" && j_pathbits&(ULong64_t(1)<<17) ) return true;
  else if( path=="HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v" && j_pathbits&(ULong64_t(1)<<18) ) return true;
  else if( path=="HLT_Mu50_v" && j_pathbits&(ULong64_t(1)<<19) ) return true;
  else if( path=="HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v" && j_pathbits&(ULong64_t(1)<<20) ) return true;
  else if( path=="HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v" && j_pathbits&(ULong64_t(1)<<21) ) return true;
  else if( path=="HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v" && j_pathbits&(ULong64_t(1)<<22) ) return true;
  else if( path=="HLT_Mu8_TrkIsoVVL_v" && j_pathbits&(ULong64_t(1)<<23) ) return true;
  else if( path=="HLT_Mu8_v" && j_pathbits&(ULong64_t(1)<<24) ) return true;
  else if( path=="HLT_OldMu100_v" && j_pathbits&(ULong64_t(1)<<25) ) return true;
  else if( path=="HLT_TkMu100_v" && j_pathbits&(ULong64_t(1)<<26) ) return true;
  else if( path=="HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v" && j_pathbits&(ULong64_t(1)<<27) ) return true;
  else if( path=="HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v" && j_pathbits&(ULong64_t(1)<<28) ) return true;
  else if( path=="HLT_TkMu50_v" && j_pathbits&(ULong64_t(1)<<29) ) return true;
  else if( path=="HLT_TripleMu_10_5_5_DZ_v" && j_pathbits&(ULong64_t(1)<<30) ) return true;
  else if( path=="HLT_TripleMu_12_10_5_v" && j_pathbits&(ULong64_t(1)<<31) ) return true;
  else if( path=="HLT_TrkMu12_DoubleTrkMu5NoFiltersNoVtx_v" && j_pathbits&(ULong64_t(1)<<32) ) return true;
  else if( path=="HLT_TrkMu15_DoubleTrkMu5NoFiltersNoVtx_v" && j_pathbits&(ULong64_t(1)<<33) ) return true;
  else{
    cout<<"[Muon::PassPath] unknown path "<<path<<endl;
    exit(ENODATA);
  }
  return false;
}
