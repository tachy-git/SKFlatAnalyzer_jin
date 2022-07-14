#include "Gen.h"

ClassImp(Gen)

Gen::Gen() : Particle() {
  j_IsEmpty = true;
  j_index = -1;
  j_PID = 0;
  j_status = 0;
  j_mother_index = -1;
  j_isPrompt = false;
  j_isPromptFinalState = false;
  j_isTauDecayProduct = false;
  j_isPromptTauDecayProduct = false;
  j_isDirectPromptTauDecayProductFinalState = false;
  j_isHardProcess = false;
  j_isLastCopy = false;
  j_isLastCopyBeforeFSR = false;
  j_isPromptDecayed = false;
  j_isDecayedLeptonHadron = false;
  j_fromHardProcessBeforeFSR = false;
  j_fromHardProcessDecayed = false;
  j_fromHardProcessFinalState = false;
  j_isMostlyLikePythia6Status3 = false;
}

Gen::~Gen(){
}

void Gen::SetIsEmpty(bool b){
  j_IsEmpty = b;
}

void Gen::SetIndexPIDStatus(int i, int p, int t){
  j_index = i;
  j_PID = p;
  j_status = t;
}

void Gen::SetMother(int mindex){
  j_mother_index = mindex;
}

void Gen::SetGenStatusFlag_isPrompt(int i){
  j_isPrompt = i;
}
void Gen::SetGenStatusFlag_isPromptFinalState(int i){
  j_isPromptFinalState = i;
}
void Gen::SetGenStatusFlag_isTauDecayProduct(int i){
  j_isTauDecayProduct = i;
}
void Gen::SetGenStatusFlag_isPromptTauDecayProduct(int i){
  j_isPromptTauDecayProduct = i;
}
void Gen::SetGenStatusFlag_isDirectPromptTauDecayProductFinalState(int i){
  j_isDirectPromptTauDecayProductFinalState = i;
}
void Gen::SetGenStatusFlag_isHardProcess(int i){
  j_isHardProcess = i;
}
void Gen::SetGenStatusFlag_isLastCopy(int i){
  j_isLastCopy = i;
}
void Gen::SetGenStatusFlag_isLastCopyBeforeFSR(int i){
  j_isLastCopyBeforeFSR = i;
}
void Gen::SetGenStatusFlag_isPromptDecayed(int i){
  j_isPromptDecayed = i;
}
void Gen::SetGenStatusFlag_isDecayedLeptonHadron(int i){
  j_isDecayedLeptonHadron = i;
}
void Gen::SetGenStatusFlag_fromHardProcessBeforeFSR(int i){
  j_fromHardProcessBeforeFSR = i;
}
void Gen::SetGenStatusFlag_fromHardProcessDecayed(int i){
  j_fromHardProcessDecayed = i;
}
void Gen::SetGenStatusFlag_fromHardProcessFinalState(int i){
  j_fromHardProcessFinalState = i;
}
void Gen::SetGenStatusFlag_isMostlyLikePythia6Status3(int i){
  j_isMostlyLikePythia6Status3 = i;
}
void Gen::Print() const {
  cout << "(Index, PID, Status, MotherIndex, Pt, Eta, Phi, M) = " << Index() << "\t" << PID() << "\t" << Status() << "\t" << MotherIndex() << "\t" << Pt() << "\t" << Eta() << "\t" << Phi() << "\t" << M() << endl;
}
