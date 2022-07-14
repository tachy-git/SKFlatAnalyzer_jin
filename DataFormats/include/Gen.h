#ifndef Gen_h
#define Gen_h

#include "Particle.h"

class Gen: public Particle {
public:

  Gen();
  ~Gen();

  //==== For empty Gen
  void SetIsEmpty(bool b);
  inline bool IsEmpty() const { return j_IsEmpty; }

  void SetIndexPIDStatus(int i, int p, int t);
  inline int Index() const { return j_index; }
  inline int PID() const { return j_PID; }
  inline int Status() const { return j_status; }
  void SetMother(int mindex);
  inline int MotherIndex() const { return j_mother_index; }

  //==== Gen Status Flags
  void SetGenStatusFlag_isPrompt(int i);
  void SetGenStatusFlag_isPromptFinalState(int i);
  void SetGenStatusFlag_isTauDecayProduct(int i);
  void SetGenStatusFlag_isPromptTauDecayProduct(int i);
  void SetGenStatusFlag_isDirectPromptTauDecayProductFinalState(int i);
  void SetGenStatusFlag_isHardProcess(int i);
  void SetGenStatusFlag_isLastCopy(int i);
  void SetGenStatusFlag_isLastCopyBeforeFSR(int i);
  void SetGenStatusFlag_isPromptDecayed(int i);
  void SetGenStatusFlag_isDecayedLeptonHadron(int i);
  void SetGenStatusFlag_fromHardProcessBeforeFSR(int i);
  void SetGenStatusFlag_fromHardProcessDecayed(int i);
  void SetGenStatusFlag_fromHardProcessFinalState(int i);
  void SetGenStatusFlag_isMostlyLikePythia6Status3(int i);

  inline bool isPrompt() const { return j_isPrompt; }
  inline bool isPromptFinalState() const { return j_isPromptFinalState; }
  inline bool isTauDecayProduct() const { return j_isTauDecayProduct; }
  inline bool isPromptTauDecayProduct() const { return j_isPromptTauDecayProduct; }
  inline bool isDirectPromptTauDecayProductFinalState() const { return j_isDirectPromptTauDecayProductFinalState; }
  inline bool isHardProcess() const { return j_isHardProcess; }
  inline bool isLastCopy() const { return j_isLastCopy; }
  inline bool isLastCopyBeforeFSR() const { return j_isLastCopyBeforeFSR; }
  inline bool isPromptDecayed() const { return j_isPromptDecayed; }
  inline bool isDecayedLeptonHadron() const { return j_isDecayedLeptonHadron; }
  inline bool fromHardProcessBeforeFSR() const { return j_fromHardProcessBeforeFSR; }
  inline bool fromHardProcessDecayed() const { return j_fromHardProcessDecayed; }
  inline bool fromHardProcessFinalState() const { return j_fromHardProcessFinalState; }
  inline bool isMostlyLikePythia6Status3() const { return j_isMostlyLikePythia6Status3; }

  virtual void Print() const;

private:

  bool j_IsEmpty;
  int j_index, j_PID, j_status;
  int j_mother_index;
  bool j_isPrompt, j_isPromptFinalState, j_isTauDecayProduct, j_isPromptTauDecayProduct, j_isDirectPromptTauDecayProductFinalState, j_isHardProcess, j_isLastCopy, j_isLastCopyBeforeFSR, j_isPromptDecayed, j_isDecayedLeptonHadron, j_fromHardProcessBeforeFSR, j_fromHardProcessDecayed, j_fromHardProcessFinalState, j_isMostlyLikePythia6Status3;

  ClassDef(Gen,1)
};

#endif
