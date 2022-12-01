#ifndef TutorialBase_h
#define TutorialBase_h

#include "AnalyzerCore.h"

class TutorialBase : public AnalyzerCore {

public:

  void initializeAnalyzer();
  
  void executeEventFromParameter(AnalyzerParameter param);
  void executeEvent();

  TString IsoMuTriggerName;
  double TriggerSafePtCut;

  vector<TString> MuonIDs, MuonIDSFKeys;
  vector<Muon> AllMuons;
  vector<Jet> AllJets;

  TutorialBase();
  ~TutorialBase();

};

#endif

