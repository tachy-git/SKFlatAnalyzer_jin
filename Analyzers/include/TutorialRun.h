#ifndef TutorialRun_h
#define TutorialRun_h

#include "AnalyzerCore.h"

class TutorialRun : public AnalyzerCore {

public:

  void initializeAnalyzer();
  
  void executeEventFromParameter(AnalyzerParameter param);
  void executeEvent();

  TString IsoMuTriggerName;
  double TriggerSafePtCut;

  vector<TString> MuonIDs, MuonIDSFKeys;
  vector<Muon> AllMuons;
  vector<Jet> AllJets;

  TutorialRun();
  ~TutorialRun();

};

#endif

