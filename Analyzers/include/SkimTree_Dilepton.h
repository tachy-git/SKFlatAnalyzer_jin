#ifndef SkimTree_Dilepton_h
#define SkimTree_Dilepton_h

#include "AnalyzerCore.h"

class SkimTree_Dilepton : public AnalyzerCore {

public:

  void initializeAnalyzer();
  void executeEventFromParameter(AnalyzerParameter param);
  void executeEvent();

  SkimTree_Dilepton();
  ~SkimTree_Dilepton();

  TTree *newtree;

  vector<TString> double_triggers;
  vector<TString> single_muon_triggers;
  vector<TString> single_electron_triggers;
  void WriteHist();

};



#endif
