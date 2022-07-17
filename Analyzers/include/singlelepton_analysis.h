#ifndef singlelepton_analysis_h
#define singlelepton_analysis_h

#include "AnalyzerCore.h"

class singlelepton_analysis : public AnalyzerCore {

public:

  void initializeAnalyzer();

  void executeEventFromParameter(AnalyzerParameter param);
  void executeEvent();

  bool RunSyst;
  bool RunNewPDF;
  bool RunXSecSyst;

  std::vector<TString> muonTightIDs, muonLooseIDs;
  std::vector<TString> electronTightIDs, electronLooseIDs;
  std::vector<TString> jetIDs, fatjetIDs;

  JetTagging::Parameters bTaggingWPMedium;
  JetTagging::Tagger pNetXbbMD, pNetXqqMD, pNetXccMD, pNetQCDMD;

  std::vector<TString> muonTriggers, electronTriggers;

  double muonPtCut, electronPtCut, jetPtCut, fatjetPtCut, leptonPtCut;

  std::vector<Muon> allMuons;
  std::vector<Electron> allElectrons;
  std::vector<Jet> allJets;
  std::vector<FatJet> allFatJets;

  double GetReconstructedNeutrinoDet(Lepton lepton, Particle missingEt);
  Particle GetReconstructedNeutrino(Lepton lepton, Particle missingEt);

  singlelepton_analysis();
  ~singlelepton_analysis();

};

#endif
