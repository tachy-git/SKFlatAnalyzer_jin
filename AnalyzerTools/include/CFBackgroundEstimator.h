#ifndef CFBackgroundEstimator_h
#define CFBackgroundEstimator_h

#include <map>
#include <vector>
#include <fstream>
#include <sstream>

#include "TFile.h"
#include "TString.h"
#include "TH1D.h"
#include "TH2D.h"

#include "Lepton.h"
#include "Muon.h"
#include "Electron.h"

#include "AnalyzerParameter.h"
#include "TDirectoryHelper.h"

class CFBackgroundEstimator{

public:

  CFBackgroundEstimator();
  ~CFBackgroundEstimator();

  TDirectory *histDir;
  void ReadHistograms();

  int DataYear;
  TString DataEra;
  void SetEra(TString era){
    DataEra=era;
    DataYear=TString(era(0,4)).Atoi();
  }
  TString GetEra() const { return DataEra; }
  int GetYear() const { return DataYear; } 

  std::map< TString, TH1D* > map_hist_Electron;
  std::map< TString, TH1D* > map_hist_Muon;

  double GetElectronCFRate(TString ID, TString key, double eta, double pt, int sys=0);
  double GetMuonCFRate(TString ID, TString key, double eta, double pt, int sys=0);

  double GetWeight(vector<Lepton *> lepptrs, AnalyzerParameter param, int sys=0);

};

#endif
