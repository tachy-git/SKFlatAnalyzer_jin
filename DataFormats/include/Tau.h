#ifndef Tau_h
#define Tau_h

#include "Lepton.h"
#include "TString.h"

class Tau: public Lepton{

public:

  Tau();
  ~Tau();



  //  inline double dr03TkSumPt() const { return j_dr03TkSumPt; }
  //tau_decaymode


  inline int DecayMode() const { return j_decaymode;}
  inline bool DecayModeNewDM () const { return j_idDecayModeNewDMs;}

  // ID functions 
  bool PassID(TString ID) const;

  // Functions to set variables
  void SetIDBit(unsigned int idbit);
  void SetDecayModeNewDM(bool dm_newdm);
  void SetDecayMode(int decay_mode);


  /// ID bit 
  static const int N_SELECTOR = 27;
  enum Selector {
    DeepTau2017v2p1VSJet_VVVLoose = 1<< 0,
    DeepTau2017v2p1VSJet_VVLoose = 1<< 1,
    DeepTau2017v2p1VSJet_VLoose = 1<< 2,
    DeepTau2017v2p1VSJet_Loose = 1<< 3,
    DeepTau2017v2p1VSJet_Medium  = 1<< 4,
    DeepTau2017v2p1VSJet_Tight = 1<< 5,
    DeepTau2017v2p1VSJet_VTight = 1<< 6,
    DeepTau2017v2p1VSJet_VVTight= 1<< 7,
    
    DeepTau2017v2p1VSEl_VVVLoose = 1<< 8,
    DeepTau2017v2p1VSEl_VVLoose = 1<< 9,
    DeepTau2017v2p1VSEl_VLoose = 1<< 10,
    DeepTau2017v2p1VSEl_Loose = 1<< 11,
    DeepTau2017v2p1VSEl_Medium  = 1<< 12,
    DeepTau2017v2p1VSEl_Tight = 1<< 13,
    DeepTau2017v2p1VSEl_VTight = 1<< 14,
    DeepTau2017v2p1VSEl_VVTight= 1<< 15,

    DeepTau2017v2p1VSMu_VVVLoose = 1<< 16,
    DeepTau2017v2p1VSMu_VVLoose = 1<< 17,
    DeepTau2017v2p1VSMu_VLoose = 1<< 18,
    DeepTau2017v2p1VSMu_Loose = 1<< 19,
    DeepTau2017v2p1VSMu_Medium  = 1<< 20,
    DeepTau2017v2p1VSMu_Tight = 1<< 21,
    DeepTau2017v2p1VSMu_VTight = 1<< 22,
    DeepTau2017v2p1VSMu_VVTight= 1<< 23,

    byLooseCombinedIsolationDeltaBetaCorr3Hits = 1 << 24,
    byMediumCombinedIsolationDeltaBetaCorr3Hits = 1 << 25,
    byTightCombinedIsolationDeltaBetaCorr3Hits = 1 << 26

  };

  inline bool PassSelector( unsigned int s ) const { return (j_IDBit & s)==s; }

  inline bool passVVLIDvJet()   const {return PassSelector(DeepTau2017v2p1VSJet_VVLoose); }
  inline bool passVLIDvJet()    const {return PassSelector(DeepTau2017v2p1VSJet_VLoose); }
  inline bool passLIDvJet()     const {return PassSelector(DeepTau2017v2p1VSJet_Loose); }
  inline bool passMIDvJet()     const {return PassSelector(DeepTau2017v2p1VSJet_Medium); }
  inline bool passTIDvJet()     const {return PassSelector(DeepTau2017v2p1VSJet_Tight); }
  inline bool passVTIDvJet()    const {return PassSelector(DeepTau2017v2p1VSJet_VTight); }
  inline bool passVVTIDvJet()   const {return PassSelector(DeepTau2017v2p1VSJet_VVTight); }
  

  inline bool passVVLIDvEl()   const {return PassSelector(DeepTau2017v2p1VSEl_VVLoose); }
  inline bool passVLIDvEl()    const {return PassSelector(DeepTau2017v2p1VSEl_VLoose); }
  inline bool passLIDvEl()     const {return PassSelector(DeepTau2017v2p1VSEl_Loose); }
  inline bool passMIDvEl()     const {return PassSelector(DeepTau2017v2p1VSEl_Medium); }
  inline bool passTIDvEl()     const {return PassSelector(DeepTau2017v2p1VSEl_Tight); }
  inline bool passVTIDvEl()    const {return PassSelector(DeepTau2017v2p1VSEl_VTight); }
  inline bool passVVTIDvEl()   const {return PassSelector(DeepTau2017v2p1VSEl_VVTight); }

  
  inline bool passVVLIDvMu()   const {return PassSelector(DeepTau2017v2p1VSMu_VVLoose); }
  inline bool passVLIDvMu()    const {return PassSelector(DeepTau2017v2p1VSMu_VLoose); }
  inline bool passLIDvMu()     const {return PassSelector(DeepTau2017v2p1VSMu_Loose); }
  inline bool passMIDvMu()     const {return PassSelector(DeepTau2017v2p1VSMu_Medium); }
  inline bool passTIDvMu()     const {return PassSelector(DeepTau2017v2p1VSMu_Tight); }
  inline bool passVTIDvMu()    const {return PassSelector(DeepTau2017v2p1VSMu_VTight); }
  inline bool passVVTIDvMu()   const {return PassSelector(DeepTau2017v2p1VSMu_VVTight); }
					  
  inline bool passCombIsoLoose()   const {return PassSelector(byLooseCombinedIsolationDeltaBetaCorr3Hits); }
  inline bool passCombIsoMedium()  const {return PassSelector(byMediumCombinedIsolationDeltaBetaCorr3Hits); }
  inline bool passCombIsoTight()   const {return PassSelector(byTightCombinedIsolationDeltaBetaCorr3Hits); }


private:

  unsigned int j_IDBit;
  
  int j_decaymode;
  bool j_idDecayModeNewDMs;

  ClassDef(Tau,1)
    
};

#endif
