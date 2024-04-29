#ifndef ElectronOptimization_h
#define ElectronOptimization_h

#include "AnalyzerCore.h"

class ElectronOptimization: public AnalyzerCore {
private:
    TTree *Events;
    unsigned int nElectrons;
    float genWeight;
    float Pt[20];
    float scEta[20];
    float MVANoIso[20];
    float MiniRelIso[20];
    float DeltaR[20];
    float SIP3D[20];
    bool  PassMVANoIsoWP90[20];
    bool  PassMVANoIsoWPLoose[20];
    int   NearestJetFlavour[20];

public:
    ElectronOptimization();
    ~ElectronOptimization();

    TString ElectronVetoID;

    void initializeAnalyzer();
    void executeEvent();

};

#endif
