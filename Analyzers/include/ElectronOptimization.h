#ifndef ElectronOptimization_h
#define ElectronOptimization_h

#include "AnalyzerCore.h"

class ElectronOptimization: public AnalyzerCore {
public:
    ElectronOptimization();
    ~ElectronOptimization();

    TString ElectronVetoID;

    void initializeAnalyzer();
    void executeEvent();

};

#endif