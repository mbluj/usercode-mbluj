#ifndef MiniAnalyzer_DiTauObj_H
#define MiniAnalyzer_DiTauObj_H

#include "Test/MiniAnalyzer/interface/ParticleObj.h"
#include <ostream>
#include <vector>

class DiTauObj : public ParticleObj {

 public:
   DiTauObj(float pt=0., float eta=0., float phi=0., float mass=0., int charge=0) 
     : ParticleObj(pt,eta,phi,mass,charge), 
    leg1Idx(-1), leg1PdgId(0),
    leg2Idx(-1), leg2PdgId(0),
    metIdx(-1),
    svFitMass(-1), svFitMassErrUp(-1), svFitMassErrDown(-1),
    svFitPt(-1), svFitPtErrUp(-1), svFitPtErrDown(-1)
    {}
  virtual ~DiTauObj(){}

 public:
  int leg1Idx, leg1PdgId;
  int leg2Idx, leg2PdgId;
  int metIdx;
  float svFitMass, svFitMassErrUp, svFitMassErrDown;
  float svFitPt,   svFitPtErrUp,   svFitPtErrDown;


 private:  
  friend ostream & operator<< (ostream &out, const DiTauObj &o);

 public:
  ClassDef(DiTauObj,2)

};

#endif
