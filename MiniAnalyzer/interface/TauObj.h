#ifndef MiniAnalyzer_TauObj_H
#define MiniAnalyzer_TauObj_H

#include "Test/MiniAnalyzer/interface/ParticleObj.h"
#include <ostream>
#include <vector>

class TauObj : public ParticleObj {

 public:
   TauObj(float pt=0., float eta=0., float phi=0., float mass=0., int charge=0) 
     : ParticleObj(pt,eta,phi,mass,charge), 
    nAllTaus(0),
    chIso(0), phIso(0), puIso(0), iso(0),
    dzPV(0),
    dzTau(0), dzTauErr(0), dxyTau(0), dxyTauErr(0), dxyTauSig(0),
    dzTrk(0), dzTrkErr(0), dxyTrk(0), dxyTrkErr(0),
    hasKft(-1), hasSecVtx(false),
    flightLength(-1), flightLengthError(-1), flightLengthSig(-1),
    genPdgId(0), nGenPart(0), leadGenPt(0), genTauJetPt(0), 
    decMode(-1), genDecMode(-1),
    visMass(mass), genVisMass(-1),
    tauIso(-1), antiE(-1), antiMu(-1) {}
  virtual ~TauObj(){}

 public:
  unsigned int nAllTaus;
  float chIso, phIso, puIso, iso;
  float dzPV;
  float dzTau, dzTauErr, dxyTau, dxyTauErr, dxyTauSig;
  float dzTrk, dzTrkErr, dxyTrk, dxyTrkErr;
  int hasKft;
  bool hasSecVtx;
  float flightLength, flightLengthError, flightLengthSig;
  int genPdgId;
  unsigned int nGenPart;
  float leadGenPt, genTauJetPt;
  int decMode, genDecMode;
  float visMass, genVisMass;
  int tauIso, antiE, antiMu;

 private:  
  friend ostream & operator<< (ostream &out, const TauObj &o);

 public:
  ClassDef(TauObj,2)

};

#endif
