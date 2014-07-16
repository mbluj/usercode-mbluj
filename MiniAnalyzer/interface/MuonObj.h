#ifndef MiniAnalyzer_MuonObj_H
#define MiniAnalyzer_MuonObj_H

#include "Test/MiniAnalyzer/interface/ParticleObj.h"
#include <ostream>
#include <vector>

class MuonObj : public ParticleObj {

 public:
    MuonObj(float pt=0., float eta=0., float phi=0., float mass=0., int charge=0) 
      : ParticleObj(pt,eta,phi,mass,charge), 
    nRPCHits(0), nDTHits(0), nCSCHits(0), nTrackerHits(0), nMatchedStations(0),
    isUnique(true), nAllMuons(0), 
    chHadIso(0), chIso(0), nHadIso(0), phIso(0), puIso(0),
    dz(0), dzErr(0), dxy(0), dxyErr(0), genPdgId(0), genPt(-1),
    theMuonBits(0) {}
  virtual ~MuonObj(){}
  void setBits(bool isGlobal, bool isTracker, bool isOuter, bool isCalo, bool isMatched, 
	       bool isPF, bool isLoose, bool isTight, bool isGood) {
    if (isGood)    theMuonBits  = 1<<8;
    if (isLoose)   theMuonBits |= 1<<7;
    if (isTight)   theMuonBits |= 1<<6;
    if (isPF)      theMuonBits |= 1<<5;
    if (isGlobal)  theMuonBits |= 1<<4;
    if (isTracker) theMuonBits |= 1<<3;
    if (isOuter)   theMuonBits |= 1<<2;
    if (isCalo)    theMuonBits |= 1<<1;
    if (isMatched) theMuonBits |= 1; 
  }
  bool isGood()    const { return  (theMuonBits>>8)&1 ;}
  bool isLoose()   const { return  (theMuonBits>>7)&1 ;}
  bool isTight()   const { return  (theMuonBits>>6)&1 ;}
  bool isPF()      const { return  (theMuonBits>>5)&1 ;}  
  bool isGlobal()  const { return  (theMuonBits>>4)&1 ;}  
  bool isTracker() const { return  (theMuonBits>>3)&1 ;}  
  bool isOuter()   const { return  (theMuonBits>>2)&1 ;}  
  bool isCalo()    const { return  (theMuonBits>>1)&1 ;}  
  bool isMatched() const { return   theMuonBits&1 ;}

 public:
  unsigned int nRPCHits, nDTHits, nCSCHits, nTrackerHits, nMatchedStations;
  bool         isUnique;
  unsigned int nAllMuons;
  float chHadIso, chIso, nHadIso, phIso, puIso;
  float dz, dzErr, dxy, dxyErr;
  int genPdgId;
  float genPt;

 private:  
  unsigned int theMuonBits; 
  friend ostream & operator<< (ostream &out, const MuonObj &o);

 public:
  ClassDef(MuonObj,2)

};

#endif
