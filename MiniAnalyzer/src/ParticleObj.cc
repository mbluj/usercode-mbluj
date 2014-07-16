#include "Test/MiniAnalyzer/interface/ParticleObj.h"

ClassImp(ParticleObj)

ostream & operator<< (ostream &out, const ParticleObj &o){
  out<<"ParticleObj: ";
  out <<" pt: "<<o.thePt<<", eta: "<<o.theEta<<", phi: "<<o.thePhi<<", mass: "<<o.theMass<<", charge: "<<o.theCharge;

  return out;
}
