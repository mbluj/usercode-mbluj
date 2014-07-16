#include "Test/MiniAnalyzer/interface/DiTauObj.h"
ClassImp(DiTauObj)

ostream & operator<< (ostream &out, const DiTauObj &o)
{
  out<<(ParticleObj)o <<" DiTauObj: ";
  out <<" svFitMass: "<<o.svFitMass;
  out <<" svFitPt: "  <<o.svFitPt;

  return out;
}
