#include "Test/MiniAnalyzer/interface/TauObj.h"
ClassImp(TauObj)

ostream & operator<< (ostream &out, const TauObj &o)
{
  out<<(ParticleObj)o  <<" TauObj: ";
  out <<" decayMode: " <<o.decMode;
  out <<" 3HitIso: "   <<o.tauIso%10;
  out <<" MVAIsoWoLt: "<<(o.tauIso%100-o.tauIso%10)/10;
  out <<" MVAIsoWLt: " <<(o.tauIso%1000-o.tauIso%100)/100;
  out <<" #AllTaus: "  <<o.nAllTaus;

  return out;
}
