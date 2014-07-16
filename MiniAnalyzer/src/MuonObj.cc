#include "Test/MiniAnalyzer/interface/MuonObj.h"
ClassImp(MuonObj)

ostream & operator<< (ostream &out, const MuonObj &o)
{
  out<<(ParticleObj)o <<" MuonObj: ";
  if (o.isTracker()) out << "_TRK";
  if (o.isOuter())   out << "_OUT";
  if (o.isGlobal())  out << "_GLB";
  if (o.isPF())      out << "_PF";
  if (o.isLoose())   out << "_Loose";
  if (o.isTight())   out << "_Tight";
  if (o.isGood())    out << "_Good";

  return out;
}
