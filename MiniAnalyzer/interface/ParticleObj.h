#ifndef MiniAnalyzer_ParticleObj_H
#define MiniAnalyzer_ParticleObj_H

#include "TObject.h"
#include <ostream>
#include <vector>

class ParticleObj : public TObject {

 public:
 
   ParticleObj(float pt=0., float eta=0., float phi=0., float mass=0., int charge=0) : 
  thePt(pt), theEta(eta), thePhi(phi), theMass(mass), theCharge(charge) {}
   void setKine(float pt, float eta, float phi, float mass, int charge) { 
    thePt=pt; theEta=eta; thePhi=phi; theMass=mass; theCharge=charge;}

   virtual ~ParticleObj(){}

   float pt() const { return thePt;}
   float eta() const { return theEta;}
   float phi() const { return thePhi;}
   float mass() const { return theMass;}
   int charge() const { return theCharge;}

 private:
   float thePt, theEta, thePhi, theMass;
   int theCharge;
   friend ostream & operator<< (ostream &out, const ParticleObj &o);

 public:
   ClassDef(ParticleObj,1)

};

#endif
