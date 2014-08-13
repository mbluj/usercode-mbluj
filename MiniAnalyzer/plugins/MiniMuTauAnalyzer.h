#ifndef MiniAnalyzer_Plugins_MiniMuTauAnalyzer_H
#define MiniAnalyzer_Plugins_MiniMuTauAnalyzer_H

// -*- C++ -*-
//
// Package:    Test/MiniAnalyzer
// Class:      MiniMuTauAnalyzer
// 
/**\class MiniMuTauAnalyzer MiniMuTauAnalyzer.cc Test/MiniAnalyzer/plugins/MiniMuTauAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Michal Bluj
//         Created:  Tue, 08 Jul 2014 10:44:22 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "Test/MiniAnalyzer/interface/EventObj.h"
#include "Test/MiniAnalyzer/interface/MuonObj.h"
#include "Test/MiniAnalyzer/interface/TauObj.h"
#include "Test/MiniAnalyzer/interface/JetObj.h"
#include "Test/MiniAnalyzer/interface/MetObj.h"
#include "Test/MiniAnalyzer/interface/DiTauObj.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"

#include <vector>
#include <string>
#include <utility>
#include <map>
#include <math.h>

//
// class declaration
//

class MiniMuTauAnalyzer : public edm::EDAnalyzer {

public:
  explicit MiniMuTauAnalyzer(const edm::ParameterSet&);
  ~MiniMuTauAnalyzer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  
  //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

  bool muonId(const pat::Muon&, bool checkLoose=false);
  int findGenDecayMode(const reco::GenJet&);
  int findDecayMode(const pat::Tau&);

  int match(const reco::Candidate *, const edm::View<reco::Candidate> *,
	    float, int&, float dPtmax=-1, float Ptmin=0, bool byPt=false);

  // ----------member data ---------------------------
private:
  // tokens
  edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
  edm::EDGetTokenT<pat::MuonCollection> muonToken_;
  edm::EDGetTokenT<pat::TauCollection> tauToken_;
  edm::EDGetTokenT<pat::METCollection> metToken_;
  edm::EDGetTokenT<pat::JetCollection> jetToken_;

  edm::EDGetTokenT<edm::View<reco::Candidate> > genTauJetsToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> prunedGenParticlesToken_;

  TFile* file_;
  TTree* tree_;

  std::vector<MuonObj> *myMuons_;
  std::vector<TauObj> *myTaus_;
  std::vector<JetObj> *myJets_;
  std::vector<MetObj> *myMets_;
  std::vector<DiTauObj> *myDiTaus_;
  EventObj theEvent_;

  bool isMC_;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

#endif
