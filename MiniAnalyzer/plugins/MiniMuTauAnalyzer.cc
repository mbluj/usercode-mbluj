
#include "MiniMuTauAnalyzer.h"

// user include files
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

//Standalone SVFit///
#include "TauAnalysis/SVfitStandalone/interface/SVfitStandaloneAlgorithm.h"

#include "PhysicsTools/JetMCUtils/interface/JetMCTag.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"

#include "CommonTools/UtilAlgos/interface/DeltaR.h"

//
// constructors and destructor
//
MiniMuTauAnalyzer::MiniMuTauAnalyzer(const edm::ParameterSet& iConfig){

  //now do what ever initialization is needed
  vtxToken_ = consumes<reco::VertexCollection>( iConfig.getParameter<edm::InputTag>("vertices") );
  muonToken_ = consumes<pat::MuonCollection>( iConfig.getParameter<edm::InputTag>("muons") );
  tauToken_ = consumes<pat::TauCollection>( iConfig.getParameter<edm::InputTag>("taus") );
  metToken_ = consumes<pat::METCollection>( iConfig.getParameter<edm::InputTag>("met") );

  isMC_ = iConfig.getParameter<bool>("isMC");
  if(isMC_){
    genTauJetsToken_ = consumes<edm::View<reco::Candidate> >( iConfig.getParameter<edm::InputTag>("genTauJets") );
    prunedGenParticlesToken_ = consumes<reco::GenParticleCollection>( iConfig.getParameter<edm::InputTag>("prunedGenParticles") );
  }

  myMuons_  = new std::vector<MuonObj>;
  myTaus_   = new std::vector<TauObj>;
  myJets_   = new std::vector<JetObj>;
  myMets_   = new std::vector<MetObj>;
  myDiTaus_ = new std::vector<DiTauObj>;

  theEvent_ = EventObj();

}


MiniMuTauAnalyzer::~MiniMuTauAnalyzer(){
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  delete myMuons_;
  delete myTaus_;
  delete myJets_;
  delete myMets_;
  delete myDiTaus_;

}


//
// member functions
//

// ------------ method called for each event  ------------
void MiniMuTauAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  using namespace edm;

  //clear (move to special method)
  myMuons_->clear();
  myTaus_->clear();
  myJets_->clear();
  myMets_->clear();
  myDiTaus_->clear();
  
  //basic event info
  theEvent_.setRunInfo(iEvent.run(), iEvent.luminosityBlock(), (iEvent.eventAuxiliary()).event(), isMC_, 1);

#ifdef THIS_IS_AN_EVENT_EXAMPLE
  Handle<ExampleData> pIn;
  iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  ESHandle<SetupData> pSetup;
  iSetup.get<SetupRecord>().get(pSetup);
#endif

  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vtxToken_, vertices);
  if( vertices->empty() ) return; // skip the event if no PV found
  const reco::Vertex &thePV = vertices->front();

  edm::Handle<pat::MuonCollection> muons;
  iEvent.getByToken(muonToken_, muons);

  for(const pat::Muon &mu : *muons) {
    if( mu.pt() < 15 || !muonId(mu,true) ) continue;
    double iso = ( mu.userIsolation(pat::PfChargedAllIso) +
		   std::max(0.0,
			    mu.userIsolation(pat::PfGammaIso) +
			    mu.userIsolation(pat::PfNeutralHadronIso) -
			    0.5*mu.userIsolation(pat::PfPUChargedHadronIso) ) );
    if(iso/mu.pt()>0.5) continue;
    /*
    printf("muon with pt %4.1f, dz(PV) %+5.3f, POG loose id %d, tight id %d\n",
	   mu.pt(), mu.muonBestTrack()->dz( thePV.position() ), mu.isLooseMuon(), mu.isTightMuon(thePV) );
    */
    //MuonObj aMuon( mu.pt(),mu.eta(),mu.phi(),mu.mass(), mu.charge() );
    MuonObj aMuon( mu.pt(),mu.eta(),mu.phi(), 0.106, mu.charge() );
    //here do the rest
    aMuon.setBits(mu.isGlobalMuon(), mu.isTrackerMuon(), mu.isStandAloneMuon(),
		  mu.isCaloMuon(), false, 
		  mu.isPFMuon(), mu.isLooseMuon(), mu.isTightMuon(thePV),
		  muonId(mu,false) );
    aMuon.nAllMuons = muons->size();
    aMuon.chHadIso = mu.userIsolation(pat::PfChargedHadronIso);
    aMuon.chIso = mu.userIsolation(pat::PfChargedAllIso);
    aMuon.nHadIso = mu.userIsolation(pat::PfNeutralHadronIso);
    aMuon.phIso = mu.userIsolation(pat::PfGammaIso);
    aMuon.puIso = mu.userIsolation(pat::PfPUChargedHadronIso);
    aMuon.dz = mu.innerTrack()->dz( thePV.position() );
    aMuon.dzErr = mu.innerTrack()->dzError();
    aMuon.dxy = mu.innerTrack()->dxy( thePV.position() );
    aMuon.dxyErr = mu.innerTrack()->dxyError();
    if(isMC_ && mu.genLepton()!=0 ){
      aMuon.genPt = mu.genLepton()->pt();
      aMuon.genPdgId = mu.genLepton()->pdgId();
    } else {
      aMuon.genPt = -1;
      aMuon.genPdgId = 0;						
    }
    //add to collection
    myMuons_->push_back(aMuon);
  }

  edm::Handle<pat::TauCollection> taus;
  iEvent.getByToken(tauToken_, taus);

  edm::Handle<edm::View<reco::Candidate> > genTauJets;
  if(isMC_){
    iEvent.getByToken(genTauJetsToken_, genTauJets);
  }

  for(const pat::Tau &tau : *taus) {
    if( tau.pt() < 20 || fabs(tau.eta() ) > 2.3 
	|| ( tau.tauID("decayModeFinding") < 0.5 &&
	     tau.tauID("decayModeFindingNewDMs") < 0.5 )
	|| ( tau.tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits") < 0.5 &&
	     tau.tauID("byVLooseIsolationMVA3oldDMwLT") < 0.5 )
	|| tau.tauID("againstMuonLoose3") < 0.5
	|| tau.tauID("againstElectronLoose") < 0.5
	) continue;
    /*
    printf("tau with pt %4.1f, dz(PV) %+5.3f, \n newDM %d, oldDM %d, \n iso3hits %f, looseMVAnewDM %d, looseMVAoldDM %d\n",
	   tau.pt(), 
	   //fabs(tau.vertex().z()-thePV.position().z() ), //MB always 0, is it OK?
	   fabs(tau.leadChargedHadrCand()->vertex().z()-thePV.position().z() ),
	   tau.tauID("decayModeFindingNewDMs")>0.5,
	   tau.tauID("decayModeFinding")>0.5,
	   tau.tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits"),
	   tau.tauID("byLooseIsolationMVA3newDMwLT")>0.5,
	   tau.tauID("byLooseIsolationMVA3oldDMwLT")>0.5);
    */
    TauObj aTau( tau.pt(),tau.eta(),tau.phi(),tau.mass(), tau.charge() );
    //here do the rest                                                         
    aTau.nAllTaus = taus->size();    
    aTau.dzPV = tau.leadChargedHadrCand()->vertex().z()-thePV.position().z();
    aTau.hasKft = -1;
    if(tau.leadChargedHadrCand().isNonnull() ){
      const reco::PFCandidate* leadPfCand = dynamic_cast<const reco::PFCandidate*>( (tau.leadChargedHadrCand()).get() ); 
      if(leadPfCand != 0) {
	if( (leadPfCand->trackRef()).isNonnull() ){
	  aTau.dxyTrk     = leadPfCand->trackRef()->dxy( thePV.position() );
	  aTau.dzTrk      = leadPfCand->trackRef()->dz( thePV.position() );
	  aTau.dxyTrkErr  = leadPfCand->trackRef()->dxyError();
	  aTau.dzTrkErr   = leadPfCand->trackRef()->dzError();
	  aTau.hasKft     = 1;
	}
	if( (leadPfCand->gsfTrackRef()).isNonnull() ){
	  aTau.dxyTrk    = leadPfCand->gsfTrackRef()->dxy( thePV.position() );
	  aTau.dzTrk     = leadPfCand->gsfTrackRef()->dz( thePV.position() );
	  aTau.dxyTrkErr = leadPfCand->gsfTrackRef()->dxyError();
	  aTau.dzTrkErr  = leadPfCand->gsfTrackRef()->dzError();
	  aTau.hasKft    = 0;
	}
      } else {
	const pat::PackedCandidate *leadCand = dynamic_cast<const pat::PackedCandidate*>( (tau.leadChargedHadrCand()).get() );
	if( leadCand !=0 ) {
	  aTau.dxyTrk     = leadCand->dxy( thePV.position() );
	  aTau.dzTrk      = leadCand->dz( thePV.position() );
	  aTau.dxyTrkErr  = leadCand->dxyError();
	  aTau.dzTrkErr   = leadCand->dzError();
	  aTau.hasKft     = -1;
	}
      }
    }
    aTau.decMode = findDecayMode(tau);
    aTau.iso = tau.tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
    aTau.chIso = tau.tauID("chargedIsoPtSum");
    aTau.phIso = tau.tauID("neutralIsoPtSum");
    aTau.puIso = tau.tauID("puCorrPtSum");

    int tightestHPSDB3HWP = 0;
    if(tau.tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits")>0.5)    tightestHPSDB3HWP=1; 
    if(tau.tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits")>0.5) tightestHPSDB3HWP=2; 
    if(tau.tauID("byTightCombinedIsolationDeltaBetaCorr3Hits")>0.5)  tightestHPSDB3HWP=3;
    int tightestHPSMVA3oldDMwoLTWP = 0;
    if(tau.tauID("byVLooseIsolationMVA3oldDMwoLT") >0.5) tightestHPSMVA3oldDMwoLTWP=1;
    if(tau.tauID("byLooseIsolationMVA3oldDMwoLT")>0.5)   tightestHPSMVA3oldDMwoLTWP=2;
    if(tau.tauID("byMediumIsolationMVA3oldDMwoLT") >0.5) tightestHPSMVA3oldDMwoLTWP=3;
    if(tau.tauID("byTightIsolationMVA3oldDMwoLT") >0.5)  tightestHPSMVA3oldDMwoLTWP=4;
    if(tau.tauID("byVTightIsolationMVA3oldDMwoLT") >0.5) tightestHPSMVA3oldDMwoLTWP=5;
    if(tau.tauID("byVVTightIsolationMVA3oldDMwoLT") >0.5)tightestHPSMVA3oldDMwoLTWP=6;
    int tightestHPSMVA3oldDMwLTWP = 0;
    if(tau.tauID("byVLooseIsolationMVA3oldDMwLT") >0.5) tightestHPSMVA3oldDMwLTWP=1;
    if(tau.tauID("byLooseIsolationMVA3oldDMwLT")>0.5)   tightestHPSMVA3oldDMwLTWP=2;
    if(tau.tauID("byMediumIsolationMVA3oldDMwLT") >0.5) tightestHPSMVA3oldDMwLTWP=3;
    if(tau.tauID("byTightIsolationMVA3oldDMwLT") >0.5)  tightestHPSMVA3oldDMwLTWP=4;
    if(tau.tauID("byVTightIsolationMVA3oldDMwLT") >0.5) tightestHPSMVA3oldDMwLTWP=5;
    if(tau.tauID("byVVTightIsolationMVA3oldDMwLT") >0.5)tightestHPSMVA3oldDMwLTWP=6;
    aTau.tauIso = (tightestHPSMVA3oldDMwLTWP*100 +
		   tightestHPSMVA3oldDMwoLTWP*10 +
		   tightestHPSDB3HWP);

    int tightestAntiMu3WP = 0;   
    if( tau.tauID("againstMuonLoose3")>0.5 )tightestAntiMu3WP = 1;   
    if( tau.tauID("againstMuonTight3")>0.5 )tightestAntiMu3WP = 2;   
    int tightestAntiMuMVAWP = 0;   
    if( tau.tauID("againstMuonLooseMVA")>0.5 )tightestAntiMuMVAWP = 1;   
    if( tau.tauID("againstMuonMediumMVA")>0.5 )tightestAntiMuMVAWP = 2;   
    if( tau.tauID("againstMuonTightMVA")>0.5 )tightestAntiMuMVAWP = 3;   
    aTau.antiMu = (tightestAntiMuMVAWP*10 +
		   tightestAntiMu3WP);    

    int tightestAntiECutWP = 0;
    if( tau.tauID("againstElectronLoose")>0.5 )tightestAntiECutWP = 1;
    if( tau.tauID("againstElectronMedium")>0.5 )tightestAntiECutWP = 2;
    if( tau.tauID("againstElectronTight")>0.5 )tightestAntiECutWP = 3;
    int tightestAntiEMVA5WP = 0;
    if( tau.tauID("againstElectronVLooseMVA5")>0.5) tightestAntiEMVA5WP  = 1;
    if( tau.tauID("againstElectronLooseMVA5")>0.5)  tightestAntiEMVA5WP  = 2;
    if( tau.tauID("againstElectronMediumMVA5")>0.5) tightestAntiEMVA5WP  = 3;
    if( tau.tauID("againstElectronTightMVA5")>0.5)  tightestAntiEMVA5WP  = 4;
    if( tau.tauID("againstElectronVTightMVA5")>0.5) tightestAntiEMVA5WP  = 5;
    aTau.antiE = (tightestAntiEMVA5WP*10 +
		  tightestAntiECutWP);

    /////NewTauID input variables
    aTau.dxyTau            = tau.dxy(); 
    aTau.dxyTauErr         = tau.dxy_error(); 
    aTau.dxyTauSig         = tau.dxy_Sig();
    aTau.hasSecVtx         = tau.hasSecondaryVertex();
    aTau.flightLength      = sqrt(tau.flightLength().Mag2());
    aTau.flightLengthError = tau.flightLengthSig()/sqrt(tau.flightLength().Mag2());
    aTau.flightLengthSig   = tau.flightLengthSig();

    if(isMC_){
      //to be computed
      aTau.genPdgId = 0;
      aTau.nGenPart = 0;
      aTau.leadGenPt = 0;

      int nMached=-1;
      int idx = match(&tau,genTauJets.product(),0.15,nMached,3.0);
      if(idx>-1){
	aTau.genTauJetPt = (*genTauJets)[idx].pt();
	aTau.genVisMass = (*genTauJets)[idx].mass();	
	const reco::GenJet *aGenTauJet = dynamic_cast<const reco::GenJet*>(&(*genTauJets)[idx]);
	aTau.genDecMode = findGenDecayMode(*aGenTauJet);
      }
      /*
      if(tau.genJet()!=0 && idx>-1){
	if( fabs(tau.genJet()->pt()/aTau.genTauJetPt-1.)>1.1)
	  std::cout<<"Different match:"
		   <<" genJet(PAT).pt="<<tau.genJet()->pt()
		   <<" genJet(custom).pt="<<aTau.genTauJetPt
		   <<std::endl;
      }
      else{
	std::cout<<"Different match:"
		 <<" *genJet(PAT)="<<tau.genJet()
		 <<" genJet(custom).idx="<<idx
	         <<" decMode: "<< (idx>-1 ? aTau.genDecMode : -1)
		 <<std::endl;
      }
      */
    }

    myTaus_->push_back(aTau);
  }

  edm::Handle<pat::METCollection> mets;
  iEvent.getByToken(metToken_, mets);
  for(const pat::MET &met : *mets) {
    MetObj aMet(met.pt(),met.phi());
    aMet.nAllMets = mets->size();
    aMet.sumEt = met.sumEt();
    aMet.sig = met.significance();
    const TMatrixD cov = met.getSignificanceMatrix();
    const double* elements;
    aMet.setSigMatrix(elements[0],elements[3],elements[1],elements[2]);
    if(isMC_ && met.genMET()!=0 ){
      aMet.genMetPt = met.genMET()->pt();
      aMet.genMetPhi = met.genMET()->phi();
    }
    myMets_->push_back(aMet);
  }

  //Build di-taus
  for(unsigned int iMu=0; iMu<myMuons_->size(); ++iMu){
    for(unsigned int iTau=0; iTau<myTaus_->size(); ++iTau){
      if(deltaR((*myMuons_)[iMu],(*myTaus_)[iTau]) < 0.5 ) continue;
      reco::Candidate::PolarLorentzVector muP4, tauP4, diTauP4;
      muP4.SetCoordinates( (*myMuons_)[iMu].pt(), 
			   (*myMuons_)[iMu].eta(),
			   (*myMuons_)[iMu].phi(),
			   (*myMuons_)[iMu].mass() );
      tauP4.SetCoordinates( (*myTaus_)[iTau].pt(), 
			    (*myTaus_)[iTau].eta(),
			    (*myTaus_)[iTau].phi(),
			    (*myTaus_)[iTau].mass() );
      diTauP4 = muP4 + tauP4; 
      DiTauObj aDiTau(diTauP4.pt(),
		      diTauP4.eta(),
		      diTauP4.phi(),
		      diTauP4.mass(),
		      (*myMuons_)[iMu].charge()+(*myTaus_)[iTau].charge() );
      aDiTau.leg1Idx = iMu;
      aDiTau.leg1PdgId = -13*(*myMuons_)[iMu].charge();
      aDiTau.leg2Idx = iTau;
      aDiTau.leg2PdgId = -15*(*myTaus_)[iTau].charge();
      aDiTau.metIdx = 0;
      const pat::MET &met = (*mets)[aDiTau.metIdx];
      const TMatrixD cov = met.getSignificanceMatrix();

      /*
      const double* elements;
      elements = cov.GetMatrixArray();
      std::cout<<" met cov modified "
	       <<" [0] "<<elements[0]
	       <<" [1] "<<elements[1]
	       <<" [2] "<<elements[2]
	       <<" [3] "<<elements[3]
	       <<std::endl;
      */
      /*
      metSgnMatrix_->push_back( elements[0] ); //sigma_E
      metSgnMatrix_->push_back( elements[1] ); //sigma_Ephi
      metSgnMatrix_->push_back( elements[3] ); //sigma_phi 
      */
      //Run Standalone SVFit  
      std::vector<svFitStandalone::MeasuredTauLepton> measuredTauLeptons; 
      svFitStandalone::Vector measuredMET( met.p4().Px(), met.p4().Py(), 0); 
      reco::Candidate::LorentzVector muP4_xyzt, tauP4_xyzt;
      muP4_xyzt.SetPxPyPzE(muP4.Px(),muP4.Py(),muP4.Pz(),muP4.E());
      tauP4_xyzt.SetPxPyPzE(tauP4.Px(),tauP4.Py(),tauP4.Pz(),tauP4.E());
      measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(svFitStandalone::kTauToMuDecay, muP4_xyzt)); //kTauToElecDecay
      measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(svFitStandalone::kTauToHadDecay, tauP4_xyzt)); 
      SVfitStandaloneAlgorithm algo(measuredTauLeptons, measuredMET, cov, 0); 
      algo.addLogM(false); 
      //algo.integrateMarkovChain(); 
      algo.integrateVEGAS(); 
      aDiTau.svFitMass = algo.getMass(); 
      aDiTau.svFitMassErrUp = algo.massUncert(); 
      aDiTau.svFitMassErrDown = aDiTau.svFitMassErrUp;
      aDiTau.svFitPt = algo.pt();  
      aDiTau.svFitPtErrUp = algo.ptUncert();
      aDiTau.svFitPtErrDown = aDiTau.svFitPtErrUp;

      myDiTaus_->push_back(aDiTau);
    }
  }
  /*
  std::cout<<"No. of muons / taus / diTaus : "
	   <<myMuons_->size()<<" / "
	   <<myTaus_->size()<<" / "
	   <<myDiTaus_->size()<<" "
	   <<std::endl;
  */
  tree_->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void MiniMuTauAnalyzer::beginJob() {
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("tree","muTau tree");

  tree_->Branch("event","EventObj",&theEvent_);
  tree_->Branch("muons","std::vector<MuonObj>",&myMuons_);
  tree_->Branch("taus","std::vector<TauObj>",&myTaus_);
  tree_->Branch("jets","std::vector<JetObj>",&myJets_);
  tree_->Branch("mets","std::vector<MetObj>",&myMets_);
  tree_->Branch("diTaus","std::vector<DiTauObj>",&myDiTaus_);

}

// ------------ method called once each job just after ending the event loop  ------------
void MiniMuTauAnalyzer::endJob() {
}

// ------------ method called when starting to processes a run  ------------
/*
void MiniMuTauAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&) {
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void MiniMuTauAnalyzer::endRun(edm::Run const&, edm::EventSetup const&) {
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void MiniMuTauAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void MiniMuTauAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void MiniMuTauAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

bool MiniMuTauAnalyzer::muonId(const pat::Muon& aMuon, bool checkLoose) {

  if(checkLoose)
    return aMuon.isGlobalMuon();
  else
    return (aMuon.isGlobalMuon() &&
	    aMuon.globalTrack().isNonnull() &&
	    aMuon.globalTrack()->normalizedChi2()<10 &&
	    (aMuon.globalTrack()->hitPattern()).numberOfValidMuonHits()>0 &&
	    aMuon.numberOfMatchedStations()>1 &&
	    (aMuon.innerTrack()->hitPattern()).numberOfValidPixelHits()>0 &&
	    (aMuon.track()->hitPattern()).trackerLayersWithMeasurement()>5);
  
  return false;
}

int MiniMuTauAnalyzer::findGenDecayMode(const reco::GenJet& genTauJet){

  int genDecayMode = -99;
  std::string genTauDecay = JetMCTagUtils::genTauDecayMode(genTauJet);
  if( genTauDecay.find("oneProng0Pi0")!=std::string::npos ) 
    genDecayMode = 0;
  else if( genTauDecay.find("oneProng1Pi0")!=std::string::npos )
    genDecayMode = 1;
  else if( genTauDecay.find("oneProng2Pi0")!=std::string::npos )
    genDecayMode = 2;
  else if( genTauDecay.find("oneProngOther")!=std::string::npos )
    genDecayMode = 3;
  if( genTauDecay.find("twoProng0Pi0")!=std::string::npos ) 
    genDecayMode = 4;
  else if( genTauDecay.find("twoProng1Pi0")!=std::string::npos )
    genDecayMode = 5;
  else if( genTauDecay.find("twoProng2Pi0")!=std::string::npos )
    genDecayMode = 6;
  else if( genTauDecay.find("twoProngOther")!=std::string::npos )
    genDecayMode = 7;
  else if( genTauDecay.find("threeProng0Pi0")!=std::string::npos )
    genDecayMode = 8;
  else if( genTauDecay.find("threeProng1Pi0")!=std::string::npos )
    genDecayMode = 9;
  else if( genTauDecay.find("threeProngOther")!=std::string::npos )
    genDecayMode = 10;
  else if( genTauDecay.find("rare")!=std::string::npos )
    genDecayMode = 11;
  else if( genTauDecay.find("electron")!=std::string::npos )
    genDecayMode = -11;
  else if( genTauDecay.find("muon")!=std::string::npos )
    genDecayMode = -13;
  else
    genDecayMode = -99;

  return genDecayMode;
}

int MiniMuTauAnalyzer::findDecayMode(const pat::Tau& aTau){

  int decayMode = -99;
  if((aTau.signalChargedHadrCands()).size()==1 && (aTau.signalGammaCands()).size()==0) decayMode = 0; 
  else if((aTau.signalChargedHadrCands()).size()==1 && (aTau.signalGammaCands()).size()>0)  decayMode = 1; 
  else if((aTau.signalChargedHadrCands()).size()==2 && (aTau.signalGammaCands()).size()==0)  decayMode = 2; 
  else if((aTau.signalChargedHadrCands()).size()==2 && (aTau.signalGammaCands()).size()>0)  decayMode = 3; 
  else if((aTau.signalChargedHadrCands()).size()==3) decayMode = 4; 
  else  decayMode = -99;

  return decayMode;
}

int MiniMuTauAnalyzer::match(const reco::Candidate *cand, const edm::View<reco::Candidate>  *collection,
			     float dRmax, int &nMached, 
			     float dPtmax, float Ptmin, bool byPt){

  int idx = -1;
  nMached = 0;
  float dR = dRmax;
  float pTmax = Ptmin;
  for(unsigned int i=0; i<collection->size(); ++i){
    if((*collection)[i].pt()<Ptmin) continue;
    if(dPtmax>0 && fabs( (*collection)[i].pt()/cand->pt()-1.)>dPtmax) continue;
    float dRtmp = reco::deltaR(cand->p4(), (*collection)[i].p4() );
    if(dRtmp<dRmax){
      nMached++;
      if(!byPt && dRtmp<dR){
	idx = i;
	dR = dRtmp;
      }
      if(byPt && (*collection)[i].pt()>pTmax){
	idx = i;
	pTmax = (*collection)[i].pt();
      }
    }
  }

  return idx;
}

//define this as a plug-in
DEFINE_FWK_MODULE(MiniMuTauAnalyzer);
