import FWCore.ParameterSet.Config as cms

process = cms.Process("ANA")

## Switches 
runOnMC = True

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
#process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

## No of events to process
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
    #input = cms.untracked.int32(10)
)

## Input source
process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(),
    #fileNames = cms.untracked.vstring('file:./miniAOD-prod_PAT.root'),
    fileNames = cms.untracked.vstring('file:./miniAOD-prod_PAT-QQH135.root'),
)

# Global Tag
from Configuration.AlCa.GlobalTag import GlobalTag
#PLS170_V7AN1 GT for 25ns (asymptotic alignment and calibration scenario), as POSTLS170_V7 + updated JEC
#PLS170_V6AN1 GT for 50ns (more pessimistic alignment and calibration scenario), as POSTLS170_V6 + updated JEC
#GR_70_V2_AN1 GT for data in 70X, as GR_R_70_V2 + updated JEC
aGT = 'PLS170_V7AN1::All'
if not runOnMC:
    aGT = 'GR_70_V2_AN1'
process.GlobalTag = GlobalTag(process.GlobalTag, aGT, '')

## Output
process.outModule = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring( 'drop *'),
    fileName = cms.untracked.string('outfile.root'),
)
process.TFileService = cms.Service("TFileService",
    fileName = cms.string("OutTree.root")
)

## Modules and sequences
# process with gen-truth
process.printTree1 = cms.EDAnalyzer(
    "ParticleListDrawer",
    src = cms.InputTag("prunedGenParticles"),
    maxEventsToPrint  = cms.untracked.int32(1)
    )

# Jet antiPU-Id
# already embedded
#process.load("RecoJets.JetProducers.PileupJetID_cfi")
#process.pileupJetIdProducer.jets="slimmedJets"
#process.pileupJetIdProducer.jec="AK4PF"
#process.puJetIdSequence = cms.Sequence(process.pileupJetIdProducer)


# Rescale Mu/El/Tau
# Producer z LLR analysis needed?

# DiTau mass per-lepron pair (MVAMet + SVFit) + up/down
#not present for 70X - requires porting
#process.load("JetMETCorrections.METPUSubtraction.mvaPFMET_cff")

process.standardSequence = cms.Sequence()
if runOnMC:
    process.standardSequence += process.printTree1
    #reproduce genTauJets to look for genDecay mode
    process.load("PhysicsTools.JetMCAlgos.TauGenJets_cfi")
    
    process.genTauJetsMini = process.tauGenJets.clone(GenParticles =  "prunedGenParticles")
    process.standardSequence += process.genTauJetsMini
# reproduce MEt (w/o correction) does not work with packedPFCandidates which are not PFCandidates
#process.load("RecoMET.Configuration.RecoPFMET_cff")
#process.pfMetMini = process.pfMet.clone(
#    src="packedPFCandidates",
#    jets="ak5PFJetsMini")
#process.load("RecoJets.JetProducers.ak5PFJets_cfi")
#process.ak5PFJetsMini = process.ak5PFJets.clone(src = 'packedPFCandidates') 
#process.standardSequence += process.ak5PFJetsMini
#process.standardSequence += process.pfMetMini
#if runOnMC:
#    process.load("RecoMET.METProducers.genMetTrue_cfi")
#    process.genMetTrueMini = process.genMetTrue.clone(src='genParticlesForMETAllVisibleMini')
#    process.load("RecoMET.Configuration.GenMETParticles_cff")
#    process.genParticlesForMETAllVisibleMini = process.genParticlesForMETAllVisible.clone(src='packedGenParticles')
#    process.standardSequence += process.genParticlesForMETAllVisibleMini
#    process.standardSequence += process.genMetTrueMini
#process.load("PhysicsTools.PatAlgos.producersLayer1.metProducer_cfi")
#process.patMETsMini = process.patMETs.clone(
#    metSource='pfMetMini',
#    addGenMET=runOnMC,
#    genMETSource='genMetTrueMini')
#process.standardSequence += process.patMETsMini

# Analyzer
process.myMuTauAna = cms.EDAnalyzer(
    "MiniMuTauAnalyzer",
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    muons    = cms.InputTag("slimmedMuons"),
    taus     = cms.InputTag("slimmedTaus"),
    met      = cms.InputTag("slimmedMETs"),
    #met      = cms.InputTag("patMETsMini"),
    isMC     = cms.bool(runOnMC),
    prunedGenParticles =  cms.InputTag("prunedGenParticles"),
    genTauJets = cms.InputTag("genTauJetsMini")
)

## Paths
process.muTauStream = cms.Path(
    process.standardSequence +
    process.myMuTauAna
)

process.endjob = cms.EndPath(process.endOfProcess)
#process.output = cms.EndPath(process.outModule)

#process.outModule.outputCommands.append('keep *_someProduct_*_*')

# stdout/err
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.options = cms.untracked.PSet(
    #process.options, #uncomment if options definied eariler
    wantSummary = cms.untracked.bool(True)
    )
