# Info at https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD
# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: miniAOD-prod -s PAT --eventcontent MINIAODSIM --runUnscheduled --mc --filein /store/relval/CMSSW_7_0_0/RelValTTbar_13/GEN-SIM-RECO/PU25ns_POSTLS170_V5_AlcaCSA14-v1/00000/1A0A1B22-D2AA-E311-8197-02163E00E8DA.root --conditions PLS170_V7AN1::All --no_exec
import FWCore.ParameterSet.Config as cms

process = cms.Process('PAT')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.PATMC_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(),
    fileNames = cms.untracked.vstring(
        #'/store/relval/CMSSW_7_0_0/RelValZTT_13/GEN-SIM-RECO/PU25ns_POSTLS170_V5_AlcaCSA14-v1/00000/08E8549F-0DAB-E311-B423-02163E00E6D6.root',
        '/store/relval/CMSSW_7_0_0/RelValQQH1352T_Tauola_13/GEN-SIM-RECO/PU25ns_POSTLS170_V5_AlcaCSA14-v1/00000/06307AD4-CCAA-E311-9F24-02163E00E99E.root',
        '/store/relval/CMSSW_7_0_0/RelValQQH1352T_Tauola_13/GEN-SIM-RECO/PU25ns_POSTLS170_V5_AlcaCSA14-v1/00000/165EE938-DEAA-E311-9C88-02163E00A0DB.root',
        '/store/relval/CMSSW_7_0_0/RelValQQH1352T_Tauola_13/GEN-SIM-RECO/PU25ns_POSTLS170_V5_AlcaCSA14-v1/00000/602EF95D-D0AA-E311-B331-002590495240.root',
        '/store/relval/CMSSW_7_0_0/RelValQQH1352T_Tauola_13/GEN-SIM-RECO/PU25ns_POSTLS170_V5_AlcaCSA14-v1/00000/62F3C1C7-D6AA-E311-8297-02163E00E6DE.root',
        '/store/relval/CMSSW_7_0_0/RelValQQH1352T_Tauola_13/GEN-SIM-RECO/PU25ns_POSTLS170_V5_AlcaCSA14-v1/00000/645D007C-EBAA-E311-AB4C-02163E00EA73.root',
    )
)

process.options = cms.untracked.PSet(
    allowUnscheduled = cms.untracked.bool(True)
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.19 $'),
    annotation = cms.untracked.string('miniAOD-prod nevts:1'),
    name = cms.untracked.string('Applications')
)

# Output definition

process.MINIAODSIMoutput = cms.OutputModule("PoolOutputModule",
    compressionLevel = cms.untracked.int32(4),
    compressionAlgorithm = cms.untracked.string('LZMA'),
    eventAutoFlushCompressedSize = cms.untracked.int32(15728640),
    outputCommands = process.MINIAODSIMEventContent.outputCommands,
    fileName = cms.untracked.string('miniAOD-prod_PAT-QQH135.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('')
    ),
    dropMetaData = cms.untracked.string('ALL'),
    fastCloning = cms.untracked.bool(False),
    overrideInputFileSplitLevels = cms.untracked.bool(True)
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
#PLS170_V7AN1 GT for 25ns (asymptotic alignment and calibration scenario), as POSTLS170_V7 + updated JEC
#PLS170_V6AN1 GT for 50ns (more pessimistic alignment and calibration scenario), as POSTLS170_V6 + updated JEC
#GR_70_V2_AN1 GT for data in 70X, as GR_R_70_V2 + updated JEC
process.GlobalTag = GlobalTag(process.GlobalTag, 'PLS170_V7AN1::All', '')

# Path and EndPath definitions
process.Flag_trackingFailureFilter = cms.Path(process.goodVertices+process.trackingFailureFilter)
process.Flag_goodVertices = cms.Path(process.goodVertices)
process.Flag_CSCTightHaloFilter = cms.Path(process.CSCTightHaloFilter)
process.Flag_trkPOGFilters = cms.Path(process.trkPOGFilters)
process.Flag_trkPOG_logErrorTooManyClusters = cms.Path(~process.logErrorTooManyClusters)
process.Flag_EcalDeadCellTriggerPrimitiveFilter = cms.Path(process.EcalDeadCellTriggerPrimitiveFilter)
process.Flag_ecalLaserCorrFilter = cms.Path(process.ecalLaserCorrFilter)
process.Flag_trkPOG_manystripclus53X = cms.Path(~process.manystripclus53X)
process.Flag_eeBadScFilter = cms.Path(process.eeBadScFilter)
process.Flag_METFilters = cms.Path(process.metFilters)
process.Flag_HBHENoiseFilter = cms.Path(process.HBHENoiseFilter)
process.Flag_trkPOG_toomanystripclus53X = cms.Path(~process.toomanystripclus53X)
process.Flag_hcalLaserEventFilter = cms.Path(process.hcalLaserEventFilter)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.MINIAODSIMoutput_step = cms.EndPath(process.MINIAODSIMoutput)

# customisation of the process.

# Automatic addition of the customisation function from PhysicsTools.PatAlgos.slimming.miniAOD_tools
from PhysicsTools.PatAlgos.slimming.miniAOD_tools import miniAOD_customizeAllMC 

#call to customisation function miniAOD_customizeAllMC imported from PhysicsTools.PatAlgos.slimming.miniAOD_tools
process = miniAOD_customizeAllMC(process)

# End of customisation functions

## MB:
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.options = cms.untracked.PSet(
    process.options,
    wantSummary = cms.untracked.bool(True)
    )
