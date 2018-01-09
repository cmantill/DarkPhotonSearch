import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


#process.source = cms.Source("PoolSource",
#                            fileNames = cms.untracked.vstring(
#        'file:/eos/user/a/amarini/ZeroBiasScouting_Run2017_v2/JetHT/JetHT_1/170624_222405/0000/outputScoutingCalo_9.root',
#        )
#)


import FWCore.Utilities.FileUtils as FileUtils
mylist = FileUtils.loadListFromFile ('infile5.txt') 
readFiles = cms.untracked.vstring( *mylist)

myParentlist = FileUtils.loadListFromFile ('infile5_parent.txt')
readParentFiles = cms.untracked.vstring( *myParentlist)

process.source = cms.Source('PoolSource', 
                            fileNames = readFiles,
                            secondaryFileNames = readParentFiles,
                            skipEvents = cms.untracked.uint32(310),

                            )

import FWCore.PythonUtilities.LumiList as LumiList
process.source.lumisToProcess = LumiList.LumiList(filename = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PromptReco/Cert_294927-306126_13TeV_PromptReco_Collisions17_JSON.txt').getVLuminosityBlockRange() 

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('AnalysisOutput.root'),
 #                                  closeFileFast = cms.untracked.bool(True)
)

## process.demo = cms.EDAnalyzer('DimuonScoutingAnalyzer',
##      triggerResults  = cms.InputTag("TriggerResults", "", "TEST"),
##      muons           = cms.InputTag("hltScoutingMuonPackerCalo", "", "TEST"),
## )

process.demo = cms.EDAnalyzer('HTScoutingAnalyzer',
     triggerResults  = cms.InputTag("TriggerResults", "", "HLT"),
     caloJets        = cms.InputTag("hltScoutingCaloPacker", "", "HLT"),
     pfJets          = cms.InputTag("hltScoutingPFPacker", "", "HLT"),
     pfMetPt         = cms.InputTag("hltScoutingPFPacker:pfMetPt"),
     candidates      = cms.InputTag("hltScoutingPFPacker", "", "HLT"),
     recoJets        = cms.InputTag("slimmedJets", "", "RECO"),
     metreco         = cms.InputTag("slimmedMETs"),
     recoAK8Jets     = cms.InputTag("slimmedJetsAK8", "", "RECO"),
     vtxreco         = cms.InputTag('offlineSlimmedPrimaryVertices', "", "RECO"),
     #softdropAK8Jets = cms.InputTag("AK8caPFJetsSoftDropPuppi", "", "RECO"),
     muons           = cms.InputTag("hltScoutingMuonPacker", "", "HLT"),
     recoMuons       = cms.InputTag("slimmedMuons", "", "RECO")
)


process.p = cms.Path(process.demo)
