import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

### For Puppi ##
puppiOn = 1

###



process = cms.Process("Tupelizer")
options = VarParsing.VarParsing ('standard')

options.register ('mode','mc',
          VarParsing.VarParsing.multiplicity.singleton,
          VarParsing.VarParsing.varType.string,
          "Switch between MC, data and SMS")

options.register ('hltName','HLT',
          VarParsing.VarParsing.multiplicity.singleton,
          VarParsing.VarParsing.varType.string,
          "HLT Trigger collection")

options.register ('GT','POSTLS170_V5::All',
          VarParsing.VarParsing.multiplicity.singleton,
          VarParsing.VarParsing.varType.string,
          "Global Tag")

options.register ('outfile','histo.root',
          VarParsing.VarParsing.multiplicity.singleton,
          VarParsing.VarParsing.varType.string,
          "outfile")

options.register ('triggers','*',
          VarParsing.VarParsing.multiplicity.list,
          VarParsing.VarParsing.varType.string,
          "Trigger requirement")

options.register ('triggersToMonitor','',
          VarParsing.VarParsing.multiplicity.list,
          VarParsing.VarParsing.varType.string,
          "Trigger list to monitor")

options.register ('verbose',False,
          VarParsing.VarParsing.multiplicity.singleton,
          VarParsing.VarParsing.varType.bool,
          "verbosity")

options.register ('keep','',
          VarParsing.VarParsing.multiplicity.singleton,
          VarParsing.VarParsing.varType.string,
          "additional keep statements")
#options.register ('addPDFWeights',False,
#          VarParsing.VarParsing.multiplicity.singleton,
#          VarParsing.VarParsing.varType.bool,
#          "whether or not to add pdfWeights")

#infiles = ['root://xrootd.unl.edu//store/mc/Spring14dr/WJetsToLNu_HT-100to200_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/00165B45-82E6-E311-B68D-002590AC4FEC.root']
#infiles = ['file:/data/schoef/local/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola_PU20bx25_POSTLS170_V5-v1_MINIAODSIM.root']
#infiles = ['root://eoscms.cern.ch//eos/cms/store/relval/CMSSW_7_0_5/RelValTTbar_13/GEN-SIM-RECO/POSTLS170_V6-v3/00000/0423767B-B5DD-E311-A1E0-02163E00E5B5.root']

#infiles = ['file:/afs/hephy.at/scratch/s/schoefbeck/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola_PU20bx25_POSTLS170_V5-v1_MINIAODSIM.root']
infiles = ['file:/afs/hephy.at/scratch/n/nrad/TTJets-53X-syncfile-AODSIM.root']
options.files=infiles

options.mode = 'mc'
options.maxEvents=1000

if not 'ipython' in VarParsing.sys.argv[0]:
  options.parseArguments()
else:
  print "No parsing of arguments!"

isMC = (options.mode.lower()=='sms' or options.mode.lower()=='mc')
if options.keep!='':
  toKeep = ['keep '+x for x in options.keep.split(',')]
else:
  toKeep=[]
print "mode",options.mode,"isMC?",isMC, ", verbose?",options.verbose,'keep?', toKeep, "GT",options.GT, ", triggers", options.triggers#, 'addPDFWeights?',options.addPDFWeights

#-- Message Logger ------------------------------------------------------------
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options = cms.untracked.PSet(
  SkipEvent = cms.untracked.vstring('ProductNotFound'),
  wantSummary = cms.untracked.bool(False),
  allowUnscheduled = cms.untracked.bool( True )
)
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.MessageLogger.categories.append('PATSummaryTables')
process.MessageLogger.cerr.PATSummaryTables = cms.untracked.PSet(
    limit = cms.untracked.int32(-1),
    reportEvery = cms.untracked.int32(100)
)

#-- Source information ------------------------------------------------------
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(options.files)
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents ) )
#process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange(
#  '190645:10-190645:110',
#)

process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

## Geometry and Detector Conditions (needed for a few patTuple production steps)
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = options.GT
process.load("Configuration.StandardSequences.MagneticField_cff")

process.out = cms.OutputModule("PoolOutputModule",
     #verbose = cms.untracked.bool(True),
     SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
     fileName = cms.untracked.string(options.outfile),
     outputCommands = cms.untracked.vstring() 
#     outputCommands = cms.untracked.vstring('keep *') 
)


################PUPPI########################
#process.load('CommonTools/ParticleFlow/PFBRECO_cff')
if puppiOn :
	process.load('Dummy/Puppi/Puppi_cff')
#	process.load('BaconProd/Ntupler/myJetExtras04Puppi_cff')
#	process.load('JetTools/AnalyzerToolbox/AnalyzerJetToolbox_cff')
	from RecoMET.METProducers.PFMET_cfi import pfMet 
	
	process.pfMetPuppi = pfMet.clone();
	process.pfMetPuppi.src = cms.InputTag('puppi','Puppi')
	process.pfMetPuppi.calculateSignificance = False
	
	
	from RecoJets.JetProducers.ak5PFJets_cfi        import ak5PFJets
	process.AK4PFJetsPuppi = ak5PFJets.clone(
	    src      = cms.InputTag('puppi','Puppi'),
	    rParam   = cms.double(0.4),
	    jetPtMin = cms.double(10)
	    )
#	process.AK4PFJetsPuppi.src = cms.InputTag('puppi','Puppi')
	process.puppiSequence = cms.Sequence(process.puppi*
					process.pfMetPuppi*
	                                     process.AK4PFJetsPuppi
	                                    )
	print "puppi is on:1"


#process.AK4PFJetsPuppi.src = cms.InputTag('puppi','Puppi')
'''
process.AK4jetsequencePuppi = cms.Sequence(
    AK4PFJetsPuppi
    )
'''       

#process.load("RecoTauTag/Configuration/RecoPFTauTag_cff")

#############################################

#import HLTrigger.HLTfilters.hltHighLevel_cfi as hlt
#process.hltFilter = hlt.hltHighLevel.clone(
#             HLTPaths = cms.vstring(options.triggers),
#             TriggerResultsTag = cms.InputTag("TriggerResults","",options.hltName),
#             throw = False
#         )
#
#process.scrapingVeto = cms.EDFilter("FilterOutScraping",
#                                             applyfilter = cms.untracked.bool(True),
#                                             debugOn = cms.untracked.bool(False),
#                                             numtrack = cms.untracked.uint32(10),
#                                             thresh = cms.untracked.double(0.25)
#                                             )
#
#process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
#                      vertexCollection = cms.InputTag('offlinePrimaryVertices'),
#                      minimumNDOF = cms.uint32(4) ,
#                      maxAbsZ = cms.double(24),
#                      maxd0 = cms.double(2))
#
#process.load('CommonTools/RecoAlgos/HBHENoiseFilter_cfi')
#
#
#process.goodVertices = cms.EDFilter(
#            "VertexSelector",
#            filter = cms.bool(False),
#            src = cms.InputTag("offlinePrimaryVertices"),
#            cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2")
#          )
#
#if not options.mode.lower()=='sms':
#  process.load("RecoMET.METFilters.hcalLaserEventFilter_cfi")
#  process.hcalLaserEventFilter.vetoByRunEventNumber=cms.untracked.bool(False)
#  process.hcalLaserEventFilter.vetoByHBHEOccupancy=cms.untracked.bool(True)
#  process.load('RecoMET.METFilters.eeBadScFilter_cfi')
#  process.load('RecoMET.METFilters.CSCTightHaloFilter_cfi')
#process.load('RecoMET.METFilters.EcalDeadCellTriggerPrimitiveFilter_cfi')
#process.load('RecoMET.METFilters.trackingFailureFilter_cfi')
#
process.load('Workspace.HEPHYCMSSWTools.EventCounter')
#
#process.EventCounterAfterHLT = process.EventCounter.clone()
#process.EventCounterAfterScraping = process.EventCounter.clone()
#process.EventCounterAfterPV = process.EventCounter.clone()
#process.EventCounterAfterHBHE = process.EventCounter.clone()
#process.EventCounterAfterTrackingFailure = process.EventCounter.clone()
#process.EventCounterAfterLaser = process.EventCounter.clone()
#process.EventCounterAfterCSC = process.EventCounter.clone()
#process.EventCounterAfterEEBadSC = process.EventCounter.clone()
#process.EventCounterAfterECALTP = process.EventCounter.clone()
#
process.filterSequence = cms.Sequence(
    process.EventCounter
#      *process.hltFilter *
#    process.EventCounterAfterHLT*
#      process.scrapingVeto *
#    process.EventCounterAfterScraping*
#      process.primaryVertexFilter*
#    process.EventCounterAfterPV
)
#process.filterSequence+= process.goodVertices
#if not options.mode.lower()=='sms':
#  process.filterSequence+= process.HBHENoiseFilter
#  process.filterSequence+= process.EventCounterAfterHBHE
#  process.filterSequence+= process.trackingFailureFilter
#  process.filterSequence+= process.EventCounterAfterTrackingFailure
#  process.filterSequence+= process.hcalLaserEventFilter
#  process.filterSequence+= process.EventCounterAfterLaser
#  process.filterSequence+= process.CSCTightHaloFilter
#  process.filterSequence+= process.EventCounterAfterCSC
#  process.filterSequence+= process.eeBadScFilter
#  process.filterSequence+= process.EventCounterAfterEEBadSC
#process.filterSequence+= process.EcalDeadCellTriggerPrimitiveFilter
#process.filterSequence+= process.EventCounterAfterECALTP
#
#if options.mode.lower()=='sms':
#  print "\nFilter List:", "HLT, scraping, PV, EcalTP\n"
#if options.mode.lower()=='mc':
#  print "\nFilter List:", "HLT, scraping, PV, HBHE, trackingFailureFilter, hcalLaser, CSCTightHalo, eeBadSC, EcalTP\n"
#
#if options.mode.lower()=='data':
#  process.load("EventFilter.HcalRawToDigi.hcallasereventfilter2012_cfi")
#  process.filterSequence+=process.hcallasereventfilter2012
#  process.EventCounterAfterHCALLaser2012 = process.EventCounter.clone()
#  process.filterSequence+=process.EventCounterAfterHCALLaser2012
#  process.load('RecoMET.METFilters.ecalLaserCorrFilter_cfi')
#  process.filterSequence+= process.ecalLaserCorrFilter
#  process.EventCounterafterECalLaserCorrectionFilter = process.EventCounter.clone()
#  process.filterSequence+=process.EventCounterafterECalLaserCorrectionFilter
#
#  print "\nFilter List:", "HLT, scraping, PV, HBHE, trackingFailureFilter, hcalLaser, CSCTightHalo, eeBadSC, EcalTP, hcalLaser2012, ecalLaserCorr\n"
from RecoMET.METProducers.PFMET_cfi import pfMet
process.pfMet = pfMet.clone(src = "packedPFCandidates")
process.pfMet.calculateSignificance = False # this can't be easily implemented on packed PF candidates at the moment

#from RecoMET.Configuration.GenMETParticles_cff import genParticlesForMETAllVisible
#vList = [abs(x) for x in list(genParticlesForMETAllVisible.ignoreParticleIDs)]
#print "status==1&("+("&".join(["abs(pdgId)!="+str(x) for x in vList]))+")"
#process.packedGenParticlesForGenMET = cms.EDFilter(
#    "CandViewRefSelector",
#    src = cms.InputTag("packedGenParticles"),
#    cut = cms.string("status==1&("+("&".join(["abs(pdgId)!="+str(x) for x in vList]))+")"),
#    noSeqChain = cms.bool(True), # Don't "chain" these sequence and update [src]
#)
#from RecoMET.METProducers.genMetTrue_cfi import genMetTrue
#process.genMetTrue = genMetTrue.clone(src = cms.InputTag('packedGenParticlesForGenMET'))
#process.metSequence = cms.Sequence( process.packedGenParticlesForGenMET * process.genMetTrue *process.pfMet)

process.metSequence = cms.Sequence(process.pfMet)

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.printTree = cms.EDAnalyzer("ParticleTreeDrawer",
                                   src = cms.InputTag("genParticles"),                                                                 
                                   printP4 = cms.untracked.bool(False),
                                   printPtEtaPhi = cms.untracked.bool(False),
                                   printVertex = cms.untracked.bool(False),
                                   printStatus = cms.untracked.bool(True),
                                   printIndex = cms.untracked.bool(False),
                                   status = cms.untracked.vint32( [1,2,3] )
                                   )
#process.printTree = cms.EDAnalyzer("ParticleListDrawer",
#    printVertex = cms.untracked.bool(False),
#    src = cms.InputTag("genParticles"),
#    maxEventsToPrint = cms.untracked.int32(-1)
#)

#-- Execution path ------------------------------------------------------------

#process.load('CommonTools.ParticleFlow.goodOfflinePrimaryVertices_cfi') #FIXME Added R.S.
##process.p += process.printTree
#if options.mode.lower()=='sms' or options.addPDFWeights:
#  process.pdfWeights = cms.EDProducer("PdfWeightProducer",
#              PdfInfoTag = cms.untracked.InputTag("generator"),
#              PdfSetNames = cms.untracked.vstring(
#    "cteq66.LHgrid"
#    , "MSTW2008nlo68cl.LHgrid"
#    , "NNPDF20_100.LHgrid"
#    ))
#  process.p += process.pdfWeights
#
#process.load("PhysicsTools.HepMCCandAlgos.flavorHistoryProducer_cfi")
#process.load("PhysicsTools.HepMCCandAlgos.flavorHistoryFilter_cfi")
#if isMC:
#  process.p +=      process.bFlavorHistoryProducer
#  process.p +=      process.cFlavorHistoryProducer
#  process.p +=      process.flavorHistoryFilter

#if options.triggersToMonitor!='':
#  options.triggersToMonitor+=options.triggers
#else:
#  options.triggersToMonitor = options.triggers
#for t in options.triggersToMonitor:
#  ts = t.replace("_v*","")
#  if ts != "*":
#    process.SUSYTupelizer.triggersToMonitor.append(ts)
#process.SUSYTupelizer.triggersToMonitor = list(set(process.SUSYTupelizer.triggersToMonitor)) #remove duplicates
#print "TriggersToMonitor:",process.SUSYTupelizer.triggersToMonitor
#from RecoJets.JetProducers.kt4PFJets_cfi import *
#process.kt6PFJetsForIsolation2011 = kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
#process.kt6PFJetsForIsolation2011.Rho_EtaMax = cms.double(2.5)
#process.p+=process.kt6PFJetsForIsolation2011

process.miniAODTupelizerSequence = cms.Sequence()
process.load('Workspace.HEPHYCMSSWTools.JetTupelizer_miniAOD_cfi')
process.JetTupelizer.useForDefaultAlias = cms.untracked.bool(True)
process.miniAODTupelizerSequence += process.JetTupelizer
process.load('Workspace.HEPHYCMSSWTools.MuonTupelizer_miniAOD_cfi')
process.MuonTupelizer.useForDefaultAlias = cms.untracked.bool(True)
process.miniAODTupelizerSequence += process.MuonTupelizer
process.load('Workspace.HEPHYCMSSWTools.ElectronTupelizer_miniAOD_cfi')
process.ElectronTupelizer.useForDefaultAlias = cms.untracked.bool(True)
process.miniAODTupelizerSequence += process.ElectronTupelizer
process.load('Workspace.HEPHYCMSSWTools.TauTupelizer_miniAOD_cfi')
process.TauTupelizer.useForDefaultAlias = cms.untracked.bool(True)
process.miniAODTupelizerSequence += process.TauTupelizer
process.load('Workspace.HEPHYCMSSWTools.TriggerTupelizer_cfi')
process.TriggerTupelizer.useForDefaultAlias = cms.untracked.bool(True)
process.miniAODTupelizerSequence += process.TriggerTupelizer
process.load('Workspace.HEPHYCMSSWTools.FilterTupelizer_cfi')
process.FilterTupelizer.useForDefaultAlias = cms.untracked.bool(True)
process.miniAODTupelizerSequence += process.FilterTupelizer
process.load('Workspace.HEPHYCMSSWTools.BasicTupelizer_miniAOD_cfi')
process.BasicTupelizer.useForDefaultAlias = cms.untracked.bool(True)
process.BasicTupelizer.addMSugraOSETInfo = cms.untracked.bool(options.mode.lower()=='sms')
process.miniAODTupelizerSequence += process.BasicTupelizer

#process.SUSYTupelizer.triggerCollection = cms.untracked.InputTag( options.hltName )
#process.SUSYTupelizer.addTriggerInfo = cms.untracked.bool(True)
#process.SUSYTupelizer.addMSugraOSETInfo = cms.untracked.bool(options.mode.lower()=='sms')
#process.SUSYTupelizer.addPDFWeights = cms.untracked.bool(options.mode.lower()=='sms' or options.addPDFWeights)
#process.SUSYTupelizer.verbose = cms.untracked.bool(options.verbose)
#process.SUSYTupelizer.addMuonVector = cms.untracked.bool(True)
#process.SUSYTupelizer.addEleVector = cms.untracked.bool(True)
#process.SUSYTupelizer.addTauiVector = cms.untracked.bool(True)
#process.SUSYTupelizer.metsToMonitor = []

if puppiOn :
	process.p = cms.Path( process.puppiSequence + process.filterSequence + process.metSequence + process.miniAODTupelizerSequence)
else : 
	process.p = cms.Path(process.filterSequence + process.metSequence + process.miniAODTupelizerSequence)
#process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

#process.printTree = cms.EDAnalyzer("ParticleListDrawer",
#  maxEventsToPrint = cms.untracked.int32(-1),
#  printVertex = cms.untracked.bool(False),
#  src = cms.InputTag("genParticles")
#)
#process.p+=process.printTree

process.out.outputCommands =  cms.untracked.vstring('drop *', 'keep *_*Tupelizer*_*_*' , 'keep *_*EventCounter*_*_*', 'keep *_*uppi*_*_*','keep *_*_*uppi*_*' ,*(toKeep))

process.outpath = cms.EndPath(process.out)
