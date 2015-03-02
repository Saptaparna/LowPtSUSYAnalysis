import FWCore.ParameterSet.Config as cms

process = cms.Process("PROD")

process.load("Configuration.StandardSequences.SimulationRandomNumberGeneratorSeeds_cff")
process.load("PhysicsTools.HepMCCandAlgos.genParticles_cfi")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.maxEvents = cms.untracked.PSet(
	    input = cms.untracked.int32(-1)
)

# The following three lines reduce the clutter of repeated printouts
# of the same exception message.
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.destinations = ['cerr']
process.MessageLogger.statistics = []
process.MessageLogger.fwkJobReports = []
process.MessageLogger.cerr.FwkReport.reportEvery = 10000


process.source = cms.Source("LHESource",
fileNames = cms.untracked.vstring('file:/uscms_data/d1/sapta/work/SUSYSearch/SignalGeneration/Pythiadecay/CMSSW_5_3_16_patch1/src/run_Stop140_pta20/PostProc_unweighted_events_Stop140_Pta20.lhe')
#fileNames = cms.untracked.vstring('file:/uscms_data/d2/sapta/work/SUSYSearch/SignalGeneration/Pythiadecay/CMSSW_5_3_16_patch1/src/run_Stop140/PostProc_unweighted_events_Stop140.lhe')
)

from Configuration.Generator.PythiaUEZ2Settings_cfi import *
#from Configuration.Generator.PythiaUED6TSettings_cfi import *


process.generator = cms.EDFilter("Pythia6HadronizerFilter",
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    maxEventsToPrint = cms.untracked.int32(1),
    pythiaPylistVerbosity = cms.untracked.int32(1),
    comEnergy = cms.double(8000.0),
    pythiaToLHE = cms.bool(True),
    PythiaParameters = cms.PSet(
        pythiaUESettingsBlock,
        processParameters = cms.vstring(
'MSEL=0','IMSS(1)=11',
'MSTP(161)=67',
'MSTP(162)=68',
'MSTP(163)=69',
),
        SLHAParameters = cms.vstring(
		'SLHAFILE = decay_PythiaCard_STOP140_WoffShell.dat'),
        parameterSets = cms.vstring('pythiaUESettings', 
            'processParameters','SLHAParameters')
    )
)

#from IOMC.RandomEngine.RandomServiceHelper import  RandomNumberServiceHelper
#randHelper =  RandomNumberServiceHelper(process.RandomNumberGeneratorService)
#randHelper.populate()
process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
    moduleSeeds = cms.PSet(
        generator = cms.untracked.uint32(123456789),
        g4SimHits = cms.untracked.uint32(123456788),
        VtxSmeared = cms.untracked.uint32(123456789)
    ),
)

process.printGenParticles = cms.EDAnalyzer("ParticleListDrawer",
     src = cms.InputTag("genParticles"),
     maxEventsToPrint = cms.untracked.int32(5)
)




process.pp = cms.Path(process.generator*process.genParticles*process.printGenParticles)

process.schedule = cms.Schedule(process.pp)

