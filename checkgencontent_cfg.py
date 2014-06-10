import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

 
process.printTree = cms.EDAnalyzer("ParticleTreeDrawer",
                                   src = cms.InputTag("genParticles"),      
                                   printP4 = cms.untracked.bool(False),
                                   printPtEtaPhi = cms.untracked.bool(False),
                                   printVertex = cms.untracked.bool(False),
                                   printStatus = cms.untracked.bool(False),
                                   printIndex = cms.untracked.bool(False),
                                   status = cms.untracked.vint32(3)
                                   )


"""
process.printTree = cms.EDAnalyzer("ParticleListDrawer",
                                    maxEventsToPrint = cms.untracked.int32(2),
                                    printVertex = cms.untracked.bool(False),
                                    src = cms.InputTag("genParticles")
)

process.printDecay = cms.EDAnalyzer("ParticleDecayDrawer",
                                   src = cms.InputTag("genParticles"),
                                   printP4 = cms.untracked.bool(False),
                                   printPtEtaPhi = cms.untracked.bool(False),
                                   printVertex = cms.untracked.bool(False)
                                   )
"""

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/uscms_data/d2/sapta/work/SUSYSearch/CMSSW_5_3_11_patch6/src/CheckGenContent/CheckGenContent/test/006B5ECF-8124-E211-9AF8-002618943982.root'#ZGammaInclusive
        #'file:/uscms_data/d2/sapta/work/SUSYSearch/CMSSW_5_3_11_patch6/src/CheckGenContent/CheckGenContent/test/003F4204-8AD4-E211-9871-E0CB4E5536EF.root'##ZGToLLG
        #fileNames = cms.untracked.vstring('root://xrootd.unl.edu//store/mc/Summer12/TTGJets_8TeV-madgraph/AODSIM/PU_S7_START52_V9-v1/0000/228F159C-AFBE-E111-8F96-003048C6617E.root')
    )
       #fileNames = cms.untracked.vstring('root://xrootd.unl.edu//store/mc/Summer12_DR53X/TTGJets_8TeV-madgraph/AODSIM/PU_S10_START53_V7A-v1/0000/02F29E4D-76D9-E111-A9DC-003048C6763A.root')
)

process.demo = cms.EDAnalyzer('CheckGenContent'
)

process.p = cms.Path(process.demo* process.printTree)
