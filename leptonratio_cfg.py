import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        '/store/user/sapta/STOP200_Publish_8TeV_FastSIM_Rerun/STOP200_Publish_8TeV_FastSIM_Rerun/ef8a184b12798caf874bdadb94795d3b/step0_1000_1_y7Y.root',
        '/store/user/sapta/STOP200_Publish_8TeV_FastSIM_Rerun/STOP200_Publish_8TeV_FastSIM_Rerun/ef8a184b12798caf874bdadb94795d3b/step0_100_1_pyA.root',
        '/store/user/sapta/STOP200_Publish_8TeV_FastSIM_Rerun/STOP200_Publish_8TeV_FastSIM_Rerun/ef8a184b12798caf874bdadb94795d3b/step0_101_1_Dq1.root',
        '/store/user/sapta/STOP200_Publish_8TeV_FastSIM_Rerun/STOP200_Publish_8TeV_FastSIM_Rerun/ef8a184b12798caf874bdadb94795d3b/step0_102_1_Mcr.root',
        '/store/user/sapta/STOP200_Publish_8TeV_FastSIM_Rerun/STOP200_Publish_8TeV_FastSIM_Rerun/ef8a184b12798caf874bdadb94795d3b/step0_103_1_J0P.root',
        '/store/user/sapta/STOP200_Publish_8TeV_FastSIM_Rerun/STOP200_Publish_8TeV_FastSIM_Rerun/ef8a184b12798caf874bdadb94795d3b/step0_104_1_KuH.root',
        '/store/user/sapta/STOP200_Publish_8TeV_FastSIM_Rerun/STOP200_Publish_8TeV_FastSIM_Rerun/ef8a184b12798caf874bdadb94795d3b/step0_105_1_mnj.root',
        '/store/user/sapta/STOP200_Publish_8TeV_FastSIM_Rerun/STOP200_Publish_8TeV_FastSIM_Rerun/ef8a184b12798caf874bdadb94795d3b/step0_106_1_TyZ.root',
        '/store/user/sapta/STOP200_Publish_8TeV_FastSIM_Rerun/STOP200_Publish_8TeV_FastSIM_Rerun/ef8a184b12798caf874bdadb94795d3b/step0_107_1_xGC.root',
        '/store/user/sapta/STOP200_Publish_8TeV_FastSIM_Rerun/STOP200_Publish_8TeV_FastSIM_Rerun/ef8a184b12798caf874bdadb94795d3b/step0_108_1_XVG.root'
        ##Stop 180##
        #'/store/user/sapta/STOP180_Publish_8TeV_FastSIM_Rerun/STOP180_Publish_8TeV_FastSIM_Rerun/ef8a184b12798caf874bdadb94795d3b/step0_1000_1_dEn.root',
        #'/store/user/sapta/STOP180_Publish_8TeV_FastSIM_Rerun/STOP180_Publish_8TeV_FastSIM_Rerun/ef8a184b12798caf874bdadb94795d3b/step0_100_1_PLD.root',
        #'/store/user/sapta/STOP180_Publish_8TeV_FastSIM_Rerun/STOP180_Publish_8TeV_FastSIM_Rerun/ef8a184b12798caf874bdadb94795d3b/step0_101_1_e3d.root',
        #'/store/user/sapta/STOP180_Publish_8TeV_FastSIM_Rerun/STOP180_Publish_8TeV_FastSIM_Rerun/ef8a184b12798caf874bdadb94795d3b/step0_102_1_vsx.root',
        #'/store/user/sapta/STOP180_Publish_8TeV_FastSIM_Rerun/STOP180_Publish_8TeV_FastSIM_Rerun/ef8a184b12798caf874bdadb94795d3b/step0_103_1_iJ9.root',
        #'/store/user/sapta/STOP180_Publish_8TeV_FastSIM_Rerun/STOP180_Publish_8TeV_FastSIM_Rerun/ef8a184b12798caf874bdadb94795d3b/step0_104_1_Rpg.root',
        #'/store/user/sapta/STOP180_Publish_8TeV_FastSIM_Rerun/STOP180_Publish_8TeV_FastSIM_Rerun/ef8a184b12798caf874bdadb94795d3b/step0_105_1_OuI.root',
        #'/store/user/sapta/STOP180_Publish_8TeV_FastSIM_Rerun/STOP180_Publish_8TeV_FastSIM_Rerun/ef8a184b12798caf874bdadb94795d3b/step0_106_1_NfV.root',
        #'/store/user/sapta/STOP180_Publish_8TeV_FastSIM_Rerun/STOP180_Publish_8TeV_FastSIM_Rerun/ef8a184b12798caf874bdadb94795d3b/step0_107_1_loV.root',
        #'/store/user/sapta/STOP180_Publish_8TeV_FastSIM_Rerun/STOP180_Publish_8TeV_FastSIM_Rerun/ef8a184b12798caf874bdadb94795d3b/step0_108_1_kyh.root' 
        ##Stop 160## 
        #'/store/user/sapta/STOP160_Publish_8TeV_FastSIM_Rerun/STOP160_Publish_8TeV_FastSIM_Rerun/ef8a184b12798caf874bdadb94795d3b/step0_1000_1_u3m.root',
        #'/store/user/sapta/STOP160_Publish_8TeV_FastSIM_Rerun/STOP160_Publish_8TeV_FastSIM_Rerun/ef8a184b12798caf874bdadb94795d3b/step0_100_1_Hkn.root',
        #'/store/user/sapta/STOP160_Publish_8TeV_FastSIM_Rerun/STOP160_Publish_8TeV_FastSIM_Rerun/ef8a184b12798caf874bdadb94795d3b/step0_101_1_laf.root',
        #'/store/user/sapta/STOP160_Publish_8TeV_FastSIM_Rerun/STOP160_Publish_8TeV_FastSIM_Rerun/ef8a184b12798caf874bdadb94795d3b/step0_102_1_EzP.root',
        #'/store/user/sapta/STOP160_Publish_8TeV_FastSIM_Rerun/STOP160_Publish_8TeV_FastSIM_Rerun/ef8a184b12798caf874bdadb94795d3b/step0_103_1_uwS.root',
        #'/store/user/sapta/STOP160_Publish_8TeV_FastSIM_Rerun/STOP160_Publish_8TeV_FastSIM_Rerun/ef8a184b12798caf874bdadb94795d3b/step0_104_1_W12.root',
        #'/store/user/sapta/STOP160_Publish_8TeV_FastSIM_Rerun/STOP160_Publish_8TeV_FastSIM_Rerun/ef8a184b12798caf874bdadb94795d3b/step0_105_1_Lee.root',
        #'/store/user/sapta/STOP160_Publish_8TeV_FastSIM_Rerun/STOP160_Publish_8TeV_FastSIM_Rerun/ef8a184b12798caf874bdadb94795d3b/step0_106_1_iZw.root',
        #'/store/user/sapta/STOP160_Publish_8TeV_FastSIM_Rerun/STOP160_Publish_8TeV_FastSIM_Rerun/ef8a184b12798caf874bdadb94795d3b/step0_107_1_XLf.root',
        #'/store/user/sapta/STOP160_Publish_8TeV_FastSIM_Rerun/STOP160_Publish_8TeV_FastSIM_Rerun/ef8a184b12798caf874bdadb94795d3b/step0_108_1_Esy.root'

   )
)

process.demo = cms.EDAnalyzer('LeptonRatio'
)


process.p = cms.Path(process.demo)
