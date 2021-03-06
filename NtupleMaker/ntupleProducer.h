// system include files
#include <memory>
#include <string>
#include <map>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/GeometryVector/interface/LocalVector.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"

// Libraries for objects
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/JetID.h"
#include "DataFormats/JetReco/interface/JPTJetCollection.h"
#include "DataFormats/JetReco/interface/JPTJet.h"

#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METCollection.h"
#include "DataFormats/METReco/interface/METFwd.h"

// PAT candidates
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/MHT.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"


#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Luminosity/interface/LumiSummary.h"
#include "DataFormats/Common/interface/ValueMap.h"

// Generator data formats
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/JetMatching/interface/JetFlavour.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourMatching.h"
#include "SimDataFormats/JetMatching/interface/MatchedPartons.h"
#include "SimDataFormats/JetMatching/interface/JetMatchedPartons.h"

// PAT 
//Not using pats no more, commentng out
//#include "DataFormats/PatCandidates/interface/Electron.h"
//#include "DataFormats/PatCandidates/interface/Muon.h"
//#include "DataFormats/PatCandidates/interface/Photon.h"
//#include "DataFormats/PatCandidates/interface/MET.h"
//#include "DataFormats/PatCandidates/interface/Jet.h"
//#include "DataFormats/PatCandidates/interface/Tau.h"

//#include "RecoVertex/PrimaryVertexProducer/interface/VertexHigherPtSquared.h"

// JEC
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
//#include "/uscms_data/d2/sapta/work/SUSYSearch/CMSSW_5_3_11_patch6/src/CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "/uscms_data/d2/sapta/work/SUSYSearch/SignalGeneration/CMSSW_5_3_16_patch1/src/CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"

// Jet associators
#include "RecoJets/JetAssociationAlgorithms/interface/JetTracksAssociationDRCalo.h"
#include "RecoJets/JetAssociationAlgorithms/interface/JetTracksAssociationDRVertex.h"

#include "RecoJets/JetProducers/interface/PileupJetIdAlgo.h"

// EGamma tools
#include "RecoEgamma/PhotonIdentification/interface/PhotonIsolationCalculator.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "EgammaAnalysis/ElectronTools/interface/ElectronEffectiveArea.h"
#include "EgammaAnalysis/ElectronTools/interface/PFIsolationEstimator.h"
//#include "EgammaAnalysis/ElectronTools/interface/ElectronEnergyRegressionEvaluate.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"

//#include "EgammaAnalysis/ElectronTools/interface/ElectronEnergyCalibrator.h"

//#include "EgammaAnalysis/ElectronTools/interface/EGammaMvaEleEstimator.h"
//#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

// Tracking tools, supposedly for track-met:
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/METReco/interface/BeamHaloSummary.h"

//#include "RecoMET/METProducers/interface/TrackMETProducer.h"
//#include "RecoMET/METProducers/interface/ParticleFlowForChargedMETProducer.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
// ntuple storage classes
#include "TCPhysObject.h"
#include "TCPrimaryVtx.h"
#include "TCJet.h"
#include "TCMET.h"
#include "TCMuon.h"
#include "TCEGamma.h"
#include "TCElectron.h"
#include "TCTau.h"
#include "TCPhoton.h"
#include "TCTriggerObject.h"
#include "TCGenJet.h"
#include "TCGenParticle.h"
//#include "TCTrack.h"

// Need for HLT trigger info:
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

//Supercluster footprint removal:
#include "PFIsolation/SuperClusterFootprintRemoval/interface/SuperClusterFootprintRemoval.h"

//Photon Lazytools and ESEff
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "Geometry/CaloTopology/interface/EcalPreshowerTopology.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "Geometry/EcalAlgo/interface/EcalPreshowerGeometry.h"
#include "RecoCaloTools/Navigation/interface/EcalPreshowerNavigator.h"

//Root  stuff
#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TString.h"
#include "TObject.h"
#include "TObjArray.h"
#include "TClonesArray.h"
#include "TRefArray.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "DataFormats/BTauReco/interface/JetTag.h"

using namespace edm;
using namespace std;
using namespace reco;

//
// class declaration
//

struct Filters {		//Filters 
  Bool_t isScraping;
  Bool_t isNoiseHcalHBHE;
  Bool_t isNoiseHcalLaser;
  Bool_t isNoiseEcalTP;
  Bool_t isNoiseEcalBE;
  Bool_t isCSCTightHalo;
  Bool_t isCSCLooseHalo;
  Bool_t isNoiseTracking;
  Bool_t isNoiseEEBadSc;
  Bool_t isNoisetrkPOG1;
  Bool_t isNoisetrkPOG2;
  Bool_t isNoisetrkPOG3;
};


class ntupleProducer : public edm::EDAnalyzer {
 public:
  explicit ntupleProducer(const edm::ParameterSet&);
  ~ntupleProducer();
  
 private:
  virtual void beginJob() ;
  virtual void beginRun(const edm::Run&, const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endLuminosityBlock(const edm::LuminosityBlock&,const edm::EventSetup&);
  virtual void endRun(const edm::Run&, const edm::EventSetup&);
  virtual void endJob() ;
  
  virtual bool  triggerDecision(edm::Handle<edm::TriggerResults>& hltR, int iTrigger);
  virtual float sumPtSquared(const Vertex& v);
  //virtual bool  associateJetToVertex(reco::PFJet inJet, Handle<reco::VertexCollection> vtxCollection, TCJet* outJet);   
  virtual void  electronMVA(const reco::GsfElectron* iElectron, TCElectron* eleCon, const edm::Event& iEvent,const edm::EventSetup& iSetup, const reco::PFCandidateCollection& PFCandidates, float Rho);
  virtual bool  isFilteredOutScraping(const edm::Event& iEvent, const edm::EventSetup& iSetup, int numtrack=10, double thresh=0.25);
  void analyzeTrigger(edm::Handle<edm::TriggerResults> &hltR, edm::Handle<trigger::TriggerEvent> &hltE, const std::string& triggerName, int* trigCount);                   
  void initJetEnergyCorrector(const edm::EventSetup &iSetup, bool isData);
  TCGenParticle* addGenParticle(const reco::GenParticle* myParticle, int& genPartCount, std::map<const reco::GenParticle*,TCGenParticle*>& genMap);

  //TCTrack::ConversionInfo CheckForConversions(const edm::Handle<reco::ConversionCollection> &convCol,
  //                                            const reco::GsfTrackRef &trk,
  //                                            const math::XYZPoint &bs, const math::XYZPoint &pv);
  vector<float> getESHits(double X, double Y, double Z, map<DetId, EcalRecHit> rechits_map, const CaloSubdetectorGeometry*& geometry_p, CaloSubdetectorTopology *topology_p, int row=0);
  vector<float> getESEffSigmaRR(vector<float> ESHits0);
  //virtual float MatchBTagsToJets(const reco::JetTagCollection, const reco::PFJet);
  virtual float MatchBTagsToJets(const reco::JetTagCollection, vector<pat::Jet>::const_iterator ipJet);
  // ----------member data ---------------------------
  
  struct JetCompare :
  public std::binary_function<reco::Jet, reco::Jet, bool> {
    inline bool operator () (const reco::Jet &j1,
                             const reco::Jet &j2) const
    { return (j1.p4().Pt() > j2.p4().Pt()); }
  };
  
  static bool EnergySortCriterium( const TCPhoton::CrystalInfo p1,const TCPhoton::CrystalInfo p2 ){
    return p1.energy > p2.energy;
  };
  
  typedef std::map<reco::Jet, unsigned int, JetCompare> flavourMap;
  typedef reco::JetTagCollection::const_iterator tag_iter;
  
  //Standard event info
  ULong64_t   eventNumber;
  UInt_t      runNumber, lumiSection, bunchCross, nEvents;
  float ptHat, qScale, evtWeight;
  float deliveredLumi, recordedLumi, lumiDeadTime;
  float rhoFactor, rho25Factor, rhoMuFactor,;
  vector<string>  savedTriggerNames;
  
  edm::Service<TFileService> fs;
  TTree* eventTree;
  TTree* jobTree;

  edm::InputTag jetTag_;
  string        jecTag_;
  edm::InputTag genJetTag_;
  edm::InputTag muonTag_;
  edm::InputTag electronTag_;
  edm::InputTag photonTag_;
  edm::InputTag primaryVtxTag_;
  edm::InputTag triggerResultsTag_;
  edm::InputTag rhoCorrTag_, rho25CorrTag_, rhoMuCorrTag_;
  edm::InputTag hcalHBHEFilterTag_;
  edm::InputTag ecalTPFilterTag_;
  edm::InputTag ecalBEFilterTag_;
  edm::InputTag hcalLaserFilterTag_;
  edm::InputTag trackingFailureTag_;
  edm::InputTag eeBadScFilterTag_;
  edm::InputTag trkPOGFiltersTag1_;
  edm::InputTag trkPOGFiltersTag2_;
  edm::InputTag trkPOGFiltersTag3_;
  edm::InputTag partFlowTag_;
  edm::ParameterSet photonIsoCalcTag_;
  edm::ParameterSet jetPUIdAlgo_;
  edm::InputTag triggerEventTag_;

  //allMet

  edm::InputTag mPhoUp;
  edm::InputTag mPhoDown;
  edm::InputTag mPhoMVAUp;
  edm::InputTag mPhoMVADown;

  edm::InputTag mNoOverlapJet;
  edm::InputTag mSelectedPatJet;
  edm::InputTag mSmearedPatJet;
  
  edm::InputTag mMetMVA;
  edm::InputTag mMVAJERup;
  edm::InputTag mMVAJERdown;
  edm::InputTag mMVAPhoup;
  edm::InputTag mMVAPhodown;
  edm::InputTag mMVAJetup;
  edm::InputTag mMVAJetdown;
  edm::InputTag mMVAUncup;
  edm::InputTag mMVAUncdown;

  edm::InputTag mMetRaw;
  edm::InputTag mMetPf;
  edm::InputTag mMetType01;
  edm::InputTag mMetJERup;
  edm::InputTag mMetJERdown;
  edm::InputTag mMetPhoup;
  edm::InputTag mMetPhodown;
  edm::InputTag mMetJetup;
  edm::InputTag mMetJetdown;
  edm::InputTag mMetUncup;
  edm::InputTag mMetUncdown;
  edm::InputTag mCaloMetMuCorr;
  edm::InputTag mCaloMet; 
  edm::InputTag ebReducedRecHitCollection_;
  edm::InputTag eeReducedRecHitCollection_;
  edm::InputTag esReducedRecHitCollection_;



  bool skimLepton_;
  bool saveMuons_;
  bool saveJets_;
  bool saveElectrons_;
  bool saveEleCrystals_;
  bool savePhotons_;
  bool savePhoCrystals_;
  bool saveMET_;
  bool saveMETExtra_;
  bool saveGenJets_;
  bool saveGenParticles_;
  bool isRealData;
  bool verboseTrigs;
  bool verboseMVAs;
  bool saveMoreEgammaVars_;
  double SCFPRemovalCone_;
  
  //Physics object containers
  TClonesArray* recoJets;
  TClonesArray* patJets;
  TClonesArray* recoJPT;
  TClonesArray* recoMuons;
  TClonesArray* recoElectrons;
  TClonesArray* recoPhotons;
  TClonesArray* triggerObjects;
  TClonesArray* genJets;
  TClonesArray* genParticles;

  TClonesArray* pho_up;
  TClonesArray* pho_down;
  TClonesArray* pho_mva_up;
  TClonesArray* pho_mva_down;

  TClonesArray* jetCon_pat;
  TClonesArray* jetCon_smear;
  

  auto_ptr<TCMET> corrMET;
  auto_ptr<TCMET> CaloMETMuCorr;
  auto_ptr<TCMET> CaloMET;
  auto_ptr<TCMET> mvaMET;

  auto_ptr<TCMET> jerUpMET;
  auto_ptr<TCMET> jerDownMET;
  auto_ptr<TCMET> jerUpMVAMET;
  auto_ptr<TCMET> jerDownMVAMET;

  auto_ptr<TCMET> phoUpMET;
  auto_ptr<TCMET> phoDownMET;
  auto_ptr<TCMET> phoUpMVAMET;
  auto_ptr<TCMET> phoDownMVAMET;

  auto_ptr<TCMET> jetupMET;
  auto_ptr<TCMET> jetdownMET;
  auto_ptr<TCMET> jetupMVAMET;
  auto_ptr<TCMET> jetdownMVAMET;

  auto_ptr<TCMET> uncUpMET;
  auto_ptr<TCMET> uncDownMET;
  auto_ptr<TCMET> uncUpMVAMET;
  auto_ptr<TCMET> uncDownMVAMET;

  //Vertex info
  TClonesArray* primaryVtx;
  TVector3*     beamSpot;
  unsigned      nPUVertices;
  float         nPUVerticesTrue;
  
  //Triggers
  HLTConfigProvider hltConfig_;
  string            hlTriggerResults_, hltProcess_, triggerName_;
  TriggerNames      triggerNames;
  vector<string>    hlNames;
  vector<string>    triggerPaths_;
  ULong64_t         triggerStatus;
  unsigned          hltPrescale[64];
  
  // Technical filters
  Filters myNoiseFilters;

  //Isolator
  PFIsolationEstimator phoIsolator;
  PFIsolationEstimator eleIsolator;

  // Histograms
  TH1F * h1_numOfEvents;

  // Electron Regression
  //auto_ptr<ElectronEnergyRegressionEvaluate> myEleReg;

  // PU Jet Id Algo
  auto_ptr<PileupJetIdAlgo> myPUJetID;
  auto_ptr<FactorizedJetCorrector> jecCor;  

  // Photon LazyTool and such
  auto_ptr<EcalClusterLazyTools> lazyTool;
  map<DetId, EcalRecHit> rechits_map_;
  auto_ptr<CaloSubdetectorTopology> topology_p;

};
