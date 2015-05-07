/* Code for accessing variables from a non-linear 
Tree and storing relevant information
Author: Saptaparna Bhattcharya

*/
#ifndef FlatTreeCreator_h
#define FlatTreeCreator_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TClonesArray.h>
#include <TVector3.h>

#include "Classes/src/TCJet.h"
#include "Classes/src/TCMET.h"
#include "Classes/src/TCElectron.h"
#include "Classes/src/TCMuon.h"
#include "Classes/src/TCTau.h"
#include "Classes/src/TCPhoton.h"
#include "Classes/src/TCGenJet.h"
#include "Classes/src/TCPrimaryVtx.h"
#include "Classes/src/TCTriggerObject.h"
#include "Classes/src/TCGenParticle.h"
#include "Classes/src/TCEGamma.h"

#include "plugins/HistManager.h"
#include "plugins/TreeManager.h"
#include "plugins/TriggerSelector.h"


// Fixed size dimensions of array or collections stored in the TTree if any.

class FlatTreeCreator : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

   // Declaration of leaf types
   TClonesArray    *patJets;
   TClonesArray    *recoElectrons;
   TClonesArray    *recoMuons;
   TClonesArray    *recoPhotons;
   TCMET           *corrMET;
   TCMET           *CaloMET;
   TCMET           *mvaMET;
   TClonesArray    *genJets;
   TClonesArray    *genParticles;
   TClonesArray    *triggerObjects;
   TClonesArray    *primaryVtx;
   TVector3        *beamSpot;
   Int_t           nPUVertices;
   Float_t         nPUVerticesTrue;
   Bool_t          isRealData;
   UInt_t          runNumber;
   ULong64_t       eventNumber;
   UInt_t          lumiSection;
   UInt_t          bunchCross;
   Float_t         ptHat;
   Float_t         qScale;
   Float_t         evtWeight;
   Float_t         rhoFactor;
   Float_t         rho25Factor;
   Float_t         rhoMuFactor;
   ULong64_t       triggerStatus;
   UInt_t          hltPrescale[64];
   Bool_t          NoiseFilters_isScraping;
   Bool_t          NoiseFilters_isNoiseHcalHBHE;
   Bool_t          NoiseFilters_isNoiseHcalLaser;
   Bool_t          NoiseFilters_isNoiseEcalTP;
   Bool_t          NoiseFilters_isNoiseEcalBE;
   Bool_t          NoiseFilters_isCSCTightHalo;
   Bool_t          NoiseFilters_isCSCLooseHalo;
   Bool_t          NoiseFilters_isNoiseTracking;
   Bool_t          NoiseFilters_isNoiseEEBadSc;
   Bool_t          NoiseFilters_isNoisetrkPOG1;
   Bool_t          NoiseFilters_isNoisetrkPOG2;
   Bool_t          NoiseFilters_isNoisetrkPOG3;

   // List of branches
   TBranch        *b_patJets;   //!
   TBranch        *b_recoElectrons;   //!
   TBranch        *b_recoMuons;   //!
   TBranch        *b_recoPhotons;   //!
   TBranch        *b_corrMET;   //!
   TBranch        *b_CaloMET;
   TBranch        *b_mvaMET;   //!
   TBranch        *b_genJets;   //!
   TBranch        *b_genParticles;   //!
   TBranch        *b_triggerObjects;   //!
   TBranch        *b_primaryVtx;   //!
   TBranch        *b_beamSpot;   //!
   TBranch        *b_nPUVertices;   //!
   TBranch        *b_nPUVerticesTrue;   //!
   TBranch        *b_isRealData;   //!
   TBranch        *b_runNumber;   //!
   TBranch        *b_eventNumber;   //!
   TBranch        *b_lumiSection;   //!
   TBranch        *b_bunchCross;   //!
   TBranch        *b_ptHat;   //!
   TBranch        *b_qScale;   //!
   TBranch        *b_evtWeight;   //!
   TBranch        *b_rhoFactor;   //!
   TBranch        *b_rho25Factor;   //!
   TBranch        *b_rhoMuFactor;   //!
   TBranch        *b_triggerStatus;   //!
   TBranch        *b_hltPrescale;   //!
   TBranch        *b_NoiseFilters;   //!

   //Variables of the output tree
   TFile          *outputFile;
   TTree          *outtree;
   int            run;
   int 		  lumi;
   int            event;            
   std::vector<float> ph_pt;
   std::vector<float> ph_eta;
   std::vector<float> ph_phi;
   std::vector<float> ph_energy;
   std::vector<float> ph_px;
   std::vector<float> ph_py;
   std::vector<float> ph_pz;
   std::vector<float> ph_SigmaIetaIeta;
   std::vector<float> ph_SigmaIetaIphi;
   std::vector<float> ph_SigmaIphiIphi;
   std::vector<int> ph_pixelVeto;
   std::vector<int> ph_conversionVeto; 
   std::vector<float> ph_HoE;
   std::vector<float> ph_preShowerOverRaw;
   std::vector<float> ph_R9;
   std::vector<bool> ph_isSTight;
   std::vector<float> ph_e1x3;
   std::vector<float> ph_e1x5;
   std::vector<float> ph_e2x2;
   std::vector<float> ph_e2x5;
   std::vector<float> ph_e5x1;
   std::vector<float> ph_e5x5;
   std::vector<float> ph_e2x5Max;
   std::vector<float> ph_e2OverE5;
   std::vector<float> ph_seedCrystalEnergy; 
   int            nPhotons;
   std::vector<float> el_pt;
   std::vector<float> el_eta;
   std::vector<float> el_phi;
   std::vector<float> el_energy;
   std::vector<int>   el_charge;
   std::vector<bool>  el_isLoose;
   std::vector<bool>  el_isLooseWG;
   std::vector<bool>  el_isTight;  
   std::vector<float> el_Dxy;
   std::vector<float> el_Dz;
   std::vector<float> el_mvaId_Hzz;
   std::vector<float> el_ip3dSig;
   int            nElectrons;
   std::vector<float> mu_pt;
   std::vector<float> mu_eta;
   std::vector<float> mu_phi;
   std::vector<float> mu_energy;
   std::vector<int>   mu_charge;
   std::vector<bool>  mu_isLoose;
   std::vector<bool>  mu_isLooseWG;
   std::vector<bool>  mu_isTight;
   std::vector<float> mu_Dxy;
   std::vector<float> mu_Dz;
   int            nMuons;
   std::vector<float> jet_pt;
   std::vector<float> jet_eta;
   std::vector<float> jet_phi;
   std::vector<float> jet_energy;
   std::vector<int> jet_mva_loose;
   std::vector<int> jet_mva_tight;
   std::vector<int> jet_mva_medium;
   std::vector<int> jet_cut_loose;
   std::vector<int> jet_cut_tight;
   std::vector<int> jet_cut_medium;
   std::vector<float> jet_btag_csv;
   std::vector<float> jet_pt_Unc;
   std::vector<float> tauPx;
   std::vector<float> tauPy;
   std::vector<float> tauPz;
   std::vector<float> tauE;
   std::vector<float> muPx;
   std::vector<float> muPy;
   std::vector<float> muPz;
   std::vector<float> muE;
   std::vector<float> ePx;
   std::vector<float> ePy;
   std::vector<float> ePz;
   std::vector<float> eE;
   std::vector<float> nuPx;
   std::vector<float> nuPy;
   std::vector<float> nuPz;
   std::vector<float> nuE;
   int            nJets;
   float          MET;
   float          MET_Phi;
   float          MET_Px;
   float          MET_Py;
   float          MET_Signxx;
   float          MET_Signxy;
   float          MET_Signyx;
   float          MET_Signyy;
   float          caloMET;
   float          caloMET_Phi;
   float          caloMET_Px;
   float          caloMET_Py;
   bool           fired_HLTPho;
   bool           fired_HLTPhoId; 
   bool           fired_HLTPhoIdMet;           
   bool           fired_HLTMET100;
   bool           fired_HLT_IsoMu24;
   bool           fired_HLT_Mu13_Mu8;
   bool           fired_HLT_Mu17_Mu8;
   bool           fired_HLT_Mu17_TkMu8; 
   bool           fired_HLT_Mu22_TkMu8; 
   bool           fired_HLT_Mu22_TkMu22;
   bool           fired_HLT_Dimuon0_Jpsi_Muon;
   bool           fired_HLT_Dimuon0_Upsilon_Muon;
   bool           fired_HLT_Dimuon11_Upsilon;
   bool           fired_HLT_Dimuon7_Upsilon;
   bool           fired_HLT_DoubleMu3_4_Dimuon5_Bs_Central;
   bool           fired_HLT_DoubleMu3p5_4_Dimuon5_Bs_Central;
   bool           fired_HLT_DoubleMu4_Dimuon7_Bs_Forward;
   bool           fired_HLT_DoubleMu4_JpsiTk_Displaced;
   bool           fired_HLT_DoubleMu4_Jpsi_Displaced;
   bool           fired_HLT_Mu5_Track2_Jpsi;
   bool           fired_HLT_Mu5_Track3p5_Jpsi;
   bool           fired_HLT_Mu7_Track7_Jpsi;
   bool           fired_HLT_Tau2Mu_ItTrack;
   bool           fired_HLT_Dimuon10_Jpsi;
   bool           fired_HLT_Dimuon5_PsiPrime;
   bool           fired_HLT_Dimuon5_Upsilon;
   bool           fired_HLT_Dimuon7_PsiPrime;
   bool           fired_HLT_Dimuon8_Jpsi;
   bool           fired_HLT_Dimuon8_Upsilon;
   bool           fired_HLT_DoubleMu3p5_LowMassNonResonant_Displaced;
   bool           fired_HLT_DoubleMu3p5_LowMass_Displaced;
   bool           fired_HLT_Mu15_TkMu5_Onia;
   bool           fired_HLT_Mu22_Photon22_CaloIdL;
   int            nVertices;
   float   trigObj1Px;
   float   trigObj1Py;
   float   trigObj1Pz;
   float   trigObj1E;
   float   trigObj2Px;
   float   trigObj2Py;
   float   trigObj2Pz;
   float   trigObj2E;
   float   trigObj3Px;
   float   trigObj3Py;
   float   trigObj3Pz;
   float   trigObj3E;
   float   trigObj4Px;
   float   trigObj4Py;
   float   trigObj4Pz;
   float   trigObj4E;
   TVector3       *pvPosition;
   std::vector<float> el_iso;
   std::vector<float> mu_iso;
   std::vector<float> ph_chIso;
   std::vector<float> ph_nuIso;
   std::vector<float> ph_phIso;
   std::vector<bool> ph_isTight;
   std::vector<bool> ph_isMedium;
   std::vector<bool> ph_isLoose;
   std::vector<bool> ph_phIsoSTight;
   std::vector<bool> ph_phIsoTight;
   std::vector<bool> ph_phIsoMedium;
   std::vector<bool> ph_phIsoLoose;
   std::vector<int> ph_Matched;
   std::vector<float> ph_MatchedPt;
   std::vector<float> ph_MatchedEta;
   std::vector<float> ph_MatchedPhi;
   std::vector<float> ph_MatchedEnergy;
   std::vector<float> ph_MatchedSigmaIetaIeta;
   std::vector<bool> ph_MatchedisTightWS;
   std::vector<bool> ph_MatchedIsoTight;
   std::vector<int>  ph_MatchedpixelVeto;
   std::vector<bool> ph_isTightWS;
   TH1F           *h1_numOfEvents;
   TFile          *inFile;
   float          unskimmedEvents;
   float          unskimmedEventsTotal;
   int            fileCount;
   float          EAEle[7];
   float          PUWeightData;
   float          PUWeightDataSys; 
   //MC info
   std::vector<int> el_Matched;
   std::vector<float> el_MatchedPt;
   std::vector<float> el_MatchedEta;
   std::vector<float> el_MatchedPhi;
   std::vector<float> el_MatchedEnergy;
   std::vector<int> el_Matched_1Mother;
   std::vector<int> el_Matched_2Mother;
   std::vector<int> el_Matched_3Mother;
   std::vector<int> el_Matched_4Mother;
   std::vector<int> el_Matched_5Mother;
   std::vector<int> el_Matched_6Mother;
   std::vector<int> el_Matched_7Mother;
   std::vector<int> el_Matched_8Mother;
   std::vector<int> el_Matched_9Mother;
   std::vector<int> el_Matched_10Mother;
   std::vector<int> mu_Matched;
   std::vector<float> mu_MatchedPt;
   std::vector<float> mu_MatchedEta;
   std::vector<float> mu_MatchedPhi;
   std::vector<float> mu_MatchedEnergy;
   std::vector<int> mu_Matched_1Mother;
   std::vector<int> mu_Matched_2Mother;
   std::vector<int> mu_Matched_3Mother;
   std::vector<int> mu_Matched_4Mother;
   std::vector<int> mu_Matched_5Mother;
   std::vector<int> mu_Matched_6Mother;
   std::vector<int> mu_Matched_7Mother;
   std::vector<int> mu_Matched_8Mother;
   std::vector<int> mu_Matched_9Mother;
   std::vector<int> mu_Matched_10Mother;
   bool Ztt;
   bool Znunu;
   bool Zmumu;
   bool Zee;
   bool Wt;
   bool isPhoton;
   bool isISRPhoton;
   FlatTreeCreator(TTree * /*tree*/ =0) : fChain(0) { }
   virtual ~FlatTreeCreator() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();
   bool isTightElectron(TCElectron *electron);
   bool isMediumElectron(TCElectron *electron);
   bool isLooseElectron(TCElectron *electron);
   bool isLooseElectronWG(TCElectron *electron);
   bool isTightMuon(TCMuon *muon);
   bool isLooseMuon(TCMuon *muon);
   bool isLooseMuonWG(TCMuon *muon);
   bool isSTightPhoton(TCPhoton *photon); 
   bool isTightPhoton(TCPhoton *photon);
   bool isMediumPhoton(TCPhoton *photon);
   bool isLoosePhoton(TCPhoton *photon);
   bool isTightWithoutSieiePhoton(TCPhoton *photon);
   float ElectronIso(TCElectron *electron);
   float MuonIso(TCMuon *muon);
   int PhotonIso(TCPhoton *photon, double &chIso, double &nuIso, double &phIso, bool &isoPassL, bool &isoPassM, bool &isoPassT, bool &isoPassST);
   double mdeltaR(double eta1, double phi1, double eta2, double phi2);
   void fillMotherInfo(TCGenParticle* mother, int i, vector <int> & momid);
 
   ClassDef(FlatTreeCreator,0);
};

#endif

#ifdef FlatTreeCreator_cxx
void FlatTreeCreator::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer

   patJets = 0;
   recoElectrons = 0;
   recoMuons = 0;
   recoPhotons = 0;
   corrMET = 0;
   CaloMET = 0;
   mvaMET = 0;
   genJets = 0;
   genParticles = 0;
   triggerObjects = 0;
   primaryVtx = 0;
   beamSpot = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("patJets", &patJets, &b_patJets);
   fChain->SetBranchAddress("recoElectrons", &recoElectrons, &b_recoElectrons);
   fChain->SetBranchAddress("recoMuons", &recoMuons, &b_recoMuons);
   fChain->SetBranchAddress("recoPhotons", &recoPhotons, &b_recoPhotons);
   fChain->SetBranchAddress("corrMET", &corrMET, &b_corrMET);
   fChain->SetBranchAddress("CaloMET", &CaloMET, &b_CaloMET);
   fChain->SetBranchAddress("mvaMET", &mvaMET, &b_mvaMET);
   fChain->SetBranchAddress("genJets", &genJets, &b_genJets);
   fChain->SetBranchAddress("genParticles", &genParticles, &b_genParticles);
   fChain->SetBranchAddress("triggerObjects", &triggerObjects, &b_triggerObjects);
   fChain->SetBranchAddress("primaryVtx", &primaryVtx, &b_primaryVtx);
   fChain->SetBranchAddress("beamSpot", &beamSpot, &b_beamSpot);
   fChain->SetBranchAddress("nPUVertices", &nPUVertices, &b_nPUVertices);
   fChain->SetBranchAddress("nPUVerticesTrue", &nPUVerticesTrue, &b_nPUVerticesTrue);
   fChain->SetBranchAddress("isRealData", &isRealData, &b_isRealData);
   fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
   fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
   fChain->SetBranchAddress("lumiSection", &lumiSection, &b_lumiSection);
   fChain->SetBranchAddress("bunchCross", &bunchCross, &b_bunchCross);
   fChain->SetBranchAddress("ptHat", &ptHat, &b_ptHat);
   fChain->SetBranchAddress("qScale", &qScale, &b_qScale);
   fChain->SetBranchAddress("evtWeight", &evtWeight, &b_evtWeight);
   fChain->SetBranchAddress("rhoFactor", &rhoFactor, &b_rhoFactor);
   fChain->SetBranchAddress("rho25Factor", &rho25Factor, &b_rho25Factor);
   fChain->SetBranchAddress("rhoMuFactor", &rhoMuFactor, &b_rhoMuFactor);
   fChain->SetBranchAddress("triggerStatus", &triggerStatus, &b_triggerStatus);
   fChain->SetBranchAddress("hltPrescale", hltPrescale, &b_hltPrescale);
   fChain->SetBranchAddress("NoiseFilters", &NoiseFilters_isScraping, &b_NoiseFilters);
}

Bool_t FlatTreeCreator::Notify()
{
   fileCount+= 1;
   inFile = fChain->GetCurrentFile();
   h1_numOfEvents = (TH1F*) inFile->Get("ntupleProducer/numOfEvents");
   unskimmedEvents = h1_numOfEvents->GetBinContent(1);
   cout<<"THIS IS FILE NUMBER: "<<fileCount<<" and it has "<<unskimmedEvents << " events " <<endl;
   unskimmedEventsTotal += unskimmedEvents;
   return kTRUE;
}

#endif // #ifdef FlatTreeCreator_cxx
