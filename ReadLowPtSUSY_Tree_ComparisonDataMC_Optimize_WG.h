#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TChain.h>
#include <TLorentzVector.h>
#include <TLegend.h>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <algorithm>
#include <TGraphAsymmErrors.h>
#include <TVector2.h>
#include <TF1.h>
#include <TMatrixD.h>
#include "/uscms_data/d2/lpcljm/sapta/SUSYSearch/CMSSW_5_3_16_patch1/src/TauAnalysis/SVfitStandalone/interface/SVfitStandaloneLikelihood.h"
#include "/uscms_data/d2/lpcljm/sapta/SUSYSearch/CMSSW_5_3_16_patch1/src/TauAnalysis/SVfitStandalone/interface/LikelihoodFunctions.h"
#include "/uscms_data/d2/lpcljm/sapta/SUSYSearch/CMSSW_5_3_16_patch1/src/TauAnalysis/SVfitStandalone/interface/SVfitStandaloneMarkovChainIntegrator.h"
#include "/uscms_data/d2/lpcljm/sapta/SUSYSearch/CMSSW_5_3_16_patch1/src/TauAnalysis/SVfitStandalone/interface/svFitStandaloneAuxFunctions.h"
using namespace svFitStandalone;
using svFitStandalone::Vector;
using svFitStandalone::LorentzVector;
using svFitStandalone::MeasuredTauLepton;
#include "/uscms_data/d2/lpcljm/sapta/SUSYSearch/CMSSW_5_3_16_patch1/src/TauAnalysis/SVfitStandalone/interface/SVfitStandaloneAlgorithm.h"

TLorentzVector fillTLorentzVector(double pT, double eta, double phi, double E)
{
  TLorentzVector object_p4;
  object_p4.SetPtEtaPhiE(pT, eta, phi, E);
  return object_p4;
}

typedef struct
{
  float pT;
  float eta;
  float phi;
  float energy;
  int charge;
  bool isTight;
  bool isLoose;
  bool isLooseWG;
  float isolation;
  float dxy;
  float dz;
  int matched;
  int mother1;
  int mother2;
  int mother3;
  int mother4;
  int mother5;
  int mother6;
  int mother7;
  int mother8;
  int mother9;
  int mother10;
} LeptonInfo;

typedef struct
{
  TLorentzVector LepLV;
  float Isolation;
  int Charge;
  float Dxy;
  float Dz;
  int Matched;
  int Mother1;
  int Mother2;
  int Mother3;
  int Mother4;
  int Mother5;
  int Mother6;
  int Mother7;
  int Mother8;
  int Mother9;
  int Mother10;
} AnalysisLeptonInfo;

typedef struct
{
  float pT;
  float eta;
  float phi;
  float energy;
  int charge;
  bool isTight;
  bool isMedium;
  bool isLoose;
  float chIsolation;
  float nuIsolation;
  float phIsolation;
  bool  phIsoTight;
  bool  phIsoMedium;
  bool  phIsoLoose;
  float phHoE;
  int   phconversionVeto;
  int   phpixelVeto;
  float phSigmaIetaIeta;
  float phSigmaIetaIphi;
  float phSigmaIphiIphi;
  float phpreShowerOverRaw;
  float phR9;
  float phe1x5;
  float phe1x3;
  float phe2x2;
  float phe2x5;
  float phe5x1;
  float phe5x5;
  float phe2x5Max;
  float phe2OverE5;
  float phseedCrystalEnergy;
  int Matched;
} PhotonInfo;

typedef struct
{
  float pT;
  float eta;
  float phi;
  float energy;
  int   PU_mva_loose;
  int   PU_mva_tight;
  int   PU_mva_medium;
  int   PU_cut_loose;
  int   PU_cut_tight;
  int   PU_cut_medium;
  float BTag_csv;
} JetInfo;

typedef struct
{
 TLorentzVector JetLV; 
 float BTag_CSV; 
} AnalysisJetInfo;

typedef struct
{
  TLorentzVector BJetLV;
  float BTag;
} BJetInfo;

typedef struct
{
  float Px;
  float Py;
  float Pz;
  float E;
} TriggerInfo;

typedef struct
{
  int matched;
  float pT;
  float eta;
  float phi;
  float energy;
} MatchedLeptonInfo;

typedef struct
{
  float Px;
  float Py;
  float Pz;
  float E;
} GenParticleInfo;


typedef struct
{
  TH1F *h_mu_pt_leading;
  TH1F *h_mu_pt_trailing; 
  TH1F *h_mu_eta_leading;
  TH1F *h_mu_eta_trailing;
  TH1F *h_mu_phi_leading;
  TH1F *h_mu_phi_trailing;
  TH1F *h_mu_energy_leading;
  TH1F *h_mu_energy_trailing;
  TH1F *h_mu_dxy_leading;
  TH1F *h_mu_dxy_trailing;
  TH1F *h_mu_dz_leading;
  TH1F *h_mu_dz_trailing;
  TH1F *h_mu_charge_leading;
  TH1F *h_mu_charge_trailing;
  TH1F *h_el_pt_leading;    
  TH1F *h_el_pt_trailing;
  TH1F *h_el_eta_leading;
  TH1F *h_el_eta_trailing;
  TH1F *h_el_phi_leading;
  TH1F *h_el_phi_trailing;
  TH1F *h_el_energy_leading;
  TH1F *h_el_energy_trailing; 
  TH1F *h_el_dxy_leading;
  TH1F *h_el_dxy_trailing;
  TH1F *h_el_dz_leading;
  TH1F *h_el_dz_trailing;
  TH1F *h_el_charge_leading;
  TH1F *h_el_charge_trailing;
  TH1F *h_TransverseMass_leadingMu;
  TH1F *h_TransverseMass_trailingMu;
  TH1F *h_TransverseMass_leadingEl;
  TH1F *h_TransverseMass_trailingEl;
  TH1F *h_InvariantMass_Ph_leadingMu;
  TH1F *h_InvariantMass_Ph_trailingMu;
  TH1F *h_InvariantMass_Ph_leadingEl;
  TH1F *h_InvariantMass_Ph_trailingEl;
  TH1F *h_InvariantMass;
  TH1F *h_InvariantMass_Ph;
  TH1F *h_Invariant_Mass_Taus;
  TH1F *h_Invariant_Mass_Mus;
  TH1F *h_Invariant_Mass_Es;
  TH1F *h_SVFit;
  TH1F *h_SVFitMass_Ph;
  TH1F *h_SVFit_Mass_Taus;
  TH1F *h_DeltaPhi_met_mu1;
  TH1F *h_DeltaPhi_met_mu2;
  TH1F *h_DeltaPhi_ph_mu1;
  TH1F *h_DeltaPhi_ph_mu2;
  TH1F *h_DeltaPhi_met_el1;
  TH1F *h_DeltaPhi_met_el2;
  TH1F *h_DeltaPhi_ph_el1;
  TH1F *h_DeltaPhi_ph_el2;
  TH1F *h_Isolation_mu1;
  TH1F *h_Isolation_mu2;
  TH1F *h_Isolation_el1;
  TH1F *h_Isolation_el2;
  TH1F *h_nVertices;
  TH1F *h_caloMET;
  TH1F *h_MET;
  TH1F *h_MET_Signxx;
  TH1F *h_MET_Signxy;
  TH1F *h_MET_Signyx;
  TH1F *h_MET_Signyy;
  TH1F *h_HT;
  TH1F *h_HTb;
  TH1F *h_nJets;
  TH1F *h_nbJets;
  TH1F *h_jet_pt_leading;
  TH1F *h_jet_pt_trailing;
  TH1F *h_jet_pt_3rd;
  TH1F *h_jet_pt_4th;
  TH1F *h_jet_pt_5th;
  TH1F *h_jet_pt_6th;
  TH1F *h_jet_eta_leading;
  TH1F *h_jet_eta_trailing;
  TH1F *h_jet_eta_3rd;
  TH1F *h_jet_eta_4th;
  TH1F *h_jet_eta_5th;
  TH1F *h_jet_eta_6th;  
  TH1F *h_jet_phi_leading;
  TH1F *h_jet_phi_trailing;
  TH1F *h_jet_phi_3rd;
  TH1F *h_jet_phi_4th;
  TH1F *h_jet_phi_5th;
  TH1F *h_jet_phi_6th;
  TH1F *h_jet_energy_leading;
  TH1F *h_jet_energy_trailing;
  TH1F *h_jet_energy_3rd;
  TH1F *h_jet_energy_4th;
  TH1F *h_jet_energy_5th;
  TH1F *h_jet_energy_6th; 
  TH1F *h_bjet_pt_leading;
  TH1F *h_bjet_pt_trailing;
  TH1F *h_bjet_eta_leading;
  TH1F *h_bjet_eta_trailing;
  TH1F *h_bjet_phi_leading;
  TH1F *h_bjet_phi_trailing;
  TH1F *h_bjet_energy_leading;
  TH1F *h_bjet_energy_trailing;
  TH1F *h_photon_pt;
  TH1F *h_photon_eta;
  TH1F *h_photon_phi;
  TH1F *h_photon_energy;

} HistCollection;

void initializeHistCollection(HistCollection &histCol, std::string suffix)
{
histCol.h_mu_pt_leading = new TH1F(("h_mu_pt_leading_"+suffix).c_str(), "Leading muon pT; pT [GeV]; Events/GeV", 1000, 0, 1000); histCol.h_mu_pt_leading->Sumw2();
histCol.h_mu_pt_trailing = new TH1F(("h_mu_pt_trailing_"+suffix).c_str(), "Trailing muon pT; pT [GeV]; Events/GeV", 1000, 0, 1000); histCol.h_mu_pt_trailing->Sumw2();
histCol.h_mu_eta_leading = new TH1F(("h_mu_eta_leading_"+suffix).c_str(), "Leading muon #eta ; #eta ; Events", 600, -3.0, 3.0); histCol.h_mu_eta_leading->Sumw2();
histCol.h_mu_eta_trailing = new TH1F(("h_mu_eta_trailing_"+suffix).c_str(), "Trailing muon #eta ; #eta ; Events", 600, -3.0, 3.0); histCol.h_mu_eta_trailing->Sumw2();
histCol.h_mu_phi_leading = new TH1F(("h_mu_phi_leading_"+suffix).c_str(), "Leading muon #phi ; #phi ; Events", 800, -4.0, 4.0); histCol.h_mu_phi_leading->Sumw2();
histCol.h_mu_phi_trailing = new TH1F(("h_mu_phi_trailing_"+suffix).c_str(), "Trailing muon #phi ; #phi ; Events", 800, -4.0, 4.0); histCol.h_mu_phi_trailing->Sumw2();
histCol.h_mu_energy_leading = new TH1F(("h_mu_energy_leading_"+suffix).c_str(), "Leading muon Energy; Energy [GeV]; Events/GeV", 1000, 0, 1000); histCol.h_mu_energy_leading->Sumw2();
histCol.h_mu_dxy_leading = new TH1F(("h_mu_dxy_leading_"+suffix).c_str(), "Leading muon dxy; dxy [cm]; Events/cm", 2000000, -10, 10); histCol.h_mu_dxy_leading->Sumw2();
histCol.h_mu_dxy_trailing = new TH1F(("h_mu_dxy_trailing_"+suffix).c_str(), "Trailing muon dxy; dxy [cm]; Events/cm", 2000000, -10, 10); histCol.h_mu_dxy_trailing->Sumw2();
histCol.h_mu_dz_leading = new TH1F(("h_mu_dz_leading_"+suffix).c_str(), "Leading muon dz; dz [cm]; Events/cm", 2000000, -10, 10); histCol.h_mu_dz_leading->Sumw2();
histCol.h_mu_dz_trailing = new TH1F(("h_mu_dz_trailing_"+suffix).c_str(), "Trailing muon dz; dz [cm]; Events/cm", 2000000, -10, 10); histCol.h_mu_dz_trailing->Sumw2();
histCol.h_Isolation_mu1 = new TH1F(("h_Isolation_mu1_"+suffix).c_str(),"Leading muon isolation; Isolation; Events", 100000, 0, 1.0);histCol.h_Isolation_mu1->Sumw2();
histCol.h_Isolation_mu2 = new TH1F(("h_Isolation_mu2_"+suffix).c_str(),"Trailing muon isolation; Isolation; Events", 100000, 0, 1.0);histCol.h_Isolation_mu2->Sumw2();
histCol.h_mu_energy_trailing = new TH1F(("h_mu_energy_trailing_"+suffix).c_str(), "Trailing muon Energy; Energy [GeV]; Events/GeV", 1000, 0, 1000); histCol.h_mu_energy_trailing->Sumw2();
histCol.h_mu_charge_leading = new TH1F(("h_mu_charge_leading_"+suffix).c_str(), "Leading muon charge; Charge; Events", 10, -5, 5); histCol.h_mu_charge_leading->Sumw2();
histCol.h_mu_charge_trailing = new TH1F(("h_mu_charge_trailing_"+suffix).c_str(), "Trailing muon charge; Charge; Events", 10, -5, 5); histCol.h_mu_charge_trailing->Sumw2();
histCol.h_el_pt_leading = new TH1F(("h_el_pt_leading_"+suffix).c_str(), "Leading electron pT; pT [GeV]; Events/GeV", 1000, 0, 1000); histCol.h_el_pt_leading->Sumw2();
histCol.h_el_pt_trailing = new TH1F(("h_el_pt_trailing_"+suffix).c_str(), "Trailing electron pT; pT [GeV]; Events/GeV", 1000, 0, 1000); histCol.h_el_pt_trailing->Sumw2();
histCol.h_el_dxy_leading = new TH1F(("h_el_dxy_leading_"+suffix).c_str(), "Leading electron dxy; dxy [cm]; Events/cm", 200000, -10, 10); histCol.h_el_dxy_leading->Sumw2();
histCol.h_el_dxy_trailing = new TH1F(("h_el_dxy_trailing_"+suffix).c_str(), "Trailing electron dxy; dxy [cm]; Events/cm", 200000, -10, 10); histCol.h_el_dxy_trailing->Sumw2();
histCol.h_el_dz_leading = new TH1F(("h_el_dz_leading_"+suffix).c_str(), "Leading electron dz; dz [cm]; Events/cm", 200000, -10, 10); histCol.h_el_dz_leading->Sumw2();
histCol.h_el_dz_trailing = new TH1F(("h_el_dz_trailing_"+suffix).c_str(), "Trailing electron dz; dz [cm]; Events/cm", 200000, 10, 10); histCol.h_el_dz_trailing->Sumw2();
histCol.h_el_eta_leading = new TH1F(("h_el_eta_leading_"+suffix).c_str(), "Leading electron #eta ; #eta ; Events", 600, -3.0, 3.0); histCol.h_el_eta_leading->Sumw2();
histCol.h_el_eta_trailing = new TH1F(("h_el_eta_trailing_"+suffix).c_str(), "Trailing electron #eta ; #eta ; Events", 600, -3.0, 3.0); histCol.h_el_eta_trailing->Sumw2();
histCol.h_el_phi_leading = new TH1F(("h_el_phi_leading_"+suffix).c_str(), "Leading electron #phi ; #phi ; Events", 800, -4.0, 4.0); histCol.h_el_phi_leading->Sumw2();
histCol.h_el_phi_trailing = new TH1F(("h_el_phi_trailing_"+suffix).c_str(), "Trailing electron #phi ; #phi ; Events", 800, -4.0, 4.0); histCol.h_el_phi_trailing->Sumw2();
histCol.h_el_energy_leading = new TH1F(("h_el_energy_leading_"+suffix).c_str(), "Leading electron Energy; Energy [GeV]; Events/GeV", 1000, 0, 1000); histCol.h_el_energy_leading->Sumw2();
histCol.h_el_energy_trailing = new TH1F(("h_el_energy_trailing_"+suffix).c_str(), "Trailing electron Energy; Energy [GeV]; Events/GeV", 1000, 0, 1000); histCol.h_el_energy_trailing->Sumw2();
histCol.h_Isolation_el1 = new TH1F(("h_Isolation_el1_"+suffix).c_str(),"Leading electron isolation; Isolation; Events", 100000, 0, 1.0);histCol.h_Isolation_el1->Sumw2();
histCol.h_Isolation_el2 = new TH1F(("h_Isolation_el2_"+suffix).c_str(),"Trailing electron isolation; Isolation; Events", 100000, 0, 1.0);histCol.h_Isolation_el2->Sumw2();
histCol.h_el_charge_leading = new TH1F(("h_el_charge_leading_"+suffix).c_str(), "Leading electron charge; Charge; Events", 10, -5, 5); histCol.h_el_charge_leading->Sumw2();
histCol.h_el_charge_trailing = new TH1F(("h_el_charge_trailing_"+suffix).c_str(), "Trailing electron charge; Charge; Events", 10, -5, 5); histCol.h_el_charge_trailing->Sumw2();
histCol.h_TransverseMass_leadingMu=new TH1F(("h_TransverseMass_leadingMu_"+suffix).c_str(), "Transverse mass (leading muon + MET); m_{T} [GeV]; Events/GeV", 9000, 0, 300); histCol.h_TransverseMass_leadingMu->Sumw2();
histCol.h_TransverseMass_trailingMu=new TH1F(("h_TransverseMass_trailingMu_"+suffix).c_str(), "Transverse mass (trailing muon + MET); m_{T} [GeV]; Events/GeV", 9000, 0, 300); histCol.h_TransverseMass_trailingMu->Sumw2();
histCol.h_TransverseMass_leadingEl=new TH1F(("h_TransverseMass_leadingEl_"+suffix).c_str(), "Transverse mass (leading electron + MET); m_{T} [GeV]; Events/GeV", 9000, 0, 300); histCol.h_TransverseMass_leadingEl->Sumw2();
histCol.h_TransverseMass_trailingEl=new TH1F(("h_TransverseMass_trailingEl_"+suffix).c_str(), "Transverse mass (trailing electron + MET); m_{T} [GeV]; Events/GeV", 9000, 0, 300); histCol.h_TransverseMass_trailingEl->Sumw2();
histCol.h_InvariantMass=new TH1F(("h_InvariantMass_"+suffix).c_str(), "Di-lepton invariant mass; m_{ll} [GeV]; Events/GeV", 9000, 0, 300); histCol.h_InvariantMass->Sumw2();
histCol.h_InvariantMass_Ph=new TH1F(("h_InvariantMass_Ph_"+suffix).c_str(), "Di-lepton and photon invariant mass; m_{ll#gamma} [GeV]; Events/GeV", 9000, 0, 300); histCol.h_InvariantMass_Ph->Sumw2();
histCol.h_InvariantMass_Ph_leadingMu =new TH1F(("h_InvariantMass_Ph_leadingMu_"+suffix).c_str(), "Invariant Mass of leading lepton and photon; m_{l#gamma} [GeV]; Events/GeV", 9000, 0, 300); histCol. h_InvariantMass_Ph_leadingMu->Sumw2();
histCol.h_InvariantMass_Ph_trailingMu =new TH1F(("h_InvariantMass_Ph_trailingMu_"+suffix).c_str(), "Invariant Mass of trailing lepton and photon; m_{l#gamma} [GeV]; Events/GeV", 9000, 0, 300); histCol.h_InvariantMass_Ph_trailingMu->Sumw2(); 
histCol.h_InvariantMass_Ph_leadingEl =new TH1F(("h_InvariantMass_Ph_leadingEl_"+suffix).c_str(), "Invariant Mass of leading lepton and photon; m_{l#gamma} [GeV]; Events/GeV", 9000, 0, 300); histCol.h_InvariantMass_Ph_leadingEl->Sumw2();
histCol.h_InvariantMass_Ph_trailingEl =new TH1F(("h_InvariantMass_Ph_trailingEl_"+suffix).c_str(), "Invariant Mass of trailing lepton and photon; m_{l#gamma} [GeV]; Events/GeV", 9000, 0, 300); histCol.h_InvariantMass_Ph_trailingEl->Sumw2();   
histCol.h_SVFitMass_Ph = new TH1F(("h_SVFitMass_Ph"+suffix).c_str(), "SVFit and photon mass; m_{ll#gamma} [GeV]; Events/GeV", 15000, 0, 500); histCol.h_SVFitMass_Ph->Sumw2();
histCol.h_Invariant_Mass_Taus=new TH1F(("h_Invariant_Mass_Taus_"+suffix).c_str(), "Invariant mass of the tau pair (with MET); m_{#tau#tau} [GeV]; Events/GeV", 9000, 0, 300); histCol.h_Invariant_Mass_Taus->Sumw2();
histCol.h_Invariant_Mass_Mus=new TH1F(("h_Invariant_Mass_Mus_"+suffix).c_str(), "Invariant mass of the muon pair; m_{#mu#mu} [GeV]; Events/GeV", 9000, 0, 300); histCol.h_Invariant_Mass_Mus->Sumw2();
histCol.h_Invariant_Mass_Es=new TH1F(("h_Invariant_Mass_Es_"+suffix).c_str(), "Invariant mass of the electron pair; m_{ee} [GeV]; Events/GeV", 9000, 0, 300); histCol.h_Invariant_Mass_Es->Sumw2();
histCol.h_SVFit = new TH1F(("h_SVFit_"+suffix).c_str(), "SVFit mass; m_{ll} [GeV]; Events/GeV", 9000, 0, 300); histCol.h_SVFit->Sumw2();
histCol.h_SVFit_Mass_Taus = new TH1F(("h_SVFit_Mass_Taus_"+suffix).c_str(), "SVFit mass with taus; m_{#tau#tau} [GeV]; Events/GeV", 9000, 0, 300); histCol.h_SVFit_Mass_Taus->Sumw2();
histCol.h_DeltaPhi_met_mu1 = new TH1F(("h_DeltaPhi_met_mu1_"+suffix).c_str(), "#Delta #phi between MET and the leading muon; #Delta #phi(MET, leading muon); Events/GeV", 3500, -3.5, 3.5); histCol.h_DeltaPhi_met_mu1->Sumw2();
histCol.h_DeltaPhi_met_mu2 = new TH1F(("h_DeltaPhi_met_mu2_"+suffix).c_str(), "#Delta #phi between MET and the trailing muon; #Delta #phi(MET, trailing muon); Events/GeV", 3500, -3.5, 3.5); histCol.h_DeltaPhi_met_mu2->Sumw2();
histCol.h_DeltaPhi_met_el1 = new TH1F(("h_DeltaPhi_met_el1_"+suffix).c_str(), "#Delta #phi between MET and the leading electron; #Delta #phi(MET, leading electron); Events/GeV", 3500, -3.5, 3.5); histCol.h_DeltaPhi_met_el1->Sumw2();
histCol.h_DeltaPhi_met_el2 = new TH1F(("h_DeltaPhi_met_el2_"+suffix).c_str(), "#Delta #phi between MET and the trailing electron; #Delta #phi(MET, trailing electron); Events/GeV", 3500, -3.5, 3.5); histCol.h_DeltaPhi_met_el2->Sumw2();
histCol.h_DeltaPhi_ph_mu1 = new TH1F(("h_DeltaPhi_ph_mu1_"+suffix).c_str(), "#Delta #phi between the photon and the leading muon; #Delta #phi(#gamma, leading muon); Events/GeV", 3500, -3.5, 3.5); histCol.h_DeltaPhi_ph_mu1->Sumw2();
histCol.h_DeltaPhi_ph_mu2 = new TH1F(("h_DeltaPhi_ph_mu2_"+suffix).c_str(), "#Delta #phi between the photon and the trailing muon; #Delta #phi(#gamma, trailing muon); Events/GeV", 3500, -3.5, 3.5); histCol.h_DeltaPhi_ph_mu2->Sumw2();
histCol.h_DeltaPhi_ph_el1 = new TH1F(("h_DeltaPhi_ph_el1_"+suffix).c_str(), "#Delta #phi between the photon and the leading electron; #Delta #phi(#gamma, leading electron); Events/GeV", 3500, -3.5, 3.5); histCol.h_DeltaPhi_ph_el1->Sumw2();
histCol.h_DeltaPhi_ph_el2 = new TH1F(("h_DeltaPhi_ph_el2_"+suffix).c_str(), "#Delta #phi between the photon and the trailing electron; #Delta #phi(#gamma, trailing electron); Events/GeV", 3500, -3.5, 3.5); histCol.h_DeltaPhi_ph_el2->Sumw2();
histCol.h_nVertices=new TH1F(("h_nVertices_"+suffix).c_str(), "Number of vertices; nvertices; Events", 50, -0.5, 49.5); histCol.h_nVertices->Sumw2();
histCol.h_caloMET=new TH1F(("h_caloMET_"+suffix).c_str(), "Missing ET; MET [GeV]; Events/GeV", 600, 0, 600); histCol.h_caloMET->Sumw2();
histCol.h_MET=new TH1F(("h_MET_"+suffix).c_str(), "Missing ET; MET [GeV]; Events/GeV", 600, 0, 600); histCol.h_MET->Sumw2();
histCol.h_MET_Signxx=new TH1F(("h_MET_Signxx_"+suffix).c_str(), "Missing ET Signxx; MET [GeV]; Events/GeV", 1000, 0, 1000); histCol.h_MET_Signxx->Sumw2();
histCol.h_MET_Signxy=new TH1F(("h_MET_Signxy_"+suffix).c_str(), "Missing ET Signxy; MET [GeV]; Events/GeV", 2000, -1000, 1000); histCol.h_MET_Signxy->Sumw2();
histCol.h_MET_Signyx=new TH1F(("h_MET_Signyx_"+suffix).c_str(), "Missing ET Signyx; MET [GeV]; Events/GeV", 2000, -1000, 1000); histCol.h_MET_Signyx->Sumw2();
histCol.h_MET_Signyy=new TH1F(("h_MET_Signyy_"+suffix).c_str(), "Missing ET Signyy; MET [GeV]; Events/GeV", 1000, 0, 1000); histCol.h_MET_Signyy->Sumw2();
histCol.h_HT = new TH1F(("h_HT_"+suffix).c_str(), "HT (scalar sum of jet pT); H_T [GeV]; Events/GeV", 5000, 0, 5000.0);histCol.h_HT->Sumw2();
histCol.h_HTb = new TH1F(("h_HTb_"+suffix).c_str(), "HTb (scalar sum of b-jet pT); H_Tb [GeV]; Events/GeV", 5000, 0, 5000.0);histCol.h_HTb->Sumw2();
histCol.h_nJets = new TH1F(("h_nJets_"+suffix).c_str(), "Number of Jets; Number of Jets; Events", 20, -0.5, 19.5);histCol.h_nJets->Sumw2();
histCol.h_nbJets = new TH1F(("h_nbJets_"+suffix).c_str(), "Number of b-Jets; Number of b-Jets; Events", 20, -0.5, 19.5);histCol.h_nbJets->Sumw2();
histCol.h_jet_pt_leading=new TH1F(("h_jet_pt_leading_"+suffix).c_str(), "Leading jet pT; pT [GeV]; Events/GeV", 10000, 0, 1000); histCol.h_jet_pt_leading->Sumw2();
histCol.h_jet_pt_trailing=new TH1F(("h_jet_pt_trailing_"+suffix).c_str(), "Trailing jet pT; pT [GeV]; Events/GeV", 10000, 0, 1000); histCol.h_jet_pt_trailing->Sumw2();
histCol.h_jet_pt_3rd=new TH1F(("h_jet_pt_3rd_"+suffix).c_str(), "3rd jet pT; pT [GeV]; Events/GeV", 10000, 0, 1000); histCol.h_jet_pt_3rd->Sumw2();
histCol.h_jet_pt_4th=new TH1F(("h_jet_pt_4th_"+suffix).c_str(), "4th jet pT; pT [GeV]; Events/GeV", 10000, 0, 1000); histCol.h_jet_pt_4th->Sumw2();
histCol.h_jet_pt_5th=new TH1F(("h_jet_pt_5th_"+suffix).c_str(), "5th jet pT; pT [GeV]; Events/GeV", 10000, 0, 1000); histCol.h_jet_pt_5th->Sumw2();
histCol.h_jet_pt_6th=new TH1F(("h_jet_pt_6th_"+suffix).c_str(), "6th jet pT; pT [GeV]; Events/GeV", 10000, 0, 1000); histCol.h_jet_pt_6th->Sumw2();
histCol.h_jet_eta_leading=new TH1F(("h_jet_eta_leading_"+suffix).c_str(), "Leading jet #eta; #eta; Events", 600, -3.0, 3.0); histCol.h_jet_eta_leading->Sumw2();
histCol.h_jet_eta_trailing=new TH1F(("h_jet_eta_trailing_"+suffix).c_str(), "Trailing jet #eta; #eta; Events", 600, -3.0, 3.0); histCol.h_jet_eta_trailing->Sumw2();
histCol.h_jet_eta_3rd=new TH1F(("h_jet_eta_3rd_"+suffix).c_str(), "3rd jet #eta; #eta; Events/GeV", 600, -3.0, 3.0); histCol.h_jet_eta_3rd->Sumw2();
histCol.h_jet_eta_4th=new TH1F(("h_jet_eta_4th_"+suffix).c_str(), "4th jet #eta; #eta; Events/GeV", 600, -3.0, 3.0); histCol.h_jet_eta_4th->Sumw2();
histCol.h_jet_eta_5th=new TH1F(("h_jet_eta_5th_"+suffix).c_str(), "5th jet #eta; #eta; Events/GeV", 600, -3.0, 3.0); histCol.h_jet_eta_5th->Sumw2();
histCol.h_jet_eta_6th=new TH1F(("h_jet_eta_6th_"+suffix).c_str(), "6th jet #eta; #eta; Events/GeV", 600, -3.0, 3.0); histCol.h_jet_eta_6th->Sumw2();
histCol.h_jet_phi_leading=new TH1F(("h_jet_phi_leading_"+suffix).c_str(), "Leading jet #phi; #phi; Events/GeV", 800, -4.0, 4.0); histCol.h_jet_phi_leading->Sumw2();
histCol.h_jet_phi_trailing=new TH1F(("h_jet_phi_trailing_"+suffix).c_str(), "Trailing jet #phi; #phi; Events/GeV", 800, -4.0, 4.0); histCol.h_jet_phi_trailing->Sumw2();
histCol.h_jet_phi_3rd=new TH1F(("h_jet_phi_3rd_"+suffix).c_str(), "3rd jet #phi; #phi; Events/GeV", 800, -4.0, 4.0); histCol.h_jet_phi_3rd->Sumw2();
histCol.h_jet_phi_4th=new TH1F(("h_jet_phi_4th_"+suffix).c_str(), "4th jet #phi; #phi; Events/GeV", 800, -4.0, 4.0); histCol.h_jet_phi_4th->Sumw2();
histCol.h_jet_phi_5th=new TH1F(("h_jet_phi_5th_"+suffix).c_str(), "5th jet #phi; #phi; Events/GeV", 800, -4.0, 4.0); histCol.h_jet_phi_5th->Sumw2();
histCol.h_jet_phi_6th=new TH1F(("h_jet_phi_6th_"+suffix).c_str(), "6th jet #phi; #phi; Events/GeV", 800, -4.0, 4.0); histCol.h_jet_phi_6th->Sumw2();
histCol.h_jet_energy_leading=new TH1F(("h_jet_energy_leading_"+suffix).c_str(), "Leading jet Energy; Energy [GeV]; Events/GeV", 10000, 0, 1000); histCol.h_jet_energy_leading->Sumw2();
histCol.h_jet_energy_trailing=new TH1F(("h_jet_energy_trailing_"+suffix).c_str(), "Trailing jet Energy; Energy [GeV]; Events/GeV", 10000, 0, 1000); histCol.h_jet_energy_trailing->Sumw2();
histCol.h_jet_energy_3rd=new TH1F(("h_jet_energy_3rd_"+suffix).c_str(), "3rd jet Energy; Energy [GeV]; Events/GeV", 10000, 0, 1000); histCol.h_jet_energy_3rd->Sumw2();
histCol.h_jet_energy_4th=new TH1F(("h_jet_energy_4th_"+suffix).c_str(), "4th jet Energy; Energy [GeV]; Events/GeV", 10000, 0, 1000); histCol.h_jet_energy_4th->Sumw2();
histCol.h_jet_energy_5th=new TH1F(("h_jet_energy_5th_"+suffix).c_str(), "5th jet Energy; Energy [GeV]; Events/GeV", 10000, 0, 1000); histCol.h_jet_energy_5th->Sumw2();
histCol.h_jet_energy_6th=new TH1F(("h_jet_energy_6th_"+suffix).c_str(), "6th jet Energy; Energy [GeV]; Events/GeV", 10000, 0, 1000); histCol.h_jet_energy_6th->Sumw2();
histCol.h_bjet_pt_leading=new TH1F(("h_bjet_pt_leading_"+suffix).c_str(), "Leading b-jet pT; pT [GeV]; Events/GeV", 10000, 0, 1000); histCol.h_bjet_pt_leading->Sumw2();
histCol.h_bjet_eta_leading=new TH1F(("h_bjet_eta_leading_"+suffix).c_str(), "Leading b-jet #eta; #eta; Events", 600, -3.0, 3.0); histCol.h_bjet_eta_leading->Sumw2();
histCol.h_bjet_phi_leading=new TH1F(("h_bjet_phi_leading_"+suffix).c_str(), "Leading b-jet #phi; #phi; Events/GeV", 800, -4.0, 4.0); histCol.h_bjet_phi_leading->Sumw2();
histCol.h_bjet_energy_leading=new TH1F(("h_bjet_energy_leading_"+suffix).c_str(), "Leading b-jet Energy; Energy [GeV]; Events/GeV", 10000, 0, 1000); histCol.h_bjet_energy_leading->Sumw2();
histCol.h_bjet_pt_trailing=new TH1F(("h_bjet_pt_trailing_"+suffix).c_str(), "Trailing b-jet pT; pT [GeV]; Events/GeV", 10000, 0, 1000); histCol.h_bjet_pt_trailing->Sumw2();
histCol.h_bjet_eta_trailing=new TH1F(("h_bjet_eta_trailing_"+suffix).c_str(), "Trailing b-jet #eta; #eta; Events", 600, -3.0, 3.0); histCol.h_bjet_eta_trailing->Sumw2();
histCol.h_bjet_phi_trailing=new TH1F(("h_bjet_phi_trailing_"+suffix).c_str(), "Trailing b-jet #phi; #phi; Events/GeV", 800, -4.0, 4.0); histCol.h_bjet_phi_trailing->Sumw2();
histCol.h_bjet_energy_trailing=new TH1F(("h_bjet_energy_trailing_"+suffix).c_str(), "Trailing b-jet Energy; Energy [GeV]; Events/GeV", 10000, 0, 1000); histCol.h_bjet_energy_trailing->Sumw2();
histCol.h_photon_pt =new TH1F(("h_photon_pt_"+suffix).c_str(), "Photon pT; pT [GeV]; Events/GeV", 1000, 0, 1000); histCol.h_photon_pt->Sumw2();
histCol.h_photon_eta =new TH1F(("h_photon_eta_"+suffix).c_str(), "Photon #eta; #eta ; Events", 600, -3.0, 3.0); histCol.h_photon_eta->Sumw2();
histCol.h_photon_phi =new TH1F(("h_photon_phi_"+suffix).c_str(), "Photon #phi; #phi ; Events", 800, -4.0, 4.0); histCol.h_photon_phi->Sumw2();
histCol.h_photon_energy =new TH1F(("h_photon_energy_"+suffix).c_str(), "Photon Energy; Energy [GeV]; Events/GeV", 1000, 0, 1000); histCol.h_photon_energy->Sumw2();

}

void writeHistCollection(HistCollection &histCol)
{
  histCol.h_mu_pt_leading->Write();
  histCol.h_mu_pt_trailing->Write();
  histCol.h_mu_eta_leading->Write();
  histCol.h_mu_eta_trailing->Write();
  histCol.h_mu_phi_leading->Write();
  histCol.h_mu_phi_trailing->Write();
  histCol.h_mu_energy_leading->Write();
  histCol.h_mu_energy_trailing->Write();
  histCol.h_mu_dxy_leading->Write();
  histCol.h_mu_dxy_trailing->Write();
  histCol.h_mu_dz_leading->Write();
  histCol.h_mu_dz_trailing->Write();
  histCol.h_el_pt_leading->Write();
  histCol.h_el_pt_trailing->Write();
  histCol.h_el_eta_leading->Write();
  histCol.h_el_eta_trailing->Write();
  histCol.h_el_phi_leading->Write();
  histCol.h_el_phi_trailing->Write();
  histCol.h_el_energy_leading->Write();
  histCol.h_el_energy_trailing->Write();
  histCol.h_el_dxy_leading->Write();
  histCol.h_el_dxy_trailing->Write();
  histCol.h_el_dz_leading->Write();
  histCol.h_el_dz_trailing->Write();
  histCol.h_el_charge_leading->Write();
  histCol.h_el_charge_trailing->Write();
  histCol.h_mu_charge_leading->Write();
  histCol.h_mu_charge_trailing->Write();
  histCol.h_TransverseMass_leadingMu->Write();
  histCol.h_TransverseMass_trailingMu->Write();
  histCol.h_TransverseMass_leadingEl->Write();
  histCol.h_TransverseMass_trailingEl->Write();
  histCol.h_InvariantMass_Ph_leadingMu->Write();
  histCol.h_InvariantMass_Ph_trailingMu->Write();
  histCol.h_InvariantMass_Ph_leadingEl->Write();
  histCol.h_InvariantMass_Ph_trailingEl->Write();
  histCol.h_InvariantMass->Write();
  histCol.h_InvariantMass_Ph->Write();
  histCol.h_SVFitMass_Ph->Write();
  histCol.h_Invariant_Mass_Taus->Write();
  histCol.h_Invariant_Mass_Mus->Write();
  histCol.h_Invariant_Mass_Es->Write();
  histCol.h_SVFit->Write();
  histCol.h_SVFit_Mass_Taus->Write();
  histCol.h_DeltaPhi_met_mu1->Write();
  histCol.h_DeltaPhi_met_mu2->Write();
  histCol.h_DeltaPhi_ph_mu1->Write();
  histCol.h_DeltaPhi_ph_mu2->Write();
  histCol.h_DeltaPhi_met_el1->Write();
  histCol.h_DeltaPhi_met_el2->Write();
  histCol.h_DeltaPhi_ph_el1->Write();
  histCol.h_DeltaPhi_ph_el2->Write();
  histCol.h_Isolation_mu1->Write();
  histCol.h_Isolation_mu2->Write();
  histCol.h_Isolation_el1->Write();
  histCol.h_Isolation_el2->Write();
  histCol.h_nVertices->Write();
  histCol.h_caloMET->Write();
  histCol.h_MET->Write();
  histCol.h_MET_Signxx->Write();
  histCol.h_MET_Signxy->Write();
  histCol.h_MET_Signyx->Write();
  histCol.h_MET_Signyy->Write();
  histCol.h_HT->Write();
  histCol.h_nJets->Write();
  histCol.h_nbJets->Write();
  histCol.h_jet_pt_leading->Write();
  histCol.h_jet_pt_trailing->Write();
  histCol.h_jet_pt_3rd->Write();
  histCol.h_jet_pt_4th->Write();
  histCol.h_jet_pt_5th->Write();
  histCol.h_jet_pt_6th->Write();
  histCol.h_jet_eta_leading->Write();
  histCol.h_jet_eta_trailing->Write();
  histCol.h_jet_eta_3rd->Write();
  histCol.h_jet_eta_4th->Write();
  histCol.h_jet_eta_5th->Write();
  histCol.h_jet_eta_6th->Write();
  histCol.h_jet_phi_leading->Write();
  histCol.h_jet_phi_trailing->Write();
  histCol.h_jet_phi_3rd->Write();
  histCol.h_jet_phi_4th->Write();
  histCol.h_jet_phi_5th->Write();
  histCol.h_jet_phi_6th->Write();
  histCol.h_jet_energy_leading->Write();
  histCol.h_jet_energy_trailing->Write();
  histCol.h_jet_energy_3rd->Write();
  histCol.h_jet_energy_4th->Write();
  histCol.h_jet_energy_5th->Write();
  histCol.h_jet_energy_6th->Write();
  histCol.h_photon_pt->Write();
  histCol.h_photon_eta->Write();
  histCol.h_photon_phi->Write();
  histCol.h_photon_energy->Write();
  histCol.h_HTb->Write();
  histCol.h_bjet_pt_leading->Write();
  histCol.h_bjet_pt_trailing->Write();
  histCol.h_bjet_eta_leading->Write();
  histCol.h_bjet_eta_trailing->Write();
  histCol.h_bjet_phi_leading->Write();
  histCol.h_bjet_phi_trailing->Write();
  histCol.h_bjet_energy_leading->Write();
  histCol.h_bjet_energy_trailing->Write();
}

double photonSF(double phPt, double phEta)
{
  double sf = 1.0;
  if(phPt > 15.00 and phPt < 20.00)
    {
    if(fabs(phEta) > 0.0 and fabs(phEta) < 0.8) sf = 0.9496;
    else if(fabs(phEta) > 0.8 and fabs(phEta) < 1.4442) sf = 0.9803;
    else if(fabs(phEta) > 1.566 and fabs(phEta) < 2.0) sf = 1.0005;
    else if(fabs(phEta) > 2.0 and fabs(phEta) < 2.5) sf = 1.0171;
    }
  else if(phPt > 20.00 and phPt < 30.00)
    {
    if(fabs(phEta) > 0.0 and fabs(phEta) < 0.8) sf = 0.9672;
    else if(fabs(phEta) > 0.8 and fabs(phEta) < 1.4442) sf = 0.9724;
    else if(fabs(phEta) > 1.566 and fabs(phEta) < 2.0) sf = 0.9867;
    else if(fabs(phEta) > 2.0 and fabs(phEta) < 2.5) sf = 1.0130;
    }
  else if(phPt > 30.00 and phPt < 40.00)
    {
    if(fabs(phEta) > 0.0 and fabs(phEta) < 0.8) sf = 0.9711;
    else if(fabs(phEta) > 0.8 and fabs(phEta) < 1.4442) sf = 0.9688;
    else if(fabs(phEta) > 1.566 and fabs(phEta) < 2.0) sf = 0.9971;
    else if(fabs(phEta) > 2.0 and fabs(phEta) < 2.5) sf = 1.0143;
    }
  else if(phPt > 40.00 and phPt < 50.00)
    {
    if(fabs(phEta) > 0.0 and fabs(phEta) < 0.8) sf = 0.9766;
    else if(fabs(phEta) > 0.8 and fabs(phEta) < 1.4442) sf = 0.9805;
    else if(fabs(phEta) > 1.566 and fabs(phEta) < 2.0) sf = 0.9996;
    else if(fabs(phEta) > 2.0 and fabs(phEta) < 2.5) sf = 1.0129;
    }
  else if(phPt > 50.00)
    {
    if(fabs(phEta) > 0.0 and fabs(phEta) < 0.8) sf = 0.9815;
    else if(fabs(phEta) > 0.8 and fabs(phEta) < 1.4442) sf = 0.9837;
    else if(fabs(phEta) > 1.566 and fabs(phEta) < 2.0) sf = 1.0034;
    else if(fabs(phEta) > 2.0 and fabs(phEta) < 2.5) sf = 1.0128;
    }
  return sf;
}

double fastsimSF_Mu(double phPt, double phEta, std::string MCSample)
{
  double sf=1.0;
  if(MCSample=="Stop100"){
    if(phPt > 0.00 and phPt < 10.00)
      {
      if(fabs(phEta) > 0.0 and fabs(phEta) < 0.9) sf = 0.955;
      if(fabs(phEta) >= 0.9 and fabs(phEta) < 1.2) sf = 0.944;
      if(fabs(phEta) >= 1.2 and fabs(phEta) < 2.1) sf = 0.940;
      if(fabs(phEta) > 2.1) sf =  0.886;
      }
    else if(phPt >= 10.00 and phPt < 20.00)
      {
      if(fabs(phEta) > 0.0 and fabs(phEta) < 0.9) sf = 0.955;
      if(fabs(phEta) >= 0.9 and fabs(phEta) < 1.2) sf = 0.944;
      if(fabs(phEta) >= 1.2 and fabs(phEta) < 2.1) sf = 0.940;
      if(fabs(phEta) > 2.1) sf =  0.886;
      } 
    else if(phPt >= 20.00 and phPt < 30.00)
      {   
      if(fabs(phEta) > 0.0 and fabs(phEta) < 0.9) sf = 0.967;
      if(fabs(phEta) >= 0.9 and fabs(phEta) < 1.2) sf = 0.962;
      if(fabs(phEta) >= 1.2 and fabs(phEta) < 2.1) sf = 0.956;
      if(fabs(phEta) > 2.1) sf =  0.917;
      }
    else if(phPt >= 30.00 and phPt < 40.00)
      {
      if(fabs(phEta) > 0.0 and fabs(phEta) < 0.9) sf = 0.971;
      if(fabs(phEta) >= 0.9 and fabs(phEta) < 1.2) sf = 0.970;
      if(fabs(phEta) >= 1.2 and fabs(phEta) < 2.1) sf = 0.962;
      if(fabs(phEta) > 2.1) sf =  0.937;
      }
    else if(phPt >= 40.00 and phPt < 50.00)
      {
      if(fabs(phEta) > 0.0 and fabs(phEta) < 0.9) sf = 0.974;
      if(fabs(phEta) >= 0.9 and fabs(phEta) < 1.2) sf = 0.971;
      if(fabs(phEta) >= 1.2 and fabs(phEta) < 2.1) sf = 0.965;
      if(fabs(phEta) > 2.1) sf =  0.949;
      }
    else if(phPt >= 50.00 and phPt < 60.00)
      {
      if(fabs(phEta) > 0.0 and fabs(phEta) < 0.9) sf = 0.975;
      if(fabs(phEta) >= 0.9 and fabs(phEta) < 1.2) sf = 0.970;
      if(fabs(phEta) >= 1.2 and fabs(phEta) < 2.1) sf = 0.965;
      if(fabs(phEta) > 2.1) sf =  0.949;
      }
    else if(phPt >= 60.00 and phPt < 90.00)
      {
      if(fabs(phEta) > 0.0 and fabs(phEta) < 0.9) sf = 0.976;
      if(fabs(phEta) >= 0.9 and fabs(phEta) < 1.2) sf = 0.972;
      if(fabs(phEta) >= 1.2 and fabs(phEta) < 2.1) sf = 0.962;
      if(fabs(phEta) > 2.1) sf =  0.940;
      }
    else if(phPt >= 90.0)
      {
      if(fabs(phEta) > 0.0 and fabs(phEta) < 0.9) sf = 0.974;
      if(fabs(phEta) >= 0.9 and fabs(phEta) < 1.2) sf = 0.970;
      if(fabs(phEta) >= 1.2 and fabs(phEta) < 2.1) sf = 0.956;
      if(fabs(phEta) > 2.1) sf =  0.940; 
      }
    } 
  else sf = 1.0;
  return sf;
}

double fastsimSF_El(double phPt, double phEta, std::string MCSample)
{
  double sf=1.0;
  if(MCSample=="Stop100"){
    if(phPt > 0.00 and phPt < 10.00)
      {
      if(fabs(phEta) > 0.0 and fabs(phEta) < 0.5) sf = 0.903;
      if(fabs(phEta) >= 0.5 and fabs(phEta) < 1.0) sf = 0.929;
      if(fabs(phEta) >= 1.0 and fabs(phEta) < 1.5) sf = 0.923;
      if(fabs(phEta) >= 1.5 and fabs(phEta) < 2.0) sf = 0.726;
      if(fabs(phEta) >= 2.1) sf =  0.711;
      }
    else if(phPt >= 10.00 and phPt < 20.00)
      {
      if(fabs(phEta) > 0.0 and fabs(phEta) < 0.5) sf = 0.903;
      if(fabs(phEta) >= 0.5 and fabs(phEta) < 1.0) sf = 0.929;
      if(fabs(phEta) >= 1.0 and fabs(phEta) < 1.5) sf = 0.923;
      if(fabs(phEta) >= 1.5 and fabs(phEta) < 2.0) sf = 0.726;
      if(fabs(phEta) >= 2.1) sf =  0.711;
      } 
    else if(phPt >= 20.00 and phPt < 30.00)
      {   
      if(fabs(phEta) > 0.0 and fabs(phEta) < 0.5) sf = 0.949;
      if(fabs(phEta) >= 0.5 and fabs(phEta) < 1.0) sf = 0.949;
      if(fabs(phEta) >= 1.0 and fabs(phEta) < 1.5) sf = 0.954;
      if(fabs(phEta) >= 1.5 and fabs(phEta) < 2.0) sf = 0.803;
      if(fabs(phEta) >= 2.1) sf =  0.810;
      }
    else if(phPt >= 30.00 and phPt < 40.00)
      {
      if(fabs(phEta) > 0.0 and fabs(phEta) < 0.5) sf = 0.967;
      if(fabs(phEta) >= 0.5 and fabs(phEta) < 1.0) sf = 0.972;
      if(fabs(phEta) >= 1.0 and fabs(phEta) < 1.5) sf = 0.991;
      if(fabs(phEta) >= 1.5 and fabs(phEta) < 2.0) sf = 0.868;
      if(fabs(phEta) >= 2.1) sf =  0.832;
      }
    else if(phPt >= 40.00 and phPt < 50.00)
      {
      if(fabs(phEta) > 0.0 and fabs(phEta) < 0.5) sf = 0.974;
      if(fabs(phEta) >= 0.5 and fabs(phEta) < 1.0) sf = 0.979;
      if(fabs(phEta) >= 1.0 and fabs(phEta) < 1.5) sf = 1.021;
      if(fabs(phEta) >= 1.5 and fabs(phEta) < 2.0) sf = 0.896;
      if(fabs(phEta) >= 2.1) sf =  0.862;
      }
    else if(phPt >= 50.00 and phPt < 60.00)
      {
      if(fabs(phEta) > 0.0 and fabs(phEta) < 0.5) sf = 0.977;
      if(fabs(phEta) >= 0.5 and fabs(phEta) < 1.0) sf = 0.983;
      if(fabs(phEta) >= 1.0 and fabs(phEta) < 1.5) sf = 1.034;
      if(fabs(phEta) >= 1.5 and fabs(phEta) < 2.0) sf = 0.896;
      if(fabs(phEta) >= 2.1) sf =  0.864;
      }
    else if(phPt >= 60.0)
      {
      if(fabs(phEta) > 0.0 and fabs(phEta) < 0.5) sf = 0.972;
      if(fabs(phEta) >= 0.5 and fabs(phEta) < 1.0) sf = 0.977;
      if(fabs(phEta) >= 1.0 and fabs(phEta) < 1.5) sf = 0.998;
      if(fabs(phEta) >= 1.5 and fabs(phEta) < 2.0) sf = 0.893;
      if(fabs(phEta) >= 2.1) sf =  0.854;
      }
    } 
  else sf = 1.0;
  return sf;
}

/*
double photonSF_WG(double phPt, double phEta)
{
  double sf = 1.0;
  if(phPt > 30.00 and phPt < 40.00)
    {
    if(fabs(phEta) > 0.0 and fabs(phEta) < 0.8) sf = 0.9954*0.9962;
    if(fabs(phEta) > 0.8 and fabs(phEta) < 1.5) sf = 0.9954*0.9967;
    }
  if(phPt > 40.00 and phPt < 50.00)
    {
    if(fabs(phEta) > 0.0 and fabs(phEta) < 0.8) sf = 0.9707*0.9966;
    if(fabs(phEta) > 0.8 and fabs(phEta) < 1.5) sf = 0.9707*0.9970;
    }
  if(phPt > 50.0)
   {
   if(fabs(phEta) > 0.0 and fabs(phEta) < 0.8) sf = 1.0122*0.9960; 
   if(fabs(phEta) > 0.8 and fabs(phEta) < 1.5) sf = 1.0122*0.9971; 
   }
  return sf;
}
*/
double photonSF_WG(double phPt, double phEta)
{
  double sf = 1.0;
  if(phPt > 30.00 and phPt < 40.00)
    {
    if(fabs(phEta) > 0.0 and fabs(phEta) < 0.8) sf = 1.0170;
    else if(fabs(phEta) > 0.8 and fabs(phEta) < 1.5) sf = 0.9837; 
    }
  else if(phPt >= 40.00 and phPt < 45.00)
    {
    if(fabs(phEta) > 0.0 and fabs(phEta) < 0.8) sf = 1.0207;
    else if(fabs(phEta) > 0.8 and fabs(phEta) < 1.5) sf = 1.0231;
    }
  else if(phPt >= 45.00 and phPt < 50.00)
    {
    if(fabs(phEta) > 0.0 and fabs(phEta) < 0.8) sf = 1.0228;
    else if(fabs(phEta) > 0.8 and fabs(phEta) < 1.5) sf = 1.0533;
    }
  else if(phPt >= 50.00 and phPt < 55.00)
    {
    if(fabs(phEta) > 0.0 and fabs(phEta) < 0.8) sf = 1.0223;    
    else if(fabs(phEta) > 0.8 and fabs(phEta) < 1.5) sf = 1.0370;
    }
  else if(phPt >= 55.00 and phPt < 60.00)
    {
    if(fabs(phEta) > 0.0 and fabs(phEta) < 0.8) sf = 1.0200; 
    else if(fabs(phEta) > 0.8 and fabs(phEta) < 1.5) sf = 1.0103; 
    }
  else if(phPt >= 60.00 and phPt < 70.00)
    {
    if(fabs(phEta) > 0.0 and fabs(phEta) < 0.8) sf = 1.0072;
    else if(fabs(phEta) > 0.8 and fabs(phEta) < 1.5) sf = 1.0144;
    }
  else if(phPt >= 70.00 and phPt < 80.00) 
    {
    if(fabs(phEta) > 0.0 and fabs(phEta) < 0.8) sf = 1.0023; 
    else if(fabs(phEta) > 0.8 and fabs(phEta) < 1.5) sf = 1.0144;
    }
  else if(phPt >= 80.00)
    {
    if(fabs(phEta) > 0.0 and fabs(phEta) < 0.8) sf = 1.0024;  
    else if(fabs(phEta) > 0.8 and fabs(phEta) < 1.5) sf = 0.9440; 
    }
  return sf;
}  
double electronSF(double elecPt, double elecEta)
{
  double sf = 1.0;
  if(elecPt > 10.00 and elecPt < 15.00)
    {
    if(fabs(elecEta) > 0.0 and fabs(elecEta) < 0.8) sf = 0.838;
    else if(fabs(elecEta) > 0.8 and fabs(elecEta) < 1.44) sf = 0.861;
    else if(fabs(elecEta) > 1.44 and fabs(elecEta) < 1.56) sf = 1.021;
    else if(fabs(elecEta) > 1.56 and fabs(elecEta) < 2.00) sf = 0.951;
    else if(fabs(elecEta) > 2.00 and fabs(elecEta) < 2.50) sf = 1.055;
    }
  else if(elecPt > 15.00 and elecPt < 20.00)
    {
    if(fabs(elecEta) > 0.0 and fabs(elecEta) < 0.8) sf = 0.942;
    else if(fabs(elecEta) > 0.8 and fabs(elecEta) < 1.44) sf = 0.925;
    else if(fabs(elecEta) > 1.44 and fabs(elecEta) < 1.56) sf = 0.889;
    else if(fabs(elecEta) > 1.56 and fabs(elecEta) < 2.00) sf = 0.932;
    else if(fabs(elecEta) > 2.00 and fabs(elecEta) < 2.50) sf = 0.984;
    }
  else if(elecPt > 20.00 and elecPt < 30.00)
    {
    if(fabs(elecEta) > 0.0 and fabs(elecEta) < 0.8) sf = 0.980;
    else if(fabs(elecEta) > 0.8 and fabs(elecEta) < 1.44) sf = 0.955;
    else if(fabs(elecEta) > 1.44 and fabs(elecEta) < 1.56) sf = 0.996;
    else if(fabs(elecEta) > 1.56 and fabs(elecEta) < 2.00) sf = 0.970;
    else if(fabs(elecEta) > 2.00 and fabs(elecEta) < 2.50) sf = 1.027;
    }
  else if(elecPt > 30.00 and elecPt < 40.00)
    {
    if(fabs(elecEta) > 0.0 and fabs(elecEta) < 0.8) sf = 0.982;
    else if(fabs(elecEta) > 0.8 and fabs(elecEta) < 1.44) sf = 0.965;
    else if(fabs(elecEta) > 1.44 and fabs(elecEta) < 1.56) sf = 0.989;
    else if(fabs(elecEta) > 1.56 and fabs(elecEta) < 2.00) sf = 0.968;
    else if(fabs(elecEta) > 2.00 and fabs(elecEta) < 2.50) sf = 1.018;
    }
  else if(elecPt > 40.00 and elecPt < 50.00)
    {
    if(fabs(elecEta) > 0.0 and fabs(elecEta) < 0.8) sf = 0.985;
    else if(fabs(elecEta) > 0.8 and fabs(elecEta) < 1.44) sf = 0.974;
    else if(fabs(elecEta) > 1.44 and fabs(elecEta) < 1.56) sf = 0.961;
    else if(fabs(elecEta) > 1.56 and fabs(elecEta) < 2.00) sf = 0.989;
    else if(fabs(elecEta) > 2.00 and fabs(elecEta) < 2.50) sf = 1.012;
    }
  else if(elecPt > 50.00 and elecPt < 200.00)
    {
    if(fabs(elecEta) > 0.0 and fabs(elecEta) < 0.8) sf = 0.984;
    else if(fabs(elecEta) > 0.8 and fabs(elecEta) < 1.44) sf = 0.978;
    else if(fabs(elecEta) > 1.44 and fabs(elecEta) < 1.56) sf = 0.982;
    else if(fabs(elecEta) > 1.56 and fabs(elecEta) < 2.00) sf = 0.991;
    else if(fabs(elecEta) > 2.00 and fabs(elecEta) < 2.50) sf = 1.008;
    }
  return sf;
}

double electronSF_WG(double elecPt, double elecEta)
{
  double sf = 1.0;
  TFile* fileElectronSF = TFile::Open("/uscms_data/d2/lpcljm/sapta/SUSYSearch/CMSSW_5_3_16_patch1/src/MPAnalyzer/CombinedMethod_ScaleFactors_IdIsoSip.root");  
  TH2F *h_ElectronSF= (TH2F*) fileElectronSF->Get("h_electronScaleFactor_IdIsoSip");
  int binx = h_ElectronSF->GetXaxis()->FindBin(elecPt);
  int biny = h_ElectronSF->GetYaxis()->FindBin(elecEta);
  sf = h_ElectronSF->GetBinContent(binx, biny);
  //int binxa = h_ElectronSF->GetXaxis()->FindBin(200);
  if(elecPt>=200) sf = 1.0;
  //if(elecPt>=200) sf = h_ElectronSF->GetBinContent(binxa, biny); 
  fileElectronSF->Close();
  return sf;
}

double muonSF(double muPt, double muEta)
{
  double sf = 1.0;
  if(muPt > 10.00 and muPt < 1000.00)
  {
    if(fabs(muEta) > 0.0 and fabs(muEta) < 0.90) sf = 0.9943;
    else if(fabs(muEta) > 0.9 and fabs(muEta) < 1.20) sf = 0.9933;
    else if(fabs(muEta) > 1.20 and fabs(muEta) < 2.50) sf = 1.0020;
  }
  return sf;
}
/*
double muonSF_WG(double muPt, double muEta)//Loose*Soft SF
{
  double sf = 1.0;
  if(muPt > 3.00 and muPt < 20.00)
  {
    if(fabs(muEta) > 0.0 and fabs(muEta) < 0.90) sf = 1.0002*1.0157;
    else if(fabs(muEta) > 0.9 and fabs(muEta) < 1.20) sf = 0.9824*0.9982;
    else if(fabs(muEta) > 1.20 and fabs(muEta) < 2.50) sf = 1.0075*1.0277;
  }
  if(muPt > 20.00)
  {
    if(fabs(muEta) > 0.0 and fabs(muEta) < 0.90) sf = 0.9984*0.9850;
    else if(fabs(muEta) > 0.9 and fabs(muEta) < 1.20) sf = 0.9990*0.9928;
    else if(fabs(muEta) > 1.2 and fabs(muEta) < 2.10) sf = 0.9986*1.0007;
    else if(fabs(muEta) > 2.1 and fabs(muEta) < 2.50) sf = 1.0000*1.0365;
  }
  return sf;
}*/

double muonSF_WG(double muPt, double muEta)//Soft SF
{
  double sf = 1.0;
  if(muPt > 3.00 and muPt < 20.00)
  {
    if(fabs(muEta) > 0.0 and fabs(muEta) < 0.90) sf = 1.0157;
    else if(fabs(muEta) > 0.9 and fabs(muEta) < 1.20) sf = 0.9982;
    else if(fabs(muEta) > 1.20 and fabs(muEta) < 2.50) sf = 1.0277;
  }
  if(muPt > 20.00)
  {
    if(fabs(muEta) > 0.0 and fabs(muEta) < 0.90) sf = 0.9850;
    else if(fabs(muEta) > 0.9 and fabs(muEta) < 1.20) sf = 0.9928;
    else if(fabs(muEta) > 1.2 and fabs(muEta) < 2.10) sf = 1.0007;
    else if(fabs(muEta) > 2.1 and fabs(muEta) < 2.50) sf = 1.0365;
  }
  return sf;
}

bool sortGenParticlesInDescendingpT(GenParticleInfo gen1, GenParticleInfo gen2)
{
  return (sqrt(gen1.Px*gen1.Px + gen1.Py*gen1.Py) > sqrt(gen2.Px*gen2.Px + gen2.Py*gen2.Py));
}

bool sortLeptonsInDescendingpT(LeptonInfo lep1, LeptonInfo lep2)
{
  return (lep1.pT > lep2.pT);
}

bool sortJetsInDescendingpT(JetInfo jet1, JetInfo jet2)
{
  return (jet1.pT > jet2.pT);
}

bool sortJetVectorsInDescendingpT(AnalysisJetInfo jet1, AnalysisJetInfo jet2)
{
  return (jet1.JetLV.Pt() > jet2.JetLV.Pt());
}

bool sortBJetVectorsInDescendingpT(BJetInfo bjet1, BJetInfo bjet2)
{
  return (bjet1.BJetLV.Pt() > bjet2.BJetLV.Pt()); 
}

bool sortPhotonsInDescendingpT(PhotonInfo pho1, PhotonInfo pho2)
{
  return (pho1.pT > pho2.pT);
}

bool sortPhotonVectorsInDescendingpT(TLorentzVector pho1, TLorentzVector pho2)
{
  return (pho1.Pt() > pho2.Pt());
}

bool sortMatchedLeptonsInDescendingpT(MatchedLeptonInfo mlep1, MatchedLeptonInfo mlep2)
{
  return (mlep1.pT > mlep2.pT);
}

bool sortVectorsInDescendingpT(AnalysisLeptonInfo lep1, AnalysisLeptonInfo lep2)
{
  return (lep1.LepLV.Pt() > lep2.LepLV.Pt());
}

bool returnAnaPhoton(PhotonInfo pho1) 
{

  if(pho1.pT > 30.0 and
     fabs(pho1.eta) < 1.4442 and
     pho1.phR9 > 0.9 and 
     pho1.phR9 < 1.0 and 
     pho1.phHoE < 0.05 and
     pho1.phSigmaIetaIeta < 0.012 and 
     pho1.phpixelVeto==0 and
     pho1.chIsolation < 2.6 and
     pho1.nuIsolation < 1.0 + 0.04*pho1.pT and
     pho1.phIsolation < 1.3 + 0.005*pho1.pT and
     pho1.phSigmaIetaIeta > 0.001 and
     std::sqrt(pho1.phSigmaIphiIphi) > 0.009) return true;
  else return false;

}
