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
  float isolation;
} LeptonInfo;

typedef struct
{
  TLorentzVector LepLV;
  float Isolation;
  int Charge;
} AnalysisLeptonInfo;

typedef struct
{
  float pT;
  float eta;
  float phi;
  float energy;
  int charge;
  bool isTight;
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
} JetInfo;


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
  TH1F *h_mu_pt_leading;
  TH1F *h_mu_pt_trailing; 
  TH1F *h_mu_eta_leading;
  TH1F *h_mu_eta_trailing;
  TH1F *h_mu_phi_leading;
  TH1F *h_mu_phi_trailing;
  TH1F *h_mu_energy_leading;
  TH1F *h_mu_energy_trailing;
  TH1F *h_el_pt_leading;    
  TH1F *h_el_pt_trailing;
  TH1F *h_el_eta_leading;
  TH1F *h_el_eta_trailing;
  TH1F *h_el_phi_leading;
  TH1F *h_el_phi_trailing;
  TH1F *h_el_energy_leading;
  TH1F *h_el_energy_trailing; 
  TH1F *h_InvariantMass;
  TH1F *h_InvariantMass_Ph;
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
  TH1F *h_HT;
  TH1F *h_nJets;
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
histCol.h_Isolation_mu1 = new TH1F(("h_Isolation_mu1_"+suffix).c_str(),"Leading muon isolation; Isolation; Events", 10000, 0, 0.15);histCol.h_Isolation_mu1->Sumw2();
histCol.h_Isolation_mu2 = new TH1F(("h_Isolation_mu2_"+suffix).c_str(),"Trailing muon isolation; Isolation; Events", 10000, 0, 0.15);histCol.h_Isolation_mu2->Sumw2();
histCol.h_mu_energy_trailing = new TH1F(("h_mu_energy_trailing_"+suffix).c_str(), "Trailing muon Energy; Energy [GeV]; Events/GeV", 1000, 0, 1000); histCol.h_mu_energy_trailing->Sumw2();
histCol.h_el_pt_leading = new TH1F(("h_el_pt_leading_"+suffix).c_str(), "Leading electron pT; pT [GeV]; Events/GeV", 1000, 0, 1000); histCol.h_el_pt_leading->Sumw2();
histCol.h_el_pt_trailing = new TH1F(("h_el_pt_trailing_"+suffix).c_str(), "Trailing electron pT; pT [GeV]; Events/GeV", 1000, 0, 1000); histCol.h_el_pt_trailing->Sumw2();
histCol.h_el_eta_leading = new TH1F(("h_el_eta_leading_"+suffix).c_str(), "Leading electron #eta ; #eta ; Events", 600, -3.0, 3.0); histCol.h_el_eta_leading->Sumw2();
histCol.h_el_eta_trailing = new TH1F(("h_el_eta_trailing_"+suffix).c_str(), "Trailing electron #eta ; #eta ; Events", 600, -3.0, 3.0); histCol.h_el_eta_trailing->Sumw2();
histCol.h_el_phi_leading = new TH1F(("h_el_phi_leading_"+suffix).c_str(), "Leading electron #phi ; #phi ; Events", 800, -4.0, 4.0); histCol.h_el_phi_leading->Sumw2();
histCol.h_el_phi_trailing = new TH1F(("h_el_phi_trailing_"+suffix).c_str(), "Trailing electron #phi ; #phi ; Events", 800, -4.0, 4.0); histCol.h_el_phi_trailing->Sumw2();
histCol.h_el_energy_leading = new TH1F(("h_el_energy_leading_"+suffix).c_str(), "Leading electron Energy; Energy [GeV]; Events/GeV", 1000, 0, 1000); histCol.h_el_energy_leading->Sumw2();
histCol.h_el_energy_trailing = new TH1F(("h_el_energy_trailing_"+suffix).c_str(), "Trailing electron Energy; Energy [GeV]; Events/GeV", 1000, 0, 1000); histCol.h_el_energy_trailing->Sumw2();
histCol.h_Isolation_el1 = new TH1F(("h_Isolation_el1_"+suffix).c_str(),"Leading electron isolation; Isolation; Events", 10000, 0, 0.15);histCol.h_Isolation_el1->Sumw2();
histCol.h_Isolation_el2 = new TH1F(("h_Isolation_el2_"+suffix).c_str(),"Trailing electron isolation; Isolation; Events", 10000, 0, 0.15);histCol.h_Isolation_el2->Sumw2();
histCol.h_InvariantMass=new TH1F(("h_InvariantMass_"+suffix).c_str(), "Di-lepton invariant mass; m_{ll} [GeV]; Events/GeV", 9000, 0, 300); histCol.h_InvariantMass->Sumw2();
histCol.h_InvariantMass_Ph=new TH1F(("h_InvariantMass_Ph_"+suffix).c_str(), "Di-lepton and photon invariant mass; m_{ll#gamma} [GeV]; Events/GeV", 9000, 0, 300); histCol.h_InvariantMass_Ph->Sumw2();
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
histCol.h_HT = new TH1F(("h_HT_"+suffix).c_str(), "HT (scalar sum of jet pT); H_T [GeV]; Events/GeV", 5000, 0, 5000.0);histCol.h_HT->Sumw2();
histCol.h_nJets = new TH1F(("h_nJets_"+suffix).c_str(), "Number of Jets; Number of Jets; Events", 20, -0.5, 19.5);histCol.h_nJets->Sumw2();
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
  histCol.h_el_pt_leading->Write();
  histCol.h_el_pt_trailing->Write();
  histCol.h_el_eta_leading->Write();
  histCol.h_el_eta_trailing->Write();
  histCol.h_el_phi_leading->Write();
  histCol.h_el_phi_trailing->Write();
  histCol.h_el_energy_leading->Write();
  histCol.h_el_energy_trailing->Write();
  histCol.h_InvariantMass->Write();
  histCol.h_InvariantMass_Ph->Write();
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
  histCol.h_HT->Write();
  histCol.h_nJets->Write();
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

bool sortLeptonsInDescendingpT(LeptonInfo lep1, LeptonInfo lep2)
{
  return (lep1.pT > lep2.pT);
}

bool sortJetsInDescendingpT(JetInfo jet1, JetInfo jet2)
{
  return (jet1.pT > jet2.pT);
}

bool sortJetVectorsInDescendingpT(TLorentzVector jet1, TLorentzVector jet2)
{
  return (jet1.Pt() > jet2.Pt());
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
