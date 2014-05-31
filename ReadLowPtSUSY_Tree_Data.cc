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

TLorentzVector fillTLorentzVector(double pT, double eta, double phi, double E)
{
  TLorentzVector object_p4;
  object_p4.SetPtEtaPhiE(pT, eta, phi, E);
  return object_p4;
}

double mdeltaR(double eta1, double phi1, double eta2, double phi2) {
  double delta_eta = fabs(eta1-eta2);
  double delta_phi = fabs(phi1-phi2);
  if(delta_phi > 3.14159265) delta_phi = delta_phi - 2*3.14159265;
  return std::sqrt((delta_eta)*(delta_eta) + (delta_phi)*(delta_phi));
}

typedef struct
{
  float pT;
  float eta;
  float phi;
  float energy;
  int charge;
  bool isTight;
  float isolation;
} LeptonInfo;

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



bool sortLeptonsInDescendingpT(LeptonInfo lep1, LeptonInfo lep2)
{
  return (lep1.pT > lep2.pT);
}

bool sortJetsInDescendingpT(JetInfo jet1, JetInfo jet2)
{
  return (jet1.pT > jet2.pT);
}

bool sortPhotonsInDescendingpT(PhotonInfo pho1, PhotonInfo pho2)
{
  return (pho1.pT > pho2.pT);
}

int ReadLowPtSUSY_Tree_Data(std::string infile, std::string outfile){
  
  std::string inputfilename=(infile+".root").c_str();
  TChain *tree=new TChain("LowPtSUSY_Tree");
  tree->Add(inputfilename.c_str());
  std::cout<<"Opened input file "<<inputfilename<<std::endl;

  Int_t           run;
  Int_t           lumi;
  Int_t           event;
  vector<float>   *ph_pt;
  vector<float>   *ph_phi;
  vector<float>   *ph_eta;
  vector<float>   *ph_energy;
  Int_t           nPhotons;
  vector<float>   *el_pt;
  vector<float>   *el_phi;
  vector<float>   *el_eta;
  vector<float>   *el_energy;
  vector<int>     *el_charge;
  vector<bool>    *el_isTight;
  Int_t           nElectrons;
  vector<float>   *mu_pt;
  vector<float>   *mu_phi;
  vector<float>   *mu_eta;
  vector<float>   *mu_energy;
  vector<int>     *mu_charge;
  vector<bool>    *mu_isTight;
  Int_t           nMuons;
  vector<float>   *jet_pt;
  vector<float>   *jet_phi;
  vector<float>   *jet_eta;
  vector<float>   *jet_energy;
  vector<int>     *jet_mva_loose;
  vector<int>     *jet_mva_tight;
  vector<int>     *jet_mva_medium;
  vector<int>     *jet_cut_loose;
  vector<int>     *jet_cut_tight;
  vector<int>     *jet_cut_medium;
  Int_t           nJets;
  Float_t         MET;
  Float_t         MET_Phi;
  Float_t         MET_Px;
  Float_t         MET_Py;
  Bool_t          fired_HLTPho;
  Bool_t          fired_HLTPhoId;
  Bool_t          fired_HLTPhoIdMet;
  Int_t           nVertices;
  vector<float>   *el_iso;
  vector<float>   *mu_iso;
  vector<float>   *ph_chIso;
  vector<float>   *ph_nuIso;
  vector<float>   *ph_phIso;
  vector<bool>    *ph_isTight;
  vector<bool>    *ph_phIsoTight;
  vector<bool>    *ph_phIsoMedium;
  vector<bool>    *ph_phIsoLoose;
  Float_t         nPUVertices;
  Float_t         nPUVerticesTrue;
  Float_t         PUWeightData;
  Float_t         PUWeightDataSys;
  float trigObj1Px;
  float trigObj1Py;
  float trigObj1Pz;
  float trigObj1E;
  float trigObj2Px;
  float trigObj2Py;
  float trigObj2Pz;
  float trigObj2E;


  // Set object pointer
  ph_pt = 0;
  ph_phi = 0;
  ph_eta = 0;
  ph_energy = 0;
  el_pt = 0;
  el_phi = 0;
  el_eta = 0;
  el_energy = 0;
  el_charge = 0;
  el_isTight = 0;
  mu_pt = 0;
  mu_phi = 0;
  mu_eta = 0;
  mu_energy = 0;
  mu_charge = 0;
  mu_isTight = 0;
  jet_pt = 0;
  jet_phi = 0;
  jet_eta = 0;
  jet_energy = 0;
  jet_mva_loose = 0;
  jet_mva_tight = 0;
  jet_mva_medium = 0;
  jet_cut_loose = 0;
  jet_cut_tight = 0;
  jet_cut_medium = 0;
  el_iso = 0;
  mu_iso = 0;
  ph_chIso = 0;
  ph_nuIso = 0;
  ph_phIso = 0;
  ph_isTight = 0;
  ph_phIsoTight = 0;
  ph_phIsoMedium = 0;
  ph_phIsoLoose = 0;
  trigObj1Px = 0;
  trigObj1Py = 0;
  trigObj1Pz = 0;
  trigObj1E = 0;
  trigObj2Px = 0;
  trigObj2Py = 0;
  trigObj2Pz = 0;
  trigObj2E = 0;

  tree->SetBranchAddress("run", &(run));
  tree->SetBranchAddress("lumi", &(lumi));
  tree->SetBranchAddress("event", &(event));
  tree->SetBranchAddress("ph_pt", &(ph_pt));
  tree->SetBranchAddress("ph_phi", &(ph_phi));
  tree->SetBranchAddress("ph_eta", &(ph_eta));
  tree->SetBranchAddress("ph_energy", &(ph_energy));
  tree->SetBranchAddress("ph_chIso", &(ph_chIso));
  tree->SetBranchAddress("ph_nuIso", &(ph_nuIso));
  tree->SetBranchAddress("ph_phIso", &(ph_phIso));
  tree->SetBranchAddress("ph_isTight", &(ph_isTight));
  tree->SetBranchAddress("nPhotons", &(nPhotons));
  tree->SetBranchAddress("el_pt", &(el_pt));
  tree->SetBranchAddress("el_eta", &(el_eta));
  tree->SetBranchAddress("el_phi", &(el_phi));
  tree->SetBranchAddress("el_energy", &(el_energy));
  tree->SetBranchAddress("el_charge", &(el_charge));
  tree->SetBranchAddress("el_isTight", &(el_isTight));  
  tree->SetBranchAddress("el_iso", &(el_iso));
  tree->SetBranchAddress("nElectrons", &(nElectrons));
  tree->SetBranchAddress("mu_pt", &(mu_pt));
  tree->SetBranchAddress("mu_eta", &(mu_eta));
  tree->SetBranchAddress("mu_phi", &(mu_phi));
  tree->SetBranchAddress("mu_energy", &(mu_energy)); 
  tree->SetBranchAddress("mu_charge", &(mu_charge)); 
  tree->SetBranchAddress("mu_isTight", &(mu_isTight)); 
  tree->SetBranchAddress("mu_iso", &(mu_iso));
  tree->SetBranchAddress("nMuons", &(nMuons));
  tree->SetBranchAddress("jet_pt", &(jet_pt));
  tree->SetBranchAddress("jet_phi", &(jet_phi));
  tree->SetBranchAddress("jet_eta", &(jet_eta));
  tree->SetBranchAddress("jet_energy", &(jet_energy));
  tree->SetBranchAddress("jet_mva_loose", &(jet_mva_loose));
  tree->SetBranchAddress("jet_mva_tight", &(jet_mva_tight));
  tree->SetBranchAddress("jet_mva_medium", &(jet_mva_medium));
  tree->SetBranchAddress("jet_cut_loose", &(jet_cut_loose));
  tree->SetBranchAddress("jet_cut_tight", &(jet_cut_tight));
  tree->SetBranchAddress("jet_cut_medium", &(jet_cut_medium));
  tree->SetBranchAddress("nJets", &(nJets));
  tree->SetBranchAddress("MET", &(MET));
  tree->SetBranchAddress("MET_Phi", &(MET_Phi));
  tree->SetBranchAddress("MET_Px", &(MET_Px));
  tree->SetBranchAddress("MET_Py", &(MET_Py));
  tree->SetBranchAddress("fired_HLTPho", &(fired_HLTPho));
  tree->SetBranchAddress("fired_HLTPhoId", &(fired_HLTPhoId));
  tree->SetBranchAddress("fired_HLTPhoIdMet", &(fired_HLTPhoIdMet));
  tree->SetBranchAddress("nVertices", &(nVertices));
  tree->SetBranchAddress("nPUVertices", &(nPUVertices));
  tree->SetBranchAddress("nPUVerticesTrue", &(nPUVerticesTrue));
  tree->SetBranchAddress("PUWeightData", &(PUWeightData));
  tree->SetBranchAddress("PUWeightDataSys", &(PUWeightDataSys)); 
  tree->SetBranchAddress("ph_phIsoTight", &(ph_phIsoTight));
  tree->SetBranchAddress("ph_phIsoMedium", &(ph_phIsoMedium));
  tree->SetBranchAddress("ph_phIsoLoose", &(ph_phIsoLoose));
  tree->SetBranchAddress("trigObj1Px", &(trigObj1Px));
  tree->SetBranchAddress("trigObj1Py", &(trigObj1Py));
  tree->SetBranchAddress("trigObj1Pz", &(trigObj1Pz));
  tree->SetBranchAddress("trigObj1E", &(trigObj1E));
  tree->SetBranchAddress("trigObj2Px", &(trigObj2Px));
  tree->SetBranchAddress("trigObj2Py", &(trigObj2Py));
  tree->SetBranchAddress("trigObj2Pz", &(trigObj2Pz));
  tree->SetBranchAddress("trigObj2E", &(trigObj2E));

  //Booking histograms:
  
  TH1F *h_jet_pt_leading=new TH1F("h_jet_pt_leading", "Leading jet pT; pT [GeV]; Events/GeV", 10000, 0, 1000); h_jet_pt_leading->Sumw2();
  TH1F *h_jet_pt_trailing=new TH1F("h_jet_pt_trailing", "Trailing jet pT; pT [GeV]; Events/GeV", 10000, 0, 1000); h_jet_pt_trailing->Sumw2();
  TH1F *h_jet_pt_3rd=new TH1F("h_jet_pt_3rd", "3rd jet pT; pT [GeV]; Events/GeV", 10000, 0, 1000); h_jet_pt_3rd->Sumw2();
  TH1F *h_jet_pt_4th=new TH1F("h_jet_pt_4th", "4th jet pT; pT [GeV]; Events/GeV", 10000, 0, 1000); h_jet_pt_4th->Sumw2();
  TH1F *h_jet_pt_5th=new TH1F("h_jet_pt_5th", "5th jet pT; pT [GeV]; Events/GeV", 10000, 0, 1000); h_jet_pt_5th->Sumw2();
  TH1F *h_jet_pt_6th=new TH1F("h_jet_pt_6th", "6th jet pT; pT [GeV]; Events/GeV", 10000, 0, 1000); h_jet_pt_6th->Sumw2();
  
  TH1F *h_jet_eta_leading=new TH1F("h_jet_eta_leading", "Leading jet #eta; #eta; Events", 600, -3.0, 3.0); h_jet_eta_leading->Sumw2();
  TH1F *h_jet_eta_trailing=new TH1F("h_jet_eta_trailing", "Trailing jet #eta; #eta; Events", 600, -3.0, 3.0); h_jet_eta_trailing->Sumw2();
  TH1F *h_jet_eta_3rd=new TH1F("h_jet_eta_3rd", "3rd jet #eta; #eta; Events/GeV", 600, -3.0, 3.0); h_jet_eta_3rd->Sumw2();
  TH1F *h_jet_eta_4th=new TH1F("h_jet_eta_4th", "4th jet #eta; #eta; Events/GeV", 600, -3.0, 3.0); h_jet_eta_4th->Sumw2();
  TH1F *h_jet_eta_5th=new TH1F("h_jet_eta_5th", "5th jet #eta; #eta; Events/GeV", 600, -3.0, 3.0); h_jet_eta_5th->Sumw2();
  TH1F *h_jet_eta_6th=new TH1F("h_jet_eta_6th", "6th jet #eta; #eta; Events/GeV", 600, -3.0, 3.0); h_jet_eta_6th->Sumw2();  
 
  TH1F *h_jet_phi_leading=new TH1F("h_jet_phi_leading", "Leading jet #phi; #phi; Events/GeV", 800, -4.0, 4.0); h_jet_phi_leading->Sumw2();
  TH1F *h_jet_phi_trailing=new TH1F("h_jet_phi_trailing", "Trailing jet #phi; #phi; Events/GeV", 800, -4.0, 4.0); h_jet_phi_trailing->Sumw2();
  TH1F *h_jet_phi_3rd=new TH1F("h_jet_phi_3rd", "3rd jet #phi; #phi; Events/GeV", 800, -4.0, 4.0); h_jet_phi_3rd->Sumw2();
  TH1F *h_jet_phi_4th=new TH1F("h_jet_phi_4th", "4th jet #phi; #phi; Events/GeV", 800, -4.0, 4.0); h_jet_phi_4th->Sumw2();
  TH1F *h_jet_phi_5th=new TH1F("h_jet_phi_5th", "5th jet #phi; #phi; Events/GeV", 800, -4.0, 4.0); h_jet_phi_5th->Sumw2();
  TH1F *h_jet_phi_6th=new TH1F("h_jet_phi_6th", "6th jet #phi; #phi; Events/GeV", 800, -4.0, 4.0); h_jet_phi_6th->Sumw2();

  TH1F *h_jet_energy_leading=new TH1F("h_jet_energy_leading", "Leading jet Energy; Energy [GeV]; Events/GeV", 10000, 0, 1000); h_jet_energy_leading->Sumw2();
  TH1F *h_jet_energy_trailing=new TH1F("h_jet_energy_trailing", "Trailing jet Energy; Energy [GeV]; Events/GeV", 10000, 0, 1000); h_jet_energy_trailing->Sumw2();
  TH1F *h_jet_energy_3rd=new TH1F("h_jet_energy_3rd", "3rd jet Energy; Energy [GeV]; Events/GeV", 10000, 0, 1000); h_jet_energy_3rd->Sumw2();
  TH1F *h_jet_energy_4th=new TH1F("h_jet_energy_4th", "4th jet Energy; Energy [GeV]; Events/GeV", 10000, 0, 1000); h_jet_energy_4th->Sumw2();
  TH1F *h_jet_energy_5th=new TH1F("h_jet_energy_5th", "5th jet Energy; Energy [GeV]; Events/GeV", 10000, 0, 1000); h_jet_energy_5th->Sumw2();
  TH1F *h_jet_energy_6th=new TH1F("h_jet_energy_6th", "6th jet Energy; Energy [GeV]; Events/GeV", 10000, 0, 1000); h_jet_energy_6th->Sumw2();

  TH1F *h_mu_pt_leading=new TH1F("h_mu_pt_leading", "Leading muon pT; pT [GeV]; Events/GeV", 1000, 0, 1000); h_mu_pt_leading->Sumw2();
  TH1F *h_mu_pt_trailing=new TH1F("h_mu_pt_trailing", "Trailing muon pT; pT [GeV]; Events/GeV", 1000, 0, 1000); h_mu_pt_trailing->Sumw2();
  TH1F *h_el_pt_leading=new TH1F("h_el_pt_leading", "Leading electron pT; pT [GeV]; Events/GeV", 1000, 0, 1000); h_el_pt_leading->Sumw2();
  TH1F *h_el_pt_trailing=new TH1F("h_el_pt_trailing", "Trailing electron pT; pT [GeV]; Events/GeV", 1000, 0, 1000); h_el_pt_trailing->Sumw2();
  TH1F *h_mu_eta_leading=new TH1F("h_mu_eta_leading", "Leading muon #eta ; #eta ; Events", 600, -3.0, 3.0); h_mu_eta_leading->Sumw2();
  TH1F *h_mu_eta_trailing=new TH1F("h_mu_eta_trailing", "Trailing muon #eta; #eta ; Events", 600, -3.0, 3.0); h_mu_eta_trailing->Sumw2();
  TH1F *h_el_eta_leading=new TH1F("h_el_eta_leading", "Leading electron #eta; #eta ; Events", 600, -3.0, 3.0); h_el_eta_leading->Sumw2();
  TH1F *h_el_eta_trailing=new TH1F("h_el_eta_trailing", "Trailing electron #eta; #eta ; Events", 600.0, -3.0, 3.0); h_el_eta_trailing->Sumw2();
  TH1F *h_mu_phi_leading=new TH1F("h_mu_phi_leading", "Leading muon #phi ; #phi ; Events", 800, -4.0, 4.0); h_mu_phi_leading->Sumw2();
  TH1F *h_mu_phi_trailing=new TH1F("h_mu_phi_trailing", "Trailing muon #phi; #phi ; Events", 800, -4.0, 4.0); h_mu_phi_trailing->Sumw2();
  TH1F *h_el_phi_leading=new TH1F("h_el_phi_leading", "Leading electron #phi; #phi ; Events", 800, -4.0, 4.0); h_el_phi_leading->Sumw2();
  TH1F *h_el_phi_trailing=new TH1F("h_el_phi_trailing", "Trailing electron #phi; #phi ; Events", 800.0, -4.0, 4.0); h_el_phi_trailing->Sumw2();
  TH1F *h_mu_energy_leading=new TH1F("h_mu_energy_leading", "Leading muon Energy; Energy [GeV]; Events/GeV", 1000, 0, 1000); h_mu_energy_leading->Sumw2();
  TH1F *h_mu_energy_trailing=new TH1F("h_mu_energy_trailing", "Trailing muon Energy; Energy [GeV]; Events/GeV", 1000, 0, 1000); h_mu_energy_trailing->Sumw2();
  TH1F *h_el_energy_leading=new TH1F("h_el_energy_leading", "Leading electron Energy; Energy [GeV]; Events/GeV", 1000, 0, 1000); h_el_energy_leading->Sumw2();
  TH1F *h_el_energy_trailing=new TH1F("h_el_energy_trailing", "Trailing electron Energy; Energy [GeV]; Events/GeV", 1000, 0, 1000); h_el_energy_trailing->Sumw2();

  TH1F *h_mu_isolation_leading=new TH1F("h_mu_isolation_leading", "Leading muon Isolation; Isolation; Events", 1000, 0, 1); h_mu_isolation_leading->Sumw2();
  TH1F *h_mu_isolation_trailing=new TH1F("h_mu_isolation_trailing", "Trailing muon Isolation; Isolation; Events", 1000, 0, 1); h_mu_isolation_trailing->Sumw2();
  TH1F *h_el_isolation_leading=new TH1F("h_el_isolation_leading", "Leading electron Isolation; Isolation; Events", 1000, 0, 1); h_el_isolation_leading->Sumw2();
  TH1F *h_el_isolation_trailing=new TH1F("h_el_isolation_trailing", "Trailing electron Isolation; Isolation; Events", 1000, 0, 1); h_el_isolation_trailing->Sumw2();

  TH1F *h_ph_chIsolation_leading=new TH1F("h_ph_chIsolation_leading", "Leading Photon Ch Isolation; Isolation; Events", 20000, -100, 100); h_ph_chIsolation_leading->Sumw2();
  TH1F *h_ph_nuIsolation_leading=new TH1F("h_ph_nuIsolation_leading", "Leading Photon Nu Isolation; Isolation; Events", 20000, -100, 100); h_ph_nuIsolation_leading->Sumw2();
  TH1F *h_ph_phIsolation_leading=new TH1F("h_ph_phIsolation_leading", "Leading Photon Ph Isolation; Isolation; Events", 20000, -100, 100); h_ph_phIsolation_leading->Sumw2();

  TH1F *h_photon_pt =new TH1F("h_photon_pt", "Photon pT; pT [GeV]; Events/GeV", 1000, 0, 1000); h_photon_pt->Sumw2();
  TH1F *h_phpt_matched = new TH1F("h_phpt_matched", "Photon pT for events: electron lies within #Delta R < 0.3 of HLT objects; pT [GeV]; Events/GeV", 1000, 0, 1000); h_phpt_matched->Sumw2();
  TH1F *h_photon_eta =new TH1F("h_photon_eta", "Photon #eta; #eta ; Events", 600, -3.0, 3.0); h_photon_eta->Sumw2();
  TH1F *h_photon_phi =new TH1F("h_photon_phi", "Photon #phi; #phi ; Events", 800, -4.0, 4.0); h_photon_phi->Sumw2(); 
  TH1F *h_photon_energy =new TH1F("h_photon_energy", "Photon Energy; Energy [GeV]; Events", 1000, 0, 1000); h_photon_energy->Sumw2();

  TH1F *h_DeltaR_el1_ph1 = new TH1F("h_DeltaR_el1_ph1", "#Delta R between leading electron and leading photon; #Delta R; Events", 3500, 0, 3.5);h_DeltaR_el1_ph1->Sumw2();
  TH1F *h_DeltaR_el2_ph1 = new TH1F("h_DeltaR_el2_ph1", "#Delta R between trailing electron and leading photon; #Delta R; Events", 3500, 0, 3.5);h_DeltaR_el2_ph1->Sumw2();

  TH1F *h_DeltaR_elphZ1 = new TH1F("h_DeltaR_elphZ1", "#Delta R between leading electron and leading photon in the Z-mass regime; #Delta R; Events", 3500, 0, 3.5);h_DeltaR_elphZ1->Sumw2();
  TH1F *h_DeltaR_elphZ2 = new TH1F("h_DeltaR_elphZ2", "#Delta R between trailing electron and leading photon in the Z-mass regime; #Delta R; Events", 3500, 0, 3.5);h_DeltaR_elphZ2->Sumw2();

  TH1F *h_minDeltaR = new TH1F("h_minDeltaR", "#Delta R between HLT objects and the leading photon; #Delta R; Events", 3500, 0, 3.5);h_minDeltaR->Sumw2();
  TH1F *h_minDeltaR_el1 = new TH1F("h_minDeltaR_el1", "#Delta R between HLT objects and the leading electron; #Delta R; Events", 3500, 0, 3.5);h_minDeltaR_el1->Sumw2();
  TH1F *h_minDeltaR_el2 = new TH1F("h_minDeltaR_el2", "#Delta R between HLT objects and the trailing electron; #Delta R; Events", 3500, 0, 3.5);h_minDeltaR_el2->Sumw2();
  TH1F *h_DeltaR_elph_unmatched = new TH1F("h_DeltaR_elph_unmatched", "#Delta R between leading electron and leading photon for unmatched electrons; #Delta R; Events", 3500, 0, 3.5);h_DeltaR_elph_unmatched->Sumw2(); 

  TH1F *h_InvariantMass_Mu=new TH1F("h_InvariantMass_Mu", "Di-muon invariant mass; m_{#mu#mu} [GeV]; Events/GeV", 9000, 0, 300); h_InvariantMass_Mu->Sumw2();
  TH1F *h_InvariantMass_El=new TH1F("h_InvariantMass_El", "Di-electron invariant mass; m_{ee} [GeV]; Events/GeV", 9000, 0, 300); h_InvariantMass_El->Sumw2();
  TH1F *h_InvariantMass_MuPh=new TH1F("h_InvariantMass_MuPh", "Di-muon and photon invariant mass; m_{#mu#mu#gamma} [GeV]; Events/GeV", 9000, 0, 300); h_InvariantMass_MuPh->Sumw2();
  TH1F *h_InvariantMass_ElPh=new TH1F("h_InvariantMass_ElPh", "Di-electron and photon invariant mass; m_{ee#gamma} [GeV]; Events/GeV", 9000, 0, 300); h_InvariantMass_ElPh->Sumw2();

  TH1F *h_InvariantMass_ElMu=new TH1F("h_InvariantMass_ElMu", "Electron-Muon invariant mass; m_{e#mu} [GeV]; Events/GeV", 9000, 0, 300); h_InvariantMass_ElMu->Sumw2();
  TH1F *h_cMt_El = new TH1F("h_cMt_El", "Cluster Transverse Mass; M_{T}(e,#gamma,MET) [GeV]; Events/GeV", 3000, 0, 150);h_cMt_El->Sumw2();
  TH1F *h_cMt_Mu = new TH1F("h_cMt_Mu", "Cluster Transverse Mass; M_{T}(#mu,#gamma,MET) [GeV]; Events/GeV", 3000, 0, 150);h_cMt_Mu->Sumw2();
  TH1F *h_Mt_El = new TH1F("h_Mt_El", "Transverse Mass; M_{T}(e, MET) [GeV]; Events/GeV", 3000, 0, 150);h_Mt_El->Sumw2();
  TH1F *h_Mt_Mu = new TH1F("h_Mt_Mu", "Transverse Mass; M_{T}(#mu, MET) [GeV]; Events/GeV", 3000, 0, 150);h_Mt_Mu->Sumw2();
  TH2F *h_Mt_cMt_El = new TH2F("h_Mt_cMt_El", "Scatter Plots of M_{T}(e,MET) versus M_{T}(e,#gamma,MET)}; M_{T}(e,MET) [GeV]; M_{T}(e,#gamma,MET) [GeV]", 4000, 0, 2000, 4000, 0, 2000);h_Mt_cMt_El->Sumw2();
  TH2F *h_Mt_cMt_Mu = new TH2F("h_Mt_cMt_Mu", "Scatter Plots of M_{T}(#mu,MET) versus M_{T}(#mu,#gamma,MET)}; M_{T}(#mu,MET) [GeV]; M_{T}(#mu,#gamma,MET) [GeV]", 4000, 0, 2000, 4000, 0, 2000);h_Mt_cMt_Mu->Sumw2();

  TH1F *h_MET=new TH1F("h_MET", "Missing ET; MET [GeV]; Events/GeV", 600, 0, 600); h_MET->Sumw2();    

  TH1F *h_Difference_Mu=new TH1F("h_Difference_Mu", "Difference pt; Difference [GeV]; Energy/GeV", 800, 0, 800); h_Difference_Mu->Sumw2();
  TH1F *h_Difference_El=new TH1F("h_Difference_El", "Difference pt; Difference [GeV]; Energy/GeV", 800, 0, 800); h_Difference_El->Sumw2();

  TH2F *h_Mmumu_MmumuGamma = new TH2F("h_Mmumu_MmumuGamma", "Scatter Plot of  M_{#mu#mu} versus M_{#mu#mu#gamma}; M_{#mu#mu} [GeV]; M_{#mu#mu#gamma} [GeV]", 4000, 0, 2000, 4000, 0, 2000);  h_Mmumu_MmumuGamma->Sumw2();
  TH2F *h_Mee_MeeGamma = new TH2F("h_Mee_MeeGamma", "Scatter Plot of  M_{ee} versus M_{ee#gamma}; M_{ee} [GeV]; M_{ee#gamma} [GeV]", 4000, 0, 2000, 4000, 0, 20000);h_Mee_MeeGamma->Sumw2();
  TH2F *h_Met_PhPt = new TH2F("h_Met_PhPt", "Scatter Plot of photon pT versus missing ET; Photon pT [GeV]; Missing ET [GeV];", 300, 0, 600, 300, 0, 600); h_Met_PhPt->Sumw2();
  TH2F *h_Met_PhEta = new TH2F("h_Met_PhEta", "Scatter Plot of photon #eta versus missing ET; Photon #eta; Missing ET [GeV];", 400, -2.0, 2.0, 300, 0, 600);h_Met_PhEta->Sumw2();
  TH1F *h_HT = new TH1F("h_HT", "HT (scalar sum of jet pT); H_{T} [GeV]; Events/GeV", 5000, 0, 5000.0);h_HT->Sumw2();
  TH1F *h_nJets = new TH1F("h_nJets", "Number of Jets; Number of Jets; Events", 20, -0.5, 19.5);h_nJets->Sumw2();
  TH2F *h_ph_el1 = new TH2F("h_ph_el1", "Scatter Plot of leading electron pT versus photon pT; Electron pT [GeV]; Photon pT [GeV];", 300, 0, 600, 300, 0, 600); h_ph_el1->Sumw2();
  TH2F *h_ph_el2 = new TH2F("h_ph_el2", "Scatter Plot of trailing electron pT versus photon pT; Electron pT [GeV]; Photon pT [GeV];", 300, 0, 600, 300, 0, 600); h_ph_el2->Sumw2();
  TH2F *h_ph_el1_Eta = new TH2F("h_ph_el1_Eta", "Scatter Plot of leading electron #eta versus photon #eta; Electron #eta; Photon #eta;", 800, -4.0, 4.0, 800, -4.0, 4.0); h_ph_el1_Eta->Sumw2();
  TH2F *h_ph_el2_Eta = new TH2F("h_ph_el2_Eta", "Scatter Plot of trailing electron #eta versus photon #eta; Electron #eta; Photon #eta;", 800, -4.0, 4.0, 800, -4.0, 4.0); h_ph_el2_Eta->Sumw2();
  TH2F *h_ph_el1_Phi = new TH2F("h_ph_el1_Phi", "Scatter Plot of leading electron #phi versus photon #phi; Electron #phi; Photon #phi", 800, -4.0, 4.0, 800, -4.0, 4.0); h_ph_el1_Phi->Sumw2();
  TH2F *h_ph_el2_Phi = new TH2F("h_ph_el2_Phi", "Scatter Plot of trailing electron #phi versus photon #phi; Electron #phi; Photon #phi;", 800, -4.0, 4.0, 800, -4.0, 4.0); h_ph_el2_Phi->Sumw2();
  TH1F *h_el1_pt50 = new TH1F("h_el1_pt50", "Leading electron pT within invariant mass range of 30-60 GeV; pT [GeV]; Events/GeV", 1000, 0, 1000); h_el1_pt50->Sumw2(); 

 int nEvents=tree->GetEntries();
 std::cout << "nEvents= " << nEvents << std::endl;
 for (int i=0; i<nEvents; ++i)
   {
   tree->GetEvent(i);
     //if(fired_HLTPho!=1 and fired_HLTPhoId!=1 and fired_HLTPhoIdMet!=1) continue;

   if(fired_HLTPhoIdMet==1){// continue;

     TriggerInfo trigger1;    
     trigger1.Px = trigObj1Px;
     trigger1.Py = trigObj1Py;
     trigger1.Pz = trigObj1Pz;
     trigger1.E  = trigObj1E;

     h_MET->Fill(MET); 

     // Filling the photon's properties into a vector of struct
     std::vector<PhotonInfo> photons;
     for (unsigned int j=0; j<ph_pt->size(); ++j)
       {
       PhotonInfo photon;
       photon.pT=ph_pt->at(j);
       photon.eta=ph_eta->at(j);
       photon.phi=ph_phi->at(j);
       photon.energy=ph_energy->at(j);
       photon.isTight=ph_isTight->at(j);
       photon.chIsolation=ph_chIso->at(j);
       photon.nuIsolation=ph_nuIso->at(j);
       photon.phIsolation=ph_phIso->at(j);
       photon.phIsoTight=ph_phIsoTight->at(j);
       photon.phIsoMedium=ph_phIsoMedium->at(j);
       photon.phIsoLoose=ph_phIsoLoose->at(j); 
       photons.push_back(photon);
      }
     // Now sorting this vector of structs
     std::sort (photons.begin(), photons.end(), sortPhotonsInDescendingpT);
      
     double phPT = 0; //photons are sorted. Don't care as far as HT is concerned.
     for(unsigned int k=0; k<photons.size(); ++k)
       {
       if(photons.at(k).pT > 30.0 and photons.at(k).isTight==1 and photons.at(k).phIsoTight==1){
         phPT += photons.at(k).pT;
       }
     }
  
     TLorentzVector ph1_p4;
     //here working with the leading photon.
     double deltaR1 = -1.0;
     if (photons.size() > 0){
       ph1_p4=fillTLorentzVector(photons.at(0).pT, photons.at(0).eta, photons.at(0).phi, photons.at(0).energy); 

       int foundHLTPhoton = 0;
       TLorentzVector trigger1_p4;
       trigger1_p4.SetPxPyPzE(trigger1.Px, trigger1.Py, trigger1.Pz, trigger1.E);
       if(ph1_p4.Pt() > 0.0) 
         {
         deltaR1 = ph1_p4.DeltaR(trigger1_p4);
         if(deltaR1>0.0) h_minDeltaR->Fill(deltaR1);
         if(deltaR1<0.3){
           foundHLTPhoton++;
         }
       }
       else{foundHLTPhoton=0;}
     }//only execute this loop if a photon exists.

     // Filling the electron's properties into a vector of struct
     std::vector<LeptonInfo> electrons;
     for (unsigned int j=0; j<el_pt->size(); ++j)
       {
        LeptonInfo electron;
        electron.pT=el_pt->at(j);
        electron.eta=el_eta->at(j);
        electron.phi=el_phi->at(j);
        electron.energy=el_energy->at(j);
        electron.charge=el_charge->at(j);
        electron.isTight=el_isTight->at(j);
        electron.isolation=el_iso->at(j);
        electrons.push_back(electron);
       }

     // Now sorting this vector of structs
     std::sort (electrons.begin(), electrons.end(), sortLeptonsInDescendingpT);
     TLorentzVector el1_p4;
     TLorentzVector el2_p4;

     int foundHLTelectron1 = 0;
     double deltaR_el1 = -1.0;
     if (electrons.size() > 0){
       el1_p4=fillTLorentzVector(electrons.at(0).pT, electrons.at(0).eta, electrons.at(0).phi, electrons.at(0).energy);
       TLorentzVector trigger1_p4;
       trigger1_p4.SetPxPyPzE(trigger1.Px, trigger1.Py, trigger1.Pz, trigger1.E);
       if(el1_p4.Pt() > 0.0)
          {
          deltaR_el1 = el1_p4.DeltaR(trigger1_p4);
          if(deltaR_el1>0.0) h_minDeltaR_el1->Fill(deltaR_el1);
          if(deltaR_el1<0.3){
            foundHLTelectron1++;
          }
       }
       else{foundHLTelectron1=0;}
     }//only execute if an electron exists.

     int foundHLTelectron2 = 0;
     double deltaR_el2 = -1.0;
     if (electrons.size() > 1){
       el2_p4=fillTLorentzVector(electrons.at(1).pT, electrons.at(1).eta, electrons.at(1).phi, electrons.at(1).energy);
       TLorentzVector trigger1_p4;
       trigger1_p4.SetPxPyPzE(trigger1.Px, trigger1.Py, trigger1.Pz, trigger1.E);
       if(el2_p4.Pt() > 0.0)
         {
         deltaR_el2 = el2_p4.DeltaR(trigger1_p4);
         if(deltaR_el2>0.0) h_minDeltaR_el2->Fill(deltaR_el2);
         if(deltaR_el2<0.3){
           foundHLTelectron2++;
        }
      }
      else{foundHLTelectron2=0;}
     }//only execute if the second electron exists.

     if(deltaR_el1 > 0.0 and deltaR_el1 < 0.3){
       if(ph1_p4.Pt()>30.0 and photons.at(0).isTight==1 and photons.at(0).phIsoTight==1){ 
       h_phpt_matched->Fill(ph1_p4.Pt());
       }
     }
     if (ph1_p4.Pt()>0.0 and photons.at(0).isTight==1 and photons.at(0).phIsoTight==1){
       h_photon_pt->Fill(ph1_p4.Pt());
       h_photon_eta->Fill(ph1_p4.Eta()); 
       h_photon_phi->Fill(ph1_p4.Phi());
       h_photon_energy->Fill(ph1_p4.E());
       h_ph_chIsolation_leading->Fill(photons.at(0).chIsolation);
       h_ph_nuIsolation_leading->Fill(photons.at(0).nuIsolation);
       h_ph_phIsolation_leading->Fill(photons.at(0).phIsolation);
       h_Met_PhPt->Fill(photons.at(0).pT, MET);
       h_Met_PhEta->Fill(photons.at(0).eta, MET);
     }  
    
    double leadingDeltaR, trailingDeltaR;
    leadingDeltaR = -1.0;
    trailingDeltaR = -1.0;

    if(el1_p4.Pt()>0.0 and ph1_p4.Pt()>0.0 and photons.at(0).isTight==1 and photons.at(0).phIsoTight==1){
      leadingDeltaR = el1_p4.DeltaR(ph1_p4); 
      if(leadingDeltaR > 0.0) h_DeltaR_el1_ph1->Fill(leadingDeltaR);
    }   

    if(el2_p4.Pt()>0.0 and ph1_p4.Pt()>0.0 and photons.at(0).isTight==1 and photons.at(0).phIsoTight==1){
      trailingDeltaR = el2_p4.DeltaR(ph1_p4); 
      if(trailingDeltaR>0.0) h_DeltaR_el2_ph1->Fill(trailingDeltaR);
      
    }
    
    double cMt2_El = -99.0;
    double Mlg_El = -99.0;
    double pTlg_El = -99.0;
    double mt2_El = -99.0;

    if(el1_p4.Pt()>0.0 and electrons.at(0).isTight==1 and electrons.at(0).isolation < 0.10 and deltaR_el1 > 0.3 and leadingDeltaR > 0.3){
      h_el_phi_leading->Fill(el1_p4.Phi());
      h_el_eta_leading->Fill(el1_p4.Eta());
      h_el_pt_leading->Fill(el1_p4.Pt());
      h_el_energy_leading->Fill(el1_p4.E());
      h_el_isolation_leading->Fill(electrons.at(0).isolation); //there will be a sharp cut at 0.10
      double dphi = fabs(MET_Phi - el1_p4.Phi());
      if (dphi > TMath::Pi()) dphi = TMath::TwoPi() - dphi;
      mt2_El = 2*el1_p4.Pt()*MET*(1 - cos(dphi));
      h_Mt_El->Fill(TMath::Sqrt(mt2_El));
      if(ph1_p4.Pt()>30.0 and photons.at(0).isTight==1 and photons.at(0).phIsoTight==1){
        h_ph_el1->Fill(el1_p4.Pt(), ph1_p4.Pt());
        h_ph_el1_Eta->Fill(el1_p4.Eta(), ph1_p4.Eta());
        h_ph_el1_Phi->Fill(el1_p4.Phi(), ph1_p4.Phi());
        h_DeltaR_elph_unmatched->Fill(leadingDeltaR); 
        //Computing the clustered mass distribution:
        TVector2 el_transverse;
        TVector2 met_transverse;
        TVector2 ph_transverse;
        el_transverse.SetMagPhi(el1_p4.Pt(), el1_p4.Phi());
        ph_transverse.SetMagPhi(ph1_p4.Pt(), ph1_p4.Phi());
        met_transverse.SetMagPhi(MET, MET_Phi);
        Mlg_El = (el1_p4+ph1_p4).M(); 
        pTlg_El = (ph_transverse + el_transverse).Mod2();
        double t1 = TMath::Sqrt(Mlg_El*Mlg_El + pTlg_El);
        cMt2_El=((t1 + MET)*(t1 + MET) - (ph_transverse+el_transverse+met_transverse).Mod2());
        h_cMt_El->Fill(TMath::Sqrt(cMt2_El));
        }
      h_Mt_cMt_El->Fill(TMath::Sqrt(mt2_El), TMath::Sqrt(cMt2_El));  
      if(el2_p4.Pt()>0.0 and electrons.at(1).isTight==1 and electrons.at(1).isolation < 0.10 and trailingDeltaR > 0.3 and ((electrons.at(0).charge*electrons.at(1).charge)==-1) and deltaR_el2 > 0.3) {
        h_el_phi_trailing->Fill(el2_p4.Phi());
        h_el_eta_trailing->Fill(el2_p4.Eta());
        h_el_pt_trailing->Fill(el2_p4.Pt());
        h_el_energy_trailing->Fill(el2_p4.E());
        h_el_isolation_trailing->Fill(electrons.at(1).isolation); //there will be a sharp cut at 0.10
        h_InvariantMass_El->Fill((el1_p4+el2_p4).M());
        h_Difference_El->Fill(el1_p4.Pt() - el2_p4.Pt());
        if(((el1_p4+el2_p4).M() > 30 and (el1_p4+el2_p4).M() < 60))h_el1_pt50->Fill(el1_p4.Pt());
        if(ph1_p4.Pt()>30.0 and photons.at(0).isTight==1 and photons.at(0).phIsoTight==1) {
          if(((el1_p4+el2_p4).M() > 60 and (el1_p4+el2_p4).M() < 120)) h_DeltaR_elphZ1->Fill(ph1_p4.DeltaR(el1_p4));
          if(((el1_p4+el2_p4).M() > 60 and (el1_p4+el2_p4).M() < 120)) h_DeltaR_elphZ2->Fill(ph1_p4.DeltaR(el2_p4)); 
          h_InvariantMass_ElPh->Fill((el1_p4+el2_p4+ph1_p4).M());
          h_Mee_MeeGamma->Fill((el1_p4+el2_p4).M(), (el1_p4+el2_p4+ph1_p4).M());
          h_ph_el2->Fill(el2_p4.Pt(), ph1_p4.Pt());
          h_ph_el2_Eta->Fill(el2_p4.Eta(), ph1_p4.Eta());
          h_ph_el2_Phi->Fill(el2_p4.Phi(), ph1_p4.Phi());
          }
        }
     }

    // Filling the muon's properties into a vector of struct

    std::vector<LeptonInfo> muons;
    for (unsigned int j=0; j<mu_pt->size(); ++j)
      {
       LeptonInfo muon;
       muon.pT=mu_pt->at(j);
       muon.eta=mu_eta->at(j);
       muon.phi=mu_phi->at(j);
       muon.energy=mu_energy->at(j);
       muon.charge=mu_charge->at(j);
       muon.isTight=mu_isTight->at(j);
       muon.isolation=mu_iso->at(j);
       muons.push_back(muon);
       } 
     // Now sorting this vector of structs
    std::sort (muons.begin(), muons.end(), sortLeptonsInDescendingpT);
    TLorentzVector mu1_p4;
    TLorentzVector mu2_p4;
    if (muons.size() > 0) mu1_p4=fillTLorentzVector(muons.at(0).pT, muons.at(0).eta, muons.at(0).phi, muons.at(0).energy);
    if (muons.size() > 1) mu2_p4=fillTLorentzVector(muons.at(1).pT, muons.at(1).eta, muons.at(1).phi, muons.at(1).energy);

    double cMt2_Mu = -99.0;
    double Mlg_Mu = -99.0;
    double pTlg_Mu = -99.0;
    double mt2_Mu = -99.0;
    if(mu1_p4.Pt()>0.0 and muons.at(0).isTight==1 and muons.at(0).isolation < 0.12) {
      h_mu_phi_leading->Fill(mu1_p4.Phi());
      h_mu_eta_leading->Fill(mu1_p4.Eta());
      h_mu_pt_leading->Fill(mu1_p4.Pt());
      h_mu_energy_leading->Fill(mu1_p4.E());
      h_mu_isolation_leading->Fill(muons.at(0).isolation);
      double dphi = fabs(MET_Phi - mu1_p4.Phi());
      if (dphi > TMath::Pi()) dphi = TMath::TwoPi() - dphi;
      mt2_Mu = 2*mu1_p4.Pt()*MET*(1 - TMath::Cos(dphi));
      h_Mt_Mu->Fill(TMath::Sqrt(mt2_Mu));
      if(ph1_p4.Pt()>30.0 and photons.at(0).isTight==1 and photons.at(0).phIsoTight==1){
        //Computing the clustered mass distribution:
        TVector2 mu_transverse;
        TVector2 met_transverse;
        TVector2 ph_transverse;
        mu_transverse.SetMagPhi(mu1_p4.Pt(), mu1_p4.Phi());
        ph_transverse.SetMagPhi(ph1_p4.Pt(), ph1_p4.Phi());
        met_transverse.SetMagPhi(MET, MET_Phi);
        Mlg_Mu = (mu1_p4+ph1_p4).M();
        pTlg_Mu = (ph_transverse + mu_transverse).Mod2();
        double t1 = TMath::Sqrt(Mlg_Mu*Mlg_Mu + pTlg_Mu);
        cMt2_Mu=((t1 + MET)*(t1 + MET) - (ph_transverse+mu_transverse+met_transverse).Mod2());
        h_cMt_Mu->Fill(TMath::Sqrt(cMt2_Mu));
      }
      h_Mt_cMt_Mu->Fill(TMath::Sqrt(mt2_Mu), TMath::Sqrt(cMt2_Mu));
      if(mu2_p4.Pt()>0.0 and muons.at(1).isTight==1 and muons.at(1).isolation < 0.12  and ((muons.at(0).charge*muons.at(1).charge)==-1)) {
        h_mu_phi_trailing->Fill(mu2_p4.Phi());
        h_mu_eta_trailing->Fill(mu2_p4.Eta());
        h_mu_pt_trailing->Fill(mu2_p4.Pt());
        h_mu_energy_trailing->Fill(mu2_p4.E());
        h_mu_isolation_trailing->Fill(muons.at(1).isolation); 
        h_InvariantMass_Mu->Fill((mu1_p4+mu2_p4).M());
        h_Difference_Mu->Fill(mu1_p4.Pt() - mu2_p4.Pt());
        if(ph1_p4.Pt()>0.0 and photons.at(0).isTight==1 and photons.at(0).phIsoTight==1) h_InvariantMass_MuPh->Fill((mu1_p4+mu2_p4+ph1_p4).M());
        if(ph1_p4.Pt()>0.0 and photons.at(0).isTight==1 and photons.at(0).phIsoTight==1) h_Mmumu_MmumuGamma->Fill((mu1_p4+mu2_p4).M(), (mu1_p4+mu2_p4+ph1_p4).M());
        }
     }
   if((mu1_p4.Pt()>0.0 and muons.at(0).isTight==1 and muons.at(0).isolation < 0.12) and (el1_p4.Pt()>0.0 and electrons.at(0).isTight==1 and electrons.at(0).isolation < 0.10 and leadingDeltaR > 0.3)) {
     h_InvariantMass_ElMu->Fill((mu1_p4+el1_p4).M());
     }//El-Mu invariant mass

   std::vector<JetInfo> jets;
   for (unsigned int j=0; j<jet_pt->size(); ++j)
     {
     JetInfo jet;
     jet.pT = jet_pt->at(j);
     jet.eta = jet_eta->at(j);
     jet.phi = jet_phi->at(j);
     jet.energy = jet_energy->at(j);
     jet.PU_mva_loose = jet_mva_loose->at(j);
     jet.PU_mva_tight = jet_mva_tight->at(j);
     jet.PU_mva_medium = jet_mva_medium->at(j);
     jet.PU_cut_loose = jet_cut_loose->at(j);
     jet.PU_cut_tight = jet_cut_tight->at(j);
     jet.PU_cut_medium = jet_cut_medium->at(j);
     jets.push_back(jet);
     }

   // Now sorting this vector of structs
   std::sort (jets.begin(), jets.end(), sortJetsInDescendingpT);

   double HT = 0.0; //jets are sorted. Don't care as far as HT is concerned.
   vector<TLorentzVector> Jet_vector;
   Jet_vector.clear();
   for(unsigned int k=0; k<jets.size(); ++k)
     {
     TLorentzVector Jet;
     if(jets.at(k).PU_mva_loose==1 and jets.at(k).pT>10.0 and fabs(jets.at(k).eta)<2.4){
       Jet.SetPtEtaPhiE(jets.at(k).pT, jets.at(k).eta, jets.at(k).phi, jets.at(k).energy);
     }
     bool isGoodJet=true;
     for(unsigned int j=0; j<electrons.size(); ++j)
       {
       TLorentzVector Electron;
       if(electrons.at(j).isTight==1 and electrons.at(j).isolation < 0.10) Electron.SetPtEtaPhiE(electrons.at(j).pT, electrons.at(j).eta, electrons.at(j).phi, electrons.at(j).energy);
       if(Electron.Pt()>0.0 and Jet.Pt() > 10.0){
         double DRjet_el = Jet.DeltaR(Electron);     
         if(DRjet_el<0.5) isGoodJet=false;
       }
     }
     for(unsigned int j=0; j<muons.size(); ++j)
       {
       TLorentzVector Muon;
       if(muons.at(j).isTight==1 and muons.at(j).isolation < 0.12) Muon.SetPtEtaPhiE(muons.at(j).pT, muons.at(j).eta, muons.at(j).phi, muons.at(j).energy);
       if(Muon.Pt()>0.0 and Jet.Pt() > 10.0){
         double DRjet_mu = Jet.DeltaR(Muon);
         if(DRjet_mu<0.5) isGoodJet=false;
       }
     }  
     for(unsigned int l=0; l<photons.size(); ++l)
       {
       TLorentzVector Photon;
       if(photons.at(l).pT>30.0 and photons.at(l).isTight==1 and photons.at(l).phIsoTight==1) Photon.SetPtEtaPhiE(photons.at(l).pT, photons.at(l).eta, photons.at(l).phi, photons.at(l).energy);
       if(Photon.Pt()>30.0 and Jet.Pt() > 10.0){
         double DRjet_ph = Jet.DeltaR(Photon);
         if(DRjet_ph<0.5) isGoodJet=false; 
       }
     }
    TLorentzVector trigger1_p4;
    trigger1_p4.SetPxPyPzE(trigger1.Px, trigger1.Py, trigger1.Pz, trigger1.E);
    if(trigger1_p4.Px() > 0.0 and Jet.Pt() > 10.0){
      double DRjet_tr = Jet.DeltaR(trigger1_p4);
      if(DRjet_tr<0.5) isGoodJet=false;
    }
    if(isGoodJet and Jet.Pt() > 10.0) Jet_vector.push_back(Jet);
   }

  for(unsigned int m=0; m<Jet_vector.size(); m++)
    {
    HT += Jet_vector.at(m).Pt();
    }

  /*
  if(HT>20 and HT<30){
    
    cout << "HT = " << HT << endl;
    if(Jet_vector.size() > 0 and Jet_vector.at(0).Pt() > 10.0) cout << "Jet_vector.at(0).Pt() = " << Jet_vector.at(0).Pt() << endl;
    if(Jet_vector.size() > 1 and Jet_vector.at(1).Pt() > 10.0) cout << "Jet_vector.at(1).Pt() = " << Jet_vector.at(1).Pt() << endl;
    if(Jet_vector.size() > 2 and Jet_vector.at(2).Pt() > 10.0) cout << "Jet_vector.at(2).Pt() = " << Jet_vector.at(2).Pt() << endl;
    if(Jet_vector.size() > 3 and Jet_vector.at(3).Pt() > 10.0) cout << "Jet_vector.at(3).Pt() = " << Jet_vector.at(3).Pt() << endl;
    if(Jet_vector.size() > 4 and Jet_vector.at(4).Pt() > 10.0) cout << "Jet_vector.at(4).Pt() = " << Jet_vector.at(4).Pt() << endl;
   }
*/
  if(Jet_vector.size()>0 and Jet_vector.at(0).Pt() > 10.0 ) h_jet_pt_leading->Fill(Jet_vector.at(0).Pt());
  if(Jet_vector.size()>1 and Jet_vector.at(1).Pt() > 10.0 ) h_jet_pt_trailing->Fill(Jet_vector.at(1).Pt());
  if(Jet_vector.size()>2 and Jet_vector.at(2).Pt() > 10.0 ) h_jet_pt_3rd->Fill(Jet_vector.at(2).Pt());
  if(Jet_vector.size()>3 and Jet_vector.at(3).Pt() > 10.0 ) h_jet_pt_4th->Fill(Jet_vector.at(3).Pt());
  if(Jet_vector.size()>4 and Jet_vector.at(4).Pt() > 10.0) h_jet_pt_5th->Fill(Jet_vector.at(4).Pt());
  if(Jet_vector.size()>5 and Jet_vector.at(5).Pt() > 10.0) h_jet_pt_6th->Fill(Jet_vector.at(5).Pt());

  if(Jet_vector.size()>0 and Jet_vector.at(0).Pt() > 10.0) h_jet_eta_leading->Fill(Jet_vector.at(0).Eta());
  if(Jet_vector.size()>1 and Jet_vector.at(1).Pt() > 10.0) h_jet_eta_trailing->Fill(Jet_vector.at(1).Eta());
  if(Jet_vector.size()>2 and Jet_vector.at(2).Pt() > 10.0) h_jet_eta_3rd->Fill(Jet_vector.at(2).Eta());
  if(Jet_vector.size()>3 and Jet_vector.at(3).Pt() > 10.0) h_jet_eta_4th->Fill(Jet_vector.at(3).Eta());
  if(Jet_vector.size()>4 and Jet_vector.at(4).Pt() > 10.0) h_jet_eta_5th->Fill(Jet_vector.at(4).Eta());
  if(Jet_vector.size()>5 and Jet_vector.at(5).Pt() > 10.0) h_jet_eta_6th->Fill(Jet_vector.at(5).Eta());

  if(Jet_vector.size()>0 and Jet_vector.at(0).Pt() > 10.0) h_jet_phi_leading->Fill(Jet_vector.at(0).Phi());
  if(Jet_vector.size()>1 and Jet_vector.at(1).Pt() > 10.0) h_jet_phi_trailing->Fill(Jet_vector.at(1).Phi());
  if(Jet_vector.size()>2 and Jet_vector.at(2).Pt() > 10.0) h_jet_phi_3rd->Fill(Jet_vector.at(2).Phi());
  if(Jet_vector.size()>3 and Jet_vector.at(3).Pt() > 10.0) h_jet_phi_4th->Fill(Jet_vector.at(3).Phi());
  if(Jet_vector.size()>4 and Jet_vector.at(4).Pt() > 10.0) h_jet_phi_5th->Fill(Jet_vector.at(4).Phi());
  if(Jet_vector.size()>5 and Jet_vector.at(5).Pt() > 10.0) h_jet_phi_6th->Fill(Jet_vector.at(5).Phi());
  
  if(Jet_vector.size()>0 and Jet_vector.at(0).Pt() > 10.0) h_jet_energy_leading->Fill(Jet_vector.at(0).E());
  if(Jet_vector.size()>1 and Jet_vector.at(1).Pt() > 10.0) h_jet_energy_trailing->Fill(Jet_vector.at(1).E());
  if(Jet_vector.size()>2 and Jet_vector.at(2).Pt() > 10.0) h_jet_energy_3rd->Fill(Jet_vector.at(2).E());
  if(Jet_vector.size()>3 and Jet_vector.at(3).Pt() > 10.0) h_jet_energy_4th->Fill(Jet_vector.at(3).E());
  if(Jet_vector.size()>4 and Jet_vector.at(4).Pt() > 10.0) h_jet_energy_5th->Fill(Jet_vector.at(4).E());
  if(Jet_vector.size()>5 and Jet_vector.at(5).Pt() > 10.0) h_jet_energy_6th->Fill(Jet_vector.at(5).E());

  h_nJets->Fill(Jet_vector.size());
  h_HT->Fill(HT);   


    }//trigger requirement
  }//event loop closed
  std::string histfilename=(outfile+".root").c_str();
  TFile *tFile=new TFile(histfilename.c_str(), "RECREATE");
  h_MET->Write();
  h_nJets->Write();
  h_jet_pt_leading->Write();
  h_jet_pt_trailing->Write();
  h_jet_pt_3rd->Write();
  h_jet_pt_4th->Write();
  h_jet_pt_5th->Write();
  h_jet_pt_6th->Write();
  
  h_jet_eta_leading->Write();
  h_jet_eta_trailing->Write();
  h_jet_eta_3rd->Write();
  h_jet_eta_4th->Write();
  h_jet_eta_5th->Write();
  h_jet_eta_6th->Write();

  h_jet_phi_leading->Write();
  h_jet_phi_trailing->Write();
  h_jet_phi_3rd->Write();
  h_jet_phi_4th->Write();
  h_jet_phi_5th->Write();
  h_jet_phi_6th->Write();

  h_jet_energy_leading->Write();
  h_jet_energy_trailing->Write();
  h_jet_energy_3rd->Write();
  h_jet_energy_4th->Write();
  h_jet_energy_5th->Write();
  h_jet_energy_6th->Write();

  h_cMt_El->Write();
  h_cMt_Mu->Write();
  h_Mt_El->Write();
  h_Mt_Mu->Write();
  h_Mt_cMt_El->Write();
  h_Mt_cMt_Mu->Write();
  h_minDeltaR->Write();
  h_minDeltaR_el1->Write();
  h_minDeltaR_el2->Write();
  h_mu_pt_leading->Write();
  h_mu_pt_trailing->Write();
  h_mu_eta_leading->Write();
  h_mu_eta_trailing->Write();
  h_mu_phi_leading->Write();
  h_mu_phi_trailing->Write();
  h_mu_energy_leading->Write();
  h_mu_energy_trailing->Write();
  h_mu_isolation_leading->Write();
  h_mu_isolation_trailing->Write();


  h_el_pt_leading->Write();
  h_el_pt_trailing->Write();
  h_el_eta_leading->Write();
  h_el_eta_trailing->Write();
  h_el_phi_leading->Write();
  h_el_phi_trailing->Write();
  h_el_energy_leading->Write();
  h_el_energy_trailing->Write();
  h_el_isolation_leading->Write();
  h_el_isolation_trailing->Write();

  h_photon_pt->Write();
  h_photon_eta->Write();
  h_photon_phi->Write();
  h_photon_energy->Write();
  h_phpt_matched->Write();

  h_ph_chIsolation_leading->Write();
  h_ph_nuIsolation_leading->Write();
  h_ph_phIsolation_leading->Write();

  h_DeltaR_el1_ph1->Write();
  h_DeltaR_el2_ph1->Write(); 

  h_DeltaR_elphZ1->Write();
  h_DeltaR_elphZ2->Write();

  h_DeltaR_elph_unmatched->Write();
  h_InvariantMass_Mu->Write();
  h_InvariantMass_El->Write();
  h_InvariantMass_ElPh->Write();
  h_InvariantMass_MuPh->Write();
  h_InvariantMass_ElMu->Write();
  h_Mmumu_MmumuGamma->Write();
  h_Mee_MeeGamma->Write();

  h_Difference_El->Write();
  h_Difference_Mu->Write();

  h_HT->Write();
  h_Met_PhPt->Write();
  h_Met_PhEta->Write();
  h_ph_el1->Write();
  h_ph_el2->Write();
  h_ph_el1_Eta->Write();
  h_ph_el1_Phi->Write();
  h_ph_el2_Eta->Write();
  h_ph_el2_Phi->Write();
  h_el1_pt50->Write();
  tFile->Close(); 
  std::cout<<"Wrote output file "<<histfilename<<std::endl;

  return 0;

}


int TriggerEfficiency(std::string infile, std::string outfile){

  std::string inputfilename=(infile+".root").c_str();
  TChain *tree=new TChain("LowPtSUSY_Tree");
  tree->Add(inputfilename.c_str());
  std::cout<<"Opened input file "<<inputfilename<<std::endl;

  Int_t           nVertices;
  vector<float>   *mu_pt;
  vector<float>   *mu_phi;
  vector<float>   *mu_eta;
  vector<float>   *mu_energy;
  vector<int>     *mu_charge;
  vector<bool>    *mu_isTight;
  vector<float>   *mu_iso;
  vector<float>   *jet_pt;
  vector<float>   *jet_phi;
  vector<float>   *jet_eta;
  vector<float>   *jet_energy;
  Bool_t          fired_HLTPho;
  Bool_t          fired_HLTPhoId;
  Bool_t          fired_HLTPhoIdMet;
  vector<float>   *ph_pt;
  vector<float>   *ph_phi;
  vector<float>   *ph_eta;
  vector<float>   *ph_energy;
  vector<bool>    *ph_phIsoTight;
  vector<bool>    *ph_phIsoMedium;
  vector<bool>    *ph_phIsoLoose;
  vector<float>   *ph_chIso;
  vector<float>   *ph_nuIso;
  vector<float>   *ph_phIso;
  vector<bool>    *ph_isTight;
  Float_t         MET;
  float trigObj1Px;
  float trigObj1Py;
  float trigObj1Pz;
  float trigObj1E;
  float trigObj2Px;
  float trigObj2Py;
  float trigObj2Pz;
  float trigObj2E;

  tree->SetBranchAddress("nVertices", &(nVertices));
  tree->SetBranchAddress("fired_HLTPho", &(fired_HLTPho));
  tree->SetBranchAddress("fired_HLTPhoId", &(fired_HLTPhoId));
  tree->SetBranchAddress("fired_HLTPhoIdMet", &(fired_HLTPhoIdMet));
  tree->SetBranchAddress("ph_pt", &(ph_pt));
  tree->SetBranchAddress("ph_phi", &(ph_phi));
  tree->SetBranchAddress("ph_eta", &(ph_eta));
  tree->SetBranchAddress("ph_energy", &(ph_energy));
  tree->SetBranchAddress("ph_chIso", &(ph_chIso));
  tree->SetBranchAddress("ph_nuIso", &(ph_nuIso));
  tree->SetBranchAddress("ph_phIso", &(ph_phIso));
  tree->SetBranchAddress("ph_isTight", &(ph_isTight));
  tree->SetBranchAddress("ph_phIsoTight", &(ph_phIsoTight));
  tree->SetBranchAddress("ph_phIsoMedium", &(ph_phIsoMedium));
  tree->SetBranchAddress("ph_phIsoLoose", &(ph_phIsoLoose));
  tree->SetBranchAddress("MET", &(MET));
  tree->SetBranchAddress("jet_pt", &(jet_pt));
  tree->SetBranchAddress("jet_phi", &(jet_phi));
  tree->SetBranchAddress("jet_eta", &(jet_eta));
  tree->SetBranchAddress("jet_energy", &(jet_energy));
  tree->SetBranchAddress("mu_pt", &(mu_pt));
  tree->SetBranchAddress("mu_eta", &(mu_eta));
  tree->SetBranchAddress("mu_phi", &(mu_phi));
  tree->SetBranchAddress("mu_energy", &(mu_energy));
  tree->SetBranchAddress("mu_charge", &(mu_charge));
  tree->SetBranchAddress("mu_isTight", &(mu_isTight));
  tree->SetBranchAddress("mu_iso", &(mu_iso));
  tree->SetBranchAddress("trigObj1Px", &(trigObj1Px));
  tree->SetBranchAddress("trigObj1Py", &(trigObj1Py));
  tree->SetBranchAddress("trigObj1Pz", &(trigObj1Pz));
  tree->SetBranchAddress("trigObj1E", &(trigObj1E));
  tree->SetBranchAddress("trigObj2Px", &(trigObj2Px));
  tree->SetBranchAddress("trigObj2Py", &(trigObj2Py));
  tree->SetBranchAddress("trigObj2Pz", &(trigObj2Pz));
  tree->SetBranchAddress("trigObj2E", &(trigObj2E));

  mu_pt = 0;
  mu_phi = 0;
  mu_eta = 0;
  mu_energy = 0;
  mu_charge = 0;
  mu_isTight = 0;
  mu_iso = 0;
  ph_pt = 0;
  ph_phi = 0;
  ph_eta = 0;
  ph_energy = 0;
  ph_chIso = 0;
  ph_nuIso = 0;
  ph_phIso = 0;
  ph_isTight = 0;
  ph_phIsoTight = 0;
  ph_phIsoMedium = 0;
  ph_phIsoLoose = 0;
  jet_pt = 0;
  jet_phi = 0;
  jet_eta = 0;
  jet_energy = 0; 
 
  TH1F *h_MET_HLTPhoId=new TH1F("h_MET_HLTPhoId", "Missing ET; MET [GeV]; Events/GeV", 600, 0, 600); h_MET_HLTPhoId->Sumw2();
  TH1F *h_MET_HLTPhoIdMet=new TH1F("h_MET_HLTPhoIdMet", "Missing ET; MET [GeV]; Events/GeV", 600, 0, 600); h_MET_HLTPhoIdMet->Sumw2();
  TH1F *h_MET_HLTPhoId_MuonVeto=new TH1F("h_MET_HLTPhoId_MuonVeto", "Missing ET; MET [GeV]; Events/GeV", 600, 0, 600); h_MET_HLTPhoId_MuonVeto->Sumw2();
  TH1F *h_MET_HLTPhoIdMet_MuonVeto=new TH1F("h_MET_HLTPhoIdMet_MuonVeto", "Missing ET; MET [GeV]; Events/GeV", 600, 0, 600); h_MET_HLTPhoIdMet_MuonVeto->Sumw2();
  TH1F *h_MET_HLTPhoId_LowHT = new TH1F("h_MET_HLTPhoId_LowHT", "Missing ET; MET [GeV]; Events/GeV", 600, 0, 600);h_MET_HLTPhoId_LowHT->Sumw2(); 
  TH1F *h_MET_HLTPhoIdMet_LowHT = new TH1F("h_MET_HLTPhoIdMet_LowHT", "Missing ET; MET [GeV]; Events/GeV", 600, 0, 600);h_MET_HLTPhoIdMet_LowHT->Sumw2(); 
  TH1F *h_MET_HLTPhoId_HighHT = new TH1F("h_MET_HLTPhoId_HighHT", "Missing ET; MET [GeV]; Events/GeV", 600, 0, 600);h_MET_HLTPhoId_HighHT->Sumw2();      
  TH1F *h_MET_HLTPhoIdMet_HighHT = new TH1F("h_MET_HLTPhoIdMet_HighHT", "Missing ET; MET [GeV]; Events/GeV", 600, 0, 600);h_MET_HLTPhoIdMet_HighHT->Sumw2();
  TH1F *h_MET_HLTPhoId_LowPU = new TH1F("h_MET_HLTPhoId_LowPU", "Missing ET; MET [GeV]; Events/GeV", 600, 0, 600);h_MET_HLTPhoId_LowPU->Sumw2();
  TH1F *h_MET_HLTPhoIdMet_LowPU = new TH1F("h_MET_HLTPhoIdMet_LowPU", "Missing ET; MET [GeV]; Events/GeV", 600, 0, 600);h_MET_HLTPhoIdMet_LowPU->Sumw2();
  TH1F *h_MET_HLTPhoId_HighPU = new TH1F("h_MET_HLTPhoId_HighPU", "Missing ET; MET [GeV]; Events/GeV", 600, 0, 600);h_MET_HLTPhoId_HighPU->Sumw2();
  TH1F *h_MET_HLTPhoIdMet_HighPU = new TH1F("h_MET_HLTPhoIdMet_HighPU", "Missing ET; MET [GeV]; Events/GeV", 600, 0, 600);h_MET_HLTPhoIdMet_HighPU->Sumw2();
  TH1F *h_PHO_HLTPho=new TH1F("h_PHO_HLTPho", "Missing ET; PHO [GeV]; Events/GeV", 600, 0, 600); h_PHO_HLTPho->Sumw2();
  TH1F *h_PHO_HLTPhoId=new TH1F("h_PHO_HLTPhoId", "Missing ET; PHO [GeV]; Events/GeV", 600, 0, 600); h_PHO_HLTPhoId->Sumw2();

  int nEvents=tree->GetEntries();
  std::cout << "nEvents= " << nEvents << std::endl;
  for (int i=0; i<nEvents; ++i)
    {
     tree->GetEvent(i);

     TriggerInfo trigger1;
     trigger1.Px = trigObj1Px;
     trigger1.Py = trigObj1Py;
     trigger1.Pz = trigObj1Pz;
     trigger1.E  = trigObj1E;

     TriggerInfo trigger2;
     trigger2.Px = trigObj2Px;
     trigger2.Py = trigObj2Py;
     trigger2.Pz = trigObj2Pz;
     trigger2.E  = trigObj2E;

     std::vector<LeptonInfo> muons;
     for (unsigned int j=0; j<mu_pt->size(); ++j)
       {
       LeptonInfo muon;
       muon.pT=mu_pt->at(j);
       muon.eta=mu_eta->at(j);
       muon.phi=mu_phi->at(j);
       muon.energy=mu_energy->at(j);
       muon.charge=mu_charge->at(j);
       muon.isTight=mu_isTight->at(j);
       muon.isolation=mu_iso->at(j);
       muons.push_back(muon);
       }

     std::vector<JetInfo> jets;
     for (unsigned int j=0; j<jet_pt->size(); ++j)
       {
       JetInfo jet;
       jet.pT = jet_pt->at(j);
       jet.eta = jet_eta->at(j);
       jet.phi = jet_phi->at(j);
       jet.energy = jet_energy->at(j);
       jets.push_back(jet);
       }

     // Now sorting this vector of structs
     std::sort (jets.begin(), jets.end(), sortJetsInDescendingpT);

     double HT = 0; //jets are sorted. Don't care as far as HT is concerned.
     for(unsigned int k=0; k<jets.size(); ++k)
       {
       HT += jets.at(k).pT;
       }

     std::vector<PhotonInfo> photons;
     for (unsigned int j=0; j<ph_pt->size(); ++j)
       {
       PhotonInfo photon;
       photon.pT=ph_pt->at(j);
       photon.eta=ph_eta->at(j);
       photon.phi=ph_phi->at(j);
       photon.energy=ph_energy->at(j);
       photon.isTight=ph_isTight->at(j);
       photon.chIsolation=ph_chIso->at(j);
       photon.nuIsolation=ph_nuIso->at(j);
       photon.phIsolation=ph_phIso->at(j);
       photon.phIsoTight=ph_phIsoTight->at(j);
       photon.phIsoMedium=ph_phIsoMedium->at(j);
       photon.phIsoLoose=ph_phIsoLoose->at(j);
       photons.push_back(photon);
      }
      // Now sorting this vector of structs
      std::sort (photons.begin(), photons.end(), sortPhotonsInDescendingpT);

     TLorentzVector ph1_p4;
     int foundHLTPhoton1 = 0;
     int foundHLTPhoton2 = 0;
     if (photons.size() > 0){
       ph1_p4=fillTLorentzVector(photons.at(0).pT, photons.at(0).eta, photons.at(0).phi, photons.at(0).energy);
       TLorentzVector trigger1_p4;
       trigger1_p4.SetPxPyPzE(trigger1.Px, trigger1.Py, trigger1.Pz, trigger1.E);
       double deltaR1 = ph1_p4.DeltaR(trigger1_p4);
         if(deltaR1 < 0.3 and fired_HLTPhoId){
            foundHLTPhoton1=1;
         }
       TLorentzVector trigger2_p4;
       trigger2_p4.SetPxPyPzE(trigger2.Px, trigger2.Py, trigger2.Pz, trigger2.E);
       double deltaR2 = ph1_p4.DeltaR(trigger2_p4);
         if(deltaR2 < 0.3 and fired_HLTPho){
            foundHLTPhoton2=1;
         } 
    }//only execute this loop if a photon exists.  
 
   if(photons.size() > 0){
     if(fired_HLTPho==1 and foundHLTPhoton2==1 and photons.at(0).pT>30.0 and photons.at(0).isTight==1 and photons.at(0).phIsoTight==1){  
       h_PHO_HLTPho->Fill(photons.at(0).pT);
       if(foundHLTPhoton1==1 and fired_HLTPhoId==1){
         h_PHO_HLTPhoId->Fill(photons.at(0).pT);
       }//main trigger if   
      }//control trigger photon 30 if
   }

     if (fired_HLTPhoId==1){  
       h_MET_HLTPhoId->Fill(MET);
       if(fired_HLTPhoIdMet==1){
         h_MET_HLTPhoIdMet->Fill(MET);
          }
        }
     
     if (fired_HLTPhoId==1 and muons.size() == 0){ //muon veto
       h_MET_HLTPhoId_MuonVeto->Fill(MET);
       if(fired_HLTPhoIdMet==1){
         h_MET_HLTPhoIdMet_MuonVeto->Fill(MET);
         }
       }

     //beginning binning in HT
     if (fired_HLTPhoId==1 and HT<100){
       h_MET_HLTPhoId_LowHT->Fill(MET);
       if(fired_HLTPhoIdMet==1){
         h_MET_HLTPhoIdMet_LowHT->Fill(MET);
         }
       }

    else if (fired_HLTPhoId==1 and HT>100){
       h_MET_HLTPhoId_HighHT->Fill(MET);
       if(fired_HLTPhoIdMet==1){
         h_MET_HLTPhoIdMet_HighHT->Fill(MET);
         }
       }//binning in HT done
      
    //binning in nVertices
    if (fired_HLTPhoId==1 and nVertices<15){
       h_MET_HLTPhoId_LowPU->Fill(MET);
       if(fired_HLTPhoIdMet==1){
         h_MET_HLTPhoIdMet_LowPU->Fill(MET);
         }
       }

    else if (fired_HLTPhoId==1 and nVertices>15){
       h_MET_HLTPhoId_HighPU->Fill(MET);
       if(fired_HLTPhoIdMet==1){
         h_MET_HLTPhoIdMet_HighPU->Fill(MET);
         }
       }   //binning in nVertices done



  }//end of event loop
  std::string histfilename=(outfile+".root").c_str();
  TFile *tFile=new TFile(histfilename.c_str(), "RECREATE");
  
  h_MET_HLTPhoIdMet->Write();
  h_MET_HLTPhoId->Write();
  h_MET_HLTPhoIdMet_MuonVeto->Write();
  h_MET_HLTPhoId_MuonVeto->Write();
  h_MET_HLTPhoIdMet_LowHT->Write();
  h_MET_HLTPhoId_LowHT->Write();
  h_MET_HLTPhoIdMet_HighHT->Write();
  h_MET_HLTPhoId_HighHT->Write();
  h_MET_HLTPhoIdMet_LowPU->Write();
  h_MET_HLTPhoId_LowPU->Write();
  h_MET_HLTPhoIdMet_HighPU->Write();
  h_MET_HLTPhoId_HighPU->Write();
  h_PHO_HLTPhoId->Write();
  h_PHO_HLTPho->Write();
  tFile->Close();
  std::cout<<"Wrote output file "<<histfilename<<std::endl;

  return 0;

}
