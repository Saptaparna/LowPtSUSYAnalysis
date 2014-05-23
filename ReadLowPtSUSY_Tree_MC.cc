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

bool sortMatchedLeptonsInDescendingpT(MatchedLeptonInfo mlep1, MatchedLeptonInfo mlep2)
{
  return (mlep1.pT > mlep2.pT);
}

int ReadLowPtSUSY_Tree_MC(std::string infile, std::string outfile){
  
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
  Int_t           nJets;
  Float_t         MET;
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
  vector<int>     *el_Matched;
  vector<float>   *el_MatchedPt;
  vector<float>   *el_MatchedEta;
  vector<float>   *el_MatchedPhi;
  vector<float>   *el_MatchedEnergy;
  vector<int>     *mu_Matched;
  vector<float>   *mu_MatchedPt;
  vector<float>   *mu_MatchedEta;
  vector<float>   *mu_MatchedPhi;
  vector<float>   *mu_MatchedEnergy;

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
  el_iso = 0;
  mu_iso = 0;
  ph_chIso = 0;
  ph_nuIso = 0;
  ph_phIso = 0;
  ph_isTight = 0;
  ph_phIsoTight = 0;
  ph_phIsoMedium = 0;
  ph_phIsoLoose = 0;
  el_Matched = 0;
  el_MatchedPt = 0;
  el_MatchedEta = 0;
  el_MatchedPhi = 0;
  el_MatchedEnergy = 0;
  mu_Matched = 0;
  mu_MatchedPt = 0;
  mu_MatchedEta = 0;
  mu_MatchedPhi = 0;
  mu_MatchedEnergy = 0;

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
  tree->SetBranchAddress("nJets", &(nJets));
  tree->SetBranchAddress("MET", &(MET));
  tree->SetBranchAddress("nVertices", &(nVertices));
  tree->SetBranchAddress("nPUVertices", &(nPUVertices));
  tree->SetBranchAddress("nPUVerticesTrue", &(nPUVerticesTrue));
  tree->SetBranchAddress("PUWeightData", &(PUWeightData));
  tree->SetBranchAddress("PUWeightDataSys", &(PUWeightDataSys)); 
  tree->SetBranchAddress("ph_phIsoTight", &(ph_phIsoTight));
  tree->SetBranchAddress("ph_phIsoMedium", &(ph_phIsoMedium));
  tree->SetBranchAddress("ph_phIsoLoose", &(ph_phIsoLoose));
  tree->SetBranchAddress("el_Matched", &(el_Matched));
  tree->SetBranchAddress("el_MatchedPt", &(el_MatchedPt));
  tree->SetBranchAddress("el_MatchedEta", &(el_MatchedEta));
  tree->SetBranchAddress("el_MatchedPhi", &(el_MatchedPhi));
  tree->SetBranchAddress("el_MatchedEnergy", &(el_MatchedEnergy));
  tree->SetBranchAddress("mu_Matched", &(mu_Matched));
  tree->SetBranchAddress("mu_MatchedPt", &(mu_MatchedPt));
  tree->SetBranchAddress("mu_MatchedEta", &(mu_MatchedEta));
  tree->SetBranchAddress("mu_MatchedPhi", &(mu_MatchedPhi));
  tree->SetBranchAddress("mu_MatchedEnergy", &(mu_MatchedEnergy));

  //Booking histograms:
  
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
  TH1F *h_photon_eta =new TH1F("h_photon_eta", "Photon #eta; #eta ; Events", 600, -3.0, 3.0); h_photon_eta->Sumw2();
  TH1F *h_photon_phi =new TH1F("h_photon_phi", "Photon #phi; #phi ; Events", 800, -4.0, 4.0); h_photon_phi->Sumw2(); 
  TH1F *h_photon_energy =new TH1F("h_photon_energy", "Photon Energy; Energy [GeV]; Events", 1000, 0, 1000); h_photon_energy->Sumw2();

  TH1F *h_DeltaR_el1_ph1 = new TH1F("h_DeltaR_el1_ph1", "#Delta R between leading electron and leading photon; #Delta R; Events", 3500, 0, 3.5);h_DeltaR_el1_ph1->Sumw2();
  TH1F *h_DeltaR_el2_ph1 = new TH1F("h_DeltaR_el2_ph1", "#Delta R between trailing electron and leading photon; #Delta R; Events", 3500, 0, 3.5);h_DeltaR_el2_ph1->Sumw2();

  TH1F *h_InvariantMass_Mu=new TH1F("h_InvariantMass_Mu", "Di-muon invariant mass; m_{#mu#mu} [GeV]; Events/GeV", 9000, 0, 300); h_InvariantMass_Mu->Sumw2();
  TH1F *h_InvariantMass_El=new TH1F("h_InvariantMass_El", "Di-electron invariant mass; m_{ee} [GeV]; Events/GeV", 9000, 0, 300); h_InvariantMass_El->Sumw2();
  TH1F *h_InvariantMass_MuPh=new TH1F("h_InvariantMass_MuPh", "Di-muon and photon invariant mass; m_{#mu#mu#gamma} [GeV]; Events/GeV", 9000, 0, 300); h_InvariantMass_MuPh->Sumw2();
  TH1F *h_InvariantMass_ElPh=new TH1F("h_InvariantMass_ElPh", "Di-electron and photon invariant mass; m_{ee#gamma} [GeV]; Events/GeV", 9000, 0, 300); h_InvariantMass_ElPh->Sumw2();

  TH1F *h_InvariantMass_ElMu=new TH1F("h_InvariantMass_ElMu", "Electron-Muon invariant mass; m_{e#mu} [GeV]; Events/GeV", 9000, 0, 300); h_InvariantMass_ElMu->Sumw2();

  TH1F *h_InvariantMass_MCMu=new TH1F("h_InvariantMass_MCMu", "Di-muon invariant mass of muons from Tau mothers; m_{#mu#mu} [GeV]; Events/GeV", 9000, 0, 300); h_InvariantMass_MCMu->Sumw2();

  TH1F *h_MET=new TH1F("h_MET", "Missing ET; MET [GeV]; Events/GeV", 600, 0, 600); h_MET->Sumw2();    

  TH1F *h_Difference_Mu=new TH1F("h_Difference_Mu", "Difference pt; Difference [GeV]; Energy/GeV", 800, 0, 800); h_Difference_Mu->Sumw2();
  TH1F *h_Difference_El=new TH1F("h_Difference_El", "Difference pt; Difference [GeV]; Energy/GeV", 800, 0, 800); h_Difference_El->Sumw2();

  TH2F *h_Mmumu_MmumuGamma = new TH2F("h_Mmumu_MmumuGamma", "Scatter Plot of  M_{#mu#mu} versus M_{#mu#mu#gamma}; M_{#mu#mu} [GeV]; M_{#mu#mu#gamma} [GeV]", 4000, 0, 2000, 4000, 0, 2000);  h_Mmumu_MmumuGamma->Sumw2();
  TH2F *h_Mee_MeeGamma = new TH2F("h_Mee_MeeGamma", "Scatter Plot of  M_{ee} versus M_{ee#gamma}; M_{ee} [GeV]; M_{ee#gamma} [GeV]", 4000, 0, 2000, 4000, 0, 20000);h_Mee_MeeGamma->Sumw2();
  TH1F *h_HT = new TH1F("h_HT", "HT (scalar sum of jet pT); H_T [GeV]; Events/GeV", 5000, 0, 5000.0);h_HT->Sumw2();

  int nEvents=tree->GetEntries();
  std::cout << "nEvents= " << nEvents << std::endl;
  for (int i=0; i<nEvents ; ++i)
    {
     tree->GetEvent(i);

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
       
     TLorentzVector ph1_p4;
     //here working with the leading photon.
     if (photons.size() > 0) ph1_p4=fillTLorentzVector(photons.at(0).pT, photons.at(0).eta, photons.at(0).phi, photons.at(0).energy); 
     
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

     if (electrons.size() > 0) el1_p4=fillTLorentzVector(electrons.at(0).pT, electrons.at(0).eta, electrons.at(0).phi, electrons.at(0).energy);
     if (electrons.size() > 1) el2_p4=fillTLorentzVector(electrons.at(1).pT, electrons.at(1).eta, electrons.at(1).phi, electrons.at(1).energy);

     if (ph1_p4.Pt()>0.0 and photons.at(0).isTight==1 and photons.at(0).phIsoTight==1){
       
       h_photon_pt->Fill(ph1_p4.Pt());
       h_photon_eta->Fill(ph1_p4.Eta()); 
       h_photon_phi->Fill(ph1_p4.Phi());
       h_photon_energy->Fill(ph1_p4.E());
       h_ph_chIsolation_leading->Fill(photons.at(0).chIsolation);
       h_ph_nuIsolation_leading->Fill(photons.at(0).nuIsolation);
       h_ph_phIsolation_leading->Fill(photons.at(0).phIsolation);
     }
    
    double leadingDeltaR, trailingDeltaR;
    leadingDeltaR = 99.0;
    trailingDeltaR = 99.0;

    if(el1_p4.Pt()>0.0 and ph1_p4.Pt()>0.0){

      leadingDeltaR = el1_p4.DeltaR(ph1_p4); 
      h_DeltaR_el1_ph1->Fill(leadingDeltaR);
    }   

    if(el2_p4.Pt()>0.0 and ph1_p4.Pt()>0.0){

      trailingDeltaR = el2_p4.DeltaR(ph1_p4); 
      h_DeltaR_el2_ph1->Fill(trailingDeltaR);
      
    }

    if(el1_p4.Pt()>0.0 and electrons.at(0).isTight==1 and electrons.at(0).isolation < 0.10 and leadingDeltaR > 0.3) {
      h_el_phi_leading->Fill(el1_p4.Phi());
      h_el_eta_leading->Fill(el1_p4.Eta());
      h_el_pt_leading->Fill(el1_p4.Pt());
      h_el_energy_leading->Fill(el1_p4.E());
      h_el_isolation_leading->Fill(electrons.at(0).isolation); //there will be a sharp cut at 0.10
      if(el2_p4.Pt()>0.0 and electrons.at(1).isTight==1 and electrons.at(1).isolation < 0.10 and trailingDeltaR > 0.3 and ((electrons.at(0).charge*electrons.at(1).charge)==-1)) {
        h_el_phi_trailing->Fill(el2_p4.Phi());
        h_el_eta_trailing->Fill(el2_p4.Eta());
        h_el_pt_trailing->Fill(el2_p4.Pt());
        h_el_energy_trailing->Fill(el2_p4.E());
        h_el_isolation_trailing->Fill(electrons.at(1).isolation); //there will be a sharp cut at 0.10
        h_InvariantMass_El->Fill((el1_p4+el2_p4).M());
        h_Difference_El->Fill(el1_p4.Pt() - el2_p4.Pt());
        if(ph1_p4.Pt()>0.0 and photons.at(0).isTight==1 and photons.at(0).phIsoTight==1) h_InvariantMass_ElPh->Fill((el1_p4+el2_p4+ph1_p4).M()); 
        if(ph1_p4.Pt()>0.0 and photons.at(0).isTight==1 and photons.at(0).phIsoTight==1) h_Mee_MeeGamma->Fill((el1_p4+el2_p4).M(), (el1_p4+el2_p4+ph1_p4).M());
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

    if(mu1_p4.Pt()>0.0 and muons.at(0).isTight==1 and muons.at(0).isolation < 0.12) {
      h_mu_phi_leading->Fill(mu1_p4.Phi());
      h_mu_eta_leading->Fill(mu1_p4.Eta());
      h_mu_pt_leading->Fill(mu1_p4.Pt());
      h_mu_energy_leading->Fill(mu1_p4.E());
      h_mu_isolation_leading->Fill(muons.at(0).isolation);
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

   //Matched muons are muons matched exclusively to tau leptons
   std::vector<MatchedLeptonInfo> matchedMuons;
   for (unsigned int j=0; j<mu_Matched->size(); ++j)
     {
     MatchedLeptonInfo matchedMuon;
     matchedMuon.pT = mu_MatchedPt->at(j);
     matchedMuon.eta = mu_MatchedEta->at(j);
     matchedMuon.phi = mu_MatchedPhi->at(j);
     matchedMuon.energy = mu_MatchedEnergy->at(j);
     matchedMuons.push_back(matchedMuon);
     }
   std::sort (matchedMuons.begin(), matchedMuons.end(), sortMatchedLeptonsInDescendingpT);
   
   TLorentzVector mu1_m;
   TLorentzVector mu2_m;
   if (matchedMuons.size() > 0) mu1_m=fillTLorentzVector(matchedMuons.at(0).pT, matchedMuons.at(0).eta, matchedMuons.at(0).phi, matchedMuons.at(0).energy);
   if (matchedMuons.size() > 1) mu2_m=fillTLorentzVector(matchedMuons.at(1).pT, matchedMuons.at(1).eta, matchedMuons.at(1).phi, matchedMuons.at(1).energy);

   if(mu1_m.Pt()>0.0) {
     if(mu2_m.Pt()>0.0){
       h_InvariantMass_MCMu->Fill((mu1_m+mu2_m).M());
        }
     }

  }//event loop closed
  std::string histfilename=(outfile+".root").c_str();
  TFile *tFile=new TFile(histfilename.c_str(), "RECREATE");
  h_MET->Write();
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

  h_ph_chIsolation_leading->Write();
  h_ph_nuIsolation_leading->Write();
  h_ph_phIsolation_leading->Write();

  h_DeltaR_el1_ph1->Write();
  h_DeltaR_el2_ph1->Write(); 

  h_InvariantMass_Mu->Write();
  h_InvariantMass_El->Write();
  h_InvariantMass_ElPh->Write();
  h_InvariantMass_MuPh->Write();
  h_InvariantMass_ElMu->Write();
  h_Mmumu_MmumuGamma->Write();
  h_Mee_MeeGamma->Write();

  h_Difference_El->Write();
  h_Difference_Mu->Write();
 
  h_InvariantMass_MCMu->Write();

  h_HT->Write();
  tFile->Close(); 
  std::cout<<"Wrote output file "<<histfilename<<std::endl;

  return 0;

}
