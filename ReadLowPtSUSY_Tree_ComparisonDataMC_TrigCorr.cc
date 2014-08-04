#include "ReadLowPtSUSY_Tree_ComparisonDataMC_TrigCorr.h"

int ReadLowPtSUSY_Tree_ComparisonDataMC_TrigCorr(std::string infile, std::string outfile, std::string type, std::string channel){
  
  std::string inputfilename=(infile+".root").c_str();
  TChain *tree=new TChain("LowPtSUSY_Tree");
  tree->Add(inputfilename.c_str());
  std::cout<<"Opened input file "<<inputfilename<<std::endl;

  Int_t           run;
  Int_t           lumi;
  Int_t           event;
  Bool_t          fired_HLTPho;
  Bool_t          fired_HLTPhoId;
  Bool_t          fired_HLTPhoIdMet;
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
  bool            Wt;
  bool            Ztt;
  bool            Zee;
  bool            Zmumu;
  bool            Znunu;
  bool            isPhoton;
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
  if(type=="MC")
    {
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
   }
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
  tree->SetBranchAddress("nVertices", &(nVertices));
  tree->SetBranchAddress("nPUVertices", &(nPUVertices));
  tree->SetBranchAddress("nPUVerticesTrue", &(nPUVerticesTrue));
  tree->SetBranchAddress("PUWeightData", &(PUWeightData));
  tree->SetBranchAddress("PUWeightDataSys", &(PUWeightDataSys)); 
  tree->SetBranchAddress("ph_phIsoTight", &(ph_phIsoTight));
  tree->SetBranchAddress("ph_phIsoMedium", &(ph_phIsoMedium));
  tree->SetBranchAddress("ph_phIsoLoose", &(ph_phIsoLoose));
  if(type=="Data"){
    tree->SetBranchAddress("fired_HLTPho", &(fired_HLTPho));
    tree->SetBranchAddress("fired_HLTPhoId", &(fired_HLTPhoId));
    tree->SetBranchAddress("fired_HLTPhoIdMet", &(fired_HLTPhoIdMet));
    tree->SetBranchAddress("trigObj1Px", &(trigObj1Px));
    tree->SetBranchAddress("trigObj1Py", &(trigObj1Py));
    tree->SetBranchAddress("trigObj1Pz", &(trigObj1Pz));
    tree->SetBranchAddress("trigObj1E", &(trigObj1E));
    tree->SetBranchAddress("trigObj2Px", &(trigObj2Px));
    tree->SetBranchAddress("trigObj2Py", &(trigObj2Py));
    tree->SetBranchAddress("trigObj2Pz", &(trigObj2Pz));
    tree->SetBranchAddress("trigObj2E", &(trigObj2E));
  }
  if(type=="MC"){
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
    tree->SetBranchAddress("Wt", &(Wt));
    tree->SetBranchAddress("Ztt", &(Ztt));
    tree->SetBranchAddress("Zee", &(Zee));
    tree->SetBranchAddress("Zmumu", &(Zmumu));
    tree->SetBranchAddress("Znunu", &(Znunu));
    tree->SetBranchAddress("isPhoton", &(isPhoton));
  }
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

  TH1F *h_photon_pt =new TH1F("h_photon_pt", "Photon pT; pT [GeV]; Events/GeV", 1000, 0, 1000); h_photon_pt->Sumw2();
  TH1F *h_photon_pt_Cut2 =new TH1F("h_photon_pt_Cut2", "Photon pT; pT [GeV]; Events/GeV", 1000, 0, 1000); h_photon_pt_Cut2->Sumw2();
  TH1F *h_photon_pt_Cut3 =new TH1F("h_photon_pt_Cut3", "Photon pT; pT [GeV]; Events/GeV", 1000, 0, 1000); h_photon_pt_Cut3->Sumw2();
  TH1F *h_photon_pt_Cut4 =new TH1F("h_photon_pt_Cut4", "Photon pT; pT [GeV]; Events/GeV", 1000, 0, 1000); h_photon_pt_Cut4->Sumw2();
  TH1F *h_photon_eta =new TH1F("h_photon_eta", "Photon #eta; #eta ; Events", 600, -3.0, 3.0); h_photon_eta->Sumw2();
  TH1F *h_photon_phi =new TH1F("h_photon_phi", "Photon #phi; #phi ; Events", 800, -4.0, 4.0); h_photon_phi->Sumw2(); 
  TH1F *h_photon_energy =new TH1F("h_photon_energy", "Photon Energy; Energy [GeV]; Events", 1000, 0, 1000); h_photon_energy->Sumw2();

  TH1F *h_InvariantMass_Mu=new TH1F("h_InvariantMass_Mu", "Di-muon invariant mass; m_{#mu#mu} [GeV]; Events/GeV", 9000, 0, 300); h_InvariantMass_Mu->Sumw2();
  TH1F *h_InvariantMass_Mu_Cut1=new TH1F("h_InvariantMass_Mu_Cut1", "Di-muon invariant mass; m_{#mu#mu} [GeV]; Events/GeV", 9000, 0, 300); h_InvariantMass_Mu_Cut1->Sumw2();
  TH1F *h_InvariantMass_Mu_Cut2=new TH1F("h_InvariantMass_Mu_Cut2", "Di-muon invariant mass; m_{#mu#mu} [GeV]; Events/GeV", 9000, 0, 300); h_InvariantMass_Mu_Cut2->Sumw2();
  TH1F *h_InvariantMass_Mu_Cut3=new TH1F("h_InvariantMass_Mu_Cut3", "Di-muon invariant mass; m_{#mu#mu} [GeV]; Events/GeV", 9000, 0, 300); h_InvariantMass_Mu_Cut3->Sumw2();
  TH1F *h_InvariantMass_Mu_Cut4=new TH1F("h_InvariantMass_Mu_Cut4", "Di-muon invariant mass; m_{#mu#mu} [GeV]; Events/GeV", 9000, 0, 300); h_InvariantMass_Mu_Cut4->Sumw2();
  TH1F *h_InvariantMass_El=new TH1F("h_InvariantMass_El", "Di-electron invariant mass; m_{ee} [GeV]; Events/GeV", 9000, 0, 300); h_InvariantMass_El->Sumw2();
  TH1F *h_InvariantMass_MuPh=new TH1F("h_InvariantMass_MuPh", "Di-muon and photon invariant mass; m_{#mu#mu#gamma} [GeV]; Events/GeV", 9000, 0, 300); h_InvariantMass_MuPh->Sumw2();
  TH1F *h_InvariantMass_MuPh_Cut2=new TH1F("h_InvariantMass_MuPh_Cut2", "Di-muon and photon invariant mass; m_{#mu#mu#gamma} [GeV]; Events/GeV", 9000, 0, 300); h_InvariantMass_MuPh_Cut2->Sumw2();
  TH1F *h_InvariantMass_MuPh_Cut3=new TH1F("h_InvariantMass_MuPh_Cut3", "Di-muon and photon invariant mass; m_{#mu#mu#gamma} [GeV]; Events/GeV", 9000, 0, 300); h_InvariantMass_MuPh_Cut3->Sumw2();
  TH1F *h_InvariantMass_MuPh_Cut4=new TH1F("h_InvariantMass_MuPh_Cut4", "Di-muon and photon invariant mass; m_{#mu#mu#gamma} [GeV]; Events/GeV", 9000, 0, 300); h_InvariantMass_MuPh_Cut4->Sumw2();
  TH1F *h_InvariantMass_ElPh=new TH1F("h_InvariantMass_ElPh", "Di-electron and photon invariant mass; m_{ee#gamma} [GeV]; Events/GeV", 9000, 0, 300); h_InvariantMass_ElPh->Sumw2();

  TH1F *h_InvariantMass_ElMu=new TH1F("h_InvariantMass_ElMu", "Electron-Muon invariant mass; m_{e#mu} [GeV]; Events/GeV", 9000, 0, 300); h_InvariantMass_ElMu->Sumw2();

  TH1F *h_MET=new TH1F("h_MET", "Missing ET; MET [GeV]; Events/GeV", 600, 0, 600); h_MET->Sumw2();    
  TH1F *h_MET_Cut0=new TH1F("h_MET_Cut0", "Missing ET; MET [GeV]; Events/GeV", 600, 0, 600); h_MET_Cut0->Sumw2();  

  TH1F *h_nVertices=new TH1F("h_nVertices", "Number of vertices; nvertices; Events", 30, -0.5, 29.5); h_nVertices->Sumw2();
  TH1F *h_nPUVertices=new TH1F("h_nPUVertices", "Number of PU vertices; nvertices; Events", 30, -0.5, 29.5); h_nPUVertices->Sumw2();
  TH1F *h_nPUVerticesTrue=new TH1F("h_nPUVerticesTrue", "Number of PU vertices; nvertices; Events", 30, -0.5, 29.5); h_nPUVerticesTrue->Sumw2();

  TH1F *h_HT = new TH1F("h_HT", "HT (scalar sum of jet pT); H_T [GeV]; Events/GeV", 5000, 0, 5000.0);h_HT->Sumw2();
  TH1F *h_cMt_Mu = new TH1F("h_cMt_Mu", "Cluster Transverse Mass; M_{T}(#mu,#gamma,MET) [GeV]; Events/GeV", 5000, 0, 500);h_cMt_Mu->Sumw2();
  TH1F *h_Mt_Mu = new TH1F("h_Mt_Mu", "Transverse Mass; M_{T}(#mu, MET) [GeV]; Events/GeV", 3000, 0, 150);h_Mt_Mu->Sumw2();

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
  TH1F *h_nJets = new TH1F("h_nJets", "Number of Jets; Number of Jets; Events", 20, -0.5, 19.5);h_nJets->Sumw2();

  TH2F *h_Mmumu_MmumuGamma = new TH2F("h_Mmumu_MmumuGamma", "Scatter Plot of  M_{#mu#mu} versus M_{#mu#mu#gamma}; M_{#mu#mu} [GeV]; M_{#mu#mu#gamma} [GeV]", 4000, 0, 2000, 4000, 0, 2000);  h_Mmumu_MmumuGamma->Sumw2();

  TH2F *h_Jet1Pt_PhPt = new TH2F("h_Jet1Pt_PhPt", "Scatter Plot of the leading Jet pT versus Photon pT; pT [GeV]; pT [GeV]", 2000, 0, 1000, 2000, 0, 1000);  h_Jet1Pt_PhPt->Sumw2();

  double Wt_events_mumu = 0;
  double Wt_events_emu = 0;
  double Ztt_events_mumu = 0;
  double Ztt_events_emu = 0;
  double MET_emu = 0;
  double MET_mumu = 0;
  double Ztt_events_LowHT_mumu = 0;
  double Ztt_events_LowHT_emu = 0;
  double LowHT_mumu = 0;
  double LowHT_emu = 0;
  double emu = 0;
  double mumu = 0;
  double Z_peakGammaEvents = 0;

  TFile *triggerCorr1=new TFile("TriggerTurnOn_METParked_Run2012D_All.root");
  TFile *triggerCorr2=new TFile("TriggerTurnOn_Data_ParkedDataFixed_All.root");

  TH2F* Eff_PHO=(TH2F*)triggerCorr1->Get("h_Eff1_2D"); 
  TH2F* Eff_MET=(TH2F*)triggerCorr2->Get("h_Eff2_2D");  

  int nEvents=tree->GetEntries();
  std::cout << "nEvents= " << nEvents << std::endl;
  for (int i=0; i<nEvents ; ++i)
    {
     tree->GetEvent(i);

     if(type=="Data")
       {
       if(fired_HLTPhoIdMet!=1) continue; 
       }

    /*if(type=="MC")
    { 
      if(Znunu==false) continue;
    }
    */
     // filling the photon's properties into a vector of struct
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
     ph1_p4 = fillTLorentzVector(0.0, 0.0, 0.0, 0.0);

     //here working with the leading photon.
     if (photons.size() > 0) ph1_p4=fillTLorentzVector(photons.at(0).pT, photons.at(0).eta, photons.at(0).phi, photons.at(0).energy); 

     double eventWeight_Trigger = 0.0;
     double eventWeight = 0.0;
     if(type=="MC")
     {
       if(ph1_p4.Pt()>30.0 and photons.at(0).isTight==1 and photons.at(0).phIsoTight==1)
       {
         double effective1_MET=MET;
         double effective1_PhpT=ph1_p4.Pt();
         if(effective1_MET > 400.0) effective1_MET=350.0;//inside the bin 
         if(effective1_PhpT > 500.0) effective1_PhpT=450.0;
         double triggerWeight1 = Eff_MET->GetBinContent(Eff_MET->FindBin(effective1_PhpT, effective1_MET));
         double effective2_MET=MET;
         double effective2_PhpT=ph1_p4.Pt();
         if(effective2_MET > 600.0) effective2_MET=550.0;//inside the bin 
         if(effective2_PhpT > 600.0) effective2_PhpT=550.0;
         double triggerWeight2 = Eff_PHO->GetBinContent(Eff_PHO->FindBin(effective2_PhpT, effective2_MET));
         eventWeight_Trigger=triggerWeight1*triggerWeight2;
         eventWeight=triggerWeight1*triggerWeight2*PUWeightData;
       }
     }
     else if(type=="Data")
     {
       if(ph1_p4.Pt()>30.0 and photons.at(0).isTight==1 and photons.at(0).phIsoTight==1) eventWeight = 1.0;
     }
     
     TriggerInfo trigger1; 
     if(type=="Data")
     {
       trigger1.Px = trigObj1Px;
       trigger1.Py = trigObj1Py;
       trigger1.Pz = trigObj1Pz;
       trigger1.E  = trigObj1E;
     } 
     // filling the electron's properties into a vector of struct
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
     el1_p4 = fillTLorentzVector(0.0, 0.0, 0.0, 0.0);
     el2_p4 = fillTLorentzVector(0.0, 0.0, 0.0, 0.0);
     double deltaR_el1 = -1.0;
     if (electrons.size() > 0)
     {
       el1_p4=fillTLorentzVector(electrons.at(0).pT, electrons.at(0).eta, electrons.at(0).phi, electrons.at(0).energy);
       if(type=="Data")
       {
         TLorentzVector trigger1_p4;
         trigger1_p4.SetPxPyPzE(trigger1.Px, trigger1.Py, trigger1.Pz, trigger1.E);
         if(el1_p4.Pt() > 0.0) deltaR_el1 = el1_p4.DeltaR(trigger1_p4);
       }
     }//only execute if an electron exists.

     double deltaR_el2 = -1.0;
     if (electrons.size() > 1)
     {
       el2_p4=fillTLorentzVector(electrons.at(1).pT, electrons.at(1).eta, electrons.at(1).phi, electrons.at(1).energy);
       if(type=="Data")
       {
         TLorentzVector trigger1_p4;
         trigger1_p4.SetPxPyPzE(trigger1.Px, trigger1.Py, trigger1.Pz, trigger1.E);
         if(el2_p4.Pt() > 0.0) deltaR_el2 = el2_p4.DeltaR(trigger1_p4);
       }
     }//only execute if the second electron exists.

    if(type=="MC") deltaR_el1 = 100.0; 
    if(type=="MC") deltaR_el2 = 100.0;

    double leadingDeltaR, trailingDeltaR;
    leadingDeltaR = -1.0;
    trailingDeltaR = -1.0;

    if(el1_p4.Pt()>0.0 and ph1_p4.Pt()>30.0 and photons.at(0).isTight==1 and photons.at(0).phIsoTight==1){
      leadingDeltaR = el1_p4.DeltaR(ph1_p4); 
    }   

    if(el2_p4.Pt()>0.0 and ph1_p4.Pt()>30.0 and photons.at(0).isTight==1 and photons.at(0).phIsoTight==1){
      trailingDeltaR = el2_p4.DeltaR(ph1_p4); 
    }
    // filling the muon's properties into a vector of struct
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
    mu1_p4 = fillTLorentzVector(0.0, 0.0, 0.0, 0.0);
    mu2_p4 = fillTLorentzVector(0.0, 0.0, 0.0, 0.0);

    if (muons.size() > 0) mu1_p4=fillTLorentzVector(muons.at(0).pT, muons.at(0).eta, muons.at(0).phi, muons.at(0).energy);
    if (muons.size() > 1) mu2_p4=fillTLorentzVector(muons.at(1).pT, muons.at(1).eta, muons.at(1).phi, muons.at(1).energy);

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
     if(fabs(jets.at(k).eta)<2.4 and jets.at(k).pT>25.0 and jets.at(k).PU_mva_loose==1)
     {
       Jet.SetPtEtaPhiE(jets.at(k).pT, jets.at(k).eta, jets.at(k).phi, jets.at(k).energy);
       bool isGoodJet=true;
       for(unsigned int j=0; j<electrons.size(); ++j)
       {
        TLorentzVector Electron;
        if(electrons.at(j).isTight==1 and electrons.at(j).isolation < 0.10)
          {
          Electron.SetPtEtaPhiE(electrons.at(j).pT, electrons.at(j).eta, electrons.at(j).phi, electrons.at(j).energy);
          double DRjet_el = Jet.DeltaR(Electron);
          if(DRjet_el<0.5) isGoodJet=false;
          }
        }
   for(unsigned int j=0; j<muons.size(); ++j)
     {
       TLorentzVector Muon;
       if(muons.at(j).isTight==1 and muons.at(j).isolation < 0.12){
       Muon.SetPtEtaPhiE(muons.at(j).pT, muons.at(j).eta, muons.at(j).phi, muons.at(j).energy);
       double DRjet_mu = Jet.DeltaR(Muon);
       if(DRjet_mu<0.5) isGoodJet=false;
     }
   }
   for(unsigned int l=0; l<photons.size(); ++l)
     {
       TLorentzVector Photon;
       if(photons.at(l).pT>30.0 and photons.at(l).isTight==1 and photons.at(l).phIsoTight==1){
         Photon.SetPtEtaPhiE(photons.at(l).pT, photons.at(l).eta, photons.at(l).phi, photons.at(l).energy);
         double DRjet_ph = Jet.DeltaR(Photon);
         if(DRjet_ph<0.5) isGoodJet=false;
     }
   }
   if(type=="Data")
     {
     TLorentzVector trigger1_p4;
     trigger1_p4.SetPxPyPzE(trigger1.Px, trigger1.Py, trigger1.Pz, trigger1.E);
     double DRjet_tr = Jet.DeltaR(trigger1_p4);
     if(DRjet_tr<0.5) isGoodJet=false;
     }

    if(isGoodJet) Jet_vector.push_back(Jet);
     }//close four vector if
  }//close jet loop

  // Now sorting this vector of structs
  std::sort (Jet_vector.begin(), Jet_vector.end(), sortJetVectorsInDescendingpT);

  for(unsigned int m=0; m<Jet_vector.size(); m++) HT += Jet_vector.at(m).Pt();

  bool elEvent = false; 
if(el1_p4.Pt()>10.0 and electrons.at(0).isTight==1 and electrons.at(0).isolation < 0.10 and leadingDeltaR > 0.4 and deltaR_el1 > 0.4){
    elEvent = true;
    if(type=="MC") eventWeight *= electronSF(el1_p4.Pt(), el1_p4.Eta());
    h_el_phi_leading->Fill(el1_p4.Phi(), eventWeight);
    h_el_eta_leading->Fill(el1_p4.Eta(), eventWeight);
    h_el_pt_leading->Fill(el1_p4.Pt(), eventWeight);
    h_el_energy_leading->Fill(el1_p4.E(), eventWeight);
    if(el2_p4.Pt()>10.0 and electrons.at(1).isTight==1 and electrons.at(1).isolation < 0.10 and trailingDeltaR > 0.4 and deltaR_el2 > 0.4 and ((electrons.at(0).charge*electrons.at(1).charge)==-1)) {
      if(type=="MC") eventWeight *= electronSF(el2_p4.Pt(), el2_p4.Eta());
      h_el_phi_trailing->Fill(el2_p4.Phi(), eventWeight);
      h_el_eta_trailing->Fill(el2_p4.Eta(), eventWeight);
      h_el_pt_trailing->Fill(el2_p4.Pt(), eventWeight);
      h_el_energy_trailing->Fill(el2_p4.E(), eventWeight);
      h_InvariantMass_El->Fill((el1_p4+el2_p4).M(), eventWeight);
      if(ph1_p4.Pt()>30.0 and photons.at(0).isTight==1 and photons.at(0).phIsoTight==1) h_InvariantMass_ElPh->Fill((el1_p4+el2_p4+ph1_p4).M(), eventWeight);
      }//trailing electron "if"
  }//closing ledinding electron "if" statement 

  double cMt2_Mu = -99.0;
  double Mlg_Mu = -99.0;
  double pTlg_Mu = -99.0;
  double mt2_Mu = -99.0;
  int mumu_event = 0;
  if(elEvent==false and mu1_p4.Pt()>5.0 and muons.at(0).isTight==1 and muons.at(0).isolation < 0.12 and eventWeight > 0.0){
    if(type=="MC") eventWeight *= muonSF(mu1_p4.Pt(), mu1_p4.Eta());
    h_mu_phi_leading->Fill(mu1_p4.Phi(), eventWeight);
    h_mu_eta_leading->Fill(mu1_p4.Eta(), eventWeight);
    h_mu_pt_leading->Fill(mu1_p4.Pt(), eventWeight);
    h_mu_energy_leading->Fill(mu1_p4.E(), eventWeight);
    double dphi = fabs(MET_Phi - mu1_p4.Phi());
    if (dphi > TMath::Pi()) dphi = TMath::TwoPi() - dphi;
    mt2_Mu = 2*mu1_p4.Pt()*MET*(1 - TMath::Cos(dphi));
    h_Mt_Mu->Fill(TMath::Sqrt(mt2_Mu), eventWeight);
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
      cMt2_Mu=((t1 + MET)*(t1 + MET) - (ph_transverse + mu_transverse + met_transverse).Mod2());
      h_cMt_Mu->Fill(TMath::Sqrt(cMt2_Mu), eventWeight);
      }
    if(mu2_p4.Pt()>5.0 and muons.at(1).isTight==1 and muons.at(1).isolation < 0.12  and ((muons.at(0).charge*muons.at(1).charge)==-1)) {
      if(MET>100.0){
        //if(((mu1_p4+mu2_p4).M()>50 and (mu1_p4+mu2_p4).M()<120)){ 
          if(type=="MC") eventWeight *= muonSF(mu2_p4.Pt(), mu2_p4.Eta());
          h_InvariantMass_Mu_Cut1->Fill((mu1_p4+mu2_p4).M(), eventWeight);
          if(ph1_p4.Pt()>30.0 and photons.at(0).isTight==1 and photons.at(0).phIsoTight==1){
            if(type=="MC") eventWeight *= photonSF(ph1_p4.Pt(), ph1_p4.Eta()); 
            h_MET_Cut0->Fill(MET, eventWeight);
            h_InvariantMass_Mu_Cut2->Fill((mu1_p4+mu2_p4).M(), eventWeight);
            h_InvariantMass_MuPh_Cut2->Fill((mu1_p4+mu2_p4+ph1_p4).M(), eventWeight);
            h_photon_pt_Cut2->Fill(ph1_p4.Pt(), eventWeight); 
            //if(((mu1_p4+mu2_p4+ph1_p4).M() < 80 or (mu1_p4+mu2_p4+ph1_p4).M() > 100)){ 
              h_InvariantMass_Mu_Cut3->Fill((mu1_p4+mu2_p4).M(), eventWeight);
              h_InvariantMass_MuPh_Cut3->Fill((mu1_p4+mu2_p4+ph1_p4).M(), eventWeight);
              h_photon_pt_Cut3->Fill(ph1_p4.Pt(), eventWeight);
              //if((mu1_p4+mu2_p4).M()+(mu1_p4+mu2_p4+ph1_p4).M() > 185.0){
              //if(type=="MC") eventWeight *= muonSF(mu2_p4.Pt(), mu2_p4.Eta())*photonSF(ph1_p4.Pt(), ph1_p4.Eta());  
                h_InvariantMass_Mu_Cut4->Fill((mu1_p4+mu2_p4).M(), eventWeight);
                h_InvariantMass_MuPh_Cut4->Fill((mu1_p4+mu2_p4+ph1_p4).M(), eventWeight);
                h_photon_pt_Cut4->Fill(ph1_p4.Pt(), eventWeight);
                mumu_event = 1;
                mumu+=eventWeight;
                h_mu_phi_trailing->Fill(mu2_p4.Phi(), eventWeight);
                h_mu_eta_trailing->Fill(mu2_p4.Eta(), eventWeight);
                h_mu_pt_trailing->Fill(mu2_p4.Pt(), eventWeight);
                h_mu_energy_trailing->Fill(mu2_p4.E(), eventWeight);
                h_InvariantMass_Mu->Fill((mu1_p4+mu2_p4).M(), eventWeight);
                h_InvariantMass_MuPh->Fill((mu1_p4+mu2_p4+ph1_p4).M(), eventWeight);
                h_Mmumu_MmumuGamma->Fill((mu1_p4+mu2_p4).M(), (mu1_p4+mu2_p4+ph1_p4).M(), eventWeight);
                if(((mu1_p4+mu2_p4+ph1_p4).M()>76 and (mu1_p4+mu2_p4+ph1_p4).M()<106)) Z_peakGammaEvents+=eventWeight;
                if(type=="MC" and Ztt==true) Ztt_events_mumu+=eventWeight;
                if(type=="MC" and Wt==true) Wt_events_mumu+=eventWeight;
                if(type=="Data" and ((mu1_p4+mu2_p4).M() > 10 and (mu1_p4+mu2_p4).M() < 60)) Ztt_events_mumu+=eventWeight;
                if(MET > 25.0) MET_mumu+=eventWeight;
                if(type=="MC" and Ztt==true and HT<300) Ztt_events_LowHT_mumu+=eventWeight; 
                if(type=="Data" and (((mu1_p4+mu2_p4).M() > 10 and (mu1_p4+mu2_p4).M() < 60) and HT<300)) Ztt_events_LowHT_mumu+=eventWeight;
                if(HT<300) LowHT_mumu+=eventWeight;
                //}//extra event selection  
             //}//Z selection 
            }//photon requirement
          //}//Zgamma event selection
        }//MET cut
     }//second muon  
  }//first muon

   int emu_event = 0;
   if((mu1_p4.Pt()>0.0 and muons.at(0).isTight==1 and muons.at(0).isolation < 0.12) and (el1_p4.Pt()>0.0 and electrons.at(0).isTight==1 and electrons.at(0).isolation < 0.10 and leadingDeltaR > 0.4 and deltaR_el1 > 0.4) and eventWeight > 0.0 and mumu_event == 0 ) {
     if(type=="MC") eventWeight *= electronSF(el1_p4.Pt(), el1_p4.Eta())*muonSF(mu1_p4.Pt(), mu1_p4.Eta());  
     emu_event = 1;
     emu+=eventWeight;
     h_InvariantMass_ElMu->Fill((mu1_p4+el1_p4).M(), eventWeight);
     if(type=="MC" and Ztt==1) Ztt_events_emu+=eventWeight;
     if(type=="MC" and Wt==true)  Wt_events_emu+=eventWeight;
     if(type=="Data" and ((mu1_p4+el1_p4).M() > 10 and (mu1_p4+el1_p4).M() < 60)) Ztt_events_emu+=eventWeight;
     if(MET > 25.0) MET_emu+=eventWeight;
     if(type=="MC" and Ztt==true and HT<300) Ztt_events_LowHT_emu+=eventWeight;
     if(type=="Data" and (((mu1_p4+el1_p4).M() > 10 and (mu1_p4+el1_p4).M() < 60) and HT<300)) Ztt_events_LowHT_emu+=eventWeight;
     if(HT<300) LowHT_emu+=eventWeight;
     
     }//El-Mu invariant mass

  if(channel=="MuMu" and mumu_event==1 and emu_event == 0){

     if (ph1_p4.Pt()>30.0 and photons.at(0).isTight==1 and photons.at(0).phIsoTight==1 and eventWeight > 0.0){
       h_photon_pt->Fill(ph1_p4.Pt(), eventWeight);
       h_photon_eta->Fill(ph1_p4.Eta(), eventWeight);
       h_photon_phi->Fill(ph1_p4.Phi(), eventWeight);
       h_photon_energy->Fill(ph1_p4.E(), eventWeight);
       if(Jet_vector.size()>0) h_Jet1Pt_PhPt->Fill(Jet_vector.at(0).Pt(), ph1_p4.Pt(), eventWeight);
    }

    if(eventWeight > 0.0) h_MET->Fill(MET, eventWeight);
    if(eventWeight > 0.0) h_nVertices->Fill(nVertices, eventWeight);
    if(eventWeight > 0.0) h_nPUVertices->Fill(nPUVertices, eventWeight);
    if(eventWeight > 0.0) h_nPUVerticesTrue->Fill(nPUVerticesTrue, eventWeight);

    if(Jet_vector.size()>0 and eventWeight > 0.0) h_jet_pt_leading->Fill(Jet_vector.at(0).Pt(), eventWeight);
    if(Jet_vector.size()>1 and eventWeight > 0.0) h_jet_pt_trailing->Fill(Jet_vector.at(1).Pt(), eventWeight);
    if(Jet_vector.size()>2 and eventWeight > 0.0) h_jet_pt_3rd->Fill(Jet_vector.at(2).Pt(), eventWeight);
    if(Jet_vector.size()>3 and eventWeight > 0.0) h_jet_pt_4th->Fill(Jet_vector.at(3).Pt(), eventWeight);
    if(Jet_vector.size()>4 and eventWeight > 0.0) h_jet_pt_5th->Fill(Jet_vector.at(4).Pt(), eventWeight);
    if(Jet_vector.size()>5 and eventWeight > 0.0) h_jet_pt_6th->Fill(Jet_vector.at(5).Pt(), eventWeight);

    if(Jet_vector.size()>0 and eventWeight > 0.0) h_jet_eta_leading->Fill(Jet_vector.at(0).Eta(), eventWeight);
    if(Jet_vector.size()>1 and eventWeight > 0.0) h_jet_eta_trailing->Fill(Jet_vector.at(1).Eta(), eventWeight);
    if(Jet_vector.size()>2 and eventWeight > 0.0) h_jet_eta_3rd->Fill(Jet_vector.at(2).Eta(), eventWeight);
    if(Jet_vector.size()>3 and eventWeight > 0.0) h_jet_eta_4th->Fill(Jet_vector.at(3).Eta(), eventWeight);
    if(Jet_vector.size()>4 and eventWeight > 0.0) h_jet_eta_5th->Fill(Jet_vector.at(4).Eta(), eventWeight);
    if(Jet_vector.size()>5 and eventWeight > 0.0) h_jet_eta_6th->Fill(Jet_vector.at(5).Eta(), eventWeight);

    if(Jet_vector.size()>0 and eventWeight > 0.0) h_jet_phi_leading->Fill(Jet_vector.at(0).Phi(), eventWeight);
    if(Jet_vector.size()>1 and eventWeight > 0.0) h_jet_phi_trailing->Fill(Jet_vector.at(1).Phi(), eventWeight);
    if(Jet_vector.size()>2 and eventWeight > 0.0) h_jet_phi_3rd->Fill(Jet_vector.at(2).Phi(), eventWeight);
    if(Jet_vector.size()>3 and eventWeight > 0.0) h_jet_phi_4th->Fill(Jet_vector.at(3).Phi(), eventWeight);
    if(Jet_vector.size()>4 and eventWeight > 0.0) h_jet_phi_5th->Fill(Jet_vector.at(4).Phi(), eventWeight);
    if(Jet_vector.size()>5 and eventWeight > 0.0) h_jet_phi_6th->Fill(Jet_vector.at(5).Phi(), eventWeight);

    if(Jet_vector.size()>0 and eventWeight > 0.0) h_jet_energy_leading->Fill(Jet_vector.at(0).E(), eventWeight);
    if(Jet_vector.size()>1 and eventWeight > 0.0) h_jet_energy_trailing->Fill(Jet_vector.at(1).E(), eventWeight);
    if(Jet_vector.size()>2 and eventWeight > 0.0) h_jet_energy_3rd->Fill(Jet_vector.at(2).E(), eventWeight);
    if(Jet_vector.size()>3 and eventWeight > 0.0) h_jet_energy_4th->Fill(Jet_vector.at(3).E(), eventWeight);
    if(Jet_vector.size()>4 and eventWeight > 0.0) h_jet_energy_5th->Fill(Jet_vector.at(4).E(), eventWeight);
    if(Jet_vector.size()>5 and eventWeight > 0.0) h_jet_energy_6th->Fill(Jet_vector.at(5).E(), eventWeight);

    if(eventWeight > 0.0) h_nJets->Fill(Jet_vector.size(), eventWeight);
    if(eventWeight > 0.0) h_HT->Fill(HT, eventWeight);
  }

  if(channel=="ElMu" and emu_event==1 and mumu_event==0){
 
    if (ph1_p4.Pt()>30.0 and photons.at(0).isTight==1 and photons.at(0).phIsoTight==1 and eventWeight > 0.0){
       h_photon_pt->Fill(ph1_p4.Pt(), eventWeight);
       h_photon_eta->Fill(ph1_p4.Eta(), eventWeight);
       h_photon_phi->Fill(ph1_p4.Phi(), eventWeight);
       h_photon_energy->Fill(ph1_p4.E(), eventWeight);
       if(Jet_vector.size()>0) h_Jet1Pt_PhPt->Fill(Jet_vector.at(0).Pt(), ph1_p4.Pt(), eventWeight); 
    }

    if(eventWeight > 0.0) h_MET->Fill(MET, eventWeight);
    if(eventWeight > 0.0) h_nVertices->Fill(nVertices, eventWeight);
    if(eventWeight > 0.0) h_nPUVertices->Fill(nPUVertices, eventWeight);
    if(eventWeight > 0.0) h_nPUVerticesTrue->Fill(nPUVerticesTrue, eventWeight);

    if(Jet_vector.size()>0 and eventWeight > 0.0) h_jet_pt_leading->Fill(Jet_vector.at(0).Pt(), eventWeight);
    if(Jet_vector.size()>1 and eventWeight > 0.0) h_jet_pt_trailing->Fill(Jet_vector.at(1).Pt(), eventWeight);
    if(Jet_vector.size()>2 and eventWeight > 0.0) h_jet_pt_3rd->Fill(Jet_vector.at(2).Pt(), eventWeight);
    if(Jet_vector.size()>3 and eventWeight > 0.0) h_jet_pt_4th->Fill(Jet_vector.at(3).Pt(), eventWeight);
    if(Jet_vector.size()>4 and eventWeight > 0.0) h_jet_pt_5th->Fill(Jet_vector.at(4).Pt(), eventWeight);
    if(Jet_vector.size()>5 and eventWeight > 0.0) h_jet_pt_6th->Fill(Jet_vector.at(5).Pt(), eventWeight);

    if(Jet_vector.size()>0 and eventWeight > 0.0) h_jet_eta_leading->Fill(Jet_vector.at(0).Eta(), eventWeight);
    if(Jet_vector.size()>1 and eventWeight > 0.0) h_jet_eta_trailing->Fill(Jet_vector.at(1).Eta(), eventWeight);
    if(Jet_vector.size()>2 and eventWeight > 0.0) h_jet_eta_3rd->Fill(Jet_vector.at(2).Eta(), eventWeight);
    if(Jet_vector.size()>3 and eventWeight > 0.0) h_jet_eta_4th->Fill(Jet_vector.at(3).Eta(), eventWeight);
    if(Jet_vector.size()>4 and eventWeight > 0.0) h_jet_eta_5th->Fill(Jet_vector.at(4).Eta(), eventWeight);
    if(Jet_vector.size()>5 and eventWeight > 0.0) h_jet_eta_6th->Fill(Jet_vector.at(5).Eta(), eventWeight);

    if(Jet_vector.size()>0 and eventWeight > 0.0) h_jet_phi_leading->Fill(Jet_vector.at(0).Phi(), eventWeight);
    if(Jet_vector.size()>1 and eventWeight > 0.0) h_jet_phi_trailing->Fill(Jet_vector.at(1).Phi(), eventWeight);
    if(Jet_vector.size()>2 and eventWeight > 0.0) h_jet_phi_3rd->Fill(Jet_vector.at(2).Phi(), eventWeight);
    if(Jet_vector.size()>3 and eventWeight > 0.0) h_jet_phi_4th->Fill(Jet_vector.at(3).Phi(), eventWeight);
    if(Jet_vector.size()>4 and eventWeight > 0.0) h_jet_phi_5th->Fill(Jet_vector.at(4).Phi(), eventWeight);
    if(Jet_vector.size()>5 and eventWeight > 0.0) h_jet_phi_6th->Fill(Jet_vector.at(5).Phi(), eventWeight);

    if(Jet_vector.size()>0 and eventWeight > 0.0) h_jet_energy_leading->Fill(Jet_vector.at(0).E(), eventWeight);
    if(Jet_vector.size()>1 and eventWeight > 0.0) h_jet_energy_trailing->Fill(Jet_vector.at(1).E(), eventWeight);
    if(Jet_vector.size()>2 and eventWeight > 0.0) h_jet_energy_3rd->Fill(Jet_vector.at(2).E(), eventWeight);
    if(Jet_vector.size()>3 and eventWeight > 0.0) h_jet_energy_4th->Fill(Jet_vector.at(3).E(), eventWeight);
    if(Jet_vector.size()>4 and eventWeight > 0.0) h_jet_energy_5th->Fill(Jet_vector.at(4).E(), eventWeight);
    if(Jet_vector.size()>5 and eventWeight > 0.0) h_jet_energy_6th->Fill(Jet_vector.at(5).E(), eventWeight);

    if(eventWeight > 0.0) h_nJets->Fill(Jet_vector.size(), eventWeight);
    if(eventWeight > 0.0) h_HT->Fill(HT, eventWeight);   
    }

  }//event loop closed

  //Cut flow table
  cout << "mumu = " << mumu << endl;
  cout << "emu = " << emu << endl;
  cout << "Ztt_events_mumu = " << Ztt_events_mumu << endl;
  cout << "Ztt_events_emu = " << Ztt_events_emu << endl;
  cout << "Wt_events_mumu = " << Wt_events_mumu << endl;
  cout << "Wt_events_emu = " << Wt_events_emu << endl;
  cout << "MET_emu = " << MET_emu << endl;
  cout << "MET_mumu = " << MET_mumu << endl;
  cout << "Ztt_events_LowHT_mumu = " << Ztt_events_LowHT_mumu << endl;
  cout << "Ztt_events_LowHT_emu = " << Ztt_events_LowHT_emu << endl;
  cout << "LowHT_mumu = " << LowHT_mumu << endl;
  cout << "LowHT_emu = " << LowHT_emu << endl;
  cout << "Z_peakGammaEvents = " << Z_peakGammaEvents << endl;
  std::string histfilename=(outfile+".root").c_str();
  TFile *tFile=new TFile(histfilename.c_str(), "RECREATE");
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
  h_HT->Write();
  h_nVertices->Write();
  h_nPUVertices->Write();
  h_nPUVerticesTrue->Write();

  h_cMt_Mu->Write();
  h_Mt_Mu->Write();
  h_MET->Write();
  h_MET_Cut0->Write();
  h_mu_pt_leading->Write();
  h_mu_pt_trailing->Write();
  h_mu_eta_leading->Write();
  h_mu_eta_trailing->Write();
  h_mu_phi_leading->Write();
  h_mu_phi_trailing->Write();
  h_mu_energy_leading->Write();
  h_mu_energy_trailing->Write();
 
  h_el_pt_leading->Write();
  h_el_pt_trailing->Write();
  h_el_eta_leading->Write();
  h_el_eta_trailing->Write();
  h_el_phi_leading->Write();
  h_el_phi_trailing->Write();
  h_el_energy_leading->Write();
  h_el_energy_trailing->Write();

  h_photon_pt->Write();
  h_photon_pt_Cut2->Write();
  h_photon_pt_Cut3->Write();
  h_photon_pt_Cut4->Write();
  h_photon_eta->Write();
  h_photon_phi->Write();
  h_photon_energy->Write();

  h_InvariantMass_Mu->Write();
  h_InvariantMass_Mu_Cut1->Write();
  h_InvariantMass_Mu_Cut2->Write();
  h_InvariantMass_Mu_Cut3->Write();
  h_InvariantMass_Mu_Cut4->Write();
  h_InvariantMass_El->Write();
  h_InvariantMass_ElPh->Write();
  h_InvariantMass_MuPh->Write();
  h_InvariantMass_MuPh_Cut2->Write();
  h_InvariantMass_MuPh_Cut3->Write();
  h_InvariantMass_MuPh_Cut4->Write();
  h_InvariantMass_ElMu->Write();
  h_Mmumu_MmumuGamma->Write();
  h_Jet1Pt_PhPt->Write();

  h_HT->Write();
  tFile->Close(); 
  std::cout<<"Wrote output file "<<histfilename<<std::endl;

  return 0;

}
