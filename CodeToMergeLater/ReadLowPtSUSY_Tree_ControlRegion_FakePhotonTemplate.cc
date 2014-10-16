#include "ReadLowPtSUSY_Tree_ControlRegion_FakePhotonTemplate.h"

int ReadLowPtSUSY_Tree_ControlRegion_FakePhotonTemplate(std::string infile,  std::string outfile,  std::string type){
  
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
  vector<float>   *ph_HoE;
  vector<int>     *ph_conversionVeto;
  vector<int>     *ph_pixelVeto;
  vector<float>   *ph_SigmaIetaIeta;
  vector<float>   *ph_SigmaIetaIphi;
  vector<float>   *ph_SigmaIphiIphi;
  vector<float>   *ph_preShowerOverRaw;
  vector<float>   *ph_R9;
  Int_t           nPhotons;
  vector<float>   *el_pt;
  vector<float>   *el_phi;
  vector<float>   *el_eta;
  vector<float>   *el_energy;
  vector<int>     *el_charge;
  vector<bool>    *el_isTight;
  vector<bool>    *el_isLoose;
  Int_t           nElectrons;
  vector<float>   *mu_pt;
  vector<float>   *mu_phi;
  vector<float>   *mu_eta;
  vector<float>   *mu_energy;
  vector<int>     *mu_charge;
  vector<bool>    *mu_isTight;
  vector<bool>    *mu_isLoose;
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
  Float_t         caloMET;
  Float_t         caloMET_Phi;
  Float_t         caloMET_Px;
  Float_t         caloMET_Py;
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
  vector<float>   *ph_e1x5;
  vector<float>   *ph_e1x3;
  vector<float>   *ph_e2x2;
  vector<float>   *ph_e2x5;
  vector<float>   *ph_e5x1;
  vector<float>   *ph_e5x5;
  vector<float>   *ph_e2x5Max;
  vector<float>   *ph_e2OverE5;
  vector<float>   *ph_seedCrystalEnergy;
  vector<int>     *ph_Matched;
  vector<float>   *ph_MatchedPt;
  vector<float>   *ph_MatchedEta;
  vector<float>   *ph_MatchedPhi;
  vector<float>   *ph_MatchedEnergy;
  vector<float>   *ph_MatchedSigmaIetaIeta;
  vector<bool>   *ph_MatchedisTightWS;
  vector<bool>   *ph_MatchedIsoTight;
  vector<int>   *ph_MatchedpixelVeto;
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
  float trigObj3Px;
  float trigObj3Py;
  float trigObj3Pz;
  float trigObj3E;

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
  el_isLoose = 0;
  mu_pt = 0;
  mu_phi = 0;
  mu_eta = 0;
  mu_energy = 0;
  mu_charge = 0;
  mu_isTight = 0;
  mu_isLoose = 0;
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
  ph_HoE = 0;
  ph_conversionVeto = 0;
  ph_pixelVeto = 0;
  ph_SigmaIetaIeta = 0;
  ph_SigmaIetaIphi = 0;
  ph_SigmaIphiIphi = 0;
  ph_preShowerOverRaw = 0;
  ph_R9 = 0; 
  ph_R9 = 0;
  ph_e1x5 = 0;
  ph_e1x3 = 0;
  ph_e2x2 = 0;
  ph_e2x5 = 0;
  ph_e5x1 = 0;
  ph_e5x5 = 0;
  ph_e2x5Max = 0;
  ph_e2OverE5 = 0;
  ph_seedCrystalEnergy = 0;
  ph_Matched = 0;
  ph_MatchedPt = 0;
  ph_MatchedEta = 0;
  ph_MatchedPhi = 0;
  ph_MatchedEnergy = 0;
  ph_MatchedSigmaIetaIeta = 0;
  ph_MatchedisTightWS = 0;
  ph_MatchedIsoTight = 0;
  ph_MatchedpixelVeto = 0;
 
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
  tree->SetBranchAddress("ph_HoE", &(ph_HoE));
  tree->SetBranchAddress("ph_conversionVeto", &(ph_conversionVeto));
  tree->SetBranchAddress("ph_pixelVeto", &(ph_pixelVeto));
  tree->SetBranchAddress("ph_SigmaIetaIeta", &(ph_SigmaIetaIeta));
  tree->SetBranchAddress("ph_SigmaIetaIphi", &(ph_SigmaIetaIphi));
  tree->SetBranchAddress("ph_SigmaIphiIphi", &(ph_SigmaIphiIphi));
  tree->SetBranchAddress("ph_preShowerOverRaw", &(ph_preShowerOverRaw));
  tree->SetBranchAddress("ph_R9", &(ph_R9));
  tree->SetBranchAddress("ph_e1x5", &(ph_e1x5));
  tree->SetBranchAddress("ph_e1x3", &(ph_e1x3));
  tree->SetBranchAddress("ph_e2x2", &(ph_e2x2));
  tree->SetBranchAddress("ph_e2x5", &(ph_e2x5));
  tree->SetBranchAddress("ph_e5x1", &(ph_e5x1));
  tree->SetBranchAddress("ph_e5x5", &(ph_e5x5));
  tree->SetBranchAddress("ph_e2x5Max", &(ph_e2x5Max));
  tree->SetBranchAddress("ph_e2OverE5", &(ph_e2OverE5));
  tree->SetBranchAddress("ph_seedCrystalEnergy", &(ph_seedCrystalEnergy));
  tree->SetBranchAddress("el_pt", &(el_pt));
  tree->SetBranchAddress("el_eta", &(el_eta));
  tree->SetBranchAddress("el_phi", &(el_phi));
  tree->SetBranchAddress("el_energy", &(el_energy));
  tree->SetBranchAddress("el_charge", &(el_charge));
  tree->SetBranchAddress("el_isTight", &(el_isTight)); 
  tree->SetBranchAddress("el_isLoose", &(el_isLoose)); 
  tree->SetBranchAddress("el_iso", &(el_iso));
  tree->SetBranchAddress("nElectrons", &(nElectrons));
  tree->SetBranchAddress("mu_pt", &(mu_pt));
  tree->SetBranchAddress("mu_eta", &(mu_eta));
  tree->SetBranchAddress("mu_phi", &(mu_phi));
  tree->SetBranchAddress("mu_energy", &(mu_energy)); 
  tree->SetBranchAddress("mu_charge", &(mu_charge)); 
  tree->SetBranchAddress("mu_isTight", &(mu_isTight));
  tree->SetBranchAddress("mu_isLoose", &(mu_isLoose)); 
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
  tree->SetBranchAddress("caloMET", &(caloMET));
  tree->SetBranchAddress("caloMET_Phi", &(caloMET_Phi));
  tree->SetBranchAddress("caloMET_Px", &(caloMET_Px));
  tree->SetBranchAddress("caloMET_Py", &(caloMET_Py));
  tree->SetBranchAddress("nVertices", &(nVertices));
  tree->SetBranchAddress("nPUVertices", &(nPUVertices));
  tree->SetBranchAddress("nPUVerticesTrue", &(nPUVerticesTrue));
  tree->SetBranchAddress("PUWeightData", &(PUWeightData));
  tree->SetBranchAddress("PUWeightDataSys", &(PUWeightDataSys)); 
  tree->SetBranchAddress("ph_phIsoTight", &(ph_phIsoTight));
  tree->SetBranchAddress("ph_phIsoMedium", &(ph_phIsoMedium));
  tree->SetBranchAddress("ph_phIsoLoose", &(ph_phIsoLoose));
  tree->SetBranchAddress("fired_HLTPho", &(fired_HLTPho));
  tree->SetBranchAddress("fired_HLTPhoId", &(fired_HLTPhoId));
  tree->SetBranchAddress("fired_HLTPhoIdMet", &(fired_HLTPhoIdMet));
  if(type=="Data"){
    tree->SetBranchAddress("trigObj1Px", &(trigObj1Px));
    tree->SetBranchAddress("trigObj1Py", &(trigObj1Py));
    tree->SetBranchAddress("trigObj1Pz", &(trigObj1Pz));
    tree->SetBranchAddress("trigObj1E", &(trigObj1E));
    tree->SetBranchAddress("trigObj2Px", &(trigObj2Px));
    tree->SetBranchAddress("trigObj2Py", &(trigObj2Py));
    tree->SetBranchAddress("trigObj2Pz", &(trigObj2Pz));
    tree->SetBranchAddress("trigObj2E", &(trigObj2E));
    tree->SetBranchAddress("trigObj3Px", &(trigObj3Px));
    tree->SetBranchAddress("trigObj3Py", &(trigObj3Py));
    tree->SetBranchAddress("trigObj3Pz", &(trigObj3Pz));
    tree->SetBranchAddress("trigObj3E", &(trigObj3E));
  }
  if(type=="MC"){
    tree->SetBranchAddress("ph_Matched", &(ph_Matched));
    tree->SetBranchAddress("ph_MatchedPt", &(ph_MatchedPt));
    tree->SetBranchAddress("ph_MatchedEta", &(ph_MatchedEta));
    tree->SetBranchAddress("ph_MatchedPhi", &(ph_MatchedPhi));
    tree->SetBranchAddress("ph_MatchedEnergy", &(ph_MatchedEnergy));
    tree->SetBranchAddress("ph_MatchedSigmaIetaIeta", &(ph_MatchedSigmaIetaIeta));
    tree->SetBranchAddress("ph_MatchedisTightWS", &(ph_MatchedisTightWS));
    tree->SetBranchAddress("ph_MatchedIsoTight", &(ph_MatchedIsoTight));
    tree->SetBranchAddress("ph_MatchedpixelVeto", &(ph_MatchedpixelVeto));
  }

  TH1F *h_SigmaIetaIeta_Num_Ph30To40 =new TH1F("h_SigmaIetaIeta_Num_Ph30To40", "Photon SigmaIetaIeta; SigmaIetaIeta; Events", 3000, 0, 0.03); h_SigmaIetaIeta_Num_Ph30To40->Sumw2();
  TH1F *h_SigmaIetaIeta_SB_Ph30To40 =new TH1F("h_SigmaIetaIeta_SB_Ph30To40", "Photon SigmaIetaIeta; SigmaIetaIeta; Events", 3000, 0, 0.03); h_SigmaIetaIeta_SB_Ph30To40->Sumw2();
  TH1F *h_SigmaIetaIeta_Matched_Ph30To40 =new TH1F("h_SigmaIetaIeta_Matched_Ph30To40", "Photon SigmaIetaIeta; SigmaIetaIeta; Events", 3000, 0, 0.03); h_SigmaIetaIeta_Matched_Ph30To40->Sumw2();
  TH1F *h_SigmaIetaIeta_Num_Ph40To50 =new TH1F("h_SigmaIetaIeta_Num_Ph40To50", "Photon SigmaIetaIeta; SigmaIetaIeta; Events", 3000, 0, 0.03); h_SigmaIetaIeta_Num_Ph40To50->Sumw2();
  TH1F *h_SigmaIetaIeta_SB_Ph40To50 =new TH1F("h_SigmaIetaIeta_SB_Ph40To50", "Photon SigmaIetaIeta; SigmaIetaIeta; Events", 3000, 0, 0.03); h_SigmaIetaIeta_SB_Ph40To50->Sumw2();
  TH1F *h_SigmaIetaIeta_Matched_Ph40To50 =new TH1F("h_SigmaIetaIeta_Matched_Ph40To50", "Photon SigmaIetaIeta; SigmaIetaIeta; Events", 3000, 0, 0.03); h_SigmaIetaIeta_Matched_Ph40To50->Sumw2();
  TH1F *h_SigmaIetaIeta_Num_Ph50To60 =new TH1F("h_SigmaIetaIeta_Num_Ph50To60", "Photon SigmaIetaIeta; SigmaIetaIeta; Events", 3000, 0, 0.03); h_SigmaIetaIeta_Num_Ph50To60->Sumw2();
  TH1F *h_SigmaIetaIeta_SB_Ph50To60 =new TH1F("h_SigmaIetaIeta_SB_Ph50To60", "Photon SigmaIetaIeta; SigmaIetaIeta; Events", 3000, 0, 0.03); h_SigmaIetaIeta_SB_Ph50To60->Sumw2();
  TH1F *h_SigmaIetaIeta_Matched_Ph50To60 =new TH1F("h_SigmaIetaIeta_Matched_Ph50To60", "Photon SigmaIetaIeta; SigmaIetaIeta; Events", 3000, 0, 0.03); h_SigmaIetaIeta_Matched_Ph50To60->Sumw2();
  TH1F *h_SigmaIetaIeta_Num_Ph60To70 =new TH1F("h_SigmaIetaIeta_Num_Ph60To70", "Photon SigmaIetaIeta; SigmaIetaIeta; Events", 3000, 0, 0.03); h_SigmaIetaIeta_Num_Ph60To70->Sumw2();
  TH1F *h_SigmaIetaIeta_SB_Ph60To70 =new TH1F("h_SigmaIetaIeta_SB_Ph60To70", "Photon SigmaIetaIeta; SigmaIetaIeta; Events", 3000, 0, 0.03); h_SigmaIetaIeta_SB_Ph60To70->Sumw2();
  TH1F *h_SigmaIetaIeta_Matched_Ph60To70 =new TH1F("h_SigmaIetaIeta_Matched_Ph60To70", "Photon SigmaIetaIeta; SigmaIetaIeta; Events", 3000, 0, 0.03); h_SigmaIetaIeta_Matched_Ph60To70->Sumw2();
  TH1F *h_SigmaIetaIeta_Num_Ph70To80 =new TH1F("h_SigmaIetaIeta_Num_Ph70To80", "Photon SigmaIetaIeta; SigmaIetaIeta; Events", 3000, 0, 0.03); h_SigmaIetaIeta_Num_Ph70To80->Sumw2();
  TH1F *h_SigmaIetaIeta_SB_Ph70To80 =new TH1F("h_SigmaIetaIeta_SB_Ph70To80", "Photon SigmaIetaIeta; SigmaIetaIeta; Events", 3000, 0, 0.03); h_SigmaIetaIeta_SB_Ph70To80->Sumw2();
  TH1F *h_SigmaIetaIeta_Matched_Ph70To80 =new TH1F("h_SigmaIetaIeta_Matched_Ph70To80", "Photon SigmaIetaIeta; SigmaIetaIeta; Events", 3000, 0, 0.03); h_SigmaIetaIeta_Matched_Ph70To80->Sumw2();
  TH1F *h_SigmaIetaIeta_Num_Ph80To100 =new TH1F("h_SigmaIetaIeta_Num_Ph80To100", "Photon SigmaIetaIeta; SigmaIetaIeta; Events", 3000, 0, 0.03); h_SigmaIetaIeta_Num_Ph80To100->Sumw2();
  TH1F *h_SigmaIetaIeta_SB_Ph80To100 =new TH1F("h_SigmaIetaIeta_SB_Ph80To100", "Photon SigmaIetaIeta; SigmaIetaIeta; Events", 3000, 0, 0.03); h_SigmaIetaIeta_SB_Ph80To100->Sumw2();
  TH1F *h_SigmaIetaIeta_Matched_Ph80To100 =new TH1F("h_SigmaIetaIeta_Matched_Ph80To100", "Photon SigmaIetaIeta; SigmaIetaIeta; Events", 3000, 0, 0.03); h_SigmaIetaIeta_Matched_Ph80To100->Sumw2();
  TH1F *h_SigmaIetaIeta_Num_Ph100To130 =new TH1F("h_SigmaIetaIeta_Num_Ph100To130", "Photon SigmaIetaIeta; SigmaIetaIeta; Events", 3000, 0, 0.03); h_SigmaIetaIeta_Num_Ph100To130->Sumw2();
  TH1F *h_SigmaIetaIeta_SB_Ph100To130 =new TH1F("h_SigmaIetaIeta_SB_Ph100To130", "Photon SigmaIetaIeta; SigmaIetaIeta; Events", 3000, 0, 0.03); h_SigmaIetaIeta_SB_Ph100To130->Sumw2();
  TH1F *h_SigmaIetaIeta_Matched_Ph100To130 =new TH1F("h_SigmaIetaIeta_Matched_Ph100To130", "Photon SigmaIetaIeta; SigmaIetaIeta; Events", 3000, 0, 0.03); h_SigmaIetaIeta_Matched_Ph100To130->Sumw2();
  
  int nEvents=tree->GetEntries();
  std::cout << "nEvents= " << nEvents << std::endl;
  for (int i=0; i<nEvents ; ++i)
    {
     tree->GetEvent(i);

     if(type=="Data")
       {
       if(fired_HLTPhoIdMet!=1) continue;
       }
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
       photon.phHoE = ph_HoE->at(j);
       photon.phconversionVeto = ph_conversionVeto->at(j);
       photon.phpixelVeto = ph_pixelVeto->at(j);
       photon.phSigmaIetaIeta = ph_SigmaIetaIeta->at(j);
       photon.phSigmaIetaIphi = ph_SigmaIetaIphi->at(j);
       photon.phSigmaIphiIphi = ph_SigmaIphiIphi->at(j);
       photon.phpreShowerOverRaw = ph_preShowerOverRaw->at(j);
       photon.phR9 = ph_R9->at(j);    
       photon.phe1x5 = ph_e1x5->at(j);
       photon.phe1x3 = ph_e1x3->at(j);
       photon.phe2x2 = ph_e2x2->at(j);
       photon.phe2x5 = ph_e2x5->at(j);
       photon.phe5x1 = ph_e5x1->at(j);
       photon.phe5x5 = ph_e5x5->at(j);
       photon.phe2x5Max = ph_e2x5Max->at(j);
       photon.phe2OverE5 = ph_e2OverE5->at(j);
       photon.phseedCrystalEnergy = ph_seedCrystalEnergy->at(j);
       if(type=="MC") photon.Matched = ph_Matched->at(j);
       photons.push_back(photon);
       }
     // Now sorting this vector of structs
       
     TriggerInfo trigger1; 
     TriggerInfo trigger3;
     
     if(type=="Data")
     {
       trigger1.Px = trigObj1Px;
       trigger1.Py = trigObj1Py;
       trigger1.Pz = trigObj1Pz;
       trigger1.E  = trigObj1E;
       trigger3.Px = trigObj3Px;
       trigger3.Py = trigObj3Py;
       trigger3.Pz = trigObj3Pz;
       trigger3.E  = trigObj3E;
     }
     
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
       electron.isLoose=el_isLoose->at(j);
       electron.isolation=el_iso->at(j);
       electrons.push_back(electron);
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
       muon.isLoose=mu_isLoose->at(j);
       muon.isolation=mu_iso->at(j);
       muons.push_back(muon);
    } 
  
   std::vector<AnalysisPhotonInfo> anaPhotons;
   std::vector<float> swissCross;
   swissCross.clear();
   for(unsigned int k=0; k<photons.size(); ++k)
   {
     AnalysisPhotonInfo anaPhoton; //analysis photon without sigma ietaieta
 
     swissCross.push_back((photons.at(k).phseedCrystalEnergy)/(photons.at(k).phe1x5 + photons.at(k).phe5x1 - photons.at(k).phseedCrystalEnergy));

     if(fabs(photons.at(k).eta)<1.444 and photons.at(k).pT>30.0 and photons.at(k).phHoE < 0.05 and photons.at(0).phR9 > 0.9 and photons.at(k).phIsoTight==1 and photons.at(k).phpixelVeto==0 and swissCross.at(k)<0.90 and photons.at(k).phR9 < 1.0 and photons.at(k).phSigmaIetaIeta > 0.001 and photons.at(k).phSigmaIphiIphi>0.0001) //tight ID except sieie cut and swiss cross cut and additional cuts for beam halo rejection
     {
       anaPhoton.PhLV.SetPtEtaPhiE(photons.at(k).pT, photons.at(k).eta, photons.at(k).phi, photons.at(k).energy);
       anaPhoton.SigmaIetaIeta = ph_SigmaIetaIeta->at(k);
       if(type=="Data") anaPhoton.matched = 0;
       if(type=="MC") anaPhoton.matched = ph_Matched->at(k);
       bool isGoodPhoton=true;
       for(unsigned int j=0; j<electrons.size(); ++j)
         {
         TLorentzVector Electron;
         if(electrons.at(j).isTight==1 and electrons.at(j).isolation < 0.10 and electrons.at(j).pT > 10.0)
           {
           Electron.SetPtEtaPhiE(electrons.at(j).pT, electrons.at(j).eta, electrons.at(j).phi, electrons.at(j).energy);
           double DRph_el = anaPhoton.PhLV.DeltaR(Electron);
           if(DRph_el<0.5) isGoodPhoton=false;
           }
         }
       for(unsigned int j=0; j<muons.size(); ++j)
         {
         TLorentzVector Muon;
         if(muons.at(j).isTight==1 and muons.at(j).isolation < 0.12){
         Muon.SetPtEtaPhiE(muons.at(j).pT, muons.at(j).eta, muons.at(j).phi, muons.at(j).energy);
         double DRph_mu = anaPhoton.PhLV.DeltaR(Muon);
         if(DRph_mu<0.5) isGoodPhoton=false;
         }
       }
       if(type=="Data")
       {
       TLorentzVector trigger1_p4;
       trigger1_p4.SetPxPyPzE(trigger1.Px, trigger1.Py, trigger1.Pz, trigger1.E);
       double DRph_tr = anaPhoton.PhLV.DeltaR(trigger1_p4);
       if(DRph_tr>0.3) isGoodPhoton=false;//matching to trigger photon object.
       }
     if(isGoodPhoton) 
       {
       anaPhotons.push_back(anaPhoton);
       }//filling the anaPhotons container
     }//close four vector if
   }//close photon loop
   // Now sorting this vector of structs
   std::sort(anaPhotons.begin(), anaPhotons.end(), sortPhotonVectorsInDescendingpT);   

   //photons in the sideband region
   std::vector<AnalysisPhotonInfo> SBPhotons;
   swissCross.clear();
   for(unsigned int k=0; k<photons.size(); ++k)
   {
     AnalysisPhotonInfo SBPhoton; //side band photon without sigma ietaieta

     swissCross.push_back((photons.at(k).phseedCrystalEnergy)/(photons.at(k).phe1x5 + photons.at(k).phe5x1 - photons.at(k).phseedCrystalEnergy));

     if(not(photons.at(0).chIsolation > 2.0 and photons.at(0).chIsolation < 6.0)) continue;

     if(fabs(photons.at(k).eta)<1.444 and photons.at(k).pT>30.0 and photons.at(k).phHoE < 0.05 and photons.at(0).phR9 > 0.9 and photons.at(k).phpixelVeto==0 and photons.at(0).nuIsolation < (0.4 + 0.04*photons.at(k).pT) and photons.at(0).phIsolation < (0.5 + 0.005*photons.at(k).pT) and swissCross.at(k)<0.90 and photons.at(k).phR9 < 1.0 and photons.at(k).phSigmaIetaIeta > 0.001 and photons.at(k).phSigmaIphiIphi>0.0001) //sideband region cut without sieie and swiss cross cut and additional cuts for beam halo rejection
     {
       SBPhoton.PhLV.SetPtEtaPhiE(photons.at(k).pT, photons.at(k).eta, photons.at(k).phi, photons.at(k).energy);
       SBPhoton.SigmaIetaIeta = ph_SigmaIetaIeta->at(k);
       bool isGoodPhoton=true;
       for(unsigned int j=0; j<electrons.size(); ++j)
         {
         TLorentzVector Electron;
         if(electrons.at(j).isTight==1 and electrons.at(j).isolation < 0.10 and electrons.at(j).pT > 10.0)
           {
           Electron.SetPtEtaPhiE(electrons.at(j).pT, electrons.at(j).eta, electrons.at(j).phi, electrons.at(j).energy);
           double DRph_el = SBPhoton.PhLV.DeltaR(Electron);
           if(DRph_el<0.5) isGoodPhoton=false;
           }
         }
       for(unsigned int j=0; j<muons.size(); ++j)
         {
         TLorentzVector Muon;
         if(muons.at(j).isTight==1 and muons.at(j).isolation < 0.12){
         Muon.SetPtEtaPhiE(muons.at(j).pT, muons.at(j).eta, muons.at(j).phi, muons.at(j).energy);
         double DRph_mu = SBPhoton.PhLV.DeltaR(Muon);
         if(DRph_mu<0.5) isGoodPhoton=false;
         }
       }
       if(type=="Data")
       {
       TLorentzVector trigger1_p4;
       trigger1_p4.SetPxPyPzE(trigger1.Px, trigger1.Py, trigger1.Pz, trigger1.E);
       double DRph_tr = SBPhoton.PhLV.DeltaR(trigger1_p4);
       if(DRph_tr>0.3) isGoodPhoton=false;//matching to trigger photon object.
       }
     if(isGoodPhoton)
       {
       SBPhotons.push_back(SBPhoton);
       }//filling the anaPhotons container
     }//close four vector if
   }//close photon loop
   // Now sorting this vector of structs
   std::sort(SBPhotons.begin(), SBPhotons.end(), sortPhotonVectorsInDescendingpT);
   //check lepton veto
   std::vector<TLorentzVector> Electron_vector;
   Electron_vector.clear();
   for (unsigned int j=0; j<electrons.size(); ++j) 
     {
     TLorentzVector Electron;
     
     if(electrons.at(j).pT > 10.0 and electrons.at(j).isTight==1 and electrons.at(j).isolation < 0.10)
       {
       Electron.SetPtEtaPhiE(electrons.at(j).pT, electrons.at(j).eta, electrons.at(j).phi, electrons.at(j).energy);
       bool isGoodElectron = true;
       if(type=="Data")
       {
         TLorentzVector trigger1_p4;
         trigger1_p4.SetPxPyPzE(trigger1.Px, trigger1.Py, trigger1.Pz, trigger1.E);
         double DRel_tr = Electron.DeltaR(trigger1_p4);
         if(DRel_tr<0.3) isGoodElectron=false;//anti-matching to trigger photon object.
       }
       if(isGoodElectron)
         {
         Electron_vector.push_back(Electron);
         }
       }//close four vector if
     }//close electron loop
   
   std::vector<TLorentzVector> Muon_vector;
   Muon_vector.clear();
   for(unsigned int j=0; j<muons.size(); ++j) 
     {
     TLorentzVector Muon;
     if(muons.at(j).pT > 5.0 and muons.at(j).isTight==1 and muons.at(j).isolation < 0.12)
       {
       Muon.SetPtEtaPhiE(muons.at(j).pT, muons.at(j).eta, muons.at(j).phi, muons.at(j).energy);
       Muon_vector.push_back(Muon);
       }//close four vector if
     }//close muon loop
   
   int nleptons =  Electron_vector.size() + Muon_vector.size();
   if (nleptons!=0) continue; 

   if(anaPhotons.size() > 0)
     {
     if(anaPhotons.at(0).PhLV.Pt() > 30.0 and anaPhotons.at(0).PhLV.Pt() < 40.0) h_SigmaIetaIeta_Num_Ph30To40->Fill(anaPhotons.at(0).SigmaIetaIeta);
     if(anaPhotons.at(0).PhLV.Pt() > 40.0 and anaPhotons.at(0).PhLV.Pt() < 50.0) h_SigmaIetaIeta_Num_Ph40To50->Fill(anaPhotons.at(0).SigmaIetaIeta);
     if(anaPhotons.at(0).PhLV.Pt() > 50.0 and anaPhotons.at(0).PhLV.Pt() < 60.0) h_SigmaIetaIeta_Num_Ph50To60->Fill(anaPhotons.at(0).SigmaIetaIeta);
     if(anaPhotons.at(0).PhLV.Pt() > 60.0 and anaPhotons.at(0).PhLV.Pt() < 70.0) h_SigmaIetaIeta_Num_Ph60To70->Fill(anaPhotons.at(0).SigmaIetaIeta);
     if(anaPhotons.at(0).PhLV.Pt() > 70.0 and anaPhotons.at(0).PhLV.Pt() < 80.0) h_SigmaIetaIeta_Num_Ph70To80->Fill(anaPhotons.at(0).SigmaIetaIeta);
     if(anaPhotons.at(0).PhLV.Pt() > 80.0 and anaPhotons.at(0).PhLV.Pt() < 100.0) h_SigmaIetaIeta_Num_Ph80To100->Fill(anaPhotons.at(0).SigmaIetaIeta);
     if(anaPhotons.at(0).PhLV.Pt() > 100.0 and anaPhotons.at(0).PhLV.Pt() < 130.0) h_SigmaIetaIeta_Num_Ph100To130->Fill(anaPhotons.at(0).SigmaIetaIeta);
   
     if(type=="MC"){
       double scaledSigmaIetaIeta = 0.891832*anaPhotons.at(0).SigmaIetaIeta + 0.0009133; //H->gg requirement 
       if(anaPhotons.at(0).PhLV.Pt() > 30.0 and anaPhotons.at(0).PhLV.Pt() < 40.0 and anaPhotons.at(0).matched==1) h_SigmaIetaIeta_Matched_Ph30To40->Fill(scaledSigmaIetaIeta);
       if(anaPhotons.at(0).PhLV.Pt() > 40.0 and anaPhotons.at(0).PhLV.Pt() < 50.0 and anaPhotons.at(0).matched==1) h_SigmaIetaIeta_Matched_Ph40To50->Fill(scaledSigmaIetaIeta);
       if(anaPhotons.at(0).PhLV.Pt() > 50.0 and anaPhotons.at(0).PhLV.Pt() < 60.0 and anaPhotons.at(0).matched==1) h_SigmaIetaIeta_Matched_Ph50To60->Fill(scaledSigmaIetaIeta);
       if(anaPhotons.at(0).PhLV.Pt() > 60.0 and anaPhotons.at(0).PhLV.Pt() < 70.0 and anaPhotons.at(0).matched==1) h_SigmaIetaIeta_Matched_Ph60To70->Fill(scaledSigmaIetaIeta);
       if(anaPhotons.at(0).PhLV.Pt() > 70.0 and anaPhotons.at(0).PhLV.Pt() < 80.0 and anaPhotons.at(0).matched==1) h_SigmaIetaIeta_Matched_Ph70To80->Fill(scaledSigmaIetaIeta);
       if(anaPhotons.at(0).PhLV.Pt() > 80.0 and anaPhotons.at(0).PhLV.Pt() < 100.0 and anaPhotons.at(0).matched==1) h_SigmaIetaIeta_Matched_Ph80To100->Fill(scaledSigmaIetaIeta);
       if(anaPhotons.at(0).PhLV.Pt() > 100.0 and anaPhotons.at(0).PhLV.Pt() < 130.0 and anaPhotons.at(0).matched==1) h_SigmaIetaIeta_Matched_Ph100To130->Fill(scaledSigmaIetaIeta);
      }//for MC matched photons
    }//anaPhotons size check
   
   if(SBPhotons.size() > 0)
     {
     if(SBPhotons.at(0).PhLV.Pt() > 30.0 and SBPhotons.at(0).PhLV.Pt() < 40.0) h_SigmaIetaIeta_SB_Ph30To40->Fill(SBPhotons.at(0).SigmaIetaIeta);
     if(SBPhotons.at(0).PhLV.Pt() > 40.0 and SBPhotons.at(0).PhLV.Pt() < 50.0) h_SigmaIetaIeta_SB_Ph40To50->Fill(SBPhotons.at(0).SigmaIetaIeta);
     if(SBPhotons.at(0).PhLV.Pt() > 50.0 and SBPhotons.at(0).PhLV.Pt() < 60.0) h_SigmaIetaIeta_SB_Ph50To60->Fill(SBPhotons.at(0).SigmaIetaIeta);
     if(SBPhotons.at(0).PhLV.Pt() > 60.0 and SBPhotons.at(0).PhLV.Pt() < 70.0) h_SigmaIetaIeta_SB_Ph60To70->Fill(SBPhotons.at(0).SigmaIetaIeta);
     if(SBPhotons.at(0).PhLV.Pt() > 70.0 and SBPhotons.at(0).PhLV.Pt() < 80.0) h_SigmaIetaIeta_SB_Ph70To80->Fill(SBPhotons.at(0).SigmaIetaIeta);
     if(SBPhotons.at(0).PhLV.Pt() > 80.0 and SBPhotons.at(0).PhLV.Pt() < 100.0) h_SigmaIetaIeta_SB_Ph80To100->Fill(SBPhotons.at(0).SigmaIetaIeta);
     if(SBPhotons.at(0).PhLV.Pt() > 100.0 and SBPhotons.at(0).PhLV.Pt() < 130.0) h_SigmaIetaIeta_SB_Ph100To130->Fill(SBPhotons.at(0).SigmaIetaIeta);
     }

  }//event loop closed

  //Cut flow table
  std::string histfilename=(outfile+".root").c_str();
  TFile *tFile=new TFile(histfilename.c_str(), "RECREATE");
  h_SigmaIetaIeta_Num_Ph30To40->Write();
  h_SigmaIetaIeta_SB_Ph30To40->Write();
  h_SigmaIetaIeta_Num_Ph40To50->Write();
  h_SigmaIetaIeta_SB_Ph40To50->Write();
  h_SigmaIetaIeta_Num_Ph50To60->Write();
  h_SigmaIetaIeta_SB_Ph50To60->Write();
  h_SigmaIetaIeta_Num_Ph60To70->Write();
  h_SigmaIetaIeta_SB_Ph60To70->Write();
  h_SigmaIetaIeta_Num_Ph70To80->Write();
  h_SigmaIetaIeta_SB_Ph70To80->Write();
  h_SigmaIetaIeta_Num_Ph80To100->Write();
  h_SigmaIetaIeta_SB_Ph80To100->Write();
  h_SigmaIetaIeta_Num_Ph100To130->Write();
  h_SigmaIetaIeta_SB_Ph100To130->Write();
  h_SigmaIetaIeta_Matched_Ph30To40->Write();
  h_SigmaIetaIeta_Matched_Ph40To50->Write();
  h_SigmaIetaIeta_Matched_Ph50To60->Write();
  h_SigmaIetaIeta_Matched_Ph60To70->Write();
  h_SigmaIetaIeta_Matched_Ph70To80->Write();
  h_SigmaIetaIeta_Matched_Ph80To100->Write();
  h_SigmaIetaIeta_Matched_Ph100To130->Write();
  tFile->Close(); 
  std::cout<<"Wrote output file "<<histfilename<<std::endl;

  return 0;

  
}
