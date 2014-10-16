#include "ReadLowPtSUSY_Tree_ComparisonDataMC_FixingNuances.h"

int ReadLowPtSUSY_Tree_ComparisonDataMC_FixingNuances(std::string infile, std::string outfile, std::string type, std::string signSelection, std::string MCSample){
  
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
  vector<int>     *ph_pixelVeto;
  vector<float>   *ph_e1x5;
  vector<float>   *ph_e1x3;
  vector<float>   *ph_e2x2;
  vector<float>   *ph_e2x5;
  vector<float>   *ph_e5x1;
  vector<float>   *ph_e5x5;
  vector<float>   *ph_e2x5Max;
  vector<float>   *ph_e2OverE5;
  vector<float>   *ph_seedCrystalEnergy;
  vector<float>   *ph_HoE;
  vector<int>     *ph_conversionVeto;
  vector<float>   *ph_SigmaIetaIeta;
  vector<float>   *ph_SigmaIetaIphi;
  vector<float>   *ph_SigmaIphiIphi;
  vector<float>   *ph_preShowerOverRaw;
  vector<float>   *ph_R9;
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
  ph_pixelVeto = 0;
  ph_HoE = 0;
  ph_conversionVeto = 0;
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
  tree->SetBranchAddress("ph_pixelVeto", &(ph_pixelVeto));
  tree->SetBranchAddress("ph_SigmaIetaIeta", &(ph_SigmaIetaIeta));
  tree->SetBranchAddress("ph_SigmaIetaIphi", &(ph_SigmaIetaIphi));
  tree->SetBranchAddress("ph_SigmaIphiIphi", &(ph_SigmaIphiIphi));
  tree->SetBranchAddress("ph_preShowerOverRaw", &(ph_preShowerOverRaw));
  tree->SetBranchAddress("ph_HoE", &(ph_HoE));
  tree->SetBranchAddress("ph_conversionVeto", &(ph_conversionVeto));
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
    tree->SetBranchAddress("trigObj3Px", &(trigObj3Px));
    tree->SetBranchAddress("trigObj3Py", &(trigObj3Py));
    tree->SetBranchAddress("trigObj3Pz", &(trigObj3Pz));
    tree->SetBranchAddress("trigObj3E", &(trigObj3E));
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

  HistCollection mumuHist;
  HistCollection mumuHistNT01;
  HistCollection mumuHistNT10;
  HistCollection mumuHistNT00;
  initializeHistCollection(mumuHist, "MuMu"); 
  initializeHistCollection(mumuHistNT01, "MuMu_NT01");
  initializeHistCollection(mumuHistNT10, "MuMu_NT10");
  initializeHistCollection(mumuHistNT00, "MuMu_NT00");
  HistCollection emuHist;
  HistCollection emuHistNT01;
  HistCollection emuHistNT10;
  HistCollection emuHistNT00;
  initializeHistCollection(emuHist, "ElMu");
  initializeHistCollection(emuHistNT01, "ElMu_NT01");
  initializeHistCollection(emuHistNT10, "ElMu_NT10");
  initializeHistCollection(emuHistNT00, "ElMu_NT00");
  HistCollection eeHist;
  HistCollection eeHistNT01;
  HistCollection eeHistNT10;
  HistCollection eeHistNT00;
  initializeHistCollection(eeHist, "ElEl"); 
  initializeHistCollection(eeHistNT01, "ElEl_NT01");
  initializeHistCollection(eeHistNT10, "ElEl_NT10");
  initializeHistCollection(eeHistNT00, "ElEl_NT00"); 
 
  TH2F *h_CaloMETTrue_CaloMET = new TH2F("h_CaloMETTrue_CaloMET", "Correlation between CaloMET_True and CaloMET; CaloMET_True [GeV]; CaloMet [GeV]", 600, 0, 600, 600, 0, 600); h_CaloMETTrue_CaloMET->Sumw2();
  TH2F *h_HLTMET_CaloMET = new TH2F("h_HLTMET_CaloMET", "Correlation between HLT MET and CaloMet Reco; HLT MET [GeV]; CaloMet Reco [GeV]", 600, 0, 600, 600, 0, 600); h_HLTMET_CaloMET->Sumw2();

  TH2F *h_Mmumu_MmumuGamma = new TH2F("h_Mmumu_MmumuGamma", "Scatter Plot of  M_{#mu#mu} versus M_{#mu#mu#gamma}; M_{#mu#mu} [GeV]; M_{#mu#mu#gamma} [GeV]", 4000, 0, 2000, 4000, 0, 2000);  h_Mmumu_MmumuGamma->Sumw2();

  double mumu_SS_NT01 = 0;
  double mumu_SS_NT10 = 0;
  double mumu_SS_NT00 = 0;
  double MET_emu_SS = 0;
  double MET_emu_OS = 0;
  double MET_mumu_SS = 0;
  double MET_mumu_OS = 0;
  double LowHT_mumu_OS = 0;
  double LowHT_mumu_SS = 0;
  double LowHT_emu_OS = 0;
  double LowHT_emu_SS = 0;
  double emu_OS = 0;
  double emu_SS = 0;
  double emu_SS_NT01 = 0;
  double emu_SS_NT10 = 0;
  double emu_SS_NT00 = 0;
  double mumu_OS = 0;
  double mumu_SS = 0;
  double ee_OS = 0;
  double ee_SS = 0;
  double ee_SS_NT01 = 0;
  double ee_SS_NT10 = 0;
  double ee_SS_NT00 = 0;

  TFile *trigger1D1=new TFile("TurnOn_SinglePhotonParked_Run2012D_22Jan2013_All.root");
  //TFile *trigger1D2=new TFile("TurnOn_METParked_Run2012D_All.root");//evaluated at the analysis photon ID point
  TFile *trigger1D2=new TFile("Test_METParked_Run2012D_All_OLD.root");

  TF1* fit_curve1=(TF1*)trigger1D1->Get("fit_caloMET");
  TF1* fit_curve2=(TF1*)trigger1D2->Get("fit_PHO");

  TFile* fileData = TFile::Open("PU_Data/PU_Data_All.root");
  TH1F *h_nVertices_Data = (TH1F*) fileData->Get("h_nVertices_Data");
  h_nVertices_Data->Scale(1.0/h_nVertices_Data->Integral()); 
  TFile* fileMC;
  TH1F *h_nWeights;
  if (type=="MC"){ 
    if(MCSample=="ZGToLLG") fileMC = TFile::Open("PU_MC/Output_LowPtSUSY_Tree_ZGToLLG_8TeV_CaloMET_All.root");
    if(MCSample=="TTG") fileMC = TFile::Open("PU_MC/Output_LowPtSUSY_Tree_TTGJets_8TeV_madgraph_CaloMET_All.root");
    if(MCSample=="ZZJetsTo4L") fileMC = TFile::Open("PU_MC/Output_LowPtSUSY_Tree_ZZJetsTo4L_TuneZ2star_CaloMET_All.root");
    if(MCSample=="Tbar_tW") fileMC = TFile::Open("PU_MC/Output_LowPtSUSY_Tree_Tbar_tW_CaloMET_All.root");
    if(MCSample=="T_tW")  fileMC = TFile::Open("PU_MC/Output_LowPtSUSY_Tree_T_tW_CaloMET_All.root");
    if(MCSample=="WGstarToLNu2Mu") fileMC = TFile::Open("PU_MC/Output_LowPtSUSY_Tree_WGstarToLNu2Mu_TuneZ2star_CaloMET_All.root");
    if(MCSample=="WmWmqq") fileMC = TFile::Open("PU_MC/Output_LowPtSUSY_Tree_WmWmqq_8TeV_madgraph_CaloMET_All.root");
    if(MCSample=="WpWpqq") fileMC = TFile::Open("PU_MC/Output_LowPtSUSY_Tree_WpWpqq_8TeV_madgraph_CaloMET_All.root");
    if(MCSample=="WZJetsTo3LNu") fileMC = TFile::Open("PU_MC/Output_LowPtSUSY_Tree_WZJetsTo3LNu_CaloMET_All.root");
    if(MCSample=="WWGJets") fileMC = TFile::Open("PU_MC/Output_LowPtSUSY_Tree_WWGJets_8TeV_madgraph_CaloMET_All.root");
    TH1F *h_nVertices_MC = (TH1F*) fileMC->Get("h_nVertices_MC");
    h_nVertices_MC->Scale(1.0/h_nVertices_MC->Integral());
    h_nWeights = (TH1F*) h_nVertices_Data->Clone("h_nWeights");
    h_nWeights->Divide(h_nVertices_MC);
  }
  int nEvents=tree->GetEntries();
  std::cout << "nEvents= " << nEvents << std::endl;
  for (int i=0; i<nEvents ; ++i)
    {
     tree->GetEvent(i);

     if(type=="Data")
       {
       if(fired_HLTPhoIdMet!=1) continue; 
       }
/*
     if(type=="MC")
       {
       if(Zmumu==false) continue;
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
       photons.push_back(photon);
       }

     // Now sorting this vector of structs
     std::sort (photons.begin(), photons.end(), sortPhotonsInDescendingpT);
       
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

   vector<TLorentzVector> Photon_vector;
   Photon_vector.clear();
   std::vector<float> swissCross;
   swissCross.clear();
   for(unsigned int k=0; k<photons.size(); ++k)
   {
     TLorentzVector Photon;
     swissCross.push_back((photons.at(k).phseedCrystalEnergy)/(photons.at(k).phe1x5 + photons.at(k).phe5x1 - photons.at(k).phseedCrystalEnergy));
     
     if(fabs(photons.at(k).eta)<1.444 and photons.at(k).pT>30.0 and photons.at(k).isTight==1 and photons.at(k).phIsoTight==1 and photons.at(k).phpixelVeto==0 and swissCross.at(k)<0.90 and photons.at(k).phR9 < 1.0 and photons.at(k).phSigmaIetaIeta > 0.001 and photons.at(k).phSigmaIphiIphi>0.0)//0.0001
    {
       Photon.SetPtEtaPhiE(photons.at(k).pT, photons.at(k).eta, photons.at(k).phi, photons.at(k).energy);
       bool isGoodPhoton=true;
       for(unsigned int j=0; j<electrons.size(); ++j)
         {
         TLorentzVector Electron;
         if(electrons.at(j).isTight==1 and electrons.at(j).isolation < 0.10 and electrons.at(j).pT > 10.0)
           {
           Electron.SetPtEtaPhiE(electrons.at(j).pT, electrons.at(j).eta, electrons.at(j).phi, electrons.at(j).energy);
           double DRph_el = Photon.DeltaR(Electron);
           if(DRph_el<0.5) isGoodPhoton=false;
           }
         }

       for(unsigned int j=0; j<electrons.size(); ++j)
         {
         TLorentzVector LNTElectron; 
         if(electrons.at(j).pT > 10.0 and electrons.at(j).isLoose==1 and electrons.at(j).isolation < 0.60 and not(electrons.at(j).isTight==1 and electrons.at(j).isolation < 0.10)){
           LNTElectron.SetPtEtaPhiE(electrons.at(j).pT, electrons.at(j).eta, electrons.at(j).phi, electrons.at(j).energy);
           double DRph_lntel = Photon.DeltaR(LNTElectron);
           if(DRph_lntel<0.5) isGoodPhoton=false;
          }
        }
       for(unsigned int j=0; j<muons.size(); ++j)
         {
         TLorentzVector Muon;
         if(muons.at(j).pT > 5.0 and muons.at(j).isTight==1 and muons.at(j).isolation < 0.12){
         Muon.SetPtEtaPhiE(muons.at(j).pT, muons.at(j).eta, muons.at(j).phi, muons.at(j).energy);
         double DRph_mu = Photon.DeltaR(Muon);
         if(DRph_mu<0.5) isGoodPhoton=false;
        }
      }

      for(unsigned int j=0; j<muons.size(); ++j)
         {
         TLorentzVector LNTMuon;
         if(muons.at(j).pT > 5.0 and muons.at(j).isLoose==1 and muons.at(j).isolation < 0.40 and not(muons.at(j).isTight==1 and muons.at(j).isolation < 0.12)){
           LNTMuon.SetPtEtaPhiE(muons.at(j).pT, muons.at(j).eta, muons.at(j).phi, muons.at(j).energy);
           double DRph_lntmu = Photon.DeltaR(LNTMuon);
           if(DRph_lntmu<0.5) isGoodPhoton=false;
        }
      }

      if(type=="Data")
       {
       TLorentzVector trigger1_p4;
       trigger1_p4.SetPxPyPzE(trigger1.Px, trigger1.Py, trigger1.Pz, trigger1.E);
       double DRph_tr = Photon.DeltaR(trigger1_p4);
       if(DRph_tr>0.3) isGoodPhoton=false;//matching to trigger photon object.
       }
      if(isGoodPhoton) Photon_vector.push_back(Photon);
    }//close four vector if
  }//close photon loop

// Now sorting this vector of structs
   std::sort (Photon_vector.begin(), Photon_vector.end(), sortPhotonVectorsInDescendingpT);
     
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
        if(electrons.at(j).pT > 10.0 and electrons.at(j).isTight==1 and electrons.at(j).isolation < 0.10)
          {
          Electron.SetPtEtaPhiE(electrons.at(j).pT, electrons.at(j).eta, electrons.at(j).phi, electrons.at(j).energy);
          double DRjet_el = Jet.DeltaR(Electron);
          if(DRjet_el<0.5) isGoodJet=false;
          }
        }
   for(unsigned int m=0; m<muons.size(); ++m)
     {
       TLorentzVector Muon;
       if(muons.at(m).pT > 5.0  and muons.at(m).isTight==1 and muons.at(m).isolation < 0.12){
       Muon.SetPtEtaPhiE(muons.at(m).pT, muons.at(m).eta, muons.at(m).phi, muons.at(m).energy);
       double DRjet_mu = Jet.DeltaR(Muon);
       if(DRjet_mu<0.5) isGoodJet=false;
     }
   }
   for(unsigned int l=0; l<photons.size(); ++l)
     {
       TLorentzVector Photon;
       swissCross.push_back((photons.at(l).phseedCrystalEnergy)/(photons.at(l).phe1x5 + photons.at(l).phe5x1 - photons.at(l).phseedCrystalEnergy));

       if(photons.at(l).pT>30.0 and photons.at(l).isTight==1 and photons.at(l).phIsoTight==1 and photons.at(l).phpixelVeto==0 and swissCross.at(l)<0.90 and photons.at(l).phR9 < 1.0 and photons.at(l).phSigmaIetaIeta > 0.001 and photons.at(l).phSigmaIphiIphi>0.0){//>0.0001
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

  std::vector<AnalysisLeptonInfo> Electrons;
  Electrons.clear();
  for (unsigned int j=0; j<electrons.size(); ++j)
     {
     AnalysisLeptonInfo Electron;
     if(electrons.at(j).pT > 10.0 and electrons.at(j).isTight==1 and electrons.at(j).isolation < 0.10)
       {
       Electron.LepLV.SetPtEtaPhiE(electrons.at(j).pT, electrons.at(j).eta, electrons.at(j).phi, electrons.at(j).energy);
       Electron.Charge=electrons.at(j).charge;
       Electron.Isolation=electrons.at(j).isolation;
       bool isGoodElectron = true;
       if(type=="Data")
       {
         TLorentzVector trigger1_p4;
         trigger1_p4.SetPxPyPzE(trigger1.Px, trigger1.Py, trigger1.Pz, trigger1.E);
         double DRel_tr = Electron.LepLV.DeltaR(trigger1_p4);
         if(DRel_tr<0.3) isGoodElectron=false;//anti-matching to trigger photon object.
       }
       if(isGoodElectron) Electrons.push_back(Electron);
       }//close four vector if
     }//close electron loop

  std::sort(Electrons.begin(), Electrons.end(), sortVectorsInDescendingpT);

  std::vector<AnalysisLeptonInfo> LNTElectrons;
  LNTElectrons.clear();
  for (unsigned int j=0; j<electrons.size(); ++j)
     {
     AnalysisLeptonInfo LNTElectron;
     if(electrons.at(j).pT > 10.0 and electrons.at(j).isLoose==1 and electrons.at(j).isolation < 0.60 and not(electrons.at(j).isTight==1 and electrons.at(j).isolation < 0.10))
       {
       LNTElectron.LepLV.SetPtEtaPhiE(electrons.at(j).pT, electrons.at(j).eta, electrons.at(j).phi, electrons.at(j).energy);
       LNTElectron.Charge=electrons.at(j).charge;
       LNTElectron.Isolation=electrons.at(j).isolation;
       bool isLNTElectron = true;
       if(type=="Data")
       {
         TLorentzVector trigger1_p4;
         trigger1_p4.SetPxPyPzE(trigger1.Px, trigger1.Py, trigger1.Pz, trigger1.E);
         double DRel_tr = LNTElectron.LepLV.DeltaR(trigger1_p4);
         if(DRel_tr<0.3) isLNTElectron=false;//anti-matching to trigger photon object.
       }
       if(isLNTElectron) LNTElectrons.push_back(LNTElectron);
       }//close four vector if
     }//close electron loop

   std::sort(LNTElectrons.begin(), LNTElectrons.end(), sortVectorsInDescendingpT);

   std::vector<AnalysisLeptonInfo> Muons;
   Muons.clear();
   for(unsigned int j=0; j<muons.size(); ++j)
     {
     AnalysisLeptonInfo Muon;
     if(muons.at(j).pT > 5.0 and muons.at(j).isTight==1 and muons.at(j).isolation < 0.12)
       {
       Muon.LepLV.SetPtEtaPhiE(muons.at(j).pT, muons.at(j).eta, muons.at(j).phi, muons.at(j).energy);
       Muon.Charge=muons.at(j).charge;
       Muon.Isolation=muons.at(j).isolation;
       Muons.push_back(Muon);
       }//close four vector if
     }//close muon loop

   std::sort(Muons.begin(), Muons.end(), sortVectorsInDescendingpT);
   std::vector<AnalysisLeptonInfo> LNTMuons;
   LNTMuons.clear();
   for(unsigned int j=0; j<muons.size(); ++j)
     {
     AnalysisLeptonInfo LNTMuon;
     if(muons.at(j).pT > 5.0 and muons.at(j).isLoose==1 and muons.at(j).isolation < 0.40 and not (muons.at(j).isTight==1 and muons.at(j).isolation < 0.12))
       {
       LNTMuon.LepLV.SetPtEtaPhiE(muons.at(j).pT, muons.at(j).eta, muons.at(j).phi, muons.at(j).energy);
       LNTMuon.Charge=muons.at(j).charge;
       LNTMuon.Isolation=muons.at(j).isolation;
       LNTMuons.push_back(LNTMuon);
       }//close four vector if
     }//close muon loop

  std::sort(LNTMuons.begin(), LNTMuons.end(), sortVectorsInDescendingpT);

  double caloMET_True = caloMET;
  if(Muons.size() > 1){
    double MEx = MET*cos(MET_Phi) + Muons.at(0).LepLV.Px() + Muons.at(1).LepLV.Px();
    double MEy = MET*sin(MET_Phi) + Muons.at(0).LepLV.Py() + Muons.at(1).LepLV.Py();
    double caloMet = sqrt(MEx*MEx + MEy*MEy);

    h_CaloMETTrue_CaloMET->Fill(caloMET_True, caloMet);
     if(type=="Data"){
       h_HLTMET_CaloMET->Fill(sqrt(trigger3.Px*trigger3.Px + trigger3.Py*trigger3.Py), caloMET);
    }
  }

  //Preliminary event cuts guided by the trigger:
  if (not(Photon_vector.size() > 0 and Photon_vector.at(0).Pt() > 35.0 and caloMET > 30.0)) continue;

  //application of trigger and pileup weights
  double eventWeight_Trigger = 0.0;
  double eventWeight = 0.0;
    if(type=="MC")
    {
      double triggerWeight1 = 0.0;
      if(caloMET>25.0 and caloMET < 200.0)triggerWeight1 = fit_curve1->Eval(caloMET);
      else if(caloMET > 200.0) triggerWeight1 = fit_curve1->Eval(200.0);
      double triggerWeight2 = 0.0;
      if(Photon_vector.at(0).Pt()>30.0 and Photon_vector.at(0).Pt()<100.0) triggerWeight2 = fit_curve2->Eval(Photon_vector.at(0).Pt());
      else if(Photon_vector.at(0).Pt()>100.0) triggerWeight2 = fit_curve2->Eval(100.0);
      eventWeight_Trigger=triggerWeight1*triggerWeight2;
      double PU_weights = h_nWeights->GetBinContent(nVertices);
      //cout << "PU_weights = " << PU_weights << endl;
      eventWeight=triggerWeight1*triggerWeight2*PU_weights;    
    }
    else if(type=="Data")
    {
      eventWeight = 1.0;
    }

  TVector2 mu1_transverse;
  TVector2 mu2_transverse;
  TVector2 met_transverse;
  TVector2 ph_transverse;
  TVector2 el1_transverse;
  TVector2 el2_transverse;
  ph_transverse.SetMagPhi(Photon_vector.at(0).Pt(), Photon_vector.at(0).Phi());
  met_transverse.SetMagPhi(MET, MET_Phi); 
  /*
  if(Muons.size()==2 and eventWeight > 0.0){
    if(Muons.at(0).Charge*Muons.at(1).Charge==-1) cout << "OS Dimuon event found" << endl;
    if(Muons.at(0).Charge*Muons.at(1).Charge==+1) cout << "SS Dimuon event found" << endl;
  }
  if(Electrons.size()==2 and eventWeight > 0.0){
    if(Electrons.at(0).Charge*Electrons.at(1).Charge==-1) cout << "OS Dielectron event found" << endl;
    if(Electrons.at(0).Charge*Electrons.at(1).Charge==+1) cout << "SS Dielectron event found" << endl;
  }
  if(Electrons.size()==1 and Muons.size()==1 and eventWeight > 0.0){
    if(Muons.at(0).Charge*Electrons.at(0).Charge==-1) cout << "OS Dilepton event found" << endl; 
    if(Muons.at(0).Charge*Electrons.at(0).Charge==+1) cout << "SS Dilepton event found" << endl;
  }
  */
  if(eventWeight > 0.0 and Muons.size()==2){ 
      if(type=="MC") eventWeight *= muonSF(Muons.at(0).LepLV.Pt(), Muons.at(0).LepLV.Eta())*muonSF(Muons.at(1).LepLV.Pt(), Muons.at(1).LepLV.Eta())*photonSF(Photon_vector.at(0).Pt(), Photon_vector.at(0).Eta()); 
      mu1_transverse.SetMagPhi(Muons.at(0).LepLV.Pt(), Muons.at(0).LepLV.Phi());
      mu2_transverse.SetMagPhi(Muons.at(1).LepLV.Pt(), Muons.at(1).LepLV.Phi()); 
      if(signSelection=="OS" and (Muons.at(0).Charge*Muons.at(1).Charge)==-1){
        mumu_OS+=eventWeight;
        h_Mmumu_MmumuGamma->Fill((Muons.at(0).LepLV+Muons.at(1).LepLV).M(), (Muons.at(0).LepLV+Muons.at(1).LepLV+Photon_vector.at(0)).M(), eventWeight);
        mumuHist.h_mu_pt_leading->Fill(Muons.at(0).LepLV.Pt(), eventWeight);
        mumuHist.h_mu_phi_leading->Fill(Muons.at(0).LepLV.Phi(), eventWeight);
        mumuHist.h_mu_eta_leading->Fill(Muons.at(0).LepLV.Eta(), eventWeight);
        mumuHist.h_mu_energy_leading->Fill(Muons.at(0).LepLV.E(), eventWeight);
        mumuHist.h_mu_phi_trailing->Fill(Muons.at(1).LepLV.Phi(), eventWeight);
        mumuHist.h_mu_eta_trailing->Fill(Muons.at(1).LepLV.Eta(), eventWeight);
        mumuHist.h_mu_pt_trailing->Fill(Muons.at(1).LepLV.Pt(), eventWeight);
        mumuHist.h_mu_energy_trailing->Fill(Muons.at(1).LepLV.E(), eventWeight);
        mumuHist.h_InvariantMass->Fill((Muons.at(0).LepLV+Muons.at(1).LepLV).M(), eventWeight);
        mumuHist.h_InvariantMass_Ph->Fill((Muons.at(0).LepLV+Muons.at(1).LepLV+Photon_vector.at(0)).M(), eventWeight);
        mumuHist.h_DeltaPhi_met_mu1->Fill(met_transverse.DeltaPhi(mu1_transverse), eventWeight);
        mumuHist.h_DeltaPhi_met_mu2->Fill(met_transverse.DeltaPhi(mu2_transverse), eventWeight);
        mumuHist.h_DeltaPhi_ph_mu1->Fill(ph_transverse.DeltaPhi(mu1_transverse), eventWeight);
        mumuHist.h_DeltaPhi_ph_mu2->Fill(ph_transverse.DeltaPhi(mu2_transverse), eventWeight);
        mumuHist.h_Isolation_mu1->Fill(Muons.at(0).Isolation, eventWeight);
        mumuHist.h_Isolation_mu2->Fill(Muons.at(1).Isolation, eventWeight);
        mumuHist.h_nVertices->Fill(nVertices, eventWeight);
        mumuHist.h_MET->Fill(MET, eventWeight);
        mumuHist.h_caloMET->Fill(caloMET, eventWeight);
        mumuHist.h_nJets->Fill(Jet_vector.size(), eventWeight);
        mumuHist.h_HT->Fill(HT, eventWeight);
        if(Jet_vector.size()>0) mumuHist.h_jet_pt_leading->Fill(Jet_vector.at(0).Pt(), eventWeight);
        if(Jet_vector.size()>1) mumuHist.h_jet_pt_trailing->Fill(Jet_vector.at(1).Pt(), eventWeight);
        if(Jet_vector.size()>2) mumuHist.h_jet_pt_3rd->Fill(Jet_vector.at(2).Pt(), eventWeight);
        if(Jet_vector.size()>3) mumuHist.h_jet_pt_4th->Fill(Jet_vector.at(3).Pt(), eventWeight);
        if(Jet_vector.size()>4) mumuHist.h_jet_pt_5th->Fill(Jet_vector.at(4).Pt(), eventWeight);
        if(Jet_vector.size()>5) mumuHist.h_jet_pt_6th->Fill(Jet_vector.at(5).Pt(), eventWeight);
        if(Jet_vector.size()>0) mumuHist.h_jet_eta_leading->Fill(Jet_vector.at(0).Eta(), eventWeight);
        if(Jet_vector.size()>1) mumuHist.h_jet_eta_trailing->Fill(Jet_vector.at(1).Eta(), eventWeight);
        if(Jet_vector.size()>2) mumuHist.h_jet_eta_3rd->Fill(Jet_vector.at(2).Eta(), eventWeight);
        if(Jet_vector.size()>3) mumuHist.h_jet_eta_4th->Fill(Jet_vector.at(3).Eta(), eventWeight);
        if(Jet_vector.size()>4) mumuHist.h_jet_eta_5th->Fill(Jet_vector.at(4).Eta(), eventWeight);
        if(Jet_vector.size()>5) mumuHist.h_jet_eta_6th->Fill(Jet_vector.at(5).Eta(), eventWeight);
        if(Jet_vector.size()>0) mumuHist.h_jet_phi_leading->Fill(Jet_vector.at(0).Phi(), eventWeight);
        if(Jet_vector.size()>1) mumuHist.h_jet_phi_trailing->Fill(Jet_vector.at(1).Phi(), eventWeight);
        if(Jet_vector.size()>2) mumuHist.h_jet_phi_3rd->Fill(Jet_vector.at(2).Phi(), eventWeight);
        if(Jet_vector.size()>3) mumuHist.h_jet_phi_4th->Fill(Jet_vector.at(3).Phi(), eventWeight);
        if(Jet_vector.size()>4) mumuHist.h_jet_phi_5th->Fill(Jet_vector.at(4).Phi(), eventWeight);
        if(Jet_vector.size()>5) mumuHist.h_jet_phi_6th->Fill(Jet_vector.at(5).Phi(), eventWeight);
        if(Jet_vector.size()>0) mumuHist.h_jet_energy_leading->Fill(Jet_vector.at(0).E(), eventWeight);
        if(Jet_vector.size()>1) mumuHist.h_jet_energy_trailing->Fill(Jet_vector.at(1).E(), eventWeight);
        if(Jet_vector.size()>2) mumuHist.h_jet_energy_3rd->Fill(Jet_vector.at(2).E(), eventWeight);
        if(Jet_vector.size()>3) mumuHist.h_jet_energy_4th->Fill(Jet_vector.at(3).E(), eventWeight);
        if(Jet_vector.size()>4) mumuHist.h_jet_energy_5th->Fill(Jet_vector.at(4).E(), eventWeight);
        if(Jet_vector.size()>5) mumuHist.h_jet_energy_6th->Fill(Jet_vector.at(5).E(), eventWeight);
        mumuHist.h_photon_pt->Fill(Photon_vector.at(0).Pt(), eventWeight);
        mumuHist.h_photon_eta->Fill(Photon_vector.at(0).Eta(), eventWeight);
        mumuHist.h_photon_phi->Fill(Photon_vector.at(0).Phi(), eventWeight);
        mumuHist.h_photon_energy->Fill(Photon_vector.at(0).E(), eventWeight);
        if(MET > 25.0) MET_mumu_OS+=eventWeight;
        if(HT<300) LowHT_mumu_OS+=eventWeight;
      }//OS mode
     else if(signSelection=="SS" and (Muons.at(0).Charge*Muons.at(1).Charge)==+1){
        mumu_SS+=eventWeight;
        mumuHist.h_mu_pt_leading->Fill(Muons.at(0).LepLV.Pt(), eventWeight);
        mumuHist.h_mu_phi_leading->Fill(Muons.at(0).LepLV.Phi(), eventWeight);
        mumuHist.h_mu_eta_leading->Fill(Muons.at(0).LepLV.Eta(), eventWeight);
        mumuHist.h_mu_energy_leading->Fill(Muons.at(0).LepLV.E(), eventWeight);
        mumuHist.h_mu_phi_trailing->Fill(Muons.at(1).LepLV.Phi(), eventWeight);
        mumuHist.h_mu_eta_trailing->Fill(Muons.at(1).LepLV.Eta(), eventWeight);
        mumuHist.h_mu_pt_trailing->Fill(Muons.at(1).LepLV.Pt(), eventWeight);
        mumuHist.h_mu_energy_trailing->Fill(Muons.at(1).LepLV.E(), eventWeight);
        mumuHist.h_InvariantMass->Fill((Muons.at(0).LepLV+Muons.at(1).LepLV).M(), eventWeight);
        mumuHist.h_InvariantMass_Ph->Fill((Muons.at(0).LepLV+Muons.at(1).LepLV+Photon_vector.at(0)).M(), eventWeight);
        mumuHist.h_DeltaPhi_met_mu1->Fill(met_transverse.DeltaPhi(mu1_transverse), eventWeight);
        mumuHist.h_DeltaPhi_met_mu2->Fill(met_transverse.DeltaPhi(mu2_transverse), eventWeight);
        mumuHist.h_DeltaPhi_ph_mu1->Fill(ph_transverse.DeltaPhi(mu1_transverse), eventWeight);
        mumuHist.h_DeltaPhi_ph_mu2->Fill(ph_transverse.DeltaPhi(mu2_transverse), eventWeight);
        mumuHist.h_Isolation_mu1->Fill(Muons.at(0).Isolation, eventWeight);
        mumuHist.h_Isolation_mu2->Fill(Muons.at(1).Isolation, eventWeight);
        mumuHist.h_nVertices->Fill(nVertices, eventWeight);
        mumuHist.h_MET->Fill(MET, eventWeight);
        mumuHist.h_caloMET->Fill(caloMET, eventWeight);
        mumuHist.h_nJets->Fill(Jet_vector.size(), eventWeight);
        mumuHist.h_HT->Fill(HT, eventWeight);
        if(Jet_vector.size()>0) mumuHist.h_jet_pt_leading->Fill(Jet_vector.at(0).Pt(), eventWeight);
        if(Jet_vector.size()>1) mumuHist.h_jet_pt_trailing->Fill(Jet_vector.at(1).Pt(), eventWeight);
        if(Jet_vector.size()>2) mumuHist.h_jet_pt_3rd->Fill(Jet_vector.at(2).Pt(), eventWeight);
        if(Jet_vector.size()>3) mumuHist.h_jet_pt_4th->Fill(Jet_vector.at(3).Pt(), eventWeight);
        if(Jet_vector.size()>4) mumuHist.h_jet_pt_5th->Fill(Jet_vector.at(4).Pt(), eventWeight);
        if(Jet_vector.size()>5) mumuHist.h_jet_pt_6th->Fill(Jet_vector.at(5).Pt(), eventWeight);
        if(Jet_vector.size()>0) mumuHist.h_jet_eta_leading->Fill(Jet_vector.at(0).Eta(), eventWeight);
        if(Jet_vector.size()>1) mumuHist.h_jet_eta_trailing->Fill(Jet_vector.at(1).Eta(), eventWeight);
        if(Jet_vector.size()>2) mumuHist.h_jet_eta_3rd->Fill(Jet_vector.at(2).Eta(), eventWeight);
        if(Jet_vector.size()>3) mumuHist.h_jet_eta_4th->Fill(Jet_vector.at(3).Eta(), eventWeight);
        if(Jet_vector.size()>4) mumuHist.h_jet_eta_5th->Fill(Jet_vector.at(4).Eta(), eventWeight);
        if(Jet_vector.size()>5) mumuHist.h_jet_eta_6th->Fill(Jet_vector.at(5).Eta(), eventWeight);
        if(Jet_vector.size()>0) mumuHist.h_jet_phi_leading->Fill(Jet_vector.at(0).Phi(), eventWeight);
        if(Jet_vector.size()>1) mumuHist.h_jet_phi_trailing->Fill(Jet_vector.at(1).Phi(), eventWeight);
        if(Jet_vector.size()>2) mumuHist.h_jet_phi_3rd->Fill(Jet_vector.at(2).Phi(), eventWeight);
        if(Jet_vector.size()>3) mumuHist.h_jet_phi_4th->Fill(Jet_vector.at(3).Phi(), eventWeight);
        if(Jet_vector.size()>4) mumuHist.h_jet_phi_5th->Fill(Jet_vector.at(4).Phi(), eventWeight);
        if(Jet_vector.size()>5) mumuHist.h_jet_phi_6th->Fill(Jet_vector.at(5).Phi(), eventWeight);
        if(Jet_vector.size()>0) mumuHist.h_jet_energy_leading->Fill(Jet_vector.at(0).E(), eventWeight);
        if(Jet_vector.size()>1) mumuHist.h_jet_energy_trailing->Fill(Jet_vector.at(1).E(), eventWeight);
        if(Jet_vector.size()>2) mumuHist.h_jet_energy_3rd->Fill(Jet_vector.at(2).E(), eventWeight);
        if(Jet_vector.size()>3) mumuHist.h_jet_energy_4th->Fill(Jet_vector.at(3).E(), eventWeight);
        if(Jet_vector.size()>4) mumuHist.h_jet_energy_5th->Fill(Jet_vector.at(4).E(), eventWeight);
        if(Jet_vector.size()>5) mumuHist.h_jet_energy_6th->Fill(Jet_vector.at(5).E(), eventWeight);
        mumuHist.h_photon_pt->Fill(Photon_vector.at(0).Pt(), eventWeight);
        mumuHist.h_photon_eta->Fill(Photon_vector.at(0).Eta(), eventWeight);
        mumuHist.h_photon_phi->Fill(Photon_vector.at(0).Phi(), eventWeight);
        mumuHist.h_photon_energy->Fill(Photon_vector.at(0).E(), eventWeight);
        if(MET > 25.0) MET_mumu_SS+=eventWeight;
        if(HT<300) LowHT_mumu_SS+=eventWeight;
       }//SS mode
   }//first+second tight muon
   else if(Muons.size() > 0.0 and LNTMuons.size() > 0.0 and eventWeight > 0.0 ) {
      if(type=="Data" and signSelection=="SS" and (Muons.at(0).Charge)*(LNTMuons.at(0).Charge)==1){
        mumu_SS_NT01+=eventWeight;
        mu1_transverse.SetMagPhi(Muons.at(0).LepLV.Pt(), Muons.at(0).LepLV.Phi());
        mu2_transverse.SetMagPhi(LNTMuons.at(0).LepLV.Pt(), LNTMuons.at(0).LepLV.Phi());  
        mumuHistNT01.h_mu_pt_leading->Fill(Muons.at(0).LepLV.Pt(), eventWeight);
        mumuHistNT01.h_mu_phi_leading->Fill(Muons.at(0).LepLV.Phi(), eventWeight);
        mumuHistNT01.h_mu_eta_leading->Fill(Muons.at(0).LepLV.Eta(), eventWeight);
        mumuHistNT01.h_mu_energy_leading->Fill(Muons.at(0).LepLV.E(), eventWeight);
        mumuHistNT01.h_mu_pt_trailing->Fill(LNTMuons.at(0).LepLV.Pt(), eventWeight);
        mumuHistNT01.h_mu_phi_trailing->Fill(LNTMuons.at(0).LepLV.Phi(), eventWeight);
        mumuHistNT01.h_mu_eta_trailing->Fill(LNTMuons.at(0).LepLV.Eta(), eventWeight);
        mumuHistNT01.h_mu_energy_trailing->Fill(LNTMuons.at(0).LepLV.E(), eventWeight);
        mumuHistNT01.h_InvariantMass->Fill((Muons.at(0).LepLV+LNTMuons.at(0).LepLV).M(), eventWeight);
        mumuHistNT01.h_InvariantMass_Ph->Fill((Muons.at(0).LepLV+LNTMuons.at(0).LepLV+Photon_vector.at(0)).M(), eventWeight);
        mumuHistNT01.h_DeltaPhi_met_mu1->Fill(met_transverse.DeltaPhi(mu1_transverse), eventWeight);
        mumuHistNT01.h_DeltaPhi_met_mu2->Fill(met_transverse.DeltaPhi(mu2_transverse), eventWeight);
        mumuHistNT01.h_DeltaPhi_ph_mu1->Fill(ph_transverse.DeltaPhi(mu1_transverse), eventWeight);
        mumuHistNT01.h_DeltaPhi_ph_mu2->Fill(ph_transverse.DeltaPhi(mu2_transverse), eventWeight);
        mumuHistNT01.h_Isolation_mu1->Fill(Muons.at(0).Isolation, eventWeight);
        mumuHistNT01.h_Isolation_mu2->Fill(LNTMuons.at(0).Isolation, eventWeight);
        mumuHistNT01.h_nVertices->Fill(nVertices, eventWeight);
        mumuHistNT01.h_MET->Fill(MET, eventWeight);
        mumuHistNT01.h_caloMET->Fill(caloMET, eventWeight);
        mumuHistNT01.h_nJets->Fill(Jet_vector.size(), eventWeight);
        mumuHistNT01.h_HT->Fill(HT, eventWeight); 
        if(Jet_vector.size()>0) mumuHistNT01.h_jet_pt_leading->Fill(Jet_vector.at(0).Pt(), eventWeight);
        if(Jet_vector.size()>1) mumuHistNT01.h_jet_pt_trailing->Fill(Jet_vector.at(1).Pt(), eventWeight);
        if(Jet_vector.size()>2) mumuHistNT01.h_jet_pt_3rd->Fill(Jet_vector.at(2).Pt(), eventWeight);
        if(Jet_vector.size()>3) mumuHistNT01.h_jet_pt_4th->Fill(Jet_vector.at(3).Pt(), eventWeight);
        if(Jet_vector.size()>4) mumuHistNT01.h_jet_pt_5th->Fill(Jet_vector.at(4).Pt(), eventWeight);
        if(Jet_vector.size()>5) mumuHistNT01.h_jet_pt_6th->Fill(Jet_vector.at(5).Pt(), eventWeight);
        if(Jet_vector.size()>0) mumuHistNT01.h_jet_eta_leading->Fill(Jet_vector.at(0).Eta(), eventWeight);
        if(Jet_vector.size()>1) mumuHistNT01.h_jet_eta_trailing->Fill(Jet_vector.at(1).Eta(), eventWeight);
        if(Jet_vector.size()>2) mumuHistNT01.h_jet_eta_3rd->Fill(Jet_vector.at(2).Eta(), eventWeight);
        if(Jet_vector.size()>3) mumuHistNT01.h_jet_eta_4th->Fill(Jet_vector.at(3).Eta(), eventWeight);
        if(Jet_vector.size()>4) mumuHistNT01.h_jet_eta_5th->Fill(Jet_vector.at(4).Eta(), eventWeight);
        if(Jet_vector.size()>5) mumuHistNT01.h_jet_eta_6th->Fill(Jet_vector.at(5).Eta(), eventWeight);
        if(Jet_vector.size()>0) mumuHistNT01.h_jet_phi_leading->Fill(Jet_vector.at(0).Phi(), eventWeight);
        if(Jet_vector.size()>1) mumuHistNT01.h_jet_phi_trailing->Fill(Jet_vector.at(1).Phi(), eventWeight);
        if(Jet_vector.size()>2) mumuHistNT01.h_jet_phi_3rd->Fill(Jet_vector.at(2).Phi(), eventWeight);
        if(Jet_vector.size()>3) mumuHistNT01.h_jet_phi_4th->Fill(Jet_vector.at(3).Phi(), eventWeight);
        if(Jet_vector.size()>4) mumuHistNT01.h_jet_phi_5th->Fill(Jet_vector.at(4).Phi(), eventWeight);
        if(Jet_vector.size()>5) mumuHistNT01.h_jet_phi_6th->Fill(Jet_vector.at(5).Phi(), eventWeight);
        if(Jet_vector.size()>0) mumuHistNT01.h_jet_energy_leading->Fill(Jet_vector.at(0).E(), eventWeight);
        if(Jet_vector.size()>1) mumuHistNT01.h_jet_energy_trailing->Fill(Jet_vector.at(1).E(), eventWeight);
        if(Jet_vector.size()>2) mumuHistNT01.h_jet_energy_3rd->Fill(Jet_vector.at(2).E(), eventWeight);
        if(Jet_vector.size()>3) mumuHistNT01.h_jet_energy_4th->Fill(Jet_vector.at(3).E(), eventWeight);
        if(Jet_vector.size()>4) mumuHistNT01.h_jet_energy_5th->Fill(Jet_vector.at(4).E(), eventWeight);
        if(Jet_vector.size()>5) mumuHistNT01.h_jet_energy_6th->Fill(Jet_vector.at(5).E(), eventWeight);
        mumuHistNT01.h_photon_pt->Fill(Photon_vector.at(0).Pt(), eventWeight);
        mumuHistNT01.h_photon_eta->Fill(Photon_vector.at(0).Eta(), eventWeight);
        mumuHistNT01.h_photon_phi->Fill(Photon_vector.at(0).Phi(), eventWeight);
        mumuHistNT01.h_photon_energy->Fill(Photon_vector.at(0).E(), eventWeight);
      }//SS mode
    }//one tight muon and loose-not-tight muon
   
  else if(LNTMuons.size()>1.0 and eventWeight > 0.0) {
     if(type=="Data" and signSelection=="SS" and ((LNTMuons.at(0).Charge)*(LNTMuons.at(1).Charge)==1)){
       mumu_SS_NT00+=eventWeight;
       mu1_transverse.SetMagPhi(LNTMuons.at(0).LepLV.Pt(), LNTMuons.at(0).LepLV.Phi());
       mu2_transverse.SetMagPhi(LNTMuons.at(1).LepLV.Pt(), LNTMuons.at(1).LepLV.Phi());
       mumuHistNT00.h_mu_pt_leading->Fill(LNTMuons.at(0).LepLV.Pt(), eventWeight);
       mumuHistNT00.h_mu_phi_leading->Fill(LNTMuons.at(0).LepLV.Phi(), eventWeight);
       mumuHistNT00.h_mu_eta_leading->Fill(LNTMuons.at(0).LepLV.Eta(), eventWeight);
       mumuHistNT00.h_mu_energy_leading->Fill(LNTMuons.at(0).LepLV.E(), eventWeight);
       mumuHistNT00.h_mu_pt_trailing->Fill(LNTMuons.at(1).LepLV.Pt(), eventWeight);
       mumuHistNT00.h_mu_phi_trailing->Fill(LNTMuons.at(1).LepLV.Phi(), eventWeight);
       mumuHistNT00.h_mu_eta_trailing->Fill(LNTMuons.at(1).LepLV.Eta(), eventWeight);
       mumuHistNT00.h_mu_energy_trailing->Fill(LNTMuons.at(1).LepLV.E(), eventWeight);
       mumuHistNT00.h_InvariantMass->Fill((LNTMuons.at(0).LepLV+LNTMuons.at(1).LepLV).M(), eventWeight);
       mumuHistNT00.h_InvariantMass_Ph->Fill((LNTMuons.at(0).LepLV+LNTMuons.at(1).LepLV+Photon_vector.at(0)).M(), eventWeight);
       mumuHistNT00.h_DeltaPhi_met_mu1->Fill(met_transverse.DeltaPhi(mu1_transverse), eventWeight);
       mumuHistNT00.h_DeltaPhi_met_mu2->Fill(met_transverse.DeltaPhi(mu2_transverse), eventWeight);
       mumuHistNT00.h_DeltaPhi_ph_mu1->Fill(ph_transverse.DeltaPhi(mu1_transverse), eventWeight);
       mumuHistNT00.h_DeltaPhi_ph_mu2->Fill(ph_transverse.DeltaPhi(mu2_transverse), eventWeight);
       mumuHistNT00.h_Isolation_mu1->Fill(LNTMuons.at(0).Isolation, eventWeight);
       mumuHistNT00.h_Isolation_mu2->Fill(LNTMuons.at(1).Isolation, eventWeight); 
       mumuHistNT00.h_nVertices->Fill(nVertices, eventWeight);
       mumuHistNT00.h_MET->Fill(MET, eventWeight);
       mumuHistNT00.h_caloMET->Fill(caloMET, eventWeight);
       mumuHistNT00.h_nJets->Fill(Jet_vector.size(), eventWeight);
       mumuHistNT00.h_HT->Fill(HT, eventWeight);
       if(Jet_vector.size()>0) mumuHistNT00.h_jet_pt_leading->Fill(Jet_vector.at(0).Pt(), eventWeight);
       if(Jet_vector.size()>1) mumuHistNT00.h_jet_pt_trailing->Fill(Jet_vector.at(1).Pt(), eventWeight);
       if(Jet_vector.size()>2) mumuHistNT00.h_jet_pt_3rd->Fill(Jet_vector.at(2).Pt(), eventWeight);
       if(Jet_vector.size()>3) mumuHistNT00.h_jet_pt_4th->Fill(Jet_vector.at(3).Pt(), eventWeight);
       if(Jet_vector.size()>4) mumuHistNT00.h_jet_pt_5th->Fill(Jet_vector.at(4).Pt(), eventWeight);
       if(Jet_vector.size()>5) mumuHistNT00.h_jet_pt_6th->Fill(Jet_vector.at(5).Pt(), eventWeight);
       if(Jet_vector.size()>0) mumuHistNT00.h_jet_eta_leading->Fill(Jet_vector.at(0).Eta(), eventWeight);
       if(Jet_vector.size()>1) mumuHistNT00.h_jet_eta_trailing->Fill(Jet_vector.at(1).Eta(), eventWeight);
       if(Jet_vector.size()>2) mumuHistNT00.h_jet_eta_3rd->Fill(Jet_vector.at(2).Eta(), eventWeight);
       if(Jet_vector.size()>3) mumuHistNT00.h_jet_eta_4th->Fill(Jet_vector.at(3).Eta(), eventWeight);
       if(Jet_vector.size()>4) mumuHistNT00.h_jet_eta_5th->Fill(Jet_vector.at(4).Eta(), eventWeight);
       if(Jet_vector.size()>5) mumuHistNT00.h_jet_eta_6th->Fill(Jet_vector.at(5).Eta(), eventWeight);
       if(Jet_vector.size()>0) mumuHistNT00.h_jet_phi_leading->Fill(Jet_vector.at(0).Phi(), eventWeight);
       if(Jet_vector.size()>1) mumuHistNT00.h_jet_phi_trailing->Fill(Jet_vector.at(1).Phi(), eventWeight);
       if(Jet_vector.size()>2) mumuHistNT00.h_jet_phi_3rd->Fill(Jet_vector.at(2).Phi(), eventWeight);
       if(Jet_vector.size()>3) mumuHistNT00.h_jet_phi_4th->Fill(Jet_vector.at(3).Phi(), eventWeight);
       if(Jet_vector.size()>4) mumuHistNT00.h_jet_phi_5th->Fill(Jet_vector.at(4).Phi(), eventWeight);
       if(Jet_vector.size()>5) mumuHistNT00.h_jet_phi_6th->Fill(Jet_vector.at(5).Phi(), eventWeight);
       if(Jet_vector.size()>0) mumuHistNT00.h_jet_energy_leading->Fill(Jet_vector.at(0).E(), eventWeight);
       if(Jet_vector.size()>1) mumuHistNT00.h_jet_energy_trailing->Fill(Jet_vector.at(1).E(), eventWeight);
       if(Jet_vector.size()>2) mumuHistNT00.h_jet_energy_3rd->Fill(Jet_vector.at(2).E(), eventWeight);
       if(Jet_vector.size()>3) mumuHistNT00.h_jet_energy_4th->Fill(Jet_vector.at(3).E(), eventWeight);
       if(Jet_vector.size()>4) mumuHistNT00.h_jet_energy_5th->Fill(Jet_vector.at(4).E(), eventWeight);
       if(Jet_vector.size()>5) mumuHistNT00.h_jet_energy_6th->Fill(Jet_vector.at(5).E(), eventWeight);
       mumuHistNT00.h_photon_pt->Fill(Photon_vector.at(0).Pt(), eventWeight);
       mumuHistNT00.h_photon_eta->Fill(Photon_vector.at(0).Eta(), eventWeight);
       mumuHistNT00.h_photon_phi->Fill(Photon_vector.at(0).Phi(), eventWeight);
       mumuHistNT00.h_photon_energy->Fill(Photon_vector.at(0).E(), eventWeight);
     }//SS mode
   }//one loose-not-tight muon and loose-not-tight muon
   

 //ElMu channel
   if( eventWeight > 0.0 and Muons.size()==1 and Electrons.size()==1) {
     if(type=="MC") eventWeight *= electronSF(Electrons.at(0).LepLV.Pt(), Electrons.at(0).LepLV.Eta())*muonSF(Muons.at(0).LepLV.Pt(), Muons.at(0).LepLV.Eta())*photonSF(Photon_vector.at(0).Pt(), Photon_vector.at(0).Eta());
     mu1_transverse.SetMagPhi(Muons.at(0).LepLV.Pt(), Muons.at(0).LepLV.Phi());
     el1_transverse.SetMagPhi(Electrons.at(0).LepLV.Pt(), Electrons.at(0).LepLV.Phi());
     if(signSelection=="OS" and ((Muons.at(0).Charge)*(Electrons.at(0).Charge)==-1)){
       if(type=="Data") cout << "Event = " << event << " Run = "<< run << " Lumi = " << lumi << endl;  
       emu_OS+=eventWeight;
       emuHist.h_mu_pt_leading->Fill(Muons.at(0).LepLV.Pt(), eventWeight);
       emuHist.h_mu_phi_leading->Fill(Muons.at(0).LepLV.Phi(), eventWeight);
       emuHist.h_mu_eta_leading->Fill(Muons.at(0).LepLV.Eta(), eventWeight);
       emuHist.h_mu_energy_leading->Fill(Muons.at(0).LepLV.E(), eventWeight);
       emuHist.h_el_pt_leading->Fill(Electrons.at(0).LepLV.Pt(), eventWeight);
       emuHist.h_el_phi_leading->Fill(Electrons.at(0).LepLV.Phi(), eventWeight);
       emuHist.h_el_eta_leading->Fill(Electrons.at(0).LepLV.Eta(), eventWeight);
       emuHist.h_el_energy_leading->Fill(Electrons.at(0).LepLV.E(), eventWeight);
       emuHist.h_InvariantMass->Fill((Muons.at(0).LepLV+Electrons.at(0).LepLV).M(), eventWeight);
       emuHist.h_InvariantMass_Ph->Fill((Muons.at(0).LepLV+Electrons.at(0).LepLV+Photon_vector.at(0)).M(), eventWeight);
       emuHist.h_DeltaPhi_met_mu1->Fill(met_transverse.DeltaPhi(mu1_transverse), eventWeight);
       emuHist.h_DeltaPhi_met_el1->Fill(met_transverse.DeltaPhi(el1_transverse), eventWeight);
       emuHist.h_DeltaPhi_ph_mu1->Fill(ph_transverse.DeltaPhi(mu1_transverse), eventWeight);
       emuHist.h_DeltaPhi_ph_el1->Fill(ph_transverse.DeltaPhi(el1_transverse), eventWeight);
       emuHist.h_Isolation_mu1->Fill(Muons.at(0).Isolation, eventWeight);
       emuHist.h_Isolation_el1->Fill(Electrons.at(0).Isolation, eventWeight);
       emuHist.h_nVertices->Fill(nVertices, eventWeight);
       emuHist.h_MET->Fill(MET, eventWeight);
       emuHist.h_caloMET->Fill(caloMET, eventWeight);
       emuHist.h_nJets->Fill(Jet_vector.size(), eventWeight);
       emuHist.h_HT->Fill(HT, eventWeight);
       if(Jet_vector.size()>0) emuHist.h_jet_pt_leading->Fill(Jet_vector.at(0).Pt(), eventWeight);
       if(Jet_vector.size()>1) emuHist.h_jet_pt_trailing->Fill(Jet_vector.at(1).Pt(), eventWeight);
       if(Jet_vector.size()>2) emuHist.h_jet_pt_3rd->Fill(Jet_vector.at(2).Pt(), eventWeight);
       if(Jet_vector.size()>3) emuHist.h_jet_pt_4th->Fill(Jet_vector.at(3).Pt(), eventWeight);
       if(Jet_vector.size()>4) emuHist.h_jet_pt_5th->Fill(Jet_vector.at(4).Pt(), eventWeight);
       if(Jet_vector.size()>5) emuHist.h_jet_pt_6th->Fill(Jet_vector.at(5).Pt(), eventWeight);
       if(Jet_vector.size()>0) emuHist.h_jet_eta_leading->Fill(Jet_vector.at(0).Eta(), eventWeight);
       if(Jet_vector.size()>1) emuHist.h_jet_eta_trailing->Fill(Jet_vector.at(1).Eta(), eventWeight);
       if(Jet_vector.size()>2) emuHist.h_jet_eta_3rd->Fill(Jet_vector.at(2).Eta(), eventWeight);
       if(Jet_vector.size()>3) emuHist.h_jet_eta_4th->Fill(Jet_vector.at(3).Eta(), eventWeight);
       if(Jet_vector.size()>4) emuHist.h_jet_eta_5th->Fill(Jet_vector.at(4).Eta(), eventWeight);
       if(Jet_vector.size()>5) emuHist.h_jet_eta_6th->Fill(Jet_vector.at(5).Eta(), eventWeight);
       if(Jet_vector.size()>0) emuHist.h_jet_phi_leading->Fill(Jet_vector.at(0).Phi(), eventWeight);
       if(Jet_vector.size()>1) emuHist.h_jet_phi_trailing->Fill(Jet_vector.at(1).Phi(), eventWeight);
       if(Jet_vector.size()>2) emuHist.h_jet_phi_3rd->Fill(Jet_vector.at(2).Phi(), eventWeight);
       if(Jet_vector.size()>3) emuHist.h_jet_phi_4th->Fill(Jet_vector.at(3).Phi(), eventWeight);
       if(Jet_vector.size()>4) emuHist.h_jet_phi_5th->Fill(Jet_vector.at(4).Phi(), eventWeight);
       if(Jet_vector.size()>5) emuHist.h_jet_phi_6th->Fill(Jet_vector.at(5).Phi(), eventWeight);
       if(Jet_vector.size()>0) emuHist.h_jet_energy_leading->Fill(Jet_vector.at(0).E(), eventWeight);
       if(Jet_vector.size()>1) emuHist.h_jet_energy_trailing->Fill(Jet_vector.at(1).E(), eventWeight);
       if(Jet_vector.size()>2) emuHist.h_jet_energy_3rd->Fill(Jet_vector.at(2).E(), eventWeight);
       if(Jet_vector.size()>3) emuHist.h_jet_energy_4th->Fill(Jet_vector.at(3).E(), eventWeight);
       if(Jet_vector.size()>4) emuHist.h_jet_energy_5th->Fill(Jet_vector.at(4).E(), eventWeight);
       if(Jet_vector.size()>5) emuHist.h_jet_energy_6th->Fill(Jet_vector.at(5).E(), eventWeight);
       emuHist.h_photon_pt->Fill(Photon_vector.at(0).Pt(), eventWeight);
       emuHist.h_photon_eta->Fill(Photon_vector.at(0).Eta(), eventWeight);
       emuHist.h_photon_phi->Fill(Photon_vector.at(0).Phi(), eventWeight);
       emuHist.h_photon_energy->Fill(Photon_vector.at(0).E(), eventWeight);
       if(MET > 25.0) MET_emu_OS+=eventWeight;
       if(HT<300) LowHT_emu_OS+=eventWeight;
     }//OS mode
    else if(signSelection=="SS" and ((Muons.at(0).Charge)*(Electrons.at(0).Charge)==+1)){
       if(type=="Data") cout << "Event = " << event << " Run = "<< run << " Lumi = " << lumi << endl;
       emu_SS+=eventWeight;
       emuHist.h_mu_pt_leading->Fill(Muons.at(0).LepLV.Pt(), eventWeight);
       emuHist.h_mu_phi_leading->Fill(Muons.at(0).LepLV.Phi(), eventWeight);
       emuHist.h_mu_eta_leading->Fill(Muons.at(0).LepLV.Eta(), eventWeight);
       emuHist.h_mu_energy_leading->Fill(Muons.at(0).LepLV.E(), eventWeight);
       emuHist.h_el_pt_leading->Fill(Electrons.at(0).LepLV.Pt(), eventWeight);
       emuHist.h_el_phi_leading->Fill(Electrons.at(0).LepLV.Phi(), eventWeight);
       emuHist.h_el_eta_leading->Fill(Electrons.at(0).LepLV.Eta(), eventWeight);
       emuHist.h_el_energy_leading->Fill(Electrons.at(0).LepLV.E(), eventWeight);
       emuHist.h_InvariantMass->Fill((Muons.at(0).LepLV+Electrons.at(0).LepLV).M(), eventWeight);
       emuHist.h_InvariantMass_Ph->Fill((Muons.at(0).LepLV+Electrons.at(0).LepLV+Photon_vector.at(0)).M(), eventWeight);
       emuHist.h_DeltaPhi_met_mu1->Fill(met_transverse.DeltaPhi(mu1_transverse), eventWeight);
       emuHist.h_DeltaPhi_met_el1->Fill(met_transverse.DeltaPhi(el1_transverse), eventWeight);
       emuHist.h_DeltaPhi_ph_mu1->Fill(ph_transverse.DeltaPhi(mu1_transverse), eventWeight);
       emuHist.h_DeltaPhi_ph_el1->Fill(ph_transverse.DeltaPhi(el1_transverse), eventWeight);
       emuHist.h_Isolation_mu1->Fill(Muons.at(0).Isolation, eventWeight);
       emuHist.h_Isolation_el1->Fill(Electrons.at(0).Isolation, eventWeight);
       emuHist.h_nVertices->Fill(nVertices, eventWeight);
       emuHist.h_MET->Fill(MET, eventWeight);
       emuHist.h_caloMET->Fill(caloMET, eventWeight);
       emuHist.h_nJets->Fill(Jet_vector.size(), eventWeight);
       emuHist.h_HT->Fill(HT, eventWeight);
       if(Jet_vector.size()>0) emuHist.h_jet_pt_leading->Fill(Jet_vector.at(0).Pt(), eventWeight);
       if(Jet_vector.size()>1) emuHist.h_jet_pt_trailing->Fill(Jet_vector.at(1).Pt(), eventWeight);
       if(Jet_vector.size()>2) emuHist.h_jet_pt_3rd->Fill(Jet_vector.at(2).Pt(), eventWeight);
       if(Jet_vector.size()>3) emuHist.h_jet_pt_4th->Fill(Jet_vector.at(3).Pt(), eventWeight);
       if(Jet_vector.size()>4) emuHist.h_jet_pt_5th->Fill(Jet_vector.at(4).Pt(), eventWeight);
       if(Jet_vector.size()>5) emuHist.h_jet_pt_6th->Fill(Jet_vector.at(5).Pt(), eventWeight);
       if(Jet_vector.size()>0) emuHist.h_jet_eta_leading->Fill(Jet_vector.at(0).Eta(), eventWeight);
       if(Jet_vector.size()>1) emuHist.h_jet_eta_trailing->Fill(Jet_vector.at(1).Eta(), eventWeight);
       if(Jet_vector.size()>2) emuHist.h_jet_eta_3rd->Fill(Jet_vector.at(2).Eta(), eventWeight);
       if(Jet_vector.size()>3) emuHist.h_jet_eta_4th->Fill(Jet_vector.at(3).Eta(), eventWeight);
       if(Jet_vector.size()>4) emuHist.h_jet_eta_5th->Fill(Jet_vector.at(4).Eta(), eventWeight);
       if(Jet_vector.size()>5) emuHist.h_jet_eta_6th->Fill(Jet_vector.at(5).Eta(), eventWeight);
       if(Jet_vector.size()>0) emuHist.h_jet_phi_leading->Fill(Jet_vector.at(0).Phi(), eventWeight);
       if(Jet_vector.size()>1) emuHist.h_jet_phi_trailing->Fill(Jet_vector.at(1).Phi(), eventWeight);
       if(Jet_vector.size()>2) emuHist.h_jet_phi_3rd->Fill(Jet_vector.at(2).Phi(), eventWeight);
       if(Jet_vector.size()>3) emuHist.h_jet_phi_4th->Fill(Jet_vector.at(3).Phi(), eventWeight);
       if(Jet_vector.size()>4) emuHist.h_jet_phi_5th->Fill(Jet_vector.at(4).Phi(), eventWeight);
       if(Jet_vector.size()>5) emuHist.h_jet_phi_6th->Fill(Jet_vector.at(5).Phi(), eventWeight);
       if(Jet_vector.size()>0) emuHist.h_jet_energy_leading->Fill(Jet_vector.at(0).E(), eventWeight);
       if(Jet_vector.size()>1) emuHist.h_jet_energy_trailing->Fill(Jet_vector.at(1).E(), eventWeight);
       if(Jet_vector.size()>2) emuHist.h_jet_energy_3rd->Fill(Jet_vector.at(2).E(), eventWeight);
       if(Jet_vector.size()>3) emuHist.h_jet_energy_4th->Fill(Jet_vector.at(3).E(), eventWeight);
       if(Jet_vector.size()>4) emuHist.h_jet_energy_5th->Fill(Jet_vector.at(4).E(), eventWeight);
       if(Jet_vector.size()>5) emuHist.h_jet_energy_6th->Fill(Jet_vector.at(5).E(), eventWeight);
       emuHist.h_photon_pt->Fill(Photon_vector.at(0).Pt(), eventWeight);
       emuHist.h_photon_eta->Fill(Photon_vector.at(0).Eta(), eventWeight);
       emuHist.h_photon_phi->Fill(Photon_vector.at(0).Phi(), eventWeight);
       emuHist.h_photon_energy->Fill(Photon_vector.at(0).E(), eventWeight);
       if(MET > 25.0) MET_emu_SS+=eventWeight;
       if(HT<300) LowHT_emu_SS+=eventWeight;
     }//SS mode     
   }//El-Mu channel
//Non-prompt background computation:
   else if(Muons.size() > 0.0 and LNTElectrons.size() > 0.0 and eventWeight > 0.0) {
     if(type=="Data" and signSelection=="SS" and ((Muons.at(0).Charge)*(LNTElectrons.at(0).Charge)==1)){
       emu_SS_NT01+=eventWeight;
       mu1_transverse.SetMagPhi(Muons.at(0).LepLV.Pt(), Muons.at(0).LepLV.Phi());
       el1_transverse.SetMagPhi(LNTElectrons.at(0).LepLV.Pt(), LNTElectrons.at(0).LepLV.Phi());
       emuHistNT01.h_mu_pt_leading->Fill(Muons.at(0).LepLV.Pt(), eventWeight);
       emuHistNT01.h_mu_phi_leading->Fill(Muons.at(0).LepLV.Phi(), eventWeight);
       emuHistNT01.h_mu_eta_leading->Fill(Muons.at(0).LepLV.Eta(), eventWeight);
       emuHistNT01.h_mu_energy_leading->Fill(Muons.at(0).LepLV.E(), eventWeight);
       emuHistNT01.h_el_pt_leading->Fill(LNTElectrons.at(0).LepLV.Pt(), eventWeight);
       emuHistNT01.h_el_phi_leading->Fill(LNTElectrons.at(0).LepLV.Phi(), eventWeight);
       emuHistNT01.h_el_eta_leading->Fill(LNTElectrons.at(0).LepLV.Eta(), eventWeight);
       emuHistNT01.h_el_energy_leading->Fill(LNTElectrons.at(0).LepLV.E(), eventWeight);
       emuHistNT01.h_InvariantMass->Fill((Muons.at(0).LepLV+LNTElectrons.at(0).LepLV).M(), eventWeight);
       emuHistNT01.h_InvariantMass_Ph->Fill((Muons.at(0).LepLV+LNTElectrons.at(0).LepLV+Photon_vector.at(0)).M(), eventWeight);
       emuHistNT01.h_DeltaPhi_met_mu1->Fill(met_transverse.DeltaPhi(mu1_transverse), eventWeight);
       emuHistNT01.h_DeltaPhi_met_el1->Fill(met_transverse.DeltaPhi(el1_transverse), eventWeight);
       emuHistNT01.h_DeltaPhi_ph_mu1->Fill(ph_transverse.DeltaPhi(mu1_transverse), eventWeight);
       emuHistNT01.h_DeltaPhi_ph_el1->Fill(ph_transverse.DeltaPhi(el1_transverse), eventWeight);
       emuHistNT01.h_Isolation_mu1->Fill(Muons.at(0).Isolation, eventWeight);
       emuHistNT01.h_Isolation_el1->Fill(LNTElectrons.at(0).Isolation, eventWeight);
       emuHistNT01.h_nVertices->Fill(nVertices, eventWeight);
       emuHistNT01.h_MET->Fill(MET, eventWeight);
       emuHistNT01.h_caloMET->Fill(caloMET, eventWeight);
       emuHistNT01.h_nJets->Fill(Jet_vector.size(), eventWeight);
       emuHistNT01.h_HT->Fill(HT, eventWeight);
       if(Jet_vector.size()>0) emuHistNT01.h_jet_pt_leading->Fill(Jet_vector.at(0).Pt(), eventWeight);
       if(Jet_vector.size()>1) emuHistNT01.h_jet_pt_trailing->Fill(Jet_vector.at(1).Pt(), eventWeight);
       if(Jet_vector.size()>2) emuHistNT01.h_jet_pt_3rd->Fill(Jet_vector.at(2).Pt(), eventWeight);
       if(Jet_vector.size()>3) emuHistNT01.h_jet_pt_4th->Fill(Jet_vector.at(3).Pt(), eventWeight);
       if(Jet_vector.size()>4) emuHistNT01.h_jet_pt_5th->Fill(Jet_vector.at(4).Pt(), eventWeight);
       if(Jet_vector.size()>5) emuHistNT01.h_jet_pt_6th->Fill(Jet_vector.at(5).Pt(), eventWeight);
       if(Jet_vector.size()>0) emuHistNT01.h_jet_eta_leading->Fill(Jet_vector.at(0).Eta(), eventWeight);
       if(Jet_vector.size()>1) emuHistNT01.h_jet_eta_trailing->Fill(Jet_vector.at(1).Eta(), eventWeight);
       if(Jet_vector.size()>2) emuHistNT01.h_jet_eta_3rd->Fill(Jet_vector.at(2).Eta(), eventWeight);
       if(Jet_vector.size()>3) emuHistNT01.h_jet_eta_4th->Fill(Jet_vector.at(3).Eta(), eventWeight);
       if(Jet_vector.size()>4) emuHistNT01.h_jet_eta_5th->Fill(Jet_vector.at(4).Eta(), eventWeight);
       if(Jet_vector.size()>5) emuHistNT01.h_jet_eta_6th->Fill(Jet_vector.at(5).Eta(), eventWeight);
       if(Jet_vector.size()>0) emuHistNT01.h_jet_phi_leading->Fill(Jet_vector.at(0).Phi(), eventWeight);
       if(Jet_vector.size()>1) emuHistNT01.h_jet_phi_trailing->Fill(Jet_vector.at(1).Phi(), eventWeight);
       if(Jet_vector.size()>2) emuHistNT01.h_jet_phi_3rd->Fill(Jet_vector.at(2).Phi(), eventWeight);
       if(Jet_vector.size()>3) emuHistNT01.h_jet_phi_4th->Fill(Jet_vector.at(3).Phi(), eventWeight);
       if(Jet_vector.size()>4) emuHistNT01.h_jet_phi_5th->Fill(Jet_vector.at(4).Phi(), eventWeight);
       if(Jet_vector.size()>5) emuHistNT01.h_jet_phi_6th->Fill(Jet_vector.at(5).Phi(), eventWeight);    
       if(Jet_vector.size()>0) emuHistNT01.h_jet_energy_leading->Fill(Jet_vector.at(0).E(), eventWeight);
       if(Jet_vector.size()>1) emuHistNT01.h_jet_energy_trailing->Fill(Jet_vector.at(1).E(), eventWeight);
       if(Jet_vector.size()>2) emuHistNT01.h_jet_energy_3rd->Fill(Jet_vector.at(2).E(), eventWeight);
       if(Jet_vector.size()>3) emuHistNT01.h_jet_energy_4th->Fill(Jet_vector.at(3).E(), eventWeight);
       if(Jet_vector.size()>4) emuHistNT01.h_jet_energy_5th->Fill(Jet_vector.at(4).E(), eventWeight);
       if(Jet_vector.size()>5) emuHistNT01.h_jet_energy_6th->Fill(Jet_vector.at(5).E(), eventWeight);
       emuHistNT01.h_photon_pt->Fill(Photon_vector.at(0).Pt(), eventWeight);
       emuHistNT01.h_photon_eta->Fill(Photon_vector.at(0).Eta(), eventWeight);
       emuHistNT01.h_photon_phi->Fill(Photon_vector.at(0).Phi(), eventWeight);
       emuHistNT01.h_photon_energy->Fill(Photon_vector.at(0).E(), eventWeight);
     }//SS mode
   }//construction of tight muon, loose electron
  else if(LNTMuons.size() > 0.0 and Electrons.size() > 0.0 and eventWeight > 0.0) {
     if(type=="Data" and signSelection=="SS" and ((LNTMuons.at(0).Charge)*(Electrons.at(0).Charge)==1)){
       emu_SS_NT10+=eventWeight;
       mu1_transverse.SetMagPhi(LNTMuons.at(0).LepLV.Pt(), LNTMuons.at(0).LepLV.Phi());
       el1_transverse.SetMagPhi(Electrons.at(0).LepLV.Pt(), Electrons.at(0).LepLV.Phi());
       emuHistNT10.h_mu_pt_leading->Fill(LNTMuons.at(0).LepLV.Pt(), eventWeight);
       emuHistNT10.h_mu_phi_leading->Fill(LNTMuons.at(0).LepLV.Phi(), eventWeight);
       emuHistNT10.h_mu_eta_leading->Fill(LNTMuons.at(0).LepLV.Eta(), eventWeight);
       emuHistNT10.h_mu_energy_leading->Fill(LNTMuons.at(0).LepLV.E(), eventWeight);
       emuHistNT10.h_el_pt_leading->Fill(Electrons.at(0).LepLV.Pt(), eventWeight);
       emuHistNT10.h_el_phi_leading->Fill(Electrons.at(0).LepLV.Phi(), eventWeight);
       emuHistNT10.h_el_eta_leading->Fill(Electrons.at(0).LepLV.Eta(), eventWeight);
       emuHistNT10.h_el_energy_leading->Fill(Electrons.at(0).LepLV.E(), eventWeight);
       emuHistNT10.h_InvariantMass->Fill((LNTMuons.at(0).LepLV+Electrons.at(0).LepLV).M(), eventWeight);
       emuHistNT10.h_InvariantMass_Ph->Fill((LNTMuons.at(0).LepLV+Electrons.at(0).LepLV+Photon_vector.at(0)).M(), eventWeight);
       emuHistNT10.h_DeltaPhi_met_mu1->Fill(met_transverse.DeltaPhi(mu1_transverse), eventWeight);
       emuHistNT10.h_DeltaPhi_met_el1->Fill(met_transverse.DeltaPhi(el1_transverse), eventWeight);
       emuHistNT10.h_DeltaPhi_ph_mu1->Fill(ph_transverse.DeltaPhi(mu1_transverse), eventWeight);
       emuHistNT10.h_DeltaPhi_ph_el1->Fill(ph_transverse.DeltaPhi(el1_transverse), eventWeight);
       emuHistNT10.h_Isolation_mu1->Fill(LNTMuons.at(0).Isolation, eventWeight);
       emuHistNT10.h_Isolation_el1->Fill(Electrons.at(0).Isolation, eventWeight);
       emuHistNT10.h_nVertices->Fill(nVertices, eventWeight);
       emuHistNT10.h_MET->Fill(MET, eventWeight);
       emuHistNT10.h_caloMET->Fill(caloMET, eventWeight);
       emuHistNT10.h_nJets->Fill(Jet_vector.size(), eventWeight);
       emuHistNT10.h_HT->Fill(HT, eventWeight);
       if(Jet_vector.size()>0) emuHistNT10.h_jet_pt_leading->Fill(Jet_vector.at(0).Pt(), eventWeight);
       if(Jet_vector.size()>1) emuHistNT10.h_jet_pt_trailing->Fill(Jet_vector.at(1).Pt(), eventWeight);
       if(Jet_vector.size()>2) emuHistNT10.h_jet_pt_3rd->Fill(Jet_vector.at(2).Pt(), eventWeight);
       if(Jet_vector.size()>3) emuHistNT10.h_jet_pt_4th->Fill(Jet_vector.at(3).Pt(), eventWeight);
       if(Jet_vector.size()>4) emuHistNT10.h_jet_pt_5th->Fill(Jet_vector.at(4).Pt(), eventWeight);
       if(Jet_vector.size()>5) emuHistNT10.h_jet_pt_6th->Fill(Jet_vector.at(5).Pt(), eventWeight);
       if(Jet_vector.size()>0) emuHistNT10.h_jet_eta_leading->Fill(Jet_vector.at(0).Eta(), eventWeight);
       if(Jet_vector.size()>1) emuHistNT10.h_jet_eta_trailing->Fill(Jet_vector.at(1).Eta(), eventWeight);
       if(Jet_vector.size()>2) emuHistNT10.h_jet_eta_3rd->Fill(Jet_vector.at(2).Eta(), eventWeight);
       if(Jet_vector.size()>3) emuHistNT10.h_jet_eta_4th->Fill(Jet_vector.at(3).Eta(), eventWeight);
       if(Jet_vector.size()>4) emuHistNT10.h_jet_eta_5th->Fill(Jet_vector.at(4).Eta(), eventWeight);
       if(Jet_vector.size()>5) emuHistNT10.h_jet_eta_6th->Fill(Jet_vector.at(5).Eta(), eventWeight);
       if(Jet_vector.size()>0) emuHistNT10.h_jet_phi_leading->Fill(Jet_vector.at(0).Phi(), eventWeight);
       if(Jet_vector.size()>1) emuHistNT10.h_jet_phi_trailing->Fill(Jet_vector.at(1).Phi(), eventWeight);
       if(Jet_vector.size()>2) emuHistNT10.h_jet_phi_3rd->Fill(Jet_vector.at(2).Phi(), eventWeight);
       if(Jet_vector.size()>3) emuHistNT10.h_jet_phi_4th->Fill(Jet_vector.at(3).Phi(), eventWeight);
       if(Jet_vector.size()>4) emuHistNT10.h_jet_phi_5th->Fill(Jet_vector.at(4).Phi(), eventWeight);
       if(Jet_vector.size()>5) emuHistNT10.h_jet_phi_6th->Fill(Jet_vector.at(5).Phi(), eventWeight);
       if(Jet_vector.size()>0) emuHistNT10.h_jet_energy_leading->Fill(Jet_vector.at(0).E(), eventWeight);
       if(Jet_vector.size()>1) emuHistNT10.h_jet_energy_trailing->Fill(Jet_vector.at(1).E(), eventWeight);
       if(Jet_vector.size()>2) emuHistNT10.h_jet_energy_3rd->Fill(Jet_vector.at(2).E(), eventWeight);
       if(Jet_vector.size()>3) emuHistNT10.h_jet_energy_4th->Fill(Jet_vector.at(3).E(), eventWeight);
       if(Jet_vector.size()>4) emuHistNT10.h_jet_energy_5th->Fill(Jet_vector.at(4).E(), eventWeight);
       if(Jet_vector.size()>5) emuHistNT10.h_jet_energy_6th->Fill(Jet_vector.at(5).E(), eventWeight);
       emuHistNT10.h_photon_pt->Fill(Photon_vector.at(0).Pt(), eventWeight);
       emuHistNT10.h_photon_eta->Fill(Photon_vector.at(0).Eta(), eventWeight);
       emuHistNT10.h_photon_phi->Fill(Photon_vector.at(0).Phi(), eventWeight);
       emuHistNT10.h_photon_energy->Fill(Photon_vector.at(0).E(), eventWeight);
     }//SS mode
   }//construction of loose muon, tight electron
      
else if(LNTMuons.size() > 0.0 and LNTElectrons.size() > 0.0 and eventWeight > 0.0) {
     if(type=="Data" and signSelection=="SS" and ((LNTMuons.at(0).Charge)*(LNTElectrons.at(0).Charge)==1)){
       emu_SS_NT00+=eventWeight;
       mu1_transverse.SetMagPhi(LNTMuons.at(0).LepLV.Pt(), LNTMuons.at(0).LepLV.Phi());
       el1_transverse.SetMagPhi(LNTElectrons.at(0).LepLV.Pt(), LNTElectrons.at(0).LepLV.Phi());
       emuHistNT00.h_mu_pt_leading->Fill(LNTMuons.at(0).LepLV.Pt(), eventWeight);
       emuHistNT00.h_mu_phi_leading->Fill(LNTMuons.at(0).LepLV.Phi(), eventWeight);
       emuHistNT00.h_mu_eta_leading->Fill(LNTMuons.at(0).LepLV.Eta(), eventWeight);
       emuHistNT00.h_mu_energy_leading->Fill(LNTMuons.at(0).LepLV.E(), eventWeight);
       emuHistNT00.h_el_pt_leading->Fill(LNTElectrons.at(0).LepLV.Pt(), eventWeight);
       emuHistNT00.h_el_phi_leading->Fill(LNTElectrons.at(0).LepLV.Phi(), eventWeight);
       emuHistNT00.h_el_eta_leading->Fill(LNTElectrons.at(0).LepLV.Eta(), eventWeight);
       emuHistNT00.h_el_energy_leading->Fill(LNTElectrons.at(0).LepLV.E(), eventWeight);
       emuHistNT00.h_InvariantMass->Fill((LNTMuons.at(0).LepLV+LNTElectrons.at(0).LepLV).M(), eventWeight);
       emuHistNT00.h_InvariantMass_Ph->Fill((LNTMuons.at(0).LepLV+LNTElectrons.at(0).LepLV+Photon_vector.at(0)).M(), eventWeight);
       emuHistNT00.h_DeltaPhi_met_mu1->Fill(met_transverse.DeltaPhi(mu1_transverse), eventWeight);
       emuHistNT00.h_DeltaPhi_met_el1->Fill(met_transverse.DeltaPhi(el1_transverse), eventWeight);
       emuHistNT00.h_DeltaPhi_ph_mu1->Fill(ph_transverse.DeltaPhi(mu1_transverse), eventWeight);
       emuHistNT00.h_DeltaPhi_ph_el1->Fill(ph_transverse.DeltaPhi(el1_transverse), eventWeight);
       emuHistNT00.h_Isolation_mu1->Fill(LNTMuons.at(0).Isolation, eventWeight);
       emuHistNT00.h_Isolation_el1->Fill(LNTElectrons.at(0).Isolation, eventWeight);
       emuHistNT00.h_nVertices->Fill(nVertices, eventWeight);
       emuHistNT00.h_MET->Fill(MET, eventWeight);
       emuHistNT00.h_caloMET->Fill(caloMET, eventWeight);
       emuHistNT00.h_nJets->Fill(Jet_vector.size(), eventWeight);
       emuHistNT00.h_HT->Fill(HT, eventWeight);
       if(Jet_vector.size()>0) emuHistNT00.h_jet_pt_leading->Fill(Jet_vector.at(0).Pt(), eventWeight);
       if(Jet_vector.size()>1) emuHistNT00.h_jet_pt_trailing->Fill(Jet_vector.at(1).Pt(), eventWeight);
       if(Jet_vector.size()>2) emuHistNT00.h_jet_pt_3rd->Fill(Jet_vector.at(2).Pt(), eventWeight);
       if(Jet_vector.size()>3) emuHistNT00.h_jet_pt_4th->Fill(Jet_vector.at(3).Pt(), eventWeight);
       if(Jet_vector.size()>4) emuHistNT00.h_jet_pt_5th->Fill(Jet_vector.at(4).Pt(), eventWeight);
       if(Jet_vector.size()>5) emuHistNT00.h_jet_pt_6th->Fill(Jet_vector.at(5).Pt(), eventWeight);
       if(Jet_vector.size()>0) emuHistNT00.h_jet_eta_leading->Fill(Jet_vector.at(0).Eta(), eventWeight);
       if(Jet_vector.size()>1) emuHistNT00.h_jet_eta_trailing->Fill(Jet_vector.at(1).Eta(), eventWeight);
       if(Jet_vector.size()>2) emuHistNT00.h_jet_eta_3rd->Fill(Jet_vector.at(2).Eta(), eventWeight);
       if(Jet_vector.size()>3) emuHistNT00.h_jet_eta_4th->Fill(Jet_vector.at(3).Eta(), eventWeight);
       if(Jet_vector.size()>4) emuHistNT00.h_jet_eta_5th->Fill(Jet_vector.at(4).Eta(), eventWeight);
       if(Jet_vector.size()>5) emuHistNT00.h_jet_eta_6th->Fill(Jet_vector.at(5).Eta(), eventWeight);
       if(Jet_vector.size()>0) emuHistNT00.h_jet_phi_leading->Fill(Jet_vector.at(0).Phi(), eventWeight);
       if(Jet_vector.size()>1) emuHistNT00.h_jet_phi_trailing->Fill(Jet_vector.at(1).Phi(), eventWeight);
       if(Jet_vector.size()>2) emuHistNT00.h_jet_phi_3rd->Fill(Jet_vector.at(2).Phi(), eventWeight);
       if(Jet_vector.size()>3) emuHistNT00.h_jet_phi_4th->Fill(Jet_vector.at(3).Phi(), eventWeight);
       if(Jet_vector.size()>4) emuHistNT00.h_jet_phi_5th->Fill(Jet_vector.at(4).Phi(), eventWeight);
       if(Jet_vector.size()>5) emuHistNT00.h_jet_phi_6th->Fill(Jet_vector.at(5).Phi(), eventWeight);
       if(Jet_vector.size()>0) emuHistNT00.h_jet_energy_leading->Fill(Jet_vector.at(0).E(), eventWeight);
       if(Jet_vector.size()>1) emuHistNT00.h_jet_energy_trailing->Fill(Jet_vector.at(1).E(), eventWeight);
       if(Jet_vector.size()>2) emuHistNT00.h_jet_energy_3rd->Fill(Jet_vector.at(2).E(), eventWeight);
       if(Jet_vector.size()>3) emuHistNT00.h_jet_energy_4th->Fill(Jet_vector.at(3).E(), eventWeight);
       if(Jet_vector.size()>4) emuHistNT00.h_jet_energy_5th->Fill(Jet_vector.at(4).E(), eventWeight);
       if(Jet_vector.size()>5) emuHistNT00.h_jet_energy_6th->Fill(Jet_vector.at(5).E(), eventWeight);
       emuHistNT00.h_photon_pt->Fill(Photon_vector.at(0).Pt(), eventWeight);
       emuHistNT00.h_photon_eta->Fill(Photon_vector.at(0).Eta(), eventWeight);
       emuHistNT00.h_photon_phi->Fill(Photon_vector.at(0).Phi(), eventWeight);
       emuHistNT00.h_photon_energy->Fill(Photon_vector.at(0).E(), eventWeight);       
     }//SS mode
   }//construction of loose muon, loose electron

  if(eventWeight > 0.0 and Electrons.size()==2){
    if(type=="MC") eventWeight *= electronSF(Electrons.at(0).LepLV.Pt(), Electrons.at(0).LepLV.Eta())*electronSF(Electrons.at(1).LepLV.Pt(), Electrons.at(1).LepLV.Eta())*photonSF(Photon_vector.at(0).Pt(), Photon_vector.at(0).Eta());
    el1_transverse.SetMagPhi(Electrons.at(0).LepLV.Pt(), Electrons.at(0).LepLV.Phi());
    el2_transverse.SetMagPhi(Electrons.at(1).LepLV.Pt(), Electrons.at(1).LepLV.Phi());
    if(signSelection=="OS" and (Electrons.at(0).Charge*Electrons.at(1).Charge)==-1) {
      ee_OS+=eventWeight;
      eeHist.h_el_pt_leading->Fill(Electrons.at(0).LepLV.Pt(), eventWeight);
      eeHist.h_el_phi_leading->Fill(Electrons.at(0).LepLV.Phi(), eventWeight);
      eeHist.h_el_eta_leading->Fill(Electrons.at(0).LepLV.Eta(), eventWeight);
      eeHist.h_el_energy_leading->Fill(Electrons.at(0).LepLV.E(), eventWeight);
      eeHist.h_el_pt_trailing->Fill(Electrons.at(1).LepLV.Pt(), eventWeight);
      eeHist.h_el_phi_trailing->Fill(Electrons.at(1).LepLV.Phi(), eventWeight);
      eeHist.h_el_eta_trailing->Fill(Electrons.at(1).LepLV.Eta(), eventWeight);
      eeHist.h_el_energy_trailing->Fill(Electrons.at(1).LepLV.E(), eventWeight);
      eeHist.h_InvariantMass->Fill((Electrons.at(0).LepLV+Electrons.at(1).LepLV).M(), eventWeight);
      eeHist.h_InvariantMass_Ph->Fill((Electrons.at(0).LepLV+Electrons.at(1).LepLV+Photon_vector.at(0)).M(), eventWeight);
      eeHist.h_DeltaPhi_met_el1->Fill(met_transverse.DeltaPhi(el1_transverse), eventWeight);
      eeHist.h_DeltaPhi_met_el2->Fill(met_transverse.DeltaPhi(el2_transverse), eventWeight); 
      eeHist.h_DeltaPhi_ph_el1->Fill(ph_transverse.DeltaPhi(el1_transverse), eventWeight);
      eeHist.h_DeltaPhi_ph_el2->Fill(ph_transverse.DeltaPhi(el2_transverse), eventWeight);
      eeHist.h_Isolation_el1->Fill(Electrons.at(0).Isolation, eventWeight);
      eeHist.h_Isolation_el2->Fill(Electrons.at(1).Isolation, eventWeight);
      eeHist.h_nVertices->Fill(nVertices, eventWeight);
      eeHist.h_MET->Fill(MET, eventWeight);
      eeHist.h_caloMET->Fill(caloMET, eventWeight);
      eeHist.h_nJets->Fill(Jet_vector.size(), eventWeight);
      eeHist.h_HT->Fill(HT, eventWeight);
      if(Jet_vector.size()>0) eeHist.h_jet_pt_leading->Fill(Jet_vector.at(0).Pt(), eventWeight);
      if(Jet_vector.size()>1) eeHist.h_jet_pt_trailing->Fill(Jet_vector.at(1).Pt(), eventWeight);
      if(Jet_vector.size()>2) eeHist.h_jet_pt_3rd->Fill(Jet_vector.at(2).Pt(), eventWeight);
      if(Jet_vector.size()>3) eeHist.h_jet_pt_4th->Fill(Jet_vector.at(3).Pt(), eventWeight);
      if(Jet_vector.size()>4) eeHist.h_jet_pt_5th->Fill(Jet_vector.at(4).Pt(), eventWeight);
      if(Jet_vector.size()>5) eeHist.h_jet_pt_6th->Fill(Jet_vector.at(5).Pt(), eventWeight);
      if(Jet_vector.size()>0) eeHist.h_jet_eta_leading->Fill(Jet_vector.at(0).Eta(), eventWeight);
      if(Jet_vector.size()>1) eeHist.h_jet_eta_trailing->Fill(Jet_vector.at(1).Eta(), eventWeight);
      if(Jet_vector.size()>2) eeHist.h_jet_eta_3rd->Fill(Jet_vector.at(2).Eta(), eventWeight);
      if(Jet_vector.size()>3) eeHist.h_jet_eta_4th->Fill(Jet_vector.at(3).Eta(), eventWeight);
      if(Jet_vector.size()>4) eeHist.h_jet_eta_5th->Fill(Jet_vector.at(4).Eta(), eventWeight);
      if(Jet_vector.size()>5) eeHist.h_jet_eta_6th->Fill(Jet_vector.at(5).Eta(), eventWeight);
      if(Jet_vector.size()>0) eeHist.h_jet_phi_leading->Fill(Jet_vector.at(0).Phi(), eventWeight);
      if(Jet_vector.size()>1) eeHist.h_jet_phi_trailing->Fill(Jet_vector.at(1).Phi(), eventWeight);
      if(Jet_vector.size()>2) eeHist.h_jet_phi_3rd->Fill(Jet_vector.at(2).Phi(), eventWeight);
      if(Jet_vector.size()>3) eeHist.h_jet_phi_4th->Fill(Jet_vector.at(3).Phi(), eventWeight);
      if(Jet_vector.size()>4) eeHist.h_jet_phi_5th->Fill(Jet_vector.at(4).Phi(), eventWeight);
      if(Jet_vector.size()>5) eeHist.h_jet_phi_6th->Fill(Jet_vector.at(5).Phi(), eventWeight);
      if(Jet_vector.size()>0) eeHist.h_jet_energy_leading->Fill(Jet_vector.at(0).E(), eventWeight);
      if(Jet_vector.size()>1) eeHist.h_jet_energy_trailing->Fill(Jet_vector.at(1).E(), eventWeight);
      if(Jet_vector.size()>2) eeHist.h_jet_energy_3rd->Fill(Jet_vector.at(2).E(), eventWeight);
      if(Jet_vector.size()>3) eeHist.h_jet_energy_4th->Fill(Jet_vector.at(3).E(), eventWeight);
      if(Jet_vector.size()>4) eeHist.h_jet_energy_5th->Fill(Jet_vector.at(4).E(), eventWeight);
      if(Jet_vector.size()>5) eeHist.h_jet_energy_6th->Fill(Jet_vector.at(5).E(), eventWeight);
      eeHist.h_photon_pt->Fill(Photon_vector.at(0).Pt(), eventWeight);
      eeHist.h_photon_eta->Fill(Photon_vector.at(0).Eta(), eventWeight);
      eeHist.h_photon_phi->Fill(Photon_vector.at(0).Phi(), eventWeight);
      eeHist.h_photon_energy->Fill(Photon_vector.at(0).E(), eventWeight);
    }//OS mode
    if(signSelection=="SS" and (Electrons.at(0).Charge*Electrons.at(1).Charge)==1) {
      ee_SS+=eventWeight;
      eeHist.h_el_pt_leading->Fill(Electrons.at(0).LepLV.Pt(), eventWeight);
      eeHist.h_el_phi_leading->Fill(Electrons.at(0).LepLV.Phi(), eventWeight);
      eeHist.h_el_eta_leading->Fill(Electrons.at(0).LepLV.Eta(), eventWeight);
      eeHist.h_el_energy_leading->Fill(Electrons.at(0).LepLV.E(), eventWeight);
      eeHist.h_el_pt_trailing->Fill(Electrons.at(1).LepLV.Pt(), eventWeight);
      eeHist.h_el_phi_trailing->Fill(Electrons.at(1).LepLV.Phi(), eventWeight);
      eeHist.h_el_eta_trailing->Fill(Electrons.at(1).LepLV.Eta(), eventWeight);
      eeHist.h_el_energy_trailing->Fill(Electrons.at(1).LepLV.E(), eventWeight);
      eeHist.h_InvariantMass->Fill((Electrons.at(0).LepLV+Electrons.at(1).LepLV).M(), eventWeight);
      eeHist.h_InvariantMass_Ph->Fill((Electrons.at(0).LepLV+Electrons.at(0).LepLV+Photon_vector.at(0)).M(), eventWeight);
      eeHist.h_DeltaPhi_met_el1->Fill(met_transverse.DeltaPhi(el1_transverse), eventWeight);
      eeHist.h_DeltaPhi_met_el2->Fill(met_transverse.DeltaPhi(el2_transverse), eventWeight);
      eeHist.h_DeltaPhi_ph_el1->Fill(ph_transverse.DeltaPhi(el1_transverse), eventWeight);
      eeHist.h_DeltaPhi_ph_el2->Fill(ph_transverse.DeltaPhi(el2_transverse), eventWeight);
      eeHist.h_Isolation_el1->Fill(Electrons.at(0).Isolation, eventWeight);
      eeHist.h_Isolation_el2->Fill(Electrons.at(1).Isolation, eventWeight);
      eeHist.h_nVertices->Fill(nVertices, eventWeight);
      eeHist.h_MET->Fill(MET, eventWeight);
      eeHist.h_caloMET->Fill(caloMET, eventWeight);
      eeHist.h_nJets->Fill(Jet_vector.size(), eventWeight);
      eeHist.h_HT->Fill(HT, eventWeight);
      if(Jet_vector.size()>0) eeHist.h_jet_pt_leading->Fill(Jet_vector.at(0).Pt(), eventWeight);
      if(Jet_vector.size()>1) eeHist.h_jet_pt_trailing->Fill(Jet_vector.at(1).Pt(), eventWeight);
      if(Jet_vector.size()>2) eeHist.h_jet_pt_3rd->Fill(Jet_vector.at(2).Pt(), eventWeight);
      if(Jet_vector.size()>3) eeHist.h_jet_pt_4th->Fill(Jet_vector.at(3).Pt(), eventWeight);
      if(Jet_vector.size()>4) eeHist.h_jet_pt_5th->Fill(Jet_vector.at(4).Pt(), eventWeight);
      if(Jet_vector.size()>5) eeHist.h_jet_pt_6th->Fill(Jet_vector.at(5).Pt(), eventWeight);
      if(Jet_vector.size()>0) eeHist.h_jet_eta_leading->Fill(Jet_vector.at(0).Eta(), eventWeight);
      if(Jet_vector.size()>1) eeHist.h_jet_eta_trailing->Fill(Jet_vector.at(1).Eta(), eventWeight);
      if(Jet_vector.size()>2) eeHist.h_jet_eta_3rd->Fill(Jet_vector.at(2).Eta(), eventWeight);
      if(Jet_vector.size()>3) eeHist.h_jet_eta_4th->Fill(Jet_vector.at(3).Eta(), eventWeight);
      if(Jet_vector.size()>4) eeHist.h_jet_eta_5th->Fill(Jet_vector.at(4).Eta(), eventWeight);
      if(Jet_vector.size()>5) eeHist.h_jet_eta_6th->Fill(Jet_vector.at(5).Eta(), eventWeight);
      if(Jet_vector.size()>0) eeHist.h_jet_phi_leading->Fill(Jet_vector.at(0).Phi(), eventWeight);
      if(Jet_vector.size()>1) eeHist.h_jet_phi_trailing->Fill(Jet_vector.at(1).Phi(), eventWeight);
      if(Jet_vector.size()>2) eeHist.h_jet_phi_3rd->Fill(Jet_vector.at(2).Phi(), eventWeight);
      if(Jet_vector.size()>3) eeHist.h_jet_phi_4th->Fill(Jet_vector.at(3).Phi(), eventWeight);
      if(Jet_vector.size()>4) eeHist.h_jet_phi_5th->Fill(Jet_vector.at(4).Phi(), eventWeight);
      if(Jet_vector.size()>5) eeHist.h_jet_phi_6th->Fill(Jet_vector.at(5).Phi(), eventWeight);
      if(Jet_vector.size()>0) eeHist.h_jet_energy_leading->Fill(Jet_vector.at(0).E(), eventWeight);
      if(Jet_vector.size()>1) eeHist.h_jet_energy_trailing->Fill(Jet_vector.at(1).E(), eventWeight);
      if(Jet_vector.size()>2) eeHist.h_jet_energy_3rd->Fill(Jet_vector.at(2).E(), eventWeight);
      if(Jet_vector.size()>3) eeHist.h_jet_energy_4th->Fill(Jet_vector.at(3).E(), eventWeight);
      if(Jet_vector.size()>4) eeHist.h_jet_energy_5th->Fill(Jet_vector.at(4).E(), eventWeight);
      if(Jet_vector.size()>5) eeHist.h_jet_energy_6th->Fill(Jet_vector.at(5).E(), eventWeight);
      eeHist.h_photon_pt->Fill(Photon_vector.at(0).Pt(), eventWeight);
      eeHist.h_photon_eta->Fill(Photon_vector.at(0).Eta(), eventWeight);
      eeHist.h_photon_phi->Fill(Photon_vector.at(0).Phi(), eventWeight);
      eeHist.h_photon_energy->Fill(Photon_vector.at(0).E(), eventWeight);    
    }//SS mode 
  }//first+second electron "if" statement 

//Non-prompt background computation:
  else if( Electrons.size() > 0.0 and LNTElectrons.size() > 0.0 and eventWeight > 0.0 ) {
  if(type=="Data" and signSelection=="SS" and (Electrons.at(0).Charge*LNTElectrons.at(0).Charge)==1) {
    ee_SS_NT01+=eventWeight;
    el1_transverse.SetMagPhi(Electrons.at(0).LepLV.Pt(), Electrons.at(0).LepLV.Phi());
    el2_transverse.SetMagPhi(LNTElectrons.at(0).LepLV.Pt(), LNTElectrons.at(0).LepLV.Phi());
    eeHistNT01.h_el_pt_leading->Fill(Electrons.at(0).LepLV.Pt(), eventWeight);
    eeHistNT01.h_el_phi_leading->Fill(Electrons.at(0).LepLV.Phi(), eventWeight);
    eeHistNT01.h_el_eta_leading->Fill(Electrons.at(0).LepLV.Eta(), eventWeight);
    eeHistNT01.h_el_energy_leading->Fill(Electrons.at(0).LepLV.E(), eventWeight);
    eeHistNT01.h_el_pt_trailing->Fill(LNTElectrons.at(0).LepLV.Pt(), eventWeight);
    eeHistNT01.h_el_phi_trailing->Fill(LNTElectrons.at(0).LepLV.Phi(), eventWeight);
    eeHistNT01.h_el_eta_trailing->Fill(LNTElectrons.at(0).LepLV.Eta(), eventWeight);
    eeHistNT01.h_el_energy_trailing->Fill(LNTElectrons.at(0).LepLV.E(), eventWeight);
    eeHistNT01.h_InvariantMass->Fill((LNTElectrons.at(0).LepLV+Electrons.at(0).LepLV).M(), eventWeight);
    eeHistNT01.h_InvariantMass_Ph->Fill((LNTElectrons.at(0).LepLV+Electrons.at(0).LepLV+Photon_vector.at(0)).M(), eventWeight);
    eeHistNT01.h_DeltaPhi_met_el1->Fill(met_transverse.DeltaPhi(el1_transverse), eventWeight);
    eeHistNT01.h_DeltaPhi_met_el2->Fill(met_transverse.DeltaPhi(el2_transverse), eventWeight);
    eeHistNT01.h_DeltaPhi_ph_el1->Fill(ph_transverse.DeltaPhi(el1_transverse), eventWeight);
    eeHistNT01.h_DeltaPhi_ph_el2->Fill(ph_transverse.DeltaPhi(el2_transverse), eventWeight);
    eeHistNT01.h_Isolation_el1->Fill(Electrons.at(0).Isolation, eventWeight);
    eeHistNT01.h_Isolation_el2->Fill(LNTElectrons.at(0).Isolation, eventWeight);
    eeHistNT01.h_nVertices->Fill(nVertices, eventWeight);
    eeHistNT01.h_MET->Fill(MET, eventWeight);
    eeHistNT01.h_caloMET->Fill(caloMET, eventWeight);
    eeHistNT01.h_nJets->Fill(Jet_vector.size(), eventWeight);
    eeHistNT01.h_HT->Fill(HT, eventWeight);
    if(Jet_vector.size()>0) eeHistNT01.h_jet_pt_leading->Fill(Jet_vector.at(0).Pt(), eventWeight);
    if(Jet_vector.size()>1) eeHistNT01.h_jet_pt_trailing->Fill(Jet_vector.at(1).Pt(), eventWeight);
    if(Jet_vector.size()>2) eeHistNT01.h_jet_pt_3rd->Fill(Jet_vector.at(2).Pt(), eventWeight);
    if(Jet_vector.size()>3) eeHistNT01.h_jet_pt_4th->Fill(Jet_vector.at(3).Pt(), eventWeight);
    if(Jet_vector.size()>4) eeHistNT01.h_jet_pt_5th->Fill(Jet_vector.at(4).Pt(), eventWeight);
    if(Jet_vector.size()>5) eeHistNT01.h_jet_pt_6th->Fill(Jet_vector.at(5).Pt(), eventWeight);
    if(Jet_vector.size()>0) eeHistNT01.h_jet_eta_leading->Fill(Jet_vector.at(0).Eta(), eventWeight);
    if(Jet_vector.size()>1) eeHistNT01.h_jet_eta_trailing->Fill(Jet_vector.at(1).Eta(), eventWeight);
    if(Jet_vector.size()>2) eeHistNT01.h_jet_eta_3rd->Fill(Jet_vector.at(2).Eta(), eventWeight);
    if(Jet_vector.size()>3) eeHistNT01.h_jet_eta_4th->Fill(Jet_vector.at(3).Eta(), eventWeight);
    if(Jet_vector.size()>4) eeHistNT01.h_jet_eta_5th->Fill(Jet_vector.at(4).Eta(), eventWeight);
    if(Jet_vector.size()>5) eeHistNT01.h_jet_eta_6th->Fill(Jet_vector.at(5).Eta(), eventWeight);
    if(Jet_vector.size()>0) eeHistNT01.h_jet_phi_leading->Fill(Jet_vector.at(0).Phi(), eventWeight);
    if(Jet_vector.size()>1) eeHistNT01.h_jet_phi_trailing->Fill(Jet_vector.at(1).Phi(), eventWeight);
    if(Jet_vector.size()>2) eeHistNT01.h_jet_phi_3rd->Fill(Jet_vector.at(2).Phi(), eventWeight);
    if(Jet_vector.size()>3) eeHistNT01.h_jet_phi_4th->Fill(Jet_vector.at(3).Phi(), eventWeight);
    if(Jet_vector.size()>4) eeHistNT01.h_jet_phi_5th->Fill(Jet_vector.at(4).Phi(), eventWeight);
    if(Jet_vector.size()>5) eeHistNT01.h_jet_phi_6th->Fill(Jet_vector.at(5).Phi(), eventWeight);
    if(Jet_vector.size()>0) eeHistNT01.h_jet_energy_leading->Fill(Jet_vector.at(0).E(), eventWeight);
    if(Jet_vector.size()>1) eeHistNT01.h_jet_energy_trailing->Fill(Jet_vector.at(1).E(), eventWeight);
    if(Jet_vector.size()>2) eeHistNT01.h_jet_energy_3rd->Fill(Jet_vector.at(2).E(), eventWeight);
    if(Jet_vector.size()>3) eeHistNT01.h_jet_energy_4th->Fill(Jet_vector.at(3).E(), eventWeight);
    if(Jet_vector.size()>4) eeHistNT01.h_jet_energy_5th->Fill(Jet_vector.at(4).E(), eventWeight);
    if(Jet_vector.size()>5) eeHistNT01.h_jet_energy_6th->Fill(Jet_vector.at(5).E(), eventWeight);
    eeHistNT01.h_photon_pt->Fill(Photon_vector.at(0).Pt(), eventWeight);
    eeHistNT01.h_photon_eta->Fill(Photon_vector.at(0).Eta(), eventWeight);
    eeHistNT01.h_photon_phi->Fill(Photon_vector.at(0).Phi(), eventWeight);
    eeHistNT01.h_photon_energy->Fill(Photon_vector.at(0).E(), eventWeight); 
    }//SS mode
  }//loose+tight electron

else if( LNTElectrons.size() > 1.0 and eventWeight > 0.0 ) {
 if(type=="Data" and signSelection=="SS" and (LNTElectrons.at(0).Charge*LNTElectrons.at(1).Charge)==1) {
    ee_SS_NT00+=eventWeight;
    el1_transverse.SetMagPhi(LNTElectrons.at(0).LepLV.Pt(), LNTElectrons.at(0).LepLV.Phi());
    el2_transverse.SetMagPhi(LNTElectrons.at(1).LepLV.Pt(), LNTElectrons.at(1).LepLV.Phi());
    eeHistNT00.h_el_pt_leading->Fill(LNTElectrons.at(0).LepLV.Pt(), eventWeight);
    eeHistNT00.h_el_phi_leading->Fill(LNTElectrons.at(0).LepLV.Phi(), eventWeight);
    eeHistNT00.h_el_eta_leading->Fill(LNTElectrons.at(0).LepLV.Eta(), eventWeight);
    eeHistNT00.h_el_energy_leading->Fill(LNTElectrons.at(0).LepLV.E(), eventWeight);
    eeHistNT00.h_el_pt_trailing->Fill(LNTElectrons.at(1).LepLV.Pt(), eventWeight);
    eeHistNT00.h_el_phi_trailing->Fill(LNTElectrons.at(1).LepLV.Phi(), eventWeight);
    eeHistNT00.h_el_eta_trailing->Fill(LNTElectrons.at(1).LepLV.Eta(), eventWeight);
    eeHistNT00.h_el_energy_trailing->Fill(LNTElectrons.at(1).LepLV.E(), eventWeight);
    eeHistNT00.h_InvariantMass->Fill((LNTElectrons.at(0).LepLV+LNTElectrons.at(1).LepLV).M(), eventWeight);
    eeHistNT00.h_InvariantMass_Ph->Fill((LNTElectrons.at(0).LepLV+LNTElectrons.at(1).LepLV+Photon_vector.at(0)).M(), eventWeight);
    eeHistNT00.h_DeltaPhi_met_el1->Fill(met_transverse.DeltaPhi(el1_transverse), eventWeight);
    eeHistNT00.h_DeltaPhi_met_el2->Fill(met_transverse.DeltaPhi(el2_transverse), eventWeight);
    eeHistNT00.h_DeltaPhi_ph_el1->Fill(ph_transverse.DeltaPhi(el1_transverse), eventWeight);
    eeHistNT00.h_DeltaPhi_ph_el2->Fill(ph_transverse.DeltaPhi(el2_transverse), eventWeight);
    eeHistNT00.h_Isolation_el1->Fill(LNTElectrons.at(0).Isolation, eventWeight);
    eeHistNT00.h_Isolation_el2->Fill(LNTElectrons.at(1).Isolation, eventWeight); 
    eeHistNT00.h_nVertices->Fill(nVertices, eventWeight);
    eeHistNT00.h_MET->Fill(MET, eventWeight);
    eeHistNT00.h_caloMET->Fill(caloMET, eventWeight);
    eeHistNT00.h_nJets->Fill(Jet_vector.size(), eventWeight);
    eeHistNT00.h_HT->Fill(HT, eventWeight);
    if(Jet_vector.size()>0) eeHistNT00.h_jet_pt_leading->Fill(Jet_vector.at(0).Pt(), eventWeight);
    if(Jet_vector.size()>1) eeHistNT00.h_jet_pt_trailing->Fill(Jet_vector.at(1).Pt(), eventWeight);
    if(Jet_vector.size()>2) eeHistNT00.h_jet_pt_3rd->Fill(Jet_vector.at(2).Pt(), eventWeight);
    if(Jet_vector.size()>3) eeHistNT00.h_jet_pt_4th->Fill(Jet_vector.at(3).Pt(), eventWeight);
    if(Jet_vector.size()>4) eeHistNT00.h_jet_pt_5th->Fill(Jet_vector.at(4).Pt(), eventWeight);
    if(Jet_vector.size()>5) eeHistNT00.h_jet_pt_6th->Fill(Jet_vector.at(5).Pt(), eventWeight);
    if(Jet_vector.size()>0) eeHistNT00.h_jet_eta_leading->Fill(Jet_vector.at(0).Eta(), eventWeight);
    if(Jet_vector.size()>1) eeHistNT00.h_jet_eta_trailing->Fill(Jet_vector.at(1).Eta(), eventWeight);
    if(Jet_vector.size()>2) eeHistNT00.h_jet_eta_3rd->Fill(Jet_vector.at(2).Eta(), eventWeight);
    if(Jet_vector.size()>3) eeHistNT00.h_jet_eta_4th->Fill(Jet_vector.at(3).Eta(), eventWeight);
    if(Jet_vector.size()>4) eeHistNT00.h_jet_eta_5th->Fill(Jet_vector.at(4).Eta(), eventWeight);
    if(Jet_vector.size()>5) eeHistNT00.h_jet_eta_6th->Fill(Jet_vector.at(5).Eta(), eventWeight);
    if(Jet_vector.size()>0) eeHistNT00.h_jet_phi_leading->Fill(Jet_vector.at(0).Phi(), eventWeight);
    if(Jet_vector.size()>1) eeHistNT00.h_jet_phi_trailing->Fill(Jet_vector.at(1).Phi(), eventWeight);
    if(Jet_vector.size()>2) eeHistNT00.h_jet_phi_3rd->Fill(Jet_vector.at(2).Phi(), eventWeight);
    if(Jet_vector.size()>3) eeHistNT00.h_jet_phi_4th->Fill(Jet_vector.at(3).Phi(), eventWeight);
    if(Jet_vector.size()>4) eeHistNT00.h_jet_phi_5th->Fill(Jet_vector.at(4).Phi(), eventWeight);
    if(Jet_vector.size()>5) eeHistNT00.h_jet_phi_6th->Fill(Jet_vector.at(5).Phi(), eventWeight);
    if(Jet_vector.size()>0) eeHistNT00.h_jet_energy_leading->Fill(Jet_vector.at(0).E(), eventWeight);
    if(Jet_vector.size()>1) eeHistNT00.h_jet_energy_trailing->Fill(Jet_vector.at(1).E(), eventWeight);
    if(Jet_vector.size()>2) eeHistNT00.h_jet_energy_3rd->Fill(Jet_vector.at(2).E(), eventWeight);
    if(Jet_vector.size()>3) eeHistNT00.h_jet_energy_4th->Fill(Jet_vector.at(3).E(), eventWeight);
    if(Jet_vector.size()>4) eeHistNT00.h_jet_energy_5th->Fill(Jet_vector.at(4).E(), eventWeight);
    if(Jet_vector.size()>5) eeHistNT00.h_jet_energy_6th->Fill(Jet_vector.at(5).E(), eventWeight);
    eeHistNT00.h_photon_pt->Fill(Photon_vector.at(0).Pt(), eventWeight);
    eeHistNT00.h_photon_eta->Fill(Photon_vector.at(0).Eta(), eventWeight);
    eeHistNT00.h_photon_phi->Fill(Photon_vector.at(0).Phi(), eventWeight);
    eeHistNT00.h_photon_energy->Fill(Photon_vector.at(0).E(), eventWeight);
    }//SS mode
  }//loose not tight+loose not tight electron

}//event loop closed

  //Cut flow table
  cout << "mumu_SS_NT01 = " << mumu_SS_NT01 << endl; 
  cout << "mumu_SS_NT10 = " << mumu_SS_NT10 << endl;
  cout << "mumu_SS_NT00 = " << mumu_SS_NT00 << endl;
  cout << "mumu_OS = " << mumu_OS << endl;
  cout << "mumu_SS = " << mumu_SS << endl;
  cout << "emu_OS = " << emu_OS << endl;
  cout << "emu_SS = " << emu_SS << endl;
  cout << "emu_SS_NT01 = " << emu_SS_NT01 << endl;
  cout << "emu_SS_NT10 = " << emu_SS_NT10 << endl;
  cout << "emu_SS_NT00 = " << emu_SS_NT00 << endl;
  cout << "ee_OS = " << ee_OS << endl;
  cout << "ee_SS = " << ee_SS << endl;
  cout << "ee_SS_NT01 = " << ee_SS_NT01 << endl;
  cout << "ee_SS_NT10 = " << ee_SS_NT10 << endl;
  cout << "ee_SS_NT00 = " << ee_SS_NT00 << endl;
  cout << "MET_emu_OS = " << MET_emu_OS << endl;
  cout << "MET_emu_SS = " << MET_emu_SS << endl;
  cout << "MET_mumu_OS = " << MET_mumu_OS << endl;
  cout << "MET_mumu_SS = " << MET_mumu_SS << endl;
  cout << "LowHT_mumu_OS = " << LowHT_mumu_OS << endl;
  cout << "LowHT_emu_SS = " << LowHT_emu_SS << endl;
  std::string histfilename=(outfile+".root").c_str();
  TFile *tFile=new TFile(histfilename.c_str(), "RECREATE");
  h_CaloMETTrue_CaloMET->Write();
  h_HLTMET_CaloMET->Write();
  h_Mmumu_MmumuGamma->Write();
  writeHistCollection(mumuHist); 
  writeHistCollection(mumuHistNT01);
  writeHistCollection(mumuHistNT10);
  writeHistCollection(mumuHistNT00);
  writeHistCollection(emuHist);
  writeHistCollection(emuHistNT01);
  writeHistCollection(emuHistNT10);
  writeHistCollection(emuHistNT00);
  writeHistCollection(eeHist);
  writeHistCollection(eeHistNT01);
  writeHistCollection(eeHistNT10);
  writeHistCollection(eeHistNT00);
  tFile->Close(); 
  std::cout<<"Wrote output file "<<histfilename<<std::endl;

  return 0;

}
