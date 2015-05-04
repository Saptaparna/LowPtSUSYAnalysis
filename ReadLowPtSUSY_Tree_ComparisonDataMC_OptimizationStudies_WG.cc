#include "ReadLowPtSUSY_Tree_ComparisonDataMC_OptimizationStudies_WG.h"

int ReadLowPtSUSY_Tree_ComparisonDataMC_OptimizationStudies_WG(std::string infile, std::string outfile, std::string type, std::string signSelection, std::string MCSample){
  
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
  vector<bool>    *el_isLooseWG;
  vector<float>   *el_Dxy;
  vector<float>   *el_Dz;
  Int_t           nElectrons;
  vector<float>   *mu_pt;
  vector<float>   *mu_phi;
  vector<float>   *mu_eta;
  vector<float>   *mu_energy;
  vector<int>     *mu_charge;
  vector<bool>    *mu_isTight;
  vector<bool>    *mu_isLoose;
  vector<bool>    *mu_isLooseWG;
  vector<float>   *mu_Dxy;
  vector<float>   *mu_Dz; 
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
  vector<float>   *jet_btag_csv; 
  Int_t           nJets;
  Float_t         MET;
  Float_t         MET_Phi;
  Float_t         MET_Px;
  Float_t         MET_Py;
  Float_t         MET_Signxx;
  Float_t         MET_Signxy;
  Float_t         MET_Signyx;
  Float_t         MET_Signyy;
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
  vector<bool>    *ph_isMedium;
  vector<bool>    *ph_isLoose;
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
  vector<int>     *mu_Matched_1Mother;
  vector<int>     *mu_Matched_2Mother;
  vector<int>     *mu_Matched_3Mother;
  vector<int>     *mu_Matched_4Mother;
  vector<int>     *mu_Matched_5Mother;
  vector<int>     *mu_Matched_6Mother;
  vector<int>     *mu_Matched_7Mother;
  vector<int>     *mu_Matched_8Mother;
  vector<int>     *mu_Matched_9Mother;
  vector<int>     *mu_Matched_10Mother;
  vector<int>     *el_Matched_1Mother;
  vector<int>     *el_Matched_2Mother;
  vector<int>     *el_Matched_3Mother;
  vector<int>     *el_Matched_4Mother;
  vector<int>     *el_Matched_5Mother;
  vector<int>     *el_Matched_6Mother;
  vector<int>     *el_Matched_7Mother;
  vector<int>     *el_Matched_8Mother;
  vector<int>     *el_Matched_9Mother;
  vector<int>     *el_Matched_10Mother;
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
  vector<float>   *tauPx;
  vector<float>   *tauPy;
  vector<float>   *tauPz;
  vector<float>   *tauE;
  vector<float>   *muPx;
  vector<float>   *muPy;
  vector<float>   *muPz;
  vector<float>   *muE;
  vector<float>   *ePx;
  vector<float>   *ePy;
  vector<float>   *ePz;
  vector<float>   *eE;
  vector<float>   *nuPx;
  vector<float>   *nuPy;
  vector<float>   *nuPz;
  vector<float>   *nuE;


  cout << "Variables defined" << endl;

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
  el_isLooseWG = 0;
  el_Dxy = 0;
  el_Dz = 0;
  mu_pt = 0;
  mu_phi = 0;
  mu_eta = 0;
  mu_energy = 0;
  mu_charge = 0;
  mu_isTight = 0;
  mu_isLoose = 0;
  mu_isLooseWG = 0;
  mu_Dxy = 0;
  mu_Dz = 0;
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
  jet_btag_csv = 0; 
  el_iso = 0;
  mu_iso = 0;
  ph_chIso = 0;
  ph_nuIso = 0;
  ph_phIso = 0;
  ph_isTight = 0;
  ph_isMedium = 0;
  ph_isLoose = 0;
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
  if(type=="MC" and not (MCSample=="Signal"))
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
    if(MCSample=="ZGToLLG") mu_Matched_1Mother = 0;
    if(MCSample=="ZGToLLG") mu_Matched_2Mother = 0;
    if(MCSample=="ZGToLLG") mu_Matched_3Mother = 0;
    if(MCSample=="ZGToLLG") mu_Matched_4Mother = 0;
    if(MCSample=="ZGToLLG") mu_Matched_5Mother = 0;
    if(MCSample=="ZGToLLG") mu_Matched_6Mother = 0;
    if(MCSample=="ZGToLLG") mu_Matched_7Mother = 0;
    if(MCSample=="ZGToLLG") mu_Matched_8Mother = 0;
    if(MCSample=="ZGToLLG") mu_Matched_9Mother = 0;
    if(MCSample=="ZGToLLG") mu_Matched_10Mother = 0;
    if(MCSample=="ZGToLLG") el_Matched_1Mother = 0;
    if(MCSample=="ZGToLLG") el_Matched_2Mother = 0;
    if(MCSample=="ZGToLLG") el_Matched_3Mother = 0;
    if(MCSample=="ZGToLLG") el_Matched_4Mother = 0;
    if(MCSample=="ZGToLLG") el_Matched_5Mother = 0;
    if(MCSample=="ZGToLLG") el_Matched_6Mother = 0;
    if(MCSample=="ZGToLLG") el_Matched_7Mother = 0;
    if(MCSample=="ZGToLLG") el_Matched_8Mother = 0;
    if(MCSample=="ZGToLLG") el_Matched_9Mother = 0;
    if(MCSample=="ZGToLLG") el_Matched_10Mother = 0;
    if(MCSample=="ZGToLLG") tauPx = 0;
    if(MCSample=="ZGToLLG") tauPy = 0;
    if(MCSample=="ZGToLLG") tauPz = 0;
    if(MCSample=="ZGToLLG") tauE = 0;
    if(MCSample=="ZGToLLG") muPx = 0;
    if(MCSample=="ZGToLLG") muPy = 0;
    if(MCSample=="ZGToLLG") muPz = 0;
    if(MCSample=="ZGToLLG") muE = 0;
    if(MCSample=="ZGToLLG") ePx = 0;
    if(MCSample=="ZGToLLG") ePy = 0;
    if(MCSample=="ZGToLLG") ePz = 0;
    if(MCSample=="ZGToLLG") eE = 0;
    if(MCSample=="ZGToLLG") nuPx = 0;
    if(MCSample=="ZGToLLG") nuPy = 0;
    if(MCSample=="ZGToLLG") nuPz = 0;
    if(MCSample=="ZGToLLG") nuE = 0;
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
  tree->SetBranchAddress("ph_isMedium", &(ph_isMedium)); 
  tree->SetBranchAddress("ph_isTight", &(ph_isTight));
  tree->SetBranchAddress("ph_isLoose", &(ph_isLoose));
  tree->SetBranchAddress("nPhotons", &(nPhotons));
  tree->SetBranchAddress("el_pt", &(el_pt));
  tree->SetBranchAddress("el_eta", &(el_eta));
  tree->SetBranchAddress("el_phi", &(el_phi));
  tree->SetBranchAddress("el_energy", &(el_energy));
  tree->SetBranchAddress("el_charge", &(el_charge));
  tree->SetBranchAddress("el_isTight", &(el_isTight)); 
  tree->SetBranchAddress("el_isLoose", &(el_isLoose)); 
  tree->SetBranchAddress("el_isLooseWG", &(el_isLooseWG));
  tree->SetBranchAddress("el_Dxy", &(el_Dxy));
  tree->SetBranchAddress("el_Dz", &(el_Dz));
  tree->SetBranchAddress("el_iso", &(el_iso));
  tree->SetBranchAddress("nElectrons", &(nElectrons));
  tree->SetBranchAddress("mu_pt", &(mu_pt));
  tree->SetBranchAddress("mu_eta", &(mu_eta));
  tree->SetBranchAddress("mu_phi", &(mu_phi));
  tree->SetBranchAddress("mu_energy", &(mu_energy)); 
  tree->SetBranchAddress("mu_charge", &(mu_charge)); 
  tree->SetBranchAddress("mu_isTight", &(mu_isTight));
  tree->SetBranchAddress("mu_isLoose", &(mu_isLoose)); 
  tree->SetBranchAddress("mu_isLooseWG", &(mu_isLooseWG));
  tree->SetBranchAddress("mu_Dxy", &(mu_Dxy));
  tree->SetBranchAddress("mu_Dz", &(mu_Dz));
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
  tree->SetBranchAddress("jet_btag_csv", &(jet_btag_csv));
  tree->SetBranchAddress("nJets", &(nJets));
  tree->SetBranchAddress("MET", &(MET));
  tree->SetBranchAddress("MET_Phi", &(MET_Phi));
  tree->SetBranchAddress("MET_Px", &(MET_Px));
  tree->SetBranchAddress("MET_Py", &(MET_Py));
  tree->SetBranchAddress("MET_Signxx", &(MET_Signxx));
  tree->SetBranchAddress("MET_Signxy", &(MET_Signxy));
  tree->SetBranchAddress("MET_Signyx", &(MET_Signyx));
  tree->SetBranchAddress("MET_Signyy", &(MET_Signyy));
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
  
  if(type=="MC" and not (MCSample=="Signal")){
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
    if(MCSample=="ZGToLLG")tree->SetBranchAddress("mu_Matched_1Mother", &(mu_Matched_1Mother));
    if(MCSample=="ZGToLLG")tree->SetBranchAddress("mu_Matched_2Mother", &(mu_Matched_2Mother));
    if(MCSample=="ZGToLLG")tree->SetBranchAddress("mu_Matched_3Mother", &(mu_Matched_3Mother));
    if(MCSample=="ZGToLLG")tree->SetBranchAddress("mu_Matched_4Mother", &(mu_Matched_4Mother));
    if(MCSample=="ZGToLLG")tree->SetBranchAddress("mu_Matched_5Mother", &(mu_Matched_5Mother));
    if(MCSample=="ZGToLLG")tree->SetBranchAddress("mu_Matched_6Mother", &(mu_Matched_6Mother));
    if(MCSample=="ZGToLLG")tree->SetBranchAddress("mu_Matched_7Mother", &(mu_Matched_7Mother)); 
    if(MCSample=="ZGToLLG")tree->SetBranchAddress("mu_Matched_8Mother", &(mu_Matched_8Mother));
    if(MCSample=="ZGToLLG")tree->SetBranchAddress("mu_Matched_9Mother", &(mu_Matched_9Mother));
    if(MCSample=="ZGToLLG")tree->SetBranchAddress("mu_Matched_10Mother", &(mu_Matched_10Mother));
    if(MCSample=="ZGToLLG")tree->SetBranchAddress("el_Matched_1Mother", &(el_Matched_1Mother));
    if(MCSample=="ZGToLLG")tree->SetBranchAddress("el_Matched_2Mother", &(el_Matched_2Mother));
    if(MCSample=="ZGToLLG")tree->SetBranchAddress("el_Matched_3Mother", &(el_Matched_3Mother));
    if(MCSample=="ZGToLLG")tree->SetBranchAddress("el_Matched_4Mother", &(el_Matched_4Mother));
    if(MCSample=="ZGToLLG")tree->SetBranchAddress("el_Matched_5Mother", &(el_Matched_5Mother));
    if(MCSample=="ZGToLLG")tree->SetBranchAddress("el_Matched_6Mother", &(el_Matched_6Mother));
    if(MCSample=="ZGToLLG")tree->SetBranchAddress("el_Matched_7Mother", &(el_Matched_7Mother));
    if(MCSample=="ZGToLLG")tree->SetBranchAddress("el_Matched_8Mother", &(el_Matched_8Mother));
    if(MCSample=="ZGToLLG")tree->SetBranchAddress("el_Matched_9Mother", &(el_Matched_9Mother));
    if(MCSample=="ZGToLLG")tree->SetBranchAddress("el_Matched_10Mother", &(el_Matched_10Mother));
    if(MCSample=="ZGToLLG")tree->SetBranchAddress("tauPx", &(tauPx));
    if(MCSample=="ZGToLLG")tree->SetBranchAddress("tauPy", &(tauPy));
    if(MCSample=="ZGToLLG")tree->SetBranchAddress("tauPz", &(tauPz));
    if(MCSample=="ZGToLLG")tree->SetBranchAddress("tauE", &(tauE)); 
    if(MCSample=="ZGToLLG")tree->SetBranchAddress("muPx", &(muPx));
    if(MCSample=="ZGToLLG")tree->SetBranchAddress("muPy", &(muPy));
    if(MCSample=="ZGToLLG")tree->SetBranchAddress("muPz", &(muPz));
    if(MCSample=="ZGToLLG")tree->SetBranchAddress("muE", &(muE));  
    if(MCSample=="ZGToLLG")tree->SetBranchAddress("ePx", &(ePx));
    if(MCSample=="ZGToLLG")tree->SetBranchAddress("ePy", &(ePy));
    if(MCSample=="ZGToLLG")tree->SetBranchAddress("ePz", &(ePz));
    if(MCSample=="ZGToLLG")tree->SetBranchAddress("eE", &(eE));  
    if(MCSample=="ZGToLLG")tree->SetBranchAddress("nuPx", &(nuPx));
    if(MCSample=="ZGToLLG")tree->SetBranchAddress("nuPy", &(nuPy));
    if(MCSample=="ZGToLLG")tree->SetBranchAddress("nuPz", &(nuPz));
    if(MCSample=="ZGToLLG")tree->SetBranchAddress("nuE", &(nuE));  
    tree->SetBranchAddress("Wt", &(Wt));
    tree->SetBranchAddress("Ztt", &(Ztt));
    tree->SetBranchAddress("Zee", &(Zee));
    tree->SetBranchAddress("Zmumu", &(Zmumu));
    tree->SetBranchAddress("Znunu", &(Znunu));
    tree->SetBranchAddress("isPhoton", &(isPhoton));
  }
  //Booking histograms:
  cout << "Branch addresses set" << endl;


  double emu_OS = 0;
  double emu_SS = 0;
  double mumu_OS = 0;
  double mumu_SS = 0;
  double ee_OS = 0;
  double ee_SS = 0;
  double mumu_OS_Z = 0;

  cout << "Histograms defined" << endl;

  TFile *trigger1D1=new TFile("CaloMET_ForTalk.root");
  //TFile *trigger1D1=new TFile("TurnOn_SinglePhotonParked_Run2012D_22Jan2013_All.root");
  //TFile *trigger1D2=new TFile("TurnOn_METParked_Run2012D_All.root");//evaluated at the analysis photon ID point
  //TFile *trigger1D2=new TFile("Test_METParked_Run2012D_All_OLD.root");
  //TFile *trigger1D2=new TFile("Test_TriggerTurn_METParked_MWP.root");
  //TFile *trigger1D2=new TFile("Test_LWP_SpikeCleaning.root");
  //TFile *trigger1D2=new TFile("TEST_TEST_LWP.root");
  //TFile *trigger1D2=new TFile("Test_TriggerLoose_WithoutBeamHalo.root");
  //TFile *trigger1D2=new TFile("TEST_MWP.root");
  //TFile *trigger1D2=new TFile("PhotonPt_TightWP_ForTalk.root");
  TFile *trigger1D2=new TFile("Trigger_Output_LowPtSUSY_Tree_METParked_Run2012D_All.root");

  TF1* fit_curve1=(TF1*)trigger1D1->Get("fit_caloMET");
  TF1* fit_curve2=(TF1*)trigger1D2->Get("fit_PHO");

  TFile* fileData = TFile::Open("/uscms_data/d2/lpcljm/sapta/SUSYSearch/CMSSW_5_3_16_patch1/src/MPAnalyzer/PU_Data/PU_Data_All.root");
  TH1F *h_nVertices_Data = (TH1F*) fileData->Get("h_nVertices_Data");
  h_nVertices_Data->Scale(1.0/h_nVertices_Data->Integral()); 
  TFile* fileMC;
  TH1F *h_nWeights;
  if (type=="MC"){ 
    if(MCSample=="Stop100") fileMC = TFile::Open("/uscms_data/d2/lpcljm/sapta/SUSYSearch/CMSSW_5_3_16_patch1/src/MPAnalyzer/PU_MC_Signal/Output_LowPtSUSY_Tree_Stop100_Chargino75_Neutralino30_Ntuple_Pta20.root");
    if(MCSample=="ZGToLLG") fileMC = TFile::Open("/uscms_data/d2/lpcljm/sapta/SUSYSearch/CMSSW_5_3_16_patch1/src/MPAnalyzer/PU_MC/Output_LowPtSUSY_Tree_ZGToLLG_8TeV_CaloMET_All.root");
    if(MCSample=="TTG") fileMC = TFile::Open("/uscms_data/d2/lpcljm/sapta/SUSYSearch/CMSSW_5_3_16_patch1/src/MPAnalyzer/PU_MC/Output_LowPtSUSY_Tree_TTGJets_8TeV_madgraph_CaloMET_All.root");
    if(MCSample=="ZZJetsTo4L") fileMC = TFile::Open("/uscms_data/d2/lpcljm/sapta/SUSYSearch/CMSSW_5_3_16_patch1/src/MPAnalyzer/PU_MC/Output_LowPtSUSY_Tree_ZZJetsTo4L_TuneZ2star_CaloMET_All.root");
    if(MCSample=="Tbar_tW") fileMC = TFile::Open("/uscms_data/d2/lpcljm/sapta/SUSYSearch/CMSSW_5_3_16_patch1/src/MPAnalyzer/PU_MC/Output_LowPtSUSY_Tree_Tbar_tW_CaloMET_All.root");
    if(MCSample=="T_tW")  fileMC = TFile::Open("/uscms_data/d2/lpcljm/sapta/SUSYSearch/CMSSW_5_3_16_patch1/src/MPAnalyzer/PU_MC/Output_LowPtSUSY_Tree_T_tW_CaloMET_All.root");
    if(MCSample=="WGstarToLNu2Mu") fileMC = TFile::Open("/uscms_data/d2/lpcljm/sapta/SUSYSearch/CMSSW_5_3_16_patch1/src/MPAnalyzer/PU_MC/Output_LowPtSUSY_Tree_WGstarToLNu2Mu_TuneZ2star_CaloMET_All.root");
    if(MCSample=="WmWmqq") fileMC = TFile::Open("/uscms_data/d2/lpcljm/sapta/SUSYSearch/CMSSW_5_3_16_patch1/src/MPAnalyzer/PU_MC/Output_LowPtSUSY_Tree_WmWmqq_8TeV_madgraph_CaloMET_All.root");
    if(MCSample=="WpWpqq") fileMC = TFile::Open("/uscms_data/d2/lpcljm/sapta/SUSYSearch/CMSSW_5_3_16_patch1/src/MPAnalyzer/PU_MC/Output_LowPtSUSY_Tree_WpWpqq_8TeV_madgraph_CaloMET_All.root");
    if(MCSample=="WZJetsTo3LNu") fileMC = TFile::Open("/uscms_data/d2/lpcljm/sapta/SUSYSearch/CMSSW_5_3_16_patch1/src/MPAnalyzer/PU_MC/Output_LowPtSUSY_Tree_WZJetsTo3LNu_CaloMET_All.root");
    if(MCSample=="WWGJets") fileMC = TFile::Open("/uscms_data/d2/lpcljm/sapta/SUSYSearch/CMSSW_5_3_16_patch1/src/MPAnalyzer/PU_MC/Output_LowPtSUSY_Tree_WWGJets_8TeV_madgraph_CaloMET_All.root");
    if(MCSample=="JPsiToMuMu") fileMC = TFile::Open("/uscms_data/d2/lpcljm/sapta/SUSYSearch/CMSSW_5_3_16_patch1/src/MPAnalyzer/PU_MC/Output_LowPtSUSY_Tree_JPsiToMuMu_All.root");
    TH1F *h_nVertices_MC = (TH1F*) fileMC->Get("h_nVertices_MC");
    h_nVertices_MC->Scale(1.0/h_nVertices_MC->Integral());
    h_nWeights = (TH1F*) h_nVertices_Data->Clone("h_nWeights");
    h_nWeights->Divide(h_nVertices_MC);
  }

  TH1F *h_HT_LowCutValue_mumu = new TH1F("h_HT_LowCutValue_mumu", "Scatter Plot of HT (MuMu Channel); HT [GeV]", 34, 30.0, 710.0);  h_HT_LowCutValue_mumu->Sumw2();
  TH1F *h_HTb_LowCutValue_mumu = new TH1F("h_HTb_LowCutValue_mumu", "Scatter Plot of HTb (MuMu Channel); HTb [GeV]", 34, 30.0, 710.0);  h_HTb_LowCutValue_mumu->Sumw2();

  TH1F *h_HT_LowCutValue_emu = new TH1F("h_HT_LowCutValue_emu", "Scatter Plot of HT (MuMu Channel); HT [GeV]", 34, 30.0, 710.0);  h_HT_LowCutValue_emu->Sumw2();
  TH1F *h_HTb_LowCutValue_emu = new TH1F("h_HTb_LowCutValue_emu", "Scatter Plot of HTb (MuMu Channel); HTb [GeV]", 34, 30.0, 710.0);  h_HTb_LowCutValue_emu->Sumw2();

  TH2F *h_HT_Mmumu_LowCutValue = new TH2F("h_HT_Mmumu_LowCutValue", "Scatter Plot of HT versus M_{#mu#mu}; HT [GeV]; M_{#mu#mu} [GeV]", 34, 30.0, 710.0, 29, 0.0, 290.0);  h_HT_Mmumu_LowCutValue->Sumw2();
  TH2F *h_HT_Mmumu_HighCutValue = new TH2F("h_HT_Mmumu_HighCutValue", "Scatter Plot of HT versus M_{#mu#mu}; HT [GeV]; M_{#mu#mu} [GeV]", 34, 30.0, 710.0, 29, 0.0, 290.0);  h_HT_Mmumu_HighCutValue->Sumw2();
  TH2F *h_HT_MmumuGamma_LowCutValue = new TH2F("h_HT_MmumuGamma_LowCutValue", "Scatter Plot of HT versus M_{#mu#mu#gamma}; HT [GeV]; M_{#mu#mu#gamma} [GeV]", 34, 30.0, 710.0, 14, 0.0, 280.0);  h_HT_MmumuGamma_LowCutValue->Sumw2();

  TH2F *h_HT_Memu_LowCutValue = new TH2F("h_HT_Memu_LowCutValue", "Scatter Plot of HT versus M_{e#mu}; HT [GeV]; M_{e#mu} [GeV]", 34, 30.0, 710.0, 29, 0.0, 290.0);  h_HT_Memu_LowCutValue->Sumw2();

  TH2F *h_HT_Memu_HighCutValue = new TH2F("h_HT_Memu_HighCutValue", "Scatter Plot of HT versus M_{e#mu}; HT [GeV]; M_{e#mu} [GeV]", 34, 30.0, 710.0, 29, 0.0, 290.0);  h_HT_Memu_HighCutValue->Sumw2();
 
  TH2F *h_HT_MemuGamma_LowCutValue = new TH2F("h_HT_MemuGamma_LowCutValue", "Scatter Plot of HT versus M_{e#mu#gamma}; HT [GeV]; M_{e#mu#gamma} [GeV]", 34, 30.0, 710.0, 14, 0.0, 280.0);  h_HT_MemuGamma_LowCutValue->Sumw2();

  TH2F *h_HTb_Mmumu_LowCutValue = new TH2F("h_HTb_Mmumu_LowCutValue", "Scatter Plot of HTb versus M_{#mu#mu}; HTb [GeV]; M_{#mu#mu} [GeV]", 34, 30.0, 710.0, 29, 0.0, 290.0);  h_HTb_Mmumu_LowCutValue->Sumw2();

  TH2F *h_HTb_Mmumu_HighCutValue = new TH2F("h_HTb_Mmumu_HighCutValue", "Scatter Plot of HTb versus M_{#mu#mu}; HTb [GeV]; M_{#mu#mu} [GeV]", 34, 30.0, 710.0, 29, 0.0, 290.0);  h_HTb_Mmumu_HighCutValue->Sumw2();
  TH2F *h_HTb_MmumuGamma_LowCutValue = new TH2F("h_HTb_MmumuGamma_LowCutValue", "Scatter Plot of HTb versus M_{#mu#mu#gamma}; HTb [GeV]; M_{#mu#mu#gamma} [GeV]", 34, 30.0, 710.0, 14, 0.0, 280.0);  h_HTb_MmumuGamma_LowCutValue->Sumw2();

  TH2F *h_HTb_Memu_LowCutValue = new TH2F("h_HTb_Memu_LowCutValue", "Scatter Plot of HTb versus M_{e#mu}; HTb [GeV]; M_{e#mu} [GeV]", 34, 30.0, 710.0, 29, 0.0, 290.0);  h_HTb_Memu_LowCutValue->Sumw2();
  TH2F *h_HTb_Memu_HighCutValue = new TH2F("h_HTb_Memu_HighCutValue", "Scatter Plot of HTb versus M_{e#mu}; HTb [GeV]; M_{e#mu} [GeV]", 34, 30.0, 710.0, 29, 0.0, 290.0);  h_HTb_Memu_HighCutValue->Sumw2();

  TH2F *h_HTb_MemuGamma_LowCutValue = new TH2F("h_HTb_MemuGamma_LowCutValue", "Scatter Plot of HTb versus M_{e#mu#gamma}; HTb [GeV]; M_{e#mu#gamma} [GeV]", 34, 30.0, 710.0, 14, 0.0, 280.0);  h_HTb_MemuGamma_LowCutValue->Sumw2();

  TH2F *h_HT_SVFitmumu_LowCutValue = new TH2F("h_HT_SVFitmumu_LowCutValue", "Scatter Plot of HT versus M_{#mu#mu}; HT [GeV]; M_{#mu#mu} [GeV]", 34, 30.0, 710.0, 24, 0.0, 480.0);  h_HT_SVFitmumu_LowCutValue->Sumw2();
  TH2F *h_HT_SVFitmumuGamma_LowCutValue = new TH2F("h_HT_SVFitmumuGamma_LowCutValue", "Scatter Plot of HT versus M_{#mu#mu#gamma}; HT [GeV]; M_{#mu#mu#gamma} [GeV]",  34, 30.0, 710.0, 24, 0.0, 480.0);  h_HT_SVFitmumuGamma_LowCutValue->Sumw2();

  TH2F *h_HT_SVFitemu_LowCutValue = new TH2F("h_HT_SVFitemu_LowCutValue", "Scatter Plot of HT versus M_{e#mu}; HT [GeV]; M_{e#mu} [GeV]", 34, 30.0, 710.0, 24, 0.0, 480.0);  h_HT_SVFitemu_LowCutValue->Sumw2();
  TH2F *h_HT_SVFitemuGamma_LowCutValue = new TH2F("h_HT_SVFitemuGamma_LowCutValue", "Scatter Plot of HT versus M_{e#mu#gamma}; HT [GeV]; M_{e#mu#gamma} [GeV]", 34, 30.0, 710.0, 24, 0.0, 480.0);  h_HT_SVFitemuGamma_LowCutValue->Sumw2();

  TH2F *h_HTb_SVFitmumu_LowCutValue = new TH2F("h_HTb_SVFitmumu_LowCutValue", "Scatter Plot of HTb versus M_{#mu#mu}; HTb [GeV]; M_{#mu#mu} [GeV]", 34, 30.0, 710.0, 24, 0.0, 480.0);  h_HTb_SVFitmumu_LowCutValue->Sumw2();
  TH2F *h_HTb_SVFitmumuGamma_LowCutValue = new TH2F("h_HTb_SVFitmumuGamma_LowCutValue", "Scatter Plot of HTb versus M_{#mu#mu#gamma}; HTb [GeV]; M_{#mu#mu#gamma} [GeV]", 34, 30.0, 710.0, 24, 0.0, 480.0);  h_HTb_SVFitmumuGamma_LowCutValue->Sumw2();

  TH2F *h_HTb_SVFitemu_LowCutValue = new TH2F("h_HTb_SVFitemu_LowCutValue", "Scatter Plot of HTb versus M_{e#mu}; HTb [GeV]; M_{e#mu} [GeV]", 34, 30.0, 710.0, 24, 0.0, 480.0);  h_HTb_SVFitemu_LowCutValue->Sumw2();
  TH2F *h_HTb_SVFitemuGamma_LowCutValue = new TH2F("h_HTb_SVFitemuGamma_LowCutValue", "Scatter Plot of HTb versus M_{e#mu#gamma}; HTb [GeV]; M_{e#mu#gamma} [GeV]", 34, 30.0, 710.0, 24, 0.0, 480.0);  h_HTb_SVFitemuGamma_LowCutValue->Sumw2();


  double vetoCounter = 0;
  double HT_Value_MuMu[35], InvMass_Value_MuMu[35], MET_Value_MuMu[35], HTb_Value_MuMu[35];
  double HT_Value_ElMu[35], InvMass_Value_ElMu[35], MET_Value_ElMu[35], HTb_Value_ElMu[35];

  double HT_Yield_MuMu[35], InvMass_Yield_MuMu[35], HTb_Yield_MuMu[35];
  double HT_Yield_ElMu[35], InvMass_Yield_ElMu[35], HTb_Yield_ElMu[35];
  
  double SVFitMass_Value_ElMu[35];
  double SVFitMass_Value_MuMu[35];

  double SVFitMassGamma_Value_ElMu[35];
  double SVFitMassGamma_Value_MuMu[35];

  double InvMassGamma_Value_ElMu[35];
  double InvMassGamma_Value_MuMu[35];

  double HT_InvMass_Yield_MuMu[35][35];
  double HT_InvMass_Yield_ElMu[35][35];

  double HT_InvMassGamma_Yield_MuMu[35][35];
  double HT_InvMassGamma_Yield_ElMu[35][35];

  double HTb_InvMassGamma_Yield_MuMu[35][35];
  double HTb_InvMassGamma_Yield_ElMu[35][35];

  double HT_InvMass_MET_Yield_MuMu[35][35][35];
  double HT_InvMass_MET_Yield_ElMu[35][35][35];

  double HTb_InvMass_Yield_MuMu[35][35];
  double HTb_InvMass_Yield_ElMu[35][35];

  double HTb_InvMass_MET_Yield_MuMu[35][35][35];
  double HTb_InvMass_MET_Yield_ElMu[35][35][35];

  double HT_svfitMass_Yield_ElMu[35][35];
  double HT_svfitMass_Yield_MuMu[35][35];

  double HTb_svfitMass_Yield_ElMu[35][35];
  double HTb_svfitMass_Yield_MuMu[35][35];

  double HT_svfitMassGamma_Yield_ElMu[35][35];
  double HT_svfitMassGamma_Yield_MuMu[35][35];
  
  double HTb_svfitMassGamma_Yield_ElMu[35][35];
  double HTb_svfitMassGamma_Yield_MuMu[35][35];

  for(unsigned int i=0; i<35; ++i){
    HT_Value_MuMu[i] = 0.0;
    InvMass_Value_MuMu[i] = 0.0;
    MET_Value_MuMu[i] = 0.0;
    HTb_Value_MuMu[i] = 0.0;
    HT_Yield_MuMu[i] = 0.0;
    HTb_Yield_MuMu[i] = 0.0;
    InvMass_Yield_MuMu[i] = 0.0; 
   for(unsigned int j=0; j<35; ++j){
      HT_InvMass_Yield_MuMu[i][j] = 0.0;
      for(unsigned int k=0; k<35; ++k){
        HT_InvMass_MET_Yield_MuMu[i][j][k] = 0.0;
      }
    }

    for(unsigned int s=0; s<35; ++s){
      HTb_InvMass_Yield_MuMu[i][s] = 0.0;
      for(unsigned int t=0; t<35; ++t){
        HTb_InvMass_MET_Yield_MuMu[i][s][t] = 0.0;
      }
    }
  }

  for(unsigned int i=0; i<35; ++i){
    HT_Value_ElMu[i] = 0.0;
    InvMass_Value_ElMu[i] = 0.0;
    MET_Value_ElMu[i] = 0.0;
    HTb_Value_ElMu[i] = 0.0;
    HT_Yield_ElMu[i] = 0.0;
    HTb_Yield_ElMu[i] = 0.0;
    InvMass_Yield_ElMu[i] = 0.0;
    for(unsigned int j=0; j<35; ++j){
      HT_InvMass_Yield_ElMu[i][j] = 0.0;
      for(unsigned int k=0; k<35; ++k){
        HT_InvMass_MET_Yield_ElMu[i][j][k] = 0.0;
      }
    }

    for(unsigned int s=0; s<35; ++s){
      HTb_InvMass_Yield_ElMu[i][s] = 0.0;
      for(unsigned int t=0; t<35; ++t){
        HTb_InvMass_MET_Yield_ElMu[i][s][t] = 0.0;
      }
    }
  }

  int nEvents=tree->GetEntries();
  int nEvents_Incoming  = 0.0;
  int nEvents_Trigger = 0.0;
  int nEvents_Dimuon = 0.0;
  int nEvents_ElMu = 0.0;
  int nEvents_ElEl = 0.0;

  std::cout << "nEvents= " << nEvents << std::endl;
  for (int i=0; i<nEvents; ++i)
    {
     tree->GetEvent(i);
     if(type=="Data")
       {
       if(fired_HLTPhoIdMet!=1) continue; 
       }
     nEvents_Incoming++; 

/*
     if(type=="MC")
       {
       if(Ztt==false) continue;
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
       photon.isMedium=ph_isMedium->at(j);
       photon.isLoose=ph_isLoose->at(j);
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
       electron.isLooseWG=el_isLooseWG->at(j);
       electron.isolation=el_iso->at(j);
       electron.dxy=el_Dxy->at(j);
       electron.dz=el_Dz->at(j);
       if(MCSample!="Signal"){
         if(type=="MC" and MCSample=="ZGToLLG") electron.matched = el_Matched->at(j);
         if(type=="MC" and MCSample=="ZGToLLG") electron.mother1 = el_Matched_1Mother->at(j);
         if(type=="MC" and MCSample=="ZGToLLG") electron.mother2 = el_Matched_2Mother->at(j);
         if(type=="MC" and MCSample=="ZGToLLG") electron.mother3 = el_Matched_3Mother->at(j);
         if(type=="MC" and MCSample=="ZGToLLG") electron.mother4 = el_Matched_4Mother->at(j);
         if(type=="MC" and MCSample=="ZGToLLG") electron.mother5 = el_Matched_5Mother->at(j);
         if(type=="MC" and MCSample=="ZGToLLG") electron.mother6 = el_Matched_6Mother->at(j);
         if(type=="MC" and MCSample=="ZGToLLG") electron.mother7 = el_Matched_7Mother->at(j);
         if(type=="MC" and MCSample=="ZGToLLG") electron.mother8 = el_Matched_8Mother->at(j);
         if(type=="MC" and MCSample=="ZGToLLG") electron.mother9 = el_Matched_9Mother->at(j);
         if(type=="MC" and MCSample=="ZGToLLG") electron.mother10 = el_Matched_10Mother->at(j);
       }
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
       muon.isLooseWG=mu_isLooseWG->at(j);
       muon.isolation=mu_iso->at(j);
       muon.dxy=mu_Dxy->at(j);
       muon.dz=mu_Dz->at(j);
       if(MCSample!="Signal"){
         if(type=="MC" and MCSample=="ZGToLLG") muon.matched = mu_Matched->at(j);
         if(type=="MC" and MCSample=="ZGToLLG") muon.mother1 = mu_Matched_1Mother->at(j);
         if(type=="MC" and MCSample=="ZGToLLG") muon.mother2 = mu_Matched_2Mother->at(j);
         if(type=="MC" and MCSample=="ZGToLLG") muon.mother3 = mu_Matched_3Mother->at(j);
         if(type=="MC" and MCSample=="ZGToLLG") muon.mother4 = mu_Matched_4Mother->at(j);
         if(type=="MC" and MCSample=="ZGToLLG") muon.mother5 = mu_Matched_5Mother->at(j);
         if(type=="MC" and MCSample=="ZGToLLG") muon.mother6 = mu_Matched_6Mother->at(j);
         if(type=="MC" and MCSample=="ZGToLLG") muon.mother7 = mu_Matched_7Mother->at(j); 
         if(type=="MC" and MCSample=="ZGToLLG") muon.mother8 = mu_Matched_8Mother->at(j);
         if(type=="MC" and MCSample=="ZGToLLG") muon.mother9 = mu_Matched_9Mother->at(j);
         if(type=="MC" and MCSample=="ZGToLLG") muon.mother10 = mu_Matched_10Mother->at(j);
       } 
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
     if(returnAnaPhoton(photons.at(k))==true and swissCross.at(k)<0.90)
     //if(fabs(photons.at(k).eta)<1.444 and photons.at(k).pT>30.0 and photons.at(k).isLoose==1 and photons.at(k).phIsoLoose==1 and photons.at(k).phpixelVeto==0 and photons.at(k).phSigmaIetaIeta > 0.001 and swissCross.at(k)<0.90 and std::sqrt(photons.at(k).phSigmaIphiIphi)>0.009 and photons.at(k).phR9 < 1.0)
     //if(fabs(photons.at(k).eta)<1.444 and photons.at(k).pT>30.0 and photons.at(k).isLoose==1 and photons.at(k).phIsoLoose==1 and photons.at(k).phpixelVeto==0)// and photons.at(k).phSigmaIetaIeta > 0.001 and photons.at(k).phSigmaIphiIphi>0.0001 and swissCross.at(k)<0.90)// and photons.at(k).phSigmaIphiIphi>0.0001)// and swissCross.at(k)<0.90)
    {
       Photon.SetPtEtaPhiE(photons.at(k).pT, photons.at(k).eta, photons.at(k).phi, photons.at(k).energy);
       bool isGoodPhoton=true;
       for(unsigned int j=0; j<electrons.size(); ++j)
         {
         TLorentzVector Electron;
         if(electrons.at(j).pT > 7.0 and electrons.at(j).isLooseWG==1 and electrons.at(j).isolation < 0.50 and electrons.at(j).isolation*electrons.at(j).pT < 5.0 and fabs(electrons.at(j).dxy)<0.01 and fabs(electrons.at(j).dz)<0.1)
           {
           Electron.SetPtEtaPhiE(electrons.at(j).pT, electrons.at(j).eta, electrons.at(j).phi, electrons.at(j).energy);
           double DRph_el = Photon.DeltaR(Electron);
           if(DRph_el<0.5) isGoodPhoton=false;
           }
         }

       for(unsigned int j=0; j<muons.size(); ++j)
         {
         TLorentzVector Muon;
         if(muons.at(j).pT > 3.0 and muons.at(j).isLoose==1 and muons.at(j).isLooseWG==1 and muons.at(j).isolation < 0.50 and muons.at(j).isolation*muons.at(j).pT < 5.0 and fabs(muons.at(j).dxy)<0.01 and fabs(muons.at(j).dz)<0.1){
         Muon.SetPtEtaPhiE(muons.at(j).pT, muons.at(j).eta, muons.at(j).phi, muons.at(j).energy);
         double DRph_mu = Photon.DeltaR(Muon);
         if(DRph_mu<0.5) isGoodPhoton=false;
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
     jet.BTag_csv = jet_btag_csv->at(j);
     jets.push_back(jet);
   }
   // Now sorting this vector of structs
   std::sort (jets.begin(), jets.end(), sortJetsInDescendingpT);
   double HT = 0.0; //jets are sorted. Don't care as far as HT is concerned.
   std::vector<AnalysisJetInfo> Jets;
   Jets.clear();
   for(unsigned int k=0; k<jets.size(); ++k)
   {
     AnalysisJetInfo Jet;
     if(fabs(jets.at(k).eta)<2.4 and jets.at(k).pT>25.0 and jets.at(k).PU_mva_loose==1)
     {
       Jet.JetLV.SetPtEtaPhiE(jets.at(k).pT, jets.at(k).eta, jets.at(k).phi, jets.at(k).energy);
       Jet.BTag_CSV = jets.at(k).BTag_csv;
       bool isGoodJet=true;
       for(unsigned int j=0; j<electrons.size(); ++j)
       {
        TLorentzVector Electron;
        if(electrons.at(j).pT > 7.0 and electrons.at(j).isLooseWG==1 and electrons.at(j).isolation < 0.50 and electrons.at(j).isolation*electrons.at(j).pT < 5.0 and fabs(electrons.at(j).dxy)<0.01 and fabs(electrons.at(j).dz)<0.1)
          {
          Electron.SetPtEtaPhiE(electrons.at(j).pT, electrons.at(j).eta, electrons.at(j).phi, electrons.at(j).energy);
          double DRjet_el = Jet.JetLV.DeltaR(Electron);
          if(DRjet_el<0.5) isGoodJet=false;
          }
        }
   for(unsigned int m=0; m<muons.size(); ++m)
     {
       TLorentzVector Muon;
       if(muons.at(m).pT > 3.0 and muons.at(m).isLoose==1 and muons.at(m).isLooseWG==1 and muons.at(m).isolation < 0.50 and muons.at(m).isolation*muons.at(m).pT < 5.0 and fabs(muons.at(m).dxy)<0.01 and fabs(muons.at(m).dz)<0.1){
       Muon.SetPtEtaPhiE(muons.at(m).pT, muons.at(m).eta, muons.at(m).phi, muons.at(m).energy);
       double DRjet_mu = Jet.JetLV.DeltaR(Muon);
       if(DRjet_mu<0.5) isGoodJet=false;
     }
   }
   for(unsigned int l=0; l<photons.size(); ++l)
     {
       TLorentzVector Photon;
       swissCross.push_back((photons.at(l).phseedCrystalEnergy)/(photons.at(l).phe1x5 + photons.at(l).phe5x1 - photons.at(l).phseedCrystalEnergy));
       if(returnAnaPhoton(photons.at(l))==true and swissCross.at(l)<0.90)
       //if(fabs(photons.at(l).eta)<1.444 and photons.at(l).pT>30.0 and photons.at(l).isLoose==1 and photons.at(l).phIsoLoose==1 and photons.at(l).phpixelVeto==0 and photons.at(l).phSigmaIetaIeta > 0.001 and swissCross.at(l)<0.90 and std::sqrt(photons.at(l).phSigmaIphiIphi)>0.009 and photons.at(l).phR9 < 1.0)
       {
       //if(fabs(photons.at(l).eta)<1.444 and photons.at(l).pT>30.0 and photons.at(l).isLoose==1 and photons.at(l).phIsoLoose==1 and photons.at(l).phpixelVeto==0){// and photons.at(l).phSigmaIetaIeta > 0.001 and photons.at(l).phSigmaIphiIphi>0.0001 and swissCross.at(l)<0.90){ //photons.at(l).phSigmaIphiIphi>0.0001 and swissCross.at(l)<0.90){
         Photon.SetPtEtaPhiE(photons.at(l).pT, photons.at(l).eta, photons.at(l).phi, photons.at(l).energy);
         double DRjet_ph = Jet.JetLV.DeltaR(Photon);
         if(DRjet_ph<0.5) isGoodJet=false;
     }
   }
   if(type=="Data")
     {
     TLorentzVector trigger1_p4;
     trigger1_p4.SetPxPyPzE(trigger1.Px, trigger1.Py, trigger1.Pz, trigger1.E);
     double DRjet_tr = Jet.JetLV.DeltaR(trigger1_p4);
     if(DRjet_tr<0.5) isGoodJet=false;
     }

    if(isGoodJet) Jets.push_back(Jet);
     }//close four vector if
  }//close jet loop

  // Now sorting this vector of structs
  std::sort (Jets.begin(), Jets.end(), sortJetVectorsInDescendingpT);
  
  std::vector<BJetInfo> Bjets;
  Bjets.clear();
  for(unsigned int m=0; m<Jets.size(); m++){
    HT += Jets.at(m).JetLV.Pt();
    BJetInfo Bjet;
    if(Jets.at(m).JetLV.Pt() > 25.0 and Jets.at(m).BTag_CSV > 0.679){
      Bjet.BJetLV.SetPtEtaPhiE(Jets.at(m).JetLV.Pt(), Jets.at(m).JetLV.Eta(), Jets.at(m).JetLV.Phi(), Jets.at(m).JetLV.E());
      Bjet.BTag=Jets.at(m).BTag_CSV;
      Bjets.push_back(Bjet);
    }
  }

  // Now sorting this vector of structs
  std::sort (Bjets.begin(), Bjets.end(), sortBJetVectorsInDescendingpT);

  double HTb = 0.0;
  for(unsigned int n=0; n<Bjets.size(); n++) HTb += Bjets.at(n).BJetLV.Pt();

  std::vector<AnalysisLeptonInfo> Electrons;
  Electrons.clear();
  for (unsigned int j=0; j<electrons.size(); ++j)
     {
     AnalysisLeptonInfo Electron;
     if(electrons.at(j).pT > 7.0 and electrons.at(j).isLooseWG==1 and electrons.at(j).isolation < 0.50 and electrons.at(j).isolation*electrons.at(j).pT < 5.0 and fabs(electrons.at(j).dxy)<0.01 and fabs(electrons.at(j).dz)<0.1)  
       {
       Electron.LepLV.SetPtEtaPhiE(electrons.at(j).pT, electrons.at(j).eta, electrons.at(j).phi, electrons.at(j).energy);
       Electron.Charge=electrons.at(j).charge;
       Electron.Isolation=electrons.at(j).isolation;
       Electron.Dxy=electrons.at(j).dxy;
       Electron.Dz=electrons.at(j).dz;
       if(MCSample!="Signal"){
         if(type=="MC" and MCSample=="ZGToLLG") Electron.Matched = electrons.at(j).matched;
         if(type=="MC" and MCSample=="ZGToLLG") Electron.Mother1 = electrons.at(j).mother1;
         if(type=="MC" and MCSample=="ZGToLLG") Electron.Mother2 = electrons.at(j).mother2;
         if(type=="MC" and MCSample=="ZGToLLG") Electron.Mother3 = electrons.at(j).mother3;
         if(type=="MC" and MCSample=="ZGToLLG") Electron.Mother4 = electrons.at(j).mother4;
         if(type=="MC" and MCSample=="ZGToLLG") Electron.Mother5 = electrons.at(j).mother5;
         if(type=="MC" and MCSample=="ZGToLLG") Electron.Mother6 = electrons.at(j).mother6;
         if(type=="MC" and MCSample=="ZGToLLG") Electron.Mother7 = electrons.at(j).mother7;
         if(type=="MC" and MCSample=="ZGToLLG") Electron.Mother8 = electrons.at(j).mother8;
         if(type=="MC" and MCSample=="ZGToLLG") Electron.Mother9 = electrons.at(j).mother9;
         if(type=="MC" and MCSample=="ZGToLLG") Electron.Mother10 = electrons.at(j).mother10;
       } 
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

   std::vector<AnalysisLeptonInfo> Muons;
   Muons.clear();
   for(unsigned int j=0; j<muons.size(); ++j)
     {
     AnalysisLeptonInfo Muon;
     if(muons.at(j).pT > 3.0 and muons.at(j).isLoose==1 and muons.at(j).isLooseWG==1 and muons.at(j).isolation < 0.50 and muons.at(j).isolation*muons.at(j).pT < 5.0 and fabs(muons.at(j).dxy)<0.01 and fabs(muons.at(j).dz)<0.1)
       {
       Muon.LepLV.SetPtEtaPhiE(muons.at(j).pT, muons.at(j).eta, muons.at(j).phi, muons.at(j).energy);
       Muon.Charge=muons.at(j).charge;
       Muon.Isolation=muons.at(j).isolation;
       Muon.Dxy=muons.at(j).dxy;
       Muon.Dz=muons.at(j).dz;
       if(MCSample!="Signal"){
         if(type=="MC" and MCSample=="ZGToLLG") Muon.Matched = muons.at(j).matched;
         if(type=="MC" and MCSample=="ZGToLLG") Muon.Mother1 = muons.at(j).mother1;
         if(type=="MC" and MCSample=="ZGToLLG") Muon.Mother2 = muons.at(j).mother2;
         if(type=="MC" and MCSample=="ZGToLLG") Muon.Mother3 = muons.at(j).mother3;
         if(type=="MC" and MCSample=="ZGToLLG") Muon.Mother4 = muons.at(j).mother4;
         if(type=="MC" and MCSample=="ZGToLLG") Muon.Mother5 = muons.at(j).mother5;
         if(type=="MC" and MCSample=="ZGToLLG") Muon.Mother6 = muons.at(j).mother6;
         if(type=="MC" and MCSample=="ZGToLLG") Muon.Mother7 = muons.at(j).mother7;
         if(type=="MC" and MCSample=="ZGToLLG") Muon.Mother8 = muons.at(j).mother8;
         if(type=="MC" and MCSample=="ZGToLLG") Muon.Mother9 = muons.at(j).mother9;
         if(type=="MC" and MCSample=="ZGToLLG") Muon.Mother10 = muons.at(j).mother10;
       }
       bool isGoodMuon = true;
       for (unsigned int k=0; k<Electrons.size(); k++){
         TLorentzVector AnaElectron;
         AnaElectron.SetPtEtaPhiE(Electrons.at(k).LepLV.Pt(), Electrons.at(k).LepLV.Eta(), Electrons.at(k).LepLV.Phi(), Electrons.at(k).LepLV.E());
         double DRel_mu = Muon.LepLV.DeltaR(AnaElectron);
         if(DRel_mu<0.3) isGoodMuon=false;//anti-matching to muons object.
       }
       if(isGoodMuon) Muons.push_back(Muon);
       }//close four vector if
     }//close muon loop

  std::sort(Muons.begin(), Muons.end(), sortVectorsInDescendingpT);

  double caloMET_True = caloMET;
  if(Muons.size() > 1){
    double MEx = MET*cos(MET_Phi) + Muons.at(0).LepLV.Px() + Muons.at(1).LepLV.Px();
    double MEy = MET*sin(MET_Phi) + Muons.at(0).LepLV.Py() + Muons.at(1).LepLV.Py();
    double caloMet = sqrt(MEx*MEx + MEy*MEy);
  }

  std::vector<GenParticleInfo> taus;
  std::vector<GenParticleInfo> mus;
  std::vector<GenParticleInfo> es;
  if(type=="MC" and MCSample=="ZGToLLG"){
    for (unsigned int j=0; j<tauPx->size(); ++j)
      {
      GenParticleInfo tau;
      tau.Px = tauPx->at(j);
      tau.Py = tauPy->at(j); 
      tau.Pz = tauPz->at(j);
      tau.E = tauE->at(j);
      taus.push_back(tau);
     }
    for (unsigned int j=0; j<muPx->size(); ++j)
      {
      GenParticleInfo mu;
      mu.Px = muPx->at(j);
      mu.Py = muPy->at(j);
      mu.Pz = muPz->at(j);
      mu.E = muE->at(j);
      mus.push_back(mu);
      }
    for (unsigned int j=0; j<ePx->size(); ++j)
      {
      GenParticleInfo ee;
      ee.Px = ePx->at(j);
      ee.Py = ePy->at(j);
      ee.Pz = ePz->at(j);
      ee.E = eE->at(j);
      es.push_back(ee);
     }
  }
  std::sort(taus.begin(), taus.end(), sortGenParticlesInDescendingpT);
  std::sort(mus.begin(), mus.end(), sortGenParticlesInDescendingpT);
  std::sort(es.begin(), es.end(), sortGenParticlesInDescendingpT);

  //Preliminary event cuts guided by the trigger:
  if (not(Photon_vector.size() > 0 and Photon_vector.at(0).Pt() > 32.0 and caloMET > 35.0)) continue;
  nEvents_Trigger++;
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
      double PU_weights = 1.0;
      //if(MCSample!="Signal") 
      PU_weights = h_nWeights->GetBinContent(nVertices);
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
 
  std::vector<svFitStandalone::MeasuredTauLepton> measuredTauLeptons;
  svFitStandalone::Vector measuredMET(MET_Px, MET_Py, 0);
  svFitStandalone::LorentzVector mu1P4_xyzt, mu2P4_xyzt, el1P4_xyzt;
  svFitStandalone::LorentzVector tau1_xyzt, tau2_xyzt; 
  TMatrixD cov(2,2);
  cov[0][0] = MET_Signxx;
  cov[0][1] = MET_Signxy; 
  cov[1][0] = MET_Signyx;
  cov[1][1] = MET_Signyy;

  if(eventWeight > 0.0 and Muons.size()==2){ 
      nEvents_Dimuon++;
      if(type=="MC") eventWeight *= muonSF_WG(Muons.at(0).LepLV.Pt(), Muons.at(0).LepLV.Eta())*muonSF_WG(Muons.at(1).LepLV.Pt(), Muons.at(1).LepLV.Eta())*photonSF_WG(Photon_vector.at(0).Pt(), Photon_vector.at(0).Eta());
      mu1_transverse.SetMagPhi(Muons.at(0).LepLV.Pt(), Muons.at(0).LepLV.Phi());
      mu2_transverse.SetMagPhi(Muons.at(1).LepLV.Pt(), Muons.at(1).LepLV.Phi()); 
      if((Muons.at(0).LepLV+Muons.at(1).LepLV).M()< 12.0) continue;
      if((Muons.at(0).LepLV+Muons.at(1).LepLV).M() > 70 and (Muons.at(0).LepLV+Muons.at(1).LepLV).M() < 110) mumu_OS_Z+=eventWeight;
      /*mu1P4_xyzt.SetPxPyPzE(Muons.at(0).LepLV.Px(), Muons.at(0).LepLV.Py(), Muons.at(0).LepLV.Pz(), Muons.at(0).LepLV.E());
      mu2P4_xyzt.SetPxPyPzE(Muons.at(1).LepLV.Px(), Muons.at(1).LepLV.Py(), Muons.at(1).LepLV.Pz(), Muons.at(1).LepLV.E());
      measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(svFitStandalone::kTauToMuDecay, mu1P4_xyzt)); //kTauToMuDecay
      measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(svFitStandalone::kTauToMuDecay, mu2P4_xyzt)); //kTauToMuDecay
      SVfitStandaloneAlgorithm algo(measuredTauLeptons, measuredMET, cov, 0);
      algo.addLogM(false);
      //algo.integrateVEGAS();  
      algo.integrateMarkovChain();*/ //new method to compute SVFit mass 
      double svfit_mass = 0;
      //svfit_mass = algo.getMass();
      TLorentzVector SVFittedVector;
      //SVFittedVector.SetPtEtaPhiM(algo.pt(), algo.eta(), algo.phi(), algo.getMass());
      SVFittedVector.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
      if(signSelection=="OS" and (Muons.at(0).Charge*Muons.at(1).Charge)==-1){
        if((Muons.at(0).LepLV+Muons.at(1).LepLV+Photon_vector.at(0)).M() < 80 or (Muons.at(0).LepLV+Muons.at(1).LepLV+Photon_vector.at(0)).M() > 100){ 
        //cout << "Event MuMu = " << event << " Run = "<< run << " Lumi = " << lumi << endl;
        mumu_OS+=eventWeight;
        bool eventVeto = false;
       for(int k=0; k<Jets.size(); k++) if(Jets.at(k).JetLV.Pt() > 40 and Jets.at(k).BTag_CSV > 0.679) {eventVeto=true; break;}
       if(eventVeto==true) vetoCounter+=eventWeight;
     
       unsigned int u_ht=0;
       for(float ht=30; ht<730;ht+=20){
          if(HT<float(ht)) HT_Yield_MuMu[u_ht] += eventWeight;
          if(HT<float(ht)) h_HT_LowCutValue_mumu->Fill(ht, eventWeight);
          ++u_ht;
       }

       unsigned int t_htb=0;
       for(float htb=30; htb<730; htb+=20.0){
          if(HTb<htb) HTb_Yield_MuMu[t_htb] += eventWeight;
          if(HTb<htb) h_HTb_LowCutValue_mumu->Fill(htb, eventWeight);
          ++t_htb;
       }
       unsigned int s_invMass=0;
       for(float invMass=0.0; invMass<300; invMass+=10){
          if((Muons.at(0).LepLV+Muons.at(1).LepLV).M()<invMass) InvMass_Yield_MuMu[s_invMass] += eventWeight;
          ++s_invMass;
       }
       unsigned int i_ht=0;
       for(float ht=30.0; ht<730.0;ht+=20.0){
         HT_Value_MuMu[i_ht] = ht;
         unsigned int j_invMass=0;
         for(float invMass=0.0; invMass<300.0; invMass+=10.0){
            InvMass_Value_MuMu[j_invMass]=invMass;
            if(HT<ht and ((Muons.at(0).LepLV+Muons.at(1).LepLV).M()<invMass)) HT_InvMass_Yield_MuMu[i_ht][j_invMass]+=eventWeight;
            if(HT<ht and ((Muons.at(0).LepLV+Muons.at(1).LepLV).M()<invMass)) h_HT_Mmumu_LowCutValue->Fill(ht, invMass, eventWeight); 
            if(HT<ht and ((Muons.at(0).LepLV+Muons.at(1).LepLV).M()>invMass)) h_HT_Mmumu_HighCutValue->Fill(ht, invMass, eventWeight);
            unsigned int k_met=0;
             for(float met=0.0; met<100.0; met+=5.0){
               if(HT<ht and ((Muons.at(0).LepLV+Muons.at(1).LepLV).M()<invMass) and MET>met) HT_InvMass_MET_Yield_MuMu[i_ht][j_invMass][k_met]+=eventWeight;
               ++k_met;
             }
          ++j_invMass;
          }
         
         unsigned int j_svfitMass=0;
         for(float svfitMass=0.0; svfitMass<500.0; svfitMass+=20.0){
            SVFitMass_Value_MuMu[j_svfitMass]=svfitMass;  
            if(HT<ht and (svfit_mass<svfitMass)) HT_svfitMass_Yield_MuMu[i_ht][j_svfitMass]+=eventWeight;
            if(HT<ht and (svfit_mass<svfitMass)) h_HT_SVFitmumu_LowCutValue->Fill(ht, svfitMass, eventWeight);
            ++j_svfitMass;
          }
         unsigned int j_invMassGamma=0;
         for(float invMassGamma=0.0; invMassGamma<300.0; invMassGamma+=20.0){
            InvMassGamma_Value_MuMu[j_invMassGamma]=invMassGamma;
            if(HT<ht and ((Muons.at(0).LepLV+Muons.at(1).LepLV+Photon_vector.at(0)).M()<invMassGamma)) HT_InvMassGamma_Yield_MuMu[i_ht][j_invMassGamma]+=eventWeight;
            if(HT<ht and ((Muons.at(0).LepLV+Muons.at(1).LepLV+Photon_vector.at(0)).M()<invMassGamma)) h_HT_MmumuGamma_LowCutValue->Fill(ht, invMassGamma, eventWeight);
            ++j_invMassGamma;
          }
        unsigned int j_svfitMassGamma=0;
        for(float svfitMassGamma=0.0; svfitMassGamma<500.0; svfitMassGamma+=20.0){
            SVFitMassGamma_Value_MuMu[j_svfitMassGamma]=svfitMassGamma;
            if(HT<ht and ((Muons.at(0).LepLV+Muons.at(1).LepLV).M()<svfitMassGamma)) HT_svfitMassGamma_Yield_MuMu[i_ht][j_svfitMassGamma]+=eventWeight;
            if(HT<ht and ((SVFittedVector+Photon_vector.at(0)).M()<svfitMassGamma)) h_HT_SVFitmumuGamma_LowCutValue->Fill(ht, svfitMassGamma, eventWeight);
            ++j_svfitMassGamma;
          }
         ++i_ht;
       }

       unsigned int i_htb=0;
       for(float htb=30.0; htb<730.0;htb+=20.0){
         HTb_Value_MuMu[i_htb] = htb;
         unsigned int j_invMass=0;
         for(float invMass=0.0; invMass<300.0; invMass+=10.0){
            InvMass_Value_MuMu[j_invMass]=invMass;
            if(HTb<htb and ((Muons.at(0).LepLV+Muons.at(1).LepLV).M()<invMass)) HTb_InvMass_Yield_MuMu[i_htb][j_invMass]+=eventWeight;
            if(HTb<htb and ((Muons.at(0).LepLV+Muons.at(1).LepLV).M()<invMass)) h_HTb_Mmumu_LowCutValue->Fill(htb, invMass, eventWeight);
            if(HTb<htb and ((Muons.at(0).LepLV+Muons.at(1).LepLV).M()>invMass)) h_HTb_Mmumu_HighCutValue->Fill(htb, invMass, eventWeight);
            unsigned int k_met=0;
             for(float met=0.0; met<100.0; met+=5.0){
               if(HTb<htb and ((Muons.at(0).LepLV+Muons.at(1).LepLV).M()<invMass) and MET>met) HTb_InvMass_MET_Yield_MuMu[i_htb][j_invMass][k_met]+=eventWeight;
               ++k_met;
             }
          ++j_invMass;
          }
         unsigned int j_svfitMass=0;
         for(float svfitMass=0.0; svfitMass<500.0; svfitMass+=20.0){
            SVFitMass_Value_MuMu[j_svfitMass]=svfitMass;
            if(HTb<htb and (svfit_mass<svfitMass)) HTb_svfitMass_Yield_MuMu[i_htb][j_svfitMass]+=eventWeight;
            if(HTb<htb and (svfit_mass<svfitMass)) h_HTb_SVFitmumu_LowCutValue->Fill(htb, svfitMass, eventWeight);
            ++j_svfitMass;
          }
         unsigned int j_invMassGamma=0;
         for(float invMassGamma=0.0; invMassGamma<300.0; invMassGamma+=20.0){
            InvMassGamma_Value_MuMu[j_invMassGamma]=invMassGamma;
            if(HTb<htb and ((Muons.at(0).LepLV+Muons.at(1).LepLV+Photon_vector.at(0)).M()<invMassGamma)) HTb_InvMassGamma_Yield_MuMu[i_htb][j_invMassGamma]+=eventWeight;
            if(HTb<htb and ((Muons.at(0).LepLV+Muons.at(1).LepLV+Photon_vector.at(0)).M()<invMassGamma)) h_HTb_MmumuGamma_LowCutValue->Fill(htb, invMassGamma, eventWeight);
            ++j_invMassGamma;
          }
         unsigned int j_svfitMassGamma=0;
         for(float svfitMassGamma=0.0; svfitMassGamma<500.0; svfitMassGamma+=20.0){
            SVFitMassGamma_Value_MuMu[j_svfitMassGamma]=svfitMassGamma;
            if(HTb<htb and ((SVFittedVector+Photon_vector.at(0)).M()<svfitMassGamma)) HTb_svfitMassGamma_Yield_MuMu[i_htb][j_svfitMassGamma]+=eventWeight;
            if(HTb<htb and ((SVFittedVector+Photon_vector.at(0)).M()<svfitMassGamma)) h_HTb_SVFitmumuGamma_LowCutValue->Fill(htb, svfitMassGamma, eventWeight);
              ++j_svfitMassGamma;
              }
           ++i_htb;
            } 
        }//ZGamma 
      }//OS mode
     else if(signSelection=="SS" and (Muons.at(0).Charge*Muons.at(1).Charge)==+1){
       if((Muons.at(0).LepLV+Muons.at(1).LepLV).M()< 15.0) continue;
       if((Muons.at(0).LepLV+Muons.at(1).LepLV+Photon_vector.at(0)).M() < 80 or (Muons.at(0).LepLV+Muons.at(1).LepLV+Photon_vector.at(0)).M() > 100){
       mumu_SS+=eventWeight;
       bool eventVeto = false;
       for(int k=0; k<Jets.size(); k++) if(Jets.at(k).JetLV.Pt() > 40 and Jets.at(k).BTag_CSV > 0.679) {eventVeto=true; break;}
       if(eventVeto==true) vetoCounter+=eventWeight;
       
        unsigned int u_ht=0;
       for(float ht=30; ht<730;ht+=20){
          if(HT<ht) HT_Yield_MuMu[u_ht] += eventWeight;
          if(HT<ht) h_HT_LowCutValue_mumu->Fill(ht, eventWeight);
          ++u_ht;
       }

       unsigned int t_htb=0;
       for(float htb=30; htb<730; htb+=20.0){
          if(HTb<htb) HTb_Yield_MuMu[t_htb] += eventWeight;
          if(HTb<htb) h_HTb_LowCutValue_mumu->Fill(htb, eventWeight);
          ++t_htb;
       }
       unsigned int s_invMass=0;
       for(int invMass=0.0; invMass<300; invMass+=10){
          if((Muons.at(0).LepLV+Muons.at(1).LepLV).M()<invMass) InvMass_Yield_MuMu[s_invMass] += eventWeight;
          ++s_invMass;
       }

       unsigned int i_ht=0;
       for(float ht=30.0; ht<730.0;ht+=20.0){
         HT_Value_MuMu[i_ht] = ht;
         unsigned int j_invMass=0;
         for(float invMass=0.0; invMass<300.0; invMass+=10.0){
            InvMass_Value_MuMu[j_invMass]=invMass;
            if(HT<ht and ((Muons.at(0).LepLV+Muons.at(1).LepLV).M()<invMass)) HT_InvMass_Yield_MuMu[i_ht][j_invMass]+=eventWeight;
            if(HT<ht and ((Muons.at(0).LepLV+Muons.at(1).LepLV).M()<invMass)) h_HT_Mmumu_LowCutValue->Fill(ht, invMass, eventWeight); 
            if(HT<ht and ((Muons.at(0).LepLV+Muons.at(1).LepLV).M()>invMass)) h_HT_Mmumu_HighCutValue->Fill(ht, invMass, eventWeight);
            unsigned int k_met=0;
             for(float met=0.0; met<100.0; met+=5.0){
               if(HT<ht and ((Muons.at(0).LepLV+Muons.at(1).LepLV).M()<invMass) and MET>met) HT_InvMass_MET_Yield_MuMu[i_ht][j_invMass][k_met]+=eventWeight;
               ++k_met;
             }
          ++j_invMass;
          }
         unsigned int j_svfitMass=0;
         for(float svfitMass=0.0; svfitMass<500.0; svfitMass+=20.0){
            SVFitMass_Value_MuMu[j_svfitMass]=svfitMass;  
            if(HT<ht and (svfit_mass<svfitMass)) HT_svfitMass_Yield_MuMu[i_ht][j_svfitMass]+=eventWeight;
            if(HT<ht and (svfit_mass<svfitMass)) h_HT_SVFitmumu_LowCutValue->Fill(ht, svfitMass, eventWeight);
            ++j_svfitMass;
          }
         unsigned int j_invMassGamma=0;
         for(float invMassGamma=0.0; invMassGamma<300.0; invMassGamma+=20.0){
            InvMassGamma_Value_MuMu[j_invMassGamma]=invMassGamma;
            if(HT<ht and ((Muons.at(0).LepLV+Muons.at(1).LepLV+Photon_vector.at(0)).M()<invMassGamma)) HT_InvMassGamma_Yield_MuMu[i_ht][j_invMassGamma]+=eventWeight;
            if(HT<ht and ((Muons.at(0).LepLV+Muons.at(1).LepLV+Photon_vector.at(0)).M()<invMassGamma)) h_HT_MmumuGamma_LowCutValue->Fill(ht, invMassGamma, eventWeight);
            ++j_invMassGamma;
          }
        unsigned int j_svfitMassGamma=0;
        for(float svfitMassGamma=0.0; svfitMassGamma<500.0; svfitMassGamma+=20.0){
            SVFitMassGamma_Value_MuMu[j_svfitMassGamma]=svfitMassGamma;
            if(HT<ht and ((SVFittedVector+Photon_vector.at(0)).M()<svfitMassGamma)) HT_svfitMassGamma_Yield_MuMu[i_ht][j_svfitMassGamma]+=eventWeight;
            if(HT<ht and ((SVFittedVector+Photon_vector.at(0)).M()<svfitMassGamma)) h_HT_SVFitmumuGamma_LowCutValue->Fill(ht, svfitMassGamma, eventWeight);
            ++j_svfitMassGamma;
          }
         ++i_ht;
       }

      unsigned int i_htb=0;
       for(float htb=30.0; htb<730.0;htb+=20.0){
         HTb_Value_MuMu[i_htb] = htb;
         unsigned int j_invMass=0;
         for(float invMass=0.0; invMass<300.0; invMass+=10.0){
            InvMass_Value_MuMu[j_invMass]=invMass;
            if(HTb<htb and ((Muons.at(0).LepLV+Muons.at(1).LepLV).M()<invMass)) HTb_InvMass_Yield_MuMu[i_htb][j_invMass]+=eventWeight;
            if(HTb<htb and ((Muons.at(0).LepLV+Muons.at(1).LepLV).M()<invMass)) h_HTb_Mmumu_LowCutValue->Fill(htb, invMass, eventWeight);
            if(HTb<htb and ((Muons.at(0).LepLV+Muons.at(1).LepLV).M()>invMass)) h_HTb_Mmumu_HighCutValue->Fill(htb, invMass, eventWeight);
            unsigned int k_met=0;
             for(float met=0.0; met<100.0; met+=5.0){
               if(HTb<htb and ((Muons.at(0).LepLV+Muons.at(1).LepLV).M()<invMass) and MET>met) HTb_InvMass_MET_Yield_MuMu[i_htb][j_invMass][k_met]+=eventWeight;
               ++k_met;
             }
          ++j_invMass;
          }
         unsigned int j_svfitMass=0;
         for(float svfitMass=0.0; svfitMass<500.0; svfitMass+=20.0){
            SVFitMass_Value_MuMu[j_svfitMass]=svfitMass;
            if(HTb<htb and (svfit_mass<svfitMass)) HTb_svfitMass_Yield_MuMu[i_htb][j_svfitMass]+=eventWeight;
            if(HTb<htb and (svfit_mass<svfitMass)) h_HTb_SVFitmumu_LowCutValue->Fill(htb, svfitMass, eventWeight);
            ++j_svfitMass;
          }
         unsigned int j_invMassGamma=0;
         for(float invMassGamma=0.0; invMassGamma<300.0; invMassGamma+=20.0){
            InvMassGamma_Value_MuMu[j_invMassGamma]=invMassGamma;
            if(HTb<htb and ((Muons.at(0).LepLV+Muons.at(1).LepLV+Photon_vector.at(0)).M()<invMassGamma)) HTb_InvMassGamma_Yield_MuMu[i_htb][j_invMassGamma]+=eventWeight;
            if(HTb<htb and ((Muons.at(0).LepLV+Muons.at(1).LepLV+Photon_vector.at(0)).M()<invMassGamma)) h_HTb_MmumuGamma_LowCutValue->Fill(htb, invMassGamma, eventWeight);
            ++j_invMassGamma;
          }
         unsigned int j_svfitMassGamma=0;
         for(float svfitMassGamma=0.0; svfitMassGamma<500.0; svfitMassGamma+=20.0){
            SVFitMassGamma_Value_MuMu[j_svfitMassGamma]=svfitMassGamma;
            if(HTb<htb and ((SVFittedVector+Photon_vector.at(0)).M()<svfitMassGamma)) HTb_svfitMassGamma_Yield_MuMu[i_htb][j_svfitMassGamma]+=eventWeight;
            if(HTb<htb and ((SVFittedVector+Photon_vector.at(0)).M()<svfitMassGamma)) h_HTb_SVFitmumuGamma_LowCutValue->Fill(htb, svfitMassGamma, eventWeight);
              ++j_svfitMassGamma;
              }
           ++i_htb;
            }
          }//ZGamma
       }//SS mode
   }//first+second tight muon

 //ElMu channel
   if( eventWeight > 0.0 and Muons.size()==1 and Electrons.size()==1) {
     nEvents_ElMu++;
     if(type=="MC") eventWeight *= electronSF_WG(Electrons.at(0).LepLV.Pt(), Electrons.at(0).LepLV.Eta())*muonSF_WG(Muons.at(0).LepLV.Pt(), Muons.at(0).LepLV.Eta())*photonSF_WG(Photon_vector.at(0).Pt(), Photon_vector.at(0).Eta());
     //if(type=="MC") eventWeight *= electronSF(Electrons.at(0).LepLV.Pt(), Electrons.at(0).LepLV.Eta())*muonSF(Muons.at(0).LepLV.Pt(), Muons.at(0).LepLV.Eta())*photonSF(Photon_vector.at(0).Pt(), Photon_vector.at(0).Eta());
     mu1_transverse.SetMagPhi(Muons.at(0).LepLV.Pt(), Muons.at(0).LepLV.Phi());
     el1_transverse.SetMagPhi(Electrons.at(0).LepLV.Pt(), Electrons.at(0).LepLV.Phi());
     /*mu1P4_xyzt.SetPxPyPzE(Muons.at(0).LepLV.Px(), Muons.at(0).LepLV.Py(), Muons.at(0).LepLV.Pz(), Muons.at(0).LepLV.E());
     el1P4_xyzt.SetPxPyPzE(Electrons.at(0).LepLV.Px(), Electrons.at(0).LepLV.Py(), Electrons.at(0).LepLV.Pz(), Electrons.at(0).LepLV.E());
     measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(svFitStandalone::kTauToMuDecay, mu1P4_xyzt)); //kTauToMuDecay
     measuredTauLeptons.push_back(svFitStandalone::MeasuredTauLepton(svFitStandalone::kTauToElecDecay, el1P4_xyzt)); //kTauToElecDecay
     SVfitStandaloneAlgorithm algo(measuredTauLeptons, measuredMET, cov, 0);
     algo.addLogM(false);
     //algo.integrateVEGAS();
     algo.integrateMarkovChain();*/
     double svfit_mass = 0;
     //svfit_mass = algo.getMass();
     TLorentzVector SVFittedVector;
     //SVFittedVector.SetPtEtaPhiM(algo.pt(), algo.eta(), algo.phi(), algo.getMass()); 
     SVFittedVector.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0); 
     if(signSelection=="OS" and ((Muons.at(0).Charge)*(Electrons.at(0).Charge)==-1)){
       if(type=="Data") cout << "Event ElMu = " << event << " Run = "<< run << " Lumi = " << lumi << endl;  
       emu_OS+=eventWeight;
       bool eventVeto = false;
       for(int k=0; k<Jets.size(); k++) if(Jets.at(k).JetLV.Pt() > 40 and Jets.at(k).BTag_CSV > 0.679) {eventVeto=true; break;}
       if(eventVeto==true) vetoCounter+=eventWeight;
     
       unsigned int u_ht=0;
       for(float ht=30; ht<730;ht+=20){
          if(HT<float(ht)) HT_Yield_ElMu[u_ht] += eventWeight;
          if(HT<float(ht)) h_HT_LowCutValue_mumu->Fill(ht, eventWeight);
          ++u_ht;
       }

       unsigned int t_htb=0;
       for(float htb=30; htb<730; htb+=20.0){
          if(HTb<htb) HTb_Yield_ElMu[t_htb] += eventWeight;
          if(HTb<htb) h_HTb_LowCutValue_mumu->Fill(htb, eventWeight);
          ++t_htb;
       }
       unsigned int s_invMass=0;
       for(float invMass=0.0; invMass<300; invMass+=10){
          if((Muons.at(0).LepLV+Electrons.at(0).LepLV).M()<invMass) InvMass_Yield_ElMu[s_invMass] += eventWeight;
          ++s_invMass;
         }
 
       unsigned int i_ht=0;
       for(float ht=30.0; ht<730.0;ht+=20.0){
         HT_Value_ElMu[i_ht] = ht;
         unsigned int j_invMass=0;
         for(float invMass=0.0; invMass<300.0; invMass+=10.0){
            InvMass_Value_ElMu[j_invMass]=invMass;
            if(HT<ht and ((Muons.at(0).LepLV+Electrons.at(0).LepLV).M()<invMass)) HT_InvMass_Yield_ElMu[i_ht][j_invMass]+=eventWeight;
            if(HT<ht and ((Muons.at(0).LepLV+Electrons.at(0).LepLV).M()<invMass)) h_HT_Memu_LowCutValue->Fill(ht, invMass, eventWeight); 
            if(HT<ht and ((Muons.at(0).LepLV+Electrons.at(0).LepLV).M()>invMass)) h_HT_Memu_HighCutValue->Fill(ht, invMass, eventWeight);
            unsigned int k_met=0;
             for(float met=0.0; met<100.0; met+=5.0){
               if(HT<ht and ((Muons.at(0).LepLV+Electrons.at(0).LepLV).M()<invMass) and MET>met) HT_InvMass_MET_Yield_ElMu[i_ht][j_invMass][k_met]+=eventWeight;
               ++k_met;
             }
          ++j_invMass;
          }
         unsigned int j_svfitMass=0;
         for(float svfitMass=0.0; svfitMass<500.0; svfitMass+=20.0){
            SVFitMass_Value_ElMu[j_svfitMass]=svfitMass;  
            if(HT<ht and (svfit_mass<svfitMass)) HT_svfitMass_Yield_ElMu[i_ht][j_svfitMass]+=eventWeight;
            if(HT<ht and (svfit_mass<svfitMass)) h_HT_SVFitemu_LowCutValue->Fill(ht, svfitMass, eventWeight);
            ++j_svfitMass;
          }
         unsigned int j_invMassGamma=0;
         for(float invMassGamma=0.0; invMassGamma<300.0; invMassGamma+=20.0){
            InvMassGamma_Value_ElMu[j_invMassGamma]=invMassGamma;
            if(HT<ht and ((Muons.at(0).LepLV+Electrons.at(0).LepLV+Photon_vector.at(0)).M()<invMassGamma)) HT_InvMassGamma_Yield_ElMu[i_ht][j_invMassGamma]+=eventWeight;
            if(HT<ht and ((Muons.at(0).LepLV+Electrons.at(0).LepLV+Photon_vector.at(0)).M()<invMassGamma)) h_HT_MemuGamma_LowCutValue->Fill(ht, invMassGamma, eventWeight);
            ++j_invMassGamma;
          }
        unsigned int j_svfitMassGamma=0;
        for(float svfitMassGamma=0.0; svfitMassGamma<500.0; svfitMassGamma+=20.0){
            SVFitMassGamma_Value_ElMu[j_svfitMassGamma]=svfitMassGamma;
            if(HT<ht and ((SVFittedVector+Photon_vector.at(0)).M()<svfitMassGamma)) HT_svfitMassGamma_Yield_ElMu[i_ht][j_svfitMassGamma]+=eventWeight;
            if(HT<ht and ((SVFittedVector+Photon_vector.at(0)).M()<svfitMassGamma)) h_HT_SVFitemuGamma_LowCutValue->Fill(ht, svfitMassGamma, eventWeight);
            ++j_svfitMassGamma;
          }
         ++i_ht;
       }

      unsigned int i_htb=0;
       for(float htb=30.0; htb<730.0;htb+=20.0){
         HTb_Value_ElMu[i_htb] = htb;
         unsigned int j_invMass=0;
         for(float invMass=0.0; invMass<300.0; invMass+=10.0){
            InvMass_Value_ElMu[j_invMass]=invMass;
            if(HTb<htb and ((Muons.at(0).LepLV+Electrons.at(0).LepLV).M()<invMass)) HTb_InvMass_Yield_ElMu[i_htb][j_invMass]+=eventWeight;
            if(HTb<htb and ((Muons.at(0).LepLV+Electrons.at(0).LepLV).M()<invMass)) h_HTb_Memu_LowCutValue->Fill(htb, invMass, eventWeight);
            if(HTb<htb and ((Muons.at(0).LepLV+Electrons.at(0).LepLV).M()>invMass)) h_HTb_Memu_HighCutValue->Fill(htb, invMass, eventWeight);
            unsigned int k_met=0;
             for(float met=0.0; met<100.0; met+=5.0){
               if(HTb<htb and ((Muons.at(0).LepLV+Electrons.at(0).LepLV).M()<invMass) and MET>met) HTb_InvMass_MET_Yield_ElMu[i_htb][j_invMass][k_met]+=eventWeight;
               ++k_met;
             }
          ++j_invMass;
          }
         unsigned int j_svfitMass=0;
         for(float svfitMass=0.0; svfitMass<500.0; svfitMass+=20.0){
            SVFitMass_Value_ElMu[j_svfitMass]=svfitMass;
            if(HTb<htb and (svfit_mass<svfitMass)) HTb_svfitMass_Yield_ElMu[i_htb][j_svfitMass]+=eventWeight;
            if(HTb<htb and (svfit_mass<svfitMass)) h_HTb_SVFitemu_LowCutValue->Fill(htb, svfitMass, eventWeight);
            ++j_svfitMass;
          }
         unsigned int j_invMassGamma=0;
         for(float invMassGamma=0.0; invMassGamma<300.0; invMassGamma+=20.0){
            InvMassGamma_Value_ElMu[j_invMassGamma]=invMassGamma;
            if(HTb<htb and ((Muons.at(0).LepLV+Electrons.at(0).LepLV+Photon_vector.at(0)).M()<invMassGamma)) HTb_InvMassGamma_Yield_ElMu[i_htb][j_invMassGamma]+=eventWeight;
            if(HTb<htb and ((Muons.at(0).LepLV+Electrons.at(0).LepLV+Photon_vector.at(0)).M()<invMassGamma)) h_HTb_MemuGamma_LowCutValue->Fill(htb, invMassGamma, eventWeight);
            ++j_invMassGamma;
          }
         unsigned int j_svfitMassGamma=0;
         for(float svfitMassGamma=0.0; svfitMassGamma<500.0; svfitMassGamma+=20.0){
            SVFitMassGamma_Value_ElMu[j_svfitMassGamma]=svfitMassGamma;
            if(HTb<htb and ((SVFittedVector+Photon_vector.at(0)).M()<svfitMassGamma)) HTb_svfitMassGamma_Yield_ElMu[i_htb][j_svfitMassGamma]+=eventWeight;
            if(HTb<htb and ((SVFittedVector+Photon_vector.at(0)).M()<svfitMassGamma)) h_HTb_SVFitemuGamma_LowCutValue->Fill(htb, svfitMassGamma, eventWeight);
            ++j_svfitMassGamma;
          }
         ++i_htb;
        }
    }//OS mode
    else if(signSelection=="SS" and ((Muons.at(0).Charge)*(Electrons.at(0).Charge)==+1)){
       if(type=="Data") cout << "Event ElMu SS = " << event << " Run = "<< run << " Lumi = " << lumi << endl;
       emu_SS+=eventWeight;
       bool eventVeto = false;
       for(int k=0; k<Jets.size(); k++) if(Jets.at(k).JetLV.Pt() > 40 and Jets.at(k).BTag_CSV > 0.679) {eventVeto=true; break;}
       if(eventVeto==true) vetoCounter+=eventWeight;
      
       unsigned int u_ht=0;
         for(float ht=30; ht<730;ht+=20){
            if(HT<float(ht)) HT_Yield_ElMu[u_ht] += eventWeight;
            if(HT<float(ht)) h_HT_LowCutValue_mumu->Fill(ht, eventWeight);
            ++u_ht;
         }

       unsigned int t_htb=0;
         for(float htb=30; htb<730; htb+=20.0){
            if(HTb<htb) HTb_Yield_ElMu[t_htb] += eventWeight;
            if(HTb<htb) h_HTb_LowCutValue_mumu->Fill(htb, eventWeight);
            ++t_htb;
         }
       unsigned int s_invMass=0;
         for(float invMass=0.0; invMass<300; invMass+=10){
           if((Muons.at(0).LepLV+Electrons.at(0).LepLV).M()<invMass) InvMass_Yield_ElMu[s_invMass] += eventWeight;
          ++s_invMass;
         }

       unsigned int i_ht=0;
       for(float ht=30.0; ht<730.0;ht+=20.0){
         HT_Value_ElMu[i_ht] = ht;
         unsigned int j_invMass=0;
         for(float invMass=0.0; invMass<300.0; invMass+=10.0){
            InvMass_Value_ElMu[j_invMass]=invMass;
            if(HT<ht and ((Muons.at(0).LepLV+Electrons.at(0).LepLV).M()<invMass)) HT_InvMass_Yield_ElMu[i_ht][j_invMass]+=eventWeight;
            if(HT<ht and ((Muons.at(0).LepLV+Electrons.at(0).LepLV).M()<invMass)) h_HT_Memu_LowCutValue->Fill(ht, invMass, eventWeight); 
            if(HT<ht and ((Muons.at(0).LepLV+Electrons.at(0).LepLV).M()>invMass)) h_HT_Memu_HighCutValue->Fill(ht, invMass, eventWeight);
            unsigned int k_met=0;
             for(float met=0.0; met<100.0; met+=5.0){
               if(HT<ht and ((Muons.at(0).LepLV+Electrons.at(0).LepLV).M()<invMass) and MET>met) HT_InvMass_MET_Yield_ElMu[i_ht][j_invMass][k_met]+=eventWeight;
               ++k_met;
             }
          ++j_invMass;
          }
         unsigned int j_svfitMass=0;
         for(float svfitMass=0.0; svfitMass<500.0; svfitMass+=20.0){
            SVFitMass_Value_ElMu[j_svfitMass]=svfitMass;  
            if(HT<ht and (svfit_mass<svfitMass)) HT_svfitMass_Yield_ElMu[i_ht][j_svfitMass]+=eventWeight;
            if(HT<ht and (svfit_mass<svfitMass)) h_HT_SVFitemu_LowCutValue->Fill(ht, svfitMass, eventWeight);
            ++j_svfitMass;
          }
         unsigned int j_invMassGamma=0;
         for(float invMassGamma=0.0; invMassGamma<300.0; invMassGamma+=20.0){
            InvMassGamma_Value_ElMu[j_invMassGamma]=invMassGamma;
            if(HT<ht and ((Muons.at(0).LepLV+Electrons.at(0).LepLV+Photon_vector.at(0)).M()<invMassGamma)) HT_InvMassGamma_Yield_ElMu[i_ht][j_invMassGamma]+=eventWeight;
            if(HT<ht and ((Muons.at(0).LepLV+Electrons.at(0).LepLV+Photon_vector.at(0)).M()<invMassGamma)) h_HT_MemuGamma_LowCutValue->Fill(ht, invMassGamma, eventWeight);
            ++j_invMassGamma;
          }
        unsigned int j_svfitMassGamma=0;
        for(float svfitMassGamma=0.0; svfitMassGamma<500.0; svfitMassGamma+=20.0){
            SVFitMassGamma_Value_ElMu[j_svfitMassGamma]=svfitMassGamma;
            if(HT<ht and ((SVFittedVector+Photon_vector.at(0)).M()<svfitMassGamma)) HT_svfitMassGamma_Yield_ElMu[i_ht][j_svfitMassGamma]+=eventWeight;
            if(HT<ht and ((SVFittedVector+Photon_vector.at(0)).M()<svfitMassGamma)) h_HT_SVFitemuGamma_LowCutValue->Fill(ht, svfitMassGamma, eventWeight);
            ++j_svfitMassGamma;
          }
         ++i_ht;
       }

       unsigned int i_htb=0;
       for(float htb=30.0; htb<730.0;htb+=20.0){
         HTb_Value_ElMu[i_htb] = htb;
         unsigned int j_invMass=0;
         for(float invMass=0.0; invMass<300.0; invMass+=10.0){
            InvMass_Value_ElMu[j_invMass]=invMass;
            if(HTb<htb and ((Muons.at(0).LepLV+Electrons.at(0).LepLV).M()<invMass)) HTb_InvMass_Yield_ElMu[i_htb][j_invMass]+=eventWeight;
            if(HTb<htb and ((Muons.at(0).LepLV+Electrons.at(0).LepLV).M()<invMass)) h_HTb_Memu_LowCutValue->Fill(htb, invMass, eventWeight);
            if(HTb<htb and ((Muons.at(0).LepLV+Electrons.at(0).LepLV).M()>invMass)) h_HTb_Memu_HighCutValue->Fill(htb, invMass, eventWeight);
            unsigned int k_met=0;
             for(float met=0.0; met<100.0; met+=5.0){
               if(HTb<htb and ((Muons.at(0).LepLV+Electrons.at(0).LepLV).M()<invMass) and MET>met) HTb_InvMass_MET_Yield_ElMu[i_htb][j_invMass][k_met]+=eventWeight;
               ++k_met;
             }
          ++j_invMass;
          }
         unsigned int j_svfitMass=0;
         for(float svfitMass=0.0; svfitMass<500.0; svfitMass+=20.0){
            SVFitMass_Value_ElMu[j_svfitMass]=svfitMass;
            if(HTb<htb and (svfit_mass<svfitMass)) HTb_svfitMass_Yield_ElMu[i_htb][j_svfitMass]+=eventWeight;
            if(HTb<htb and (svfit_mass<svfitMass)) h_HTb_SVFitemu_LowCutValue->Fill(htb, svfitMass, eventWeight);
            ++j_svfitMass;
          }
         unsigned int j_invMassGamma=0;
         for(float invMassGamma=0.0; invMassGamma<300.0; invMassGamma+=20.0){
            InvMassGamma_Value_ElMu[j_invMassGamma]=invMassGamma;
            if(HTb<htb and ((Muons.at(0).LepLV+Electrons.at(0).LepLV+Photon_vector.at(0)).M()<invMassGamma)) HTb_InvMassGamma_Yield_ElMu[i_htb][j_invMassGamma]+=eventWeight;
            if(HTb<htb and ((Muons.at(0).LepLV+Electrons.at(0).LepLV+Photon_vector.at(0)).M()<invMassGamma)) h_HTb_MemuGamma_LowCutValue->Fill(htb, invMassGamma, eventWeight);
            ++j_invMassGamma;
          }
         unsigned int j_svfitMassGamma=0;
         for(float svfitMassGamma=0.0; svfitMassGamma<500.0; svfitMassGamma+=20.0){
            SVFitMassGamma_Value_ElMu[j_svfitMassGamma]=svfitMassGamma;
            if(HTb<htb and ((SVFittedVector+Photon_vector.at(0)).M()<svfitMassGamma)) HTb_svfitMassGamma_Yield_ElMu[i_htb][j_svfitMassGamma]+=eventWeight;
            if(HTb<htb and ((SVFittedVector+Photon_vector.at(0)).M()<svfitMassGamma)) h_HTb_SVFitemuGamma_LowCutValue->Fill(htb, svfitMassGamma, eventWeight);
            ++j_svfitMassGamma;
          }
         ++i_htb;
        }
     }//SS mode     
   }//El-Mu channel
  
 //ElEl channel
   if( eventWeight > 0.0 and Electrons.size()==2) {
       nEvents_ElEl++;
       if(type=="MC") eventWeight *= electronSF_WG(Electrons.at(0).LepLV.Pt(), Electrons.at(0).LepLV.Eta())*electronSF_WG(Electrons.at(1).LepLV.Pt(), Electrons.at(1).LepLV.Eta())*photonSF_WG(Photon_vector.at(0).Pt(), Photon_vector.at(0).Eta());
     //if(type=="MC") eventWeight *= electronSF(Electrons.at(0).LepLV.Pt(), Electrons.at(0).LepLV.Eta())*electronSF(Electrons.at(1).LepLV.Pt(), Electrons.at(1).LepLV.Eta())*photonSF(Photon_vector.at(0).Pt(), Photon_vector.at(0).Eta());
     el1_transverse.SetMagPhi(Electrons.at(0).LepLV.Pt(), Electrons.at(0).LepLV.Phi());
     el2_transverse.SetMagPhi(Electrons.at(1).LepLV.Pt(), Electrons.at(1).LepLV.Phi());
     if(signSelection=="OS" and ((Electrons.at(0).Charge)*(Electrons.at(1).Charge)==-1)){
       ee_OS+=eventWeight;
    }//OS mode
    else if(signSelection=="SS" and ((Electrons.at(0).Charge)*(Electrons.at(1).Charge)==+1)){
       if(type=="Data") cout << "Event ElEl = " << event << " Run = "<< run << " Lumi = " << lumi << endl;
       ee_SS+=eventWeight;
     }//SS mode     
   }//El-El channel
}//event loop closed

  cout << "nEvents_Incoming = " << nEvents_Incoming << endl;
  cout << "nEvents_Trigger = " << nEvents_Trigger << endl;
  cout << "nEvents_Dimuon = " << nEvents_Dimuon << endl;
  cout << "nEvents_ElMu = " << nEvents_ElMu << endl;
  cout << "nEvents_ElEl = " << nEvents_ElEl << endl;
  //Cut flow table
  cout << "mumu_OS_Z = " << mumu_OS_Z << endl;
  cout << "mumu_OS = " << mumu_OS << endl;
  cout << "mumu_SS = " << mumu_SS << endl;
  cout << "emu_OS = " << emu_OS << endl;
  cout << "emu_SS = " << emu_SS << endl;
  cout << "ee_OS = " << ee_OS << endl;
  cout << "ee_SS = " << ee_SS << endl;

  cout << "Only HT optimization with HTb < 75"<< endl;

  for(int s=0; s<35.00; s++){
    std::cout << "HT = " << s*20+30 << " and HT_Yield ElMu = " << HT_Yield_ElMu[s] << std::endl;
  }

  cout << "Only Invariant mass optimization with HTb < 75" << endl;

  for(int t=0; t<35.00; t++){
    std::cout << "Inv mass = " << t*10 << " and InvMass_Yield ElMu = " << InvMass_Yield_ElMu[t] << std::endl;
  }

  cout << "2D HT and invariant mass optimization" << endl;

  for(int l=0; l<35.0; l++){
    for(int m=0; m<35.0; m++){
      cout <<  "HT = " << l*20+30 << " and " << " inv mass = " << m*10 << " and HT_InvMass_Yield ElMu = " << HT_InvMass_Yield_ElMu[l][m] << std::endl;
    }
  }


 for(int l=0; l<35.0; l++){
    for(int m=0; m<35.0; m++){
      cout << "HT = " << l*20+30 << " and " << " inv mass = " << m*20 << " and HT_InvMassGamma_Yield ElMu = " << HT_InvMassGamma_Yield_ElMu[l][m] << std::endl;
    }
  }


  for(int l=0; l<35.0; l++){
    for(int m=0; m<35.0; m++){
      cout << "HT = " << l*20+30 << " and " << " inv mass = " << m*20 << " and HT_svfitMass_Yield ElMu = " << HT_svfitMass_Yield_ElMu[l][m] << std::endl;
    }
  }


  for(int l=0; l<35.0; l++){
    for(int m=0; m<35.0; m++){
      cout << "HT = " << l*20+30 << " and " << " inv mass = " << m*20 << " and HT_svfitMassGamma_Yield ElMu = " << HT_svfitMassGamma_Yield_ElMu[l][m] << std::endl;
    }
  }

  cout << "3D HT, invariant mass and met optimization" << endl;

  for(int l=0; l<35.0; l++){
    for(int m=0; m<35.0; m++){
      for(int n=0; n<35.0; n++){
        cout << "HT = " << l*20+30 << " and " << " inv mass = " << m*10 << " met = " << n*5 << " and HT_InvMass_MET_Yield ElMu = " << HT_InvMass_MET_Yield_ElMu[l][m][n] << std::endl;
      }
    }
  }

   cout << "Only HTb optimization with HTb < 75" << endl;

   for(int s=0; s<35.00; s++){
     std::cout << "HTb = " << s*20+30 << " and HTb_Yield ElMu = " << HTb_Yield_ElMu[s] << std::endl;
   }

   cout << "2D HTb and invariant mass optimization" << endl;

  for(int l=0; l<35.0; l++){
    for(int m=0; m<35.0; m++){
      cout <<  "HTb = " << l*20+30 << " and " << " inv mass = " << m*10 << " and HTb_InvMass_Yield ElMu = " << HTb_InvMass_Yield_ElMu[l][m] << std::endl;
    }
  }

  for(int l=0; l<35.0; l++){
    for(int m=0; m<35.0; m++){
      cout << "HTb = " << l*20+30 << " and " << " inv mass = " << m*20 << " and HTb_InvMassGamma_Yield ElMu = " << HTb_InvMassGamma_Yield_ElMu[l][m] << std::endl;
    }
  }


  for(int l=0; l<35.0; l++){
    for(int m=0; m<35.0; m++){
      cout << "HTb = " << l*20+30 << " and " << " inv mass = " << m*20 << " and HTb_svfitMass_Yield ElMu = " << HTb_svfitMass_Yield_ElMu[l][m] << std::endl;
    }
  }


  for(int l=0; l<35.0; l++){
    for(int m=0; m<35.0; m++){
      cout << "HTb = " << l*20+30 << " and " << " inv mass = " << m*20 << " and HTb_svfitMassGamma_Yield ElMu = " << HTb_svfitMassGamma_Yield_ElMu[l][m] << std::endl;
    }
  }


  cout << "3D HTb, invariant mass and met optimization" << endl;

  for(int l=0; l<35.0; l++){
    for(int m=0; m<35.0; m++){
      for(int n=0; n<35.0; n++){
        cout << "HTb = " << l*20+30 << " and " << " inv mass = " << m*10 << " met = " << n*5 << " and HTb_InvMass_MET_Yield ElMu = " << HTb_InvMass_MET_Yield_ElMu[l][m][n] << std::endl;
      }
    }
  }

  cout << "Only HT optimization with HTb < 75" << endl;

  for(int s=0; s<35.00; s++){
    std::cout << "HT = " << s*20+30 << " and HT_Yield MuMu = " << HT_Yield_MuMu[s] << std::endl;
  }

  cout << "Only Invariant mass optimization with HTb < 75" << endl;

  for(int t=0; t<35.00; t++){
    std::cout << "Inv mass = " << t*10 << " and InvMass_Yield MuMu = " << InvMass_Yield_MuMu[t] << std::endl;
  }

  cout << "2D HT and invariant mass optimization" << endl;

  for(int l=0; l<3.0; l++){
    for(int m=0; m<3.0; m++){
      cout << "HT = " << l*20+30 << " and " << " inv mass = " << m*10 << " and HT_InvMass_Yield MuMu = " << HT_InvMass_Yield_MuMu[l][m] << std::endl;
    }
  }

 cout << "3D HT, invariant mass and met optimization" << endl;

  for(int l=0; l<35.0; l++){
    for(int m=0; m<35.0; m++){
      for(int n=0; n<35.0; n++){
        cout << "HT = " << l*20+30 << " and " << " inv mass = " << m*10 << " met = " << n*5 << " and HT_InvMass_MET_Yield MuMu = " << HT_InvMass_MET_Yield_MuMu[l][m][n] << std::endl;
      }
    }
  }

  cout << "Only HTb optimization with HTb < 75" << endl;

  for(int s=0; s<35.00; s++){
    std::cout << "HTb = " << s*20+30 << " and HTb_Yield MuMu = " << HTb_Yield_MuMu[s] << std::endl;
  }
 
  cout << "2D HTb and invariant mass optimization" << endl;

  for(int l=0; l<35.0; l++){
    for(int m=0; m<35.0; m++){
      cout << "HTb = " << l*20+30 << " and " << " inv mass = " << m*20 << " and HTb_InvMass_Yield MuMu = " << HTb_InvMass_Yield_MuMu[l][m] << std::endl;
    }
  }

  cout << "3D HTb, invariant mass and met optimization" << endl;

  for(int l=0; l<35.0; l++){
    for(int m=0; m<35.0; m++){
      for(int n=0; n<35.0; n++){
        cout << "HTb = " << l*20+30 << " and " << " inv mass = " << m*10 << " met = " << n*5 << " and HTb_InvMass_MET_Yield MuMu = " << HTb_InvMass_MET_Yield_MuMu[l][m][n] << std::endl;
      }
    }
  }

  std::string histfilename=(outfile+".root").c_str();
  TFile *tFile=new TFile(histfilename.c_str(), "RECREATE");
  h_HT_LowCutValue_mumu->Write();
  h_HT_LowCutValue_emu->Write();
  h_HTb_LowCutValue_mumu->Write();
  h_HTb_LowCutValue_emu->Write();
  h_HT_Mmumu_LowCutValue->Write();
  h_HT_Mmumu_HighCutValue->Write();
  h_HT_MmumuGamma_LowCutValue->Write();
  h_HT_Memu_LowCutValue->Write(); 
  h_HT_Memu_HighCutValue->Write();
  h_HT_MemuGamma_LowCutValue->Write();
  h_HTb_Mmumu_LowCutValue->Write();
  h_HTb_Mmumu_HighCutValue->Write();
  h_HTb_MmumuGamma_LowCutValue->Write();
  h_HTb_Memu_LowCutValue->Write();
  h_HTb_Memu_HighCutValue->Write();
  h_HTb_MemuGamma_LowCutValue->Write();
  h_HT_SVFitmumu_LowCutValue->Write();
  h_HT_SVFitmumuGamma_LowCutValue->Write();
  h_HT_SVFitemu_LowCutValue->Write();
  h_HT_SVFitemuGamma_LowCutValue->Write(); 
  h_HTb_SVFitmumu_LowCutValue->Write();
  h_HTb_SVFitmumuGamma_LowCutValue->Write();
  h_HTb_SVFitemu_LowCutValue->Write();
  h_HTb_SVFitemuGamma_LowCutValue->Write(); 
  tFile->Close(); 
  std::cout<<"Wrote output file "<<histfilename<<std::endl;

  return 0;

}
