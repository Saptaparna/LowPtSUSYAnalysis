#include "ReadLowPtSUSY_Tree_ControlRegion_CMID.h"

int ReadLowPtSUSY_Tree_ControlRegion_CMID(std::string infile, std::string type, std::string signSelection){
  
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

  double ee_OS = 0;
  double ee_SS = 0;
  double mumu_OS = 0;
  double mumu_SS = 0;

  TFile *trigger1D1=new TFile("TurnOn_SinglePhotonParked_Run2012D_22Jan2013_All.root");
  TFile *trigger1D2=new TFile("TurnOn_METParked_Run2012D_All.root");

  TF1* fit_curve1=(TF1*)trigger1D1->Get("fit_caloMET");
  TF1* fit_curve2=(TF1*)trigger1D2->Get("fit_PHO");

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
       photons.push_back(photon);
      }
     // Now sorting this vector of structs
     std::sort (photons.begin(), photons.end(), sortPhotonsInDescendingpT);
       
     TLorentzVector ph1_p4;
     ph1_p4 = fillTLorentzVector(0.0, 0.0, 0.0, 0.0);

     //here working with the leading photon.
     if (photons.size() > 0) ph1_p4=fillTLorentzVector(photons.at(0).pT, photons.at(0).eta, photons.at(0).phi, photons.at(0).energy); 
     
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
         if(el1_p4.Pt() > 10.0) deltaR_el1 = el1_p4.DeltaR(trigger1_p4);
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
         if(el2_p4.Pt() > 10.0) deltaR_el2 = el2_p4.DeltaR(trigger1_p4);
       }
     }//only execute if the second electron exists.

    if(type=="MC") deltaR_el1 = 100.0; 
    if(type=="MC") deltaR_el2 = 100.0;

    double leadingDeltaR_el, trailingDeltaR_el;
    leadingDeltaR_el = -1.0;
    trailingDeltaR_el = -1.0;

    if(el1_p4.Pt()>10.0 and ph1_p4.Pt()>30.0 and photons.at(0).isTight==1 and photons.at(0).phIsoTight==1){
      leadingDeltaR_el = el1_p4.DeltaR(ph1_p4); 
    }   

    if(el2_p4.Pt()>10.0 and ph1_p4.Pt()>30.0 and photons.at(0).isTight==1 and photons.at(0).phIsoTight==1){
      trailingDeltaR_el = el2_p4.DeltaR(ph1_p4); 
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
     // Now sorting this vector of structs
    std::sort (muons.begin(), muons.end(), sortLeptonsInDescendingpT);
    TLorentzVector mu1_p4;
    TLorentzVector mu2_p4;
    mu1_p4 = fillTLorentzVector(0.0, 0.0, 0.0, 0.0);
    mu2_p4 = fillTLorentzVector(0.0, 0.0, 0.0, 0.0);

   if (muons.size() > 0) mu1_p4=fillTLorentzVector(muons.at(0).pT, muons.at(0).eta, muons.at(0).phi, muons.at(0).energy);
   if (muons.size() > 1) mu2_p4=fillTLorentzVector(muons.at(1).pT, muons.at(1).eta, muons.at(1).phi, muons.at(1).energy);

   double leadingDeltaR_mu, trailingDeltaR_mu;
    leadingDeltaR_mu = -1.0;
    trailingDeltaR_mu = -1.0;

    if(mu1_p4.Pt()>5.0 and ph1_p4.Pt()>30.0 and photons.at(0).isTight==1 and photons.at(0).phIsoTight==1){
      leadingDeltaR_mu = mu1_p4.DeltaR(ph1_p4);
    }

    if(mu2_p4.Pt()>5.0 and ph1_p4.Pt()>30.0 and photons.at(0).isTight==1 and photons.at(0).phIsoTight==1){
      trailingDeltaR_mu = mu2_p4.DeltaR(ph1_p4);
    }

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

  //Preliminary event cuts guided by the trigger:
  if (not(ph1_p4.Pt()>30.0 and photons.at(0).isTight==1 and photons.at(0).phIsoTight==1 and caloMET>25.0)) continue;

  //application of trigger weights
  double eventWeight_Trigger = 0.0;
  double eventWeight = 0.0;
    if(type=="MC")
    {
      double triggerWeight1 = 0.0;
      if(caloMET>25.0 and caloMET < 200.0)triggerWeight1 = fit_curve1->Eval(caloMET);
      else if(caloMET > 200.0) triggerWeight1 = fit_curve1->Eval(200.0);
      double triggerWeight2 = 0.0;
      if(ph1_p4.Pt()>30.0 and ph1_p4.Pt()<100.0 and photons.at(0).isTight==1 and photons.at(0).phIsoTight==1) triggerWeight2 = fit_curve2->Eval(ph1_p4.Pt());
      else if(ph1_p4.Pt()>100.0 and photons.at(0).isTight==1 and photons.at(0).phIsoTight==1) triggerWeight2 = fit_curve2->Eval(100.0);
      eventWeight_Trigger=triggerWeight1*triggerWeight2;
      eventWeight=triggerWeight1*triggerWeight2*PUWeightData;    
    }
    else if(type=="Data")
    {
      eventWeight = 1.0;
    }

  int elel_event = 0;
  int mumu_event = 0;
  int elmu_event = 0;

  int nmuons = 0;
  for(unsigned int j=0; j<muons.size(); ++j)
    {
    if(muons.at(j).pT > 5.0 and muons.at(j).isTight==1 and muons.at(j).isolation < 0.12){
       nmuons++;
      } 
    }

  int nelectrons = 0;
  for(unsigned int j=0; j<electrons.size(); ++j)
    {
    if(electrons.at(j).pT > 10.0 and electrons.at(j).isTight==1 and electrons.at(j).isolation < 0.10){
       nelectrons++;
      }
    }

 
  if(mu1_p4.Pt()>5.0 and muons.at(0).isTight==1 and muons.at(0).isolation < 0.12 and leadingDeltaR_mu > 0.4 and mu2_p4.Pt()>5.0 and muons.at(1).isTight==1 and muons.at(1).isolation < 0.12 and trailingDeltaR_mu > 0.4 and eventWeight > 0.0) {
    if(type=="MC") eventWeight *= muonSF(mu1_p4.Pt(), mu1_p4.Eta())*muonSF(mu2_p4.Pt(), mu2_p4.Eta())*photonSF(ph1_p4.Pt(), ph1_p4.Eta());  
    if(nmuons!=2) continue;  
    mumu_event = 1;
    if (not((mu1_p4+mu2_p4).M() > 60 and (mu1_p4+mu2_p4).M() < 120)) continue;
    if(signSelection=="OS" and (muons.at(0).charge*muons.at(1).charge)==-1){
      mumu_OS+=eventWeight;
    }//OS mode
    if(signSelection=="SS" and (muons.at(0).charge*muons.at(1).charge)==+1){
      mumu_SS+=eventWeight;
    }//SS mode
  }//tight + tight muon

 if(el1_p4.Pt()>10.0 and electrons.at(0).isTight==1 and electrons.at(0).isolation < 0.10 and leadingDeltaR_el > 0.4 and deltaR_el1 > 0.4 and el2_p4.Pt()>10.0 and electrons.at(1).isTight==1 and electrons.at(1).isolation < 0.10 and trailingDeltaR_el > 0.4 and deltaR_el2 > 0.4 and eventWeight > 0.0 and mumu_event==0){
    if(type=="MC") eventWeight *= electronSF(el1_p4.Pt(), el1_p4.Eta())*electronSF(el2_p4.Pt(), el2_p4.Eta())*photonSF(ph1_p4.Pt(), ph1_p4.Eta());
    if (nelectrons!=2) continue; 
    elel_event = 1;
    //if (not((el1_p4+el2_p4).M() > 60 and (el1_p4+el2_p4).M() < 120)) continue;
    if(signSelection=="OS" and (electrons.at(0).charge*electrons.at(1).charge)==-1) {
      ee_OS+=eventWeight;
    }//OS mode
    if(signSelection=="SS" and (electrons.at(0).charge*electrons.at(1).charge)==+1){
      ee_SS+=eventWeight;
    }//SS mode
  }//tight + tight electron
}//event loop closed

  //Cut flow table
  cout << "mumu_OS = " << mumu_OS << endl;
  cout << "mumu_SS = " << mumu_SS << endl;
  cout << "ee_OS = " << ee_OS << endl;
  cout << "ee_SS = " << ee_SS << endl;

}
