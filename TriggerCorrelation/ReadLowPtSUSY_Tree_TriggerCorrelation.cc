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

double mdeltaR(double eta1, double phi1, double eta2, double phi2) { //this function only works if both phi1 and phi2 are defined from -pi to pi or 0 to 2*pi
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

bool sortJetVectorsInDescendingpT(TLorentzVector jet1, TLorentzVector jet2)
{
  return (jet1.Pt() > jet2.Pt());
}

bool sortPhotonsInDescendingpT(PhotonInfo pho1, PhotonInfo pho2)
{
  return (pho1.pT > pho2.pT);
}

int ReadLowPtSUSY_Tree_TriggerCorrelation(std::string infile, std::string outfile){

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

  TH2F *h_ControlTrigger1 = new TH2F("h_ControlTrigger1", "Scatter Plot of Photon pT versus MET; Photon pT [GeV]; MET [GeV]", 1000, 0, 1000.0, 1000.0, 0, 1000.0);  h_ControlTrigger1->Sumw2(); 
  TH2F *h_AnalysisTrigger1 = new TH2F("h_AnalysisTrigger1", "Scatter Plot of Photon pT versus MET; Photon pT [GeV]; MET [GeV]", 1000, 0, 1000.0, 1000.0, 0, 1000.0);  h_AnalysisTrigger1->Sumw2();
  TH2F *h_ControlTrigger2 = new TH2F("h_ControlTrigger2", "Scatter Plot of Photon pT versus MET; Photon pT [GeV]; MET [GeV]", 1000, 0, 1000.0, 1000.0, 0, 1000.0);  h_ControlTrigger2->Sumw2();
  TH2F *h_AnalysisTrigger2 = new TH2F("h_AnalysisTrigger2", "Scatter Plot of Photon pT versus MET; Photon pT [GeV]; MET [GeV]", 1000, 0, 1000.0, 1000.0, 0, 1000.0);  h_AnalysisTrigger2->Sumw2();

  
  int nEvents=tree->GetEntries();
  std::cout << "nEvents= " << nEvents << std::endl;
  for (int i=0; i<nEvents; ++i)
    {
     tree->GetEvent(i);
 
     TriggerInfo trigger1;//Corresponds to: hltPhoton30R9Id90CaloIdHE10Iso40EBOnlyTrackIsoLastFilter
     //fired for both
     //HLT_Photon30_R9Id90_CaloId_HE10_Iso40_EBOnly_v* (Prescaled by 20)
     //HLT_Photon30_R9Id90_CaloId_HE10_Iso40_EBOnly_Met25_HBHENoiseCleaned_v* //main analysis trigger 
     trigger1.Px = trigObj1Px;
     trigger1.Py = trigObj1Py;
     trigger1.Pz = trigObj1Pz;
     trigger1.E  = trigObj1E;

     TriggerInfo trigger2; //Corresponds to: 
     //fired for hltPhoton30HEFilter
     //HLT_Photon30_v*  
     trigger2.Px = trigObj2Px;
     trigger2.Py = trigObj2Py;
     trigger2.Pz = trigObj2Pz;
     trigger2.E  = trigObj2E;

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

   //To understand this statement better: fired_HLTPho corresponds to HLT_Photon30_v* (Prescaled by 500)
   //fired_HLTPhoId correspond to  HLT_Photon30_R9Id90_CaloId_HE10_Iso40_EBOnly_v* (Prescaled by 20)
   //fired_HLTPhoIdMet corresponds to HLT_Photon30_R9Id90_CaloId_HE10_Iso40_EBOnly_Met25_HBHENoiseCleaned_v*

   if(photons.size() > 0){
     if(fired_HLTPho==1 and foundHLTPhoton2==1 and photons.at(0).pT>30.0 and photons.at(0).isTight==1 and photons.at(0).phIsoTight==1){
       h_ControlTrigger1->Fill(photons.at(0).pT, MET);
       if(foundHLTPhoton1==1 and fired_HLTPhoId==1){
         h_AnalysisTrigger1->Fill(photons.at(0).pT, MET);
       }//main trigger if   
     }//control trigger photon 30 if
   
   if (fired_HLTPhoId==1 and photons.at(0).pT>30.0 and photons.at(0).isTight==1 and photons.at(0).phIsoTight==1){
     h_ControlTrigger2->Fill(photons.at(0).pT, MET);
     if(fired_HLTPhoIdMet==1){
        h_AnalysisTrigger2->Fill(photons.at(0).pT, MET);
       }//main trigger if 
     }//control trigger photon 30 if
   }//photon check

  }//end of event loop
  std::string histfilename=(outfile+".root").c_str();
  TFile *tFile=new TFile(histfilename.c_str(), "RECREATE");
  
  h_ControlTrigger1->Write();
  h_AnalysisTrigger1->Write();
  h_ControlTrigger2->Write();
  h_AnalysisTrigger2->Write();
  tFile->Close();
  std::cout<<"Wrote output file "<<histfilename<<std::endl;

  return 0;

}
