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
  bool isLoose;
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
  bool isMedium;
  bool isLoose;
  float chIsolation;
  float nuIsolation;
  float phIsolation;
  bool  phIsoSTight;
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
    

int ReadLowPtSUSY_Tree_WithMETTrigger_FixingNuances(std::string infile, std::string outfile){

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
  vector<int>     *jet_mva_loose;
  vector<int>     *jet_mva_tight;
  vector<int>     *jet_mva_medium;
  vector<int>     *jet_cut_loose;
  vector<int>     *jet_cut_tight;
  vector<int>     *jet_cut_medium;
  Bool_t          fired_HLTPho;
  Bool_t          fired_HLTPhoId;
  Bool_t          fired_HLTPhoIdMet;
  Bool_t          fired_HLTMET100;
  vector<float>   *ph_pt;
  vector<float>   *ph_phi;
  vector<float>   *ph_eta;
  vector<float>   *ph_energy;
  vector<bool>    *ph_phIsoSTight;
  vector<bool>    *ph_phIsoTight;
  vector<bool>    *ph_phIsoMedium;
  vector<bool>    *ph_phIsoLoose;
  vector<float>   *ph_chIso;
  vector<float>   *ph_nuIso;
  vector<float>   *ph_phIso;
  vector<bool>    *ph_isTight;
  vector<bool>    *ph_isMedium;
  vector<bool>    *ph_isLoose;
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
  vector<float>   *el_pt;
  vector<float>   *el_phi;
  vector<float>   *el_eta;
  vector<float>   *el_energy;
  vector<int>     *el_charge;
  vector<bool>    *el_isTight;
  vector<bool>    *el_isLoose;
  vector<float>   *el_iso;
  Float_t         MET;
  float trigObj1Px;
  float trigObj1Py;
  float trigObj1Pz;
  float trigObj1E;
  float trigObj2Px;
  float trigObj2Py;
  float trigObj2Pz;
  float trigObj2E;

  mu_pt = 0;
  mu_phi = 0;
  mu_eta = 0;
  mu_energy = 0;
  mu_charge = 0;
  mu_isTight = 0;
  mu_iso = 0;
  el_pt = 0;
  el_phi = 0;
  el_eta = 0;
  el_energy = 0;
  el_charge = 0;
  el_isTight = 0;
  el_isLoose = 0;
  el_iso = 0;
  ph_pt = 0;
  ph_phi = 0;
  ph_eta = 0;
  ph_energy = 0;
  ph_chIso = 0;
  ph_nuIso = 0;
  ph_phIso = 0;
  ph_isTight = 0;
  ph_isMedium = 0;
  ph_isLoose = 0;
  ph_phIsoSTight = 0;
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

  tree->SetBranchAddress("el_pt", &(el_pt));
  tree->SetBranchAddress("el_eta", &(el_eta));
  tree->SetBranchAddress("el_phi", &(el_phi));
  tree->SetBranchAddress("el_energy", &(el_energy));
  tree->SetBranchAddress("el_charge", &(el_charge));
  tree->SetBranchAddress("el_isTight", &(el_isTight));
  tree->SetBranchAddress("el_isLoose", &(el_isLoose));
  tree->SetBranchAddress("el_iso", &(el_iso));
  tree->SetBranchAddress("nVertices", &(nVertices));
  tree->SetBranchAddress("fired_HLTPho", &(fired_HLTPho));
  tree->SetBranchAddress("fired_HLTPhoId", &(fired_HLTPhoId));
  tree->SetBranchAddress("fired_HLTPhoIdMet", &(fired_HLTPhoIdMet));
  tree->SetBranchAddress("fired_HLTPhoIdMet", &(fired_HLTPhoIdMet));
  tree->SetBranchAddress("fired_HLTMET100", &(fired_HLTMET100));
  tree->SetBranchAddress("ph_pt", &(ph_pt));
  tree->SetBranchAddress("ph_phi", &(ph_phi));
  tree->SetBranchAddress("ph_eta", &(ph_eta));
  tree->SetBranchAddress("ph_energy", &(ph_energy));
  tree->SetBranchAddress("ph_chIso", &(ph_chIso));
  tree->SetBranchAddress("ph_nuIso", &(ph_nuIso));
  tree->SetBranchAddress("ph_phIso", &(ph_phIso));
  tree->SetBranchAddress("ph_isTight", &(ph_isTight));
  tree->SetBranchAddress("ph_isMedium", &(ph_isMedium));
  tree->SetBranchAddress("ph_isLoose", &(ph_isLoose));
  tree->SetBranchAddress("ph_phIsoSTight", &(ph_phIsoSTight));
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
  tree->SetBranchAddress("MET", &(MET));
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

  TH1F *h_PHO_HLT_METTrig=new TH1F("h_PHO_HLT_METTrig", "Photon pT; photon pT [GeV]; Events/GeV", 600, 0, 600); h_PHO_HLT_METTrig->Sumw2();
  TH1F *h_PHO_HLTPho_METTrig=new TH1F("h_PHO_HLTPho_METTrig", "Photon pT; photon pT [GeV]; Events/GeV", 600, 0, 600); h_PHO_HLTPho_METTrig->Sumw2();  
  TH1F *h_PHO_HLT_METTrig_RemoveJets=new TH1F("h_PHO_HLT_METTrig_RemoveJets", "Photon pT; photon pT [GeV]; Events/GeV", 600, 0, 600); h_PHO_HLT_METTrig_RemoveJets->Sumw2();
  TH1F *h_PHO_HLTPho_METTrig_RemoveJets=new TH1F("h_PHO_HLTPho_METTrig_RemoveJets", "Photon pT; photon pT [GeV]; Events/GeV", 600, 0, 600); h_PHO_HLTPho_METTrig_RemoveJets->Sumw2();
  TH1F *h_DeltaR_Ph_Jet = new TH1F("h_DeltaR_Ph_Jet", "#Delta R between the photon and the jet; #Delta R; Events", 4000, 0, 4.0);h_DeltaR_Ph_Jet->Sumw2();
  TH1F *h_PHO_HLT_METTrig_ST=new TH1F("h_PHO_HLT_METTrig_ST", "Photon pT; photon pT [GeV]; Events/GeV", 600, 0, 600); h_PHO_HLT_METTrig_ST->Sumw2();
  TH1F *h_PHO_HLTPho_METTrig_ST=new TH1F("h_PHO_HLTPho_METTrig_ST", "Photon pT; photon pT [GeV]; Events/GeV", 600, 0, 600); h_PHO_HLTPho_METTrig_ST->Sumw2();
  TH1F *h_PHO_HLT_METTrig_M=new TH1F("h_PHO_HLT_METTrig_M", "Photon pT; photon pT [GeV]; Events/GeV", 600, 0, 600); h_PHO_HLT_METTrig_M->Sumw2();
  TH1F *h_PHO_HLTPho_METTrig_M=new TH1F("h_PHO_HLTPho_METTrig_M", "Photon pT; photon pT [GeV]; Events/GeV", 600, 0, 600); h_PHO_HLTPho_METTrig_M->Sumw2();
  TH1F *h_PHO_HLT_METTrig_L=new TH1F("h_PHO_HLT_METTrig_L", "Photon pT; photon pT [GeV]; Events/GeV", 600, 0, 600); h_PHO_HLT_METTrig_L->Sumw2();
  TH1F *h_PHO_HLTPho_METTrig_L=new TH1F("h_PHO_HLTPho_METTrig_L", "Photon pT; photon pT [GeV]; Events/GeV", 600, 0, 600); h_PHO_HLTPho_METTrig_L->Sumw2();

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

   vector<TLorentzVector> Photon_vector;
   Photon_vector.clear();
   std::vector<float> swissCross;
   swissCross.clear();
   for(unsigned int k=0; k<photons.size(); ++k)
   {
     TLorentzVector Photon;
     swissCross.push_back((photons.at(k).phseedCrystalEnergy)/(photons.at(k).phe1x5 + photons.at(k).phe5x1 - photons.at(k).phseedCrystalEnergy));

    if(photons.at(k).pT>0.0 and photons.at(k).isTight==1 and photons.at(k).phIsoTight==1 and photons.at(k).phpixelVeto==0 and photons.at(k).phSigmaIetaIeta > 0.001 and swissCross.at(k)<0.90 and photons.at(k).phR9 < 1.0 and photons.at(k).phSigmaIphiIphi>0.0)////0.0001
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
       for(unsigned int j=0; j<muons.size(); ++j)
         {
         TLorentzVector Muon;
         if(muons.at(j).pT > 5.0 and muons.at(j).isTight==1 and muons.at(j).isolation < 0.12){
         Muon.SetPtEtaPhiE(muons.at(j).pT, muons.at(j).eta, muons.at(j).phi, muons.at(j).energy);
         double DRph_mu = Photon.DeltaR(Muon);
        if(DRph_mu<0.5) isGoodPhoton=false;
        }
      }
      if(isGoodPhoton) Photon_vector.push_back(Photon);
    }//close four vector if
  }//close photon loop

// Now sorting this vector of structs
  std::sort (Photon_vector.begin(), Photon_vector.end(), sortPhotonVectorsInDescendingpT);
//  cout << "Photon_vector.size() = " << Photon_vector.size() << endl; 
 
  std::sort (jets.begin(), jets.end(), sortJetsInDescendingpT);

   double minDR_PhJ = 9999;
   int foundJetNearPhoton=0;
   if(Photon_vector.size() > 0)
   {
     if(Photon_vector.at(0).Pt()>0.0)
     {
       TLorentzVector Photon;
       Photon.SetPtEtaPhiE(Photon_vector.at(0).Pt(), Photon_vector.at(0).Eta(), Photon_vector.at(0).Phi(), Photon_vector.at(0).E());
       for(unsigned int k=0; k<jets.size(); ++k)
       {
         TLorentzVector Jet;
         if(fabs(jets.at(k).eta)<2.4 and jets.at(k).pT>25.0 and jets.at(k).PU_mva_loose==1)
         {
           Jet.SetPtEtaPhiE(jets.at(k).pT, jets.at(k).eta, jets.at(k).phi, jets.at(k).energy);
           double dRjet_ph=Jet.DeltaR(Photon);
           if (dRjet_ph<minDR_PhJ) minDR_PhJ=dRjet_ph;
         }
       }
       h_DeltaR_Ph_Jet->Fill(minDR_PhJ);
       if (minDR_PhJ<0.5) foundJetNearPhoton=1;
     }
   }

 
   if(Photon_vector.size() > 0){
    if(fired_HLTMET100==1){
       h_PHO_HLT_METTrig->Fill(Photon_vector.at(0).Pt());
       if(fired_HLTPhoIdMet==1){
         h_PHO_HLTPho_METTrig->Fill(Photon_vector.at(0).Pt());
       }
     }
   }
   
 if(Photon_vector.size() > 0){
     if(fired_HLTMET100==1 and foundJetNearPhoton==0){
       h_PHO_HLT_METTrig_RemoveJets->Fill(Photon_vector.at(0).Pt());
       if(fired_HLTPhoIdMet==1){
         h_PHO_HLTPho_METTrig_RemoveJets->Fill(Photon_vector.at(0).Pt());
       }
     }
   }

  }//end of event loop
  std::string histfilename=(outfile+".root").c_str();
  TFile *tFile=new TFile(histfilename.c_str(), "RECREATE");
  h_DeltaR_Ph_Jet->Write(); 
  h_PHO_HLTPho_METTrig->Write();
  h_PHO_HLT_METTrig->Write();
  h_PHO_HLT_METTrig_RemoveJets->Write();
  h_PHO_HLTPho_METTrig_RemoveJets->Write();
  tFile->Close();
  std::cout<<"Wrote output file "<<histfilename<<std::endl;

  return 0;

}
