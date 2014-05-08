#define FlatTreeCreator_cxx
/* Code for accessing variables from a non-linear 
Tree and storing relevant information
Author: Saptaparna Bhattcharya

*/


#include "FlatTreeCreator.h"
#include <TH2.h>
#include <TStyle.h>
using namespace std;

void FlatTreeCreator::Begin(TTree *tree)
{

   TString option = GetOption();
   
   outputFile = new TFile("LowPtSUSY_Tree.root","RECREATE");
   outtree=new TTree("LowPtSUSY_Tree", "LowPtSUSY_Tree"); 
   outtree->Branch("run", &run, "run/I");
   outtree->Branch("lumi", &lumi, "lumi/I");
   outtree->Branch("event", &event, "event/I");
   outtree->Branch("ph_pt", &ph_pt);
   outtree->Branch("ph_phi", &ph_phi);
   outtree->Branch("ph_eta", &ph_eta);
   outtree->Branch("ph_energy", &ph_energy);
   outtree->Branch("nPhotons", &nPhotons, "nPhotons/I");
   outtree->Branch("el_pt", &el_pt);
   outtree->Branch("el_phi", &el_phi);
   outtree->Branch("el_eta", &el_eta);
   outtree->Branch("el_energy", &el_energy);
   outtree->Branch("nElectrons", &nElectrons, "nElectrons/I");
   outtree->Branch("mu_pt", &mu_pt);
   outtree->Branch("mu_phi", &mu_phi);
   outtree->Branch("mu_eta", &mu_eta);
   outtree->Branch("mu_energy", &mu_energy);
   outtree->Branch("nMuons", &nMuons, "nMuons/I");
   outtree->Branch("jet_pt", &jet_pt);
   outtree->Branch("jet_phi", &jet_phi);
   outtree->Branch("jet_eta", &jet_eta);
   outtree->Branch("jet_energy", &jet_energy);
   outtree->Branch("nJets", &nJets, "nJets/I");

   run = 0;
   lumi = 0;
   event = 0;

   ph_pt.clear();
   ph_phi.clear();
   ph_eta.clear();
   ph_energy.clear();
   nPhotons = -1;
   el_pt.clear();
   el_phi.clear();
   el_eta.clear();
   el_energy.clear();
   nElectrons = -1;
   mu_pt.clear();
   mu_phi.clear();
   mu_eta.clear();
   mu_energy.clear(); 
   nMuons = -1;
   jet_pt.clear();
   jet_phi.clear();
   jet_eta.clear();
   jet_energy.clear(); 
   nJets = -1;
   MET = 0;
}

void FlatTreeCreator::SlaveBegin(TTree *)
{
   TString option = GetOption();
}

Bool_t FlatTreeCreator::Process(Long64_t entry)
{
  GetEntry(entry);
  if(entry % 1000 == 0) cout << "Processing event number: " << entry << endl;
  
  run = runNumber;
  lumi = lumiSection;
  event = eventNumber;
  
  vector<TCPhoton> vPhotons;

  for (Int_t i = 0; i < recoPhotons->GetSize(); ++i) {
    TCPhoton* photon = (TCPhoton*) recoPhotons->At(i);
    vPhotons.push_back(*photon);
  }  

  nPhotons = vPhotons.size();

  for (int i=0; i < nPhotons; i++){
    ph_pt.push_back(vPhotons[i].Pt());
    ph_eta.push_back(vPhotons[i].Eta());
    ph_phi.push_back(vPhotons[i].Phi());
    ph_energy.push_back(vPhotons[i].Energy());
  }
  
 vector<TCElectron> vElectrons;

 for (Int_t i = 0; i < recoElectrons->GetSize(); ++i) {
   TCElectron* electron = (TCElectron*) recoElectrons->At(i);
   vElectrons.push_back(*electron);
  }

 nElectrons = vElectrons.size();

 for (int i=0; i < nElectrons; i++){
   el_pt.push_back(vElectrons[i].Pt());
   el_eta.push_back(vElectrons[i].Eta());
   el_phi.push_back(vElectrons[i].Phi());
   el_energy.push_back(vElectrons[i].Energy());
  }

 vector<TCMuon> vMuons;

 for (Int_t i = 0; i < recoMuons->GetSize(); ++i) {
   TCMuon* muon = (TCMuon*) recoMuons->At(i);
   vMuons.push_back(*muon);
  }

 nMuons = vMuons.size();

 for (int i=0; i < nMuons; i++){
   mu_pt.push_back(vMuons[i].Pt());
   mu_eta.push_back(vMuons[i].Eta());
   mu_phi.push_back(vMuons[i].Phi());
   mu_energy.push_back(vMuons[i].Energy());
  }
 

 vector<TCJet> vJets;

 for (Int_t i = 0; i < patJets->GetSize(); ++i) {
   TCJet* jet = (TCJet*) patJets->At(i);
   vJets.push_back(*jet);
  }

 nJets = vJets.size();

 for (int i=0; i < nJets; i++){
   jet_pt.push_back(vJets[i].Pt());
   jet_eta.push_back(vJets[i].Eta());
   jet_phi.push_back(vJets[i].Phi());
   jet_energy.push_back(vJets[i].Energy());
  }


 MET = corrMET->Mod();

 outtree->Fill();
 return kTRUE;
}

void FlatTreeCreator::SlaveTerminate()
{
}

void FlatTreeCreator::Terminate()
{
 outtree->Write();
}

