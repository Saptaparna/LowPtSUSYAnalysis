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
   outtree->Branch("ph_pt_leading", &ph_pt_leading, "ph_pt_leading/F");
   outtree->Branch("ph_energy_leading", &ph_energy_leading, "ph_energy_leading/F");
   outtree->Branch("ph_eta_leading", &ph_eta_leading, "ph_eta_leading/F");
   outtree->Branch("ph_phi_leading", &ph_phi_leading, "ph_phi_leading/F");
   outtree->Branch("nPhotons", &nPhotons, "nPhotons/I");


   run = 0;
   lumi = 0;
   event = 0;

   ph_pt_leading=-99.0;
   ph_energy_leading=-99.0;
   ph_eta_leading=-99.0;
   ph_phi_leading=-99.0;

   nPhotons = -1;


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
  
  //Accessing Photons
  
  vector<TCPhoton> vPhotons;

  for (Int_t i = 0; i < recoPhotons->GetSize(); ++i) {
    TCPhoton* photon = (TCPhoton*) recoPhotons->At(i);
    vPhotons.push_back(*photon);
  }  

  nPhotons = vPhotons.size();

  if(vPhotons.size() > 1){

    ph_pt_leading = vPhotons[0].Pt();
    ph_eta_leading = vPhotons[0].Eta();
    ph_phi_leading = vPhotons[0].Phi();
    ph_energy_leading = vPhotons[0].Energy();
  }

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

