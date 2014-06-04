#define FlatTreeCreator_cxx
/* Code for accessing variables from a non-linear 
Tree and storing relevant information
Author: Saptaparna Bhattcharya

*/


#include "FlatTreeCreator.h"
#include <TH2.h>
#include <TStyle.h>
#include "PhysicsTools/Utilities/interface/LumiReweightingStandAlone.h"
using namespace std;

reweight::LumiReWeighting LumiWeightsD_;
reweight::LumiReWeighting LumiWeightsD_sys_;
bool isMC= true;
//bool isMC=false;
string  suffix = "SUFFIX";

void FlatTreeCreator::Begin(TTree *tree)
{

   TString option = GetOption();

   //PileUp Reweighting:

   Float_t DataDist_2012D[60] = {1768.77, 3651.68, 7356.92, 14793.3, 51246.8, 489173, 2.63921e+06, 6.82361e+06, 1.88317e+07, 5.19794e+07, 1.08899e+08, 1.88257e+08, 2.57316e+08, 3.01258e+08, 3.27492e+08, 3.44354e+08, 3.59374e+08, 3.74823e+08, 3.90058e+08, 4.0217e+08, 4.08643e+08, 4.11422e+08, 4.11151e+08, 4.05573e+08, 3.92856e+08, 3.72226e+08, 3.44284e+08, 3.10504e+08, 2.7226e+08, 2.30759e+08, 1.8777e+08, 1.45857e+08, 1.07763e+08, 7.56393e+07, 5.0551e+07, 3.23928e+07, 2.0178e+07, 1.25011e+07, 7.95461e+06, 5.38133e+06, 3.95572e+06, 3.15182e+06, 2.66553e+06, 2.33493e+06, 2.08e+06, 1.86352e+06, 1.669e+06, 1.48952e+06, 1.32229e+06, 1.16642e+06, 1.02175e+06, 888282, 766092, 655189, 555473, 466702, 388487, 320302, 261506, 211367};

   Float_t DataDist_2012D_sys[60] = {1632.3, 3240.95, 6296.46, 12068, 30517.4, 227305, 1.49497e+06, 4.52239e+06, 1.0476e+07, 2.96274e+07, 6.82538e+07, 1.28218e+08, 2.00114e+08, 2.55292e+08, 2.90327e+08, 3.11739e+08, 3.26279e+08, 3.39668e+08, 3.53444e+08, 3.67108e+08, 3.78444e+08, 3.85013e+08, 3.88036e+08, 3.88742e+08, 3.85622e+08, 3.76942e+08, 3.61739e+08, 3.39978e+08, 3.12601e+08, 2.80828e+08, 2.45645e+08, 2.08081e+08, 1.69728e+08, 1.32723e+08, 9.92184e+07, 7.08564e+07, 4.84379e+07, 3.18819e+07, 2.04316e+07, 1.29836e+07, 8.39661e+06, 5.69316e+06, 4.14131e+06, 3.24848e+06, 2.71199e+06, 2.35979e+06, 2.10068e+06, 1.88918e+06, 1.70364e+06, 1.53423e+06, 1.37666e+06, 1.22917e+06, 1.0912e+06, 962650, 843537, 733913, 633795, 543116, 461705, 389280};
  
   Float_t MCDist_Summer2012_S10[60] = {2.560E-06,5.239E-06,1.420E-05,5.005E-05,1.001E-04,2.705E-04,1.999E-03,6.097E-03,1.046E-02,1.383E-02,1.685E-02,2.055E-02,2.572E-02,3.262E-02,4.121E-02,4.977E-02,5.539E-02,5.725E-02,5.607E-02,5.312E-02,5.008E-02,4.763E-02,4.558E-02,4.363E-02,4.159E-02,3.933E-02,3.681E-02,3.406E-02,3.116E-02,2.818E-02,2.519E-02,2.226E-02,1.946E-02,1.682E-02,1.437E-02,1.215E-02,1.016E-02,8.400E-03,6.873E-03,5.564E-03,4.457E-03,3.533E-03,2.772E-03,2.154E-03,1.656E-03,1.261E-03,9.513E-04,7.107E-04,5.259E-04,3.856E-04,2.801E-04,2.017E-04,1.439E-04,1.017E-04,7.126E-05,4.948E-05,3.405E-05,2.322E-05,1.570E-05,5.005E-06};

   std::vector<float> DataDistD;
   std::vector<float> DataDistD_sys;
   std::vector<float> MCDist;
  
   for( int i=0; i<60; i++) {
     DataDistD.push_back(DataDist_2012D[i]);
     DataDistD_sys.push_back(DataDist_2012D_sys[i]);
     MCDist.push_back(MCDist_Summer2012_S10[i]);
   }

   LumiWeightsD_ = reweight::LumiReWeighting(MCDist, DataDistD);
   LumiWeightsD_sys_ = reweight::LumiReWeighting(MCDist, DataDistD_sys);

   fileCount = 0; 
   outputFile = new TFile("LowPtSUSY_Tree_SUFFIX.root","RECREATE");
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
   outtree->Branch("el_charge", &el_charge);
   outtree->Branch("el_isTight", &el_isTight);
   outtree->Branch("nElectrons", &nElectrons, "nElectrons/I");
   outtree->Branch("mu_pt", &mu_pt);
   outtree->Branch("mu_phi", &mu_phi);
   outtree->Branch("mu_eta", &mu_eta);
   outtree->Branch("mu_energy", &mu_energy);
   outtree->Branch("mu_charge", &mu_charge);
   outtree->Branch("mu_isTight", &mu_isTight);
   outtree->Branch("nMuons", &nMuons, "nMuons/I");
   outtree->Branch("jet_pt", &jet_pt);
   outtree->Branch("jet_phi", &jet_phi);
   outtree->Branch("jet_eta", &jet_eta);
   outtree->Branch("jet_energy", &jet_energy);
   outtree->Branch("nJets", &nJets, "nJets/I");
   outtree->Branch("jet_mva_loose", &jet_mva_loose);
   outtree->Branch("jet_mva_tight", &jet_mva_tight);
   outtree->Branch("jet_mva_medium", &jet_mva_medium);
   outtree->Branch("jet_cut_loose", &jet_cut_loose);
   outtree->Branch("jet_cut_tight", &jet_cut_tight);
   outtree->Branch("jet_cut_medium", &jet_cut_medium);
   outtree->Branch("MET", &MET, "MET/F");
   outtree->Branch("MET_Phi", &MET_Phi, "MET_Phi/F");
   outtree->Branch("MET_Px", &MET_Px, "MET_Px/F");
   outtree->Branch("MET_Py", &MET_Py, "MET_Py/F");
   outtree->Branch("fired_HLTPho", &fired_HLTPho, "fired_HLTPho/O");
   outtree->Branch("fired_HLTPhoId", &fired_HLTPhoId, "fired_HLTPhoId/O");
   outtree->Branch("fired_HLTPhoIdMet", &fired_HLTPhoIdMet, "fired_HLTPhoIdMet/O");
   outtree->Branch("nVertices", &nVertices, "nVertices/I");
   outtree->Branch("el_iso", &el_iso);
   outtree->Branch("mu_iso", &mu_iso);
   outtree->Branch("ph_chIso", &ph_chIso);
   outtree->Branch("ph_nuIso", &ph_nuIso);
   outtree->Branch("ph_phIso", &ph_phIso);
   outtree->Branch("ph_isTight", &ph_isTight);
   outtree->Branch("ph_phIsoTight", &ph_phIsoTight);
   outtree->Branch("ph_phIsoMedium", &ph_phIsoMedium);
   outtree->Branch("ph_phIsoLoose", &ph_phIsoLoose);
   //Trigger Objects saved
   outtree->Branch("trigObj1Px", &trigObj1Px, "trigObj1Px/F");
   outtree->Branch("trigObj1Py", &trigObj1Py, "trigObj1Py/F");
   outtree->Branch("trigObj1Pz", &trigObj1Pz, "trigObj1Pz/F");
   outtree->Branch("trigObj1E",  &trigObj1E, "trigObj1E/F");
   outtree->Branch("trigObj2Px", &trigObj2Px, "trigObj2Px/F");
   outtree->Branch("trigObj2Py", &trigObj2Py, "trigObj2Py/F");
   outtree->Branch("trigObj2Pz", &trigObj2Pz, "trigObj2Pz/F");
   outtree->Branch("trigObj2E",  &trigObj2E, "trigObj2E/F");

   //MC PU-related variables
   outtree->Branch("nPUVertices", &nPUVertices, "nPUVertices/F");
   outtree->Branch("nPUVerticesTrue", &nPUVerticesTrue, "nPUVerticesTrue/F");
   outtree->Branch("PUWeightData", &PUWeightData, "PUWeightData/F");
   outtree->Branch("PUWeightDataSys", &PUWeightDataSys, "PUWeightDataSys/F");
   //MC GenInfo Related Variables
   outtree->Branch("el_Matched",  &el_Matched); 
   outtree->Branch("el_MatchedPt",  &el_MatchedPt);  
   outtree->Branch("el_MatchedEta",  &el_MatchedEta);  
   outtree->Branch("el_MatchedPhi",  &el_MatchedPhi);
   outtree->Branch("el_MatchedEnergy",  &el_MatchedEnergy);

   outtree->Branch("mu_Matched",  &mu_Matched);
   outtree->Branch("mu_MatchedPt",  &mu_MatchedPt);
   outtree->Branch("mu_MatchedEta",  &mu_MatchedEta);
   outtree->Branch("mu_MatchedPhi",  &mu_MatchedPhi);
   outtree->Branch("mu_MatchedEnergy",  &mu_MatchedEnergy);

   outtree->Branch("Ztt", &Ztt, "Ztt/O");
}

void FlatTreeCreator::SlaveBegin(TTree *)
{
   TString option = GetOption();
}

Bool_t FlatTreeCreator::Process(Long64_t entry)
{
  GetEntry(entry);
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
  el_charge.clear();
  el_isTight.clear();
  nElectrons = -1;
  mu_pt.clear();
  mu_phi.clear();
  mu_eta.clear();
  mu_energy.clear();
  mu_charge.clear();
  mu_isTight.clear();
  nMuons = -1;
  jet_pt.clear();
  jet_phi.clear();
  jet_eta.clear();
  jet_energy.clear();
  nJets = -1;
  MET = 0;
  MET_Phi = -99.0;
  MET_Px = -99.0;
  MET_Py = -99.0;
  fired_HLTPho = false;
  fired_HLTPhoId = false;
  fired_HLTPhoIdMet = false; 
  nVertices = -1;
  el_iso.clear();
  mu_iso.clear();
  ph_chIso.clear();
  ph_nuIso.clear();
  ph_phIso.clear();
  ph_isTight.clear();
  ph_phIsoTight.clear();
  ph_phIsoMedium.clear();
  ph_phIsoLoose.clear();
  PUWeightData = -1.0;
  PUWeightDataSys = -1.0;

  jet_mva_loose.clear();
  jet_mva_tight.clear();
  jet_mva_medium.clear();
  jet_cut_loose.clear();
  jet_cut_tight.clear();
  jet_cut_medium.clear();

  //Trigger object 1 refers to objects passed by the hltPhoton30R9Id90CaloIdHE10Iso40EBOnlyTrackIsoLastFilter
  trigObj1Px= -99.0;
  trigObj1Py= -99.0;
  trigObj1Pz= -99.0;
  trigObj1E= -99.0;
  //Trigger object 1 refers to objects passed by the hltPhoton30HEFilter
  trigObj2Px= -99.0;
  trigObj2Py= -99.0;
  trigObj2Pz= -99.0;
  trigObj2E= -99.0;
  //MC history and genparticle info
  el_Matched.clear();
  el_MatchedPt.clear();
  el_MatchedEta.clear();
  el_MatchedPhi.clear();
  el_MatchedEnergy.clear();

  mu_Matched.clear();
  mu_MatchedPt.clear();
  mu_MatchedEta.clear();
  mu_MatchedPhi.clear();
  mu_MatchedEnergy.clear();
  Ztt = false;

  if(entry % 1000 == 0) cout << "Processing event number: " << entry << endl;

  
 run = runNumber;
 lumi = lumiSection;
 event = eventNumber;

 //Trigger Information
 if(not isMC){
   for (int i = 0; i <  triggerObjects->GetSize(); i++) {
   TCTriggerObject* thisTrigObj = (TCTriggerObject*) triggerObjects->At(i);
   if( thisTrigObj->GetModuleName() == "hltPhoton30R9Id90CaloIdHE10Iso40EBOnlyTrackIsoLastFilter"){
     trigObj1Px = thisTrigObj->Px();
     trigObj1Py = thisTrigObj->Py();
     trigObj1Pz = thisTrigObj->Pz();
     trigObj1E  = thisTrigObj->E();
     }

   if(thisTrigObj->GetModuleName() == "hltPhoton30HEFilter"){
     trigObj2Px = thisTrigObj->Px();
     trigObj2Py = thisTrigObj->Py();
     trigObj2Pz = thisTrigObj->Pz();
     trigObj2E  = thisTrigObj->E();
    }

   //cout << "thisTrigObj->GetModuleName() == " << thisTrigObj->GetModuleName() << endl; 

   if(thisTrigObj->GetHLTName().find("HLT_Photon30_v")!=std::string::npos) fired_HLTPho = true;
   if(thisTrigObj->GetHLTName().find("HLT_Photon30_R9Id90_CaloId_HE10_Iso40_EBOnly_v")!=std::string::npos) fired_HLTPhoId = true;
   if(thisTrigObj->GetHLTName().find("HLT_Photon30_R9Id90_CaloId_HE10_Iso40_EBOnly_Met25_HBHENoiseCleaned_v")!=std::string::npos) fired_HLTPhoIdMet = true;
   }
   if (NoiseFilters_isScraping) return kTRUE;
   if (NoiseFilters_isNoiseHcalHBHE) return kTRUE;
   if (NoiseFilters_isNoiseHcalLaser) return kTRUE;
   if (NoiseFilters_isNoiseEcalTP) return kTRUE;
   if (NoiseFilters_isNoiseEcalBE) return kTRUE;
   if (NoiseFilters_isCSCTightHalo) return kTRUE;
   if (NoiseFilters_isNoiseEEBadSc) return kTRUE;
   if (!NoiseFilters_isNoisetrkPOG1) return kTRUE;
   if (!NoiseFilters_isNoisetrkPOG2) return kTRUE;
   if (!NoiseFilters_isNoisetrkPOG3) return kTRUE; 

 }

  vector<TVector3> goodVertices;
  for (int i = 0; i < primaryVtx->GetSize(); i++) {
    TCPrimaryVtx* pVtx = (TCPrimaryVtx*) primaryVtx->At(i);
    if ((not pVtx->IsFake()) && pVtx->NDof() > 4. && fabs(pVtx->z()) <= 24. && fabs(pVtx->Perp()) <= 2.) {
      goodVertices.push_back(*pVtx);
    }
  }

  //Reject events that do not have a good vertex 

  if (goodVertices.size() < 1) return kTRUE;

  nVertices = goodVertices.size();
  pvPosition = new TVector3();
  *pvPosition = goodVertices[0]; 
 
  vector<TCPhoton> vPhotons;

  for (Int_t i = 0; i < recoPhotons->GetSize(); i++) {
    TCPhoton* photon = (TCPhoton*) recoPhotons->At(i);
    if (!(fabs(photon->SCEta())<1.4442)) continue;
    if (!(photon->R9()>0.9)) continue;
    vPhotons.push_back(*photon);
    ph_pt.push_back(photon->Pt());
    ph_eta.push_back(photon->Eta());
    ph_phi.push_back(photon->Phi());
    ph_energy.push_back(photon->Energy());
    double chIsoTemp, nuIsoTemp, phIsoTemp;
    bool isoPassLTemp, isoPassMTemp, isoPassTTemp;
    PhotonIso(photon, chIsoTemp, nuIsoTemp, phIsoTemp, isoPassLTemp, isoPassMTemp, isoPassTTemp);
    ph_chIso.push_back(chIsoTemp);
    ph_nuIso.push_back(nuIsoTemp);
    ph_phIso.push_back(phIsoTemp);
    ph_phIsoTight.push_back(isoPassTTemp);
    ph_phIsoMedium.push_back(isoPassMTemp);
    ph_phIsoLoose.push_back(isoPassLTemp);
    ph_isTight.push_back(isTightPhoton(photon));
  }  

  nPhotons = vPhotons.size();

 vector<TCElectron> vElectrons;

 for (Int_t i = 0; i < recoElectrons->GetSize(); i++) {
   TCElectron* electron = (TCElectron*) recoElectrons->At(i);
   vElectrons.push_back(*electron);
   el_pt.push_back(electron->Pt());
   el_eta.push_back(electron->Eta());
   el_phi.push_back(electron->Phi());
   el_energy.push_back(electron->Energy());
   el_charge.push_back(electron->Charge());
   el_iso.push_back(ElectronIso(electron));
   el_isTight.push_back(isTightElectron(electron));
   if(isMC){
     double closestDR = 0.3;
     int closestIndex=-1;
     for (int g = 0; g <  genParticles->GetSize(); g++) {
       TCGenParticle* genParticle = (TCGenParticle*) genParticles->At(g);
       if(abs(genParticle->GetPDGId())==11){
         double tmpDR = mdeltaR(electron->Eta(), electron->Phi(), genParticle->Eta(), genParticle->Phi());
         if(tmpDR<closestDR){
           closestDR = tmpDR;
           closestIndex = g;
          }
       }
    }//closing the gen particle loop
    if(closestIndex!=-1){ 
      TCGenParticle* genParticle = (TCGenParticle*) genParticles->At(closestIndex);
      if((genParticle->Mother() and abs(genParticle->Mother()->GetPDGId())==15)){ //Checking that the Tau (15) is the first mother.
         el_Matched.push_back(1);
         el_MatchedPt.push_back(electron->Pt());     
         el_MatchedEta.push_back(electron->Eta()); 
         el_MatchedPhi.push_back(electron->Phi()); 
         el_MatchedEnergy.push_back(electron->Energy()); 
        }//mother checked
      }//gen matching criteria checked
    else{
      el_Matched.push_back(0);
      el_MatchedPt.push_back(-99.0);
      el_MatchedEta.push_back(-99.0);
      el_MatchedPhi.push_back(-99.0);
      el_MatchedEnergy.push_back(-99.0);
        }//gen matching else statement closed
    }//isMC if statement ended 
 }//end reco electron loop
 
 nElectrons = vElectrons.size();

 if(isMC){
  for (int g = 0; g <  genParticles->GetSize(); g++) {
    TCGenParticle* genParticle = (TCGenParticle*) genParticles->At(g);
    if(abs(genParticle->GetPDGId())==15 and genParticle->GetStatus()==3){
      if(genParticle->Mother()->GetPDGId()==23){
       Ztt = true;
       }
     }
   }
 }


 vector<TCMuon> vMuons;

 for (Int_t i = 0; i < recoMuons->GetSize(); i++) {
   TCMuon* muon = (TCMuon*) recoMuons->At(i);
   vMuons.push_back(*muon);
   mu_pt.push_back(muon->Pt());
   mu_eta.push_back(muon->Eta());
   mu_phi.push_back(muon->Phi()); 
   mu_energy.push_back(muon->Energy());
   mu_charge.push_back(muon->Charge());
   mu_iso.push_back(MuonIso(muon));
   mu_isTight.push_back(isTightMuon(muon));
   if(isMC){
     double closestDR = 0.3;
     int closestIndex=-1;
     for (int g = 0; g <  genParticles->GetSize(); g++) {
       TCGenParticle* genParticle = (TCGenParticle*) genParticles->At(g);
       if(abs(genParticle->GetPDGId())==13 ){//and genParticle->GetStatus()==3){
        /* if((genParticle->Mother() and abs(genParticle->Mother()->GetPDGId())==15)){
            mu_Matched.push_back(1);
            mu_MatchedPt.push_back(genParticle->Pt());
            mu_MatchedEta.push_back(genParticle->Eta());
            mu_MatchedPhi.push_back(genParticle->Phi());
            mu_MatchedEnergy.push_back(genParticle->Energy());
           }
         else{
           mu_Matched.push_back(0);
           mu_MatchedPt.push_back(-99.0);
           mu_MatchedEta.push_back(-99.0);
           mu_MatchedPhi.push_back(-99.0);
           mu_MatchedEnergy.push_back(-99.0);
            }//else mother check
         }//gen muon check
      } //directly looking at the genParticles.
  */
        double tmpDR = mdeltaR(muon->Eta(), muon->Phi(), genParticle->Eta(), genParticle->Phi());
        if(tmpDR<closestDR){
           closestDR = tmpDR;
           closestIndex = g;
          }
       }
    }//closing the gen particle loop
    if(closestIndex!=-1){
      TCGenParticle* genParticle = (TCGenParticle*) genParticles->At(closestIndex);
      if((genParticle->Mother() and abs(genParticle->Mother()->GetPDGId())==15)){// and genParticle->Mother()->GetStatus()==3){ //Checking that the Tau (15) is the first mother.
         mu_Matched.push_back(1);
         mu_MatchedPt.push_back(muon->Pt());
         mu_MatchedEta.push_back(muon->Eta());
         mu_MatchedPhi.push_back(muon->Phi());
         mu_MatchedEnergy.push_back(muon->Energy());
        }//mother checked
      }//gen matching criteria checked
    else{
      mu_Matched.push_back(0);
      mu_MatchedPt.push_back(-99.0);
      mu_MatchedEta.push_back(-99.0);
      mu_MatchedPhi.push_back(-99.0);
      mu_MatchedEnergy.push_back(-99.0);
        }//gen matching else statement closed
    }//isMC if statement ended 
 }//end reco muon loop

 nMuons = vMuons.size();

 vector<TCJet> vJets;

 for (Int_t i = 0; i < patJets->GetSize(); i++) {
   TCJet* jet = (TCJet*) patJets->At(i);
   vJets.push_back(*jet);
   jet_pt.push_back(jet->Pt());
   jet_eta.push_back(jet->Eta());
   jet_phi.push_back(jet->Phi());
   jet_energy.push_back(jet->Energy());
   jet_mva_loose.push_back(jet->PuJetIdFlag_mva_loose());
   jet_mva_medium.push_back(jet->PuJetIdFlag_mva_medium());
   jet_mva_tight.push_back(jet->PuJetIdFlag_mva_tight());
   jet_cut_loose.push_back(jet->PuJetIdFlag_cut_loose());
   jet_cut_medium.push_back(jet->PuJetIdFlag_cut_medium());
   jet_cut_tight.push_back(jet->PuJetIdFlag_cut_tight());

  }

 nJets = vJets.size();
 
 MET = corrMET->Mod();
 MET_Phi = corrMET->Phi();
 MET_Px  = corrMET->Px();
 MET_Py  = corrMET->Py();

 if(isMC){ 
   PUWeightData = LumiWeightsD_.weight(nPUVerticesTrue);
   PUWeightDataSys = LumiWeightsD_sys_.weight(nPUVerticesTrue);
 }

 outtree->Fill();
 return kTRUE;
}

float FlatTreeCreator::ElectronIso(TCElectron *electron){
  
 float EA = 0;
 EAEle[0] = 0.208;
 EAEle[1] = 0.209;
 EAEle[2] = 0.115;
 EAEle[3] = 0.143;
 EAEle[4] = 0.183;
 EAEle[5] = 0.194;
 EAEle[6] = 0.261;
 if (fabs(electron->Eta())  <  1.0) EA = EAEle[0];
 else if(fabs(electron->Eta()) < 1.479) EA = EAEle[1]; 
 else if(fabs(electron->Eta()) < 2.0) EA = EAEle[2];   
 else if(fabs(electron->Eta()) < 2.2) EA = EAEle[3];
 else if(fabs(electron->Eta()) < 2.3) EA = EAEle[4]; 
 else if(fabs(electron->Eta()) < 2.4) EA = EAEle[5];
 else if(fabs(electron->Eta()) > 2.4) EA = EAEle[6];

 float combIso = (electron->IsoMap("pfChIso_R04")
    + max(0.,(double)electron->IsoMap("pfNeuIso_R04") + electron->IsoMap("pfPhoIso_R04") - rhoFactor*EA));
 return (combIso/(electron->Pt()));
}


bool FlatTreeCreator::isLooseElectron(TCElectron *electron){

  if(electron->SCEta() > 2.5) return false;
  if(fabs(electron->SCEta()) > 1.4442 and fabs(electron->SCEta()) < 1.566) return false;
  if(fabs(electron->Eta()) < 1.4442){
    if(fabs(electron->SCDeltaEta())    > 7.0e-03) return false;   
    if(fabs(electron->SCDeltaPhi())    > 1.5e-01) return false;
    if(electron->SigmaIEtaIEta()       > 1.0e-02) return false;
    if(electron->HadOverEm()           > 1.2e-01) return false;
    if(fabs(electron->Dxy(pvPosition)) > 0.02)    return false; 
    if(fabs(electron->Dz(pvPosition))  > 0.2)     return false; 
    if(electron->IdMap("fabsEPDiff")   > 0.05)    return false;
    if(electron->ConversionMissHits()  > 1)       return false;
    if(electron->PassConversionVeto() != 1)       return false;
  }

  if(fabs(electron->Eta()) > 1.566){
    if(fabs(electron->SCDeltaEta())    > 9.0e-03) return false;
    if(fabs(electron->SCDeltaPhi())    > 1.0e-01) return false;
    if(electron->SigmaIEtaIEta()       > 3.0e-02) return false;
    if(electron->HadOverEm()           > 1.0e-01) return false;
    if(fabs(electron->Dxy(pvPosition)) > 0.02)    return false;
    if(fabs(electron->Dz(pvPosition))  > 0.2)     return false;
    if(electron->IdMap("fabsEPDiff")   > 0.05)    return false;
    if(electron->ConversionMissHits()  > 1)       return false;
    if(electron->PassConversionVeto() != 1)       return false;
  }

  return true;

}

bool FlatTreeCreator::isMediumElectron(TCElectron *electron){

  if(electron->SCEta() > 2.5) return false;
  if(fabs(electron->SCEta()) > 1.4442 and fabs(electron->SCEta()) < 1.566) return false;
  if(fabs(electron->Eta()) < 1.4442){
    if(fabs(electron->SCDeltaEta())    > 4.0e-03) return false;
    if(fabs(electron->SCDeltaPhi())    > 0.6e-01) return false;
    if(electron->SigmaIEtaIEta()       > 1.0e-02) return false;
    if(electron->HadOverEm()           > 1.2e-01) return false;
    if(fabs(electron->Dxy(pvPosition)) > 0.02)    return false;
    if(fabs(electron->Dz(pvPosition))  > 0.1)     return false;
    if(electron->IdMap("fabsEPDiff")   > 0.05)    return false;
    if(electron->ConversionMissHits()  > 1)       return false;
    if(electron->PassConversionVeto() != 1)       return false;
  }

  if(fabs(electron->Eta()) > 1.566){
    if(fabs(electron->SCDeltaEta())    > 7.0e-03) return false;
    if(fabs(electron->SCDeltaPhi())    > 0.3e-01) return false;
    if(electron->SigmaIEtaIEta()       > 3.0e-02) return false;
    if(electron->HadOverEm()           > 1.0e-01) return false;
    if(fabs(electron->Dxy(pvPosition)) > 0.02)    return false;
    if(fabs(electron->Dz(pvPosition))  > 0.1)     return false;
    if(electron->IdMap("fabsEPDiff")   > 0.05)    return false;
    if(electron->ConversionMissHits()  > 1)       return false;
    if(electron->PassConversionVeto() != 1)       return false;
  }
  return true;

}

bool FlatTreeCreator::isTightElectron(TCElectron *electron){

  if(electron->SCEta() > 2.5) return false;
  if(fabs(electron->SCEta()) > 1.4442 and fabs(electron->SCEta()) < 1.566) return false;
  if(fabs(electron->Eta()) < 1.4442){
    if(fabs(electron->SCDeltaEta())    > 4.0e-03) return false;
    if(fabs(electron->SCDeltaPhi())    > 0.3e-01) return false;
    if(electron->SigmaIEtaIEta()       > 1.0e-02) return false;
    if(electron->HadOverEm()           > 1.2e-01) return false;
    if(fabs(electron->Dxy(pvPosition)) > 0.02)    return false;
    if(fabs(electron->Dz(pvPosition))  > 0.1)     return false;
    if(electron->IdMap("fabsEPDiff")   > 0.05)    return false;
    if(electron->ConversionMissHits()  > 0)       return false;
    if(electron->PassConversionVeto() != 1)       return false;
  }

  if(fabs(electron->Eta()) > 1.566){
    if(fabs(electron->SCDeltaEta())    > 5.0e-03) return false;
    if(fabs(electron->SCDeltaPhi())    > 0.2e-01) return false;
    if(electron->SigmaIEtaIEta()       > 3.0e-02) return false;
    if(electron->HadOverEm()           > 1.0e-01) return false;
    if(fabs(electron->Dxy(pvPosition)) > 0.02)    return false;
    if(fabs(electron->Dz(pvPosition))  > 0.1)     return false;
    if(electron->IdMap("fabsEPDiff")   > 0.05)    return false;
    if(electron->ConversionMissHits()  > 0)       return false;
    if(electron->PassConversionVeto() != 1)       return false;
  }
  return true;

}


float FlatTreeCreator::MuonIso(TCMuon *muon){
  
  float combIso = (muon->IsoMap("pfChargedHadronPt_R04")
    + max(0.,(double)muon->IsoMap("pfNeutralHadronEt_R04") + muon->IsoMap("pfPhotonEt_R04") - 0.5*muon->IsoMap("pfPUPt_R04")));

  return (combIso/(muon->Pt())); 

} 



bool FlatTreeCreator::isTightMuon(TCMuon *muon){

  if(fabs(muon->Eta()) > 2.4)                     return false;
  if(fabs(muon->Eta()) < 2.4){
    if(muon->IsPF() != 1)                         return false;
    if(muon->IsGLB()!= 1)                         return false;
    if(muon->NormalizedChi2() >= 10)              return false;
    if(muon->NumberOfValidMuonHits() <= 0)        return false;
    if(muon->NumberOfMatchedStations() <= 1)      return false;
    if(muon->NumberOfValidPixelHits() <= 0)       return false; 
    if(muon->TrackLayersWithMeasurement() <= 5)   return false;
    if(muon->Dxy(pvPosition) > 0.2)               return false;
    if(fabs(muon->Dz(pvPosition)) > 0.5)          return false;  
  }
  return true;

}

bool FlatTreeCreator::isLooseMuon(TCMuon *muon){

  if(fabs(muon->Eta()) > 2.4)                     return false;
  if(fabs(muon->Eta()) < 2.4){
    if(muon->IsPF() != 1)                         return false;
    if(muon->IsGLB()!= 1)                         return false;
    if(muon->NormalizedChi2() >= 50)              return false;
    if(muon->NumberOfValidMuonHits() <= 0)        return false;
    if(muon->NumberOfMatchedStations() <= 1)      return false;
    if(muon->NumberOfValidPixelHits() <= 0)       return false;
    if(muon->TrackLayersWithMeasurement() <= 5)   return false;
    if(muon->Dxy(pvPosition) > 2.0)               return false;
  }
  return true;

}

bool FlatTreeCreator::isTightPhoton(TCPhoton *photon){

  if(fabs(photon->SCEta()) > 2.5)                 return false;
  if(fabs(photon->SCEta()) > 1.4442 and fabs(photon->SCEta()) < 1.566) return false;    
  if(fabs(photon->SCEta()) < 1.4442){
    if(photon->HadOverEm() > 0.05)                return false;
    if(photon->SigmaIEtaIEta() >  0.011)          return false;
  }

  if(fabs(photon->SCEta()) > 1.566){
    if(photon->HadOverEm() > 0.05)                return false;
    if(photon->SigmaIEtaIEta() >  0.031)          return false;
  }

  return true;
}


int FlatTreeCreator::PhotonIso(TCPhoton *photon, double &chIso, double &nuIso, double &phIso, bool &isoPassL, bool &isoPassM, bool &isoPassT){
  float chEA,nhEA,phEA,chIsoCor,nhIsoCor,phIsoCor;
  float EAPho[7][3] = {
    {0.012,  0.030,   0.148}, //         eta < 1.0  
    {0.010,  0.057,   0.130}, // 1.0   < eta < 1.479   
    {0.014,  0.039,   0.112}, // 1.479 < eta < 2.0  
    {0.012,  0.015,   0.216}, // 2.0   < eta < 2.2 
    {0.016,  0.024,   0.262}, // 2.2   < eta < 2.3  
    {0.020,  0.039,   0.260}, // 2.3   < eta < 2.4 
    {0.012,  0.072,   0.266}  // 2.4   < eta       
  };

  
  if (fabs(photon->SCEta()) < 1.0){
    chEA = EAPho[0][0];
    nhEA = EAPho[0][1];
    phEA = EAPho[0][2];
  }else if (fabs(photon->SCEta()) < 1.479){
    chEA = EAPho[1][0];
    nhEA = EAPho[1][1];
    phEA = EAPho[1][2];
  }else if (fabs(photon->SCEta()) < 2.0){
    chEA = EAPho[2][0];
    nhEA = EAPho[2][1];
    phEA = EAPho[2][2];
  }else if (fabs(photon->SCEta()) < 2.2){
    chEA = EAPho[3][0];
    nhEA = EAPho[3][1];
    phEA = EAPho[3][2];
  }else if (fabs(photon->SCEta()) < 2.3){
    chEA = EAPho[4][0];
    nhEA = EAPho[4][1];
    phEA = EAPho[4][2];
  }else if (fabs(photon->SCEta()) < 2.4){
    chEA = EAPho[5][0];
    nhEA = EAPho[5][1];
    phEA = EAPho[5][2];
  }else{
    chEA = EAPho[6][0];
    nhEA = EAPho[6][1];
    phEA = EAPho[6][2];
  }

  chIsoCor = photon->IsoMap("chIso03")-rhoFactor*chEA;
  nhIsoCor = photon->IsoMap("nhIso03")-rhoFactor*nhEA;
  phIsoCor = photon->IsoMap("phIso03")-rhoFactor*phEA;

  chIso = chIsoCor;
  nuIso = nhIsoCor;
  phIso = phIsoCor;

  bool isoPassLoose = false;
  bool isoPassMedium = false;
  bool isoPassTight = false;

  if ((fabs(photon->SCEta()) < 1.4442
       and max((double)chIsoCor,0.)          < 2.6
       and max((double)nhIsoCor,0.)          < 3.5 + 0.04*photon->Pt()
       and max((double)phIsoCor,0.)          < 1.3 + 0.005*photon->Pt()
	 ) || (fabs(photon->SCEta()) > 1.566
       and max((double)chIsoCor,0.)          < 2.3
       and max((double)nhIsoCor,0.)          < 2.9 + 0.04*photon->Pt()
	 )) isoPassLoose = true;
  
  if ((fabs(photon->SCEta()) < 1.4442
       and max((double)chIsoCor,0.)          < 1.5
       and max((double)nhIsoCor,0.)          < 1.0 + 0.04*photon->Pt()
       and max((double)phIsoCor,0.)          < 0.7 + 0.005*photon->Pt()
         ) || (fabs(photon->SCEta()) > 1.566
       and max((double)chIsoCor,0.)          < 1.2
       and max((double)nhIsoCor,0.)          < 1.5 + 0.04*photon->Pt()
       and max((double)nhIsoCor,0.)          < 1.0 + 0.005*photon->Pt()
         )) isoPassMedium = true;

  if ((fabs(photon->SCEta()) < 1.4442
       and max((double)chIsoCor,0.)          < 0.7
       and max((double)nhIsoCor,0.)          < 0.4 + 0.04*photon->Pt()
       and max((double)phIsoCor,0.)          < 0.5 + 0.005*photon->Pt()
         ) || (fabs(photon->SCEta()) > 1.566
       and max((double)chIsoCor,0.)          < 0.5
       and max((double)nhIsoCor,0.)          < 1.5 + 0.04*photon->Pt()
       and max((double)nhIsoCor,0.)          < 1.0 + 0.005*photon->Pt()
         )) isoPassTight = true;

  isoPassL = isoPassLoose;
  isoPassM = isoPassMedium;
  isoPassT = isoPassTight; 
  return 0;
}


bool FlatTreeCreator::isMediumPhoton(TCPhoton *photon){

  if(fabs(photon->SCEta()) > 2.5)                 return false;
  if(fabs(photon->SCEta()) > 1.4442 and fabs(photon->SCEta()) < 1.566) return false;
  if(fabs(photon->SCEta()) < 1.4442){
    if(photon->HadOverEm() > 0.05)                return false;
    if(photon->SigmaIEtaIEta() >  0.011)          return false;
  }

  if(fabs(photon->SCEta()) > 1.566){
    if(photon->HadOverEm() > 0.05)                return false;
    if(photon->SigmaIEtaIEta() >  0.033)          return false;
  }

  return true;
}

bool FlatTreeCreator::isLoosePhoton(TCPhoton *photon){

  if(fabs(photon->SCEta()) > 2.5)                 return false;
  if(fabs(photon->SCEta()) > 1.4442 and fabs(photon->SCEta()) < 1.566) return false;
  if(fabs(photon->SCEta()) < 1.4442){
    if(photon->HadOverEm() > 0.05)                return false;
    if(photon->SigmaIEtaIEta() >  0.012)          return false;
  }

  if(fabs(photon->SCEta()) > 1.566){
    if(photon->HadOverEm() > 0.05)                return false;
    if(photon->SigmaIEtaIEta() >  0.034)          return false;
  }

  return true;
}

double FlatTreeCreator::mdeltaR(double eta1, double phi1, double eta2, double phi2) {
  double delta_eta = fabs(eta1-eta2);
  double delta_phi = fabs(phi1-phi2);
  if(delta_phi > 3.14159265) delta_phi = delta_phi - 2*3.14159265;
  double pt = 5.0; //some non zero value
  TVector3 vector1;
  TVector3 vector2;
  vector1.SetPtEtaPhi(pt, eta1, phi1);
  vector2.SetPtEtaPhi(pt, eta2, phi2);
  double deltaR1 = vector1.DeltaR(vector2);
  double deltaR2 = std::sqrt((delta_eta)*(delta_eta) + (delta_phi)*(delta_phi));
  if (fabs(deltaR1-deltaR2) > 0.0001) std::cout<<"Difference spotted "<<deltaR1<<", "<<deltaR2<<std::endl;
  return std::sqrt((delta_eta)*(delta_eta) + (delta_phi)*(delta_phi));

}

void FlatTreeCreator::SlaveTerminate()
{
}

void FlatTreeCreator::Terminate()
{
 outtree->Write();
}



