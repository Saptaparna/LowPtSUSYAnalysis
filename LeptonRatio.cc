#include"LeptonRatio/LeptonRatio/interface/LeptonRatio.h"
#include"DataFormats/HepMCCandidate/interface/GenParticle.h"
#include"RecoEgamma/Examples/plugins/GsfElectronMCAnalyzer.h"
#include"FWCore/ParameterSet/interface/ParameterSet.h"
#include"FWCore/Framework/interface/EDAnalyzer.h"
#include"FWCore/Framework/interface/Event.h"
#include"DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include"TLorentzVector.h"

LeptonRatio::LeptonRatio(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
    electrons=0;
    muons=0;
    taus=0;
    nues=0;
    numus=0;
    nutaus=0;
    n_dielectrons=0;
    n_dimuons=0;
    n_ditaus=0;
    n_dileptons=0;
    n_electronTau=0;
    n_muonTau=0;
    Wmass=0;
    NeutralinoMass=0;
    CharginoMass=0;
    Ws=0;
}


LeptonRatio::~LeptonRatio()
{
 
}

void
LeptonRatio::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   Handle<std::vector<reco::GenParticle> > genParticles;
   iEvent.getByLabel("genParticles", genParticles) ;
   int n_electrons=0;
   int n_muons=0;
   int n_taus=0;
   int n_nues=0;
   int n_numus=0;
   int n_nutaus=0;
   int n_Ws=0;
   genElectronPt.clear();
   genElectronEta.clear();
   genElectronPhi.clear();
   genElectronEnergy.clear();
   genMuonPt.clear();
   genMuonEta.clear();
   genMuonPhi.clear();
   genMuonEnergy.clear();
   genTauPt.clear();
   genTauEta.clear();
   genTauPhi.clear();
   genTauEnergy.clear();
   genElectronNeutrinoPt.clear();
   genElectronNeutrinoEta.clear();
   genElectronNeutrinoPhi.clear();
   genElectronNeutrinoEnergy.clear(); 
   genMuonNeutrinoPt.clear();
   genMuonNeutrinoEta.clear();
   genMuonNeutrinoPhi.clear();
   genMuonNeutrinoEnergy.clear();
   genTauNeutrinoPt.clear();
   genTauNeutrinoEta.clear();
   genTauNeutrinoPhi.clear();
   genTauNeutrinoEnergy.clear();
   genNeutralinoPt.clear();
   genNeutralinoEta.clear();
   genNeutralinoPhi.clear();
   genNeutralinoEnergy.clear();
   genCharginoPt.clear();
   genCharginoEta.clear();
   genCharginoPhi.clear();
   genCharginoEnergy.clear();
   genWId_El.clear();
   genWPt_El.clear();
   genWEta_El.clear();
   genWPhi_El.clear();   
   genWEnergy_El.clear();
   genWId_Mu.clear();
   genWPt_Mu.clear();
   genWEta_Mu.clear();
   genWPhi_Mu.clear();
   genWEnergy_Mu.clear();
   genWId_Tau.clear();
   genWPt_Tau.clear();
   genWEta_Tau.clear();
   genWPhi_Tau.clear();
   genWEnergy_Tau.clear();

   for(std::vector<reco::GenParticle>::const_iterator genParticle = genParticles->begin(); genParticle != genParticles->end(); genParticle++){
     int id = genParticle->pdgId();
     int st = genParticle->status(); 
     int mother_id = 0;
     if(genParticle->mother()) mother_id = genParticle->mother()->pdgId();
     
     if(abs(id)==11 and st==3){// and abs(mother_id)==24){ 
       electrons++;
       n_electrons++;
       genElectronPt.push_back(genParticle->pt()); 
       genElectronEta.push_back(genParticle->eta());
       genElectronPhi.push_back(genParticle->phi());
       genElectronEnergy.push_back(genParticle->energy());
       int nudaughter_id = -1;
       /*if(abs(genParticle->mother()->daughter(0)->pdgId()) == 12) nudaughter_id=0;
       else if(abs(genParticle->mother()->daughter(1)->pdgId()) == 12) nudaughter_id=1;
       if (nudaughter_id==-1) std::cout<<"No neutrino daughter found!"<<std::endl;
       else{
         genElectronNeutrinoPt.push_back(genParticle->mother()->daughter(nudaughter_id)->pt());
         genElectronNeutrinoEta.push_back(genParticle->mother()->daughter(nudaughter_id)->eta());
         genElectronNeutrinoPhi.push_back(genParticle->mother()->daughter(nudaughter_id)->phi());
         genElectronNeutrinoEnergy.push_back(genParticle->mother()->daughter(nudaughter_id)->energy());
       }
       if(genParticle->mother()->mother()) genWId_El.push_back(genParticle->mother()->mother()->pdgId());
       genWPt_El.push_back(genParticle->mother()->pt());
       genWEta_El.push_back(genParticle->mother()->eta());
       genWPhi_El.push_back(genParticle->mother()->phi());
       genWEnergy_El.push_back(genParticle->mother()->energy());
       std::cout<<"electron + neutrino energy = "<<genParticle->energy() + genParticle->mother()->daughter(nudaughter_id)->energy() <<"; W energy = "<<genParticle->mother()->energy()<<std::endl;*/
     }

    if(abs(id)==13 and st==3){// and abs(mother_id)==24){
       muons++;
       n_muons++;
       genMuonPt.push_back(genParticle->pt());
       genMuonEta.push_back(genParticle->eta());
       genMuonPhi.push_back(genParticle->phi());
       genMuonEnergy.push_back(genParticle->energy());
       int nudaughter_id = -1;
  /*     if(abs(genParticle->mother()->daughter(0)->pdgId()) == 14) nudaughter_id=0;
       else if(abs(genParticle->mother()->daughter(1)->pdgId()) == 14) nudaughter_id=1;
       if (nudaughter_id==-1) std::cout<<"No neutrino daughter found!"<<std::endl;
       else{
         genMuonNeutrinoPt.push_back(genParticle->mother()->daughter(nudaughter_id)->pt());
         genMuonNeutrinoEta.push_back(genParticle->mother()->daughter(nudaughter_id)->eta());
         genMuonNeutrinoPhi.push_back(genParticle->mother()->daughter(nudaughter_id)->phi());
         genMuonNeutrinoEnergy.push_back(genParticle->mother()->daughter(nudaughter_id)->energy());
       }
       if(genParticle->mother()->mother()) genWId_Mu.push_back(genParticle->mother()->mother()->pdgId());
       genWPt_Mu.push_back(genParticle->mother()->pt());
       genWEta_Mu.push_back(genParticle->mother()->eta());
       genWPhi_Mu.push_back(genParticle->mother()->phi());
       genWEnergy_Mu.push_back(genParticle->mother()->energy());
       std::cout<<"muon + neutrino energy = "<<genParticle->energy() + genParticle->mother()->daughter(nudaughter_id)->energy() <<"; W energy = "<<genParticle->mother()->energy()<<std::endl;*/
     }    

    if(abs(id)==15 and st==3){// and abs(mother_id)==24){
       taus++;
       n_taus++;
       genTauPt.push_back(genParticle->pt());
       genTauEta.push_back(genParticle->eta());
       genTauPhi.push_back(genParticle->phi());
       genTauEnergy.push_back(genParticle->energy());
       int nudaughter_id = -1;
       /*if(abs(genParticle->mother()->daughter(0)->pdgId()) == 16) nudaughter_id=0;
       else if(abs(genParticle->mother()->daughter(1)->pdgId()) == 16) nudaughter_id=1;
       if (nudaughter_id==-1) std::cout<<"No neutrino daughter found!"<<std::endl;
       else{
         genTauNeutrinoPt.push_back(genParticle->mother()->daughter(nudaughter_id)->pt());
         genTauNeutrinoEta.push_back(genParticle->mother()->daughter(nudaughter_id)->eta());
         genTauNeutrinoPhi.push_back(genParticle->mother()->daughter(nudaughter_id)->phi());
         genTauNeutrinoEnergy.push_back(genParticle->mother()->daughter(nudaughter_id)->energy());
       }
       if(genParticle->mother()->mother()) genWId_Tau.push_back(genParticle->mother()->mother()->pdgId());
       genWPt_Tau.push_back(genParticle->mother()->pt());
       genWEta_Tau.push_back(genParticle->mother()->eta());
       genWPhi_Tau.push_back(genParticle->mother()->phi());
       genWEnergy_Tau.push_back(genParticle->mother()->energy());
       std::cout<<"tau + neutrino energy = "<<genParticle->energy() + genParticle->mother()->daughter(nudaughter_id)->energy() <<"; W energy = "<<genParticle->mother()->energy()<<std::endl;
*/
    }
/*
    if(abs(id)==1000022 and st==3 and abs(mother_id)==1000024){
       genNeutralinoPt.push_back(genParticle->pt());
       genNeutralinoEta.push_back(genParticle->eta());
       genNeutralinoPhi.push_back(genParticle->phi());
       genNeutralinoEnergy.push_back(genParticle->energy());
      }

    if(abs(id)==1000024 and st==3 and abs(mother_id)==1000006){
       genCharginoPt.push_back(genParticle->pt());
       genCharginoEta.push_back(genParticle->eta());
       genCharginoPhi.push_back(genParticle->phi());
       genCharginoEnergy.push_back(genParticle->energy());
       } */
     }
   if(n_electrons==2 && n_muons==0 && n_taus==0){
       n_dielectrons++;
       }

     if(n_muons==2 && n_electrons==0 && n_taus==0){
       n_dimuons++;
       }

     if(n_electrons==1 && n_muons==1 && n_taus==0){
       n_dileptons++;
       }
    
     if(n_taus==2 && n_muons==0 && n_electrons==0){
       n_ditaus++;
       }
    
     if(n_electrons==1 && n_taus==1 && n_muons==0){
       n_electronTau++;
       } 
     
     if(n_electrons==0 && n_taus==1 && n_muons==1){
       n_muonTau++;
       }    
/*
   TLorentzVector leadingEl;
   TLorentzVector leadingElNu;
   if(genElectronPt.size() > 0) leadingEl.SetPtEtaPhiE(genElectronPt[0], genElectronEta[0], genElectronPhi[0], genElectronEnergy[0]);
   if(genElectronNeutrinoPt.size() > 0) leadingElNu.SetPtEtaPhiE(genElectronNeutrinoPt[0], genElectronNeutrinoEta[0], genElectronNeutrinoPhi[0], genElectronNeutrinoEnergy[0]);
   Wmass = (leadingEl+leadingElNu).M();
   TLorentzVector neutralino;
   if(genNeutralinoPt.size()>0.0) neutralino.SetPtEtaPhiE(genNeutralinoPt[0], genNeutralinoEta[0], genNeutralinoPhi[0], genNeutralinoEnergy[0]);
   NeutralinoMass = neutralino.M();
   std::cout << "NeutralinoMass = " << NeutralinoMass << std::endl;
   CharginoMass = (leadingEl+leadingElNu+neutralino).M();*/
    _tree->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void 
LeptonRatio::beginJob()
{
 _file = new TFile("LeptonProperties_GenLevel.root", "RECREATE");
 _tree = new TTree("LeptonProperties_GenLevel", "LeptonProperties_GenLevel", 64000000);
/* _branch = _tree->Branch("genElectronPt", &genElectronPt);
 _branch = _tree->Branch("genElectronEta", &genElectronEta);
 _branch = _tree->Branch("genElectronPhi", &genElectronPhi);
 _branch = _tree->Branch("genElectronEnergy", &genElectronEnergy);
 _branch = _tree->Branch("genElectronNeutrinoPt", &genElectronNeutrinoPt);
 _branch = _tree->Branch("genElectronNeutrinoEta", &genElectronNeutrinoEta);
 _branch = _tree->Branch("genElectronNeutrinoPhi", &genElectronNeutrinoPhi);
 _branch = _tree->Branch("genElectronNeutrinoEnergy", &genElectronNeutrinoEnergy);
 _branch = _tree->Branch("genMuonPt", &genMuonPt);
 _branch = _tree->Branch("genMuonEta", &genMuonEta);
 _branch = _tree->Branch("genMuonPhi", &genMuonPhi);
 _branch = _tree->Branch("genMuonEnergy", &genMuonEnergy);
 _branch = _tree->Branch("genMuonNeutrinoPt", &genMuonNeutrinoPt);
 _branch = _tree->Branch("genMuonNeutrinoEta", &genMuonNeutrinoEta);
 _branch = _tree->Branch("genMuonNeutrinoPhi", &genMuonNeutrinoPhi);
 _branch = _tree->Branch("genMuonNeutrinoEnergy", &genMuonNeutrinoEnergy);
 _branch = _tree->Branch("genTauPt", &genTauPt);
 _branch = _tree->Branch("genTauEta", &genTauEta);
 _branch = _tree->Branch("genTauPhi", &genTauPhi);
 _branch = _tree->Branch("genTauEnergy", &genTauEnergy);
 _branch = _tree->Branch("genTauNeutrinoPt", &genTauNeutrinoPt);
 _branch = _tree->Branch("genTauNeutrinoEta", &genTauNeutrinoEta);
 _branch = _tree->Branch("genTauNeutrinoPhi", &genTauNeutrinoPhi);
 _branch = _tree->Branch("genTauNeutrinoEnergy", &genTauNeutrinoEnergy);
 _branch = _tree->Branch("genWId_El", &genWId_El);
 _branch = _tree->Branch("genWPt_El", &genWPt_El);
 _branch = _tree->Branch("genWEta_El", &genWEta_El);
 _branch = _tree->Branch("genWPhi_El", &genWPhi_El);
 _branch = _tree->Branch("genWEnergy_El", &genWEnergy_El);
 _branch = _tree->Branch("genWId_Mu", &genWId_Mu);
 _branch = _tree->Branch("genWPt_Mu", &genWPt_Mu);
 _branch = _tree->Branch("genWEta_Mu", &genWEta_Mu);
 _branch = _tree->Branch("genWPhi_Mu", &genWPhi_Mu);
 _branch = _tree->Branch("genWEnergy_Mu", &genWEnergy_Mu);
 _branch = _tree->Branch("genWId_Tau", &genWId_Tau);
 _branch = _tree->Branch("genWPt_Tau", &genWPt_Tau);
 _branch = _tree->Branch("genWEta_Tau", &genWEta_Tau);
 _branch = _tree->Branch("genWPhi_Tau", &genWPhi_Tau);
 _branch = _tree->Branch("genWEnergy_Tau", &genWEnergy_Tau);
 _branch = _tree->Branch("Wmass", &Wmass, "Wmass/D");
 _branch = _tree->Branch("NeutralinoMass", &NeutralinoMass, "NeutralinoMass/D");
 _branch = _tree->Branch("CharginoMass", &CharginoMass, "CharginoMass/D");*/
}

// ------------ method called once each job just after ending the event loop  ------------
void 
LeptonRatio::endJob() 
{
std::cout << "electrons = " << electrons << std::endl;
std::cout << "muons = " << muons << std::endl;
std::cout << "taus = " << taus << std::endl;
std::cout << "n_dielectrons = " << n_dielectrons << std::endl;
std::cout << "n_dimuons = " << n_dimuons << std::endl;
std::cout << "n_dileptons = " << n_dileptons << std::endl;
std::cout << "n_ditaus = " << n_ditaus  << std::endl;
std::cout << "n_electronTau = " << n_electronTau << std::endl;
std::cout << "n_muonTau = " << n_muonTau << std::endl;
_file->Write();
_file->Close();

}

// ------------ method called when starting to processes a run  ------------
void 
LeptonRatio::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
LeptonRatio::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
LeptonRatio::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
LeptonRatio::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
LeptonRatio::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(LeptonRatio);
