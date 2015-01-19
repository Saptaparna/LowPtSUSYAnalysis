#include"LeptonRatio/LeptonRatio/interface/LeptonRatio.h"
#include"DataFormats/HepMCCandidate/interface/GenParticle.h"
#include"RecoEgamma/Examples/plugins/GsfElectronMCAnalyzer.h"
#include"FWCore/ParameterSet/interface/ParameterSet.h"
#include"FWCore/Framework/interface/EDAnalyzer.h"
#include"FWCore/Framework/interface/Event.h"
#include"DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

LeptonRatio::LeptonRatio(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
    electrons=0;
    muons=0;
    taus=0;
    n_dielectrons=0;
    n_dimuons=0;
    n_ditaus=0;
    n_dileptons=0;
    n_electronTau=0;
    n_muonTau=0;
}


LeptonRatio::~LeptonRatio()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
LeptonRatio::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   Handle<std::vector<reco::GenParticle> > genParticles;
   iEvent.getByLabel("genParticles", genParticles) ;

   int n_electrons=0;
   int n_muons=0;
   int n_taus=0;
   for(std::vector<reco::GenParticle>::const_iterator genParticle = genParticles->begin(); genParticle != genParticles->end(); genParticle++){
    int id = genParticle->pdgId();
    int st = genParticle->status();
    int mother_id = 0;
    if( genParticle->mother()) mother_id = genParticle->mother()->pdgId();
    double pt = genParticle->pt();
    if (abs(id)==11 && st==3 and pt>0.5 and abs(mother_id)==24){
       electrons++;
       n_electrons++;
                }

    if (abs(id)==13 && st==3 and pt>0.5 and abs(mother_id)==24){
        muons++;
        n_muons++;
                }

    if (abs(id)==15 && st==3 and pt>0.5 and abs(mother_id)==24){
       taus++;
       n_taus++;
       }
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

}


// ------------ method called once each job just before starting event loop  ------------
void 
LeptonRatio::beginJob()
{
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
