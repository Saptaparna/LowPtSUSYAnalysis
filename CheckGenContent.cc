// -*- C++ -*-
//
// Package:    CheckGenContent
// Class:      CheckGenContent
// 
/**\class CheckGenContent CheckGenContent.cc CheckGenContent/CheckGenContent/src/CheckGenContent.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Saptaparna Bhattacharya
//         Created:  Thu May 22 12:21:29 CDT 2014
// $Id$
//
//

#include"CheckGenContent/CheckGenContent/interface/CheckGenContent.h"
#include"DataFormats/HepMCCandidate/interface/GenParticle.h"
#include"RecoEgamma/Examples/plugins/GsfElectronMCAnalyzer.h"
#include"FWCore/ParameterSet/interface/ParameterSet.h"
#include"FWCore/Framework/interface/EDAnalyzer.h"
#include"FWCore/Framework/interface/Event.h"
#include"DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

CheckGenContent::CheckGenContent(const edm::ParameterSet& iConfig)

{
    electrons=0;
    muons=0;
    taus = 0;
    n_ditaus = 0;
    n_dielectrons=0;
    n_dimuons=0;

}


CheckGenContent::~CheckGenContent()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
CheckGenContent::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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
    double pt = genParticle->pt();
    if (abs(id)==11 && st==3){
       electrons++;
       n_electrons++;
                }

    if (abs(id)==13 && st==3){
        muons++;
        n_muons++;
                }
 
    if (abs(id)==15 && st==3){
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

    if(n_taus==2 && n_electrons==0 && n_muons==0){
       n_ditaus++;
        }

}


// ------------ method called once each job just before starting event loop  ------------
void 
CheckGenContent::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
CheckGenContent::endJob() 
{
  std::cout << "electrons = " << electrons << std::endl;
  std::cout << "muons = " << muons << std::endl;
  std::cout << "taus = " << taus << std::endl;  
  std::cout << "n_ditaus = " << n_ditaus << std::endl;   
  std::cout << "n_dielectrons = " << n_dielectrons << std::endl;   
  std::cout << "n_dimuons = " << n_dimuons << std::endl;    

}

// ------------ method called when starting to processes a run  ------------
void 
CheckGenContent::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
CheckGenContent::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
CheckGenContent::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
CheckGenContent::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
CheckGenContent::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(CheckGenContent);
