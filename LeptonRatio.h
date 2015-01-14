// -*- C++ -*-
//
// Package:    LeptonRatio
// Class:      LeptonRatio
// 
/**\class LeptonRatio LeptonRatio.cc LeptonRatio/LeptonRatio/src/LeptonRatio.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Saptaparna Bhattacharya
//         Created:  Tue Jan 13 15:18:16 CST 2015
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
//
// class declaration
//
class LeptonRatio : public edm::EDAnalyzer {
   public:
      explicit LeptonRatio(const edm::ParameterSet&);
      ~LeptonRatio();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      int electrons;
      int muons;
      int taus;
      int n_dielectrons;
      int n_dimuons;
      int n_dileptons;
      int n_ditaus;
      int n_muonTau;
      int n_electronTau;
      double pt;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//

