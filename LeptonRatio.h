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
#include"TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include <TChain.h>
#include <TLorentzVector.h>
#include <TLegend.h>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <algorithm>
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
      TFile * _file;
      TTree * _tree;
      TBranch * _branch;
      int electrons;
      int muons;
      int taus;
      int nues;
      int numus;
      int nutaus;
      int n_dielectrons;
      int n_dimuons;
      int n_dileptons;
      int n_ditaus;
      int n_muonTau;
      int n_electronTau;
      int Ws;
      int n_Ws;
      double pt;
      double Wmass;
      double NeutralinoMass;
      double CharginoMass;
      std::vector<double> genWId_El; 
      std::vector<double> genWPt_El;
      std::vector<double> genWEta_El;
      std::vector<double> genWPhi_El;
      std::vector<double> genWEnergy_El;
    
      std::vector<double> genWId_Mu;
      std::vector<double> genWPt_Mu;
      std::vector<double> genWEta_Mu;
      std::vector<double> genWPhi_Mu;
      std::vector<double> genWEnergy_Mu;

      std::vector<double> genWId_Tau;
      std::vector<double> genWPt_Tau;
      std::vector<double> genWEta_Tau;
      std::vector<double> genWPhi_Tau;
      std::vector<double> genWEnergy_Tau;

      std::vector<double> genElectronPt;
      std::vector<double> genElectronEta;
      std::vector<double> genElectronPhi; 
      std::vector<double> genElectronEnergy;

      std::vector<double> genMuonPt;
      std::vector<double> genMuonEta;
      std::vector<double> genMuonPhi;
      std::vector<double> genMuonEnergy;

      std::vector<double> genTauPt;
      std::vector<double> genTauEta;
      std::vector<double> genTauPhi;
      std::vector<double> genTauEnergy;

      std::vector<double> genElectronNeutrinoPt;
      std::vector<double> genElectronNeutrinoEta;
      std::vector<double> genElectronNeutrinoPhi;
      std::vector<double> genElectronNeutrinoEnergy;

      std::vector<double> genMuonNeutrinoPt;
      std::vector<double> genMuonNeutrinoEta;
      std::vector<double> genMuonNeutrinoPhi;
      std::vector<double> genMuonNeutrinoEnergy;
      
      std::vector<double> genTauNeutrinoPt;
      std::vector<double> genTauNeutrinoEta;
      std::vector<double> genTauNeutrinoPhi;
      std::vector<double> genTauNeutrinoEnergy;

      std::vector<double> genNeutralinoPt;
      std::vector<double> genNeutralinoEta;
      std::vector<double> genNeutralinoPhi;
      std::vector<double> genNeutralinoEnergy;

      std::vector<double> genCharginoPt;
      std::vector<double> genCharginoEta;
      std::vector<double> genCharginoPhi;
      std::vector<double> genCharginoEnergy;

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

