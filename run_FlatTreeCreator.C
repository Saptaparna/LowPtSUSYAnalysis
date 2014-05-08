#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <stdio.h>
#include "TChain.h"
#include "TROOT.h"
#include "TStopwatch.h" 


void run_FlatTreeCreator(){

  gROOT->LoadMacro("Classes/src/TCPhysObject_cc.so");
  gROOT->LoadMacro("Classes/src/TCJet_cc.so");
  gROOT->LoadMacro("Classes/src/TCMET_cc.so");
  gROOT->LoadMacro("Classes/src/TCEGamma_cc.so");
  gROOT->LoadMacro("Classes/src/TCElectron_cc.so");
  gROOT->LoadMacro("Classes/src/TCMuon_cc.so");
  gROOT->LoadMacro("Classes/src/TCTau_cc.so");
  gROOT->LoadMacro("Classes/src/TCPhoton_cc.so");
  gROOT->LoadMacro("Classes/src/TCGenJet_cc.so");
  gROOT->LoadMacro("Classes/src/TCGenParticle_cc.so");
  gROOT->LoadMacro("Classes/src/TCPrimaryVtx_cc.so");
  gROOT->LoadMacro("Classes/src/TCTriggerObject_cc.so");
  gROOT->LoadMacro("plugins/HistManager_cc.so");
  gROOT->LoadMacro("plugins/TreeManager_C.so");
  gROOT->LoadMacro("plugins/TriggerSelector_cc.so");
 
  TChain* fChain = new TChain("ntupleProducer/eventTree");
  ifstream sourceFiles("TEST_DATAFILE.txt");
  char line[128];
  int  count = 0;
  cout<< "Adding files from TEST_DATAFILE.txt to chain..."<< endl;
  while (sourceFiles >> line) {
    fChain->Add(line);
    count++;
    }
  cout << count << " file/files have been added" <<endl;
  sourceFiles.close();

  fChain->Process("FlatTreeCreator.C+");

}
