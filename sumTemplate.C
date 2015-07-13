#include <vector>
#include <iostream>
#include <string.h>
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TKey.h"
#include "Riostream.h"
#include "TGraphErrors.h"
#include "TCanvas.h"

using std::cout;
using std::endl;

void sumTemplate(){

TFile* file1 = TFile::Open("Output_LowPtSUSY_Tree_QCD_Pt_15to3000_TuneZ2_Flat_8TeV_pythia6_May26_All_FakePhotonTemplate.root");
TFile* file2 = TFile::Open("Output_LowPtSUSY_Tree_G_Pt_15to3000_All_FakePhotonTemplate.root");

TH1F* h1 = (TH1F*)file1->Get("h_SigmaIetaIeta_Num_Ph30To40");
h1->Scale((2.99815997E10*7.4*1000)/9904852);

TH1F* h2 = (TH1F*)file2->Get("h_SigmaIetaIeta_Num_Ph30To40");
h2->Scale((2.06583E7*7.4*1000)/9755606);

cout << "h1 + h2 entries = " << long(h1->Integral() + h2->Integral())  << endl;

h1->Add(h2);

cout << "h1 entries = " << long(h1->Integral()) << endl;
 
h1->SetName("h_SigmaIetaIeta_Num_Ph30To40");

TH1F* h1a = (TH1F*)file1->Get("h_SigmaIetaIeta_SB_Ph30To40");
h1a->Scale((2.99815997E10*7.4*1000)/9904852);

TH1F* h2a = (TH1F*)file2->Get("h_SigmaIetaIeta_SB_Ph30To40");
h2a->Scale((2.06583E7*7.4*1000)/9755606);

cout << "h1a + h2a entries = " << long(h1a->Integral() + h2a->Integral())  << endl;

h1a->Add(h2a);

cout << "h1a entries = " << long(h1a->Integral()) << endl;

h1a->SetName("h_SigmaIetaIeta_SB_Ph30To40");

TH1F* h3 = (TH1F*)file1->Get("h_SigmaIetaIeta_Num_Ph40To50");
h3->Scale((2.99815997E10*7.4*1000)/9904852);

TH1F* h4 = (TH1F*)file2->Get("h_SigmaIetaIeta_Num_Ph40To50");
h4->Scale((2.06583E7*7.4*1000)/9755606);

cout << "h3 + h4 entries = " << long(h3->Integral() + h4->Integral())  << endl;

h3->Add(h4);

cout << "h3 entries = " << long(h3->Integral()) << endl;

h3->SetName("h_SigmaIetaIeta_Num_Ph40To50");

TH1F* h3a = (TH1F*)file1->Get("h_SigmaIetaIeta_SB_Ph40To50");
h3a->Scale((2.99815997E10*7.4*1000)/9904852);

TH1F* h4a = (TH1F*)file2->Get("h_SigmaIetaIeta_SB_Ph40To50");
h4a->Scale((2.06583E7*7.4*1000)/9755606);

cout << "h3a + h4a entries = " << long(h3a->Integral() + h4a->Integral())  << endl;

h3a->Add(h4a);

cout << "h3a entries = " << long(h3a->Integral()) << endl;

h3a->SetName("h_SigmaIetaIeta_SB_Ph40To50");

TH1F* h5 = (TH1F*)file1->Get("h_SigmaIetaIeta_Num_Ph50To60");
h5->Scale((2.99815997E10*7.4*1000)/9904852);

TH1F* h6 = (TH1F*)file2->Get("h_SigmaIetaIeta_Num_Ph50To60");
h6->Scale((2.06583E7*7.4*1000)/9755606);

cout << "h5 + h6 entries = " << long(h5->Integral() + h6->Integral())  << endl;

h5->Add(h6);

cout << "h5 entries = " << long(h5->Integral()) << endl;

h5->SetName("h_SigmaIetaIeta_Num_Ph50To60");

TH1F* h5a = (TH1F*)file1->Get("h_SigmaIetaIeta_SB_Ph50To60");
h5a->Scale((2.99815997E10*7.4*1000)/9904852);

TH1F* h6a = (TH1F*)file2->Get("h_SigmaIetaIeta_SB_Ph50To60");
h6a->Scale((2.06583E7*7.4*1000)/9755606);

cout << "h5a + h6a entries = " << long(h5a->Integral() + h6a->Integral())  << endl;

h5a->Add(h6a);

cout << "h5a entries = " << long(h5a->Integral()) << endl;

h5a->SetName("h_SigmaIetaIeta_SB_Ph50To60");

TH1F* h7 = (TH1F*)file1->Get("h_SigmaIetaIeta_Num_Ph60To70");
h7->Scale((2.99815997E10*7.4*1000)/9904852);

TH1F* h8 = (TH1F*)file2->Get("h_SigmaIetaIeta_Num_Ph60To70");
h8->Scale((2.06583E7*7.4*1000)/9755606);

cout << "h7 + h8 entries = " << long(h7->Integral() + h8->Integral())  << endl;

h7->Add(h8);

cout << "h7 entries = " << long(h7->Integral()) << endl;

h7->SetName("h_SigmaIetaIeta_Num_Ph70To60");

TH1F* h7a = (TH1F*)file1->Get("h_SigmaIetaIeta_SB_Ph60To70");
h7a->Scale((2.99815997E10*7.4*1000)/9904852);

TH1F* h8a = (TH1F*)file2->Get("h_SigmaIetaIeta_SB_Ph60To70");
h8a->Scale((2.06583E7*7.4*1000)/9755606);

cout << "h7a + h8a entries = " << long(h7a->Integral() + h8a->Integral())  << endl;

h7a->Add(h8a);

cout << "h7a entries = " << long(h7a->Integral()) << endl;

h7a->SetName("h_SigmaIetaIeta_SB_Ph60To70");

TH1F* h9 = (TH1F*)file1->Get("h_SigmaIetaIeta_Num_Ph70To80");
h9->Scale((2.99815997E10*7.4*1000)/9904852);

TH1F* h10 = (TH1F*)file2->Get("h_SigmaIetaIeta_Num_Ph70To80");
h10->Scale((2.06583E7*7.4*1000)/9755606);

cout << "h9 + h10 entries = " << long(h9->Integral() + h10->Integral())  << endl;

h9->Add(h10);

cout << "h9 entries = " << long(h9->Integral()) << endl;

h9->SetName("h_SigmaIetaIeta_Num_Ph70To80");

TH1F* h9a = (TH1F*)file1->Get("h_SigmaIetaIeta_SB_Ph70To80");
h9a->Scale((2.99815997E10*7.4*1000)/9904852);

TH1F* h10a = (TH1F*)file2->Get("h_SigmaIetaIeta_SB_Ph70To80");
h10a->Scale((2.06583E7*7.4*1000)/9755606);

cout << "h9a + h10a entries = " << long(h9a->Integral() + h10a->Integral())  << endl;

h9a->Add(h10a);

cout << "h9a entries = " << long(h9a->Integral()) << endl;

h9a->SetName("h_SigmaIetaIeta_SB_Ph70To80");

TFile *f = new TFile("DataFile.root","recreate");
h1->Write();
h1a->Write();
h3->Write();
h3a->Write();
h5->Write();
h5a->Write();
h7->Write();
h7a->Write();
h9->Write();
h9a->Write();
f->Close();

}
