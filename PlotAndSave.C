#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TChain.h>
#include <TLorentzVector.h>
#include <TLegend.h>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TFrame.h>


void PlotAndSave(std::string infile){

  gROOT->SetStyle("Plain");
  
  TCanvas *c1 = new TCanvas("c1", "c1", 600, 400);
  TFile* file1 = TFile::Open((infile+".root").c_str());
  TH2F *h1 = (TH2F*)gDirectory->Get("h_InvariantMass_Mu");
  h1->Rebin(10);
  h1->Draw();
  c1->SaveAs((infile+"_InvariantMass_Mu.pdf").c_str());
  c1->SaveAs((infile+"_InvariantMass_Mu.pdf").c_str());

  TH2F *h2 = (TH2F*)gDirectory->Get("h_InvariantMass_El");
  h2->Rebin(10);
  h2->Draw();
  c1->SaveAs((infile+"_InvariantMass_El.pdf").c_str());
  c1->SaveAs((infile+"_InvariantMass_El.pdf").c_str());

  TH2F *h3 = (TH2F*)gDirectory->Get("h_mu_pt_leading");
  h3->Rebin(5);
  h3->GetXaxis()->SetRangeUser(0, 200);
  h3->Draw();
  c1->SaveAs((infile+"_mu_pt_leading.pdf").c_str());
  c1->SaveAs((infile+"_mu_pt_leading.pdf").c_str());
  
  TH2F *h4 = (TH2F*)gDirectory->Get("h_mu_pt_trailing");
  h4->Rebin(5);
  h4->GetXaxis()->SetRangeUser(0, 200);
  h4->Draw();
  c1->SaveAs((infile+"_mu_pt_trailing.pdf").c_str());
  c1->SaveAs((infile+"_mu_pt_trailing.pdf").c_str());

  TH2F *h5 = (TH2F*)gDirectory->Get("h_mu_eta_leading");
  h5->Rebin(10);
  h5->Draw();
  c1->SaveAs((infile+"_mu_eta_leading.pdf").c_str());
  c1->SaveAs((infile+"_mu_eta_leading.pdf").c_str());
  
  TH2F *h6 = (TH2F*)gDirectory->Get("h_mu_eta_trailing");
  h6->Rebin(10);
  h6->Draw();
  c1->SaveAs((infile+"_mu_eta_trailing.pdf").c_str());
  c1->SaveAs((infile+"_mu_eta_trailing.pdf").c_str());

  TH2F *h7 = (TH2F*)gDirectory->Get("h_mu_phi_leading");
  h7->Rebin(10);
  h7->Draw();
  c1->SaveAs((infile+"_mu_phi_leading.pdf").c_str());
  c1->SaveAs((infile+"_mu_phi_leading.pdf").c_str());
  
  TH2F *h8 = (TH2F*)gDirectory->Get("h_mu_phi_trailing");
  h8->Rebin(10);
  h8->Draw();
  c1->SaveAs((infile+"_mu_phi_trailing.pdf").c_str());
  c1->SaveAs((infile+"_mu_phi_trailing.pdf").c_str());
   
  TH2F *h9 = (TH2F*)gDirectory->Get("h_mu_energy_leading");
  h9->Rebin(5);
  h9->GetXaxis()->SetRangeUser(0, 500);
  h9->Draw();
  c1->SaveAs((infile+"_mu_energy_leading.pdf").c_str());
  c1->SaveAs((infile+"_mu_energy_leading.pdf").c_str());
  
  TH2F *h10 = (TH2F*)gDirectory->Get("h_mu_energy_trailing");
  h10->Rebin(5);
  h10->GetXaxis()->SetRangeUser(0, 500);
  h10->Draw();
  c1->SaveAs((infile+"_mu_energy_trailing.pdf").c_str());
  c1->SaveAs((infile+"_mu_energy_trailing.pdf").c_str());

  TH2F *h11 = (TH2F*)gDirectory->Get("h_el_pt_leading"); 
  h11->Rebin(5);
  h11->GetXaxis()->SetRangeUser(0, 200);
  h11->Draw(); 
  c1->SaveAs((infile+"_el_pt_leading.pdf").c_str()); 
  c1->SaveAs((infile+"_el_pt_leading.pdf").c_str()); 
   
  TH2F *h12 = (TH2F*)gDirectory->Get("h_el_pt_trailing"); 
  h12->Rebin(5); 
  h12->GetXaxis()->SetRangeUser(0, 200);
  h12->Draw(); 
  c1->SaveAs((infile+"_el_pt_trailing.pdf").c_str()); 
  c1->SaveAs((infile+"_el_pt_trailing.pdf").c_str()); 
 
  TH2F *h13 = (TH2F*)gDirectory->Get("h_el_eta_leading"); 
  h13->Rebin(10); 
  h13->Draw(); 
  c1->SaveAs((infile+"_el_eta_leading.pdf").c_str()); 
  c1->SaveAs((infile+"_el_eta_leading.pdf").c_str()); 
   
  TH2F *h14 = (TH2F*)gDirectory->Get("h_el_eta_trailing"); 
  h14->Rebin(10); 
  h14->Draw(); 
  c1->SaveAs((infile+"_el_eta_trailing.pdf").c_str()); 
  c1->SaveAs((infile+"_el_eta_trailing.pdf").c_str()); 
 
  TH2F *h15 = (TH2F*)gDirectory->Get("h_el_phi_leading"); 
  h15->Rebin(10); 
  h15->Draw(); 
  c1->SaveAs((infile+"_el_phi_leading.pdf").c_str()); 
  c1->SaveAs((infile+"_el_phi_leading.pdf").c_str()); 
   
  TH2F *h16 = (TH2F*)gDirectory->Get("h_el_phi_trailing"); 
  h16->Rebin(10); 
  h16->Draw(); 
  c1->SaveAs((infile+"_el_phi_trailing.pdf").c_str()); 
  c1->SaveAs((infile+"_el_phi_trailing.pdf").c_str()); 
    
  TH2F *h17 = (TH2F*)gDirectory->Get("h_el_energy_leading"); 
  h17->Rebin(5); 
  h17->GetXaxis()->SetRangeUser(0, 500);
  h17->Draw(); 
  c1->SaveAs((infile+"_el_energy_leading.pdf").c_str()); 
  c1->SaveAs((infile+"_el_energy_leading.pdf").c_str()); 
   
  TH2F *h18 = (TH2F*)gDirectory->Get("h_el_energy_trailing");
  h18->Rebin(5); 
  h18->GetXaxis()->SetRangeUser(0, 500);
  h18->Draw();
  c1->SaveAs((infile+"_el_energy_trailing.pdf").c_str()); 
  c1->SaveAs((infile+"_el_energy_trailing.pdf").c_str());  
 
  TH2F *h19 = (TH2F*)gDirectory->Get("h_photon_pt"); 
  h19->GetXaxis()->SetRangeUser(0, 500);
  h19->Rebin(5);
  h19->Draw();
  c1->SaveAs((infile+"_photon_pt.pdf").c_str()); 
  c1->SaveAs((infile+"_photon_pt.pdf").c_str()); 
   
  TH2F *h20 = (TH2F*)gDirectory->Get("h_photon_eta");          
  h20->Rebin(10);
  h20->Draw();
  c1->SaveAs((infile+"_photon_eta.pdf").c_str());         
  c1->SaveAs((infile+"_photon_eta.pdf").c_str());      
  
  TH2F *h21 = (TH2F*)gDirectory->Get("h_photon_phi");          
  h21->Rebin(10);
  h21->Draw();
  c1->SaveAs((infile+"_photon_phi.pdf").c_str());         
  c1->SaveAs((infile+"_photon_phi.pdf").c_str());    

  TH2F *h22 = (TH2F*)gDirectory->Get("h_photon_energy");          
  h22->Rebin(10);
  h22->Draw();
  c1->SaveAs((infile+"_photon_energy.pdf").c_str());         
  c1->SaveAs((infile+"_photon_energy.pdf").c_str());    

  c1->SetLogy();

  TH2F *h23 = (TH2F*)gDirectory->Get("h_DeltaR_el1_ph1");        
  h23->Rebin(10);
  h23->Draw();
  c1->SaveAs((infile+"_DeltaR_el1_ph1.pdf").c_str());       
  c1->SaveAs((infile+"_DeltaR_el1_ph1.pdf").c_str());  

  TH2F *h24 = (TH2F*)gDirectory->Get("h_DeltaR_el2_ph1");
  h24->Rebin(10);
  h24->Draw();
  c1->SaveAs((infile+"_DeltaR_el2_ph1.pdf").c_str());
  c1->SaveAs((infile+"_DeltaR_el2_ph1.pdf").c_str());



/*
  KEY: TH1F	h_el_isolation_leading;1	Leading electron Isolation
  KEY: TH1F	h_el_isolation_trailing;1	Trailing electron Isolation
  KEY: TH1F	h_photon_energy;1	Photon Energy
  KEY: TH1F	h_ph_chIsolation_leading;1	Leading Photon Ch Isolation
  KEY: TH1F	h_ph_nuIsolation_leading;1	Leading Photon Nu Isolation
  KEY: TH1F	h_ph_phIsolation_leading;1	Leading Photon Ph Isolation
  KEY: TH1F	h_InvariantMass_ElPh;1	Di-electron and photon invariant mass
  KEY: TH1F	h_InvariantMass_MuPh;1	Di-muon and photon invariant mass
  KEY: TH2F	h_Mmumu_MmumuGamma;1	Scatter Plot of  M_{#mu#mu} versus M_{#mu#mu#gamma}
  KEY: TH2F	h_Mee_MeeGamma;1	Scatter Plot of  M_{ee} versus M_{ee#gamma}
  KEY: TH1F	h_Difference_El;1	Difference pt
  KEY: TH1F	h_Difference_Mu;1	Difference pt 
*/



}
