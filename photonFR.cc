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

void photonFR(std::string infile, std::string outfile){


 TCanvas *c1 = new TCanvas("c1_fit1","The Fit Canvas",200,10,700,500);
 TFile* file = TFile::Open((infile+".root").c_str());
 
 TH1F *h_photonpT_N = (TH1F*) file->Get("h_photonpT_Num");
 TH1F *h_photonpT_D = (TH1F*) file->Get("h_photonpT_Den");
 Double_t rebin_array_X[8]={30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 100.0, 130.0};
 TH1F *h_photonpT_N_new = (TH1F*)h_photonpT_N->Rebin(7, "h_photonpT_N_new", rebin_array_X);
 h_photonpT_N_new->SetBinError(6, TMath::Sqrt(h_photonpT_N->GetBinContent(6))/2.0);
 h_photonpT_N_new->SetBinContent(6, h_photonpT_N->GetBinContent(6)/2.0);
 h_photonpT_N_new->SetBinError(7, TMath::Sqrt(h_photonpT_N->GetBinContent(7))/3.0);
 h_photonpT_N_new->SetBinContent(7, h_photonpT_N->GetBinContent(7)/3.0);
 TH1F *h_photonpT_D_new = (TH1F*)h_photonpT_D->Rebin(7, "h_photonpT_D_new", rebin_array_X);
 h_photonpT_D_new->SetBinError(6, TMath::Sqrt(h_photonpT_D->GetBinContent(6))/2.0);
 h_photonpT_D_new->SetBinContent(6, h_photonpT_D->GetBinContent(6)/2.0);
 h_photonpT_D_new->SetBinError(7, TMath::Sqrt(h_photonpT_D->GetBinContent(7))/3.0);
 h_photonpT_D_new->SetBinContent(7, h_photonpT_D->GetBinContent(7)/3.0);
 TH1F *h_photonpT=(TH1F*)h_photonpT_N_new->Clone("h_photonpT");
 h_photonpT->Divide(h_photonpT_D_new);
 h_photonpT->SetTitle("Photon pT versus fake ratio");
 h_photonpT->GetYaxis()->SetTitle("Fake Ratio");
 h_photonpT->SetStats(kFALSE);
 h_photonpT->Draw("HIST");

 std::string histfilename=(outfile+".root").c_str();
 TFile *tFile=new TFile(histfilename.c_str(), "RECREATE");
 h_photonpT->Write();
 tFile->Close();

}
