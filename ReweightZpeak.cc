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

void ReweightZpeak(std::string infileData, std::string infileMC, std::string outfile){

 TFile* fileData = TFile::Open((infileData+".root").c_str());
 TFile* fileMC = TFile::Open((infileMC+".root").c_str());

 TH1F *h_InvariantMass_MuMu_Data = (TH1F*) fileData->Get("h_InvariantMass_MuMu");
 h_InvariantMass_MuMu_Data->Rebin(10);
 TH1F *h_InvariantMass_MuMu_MC = (TH1F*) fileMC->Get("h_InvariantMass_MuMu");
 h_InvariantMass_MuMu_MC->Rebin(10);
 h_InvariantMass_MuMu_MC->Scale((1.2*132.6*7.4*1000)/6583032);

 TH1F *h_InvariantMass_MuMu_DataDist=(TH1F*)h_InvariantMass_MuMu_Data->Clone("h_InvariantMass_MuMu_DataDist");
 h_InvariantMass_MuMu_DataDist->Divide(h_InvariantMass_MuMu_MC);

 h_InvariantMass_MuMu_DataDist->SetTitle("Weights in the Z peak");
 h_InvariantMass_MuMu_DataDist->SetStats(kFALSE);

 TF1 *fit_ratio = new TF1("fit_ratio","[0]",0.0, 300.0);
 fit_ratio->SetParLimits(0,0.0, 1.5);
 fit_ratio->SetLineColor(kRed);
 fit_ratio->SetLineWidth(3);

 h_InvariantMass_MuMu_DataDist->Fit("fit_ratio", "", "", 0.0,300.0);

 std::string histfilename=(outfile+".root").c_str();
 TFile *tFile=new TFile(histfilename.c_str(), "RECREATE");
 h_InvariantMass_MuMu_DataDist->Write();
 fit_ratio->Write("fit_ratio");
 tFile->Close();

}
