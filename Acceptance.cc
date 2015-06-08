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
#include <TGraphErrors.h>

void Acceptance(std::string outfile){

  TCanvas *c1 = new TCanvas("c1_fit1","The Fit Canvas",200,10,700,500);
  //double FR[7] = {3.49656, 2.63871, 2.18338, 1.9606, 1.82274, 1.71089, 1.51592};
  //double FR[7] = {2.46412, 1.79744,  1.42782, 1.26439, 1.12921, 1.10087, 1.09296};
  //double FR[7] = {2.71755, 2.06168, 1.69469, 1.52471, 1.41677, 1.32065, 1.15989};
  //double FR[7] = {4.002, 2.710, 2.212, 2.044, 1.919, 1.827, 1.660};
  //double FR[7] = {4.92, 3.68, 3.07, 2.77, 2.60, 2.43, 2.15};
  //double FR[7] = {4.8156, 3.60663, 2.96109, 2.69444, 2.49084, 2.34777, 2.09628};//With spike cleaning
  //double phpT[7] = {35.0, 45.0, 55.0, 65.0, 75.0, 90.0, 115.0};
  //double FR[7] = {4.66094, 3.48385, 2.85518, 2.59734, 2.39063, 2.2498, 2.00732};//for Wolfgang
  double FR[11] = {0.29, 0.46, 0.57, 0.65, 0.71, 0.75, 0.79, 0.84, 0.85, 0.87, 0.87};
  double phpT[11] = {30, 50, 70, 90, 110, 130, 150, 210, 230, 370, 710};

  gStyle->SetErrorX(0.5);

  TGraph *phpT_FR = new TGraph(11, phpT, FR);
  phpT_FR->SetTitle("");
  phpT_FR->SetLineColor(kBlue+4);
  phpT_FR->SetMarkerColor(kBlue);
  phpT_FR->SetMarkerStyle(20);
  phpT_FR->SetMarkerSize(1.0);
  phpT_FR->SetLineWidth(2);
  phpT_FR->SetLineStyle(1);
  phpT_FR->SetTitle("Signal Acceptance");
  phpT_FR->GetXaxis()->SetTitle("HT < HT_{min} [GeV]");
  phpT_FR->GetYaxis()->SetTitle("Signal Acceptance in %");
  phpT_FR->GetYaxis()->SetLabelSize(0.035);
  phpT_FR->GetXaxis()->SetLabelSize(0.035);
  phpT_FR->GetYaxis()->SetTitleSize(0.035);
  phpT_FR->GetXaxis()->SetTitleSize(0.035);
  phpT_FR->GetXaxis()->SetLabelFont(62);
  phpT_FR->GetYaxis()->SetLabelFont(62);
  phpT_FR->GetXaxis()->SetTitleFont(62);
  phpT_FR->GetYaxis()->SetTitleFont(62);
  phpT_FR->GetYaxis()->SetTitleOffset(1.4);
  phpT_FR->SetMaximum(1.0);
  //phpT_FR->SetMinimum(1.0);
  phpT_FR->Draw("APL");

  std::string histfilename=(outfile+".root").c_str();
  TFile *tFile=new TFile(histfilename.c_str(), "RECREATE");
  phpT_FR->Write("phpT_FR");
  tFile->Close();

}

void Acceptance_Mll(std::string outfile){

  TCanvas *c1 = new TCanvas("c1_fit1","The Fit Canvas",200,10,700,500);
  //double FR[7] = {3.49656, 2.63871, 2.18338, 1.9606, 1.82274, 1.71089, 1.51592};
  //double FR[7] = {2.46412, 1.79744,  1.42782, 1.26439, 1.12921, 1.10087, 1.09296};
  //double FR[7] = {2.71755, 2.06168, 1.69469, 1.52471, 1.41677, 1.32065, 1.15989};
  //double FR[7] = {4.002, 2.710, 2.212, 2.044, 1.919, 1.827, 1.660};
  //double FR[7] = {4.92, 3.68, 3.07, 2.77, 2.60, 2.43, 2.15};
  //double FR[7] = {4.8156, 3.60663, 2.96109, 2.69444, 2.49084, 2.34777, 2.09628};//With spike cleaning
  //double phpT[7] = {35.0, 45.0, 55.0, 65.0, 75.0, 90.0, 115.0};
  //double FR[7] = {4.66094, 3.48385, 2.85518, 2.59734, 2.39063, 2.2498, 2.00732};//for Wolfgang
  double FR[9] = {0.04, 0.20, 0.43, 0.58, 0.69, 0.80, 0.84, 0.87, 0.87};
  double phpT[9] = {20, 40, 60, 80, 100, 140, 180, 250, 270};

  gStyle->SetErrorX(0.5);

  TGraph *phpT_FR = new TGraph(9, phpT, FR);
  phpT_FR->SetTitle("");
  phpT_FR->SetLineColor(kBlue+4);
  phpT_FR->SetMarkerColor(kBlue);
  phpT_FR->SetMarkerStyle(20);
  phpT_FR->SetMarkerSize(1.0);
  phpT_FR->SetLineWidth(2);
  phpT_FR->SetLineStyle(1);
  phpT_FR->SetTitle("Signal Acceptance");
  phpT_FR->GetXaxis()->SetTitle("Invariant Mass (M_{ll}) < Invariant Mass min (M_{llmin}) [GeV]");
  phpT_FR->GetYaxis()->SetTitle("Signal Acceptance in %");
  phpT_FR->GetYaxis()->SetLabelSize(0.035);
  phpT_FR->GetXaxis()->SetLabelSize(0.035);
  phpT_FR->GetYaxis()->SetTitleSize(0.035);
  phpT_FR->GetXaxis()->SetTitleSize(0.035);
  phpT_FR->GetXaxis()->SetLabelFont(62);
  phpT_FR->GetYaxis()->SetLabelFont(62);
  phpT_FR->GetXaxis()->SetTitleFont(62);
  phpT_FR->GetYaxis()->SetTitleFont(62);
  phpT_FR->GetYaxis()->SetTitleOffset(1.4);
  phpT_FR->SetMaximum(1.0);
  //phpT_FR->SetMinimum(1.0);
  phpT_FR->Draw("APL");

  std::string histfilename=(outfile+".root").c_str();
  TFile *tFile=new TFile(histfilename.c_str(), "RECREATE");
  phpT_FR->Write("phpT_FR");
  tFile->Close();

}

