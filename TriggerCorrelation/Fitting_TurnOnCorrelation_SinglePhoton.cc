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


void Fitting_TurnOnCorrelation_SinglePhoton(std::string infile, std::string outfile){

  gStyle->SetErrorX(0);
  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat("4.1f m");

  TFile* file = TFile::Open((infile+".root").c_str());
  TH2F *h_ControlTrigger_2 = (TH2F*) file->Get("h_ControlTrigger2");
  Double_t rebin_array_2_X[13]={30.0, 50.0, 70.0, 90.0, 110.0, 130.0, 150.0, 170.0, 190.0, 210.0, 250.0, 300.0, 600.0};
  Double_t rebin_array_2_Y[10]={0.0,  30.0, 50.0, 70.0, 90.0,  110.0, 130.0, 150.0, 200.0, 600.0};
  TH2F *h_ControlTrigger2_Rebin = new TH2F("h_ControlTrigger2_Rebin",h_ControlTrigger_2->GetTitle(), 12, rebin_array_2_X, 9, rebin_array_2_Y);
  TAxis *xaxis = h_ControlTrigger_2->GetXaxis();
  TAxis *yaxis = h_ControlTrigger_2->GetYaxis();
  for (int j=1; j<=yaxis->GetNbins();j++) {
      for (int i=1; i<=xaxis->GetNbins();i++) {
         h_ControlTrigger2_Rebin->Fill(xaxis->GetBinCenter(i),yaxis->GetBinCenter(j),h_ControlTrigger_2->GetBinContent(i,j));
      }
  }

  h_ControlTrigger2_Rebin->GetXaxis()->SetTitle("photon pT [GeV]");
  h_ControlTrigger2_Rebin->GetYaxis()->SetTitle("MET [GeV]");

  TH2F *h_AnalysisTrigger_2 = (TH2F*) file->Get("h_AnalysisTrigger2");
  TH2F *h_AnalysisTrigger2_Rebin = new TH2F("h_AnalysisTrigger2_Rebin",h_AnalysisTrigger_2->GetTitle(), 12, rebin_array_2_X, 9, rebin_array_2_Y);
  TAxis *xaxis = h_AnalysisTrigger_2->GetXaxis();
  TAxis *yaxis = h_AnalysisTrigger_2->GetYaxis();
  for (int j=1; j<=yaxis->GetNbins();j++) {
      for (int i=1; i<=xaxis->GetNbins();i++) {
         h_AnalysisTrigger2_Rebin->Fill(xaxis->GetBinCenter(i),yaxis->GetBinCenter(j),h_AnalysisTrigger_2->GetBinContent(i,j));
      }
  }

  h_AnalysisTrigger2_Rebin->GetXaxis()->SetTitle("photon pT [GeV]");
  h_AnalysisTrigger2_Rebin->GetYaxis()->SetTitle("MET [GeV]");
 
  TH2F *h_Eff2_2D = new TH2F("h_Eff2_2D", "Efficiency 2: Analysis Trigger versus Control Trigger", 12, rebin_array_2_X, 9, rebin_array_2_Y);
  h_Eff2_2D->Divide(h_AnalysisTrigger2_Rebin, h_ControlTrigger2_Rebin, 1., 1., "b");
  h_Eff2_2D->GetXaxis()->SetTitle("photon pT [GeV]");
  h_Eff2_2D->GetYaxis()->SetTitle("MET [GeV]");

  std::string histfilename=(outfile+".root").c_str();
  TFile *tFile=new TFile(histfilename.c_str(), "RECREATE");
  h_AnalysisTrigger_2->Write();
  h_ControlTrigger_2->Write();
  h_AnalysisTrigger2_Rebin->Write();
  h_ControlTrigger2_Rebin->Write();
  h_Eff2_2D->Write();
  tFile->Close();
}
