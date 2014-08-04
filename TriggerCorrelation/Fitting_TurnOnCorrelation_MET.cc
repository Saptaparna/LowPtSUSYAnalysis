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


void Fitting_TurnOnCorrelation(std::string infile, std::string outfile){

  gStyle->SetErrorX(0);
  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat("4.1f m");

  TFile* file = TFile::Open((infile+".root").c_str());
  TH2F *h_ControlTrigger_1 = (TH2F*) file->Get("h_ControlTrigger1");
  Double_t rebin_array_1_X[7]={30.0, 50.0, 70.0, 90.0, 110.0, 200.0, 600.0};
  Double_t rebin_array_1_Y[5]={0.0,  20.0, 40.0, 80.0, 100.0}; //, 200.0, 600.0};
  TH2F *h_ControlTrigger1_Rebin = new TH2F("h_ControlTrigger1_Rebin",h_ControlTrigger_1->GetTitle(), 6, rebin_array_1_X, 4, rebin_array_1_Y);
  TAxis *xaxis = h_ControlTrigger_1->GetXaxis();
  TAxis *yaxis = h_ControlTrigger_1->GetYaxis();
  for (int j=1; j<=yaxis->GetNbins();j++) {
      for (int i=1; i<=xaxis->GetNbins();i++) {
         h_ControlTrigger1_Rebin->Fill(xaxis->GetBinCenter(i),yaxis->GetBinCenter(j),h_ControlTrigger_1->GetBinContent(i,j));
      }
  }

  h_ControlTrigger1_Rebin->GetXaxis()->SetTitle("photon pT [GeV]");
  h_ControlTrigger1_Rebin->GetYaxis()->SetTitle("MET [GeV]");  
  h_ControlTrigger1_Rebin->Draw("COLZ");

  TH2F *h_AnalysisTrigger_1 = (TH2F*) file->Get("h_AnalysisTrigger1");
  TH2F *h_AnalysisTrigger1_Rebin = new TH2F("h_AnalysisTrigger1_Rebin",h_AnalysisTrigger_1->GetTitle(), 6, rebin_array_1_X, 4, rebin_array_1_Y);
  TAxis *xaxis = h_AnalysisTrigger_1->GetXaxis();
  TAxis *yaxis = h_AnalysisTrigger_1->GetYaxis();
  for (int j=1; j<=yaxis->GetNbins();j++) {
      for (int i=1; i<=xaxis->GetNbins();i++) {
         h_AnalysisTrigger1_Rebin->Fill(xaxis->GetBinCenter(i),yaxis->GetBinCenter(j),h_AnalysisTrigger_1->GetBinContent(i,j));
      }
  }

  h_AnalysisTrigger1_Rebin->GetXaxis()->SetTitle("photon pT [GeV]");
  h_AnalysisTrigger1_Rebin->GetYaxis()->SetTitle("MET [GeV]");    
  h_ControlTrigger1_Rebin->Draw("COLZ");

  TH2F *h_ControlTrigger_1_RemoveJets = (TH2F*) file->Get("h_ControlTrigger1_RemoveJets");
  TH2F *h_ControlTrigger1_RemoveJets_Rebin = new TH2F("h_ControlTrigger1_RemoveJets_Rebin",h_ControlTrigger_1_RemoveJets->GetTitle(), 6, rebin_array_1_X, 4, rebin_array_1_Y);
  TAxis *xaxis = h_ControlTrigger_1_RemoveJets->GetXaxis();
  TAxis *yaxis = h_ControlTrigger_1_RemoveJets->GetYaxis();
  for (int j=1; j<=yaxis->GetNbins();j++) {
      for (int i=1; i<=xaxis->GetNbins();i++) {
         h_ControlTrigger1_RemoveJets_Rebin->Fill(xaxis->GetBinCenter(i),yaxis->GetBinCenter(j),h_ControlTrigger_1_RemoveJets->GetBinContent(i,j));
      }
  }

  h_ControlTrigger1_RemoveJets_Rebin->GetXaxis()->SetTitle("photon pT [GeV]");
  h_ControlTrigger1_RemoveJets_Rebin->GetYaxis()->SetTitle("MET [GeV]");
  h_ControlTrigger1_RemoveJets_Rebin->Draw("COLZ");

  TH2F *h_AnalysisTrigger_1_RemoveJets = (TH2F*) file->Get("h_AnalysisTrigger1_RemoveJets");
  TH2F *h_AnalysisTrigger1_RemoveJets_Rebin = new TH2F("h_AnalysisTrigger1_RemoveJets_Rebin",h_AnalysisTrigger_1_RemoveJets->GetTitle(), 6, rebin_array_1_X, 4, rebin_array_1_Y);
  TAxis *xaxis = h_AnalysisTrigger_1_RemoveJets->GetXaxis();
  TAxis *yaxis = h_AnalysisTrigger_1_RemoveJets->GetYaxis();
  for (int j=1; j<=yaxis->GetNbins();j++) {
      for (int i=1; i<=xaxis->GetNbins();i++) {
         h_AnalysisTrigger1_RemoveJets_Rebin->Fill(xaxis->GetBinCenter(i),yaxis->GetBinCenter(j),h_AnalysisTrigger_1_RemoveJets->GetBinContent(i,j));
      }
  }

  h_AnalysisTrigger1_RemoveJets_Rebin->GetXaxis()->SetTitle("photon pT [GeV]");
  h_AnalysisTrigger1_RemoveJets_Rebin->GetYaxis()->SetTitle("MET [GeV]");
  h_ControlTrigger1_RemoveJets_Rebin->Draw("COLZ");

  TH2F *h_Eff1_2D = new TH2F("h_Eff1_2D", "Efficiency 1: Analysis Trigger versus Control Trigger", 6, rebin_array_1_X, 4, rebin_array_1_Y); 
  h_Eff1_2D->Divide(h_AnalysisTrigger1_Rebin, h_ControlTrigger1_Rebin, 1., 1., "b");
  h_Eff1_2D->GetXaxis()->SetTitle("photon pT [GeV]");
  h_Eff1_2D->GetYaxis()->SetTitle("MET [GeV]");  
  h_Eff1_2D->Draw("COLZ");

  TH2F *h_Eff1_2D_RemoveJets = new TH2F("h_Eff1_2D_RemoveJets", "Efficiency 1: Analysis Trigger versus Control Trigger", 6, rebin_array_1_X, 4, rebin_array_1_Y);
  h_Eff1_2D_RemoveJets->Divide(h_AnalysisTrigger1_RemoveJets_Rebin, h_ControlTrigger1_RemoveJets_Rebin, 1., 1., "b");
  h_Eff1_2D_RemoveJets->GetXaxis()->SetTitle("photon pT [GeV]");
  h_Eff1_2D_RemoveJets->GetYaxis()->SetTitle("MET [GeV]");
  h_Eff1_2D_RemoveJets->Draw("COLZ");

  std::string histfilename=(outfile+".root").c_str();
  TFile *tFile=new TFile(histfilename.c_str(), "RECREATE");
  h_AnalysisTrigger_1->Write();
  h_ControlTrigger_1->Write();
  h_AnalysisTrigger1_Rebin->Write();
  h_ControlTrigger1_Rebin->Write();
  h_Eff1_2D->Write();
  h_AnalysisTrigger_1_RemoveJets->Write();
  h_ControlTrigger_1_RemoveJets->Write();
  h_AnalysisTrigger1_RemoveJets_Rebin->Write();
  h_ControlTrigger1_RemoveJets_Rebin->Write();
  h_Eff1_2D_RemoveJets->Write();
  tFile->Close();
}
