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


void Fitting_TurnOn(std::string infile, std::string outfile){

  TCanvas *c1 = new TCanvas("c1_fit1","The Fit Canvas",200,10,700,500);

  gStyle->SetErrorX(0);

  TFile* file = TFile::Open((infile+".root").c_str());
  TH1F *h_MET_HLTPhoIdMet = (TH1F*) file->Get("h_MET_HLTPhoIdMet");
  h_MET_HLTPhoIdMet->Rebin(3);
  TH1F *h_MET_HLTPhoId = (TH1F*) file->Get("h_MET_HLTPhoId");
  h_MET_HLTPhoId->Rebin(3);
  TGraphAsymmErrors *Eff_MET = new TGraphAsymmErrors;
  Eff_MET->BayesDivide(h_MET_HLTPhoIdMet, h_MET_HLTPhoId, "");
  Eff_MET->SetLineColor(kBlue);
  Eff_MET->GetXaxis()->SetRangeUser(25, 800);
  Eff_MET->GetYaxis()->SetTitle("Efficiency");
  Eff_MET->GetXaxis()->SetTitle("Offline MET");
  Eff_MET->SetMinimum(0.0);
  Eff_MET->SetTitle("Efficiency: HLT_Photon30_R9Id90_Met25 w.r.t HLT_Photon30_R9Id90");

  Eff_MET->Draw("AP");
  
  TF1 *fit_MET = new TF1("fit_MET","0.5*[2]*(1+TMath::Erf((x-[0])/[1]))",25.0, 800.0); 
  fit_MET->SetParameter(0, 20);
  fit_MET->SetParameter(1, 40);
  fit_MET->SetParameter(2, 0.8);
  fit_MET->SetParLimits(2,0.0,1.0);
  Eff_MET->Fit("fit_MET", "", "", 25.0,800.0); 

  fit_MET->SetLineColor(kRed);
  fit_MET->SetLineWidth(2);

  TH1F *h_MET_HLTPho = (TH1F*) file->Get("h_PHO_HLTPho");
  h_MET_HLTPho->Rebin(3);
  TH1F *h_MET_HLTPhoId = (TH1F*) file->Get("h_PHO_HLTPhoId");
  h_MET_HLTPhoId->Rebin(3);
  TGraphAsymmErrors *Eff_PHO = new TGraphAsymmErrors;
  Eff_PHO->BayesDivide(h_PHO_HLTPhoId, h_PHO_HLTPho, "");
  Eff_PHO->SetLineColor(kBlue);
  Eff_PHO->GetXaxis()->SetRangeUser(0, 800);
  Eff_PHO->GetYaxis()->SetTitle("Efficiency");
  Eff_PHO->GetXaxis()->SetTitle("Offline Photon pT [GeV]");
  Eff_PHO->SetMinimum(0.0);
  Eff_PHO->SetTitle("Efficiency: HLT_Photon30_R9Id90 w.r.t HLT_Photon30");

  Eff_PHO->Draw("AP E");

  TF1 *fit_PHO = new TF1("fit_PHO","0.5*[2]*(1+TMath::Erf((x-[0])/[1]))",0.0, 800.0);
  fit_PHO->SetParameter(0, 20);
  fit_PHO->SetParameter(1, 40);
  fit_PHO->SetParameter(2, 0.8);
  fit_PHO->SetParLimits(2,0.0,1.0);
  Eff_PHO->Fit("fit_PHO", "", "", 0,800.0);

  fit_PHO->SetLineColor(kRed);
  fit_PHO->SetLineWidth(2);

  std::string histfilename=(outfile+".root").c_str();
  TFile *tFile=new TFile(histfilename.c_str(), "RECREATE");
  Eff_MET->Write("Eff_MET");
  fit_MET->Write("fit_MET");
  Eff_PHO->Write("Eff_PHO");
  fit_PHO->Write("fit_PHO");
  tFile->Close();
}
