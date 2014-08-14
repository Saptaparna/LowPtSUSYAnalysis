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
  TH1F *h_caloMET_HLTPhoIdMet = (TH1F*) file->Get("h_caloMET_HLTPhoIdMet");
  h_caloMET_HLTPhoIdMet->Rebin(3);
  TH1F *h_caloMET_HLTPhoId = (TH1F*) file->Get("h_caloMET_HLTPhoId");
  h_caloMET_HLTPhoId->Rebin(3);
  TGraphAsymmErrors *Eff_caloMET = new TGraphAsymmErrors;
  Eff_caloMET->BayesDivide(h_caloMET_HLTPhoIdMet, h_caloMET_HLTPhoId, "");
  Eff_caloMET->SetLineColor(kBlue);
  Eff_caloMET->GetXaxis()->SetRangeUser(25, 800);
  Eff_caloMET->GetYaxis()->SetTitle("Efficiency");
  Eff_caloMET->GetXaxis()->SetTitle("Offline caloMET");
  Eff_caloMET->SetMinimum(0.0);
  Eff_caloMET->SetTitle("Efficiency: HLT_Photon30_R9Id90_Met25 w.r.t HLT_Photon30_R9Id90");

  Eff_caloMET->Draw("AP");
  
  TF1 *fit_caloMET = new TF1("fit_caloMET","0.5*[2]*(1+TMath::Erf((x-[0])/[1]))",25.0, 800.0); 
  fit_caloMET->SetParameter(0, 18);
  fit_caloMET->SetParameter(1, 34);
  fit_caloMET->SetParameter(2, 0.9);
  fit_caloMET->SetParLimits(2,0.0,1.0);
  Eff_caloMET->Fit("fit_caloMET", "", "", 25.0,800.0); 

  fit_caloMET->SetLineColor(kRed);
  fit_caloMET->SetLineWidth(3.5);

  std::string histfilename=(outfile+".root").c_str();
  TFile *tFile=new TFile(histfilename.c_str(), "RECREATE");
  Eff_caloMET->Write("Eff_caloMET");
  fit_caloMET->Write("fit_caloMET");
  tFile->Close();
}
