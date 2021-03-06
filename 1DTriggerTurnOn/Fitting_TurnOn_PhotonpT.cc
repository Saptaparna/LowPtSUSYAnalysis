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


void Fitting_TurnOn_METTrigger(std::string infile, std::string outfile){

  TCanvas *c1 = new TCanvas("c1_fit1","The Fit Canvas",200,10,700,500);

  gStyle->SetErrorX(0);

  TFile* file = TFile::Open((infile+".root").c_str());

  TH1F *h_MET_HLTPho = (TH1F*) file->Get("h_PHO_HLT_METTrig");
  h_MET_HLTPho->Rebin(5);
  TH1F *h_MET_HLTPhoId = (TH1F*) file->Get("h_PHO_HLTPho_METTrig");
  h_MET_HLTPhoId->Rebin(5);
  TGraphAsymmErrors *Eff_PHO = new TGraphAsymmErrors;
  Eff_PHO->BayesDivide(h_MET_HLTPhoId, h_MET_HLTPho, "");
  Eff_PHO->SetLineColor(kBlue);
  Eff_PHO->GetXaxis()->SetRangeUser(0, 800);
  Eff_PHO->GetYaxis()->SetTitle("Efficiency");
  Eff_PHO->GetXaxis()->SetTitle("Offline Photon pT [GeV]");
  Eff_PHO->SetMinimum(0.0);
  Eff_PHO->SetTitle("Efficiency: HLT_Photon30_R9Id90_Met25 w.r.t HLT_Met100");

  Eff_PHO->Draw("AP E");

  TH1F *h_MET_HLTPho_RemoveJets = (TH1F*) file->Get("h_PHO_HLT_METTrig_RemoveJets");
  h_MET_HLTPho_RemoveJets->Rebin(5);
  TH1F *h_MET_HLTPhoId_RemoveJets = (TH1F*) file->Get("h_PHO_HLTPho_METTrig_RemoveJets");
  h_MET_HLTPhoId_RemoveJets->Rebin(5);
  TGraphAsymmErrors *Eff_PHO_RemoveJets = new TGraphAsymmErrors;
  Eff_PHO_RemoveJets->BayesDivide(h_MET_HLTPhoId_RemoveJets, h_MET_HLTPho_RemoveJets, "");
  Eff_PHO_RemoveJets->SetLineColor(kBlue);
  Eff_PHO_RemoveJets->GetXaxis()->SetRangeUser(0, 800);
  Eff_PHO_RemoveJets->GetYaxis()->SetTitle("Efficiency");
  Eff_PHO_RemoveJets->GetXaxis()->SetTitle("Offline Photon pT [GeV]");
  Eff_PHO_RemoveJets->SetMinimum(0.0);
  Eff_PHO_RemoveJets->SetTitle("Efficiency: HLT_Photon30_R9Id90_Met25 w.r.t HLT_Met100");

  Eff_PHO_RemoveJets->Draw("AP E");


  TF1 *fit_PHO = new TF1("fit_PHO","0.5*[2]*(1+TMath::Erf((x-[0])/[1]))",0.0, 800.0);
  fit_PHO->SetParameter(0, 33);
  fit_PHO->SetParameter(1, 11);
  fit_PHO->SetParameter(2, 0.9);
  fit_PHO->SetParLimits(2,0.0,1.0);
  Eff_PHO->Fit("fit_PHO", "", "", 0,800.0);

  fit_PHO->SetLineColor(kRed);
  fit_PHO->SetLineWidth(2);

  TF1 *fit_PHO_RemoveJets = new TF1("fit_PHO_RemoveJets","0.5*[2]*(1+TMath::Erf((x-[0])/[1]))",0.0, 800.0);
  fit_PHO_RemoveJets->SetParameter(0, 20);
  fit_PHO_RemoveJets->SetParameter(1, 40);
  fit_PHO_RemoveJets->SetParameter(2, 0.8);
  fit_PHO_RemoveJets->SetParLimits(2,0.0,1.0);
  Eff_PHO_RemoveJets->Fit("fit_PHO_RemoveJets", "", "", 0,800.0);

  fit_PHO_RemoveJets->SetLineColor(kRed);
  fit_PHO_RemoveJets->SetLineWidth(2);

  std::string histfilename=(outfile+".root").c_str();
  TFile *tFile=new TFile(histfilename.c_str(), "RECREATE");
  Eff_PHO->Write("Eff_PHO");
  fit_PHO->Write("fit_PHO");
  Eff_PHO_RemoveJets->Write("Eff_PHO_RemoveJets");
  fit_PHO_RemoveJets->Write("fit_PHO_RemoveJets");
  tFile->Close();
}
