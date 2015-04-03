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

void Optimization1D_h_HT_LowCutValue_emu(std::string outfile){

  gStyle->SetErrorX(0);
  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat("4.1f");
  gStyle->SetLineWidth(3);

  TCanvas c1("c1","Significance Plot",10,10,900,800);

  double ZGammaLL_SF = (1.27*1.2*132.6*7.4*1000)/(6583032.0);
  double ZGammaNU_SF = (1.2*123.9*7.4*1000)/3169096;
  double WGamma_SF = (1.2*461.6*7.4*1000)/4802339;
  double TTG_SF = (1.2*2.166*7.4*1000)/71328;
  double DY_SF = (1.2*3532.8*7.4*1000)/25283595;
  double WWG_SF = (0.528*1.2*7.4*1000)/214538;
  double WGSMu_SF = (1.914*1.2*7.4*1000)/299866;
  double ZZ_SF = (0.1769*1.2*7.4*1000)/4622329;
  double Tbar_tW_SF = (11.1*1.2*7.4*1000)/452907;
  double T_tW_SF = (11.1*1.2*7.4*1000)/445199;
  double WZ_SF = (1.0575*1.2*(7.4*1000))/1911473;
  double signal_stop_SF_120 = (1.66*0.1943*1000*7.4*0.33*0.33)/99986;
  double signal_stop_SF_100 = (1.66*0.3983*1000*7.4*0.33*0.33*2)/99994;
  double signal_stop_SF_140 = (1.66*0.104*1000*7.4*0.33*0.33*2)/99996;
  double signal_stop_SF_180 = (1.66*0.03592*1000*7.4*0.33*0.33*2)/99987;
  double signal_Chargino1 = (1.66*0.0004318*100*1000*7.4*0.33*0.33*2)/99993;
  double signal_SF = (0.0004318*(7.4)*1000*1000*1.66)/99990;

  TFile* file1 = TFile::Open("MC_OS_OptimizationStudies/Output_LowPtSUSY_Tree_ZGToLLG_8TeV_CaloMET_All.root");
  TH1F *h1 = (TH1F*)file1->Get("h_HT_LowCutValue_emu");
  //h1->SetLineWidth(3);
  h1->SetLineColor(kBlue);
  h1->SetLineWidth(3);
  h1->SetStats(kFALSE);
  h1->SetTitle("Significance");
  h1->GetXaxis()->SetLabelSize(0.03);
  h1->GetYaxis()->SetLabelSize(0.03);
  h1->GetXaxis()->SetTitle("HT < Cut Value [GeV]");
  h1->GetYaxis()->SetTitle("s/#sqrt{b}");
  h1->GetXaxis()->SetRangeUser(1, 500);
  h1->GetYaxis()->SetRangeUser(1, 600);
  h1->GetYaxis()->SetTitleOffset(1.2);
  //h1->Draw("BOX");

  TFile* file2 = TFile::Open("MC_OS_OptimizationStudies/Output_LowPtSUSY_Tree_TTGJets_8TeV_madgraph_CaloMET_All.root");
  TH1F *h2 = (TH1F*)file2->Get("h_HT_LowCutValue_emu");
  h2->SetLineWidth(3);
  h2->SetLineColor(kRed);
  //h2->Draw("BOX SAME");

  TFile* file3 = new TFile("MC_OS_OptimizationStudies/Output_LowPtSUSY_Tree_WGstarToLNu2Mu_TuneZ2star_CaloMET_All.root");
  TH1F *h3 = (TH1F*)file3->Get("h_HT_LowCutValue_emu");

  TFile* file4 = new TFile("MC_OS_OptimizationStudies/Output_LowPtSUSY_Tree_ZZJetsTo4L_TuneZ2star_CaloMET_All.root");
  TH1F *h4 = (TH1F*)file4->Get("h_HT_LowCutValue_emu");

  TFile* file5 = new TFile("MC_OS_OptimizationStudies/Output_LowPtSUSY_Tree_T_tW_CaloMET_All.root");
  TH1F *h5 = (TH1F*)file5->Get("h_HT_LowCutValue_emu");

  TFile* file6 = new TFile("MC_OS_OptimizationStudies/Output_LowPtSUSY_Tree_Tbar_tW_CaloMET_All.root");  
  TH1F *h6 = (TH1F*)file6->Get("h_HT_LowCutValue_emu");

  TFile* file7 = new TFile("MC_OS_OptimizationStudies/Output_LowPtSUSY_Tree_WZJetsTo3LNu_CaloMET_All.root");
  TH1F *h7 = (TH1F*)file7->Get("h_HT_LowCutValue_emu");

  TFile* file8 = new TFile("MC_OS_OptimizationStudies/Output_LowPtSUSY_Tree_WWGJets_8TeV_madgraph_CaloMET_All.root");
  TH1F *h8 = (TH1F*)file8->Get("h_HT_LowCutValue_emu");

  TFile* file9 = new TFile("Data_SS_OptimizationStudies/Output_LowPtSUSY_Tree_SinglePhotonParked_Run2012D_22Jan2013_All.root");
  TH1F *h9 = (TH1F*)file9->Get("h_HT_LowCutValue_emu");

  //TFile* file10 = TFile::Open("MC_OS_Signal_OptimizationStudies/Output_LowPtSUSY_Tree_Stop80_Chargino55_Neutralino10_Ntuple_Pta20.root");
  TFile* file10 = TFile::Open("MC_OS_Signal_OptimizationStudies/Output_LowPtSUSY_Tree_Stop100_Chargino75_Neutralino30_Ntuple_Pta20.root");
  TH1F *h10 = (TH1F*)file10->Get("h_HT_LowCutValue_emu");
  h10->SetLineWidth(3);
  h10->SetLineColor(kBlack);
  //h10->Draw("BOX SAME");

  TFile* file11 = TFile::Open("MC_OS_Signal_OptimizationStudies/Output_LowPtSUSY_Tree_Stop140_Chargino115_Neutralino70_Ntuple_Pta20.root");
  TH1F *h11 = (TH1F*)file11->Get("h_HT_LowCutValue_emu");
  h11->SetLineWidth(3);
  h11->SetLineColor(kBlack);

  TFile* file12 = TFile::Open("MC_OS_Signal_OptimizationStudies/Output_LowPtSUSY_Tree_Stop180_Chargino155_Neutralino110_Ntuple_Pta20.root");
  TH1F *h12 = (TH1F*)file12->Get("h_HT_LowCutValue_emu");
  h12->SetLineWidth(3);
  h12->SetLineColor(kBlack);

  TFile* file13 = TFile::Open("MC_OS_Signal_OptimizationStudies/Output_LowPtSUSY_Tree_Chargino1.root");
  TH1F *h13 = (TH1F*)file13->Get("h_HT_LowCutValue_emu");
  h13->SetLineWidth(3);
  h13->SetLineColor(kBlack);

  TLegend *leg1 = new TLegend(0.7787356,0.756871,0.8793103,0.8773784,NULL,"brNDC");
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.03);
  leg1->SetLineColor(1);
  leg1->SetLineStyle(0);
  leg1->SetLineWidth(1);
  leg1->SetFillColor(10);
  leg1->SetFillStyle(0);
  leg1->AddEntry(h1,"ZGamma","l");
  leg1->AddEntry(h2,"TTG","l");
  leg1->AddEntry(h3,"Stop 100","l");
  //leg1->Draw();


  TH1F *h_AllBackgrounds=(TH1F*)h1->Clone("h_AllBackgrounds");
  h_AllBackgrounds->Reset();
  h_AllBackgrounds->Add(h1);
  h_AllBackgrounds->Add(h2);
  h_AllBackgrounds->Add(h3);
  h_AllBackgrounds->Add(h4);
  h_AllBackgrounds->Add(h5);
  h_AllBackgrounds->Add(h6);
  h_AllBackgrounds->Add(h7);
  h_AllBackgrounds->Add(h8);

  TH1F *h_Significance_Stop100=(TH1F*)h_AllBackgrounds->Clone("h_Significance_Stop100");
  h_Significance_Stop100->Reset();

  TH1F *h_Significance_Stop140=(TH1F*)h_AllBackgrounds->Clone("h_Significance_Stop140");
  h_Significance_Stop140->Reset();

  TH1F *h_Significance_Stop180=(TH1F*)h_AllBackgrounds->Clone("h_Significance_Stop180");
  h_Significance_Stop180->Reset();

  TH1F *h_Significance_Chargino1=(TH1F*)h_AllBackgrounds->Clone("h_Significance_Chargino1");
  h_Significance_Chargino1->Reset();

  TAxis *xaxis = h3->GetXaxis();
  for (int i=1; i<=xaxis->GetNbins();i++) {
        double sumOfBackgrounds = h1->GetBinContent(i) + h2->GetBinContent(i) + h3->GetBinContent(i) + h4->GetBinContent(i) + h5->GetBinContent(i) + h6->GetBinContent(i) + h7->GetBinContent(i) + h8->GetBinContent(i);
         if(sumOfBackgrounds > 0.0) h_Significance_Stop100->SetBinContent(i, (h10->GetBinContent(i)*signal_stop_SF_100/
                                                    (TMath::Sqrt(h1->GetBinContent(i)*ZGammaLL_SF + 
                                                                 h2->GetBinContent(i)*TTG_SF + 
                                                                 h3->GetBinContent(i)*WGSMu_SF + 
                                                                 h4->GetBinContent(i)*ZZ_SF + 
                                                                 h5->GetBinContent(i)*T_tW_SF + 
                                                                 h6->GetBinContent(i)*Tbar_tW_SF +
                                                                 h7->GetBinContent(i)*WZ_SF +
                                                                 h8->GetBinContent(i)*WWG_SF +
                                                                 h9->GetBinContent(i)  +
                                                                 h1->GetBinContent(i)*ZGammaLL_SF*115.118/2332.79
                                                   ))));
         else h_Significance_Stop100->SetBinContent(i, 0.0);
      
  }


   for (int i=1; i<=xaxis->GetNbins();i++) {
        double sumOfBackgrounds = h1->GetBinContent(i) + h2->GetBinContent(i) + h3->GetBinContent(i) + h4->GetBinContent(i) + h5->GetBinContent(i) + h6->GetBinContent(i) + h7->GetBinContent(i) + h8->GetBinContent(i);
         if(sumOfBackgrounds > 0.0) h_Significance_Stop140->SetBinContent(i, (h11->GetBinContent(i)*signal_stop_SF_140/
                                                    (TMath::Sqrt(h1->GetBinContent(i)*ZGammaLL_SF +
                                                                 h2->GetBinContent(i)*TTG_SF +
                                                                 h3->GetBinContent(i)*WGSMu_SF +
                                                                 h4->GetBinContent(i)*ZZ_SF +
                                                                 h5->GetBinContent(i)*T_tW_SF +
                                                                 h6->GetBinContent(i)*Tbar_tW_SF +
                                                                 h7->GetBinContent(i)*WZ_SF +
                                                                 h8->GetBinContent(i)*WWG_SF +
                                                                 h9->GetBinContent(i)  +
                                                                 h1->GetBinContent(i)*ZGammaLL_SF*115.118/2332.79
                                                   ))));
         else h_Significance_Stop140->SetBinContent(i, 0.0);

  }  

   for (int i=1; i<=xaxis->GetNbins();i++) {
        double sumOfBackgrounds = h1->GetBinContent(i) + h2->GetBinContent(i) + h3->GetBinContent(i) + h4->GetBinContent(i) + h5->GetBinContent(i) + h6->GetBinContent(i) + h7->GetBinContent(i) + h8->GetBinContent(i);
         if(sumOfBackgrounds > 0.0) h_Significance_Stop180->SetBinContent(i, (h12->GetBinContent(i)*signal_stop_SF_180/
                                                    (TMath::Sqrt(h1->GetBinContent(i)*ZGammaLL_SF +
                                                                 h2->GetBinContent(i)*TTG_SF +
                                                                 h3->GetBinContent(i)*WGSMu_SF +
                                                                 h4->GetBinContent(i)*ZZ_SF +
                                                                 h5->GetBinContent(i)*T_tW_SF +
                                                                 h6->GetBinContent(i)*Tbar_tW_SF +
                                                                 h7->GetBinContent(i)*WZ_SF +
                                                                 h8->GetBinContent(i)*WWG_SF +
                                                                 h9->GetBinContent(i)  +
                                                                 h1->GetBinContent(i)*ZGammaLL_SF*115.118/2332.79
                                                   ))));
         else h_Significance_Stop180->SetBinContent(i, 0.0);

  }  

  for (int i=1; i<=xaxis->GetNbins();i++) {
        double sumOfBackgrounds = h1->GetBinContent(i) + h2->GetBinContent(i) + h3->GetBinContent(i) + h4->GetBinContent(i) + h5->GetBinContent(i) + h6->GetBinContent(i) + h7->GetBinContent(i) + h8->GetBinContent(i);
         if(sumOfBackgrounds > 0.0) h_Significance_Chargino1->SetBinContent(i, (h13->GetBinContent(i)*signal_Chargino1/
                                                    (TMath::Sqrt(h1->GetBinContent(i)*ZGammaLL_SF +
                                                                 h2->GetBinContent(i)*TTG_SF +
                                                                 h3->GetBinContent(i)*WGSMu_SF +
                                                                 h4->GetBinContent(i)*ZZ_SF +
                                                                 h5->GetBinContent(i)*T_tW_SF +
                                                                 h6->GetBinContent(i)*Tbar_tW_SF +
                                                                 h7->GetBinContent(i)*WZ_SF +
                                                                 h8->GetBinContent(i)*WWG_SF +
                                                                 h9->GetBinContent(i)  +
                                                                 h1->GetBinContent(i)*ZGammaLL_SF*115.118/2332.79
                                                   ))));
         else h_Significance_Chargino1->SetBinContent(i, 0.0);
         //cout << "h_Significance_Chargino1->GetBinContent(i) = " << h_Significance_Chargino1->GetBinContent(i) << endl;
         //cout << "h13->GetBinContent(i) = " << h13->GetBinContent(i) << endl;

  }  


  h_Significance_Stop100->Draw("");
  h_Significance_Stop140->SetLineColor(kBlack);
  h_Significance_Stop140->Draw("SAME");
  h_Significance_Stop180->SetLineColor(kRed);
  h_Significance_Stop180->Draw("SAME");
  h_Significance_Chargino1->SetLineColor(kGreen);
  h_Significance_Chargino1->Draw("SAME");

  TLegend *leg1 = new TLegend(0.52,0.30,0.90,0.55,NULL,"brNDC");
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.03);
  leg1->SetLineColor(1);
  leg1->SetLineStyle(0);
  leg1->SetLineWidth(1);
  leg1->SetFillColor(10);
  leg1->SetFillStyle(0);
  leg1->AddEntry(h_Significance_Stop100,"Stop 100 GeV","l");
  leg1->AddEntry(h_Significance_Stop140,"Stop 140 GeV","l");
  leg1->AddEntry(h_Significance_Stop180,"Stop 180 GeV","l");
  leg1->AddEntry(h_Significance_Chargino1,"Chargino 100 GeV X 100","l");
  leg1->Draw();

  c1.SaveAs("h_Significance_h_HT_LowCutValue_emu.pdf");
  c1.SaveAs("h_Significance_h_HT_LowCutValue_emu.png");

  std::string histfilename=(outfile+".root").c_str();
  TFile *tFile=new TFile(histfilename.c_str(), "RECREATE");
  h1->Write("h_ZGToLLG");
  h2->Write("h_TTG");
  h3->Write("h_WGstarToLNu2Mu");
  h4->Write("h_ZZJetsTo4L");
  h5->Write("h_T_tW");
  h6->Write("h_Tbar_tW");
  h7->Write("h_WZJetsTo3LNu");
  h8->Write("h_WWGJets");
  h9->Write("h_SS");
  h10->Write("h_STOP100");
  h11->Write("h_STOP140");
  h12->Write("h_STOP180");
  h13->Write("h_Chargino1");
  h_Significance_Stop100->Write();
  tFile->Close();


}

void Optimization1D_h_HTb_LowCutValue_emu(std::string outfile){

  gStyle->SetErrorX(0);
  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat("4.1f");
  gStyle->SetLineWidth(3);

  TCanvas c1("c1","Significance Plot",10,10,900,800);

  double ZGammaLL_SF = (1.27*1.2*132.6*7.4*1000)/(6583032.0);
  double ZGammaNU_SF = (1.2*123.9*7.4*1000)/3169096;
  double WGamma_SF = (1.2*461.6*7.4*1000)/4802339;
  double TTG_SF = (1.2*2.166*7.4*1000)/71328;
  double DY_SF = (1.2*3532.8*7.4*1000)/25283595;
  double WWG_SF = (0.528*1.2*7.4*1000)/214538;
  double WGSMu_SF = (1.914*1.2*7.4*1000)/299866;
  double ZZ_SF = (0.1769*1.2*7.4*1000)/4622329;
  double Tbar_tW_SF = (11.1*1.2*7.4*1000)/452907;
  double T_tW_SF = (11.1*1.2*7.4*1000)/445199;
  double WZ_SF = (1.0575*1.2*(7.4*1000))/1911473;

  double signal_stop_SF_120 = (1.66*0.1943*1000*7.4*0.33*0.33)/99986;
  double signal_stop_SF_100 = (1.66*0.3983*1000*7.4*0.33*0.33*2)/99994;
  double signal_stop_SF_140 = (1.66*0.104*1000*7.4*0.33*0.33*2)/99996;
  double signal_stop_SF_180 = (1.66*0.03592*1000*7.4*0.33*0.33*2)/99987;
  double signal_Chargino1 = (1.66*0.0004318*100*1000*7.4*0.33*0.33*2)/99993;

  TFile* file1 = TFile::Open("MC_OS_OptimizationStudies/Output_LowPtSUSY_Tree_ZGToLLG_8TeV_CaloMET_All.root");
  TH1F *h1 = (TH1F*)file1->Get("h_HTb_LowCutValue_emu");
  //h1->SetLineWidth(3);
  h1->SetLineColor(kBlue);
  h1->SetLineWidth(3);
  h1->SetStats(kFALSE);
  h1->SetTitle("Significance");
  h1->GetXaxis()->SetLabelSize(0.03);
  h1->GetYaxis()->SetLabelSize(0.03);
  h1->GetXaxis()->SetTitle("HTb < Cut Value [GeV]");
  h1->GetYaxis()->SetTitle("s/#sqrt{b}");
  h1->GetXaxis()->SetRangeUser(1, 500);
  h1->GetYaxis()->SetRangeUser(1, 600);
  h1->GetYaxis()->SetTitleOffset(1.2);
  //h1->Draw("BOX");

  TFile* file2 = TFile::Open("MC_OS_OptimizationStudies/Output_LowPtSUSY_Tree_TTGJets_8TeV_madgraph_CaloMET_All.root");
  TH1F *h2 = (TH1F*)file2->Get("h_HTb_LowCutValue_emu");
  h2->SetLineWidth(3);
  h2->SetLineColor(kRed);
  //h2->Draw("BOX SAME");

  TFile* file3 = new TFile("MC_OS_OptimizationStudies/Output_LowPtSUSY_Tree_WGstarToLNu2Mu_TuneZ2star_CaloMET_All.root");
  TH1F *h3 = (TH1F*)file3->Get("h_HTb_LowCutValue_emu");

  TFile* file4 = new TFile("MC_OS_OptimizationStudies/Output_LowPtSUSY_Tree_ZZJetsTo4L_TuneZ2star_CaloMET_All.root");
  TH1F *h4 = (TH1F*)file4->Get("h_HTb_LowCutValue_emu");

  TFile* file5 = new TFile("MC_OS_OptimizationStudies/Output_LowPtSUSY_Tree_T_tW_CaloMET_All.root");
  TH1F *h5 = (TH1F*)file5->Get("h_HTb_LowCutValue_emu");

  TFile* file6 = new TFile("MC_OS_OptimizationStudies/Output_LowPtSUSY_Tree_Tbar_tW_CaloMET_All.root");  
  TH1F *h6 = (TH1F*)file6->Get("h_HTb_LowCutValue_emu");

  TFile* file7 = new TFile("MC_OS_OptimizationStudies/Output_LowPtSUSY_Tree_WZJetsTo3LNu_CaloMET_All.root");
  TH1F *h7 = (TH1F*)file7->Get("h_HTb_LowCutValue_emu");

  TFile* file8 = new TFile("MC_OS_OptimizationStudies/Output_LowPtSUSY_Tree_WWGJets_8TeV_madgraph_CaloMET_All.root");
  TH1F *h8 = (TH1F*)file8->Get("h_HTb_LowCutValue_emu");

  TFile* file9 = new TFile("Data_SS_OptimizationStudies/Output_LowPtSUSY_Tree_SinglePhotonParked_Run2012D_22Jan2013_All.root");
  TH1F *h9 = (TH1F*)file9->Get("h_HTb_LowCutValue_emu");

  //TFile* file10 = TFile::Open("MC_OS_Signal_OptimizationStudies/Output_LowPtSUSY_Tree_Stop80_Chargino55_Neutralino10_Ntuple_Pta20.root");
  TFile* file10 = TFile::Open("MC_OS_Signal_OptimizationStudies/Output_LowPtSUSY_Tree_Stop100_Chargino75_Neutralino30_Ntuple_Pta20.root");
  TH1F *h10 = (TH1F*)file10->Get("h_HTb_LowCutValue_emu");
  h10->SetLineWidth(3);
  h10->SetLineColor(kBlack);
  //h10->Draw("BOX SAME");

  TFile* file11 = TFile::Open("MC_OS_Signal_OptimizationStudies/Output_LowPtSUSY_Tree_Stop140_Chargino115_Neutralino70_Ntuple_Pta20.root");
  TH1F *h11 = (TH1F*)file11->Get("h_HTb_LowCutValue_emu");
  h11->SetLineWidth(3);
  h11->SetLineColor(kBlack);

  TFile* file12 = TFile::Open("MC_OS_Signal_OptimizationStudies/Output_LowPtSUSY_Tree_Stop180_Chargino155_Neutralino110_Ntuple_Pta20.root");
  TH1F *h12 = (TH1F*)file12->Get("h_HTb_LowCutValue_emu");
  h12->SetLineWidth(3);
  h12->SetLineColor(kBlack);

  TFile* file13 = TFile::Open("MC_OS_Signal_OptimizationStudies/Output_LowPtSUSY_Tree_Chargino1.root");
  TH1F *h13 = (TH1F*)file13->Get("h_HTb_LowCutValue_emu");
  h13->SetLineWidth(3);
  h13->SetLineColor(kBlack);

  TLegend *leg1 = new TLegend(0.7787356,0.756871,0.8793103,0.8773784,NULL,"brNDC");
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.03);
  leg1->SetLineColor(1);
  leg1->SetLineStyle(0);
  leg1->SetLineWidth(1);
  leg1->SetFillColor(10);
  leg1->SetFillStyle(0);
  leg1->AddEntry(h1,"ZGamma","l");
  leg1->AddEntry(h2,"TTG","l");
  leg1->AddEntry(h3,"Stop 100","l");
  //leg1->Draw();


  TH1F *h_AllBackgrounds=(TH1F*)h1->Clone("h_AllBackgrounds");
  h_AllBackgrounds->Reset();
  h_AllBackgrounds->Add(h1);
  h_AllBackgrounds->Add(h2);
  h_AllBackgrounds->Add(h3);
  h_AllBackgrounds->Add(h4);
  h_AllBackgrounds->Add(h5);
  h_AllBackgrounds->Add(h6);
  h_AllBackgrounds->Add(h7);
  h_AllBackgrounds->Add(h8);

  TH1F *h_Significance_Stop100=(TH1F*)h_AllBackgrounds->Clone("h_Significance_Stop100");
  h_Significance_Stop100->Reset();

  TH1F *h_Significance_Stop140=(TH1F*)h_AllBackgrounds->Clone("h_Significance_Stop140");
  h_Significance_Stop140->Reset();

  TH1F *h_Significance_Stop180=(TH1F*)h_AllBackgrounds->Clone("h_Significance_Stop180");
  h_Significance_Stop180->Reset();

  TH1F *h_Significance_Chargino1=(TH1F*)h_AllBackgrounds->Clone("h_Significance_Chargino1");
  h_Significance_Chargino1->Reset();

  TAxis *xaxis = h3->GetXaxis();
  for (int i=1; i<=xaxis->GetNbins();i++) {
        double sumOfBackgrounds = h1->GetBinContent(i) + h2->GetBinContent(i) + h3->GetBinContent(i) + h4->GetBinContent(i) + h5->GetBinContent(i) + h6->GetBinContent(i) + h7->GetBinContent(i) + h8->GetBinContent(i);
         if(sumOfBackgrounds > 0.0) h_Significance_Stop100->SetBinContent(i, (h10->GetBinContent(i)*signal_stop_SF_100/
                                                    (TMath::Sqrt(h1->GetBinContent(i)*ZGammaLL_SF + 
                                                                 h2->GetBinContent(i)*TTG_SF + 
                                                                 h3->GetBinContent(i)*WGSMu_SF + 
                                                                 h4->GetBinContent(i)*ZZ_SF + 
                                                                 h5->GetBinContent(i)*T_tW_SF + 
                                                                 h6->GetBinContent(i)*Tbar_tW_SF +
                                                                 h7->GetBinContent(i)*WZ_SF +
                                                                 h8->GetBinContent(i)*WWG_SF +
                                                                 h9->GetBinContent(i)  +
                                                                 h1->GetBinContent(i)*ZGammaLL_SF*115.118/2332.79
                                                   ))));
         else h_Significance_Stop100->SetBinContent(i, 0.0);
      
  }


   for (int i=1; i<=xaxis->GetNbins();i++) {
        double sumOfBackgrounds = h1->GetBinContent(i) + h2->GetBinContent(i) + h3->GetBinContent(i) + h4->GetBinContent(i) + h5->GetBinContent(i) + h6->GetBinContent(i) + h7->GetBinContent(i) + h8->GetBinContent(i);
         if(sumOfBackgrounds > 0.0) h_Significance_Stop140->SetBinContent(i, (h11->GetBinContent(i)*signal_stop_SF_140/
                                                    (TMath::Sqrt(h1->GetBinContent(i)*ZGammaLL_SF +
                                                                 h2->GetBinContent(i)*TTG_SF +
                                                                 h3->GetBinContent(i)*WGSMu_SF +
                                                                 h4->GetBinContent(i)*ZZ_SF +
                                                                 h5->GetBinContent(i)*T_tW_SF +
                                                                 h6->GetBinContent(i)*Tbar_tW_SF +
                                                                 h7->GetBinContent(i)*WZ_SF +
                                                                 h8->GetBinContent(i)*WWG_SF +
                                                                 h9->GetBinContent(i)  +
                                                                 h1->GetBinContent(i)*ZGammaLL_SF*115.118/2332.79
                                                   ))));
         else h_Significance_Stop140->SetBinContent(i, 0.0);

  }  

   for (int i=1; i<=xaxis->GetNbins();i++) {
        double sumOfBackgrounds = h1->GetBinContent(i) + h2->GetBinContent(i) + h3->GetBinContent(i) + h4->GetBinContent(i) + h5->GetBinContent(i) + h6->GetBinContent(i) + h7->GetBinContent(i) + h8->GetBinContent(i);
         if(sumOfBackgrounds > 0.0) h_Significance_Stop180->SetBinContent(i, (h12->GetBinContent(i)*signal_stop_SF_180/
                                                    (TMath::Sqrt(h1->GetBinContent(i)*ZGammaLL_SF +
                                                                 h2->GetBinContent(i)*TTG_SF +
                                                                 h3->GetBinContent(i)*WGSMu_SF +
                                                                 h4->GetBinContent(i)*ZZ_SF +
                                                                 h5->GetBinContent(i)*T_tW_SF +
                                                                 h6->GetBinContent(i)*Tbar_tW_SF +
                                                                 h7->GetBinContent(i)*WZ_SF +
                                                                 h8->GetBinContent(i)*WWG_SF +
                                                                 h9->GetBinContent(i)  +
                                                                 h1->GetBinContent(i)*ZGammaLL_SF*115.118/2332.79
                                                   ))));
         else h_Significance_Stop180->SetBinContent(i, 0.0);

  }  

  for (int i=1; i<=xaxis->GetNbins();i++) {
        double sumOfBackgrounds = h1->GetBinContent(i) + h2->GetBinContent(i) + h3->GetBinContent(i) + h4->GetBinContent(i) + h5->GetBinContent(i) + h6->GetBinContent(i) + h7->GetBinContent(i) + h8->GetBinContent(i);
         if(sumOfBackgrounds > 0.0) h_Significance_Chargino1->SetBinContent(i, (h13->GetBinContent(i)*signal_Chargino1/
                                                    (TMath::Sqrt(h1->GetBinContent(i)*ZGammaLL_SF +
                                                                 h2->GetBinContent(i)*TTG_SF +
                                                                 h3->GetBinContent(i)*WGSMu_SF +
                                                                 h4->GetBinContent(i)*ZZ_SF +
                                                                 h5->GetBinContent(i)*T_tW_SF +
                                                                 h6->GetBinContent(i)*Tbar_tW_SF +
                                                                 h7->GetBinContent(i)*WZ_SF +
                                                                 h8->GetBinContent(i)*WWG_SF +
                                                                 h9->GetBinContent(i)  +
                                                                 h1->GetBinContent(i)*ZGammaLL_SF*115.118/2332.79
                                                   ))));
         else h_Significance_Chargino1->SetBinContent(i, 0.0);
         //cout << "h_Significance_Chargino1->GetBinContent(i) = " << h_Significance_Chargino1->GetBinContent(i) << endl;
         //cout << "h13->GetBinContent(i) = " << h13->GetBinContent(i) << endl;
  }  


  h_Significance_Stop100->Draw("");
  h_Significance_Stop140->SetLineColor(kBlack);
  h_Significance_Stop140->Draw("SAME");
  h_Significance_Stop180->SetLineColor(kRed);
  h_Significance_Stop180->Draw("SAME");
  h_Significance_Chargino1->SetLineColor(kGreen);
  h_Significance_Chargino1->Draw("SAME");

  TLegend *leg1 = new TLegend(0.52,0.30,0.90,0.55,NULL,"brNDC");
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.03);
  leg1->SetLineColor(1);
  leg1->SetLineStyle(0);
  leg1->SetLineWidth(1);
  leg1->SetFillColor(10);
  leg1->SetFillStyle(0);
  leg1->AddEntry(h_Significance_Stop100,"Stop 100 GeV","l");
  leg1->AddEntry(h_Significance_Stop140,"Stop 140 GeV","l");
  leg1->AddEntry(h_Significance_Stop180,"Stop 180 GeV","l");
  leg1->AddEntry(h_Significance_Chargino1,"Chargino 100 GeV X 100","l");
  leg1->Draw();

  c1.SaveAs("h_Significance_h_HTb_LowCutValue_emu.pdf");
  c1.SaveAs("h_Significance_h_HTb_LowCutValue_emu.png");

  std::string histfilename=(outfile+".root").c_str();
  TFile *tFile=new TFile(histfilename.c_str(), "RECREATE");
  h1->Write("h_ZGToLLG");
  h2->Write("h_TTG");
  h3->Write("h_WGstarToLNu2Mu");
  h4->Write("h_ZZJetsTo4L");
  h5->Write("h_T_tW");
  h6->Write("h_Tbar_tW");
  h7->Write("h_WZJetsTo3LNu");
  h8->Write("h_WWGJets");
  h9->Write("h_SS");
  h10->Write("h_STOP100");
  h11->Write("h_STOP140");
  h12->Write("h_STOP180");
  h13->Write("h_Chargino1");
  h_Significance_Stop100->Write();
  h_Significance_Stop140->Write();
  h_Significance_Stop180->Write();
  h_Significance_Chargino1->Write();
  tFile->Close();


}

