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

void Optimization2D_h_HT_Memu_LowCutValue(std::string outfile){

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
  double signal_stop_SF_120 = (1.66*0.1943*1000*7.4*100*0.33*0.33)/99986;
  double signal_stop_SF_140 = (1.66*0.104*1000*7.4*0.33*0.33*2)/99986;
  double signal_SF = (0.0004318*(7.4)*1000*1000*1.66)/99990;

  TFile* file1 = TFile::Open("MC_OS_OptimizationStudies/Output_LowPtSUSY_Tree_ZGToLLG_8TeV_CaloMET_All.root");
  TH2F *h1 = (TH2F*)file1->Get("h_HT_Memu_LowCutValue");
  //h1->SetLineWidth(3);
  h1->SetLineColor(kBlue);
  h1->SetLineWidth(3);
  h1->SetStats(kFALSE);
  h1->SetTitle("HT versus M_{e#mu}");
  h1->GetXaxis()->SetLabelSize(0.03);
  h1->GetYaxis()->SetLabelSize(0.03);
  h1->GetXaxis()->SetTitle("HT < Cut Value [GeV]");
  h1->GetYaxis()->SetTitle("M_{e#mu} < Cut Value [GeV]");
  h1->GetXaxis()->SetRangeUser(0, 500);
  h1->GetYaxis()->SetRangeUser(0, 600);
  h1->GetYaxis()->SetTitleOffset(1.2);
  //h1->Draw("BOX");

  TFile* file2 = TFile::Open("MC_OS_OptimizationStudies/Output_LowPtSUSY_Tree_TTGJets_8TeV_madgraph_CaloMET_All.root");
  TH2F *h2 = (TH2F*)file2->Get("h_HT_Memu_LowCutValue");
  h2->SetLineWidth(3);
  h2->SetLineColor(kRed);
  //h2->Draw("BOX SAME");

  TFile* file3 = new TFile("MC_OS_OptimizationStudies/Output_LowPtSUSY_Tree_WGstarToLNu2Mu_TuneZ2star_CaloMET_All.root");
  TH2F *h3 = (TH2F*)file3->Get("h_HT_Memu_LowCutValue");

  TFile* file4 = new TFile("MC_OS_OptimizationStudies/Output_LowPtSUSY_Tree_ZZJetsTo4L_TuneZ2star_CaloMET_All.root");
  TH2F *h4 = (TH2F*)file4->Get("h_HT_Memu_LowCutValue");

  TFile* file5 = new TFile("MC_OS_OptimizationStudies/Output_LowPtSUSY_Tree_T_tW_CaloMET_All.root");
  TH2F *h5 = (TH2F*)file5->Get("h_HT_Memu_LowCutValue");

  TFile* file6 = new TFile("MC_OS_OptimizationStudies/Output_LowPtSUSY_Tree_Tbar_tW_CaloMET_All.root");  
  TH2F *h6 = (TH2F*)file6->Get("h_HT_Memu_LowCutValue");

  TFile* file7 = new TFile("MC_OS_OptimizationStudies/Output_LowPtSUSY_Tree_WZJetsTo3LNu_CaloMET_All.root");
  TH2F *h7 = (TH2F*)file7->Get("h_HT_Memu_LowCutValue");

  TFile* file8 = new TFile("MC_OS_OptimizationStudies/Output_LowPtSUSY_Tree_WWGJets_8TeV_madgraph_CaloMET_All.root");
  TH2F *h8 = (TH2F*)file8->Get("h_HT_Memu_LowCutValue");

  TFile* file9 = new TFile("Data_SS_OptimizationStudies/Output_LowPtSUSY_Tree_SinglePhotonParked_Run2012D_22Jan2013_All.root");
  TH2F *h9 = (TH2F*)file9->Get("h_HT_Memu_LowCutValue");

  TFile* file10 = TFile::Open("MC_OS_Signal_OptimizationStudies/Output_LowPtSUSY_Tree_Stop140_Chargino115_Neutralino70_Ntuple_Pta20.root");
  TH2F *h10 = (TH2F*)file10->Get("h_HT_Memu_LowCutValue");
  h10->SetLineWidth(3);
  h10->SetLineColor(kBlack);
  //h10->Draw("BOX SAME");

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
  leg1->AddEntry(h3,"Stop 140","l");
  //leg1->Draw();


  TH2F *h_AllBackgrounds=(TH2F*)h1->Clone("h_AllBackgrounds");
  h_AllBackgrounds->Reset();
  h_AllBackgrounds->Add(h1);
  h_AllBackgrounds->Add(h2);
  h_AllBackgrounds->Add(h3);
  h_AllBackgrounds->Add(h4);
  h_AllBackgrounds->Add(h5);
  h_AllBackgrounds->Add(h6);
  h_AllBackgrounds->Add(h7);
  h_AllBackgrounds->Add(h8);

  TH2F *h_Significance=(TH2F*)h_AllBackgrounds->Clone("h_Significance");
  h_Significance->Reset();

  TAxis *xaxis = h3->GetXaxis();
  TAxis *yaxis = h3->GetYaxis();
  for (int j=1; j<=yaxis->GetNbins();j++) {
      for (int i=1; i<=xaxis->GetNbins();i++) {
         double sumOfBackgrounds = h1->GetBinContent(i, j) + h2->GetBinContent(i, j) + h3->GetBinContent(i, j) + h4->GetBinContent(i, j) + h5->GetBinContent(i, j) + h6->GetBinContent(i, j) + h7->GetBinContent(i, j) + h8->GetBinContent(i, j);
         if(sumOfBackgrounds > 0.0) h_Significance->SetBinContent(i, j, (h10->GetBinContent(i, j)*signal_stop_SF_140/
                                                    (TMath::Sqrt(h1->GetBinContent(i, j)*ZGammaLL_SF + 
                                                                 h2->GetBinContent(i, j)*TTG_SF + 
                                                                 h3->GetBinContent(i, j)*WGSMu_SF + 
                                                                 h4->GetBinContent(i, j)*ZZ_SF + 
                                                                 h5->GetBinContent(i, j)*T_tW_SF + 
                                                                 h6->GetBinContent(i, j)*Tbar_tW_SF +
                                                                 h7->GetBinContent(i, j)*WZ_SF +
                                                                 h8->GetBinContent(i, j)*WWG_SF +
                                                                 h9->GetBinContent(i, j)  +
                                                                 h1->GetBinContent(i, j)*ZGammaLL_SF*115.118/2332.79
                                                 ))));
         else h_Significance->SetBinContent(i, j, 0.0);
         cout << "h1->GetBinContent(6, 10) = " << h1->GetBinContent(6, 10) << endl;
         cout << "h1->FindBinX(145) = " << h1->GetXaxis()->FindBin(145) << endl;
         cout << "h1->FindBinY(90) = " << h1->GetYaxis()->FindBin(90) << endl;
      }
  }

  h_Significance->Draw("COLZ TEXT");
  c1.SaveAs("h_Significance_h_HT_Memu_LowCutValue.pdf");
  c1.SaveAs("h_Significance_h_HT_Memu_LowCutValue.png");
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
  h10->Write("h_STOP140");
  h_Significance->Write();
  tFile->Close();

}


void Optimization2D_h_HT_Memu_HighCutValue(std::string outfile){

  TCanvas c1("c1","Significance Plot",10,10,900,800); 

  gStyle->SetErrorX(0);
  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat("4.1f");
  gStyle->SetLineWidth(3);

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
  double signal_stop_SF_120 = (1.66*0.1943*1000*7.4*100*0.33*0.33)/99986;
  double signal_stop_SF_140 = (1.66*0.104*1000*7.4*0.33*0.33*2)/99986;
  double signal_SF = (0.0004318*(7.4)*1000*1000*1.66)/99990;

  TFile* file1 = TFile::Open("MC_OS_OptimizationStudies/Output_LowPtSUSY_Tree_ZGToLLG_8TeV_CaloMET_All.root");
  TH2F *h1 = (TH2F*)file1->Get("h_HT_Memu_HighCutValue");
  //h1->SetLineWidth(3);
  h1->SetLineColor(kBlue);
  h1->SetLineWidth(3);
  h1->SetStats(kFALSE);
  h1->SetTitle("HT versus M_{e#mu}");
  h1->GetXaxis()->SetLabelSize(0.03);
  h1->GetYaxis()->SetLabelSize(0.03);
  h1->GetXaxis()->SetTitle("HT < Cut Value [GeV]");
  h1->GetYaxis()->SetTitle("M_{e#mu} > Cut Value [GeV]");
  h1->GetXaxis()->SetRangeUser(0, 500);
  h1->GetYaxis()->SetRangeUser(0, 600);
  h1->GetYaxis()->SetTitleOffset(1.2);
  //h1->Draw("BOX");

  TFile* file2 = TFile::Open("MC_OS_OptimizationStudies/Output_LowPtSUSY_Tree_TTGJets_8TeV_madgraph_CaloMET_All.root");
  TH2F *h2 = (TH2F*)file2->Get("h_HT_Memu_HighCutValue");
  h2->SetLineWidth(3);
  h2->SetLineColor(kRed);
  //h2->Draw("BOX SAME");

  TFile* file3 = new TFile("MC_OS_OptimizationStudies/Output_LowPtSUSY_Tree_WGstarToLNu2Mu_TuneZ2star_CaloMET_All.root");
  TH2F *h3 = (TH2F*)file3->Get("h_HT_Memu_HighCutValue");

  TFile* file4 = new TFile("MC_OS_OptimizationStudies/Output_LowPtSUSY_Tree_ZZJetsTo4L_TuneZ2star_CaloMET_All.root");
  TH2F *h4 = (TH2F*)file4->Get("h_HT_Memu_HighCutValue");

  TFile* file5 = new TFile("MC_OS_OptimizationStudies/Output_LowPtSUSY_Tree_T_tW_CaloMET_All.root");
  TH2F *h5 = (TH2F*)file5->Get("h_HT_Memu_HighCutValue");

  TFile* file6 = new TFile("MC_OS_OptimizationStudies/Output_LowPtSUSY_Tree_Tbar_tW_CaloMET_All.root");  
  TH2F *h6 = (TH2F*)file6->Get("h_HT_Memu_HighCutValue");

  TFile* file7 = new TFile("MC_OS_OptimizationStudies/Output_LowPtSUSY_Tree_WZJetsTo3LNu_CaloMET_All.root");
  TH2F *h7 = (TH2F*)file7->Get("h_HT_Memu_HighCutValue");

  TFile* file8 = new TFile("MC_OS_OptimizationStudies/Output_LowPtSUSY_Tree_WWGJets_8TeV_madgraph_CaloMET_All.root");
  TH2F *h8 = (TH2F*)file8->Get("h_HT_Memu_HighCutValue");

  TFile* file9 = new TFile("Data_SS_OptimizationStudies/Output_LowPtSUSY_Tree_SinglePhotonParked_Run2012D_22Jan2013_All.root");
  TH2F *h9 = (TH2F*)file9->Get("h_HT_Memu_HighCutValue");

  TFile* file10 = TFile::Open("MC_OS_Signal_OptimizationStudies/Output_LowPtSUSY_Tree_Stop140_Chargino115_Neutralino70_Ntuple_Pta20.root");
  TH2F *h10 = (TH2F*)file10->Get("h_HT_Memu_HighCutValue");
  h10->SetLineWidth(3);
  h10->SetLineColor(kBlack);
  //h10->Draw("BOX SAME");

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
  leg1->AddEntry(h3,"Stop 140","l");
  //leg1->Draw();


  TH2F *h_AllBackgrounds=(TH2F*)h1->Clone("h_AllBackgrounds");
  h_AllBackgrounds->Reset();
  h_AllBackgrounds->Add(h1);
  h_AllBackgrounds->Add(h2);
  h_AllBackgrounds->Add(h3);
  h_AllBackgrounds->Add(h4);
  h_AllBackgrounds->Add(h5);
  h_AllBackgrounds->Add(h6);
  h_AllBackgrounds->Add(h7);
  h_AllBackgrounds->Add(h8);

  TH2F *h_Significance=(TH2F*)h_AllBackgrounds->Clone("h_Significance");
  h_Significance->Reset();

  TAxis *xaxis = h3->GetXaxis();
  TAxis *yaxis = h3->GetYaxis();
  for (int j=1; j<=yaxis->GetNbins();j++) {
      for (int i=1; i<=xaxis->GetNbins();i++) {
         double sumOfBackgrounds = h1->GetBinContent(i, j) + h2->GetBinContent(i, j) + h3->GetBinContent(i, j) + h4->GetBinContent(i, j) + h5->GetBinContent(i, j) + h6->GetBinContent(i, j) + h7->GetBinContent(i, j) + h8->GetBinContent(i, j);
         if(sumOfBackgrounds > 0.0) h_Significance->SetBinContent(i, j, (h10->GetBinContent(i, j)*signal_stop_SF_140/
                                                    (TMath::Sqrt(h1->GetBinContent(i, j)*ZGammaLL_SF + 
                                                                 h2->GetBinContent(i, j)*TTG_SF + 
                                                                 h3->GetBinContent(i, j)*WGSMu_SF + 
                                                                 h4->GetBinContent(i, j)*ZZ_SF + 
                                                                 h5->GetBinContent(i, j)*T_tW_SF + 
                                                                 h6->GetBinContent(i, j)*Tbar_tW_SF +
                                                                 h7->GetBinContent(i, j)*WZ_SF +
                                                                 h8->GetBinContent(i, j)*WWG_SF +
                                                                 h9->GetBinContent(i, j) +
                                                                 h1->GetBinContent(i, j)*ZGammaLL_SF*(115.118/2332.79) 
                                                   ))));
         else h_Significance->SetBinContent(i, j, 0.0);
      }
  }

  h_Significance->Draw("COLZ TEXT");

  c1.SaveAs("h_Significance_h_HT_Memu_HighCutValue.pdf");
  c1.SaveAs("h_Significance_h_HT_Memu_HighCutValue.png");
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
  h10->Write("h_STOP140");
  h_Significance->Write();
  tFile->Close();

}

void Optimization2D_h_HT_MemuGamma_LowCutValue(std::string outfile){

  TCanvas c1("c1","Significance Plot",10,10,900,800); 

  gStyle->SetErrorX(0);
  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat("4.1f");
  gStyle->SetLineWidth(3);

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
  double signal_stop_SF_120 = (1.66*0.1943*1000*7.4*100*0.33*0.33)/99986;
  double signal_stop_SF_140 = (1.66*0.104*1000*7.4*0.33*0.33*2)/99986;
  double signal_SF = (0.0004318*(7.4)*1000*1000*1.66)/99990;

  TFile* file1 = TFile::Open("MC_OS_OptimizationStudies/Output_LowPtSUSY_Tree_ZGToLLG_8TeV_CaloMET_All.root");
  TH2F *h1 = (TH2F*)file1->Get("h_HT_MemuGamma_LowCutValue");
  //h1->SetLineWidth(3);
  h1->SetLineColor(kBlue);
  h1->SetLineWidth(3);
  h1->SetStats(kFALSE);
  h1->SetTitle("HT versus M_{e#mu#gamma}");
  h1->GetXaxis()->SetLabelSize(0.03);
  h1->GetYaxis()->SetLabelSize(0.03);
  h1->GetXaxis()->SetTitle("HT < Cut Value [GeV]");
  h1->GetYaxis()->SetTitle("M_{e#mu#gamma} < Cut Value [GeV]");
  h1->GetXaxis()->SetRangeUser(0, 500);
  h1->GetYaxis()->SetRangeUser(0, 600);
  h1->GetYaxis()->SetTitleOffset(1.2);
  //h1->Draw("BOX");

  TFile* file2 = TFile::Open("MC_OS_OptimizationStudies/Output_LowPtSUSY_Tree_TTGJets_8TeV_madgraph_CaloMET_All.root");
  TH2F *h2 = (TH2F*)file2->Get("h_HT_MemuGamma_LowCutValue");
  h2->SetLineWidth(3);
  h2->SetLineColor(kRed);
  //h2->Draw("BOX SAME");

  TFile* file3 = new TFile("MC_OS_OptimizationStudies/Output_LowPtSUSY_Tree_WGstarToLNu2Mu_TuneZ2star_CaloMET_All.root");
  TH2F *h3 = (TH2F*)file3->Get("h_HT_MemuGamma_LowCutValue");

  TFile* file4 = new TFile("MC_OS_OptimizationStudies/Output_LowPtSUSY_Tree_ZZJetsTo4L_TuneZ2star_CaloMET_All.root");
  TH2F *h4 = (TH2F*)file4->Get("h_HT_MemuGamma_LowCutValue");

  TFile* file5 = new TFile("MC_OS_OptimizationStudies/Output_LowPtSUSY_Tree_T_tW_CaloMET_All.root");
  TH2F *h5 = (TH2F*)file5->Get("h_HT_MemuGamma_LowCutValue");

  TFile* file6 = new TFile("MC_OS_OptimizationStudies/Output_LowPtSUSY_Tree_Tbar_tW_CaloMET_All.root");  
  TH2F *h6 = (TH2F*)file6->Get("h_HT_MemuGamma_LowCutValue");

  TFile* file7 = new TFile("MC_OS_OptimizationStudies/Output_LowPtSUSY_Tree_WZJetsTo3LNu_CaloMET_All.root");
  TH2F *h7 = (TH2F*)file7->Get("h_HT_MemuGamma_LowCutValue");

  TFile* file8 = new TFile("MC_OS_OptimizationStudies/Output_LowPtSUSY_Tree_WWGJets_8TeV_madgraph_CaloMET_All.root");
  TH2F *h8 = (TH2F*)file8->Get("h_HT_MemuGamma_LowCutValue");

  TFile* file9 = new TFile("Data_SS_OptimizationStudies/Output_LowPtSUSY_Tree_SinglePhotonParked_Run2012D_22Jan2013_All.root");
  TH2F *h9 = (TH2F*)file9->Get("h_HT_MemuGamma_LowCutValue");

  TFile* file10 = TFile::Open("MC_OS_Signal_OptimizationStudies/Output_LowPtSUSY_Tree_Stop140_Chargino115_Neutralino70_Ntuple_Pta20.root");
  TH2F *h10 = (TH2F*)file10->Get("h_HT_MemuGamma_LowCutValue");
  h10->SetLineWidth(3);
  h10->SetLineColor(kBlack);
  //h10->Draw("BOX SAME");

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
  leg1->AddEntry(h3,"Stop 140","l");
  //leg1->Draw();


  TH2F *h_AllBackgrounds=(TH2F*)h1->Clone("h_AllBackgrounds");
  h_AllBackgrounds->Reset();
  h_AllBackgrounds->Add(h1);
  h_AllBackgrounds->Add(h2);
  h_AllBackgrounds->Add(h3);
  h_AllBackgrounds->Add(h4);
  h_AllBackgrounds->Add(h5);
  h_AllBackgrounds->Add(h6);
  h_AllBackgrounds->Add(h7);
  h_AllBackgrounds->Add(h8);

  TH2F *h_Significance=(TH2F*)h_AllBackgrounds->Clone("h_Significance");
  h_Significance->Reset();

  TAxis *xaxis = h3->GetXaxis();
  TAxis *yaxis = h3->GetYaxis();
  for (int j=1; j<=yaxis->GetNbins();j++) {
      for (int i=1; i<=xaxis->GetNbins();i++) {
         double sumOfBackgrounds = h1->GetBinContent(i, j) + h2->GetBinContent(i, j) + h3->GetBinContent(i, j) + h4->GetBinContent(i, j) + h5->GetBinContent(i, j) + h6->GetBinContent(i, j) + h7->GetBinContent(i, j) + h8->GetBinContent(i, j);
         if(sumOfBackgrounds > 0.0) h_Significance->SetBinContent(i, j, (h10->GetBinContent(i, j)*signal_stop_SF_140/
                                                    (TMath::Sqrt(h1->GetBinContent(i, j)*ZGammaLL_SF + 
                                                                 h2->GetBinContent(i, j)*TTG_SF + 
                                                                 h3->GetBinContent(i, j)*WGSMu_SF + 
                                                                 h4->GetBinContent(i, j)*ZZ_SF + 
                                                                 h5->GetBinContent(i, j)*T_tW_SF + 
                                                                 h6->GetBinContent(i, j)*Tbar_tW_SF +
                                                                 h7->GetBinContent(i, j)*WZ_SF +
                                                                 h8->GetBinContent(i, j)*WWG_SF +
                                                                 h9->GetBinContent(i, j) +
                                                                 h1->GetBinContent(i, j)*ZGammaLL_SF*(115.118/2332.79) 
                                                   ))));
         else h_Significance->SetBinContent(i, j, 0.0);
      }
  }

  h_Significance->Draw("COLZ TEXT");

  c1.SaveAs("h_Significance_h_HT_MemuGamma_LowCutValue.pdf");
  c1.SaveAs("h_Significance_h_HT_MemuGamma_LowCutValue.png");
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
  h10->Write("h_STOP140");
  h_Significance->Write();
  tFile->Close();

}


void Optimization2D_h_HTb_MemuGamma_LowCutValue(std::string outfile){

  TCanvas c1("c1","Significance Plot",10,10,900,800); 

  gStyle->SetErrorX(0);
  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat("4.1f");
  gStyle->SetLineWidth(3);

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
  double signal_stop_SF_120 = (1.66*0.1943*1000*7.4*100*0.33*0.33)/99986;
  double signal_stop_SF_140 = (1.66*0.104*1000*7.4*0.33*0.33*2)/99986;
  double signal_SF = (0.0004318*(7.4)*1000*1000*1.66)/99990;

  TFile* file1 = TFile::Open("MC_OS_OptimizationStudies/Output_LowPtSUSY_Tree_ZGToLLG_8TeV_CaloMET_All.root");
  TH2F *h1 = (TH2F*)file1->Get("h_HTb_MemuGamma_LowCutValue");
  //h1->SetLineWidth(3);
  h1->SetLineColor(kBlue);
  h1->SetLineWidth(3);
  h1->SetStats(kFALSE);
  h1->SetTitle("HTb versus M_{e#mu#gamma}");
  h1->GetXaxis()->SetLabelSize(0.03);
  h1->GetYaxis()->SetLabelSize(0.03);
  h1->GetXaxis()->SetTitle("HTb < Cut Value [GeV]");
  h1->GetYaxis()->SetTitle("M_{e#mu#gamma} < Cut Value [GeV]");
  h1->GetXaxis()->SetRangeUser(0, 500);
  h1->GetYaxis()->SetRangeUser(0, 600);
  h1->GetYaxis()->SetTitleOffset(1.2);
  //h1->Draw("BOX");

  TFile* file2 = TFile::Open("MC_OS_OptimizationStudies/Output_LowPtSUSY_Tree_TTGJets_8TeV_madgraph_CaloMET_All.root");
  TH2F *h2 = (TH2F*)file2->Get("h_HTb_MemuGamma_LowCutValue");
  h2->SetLineWidth(3);
  h2->SetLineColor(kRed);
  //h2->Draw("BOX SAME");

  TFile* file3 = new TFile("MC_OS_OptimizationStudies/Output_LowPtSUSY_Tree_WGstarToLNu2Mu_TuneZ2star_CaloMET_All.root");
  TH2F *h3 = (TH2F*)file3->Get("h_HTb_MemuGamma_LowCutValue");

  TFile* file4 = new TFile("MC_OS_OptimizationStudies/Output_LowPtSUSY_Tree_ZZJetsTo4L_TuneZ2star_CaloMET_All.root");
  TH2F *h4 = (TH2F*)file4->Get("h_HTb_MemuGamma_LowCutValue");

  TFile* file5 = new TFile("MC_OS_OptimizationStudies/Output_LowPtSUSY_Tree_T_tW_CaloMET_All.root");
  TH2F *h5 = (TH2F*)file5->Get("h_HTb_MemuGamma_LowCutValue");

  TFile* file6 = new TFile("MC_OS_OptimizationStudies/Output_LowPtSUSY_Tree_Tbar_tW_CaloMET_All.root");  
  TH2F *h6 = (TH2F*)file6->Get("h_HTb_MemuGamma_LowCutValue");

  TFile* file7 = new TFile("MC_OS_OptimizationStudies/Output_LowPtSUSY_Tree_WZJetsTo3LNu_CaloMET_All.root");
  TH2F *h7 = (TH2F*)file7->Get("h_HTb_MemuGamma_LowCutValue");

  TFile* file8 = new TFile("MC_OS_OptimizationStudies/Output_LowPtSUSY_Tree_WWGJets_8TeV_madgraph_CaloMET_All.root");
  TH2F *h8 = (TH2F*)file8->Get("h_HTb_MemuGamma_LowCutValue");

  TFile* file9 = new TFile("Data_SS_OptimizationStudies/Output_LowPtSUSY_Tree_SinglePhotonParked_Run2012D_22Jan2013_All.root");
  TH2F *h9 = (TH2F*)file9->Get("h_HTb_MemuGamma_LowCutValue");

  TFile* file10 = TFile::Open("MC_OS_Signal_OptimizationStudies/Output_LowPtSUSY_Tree_Stop140_Chargino115_Neutralino70_Ntuple_Pta20.root");
  TH2F *h10 = (TH2F*)file10->Get("h_HTb_MemuGamma_LowCutValue");
  h10->SetLineWidth(3);
  h10->SetLineColor(kBlack);
  //h10->Draw("BOX SAME");

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
  leg1->AddEntry(h3,"Stop 140","l");
  //leg1->Draw();


  TH2F *h_AllBackgrounds=(TH2F*)h1->Clone("h_AllBackgrounds");
  h_AllBackgrounds->Reset();
  h_AllBackgrounds->Add(h1);
  h_AllBackgrounds->Add(h2);
  h_AllBackgrounds->Add(h3);
  h_AllBackgrounds->Add(h4);
  h_AllBackgrounds->Add(h5);
  h_AllBackgrounds->Add(h6);
  h_AllBackgrounds->Add(h7);
  h_AllBackgrounds->Add(h8);

  TH2F *h_Significance=(TH2F*)h_AllBackgrounds->Clone("h_Significance");
  h_Significance->Reset();

  TAxis *xaxis = h3->GetXaxis();
  TAxis *yaxis = h3->GetYaxis();
  for (int j=1; j<=yaxis->GetNbins();j++) {
      for (int i=1; i<=xaxis->GetNbins();i++) {
         double sumOfBackgrounds = h1->GetBinContent(i, j) + h2->GetBinContent(i, j) + h3->GetBinContent(i, j) + h4->GetBinContent(i, j) + h5->GetBinContent(i, j) + h6->GetBinContent(i, j) + h7->GetBinContent(i, j) + h8->GetBinContent(i, j);
         if(sumOfBackgrounds > 0.0) h_Significance->SetBinContent(i, j, (h10->GetBinContent(i, j)*signal_stop_SF_140/
                                                    (TMath::Sqrt(h1->GetBinContent(i, j)*ZGammaLL_SF + 
                                                                 h2->GetBinContent(i, j)*TTG_SF + 
                                                                 h3->GetBinContent(i, j)*WGSMu_SF + 
                                                                 h4->GetBinContent(i, j)*ZZ_SF + 
                                                                 h5->GetBinContent(i, j)*T_tW_SF + 
                                                                 h6->GetBinContent(i, j)*Tbar_tW_SF +
                                                                 h7->GetBinContent(i, j)*WZ_SF +
                                                                 h8->GetBinContent(i, j)*WWG_SF +
                                                                 h9->GetBinContent(i, j) +
                                                                 h1->GetBinContent(i, j)*ZGammaLL_SF*(115.118/2332.79) 
                                                   ))));
         else h_Significance->SetBinContent(i, j, 0.0);
      }
  }

  h_Significance->Draw("COLZ TEXT");

  c1.SaveAs("h_Significance_h_HTb_MemuGamma_LowCutValue.pdf");
  c1.SaveAs("h_Significance_h_HTb_MemuGamma_LowCutValue.png");
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
  h10->Write("h_STOP140");
  h_Significance->Write();
  tFile->Close();

}


void Optimization2D_h_HTb_Memu_HighCutValue(std::string outfile){

  TCanvas c1("c1","Significance Plot",10,10,900,800); 

  gStyle->SetErrorX(0);
  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat("4.1f");
  gStyle->SetLineWidth(3);

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
  double signal_stop_SF_120 = (1.66*0.1943*1000*7.4*100*0.33*0.33)/99986;
  double signal_stop_SF_140 = (1.66*0.104*1000*7.4*0.33*0.33*2)/99986;
  double signal_SF = (0.0004318*(7.4)*1000*1000*1.66)/99990;

  TFile* file1 = TFile::Open("MC_OS_OptimizationStudies/Output_LowPtSUSY_Tree_ZGToLLG_8TeV_CaloMET_All.root");
  TH2F *h1 = (TH2F*)file1->Get("h_HTb_Memu_HighCutValue");
  //h1->SetLineWidth(3);
  h1->SetLineColor(kBlue);
  h1->SetLineWidth(3);
  h1->SetStats(kFALSE);
  h1->SetTitle("HTb versus M_{e#mu}");
  h1->GetXaxis()->SetLabelSize(0.03);
  h1->GetYaxis()->SetLabelSize(0.03);
  h1->GetXaxis()->SetTitle("HTb < Cut Value [GeV]");
  h1->GetYaxis()->SetTitle("M_{e#mu} > Cut Value [GeV]");
  h1->GetXaxis()->SetRangeUser(0, 500);
  h1->GetYaxis()->SetRangeUser(0, 600);
  h1->GetYaxis()->SetTitleOffset(1.2);
  //h1->Draw("BOX");

  TFile* file2 = TFile::Open("MC_OS_OptimizationStudies/Output_LowPtSUSY_Tree_TTGJets_8TeV_madgraph_CaloMET_All.root");
  TH2F *h2 = (TH2F*)file2->Get("h_HTb_Memu_HighCutValue");
  h2->SetLineWidth(3);
  h2->SetLineColor(kRed);
  //h2->Draw("BOX SAME");

  TFile* file3 = new TFile("MC_OS_OptimizationStudies/Output_LowPtSUSY_Tree_WGstarToLNu2Mu_TuneZ2star_CaloMET_All.root");
  TH2F *h3 = (TH2F*)file3->Get("h_HTb_Memu_HighCutValue");

  TFile* file4 = new TFile("MC_OS_OptimizationStudies/Output_LowPtSUSY_Tree_ZZJetsTo4L_TuneZ2star_CaloMET_All.root");
  TH2F *h4 = (TH2F*)file4->Get("h_HTb_Memu_HighCutValue");

  TFile* file5 = new TFile("MC_OS_OptimizationStudies/Output_LowPtSUSY_Tree_T_tW_CaloMET_All.root");
  TH2F *h5 = (TH2F*)file5->Get("h_HTb_Memu_HighCutValue");

  TFile* file6 = new TFile("MC_OS_OptimizationStudies/Output_LowPtSUSY_Tree_Tbar_tW_CaloMET_All.root");  
  TH2F *h6 = (TH2F*)file6->Get("h_HTb_Memu_HighCutValue");

  TFile* file7 = new TFile("MC_OS_OptimizationStudies/Output_LowPtSUSY_Tree_WZJetsTo3LNu_CaloMET_All.root");
  TH2F *h7 = (TH2F*)file7->Get("h_HTb_Memu_HighCutValue");

  TFile* file8 = new TFile("MC_OS_OptimizationStudies/Output_LowPtSUSY_Tree_WWGJets_8TeV_madgraph_CaloMET_All.root");
  TH2F *h8 = (TH2F*)file8->Get("h_HTb_Memu_HighCutValue");

  TFile* file9 = new TFile("Data_SS_OptimizationStudies/Output_LowPtSUSY_Tree_SinglePhotonParked_Run2012D_22Jan2013_All.root");
  TH2F *h9 = (TH2F*)file9->Get("h_HTb_Memu_HighCutValue");

  TFile* file10 = TFile::Open("MC_OS_Signal_OptimizationStudies/Output_LowPtSUSY_Tree_Stop140_Chargino115_Neutralino70_Ntuple_Pta20.root");
  TH2F *h10 = (TH2F*)file10->Get("h_HTb_Memu_HighCutValue");
  h10->SetLineWidth(3);
  h10->SetLineColor(kBlack);
  //h10->Draw("BOX SAME");

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
  leg1->AddEntry(h3,"Stop 140","l");
  //leg1->Draw();


  TH2F *h_AllBackgrounds=(TH2F*)h1->Clone("h_AllBackgrounds");
  h_AllBackgrounds->Reset();
  h_AllBackgrounds->Add(h1);
  h_AllBackgrounds->Add(h2);
  h_AllBackgrounds->Add(h3);
  h_AllBackgrounds->Add(h4);
  h_AllBackgrounds->Add(h5);
  h_AllBackgrounds->Add(h6);
  h_AllBackgrounds->Add(h7);
  h_AllBackgrounds->Add(h8);

  TH2F *h_Significance=(TH2F*)h_AllBackgrounds->Clone("h_Significance");
  h_Significance->Reset();

  TAxis *xaxis = h3->GetXaxis();
  TAxis *yaxis = h3->GetYaxis();
  for (int j=1; j<=yaxis->GetNbins();j++) {
      for (int i=1; i<=xaxis->GetNbins();i++) {
         double sumOfBackgrounds = h1->GetBinContent(i, j) + h2->GetBinContent(i, j) + h3->GetBinContent(i, j) + h4->GetBinContent(i, j) + h5->GetBinContent(i, j) + h6->GetBinContent(i, j) + h7->GetBinContent(i, j) + h8->GetBinContent(i, j);
         if(sumOfBackgrounds > 0.0) h_Significance->SetBinContent(i, j, (h10->GetBinContent(i, j)*signal_stop_SF_140/
                                                    (TMath::Sqrt(h1->GetBinContent(i, j)*ZGammaLL_SF + 
                                                                 h2->GetBinContent(i, j)*TTG_SF + 
                                                                 h3->GetBinContent(i, j)*WGSMu_SF + 
                                                                 h4->GetBinContent(i, j)*ZZ_SF + 
                                                                 h5->GetBinContent(i, j)*T_tW_SF + 
                                                                 h6->GetBinContent(i, j)*Tbar_tW_SF +
                                                                 h7->GetBinContent(i, j)*WZ_SF +
                                                                 h8->GetBinContent(i, j)*WWG_SF +
                                                                 h9->GetBinContent(i, j) +
                                                                 h1->GetBinContent(i, j)*ZGammaLL_SF*(115.118/2332.79) 
                                                   ))));
         else h_Significance->SetBinContent(i, j, 0.0);
      }
  }

  h_Significance->Draw("COLZ TEXT");

  c1.SaveAs("h_Significance_h_HTb_Memu_HighCutValue.pdf");
  c1.SaveAs("h_Significance_h_HTb_Memu_HighCutValue.png");
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
  h10->Write("h_STOP140");
  h_Significance->Write();
  tFile->Close();

}

void Optimization2D_h_HTb_Memu_LowCutValue(std::string outfile){

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
  double signal_stop_SF_120 = (1.66*0.1943*1000*7.4*100*0.33*0.33)/99986;
  double signal_stop_SF_140 = (1.66*0.104*1000*7.4*0.33*0.33*2)/99986;
  double signal_SF = (0.0004318*(7.4)*1000*1000*1.66)/99990;

  TFile* file1 = TFile::Open("MC_OS_OptimizationStudies/Output_LowPtSUSY_Tree_ZGToLLG_8TeV_CaloMET_All.root");
  TH2F *h1 = (TH2F*)file1->Get("h_HTb_Memu_LowCutValue");
  //h1->SetLineWidth(3);
  h1->SetLineColor(kBlue);
  h1->SetLineWidth(3);
  h1->SetStats(kFALSE);
  h1->SetTitle("HTb versus M_{e#mu}");
  h1->GetXaxis()->SetLabelSize(0.03);
  h1->GetYaxis()->SetLabelSize(0.03);
  h1->GetXaxis()->SetTitle("HTb < Cut Value [GeV]");
  h1->GetYaxis()->SetTitle("M_{e#mu} < Cut Value [GeV]");
  h1->GetXaxis()->SetRangeUser(0, 500);
  h1->GetYaxis()->SetRangeUser(0, 600);
  h1->GetYaxis()->SetTitleOffset(1.2);
  //h1->Draw("BOX");

  TFile* file2 = TFile::Open("MC_OS_OptimizationStudies/Output_LowPtSUSY_Tree_TTGJets_8TeV_madgraph_CaloMET_All.root");
  TH2F *h2 = (TH2F*)file2->Get("h_HTb_Memu_LowCutValue");
  h2->SetLineWidth(3);
  h2->SetLineColor(kRed);
  //h2->Draw("BOX SAME");

  TFile* file3 = new TFile("MC_OS_OptimizationStudies/Output_LowPtSUSY_Tree_WGstarToLNu2Mu_TuneZ2star_CaloMET_All.root");
  TH2F *h3 = (TH2F*)file3->Get("h_HTb_Memu_LowCutValue");

  TFile* file4 = new TFile("MC_OS_OptimizationStudies/Output_LowPtSUSY_Tree_ZZJetsTo4L_TuneZ2star_CaloMET_All.root");
  TH2F *h4 = (TH2F*)file4->Get("h_HTb_Memu_LowCutValue");

  TFile* file5 = new TFile("MC_OS_OptimizationStudies/Output_LowPtSUSY_Tree_T_tW_CaloMET_All.root");
  TH2F *h5 = (TH2F*)file5->Get("h_HTb_Memu_LowCutValue");

  TFile* file6 = new TFile("MC_OS_OptimizationStudies/Output_LowPtSUSY_Tree_Tbar_tW_CaloMET_All.root");  
  TH2F *h6 = (TH2F*)file6->Get("h_HTb_Memu_LowCutValue");

  TFile* file7 = new TFile("MC_OS_OptimizationStudies/Output_LowPtSUSY_Tree_WZJetsTo3LNu_CaloMET_All.root");
  TH2F *h7 = (TH2F*)file7->Get("h_HTb_Memu_LowCutValue");

  TFile* file8 = new TFile("MC_OS_OptimizationStudies/Output_LowPtSUSY_Tree_WWGJets_8TeV_madgraph_CaloMET_All.root");
  TH2F *h8 = (TH2F*)file8->Get("h_HTb_Memu_LowCutValue");

  TFile* file9 = new TFile("Data_SS_OptimizationStudies/Output_LowPtSUSY_Tree_SinglePhotonParked_Run2012D_22Jan2013_All.root");
  TH2F *h9 = (TH2F*)file9->Get("h_HTb_Memu_LowCutValue");

  TFile* file10 = TFile::Open("MC_OS_Signal_OptimizationStudies/Output_LowPtSUSY_Tree_Stop140_Chargino115_Neutralino70_Ntuple_Pta20.root");
  TH2F *h10 = (TH2F*)file10->Get("h_HTb_Memu_LowCutValue");
  h10->SetLineWidth(3);
  h10->SetLineColor(kBlack);
  //h10->Draw("BOX SAME");

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
  leg1->AddEntry(h3,"Stop 140","l");
  //leg1->Draw();


  TH2F *h_AllBackgrounds=(TH2F*)h1->Clone("h_AllBackgrounds");
  h_AllBackgrounds->Reset();
  h_AllBackgrounds->Add(h1);
  h_AllBackgrounds->Add(h2);
  h_AllBackgrounds->Add(h3);
  h_AllBackgrounds->Add(h4);
  h_AllBackgrounds->Add(h5);
  h_AllBackgrounds->Add(h6);
  h_AllBackgrounds->Add(h7);
  h_AllBackgrounds->Add(h8);

  TH2F *h_Significance=(TH2F*)h_AllBackgrounds->Clone("h_Significance");
  h_Significance->Reset();

  TAxis *xaxis = h3->GetXaxis();
  TAxis *yaxis = h3->GetYaxis();
  for (int j=1; j<=yaxis->GetNbins();j++) {
      for (int i=1; i<=xaxis->GetNbins();i++) {
         double sumOfBackgrounds = h1->GetBinContent(i, j) + h2->GetBinContent(i, j) + h3->GetBinContent(i, j) + h4->GetBinContent(i, j) + h5->GetBinContent(i, j) + h6->GetBinContent(i, j) + h7->GetBinContent(i, j) + h8->GetBinContent(i, j);
         if(sumOfBackgrounds > 0.0) h_Significance->SetBinContent(i, j, (h10->GetBinContent(i, j)*signal_stop_SF_140/
                                                    (TMath::Sqrt(h1->GetBinContent(i, j)*ZGammaLL_SF + 
                                                                 h2->GetBinContent(i, j)*TTG_SF + 
                                                                 h3->GetBinContent(i, j)*WGSMu_SF + 
                                                                 h4->GetBinContent(i, j)*ZZ_SF + 
                                                                 h5->GetBinContent(i, j)*T_tW_SF + 
                                                                 h6->GetBinContent(i, j)*Tbar_tW_SF +
                                                                 h7->GetBinContent(i, j)*WZ_SF +
                                                                 h8->GetBinContent(i, j)*WWG_SF +
                                                                 h9->GetBinContent(i, j)  +
                                                                 h1->GetBinContent(i, j)*ZGammaLL_SF*115.118/2332.79
                                                   ))));
         else h_Significance->SetBinContent(i, j, 0.0);
      }
  }

  h_Significance->Draw("COLZ TEXT");
  c1.SaveAs("h_Significance_h_HTb_Memu_LowCutValue.pdf");
  c1.SaveAs("h_Significance_h_HTb_Memu_LowCutValue.png");
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
  h10->Write("h_STOP140");
  h_Significance->Write();
  tFile->Close();

}
