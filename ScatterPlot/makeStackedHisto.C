#include "TH1.h"
#include "TTree.h"
#include "TKey.h"
#include "Riostream.h"
#include <TCanvas.h>
#include <TLatex.h>
#include "TGraphErrors.h"
#include "TLegend.h"
#include <TPad.h>
#include <sstream>
#include "TVectorD.h"
#include "TGraph.h"

using std::string;
using std::cout;
using std::endl;
using std::istringstream;


void makeStackedHisto(){

gROOT->SetStyle("Plain");

TCanvas c1("c1","Stacked Histogram",10,10,600,400);
c1.SetLogy();
double ZGammaLL_SF = (1.2*132.6*7.3*1000)/6583032;
double ZGammaNU_SF = (1.2*123.9*7.3*1000)/3169096;
double WGamma_SF = (1.2*461.6*7.3*1000)/4802339;
//double ZGammaLL_SF = (132.6*0.10*7.3*1000*1.2)/6583032; //Z->LL BR
//double ZGammaNU_SF = (123.9*0.20*7.3*1000*1.2)/3169096; //Z->nunu BR
//double WGamma_SF = (461.6*0.33*7.3*1000)/4802339;
double TTG_SF = (1.2*2.166*7.3*1000)/71328;

THStack hs("hs","Stacked Plot of HT");

TFile* file1 = new TFile("Output_ZGamma_Inclusive_Jun12.root");
TH1F *h1 = (TH1F*)file1->Get("h_HT");
h1->GetXaxis()->SetRangeUser(1, 500);
h1->Scale(ZGammaNU_SF);
h1->SetFillColor(kRed);
hs.Add(h1);

TFile* file2 = TFile::Open("Output_ZGToLLG_Jun11.root");
TH1F *h2 = (TH1F*)file2->Get("h_HT");
h2->GetXaxis()->SetRangeUser(1, 500);
h2->Scale(ZGammaLL_SF);
h2->SetFillColor(kBlue);
hs.Add(h2);

TFile* file3 = TFile::Open("Output_WGamma_Jun11.root");
TH1F *h3 = (TH1F*)file3->Get("h_HT");
h3->GetXaxis()->SetRangeUser(1, 500);
h3->Scale(WGamma_SF);
h3->SetFillColor(kGreen);
hs.Add(h3);

TFile* file4 = TFile::Open("Output_LowPtSUSY_Tree_TTG.root");
TH1F *h4 = (TH1F*)file4->Get("h_HT");
h4->GetXaxis()->SetRangeUser(1, 500);
h4->Scale(TTG_SF);
h4->SetFillColor(kCyan);
hs.Add(h4);

hs.Draw("HIST");
//hs.SetMaximum(1000000);
hs.GetXaxis()->SetRangeUser(1, 500);

TLegend *leg1 = new TLegend(0.55,0.55,0.90,0.90,NULL,"brNDC");
leg1->SetBorderSize(0);
leg1->SetTextSize(0.03);
leg1->SetLineColor(1);
leg1->SetLineStyle(0);
leg1->SetLineWidth(1);
leg1->SetFillColor(10);
leg1->SetFillStyle(0);
leg1->AddEntry(h1,"Z#gamma #rightarrow #nu#nu#gamma","f");
leg1->AddEntry(h2,"Z#gamma #rightarrow ll#gamma","f");
leg1->AddEntry(h3,"W#gamma #rightarrow l#nu#gamma","f");
leg1->AddEntry(h4,"t#bar{t}#gamma","f");
leg1->Draw();


/*TFile* file5 = TFile::Open("Output_PULooseJetID.root");
TH1F *h5 = (TH1F*)file5->Get("h_HT");
h5->GetXaxis()->SetRangeUser(1, 500);
h5->SetMarkerSize(1.3);
h5->Draw("SAME e");
*/
c1.SaveAs("h_HT_Stacked.pdf");
c1.SaveAs("h_HT_Stacked.png");
}
