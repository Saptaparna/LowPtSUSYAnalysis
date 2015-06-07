#ifndef DrawOneHistoWithSyst_cxx
#define DrawOneHistoWithSyst_cxx

#include <TStyle.h>
#include <TString.h>
#include <TFile.h>
#include <TH1F.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TPaveText.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <TMath.h>
//include "LoadData_2012.C"
#include "tdrstyle.C"

// #include "tPrimeStandardHistos.C"

using namespace std;

enum Samples_t {Data,    ChargeMisID, FakeRate,
		TTbar,   ZJets,   DY1050,
		SingleTop,    DrellYan, WJets,
		WW,      WZ,      ZZ,
		T_tW,    Tbar_tW, T_t,     Tbar_t, T_s,   Tbar_s,
		WWm,     WWp,     WWW,
		TTW,     TTZ, 	  TTWW,
		ZZZNoGs, WWZNoGs, WZZNoGs, WWds,
		WGstarToLNu2E, WGstarToLNu2Mu, WGstarToLNu2Tau, ZGToLLG,
		Tprime400_bWbW, Tprime400_bWtZ, Tprime400_bWtH, Tprime400_tHtZ, Tprime400_tZtZ, Tprime400_tHtH,
		Tprime500_bWbW, Tprime500_bWtZ, Tprime500_bWtH, Tprime500_tHtZ, Tprime500_tZtZ, Tprime500_tHtH,
		Tprime600_bWbW, Tprime600_bWtZ, Tprime600_bWtH, Tprime600_tHtZ, Tprime600_tZtZ, Tprime600_tHtH,
		Tprime700_bWbW, Tprime700_bWtZ, Tprime700_bWtH, Tprime700_tHtZ, Tprime700_tZtZ, Tprime700_tHtH,
		Tprime800_bWbW, Tprime800_bWtZ, Tprime800_bWtH, Tprime800_tHtZ, Tprime800_tZtZ, Tprime800_tHtH,
		Tprime900_bWbW, Tprime900_bWtZ, Tprime900_bWtH, Tprime900_tHtZ, Tprime900_tZtZ, Tprime900_tHtH,
		Tprime1000_bWbW, Tprime1000_bWtZ, Tprime1000_bWtH, Tprime1000_tHtZ, Tprime1000_tZtZ, Tprime1000_tHtH,
		Tprime1100_bWbW, Tprime1100_bWtZ, Tprime1100_bWtH, Tprime1100_tHtZ, Tprime1100_tZtZ, Tprime1100_tHtH,
		Tprime1200_bWbW, Tprime1200_bWtZ, Tprime1200_bWtH, Tprime1200_tHtZ, Tprime1200_tZtZ, Tprime1200_tHtH,
		Tprime1300_bWbW, Tprime1300_bWtZ, Tprime1300_bWtH, Tprime1300_tHtZ, Tprime1300_tZtZ, Tprime1300_tHtH,
		Tprime1400_bWbW, Tprime1400_bWtZ, Tprime1400_bWtH, Tprime1400_tHtZ, Tprime1400_tZtZ, Tprime1400_tHtH,
		Tprime1500_bWbW, Tprime1500_bWtZ, Tprime1500_bWtH, Tprime1500_tHtZ, Tprime1500_tZtZ, Tprime1500_tHtH
};

const unsigned int totalSamples = 104;

Samples_t allSamples[totalSamples] = {Data,    ChargeMisID, FakeRate,
		TTbar,   ZJets,   DY1050,SingleTop,    DrellYan, WJets,
		WW,      WZ,      ZZ,
		T_tW,    Tbar_tW, T_t,     Tbar_t, T_s,   Tbar_s,
		WWm,     WWp,     WWW,
		TTW,     TTZ, 	  TTWW,
		ZZZNoGs, WWZNoGs, WZZNoGs, WWds,
		WGstarToLNu2E, WGstarToLNu2Mu, WGstarToLNu2Tau, ZGToLLG,
		Tprime400_bWbW, Tprime400_bWtZ, Tprime400_bWtH, Tprime400_tHtZ, Tprime400_tZtZ, Tprime400_tHtH,
		Tprime500_bWbW, Tprime500_bWtZ, Tprime500_bWtH, Tprime500_tHtZ, Tprime500_tZtZ, Tprime500_tHtH,
		Tprime600_bWbW, Tprime600_bWtZ, Tprime600_bWtH, Tprime600_tHtZ, Tprime600_tZtZ, Tprime600_tHtH,
		Tprime700_bWbW, Tprime700_bWtZ, Tprime700_bWtH, Tprime700_tHtZ, Tprime700_tZtZ, Tprime700_tHtH,
		Tprime800_bWbW, Tprime800_bWtZ, Tprime800_bWtH, Tprime800_tHtZ, Tprime800_tZtZ, Tprime800_tHtH,
		Tprime900_bWbW, Tprime900_bWtZ, Tprime900_bWtH, Tprime900_tHtZ, Tprime900_tZtZ, Tprime900_tHtH,
		Tprime1000_bWbW, Tprime1000_bWtZ, Tprime1000_bWtH, Tprime1000_tHtZ, Tprime1000_tZtZ, Tprime1000_tHtH,
		Tprime1100_bWbW, Tprime1100_bWtZ, Tprime1100_bWtH, Tprime1100_tHtZ, Tprime1100_tZtZ, Tprime1100_tHtH,
		Tprime1200_bWbW, Tprime1200_bWtZ, Tprime1200_bWtH, Tprime1200_tHtZ, Tprime1200_tZtZ, Tprime1200_tHtH,
		Tprime1300_bWbW, Tprime1300_bWtZ, Tprime1300_bWtH, Tprime1300_tHtZ, Tprime1300_tZtZ, Tprime1300_tHtH,
		Tprime1400_bWbW, Tprime1400_bWtZ, Tprime1400_bWtH, Tprime1400_tHtZ, Tprime1400_tZtZ, Tprime1400_tHtH,
		Tprime1500_bWbW, Tprime1500_bWtZ, Tprime1500_bWtH, Tprime1500_tHtZ, Tprime1500_tZtZ, Tprime1500_tHtH
};


TString allNames[totalSamples] = {"Data",    "ChargeMisID", "FakeRate",
			    "TTbar",   "ZJets",   "DY1050","SingleTop",    "DrellYan","WJets",
			    "WW",      "WZ",      "ZZ", 
			    "T_tW",    "Tbar_tW", "T_t",     "Tbar_t", "T_s", "Tbar_s",
			    "WWm",     "WWp",     "WWW",  
			    "TTW",     "TTZ", "TTWW",
			     "ZZZNoGs", "WWZNoGs", "WZZNoGs", "WWds",
			     "WGstarToLNu2E", "WGstarToLNu2Mu", "WGstarToLNu2Tau", "ZGToLLG",
			    "Tprime400_BWBW", "Tprime400_BWTZ", "Tprime400_BWTH", "Tprime400_THTZ", "Tprime400_TZTZ", "Tprime400_THTH",
			    "Tprime500_BWBW", "Tprime500_BWTZ", "Tprime500_BWTH", "Tprime500_THTZ", "Tprime500_TZTZ", "Tprime500_THTH",
			    "Tprime600_BWBW", "Tprime600_BWTZ", "Tprime600_BWTH", "Tprime600_THTZ", "Tprime600_TZTZ", "Tprime600_THTH",
			    "Tprime700_BWBW", "Tprime700_BWTZ", "Tprime700_BWTH", "Tprime700_THTZ", "Tprime700_TZTZ", "Tprime700_THTH",
			    "Tprime800_BWBW", "Tprime800_BWTZ", "Tprime800_BWTH", "Tprime800_THTZ", "Tprime800_TZTZ", "Tprime800_THTH",
			    "Tprime900_BWBW", "Tprime900_BWTZ", "Tprime900_BWTH", "Tprime900_THTZ", "Tprime900_TZTZ", "Tprime900_THTH",
			    "Tprime1000_BWBW", "Tprime1000_BWTZ", "Tprime1000_BWTH", "Tprime1000_THTZ", "Tprime1000_TZTZ", "Tprime1000_THTH",
			    "Tprime1100_BWBW", "Tprime1100_BWTZ", "Tprime1100_BWTH", "Tprime1100_THTZ", "Tprime1100_TZTZ", "Tprime1100_THTH",
			    "Tprime1200_BWBW", "Tprime1200_BWTZ", "Tprime1200_BWTH", "Tprime1200_THTZ", "Tprime1200_TZTZ", "Tprime1200_THTH",
			    "Tprime1300_BWBW", "Tprime1300_BWTZ", "Tprime1300_BWTH", "Tprime1300_THTZ", "Tprime1300_TZTZ", "Tprime1300_THTH",
			    "Tprime1400_BWBW", "Tprime1400_BWTZ", "Tprime1400_BWTH", "Tprime1400_THTZ", "Tprime1400_TZTZ", "Tprime1400_THTH",
			    "Tprime1500_BWBW", "Tprime1500_BWTZ", "Tprime1500_BWTH", "Tprime1500_THTZ", "Tprime1500_TZTZ", "Tprime1500_THTH"
};

TString legendNames[totalSamples] = {"Data",    "Charge MisID", "Non-prompt",
			    "t#bar{t}",   "ZJets",   "DY1050","Single top",    "Drell Yan","W+Jets",
			    "WW",      "WZ",      "ZZ",
			    "T_tW",    "Tbar_tW", "T_t",     "Tbar_t", "T_s", "Tbar_s",
			    "WWm",     "WWp",     "WWW",
			    "t#bar{t}W",     "t#bar{t}Z", 	  "t#bar{t}WW",
			     "ZZZNoGs", "WWZNoGs", "WZZNoGs", "WWds",
			     "WGstarToLNu2E", "WGstarToLNu2Mu", "WGstarToLNu2Tau", "ZGToLLG",
			    "t'400 bWbW", "t'400 bWtZ", "t'400 bWtH", "t'400 tHtZ", "t'400 tZtZ", "t'400 tHtH",
			    "t'500 bWbW", "t'500 bWtZ", "t'500 bWtH", "t'500 tHtZ", "t'500 tZtZ", "t'500 tHtH",
			    "t'600 bWbW", "t'600 bWtZ", "t'600 bWtH", "t'600 tHtZ", "t'600 tZtZ", "t'600 tHtH",
			    "t'700 bWbW", "t'700 bWtZ", "t'700 bWtH", "t'700 tHtZ", "t'700 tZtZ", "t'700 tHtH",
			    "t'800 bWbW", "t'800 bWtZ", "t'800 bWtH", "t'800 tHtZ", "t'800 tZtZ", "t'800 tHtH",
			    "t'900 bWbW", "t'900 bWtZ", "t'900 bWtH", "t'900 tHtZ", "t'900 tZtZ", "t'900 tHtH",
			    "t'1000 bWbW", "t'1000 bWtZ", "t'1000 bWtH", "t'1000 tHtZ", "t'1000 tZtZ", "t'1000 tHtH",
			    "t'1100 bWbW", "t'1100 bWtZ", "t'1100 bWtH", "t'1100 tHtZ", "t'1100 tZtZ", "t'1100 tHtH",
			    "t'1200 bWbW", "t'1200 bWtZ", "t'1200 bWtH", "t'1200 tHtZ", "t'1200 tZtZ", "t'1200 tHtH",
			    "t'1300 bWbW", "t'1300 bWtZ", "t'1300 bWtH", "t'1300 tHtZ", "t'1300 tZtZ", "t'1300 tHtH",
			    "t'1400 bWbW", "t'1400 bWtZ", "t'1400 bWtH", "t'1400 tHtZ", "t'1400 tZtZ", "t'1400 tHtH",
			    "t'1500 bWbW", "t'1500 bWtZ", "t'1500 bWtH", "t'1500 tHtZ", "t'1500 tZtZ", "t'1500 tHtH"
};

Color_t color[totalSamples] = {
	kBlack, kRed, kGray,				// Data,    ChargeMisID, FakeRate,
	kRed+1,   kBlue,kBlue,  kMagenta, kBlue,  kBlue+3, 	// TTbar,   ZJets,   DY1050,Single Top, DY,WJets,
	kGreen -3,kGreen -3,kGreen -3,  //WW,      WZ,      ZZ,
	//kGreen-9,      kViolet+1,      kAzure-2,
	kMagenta,    kMagenta, kMagenta,  kMagenta, kMagenta,   kMagenta, //Single Top
	kGreen-3,     kGreen-3,     kGreen+2,  	// WWm,     WWp,     WWW,
	kOrange,     kOrange,kOrange+2,		// TTW,     TTZ,     TTWW,
kCyan, kBlue, kOrange+4, kViolet+2,  kMagenta, 0,kMagenta,
	kViolet-4, kRed-4, kGreen-4, kGray+2,  kBlue-4, kMagenta-4, // Tprime 400 
	kViolet-4, kRed-4, kGreen-4, kGray+2,  kBlue-4, kMagenta-4, // Tprime 500 
	kViolet-4, kRed-4, kGreen-4, kGray+2,  kBlue-4, kMagenta-4, // Tprime 600 
	kViolet-4, kRed-4, kGreen-4, kGray+2,  kBlue-4, kMagenta-4, // Tprime 700 
	kViolet-4, kRed-4, kGreen-4, kGray+2,  kBlue-4, kMagenta-4, // Tprime 800 
	kViolet-4, kRed-4, kGreen-4, kGray+2,  kBlue-4, kMagenta-4, // Tprime 900 
	kViolet-4, kRed-4, kGreen-4, kGray+2,  kBlue-4, kMagenta-4, // Tprime 1000
	kViolet-4, kRed-4, kGreen-4, kGray+2,  kBlue-4, kMagenta-4, // Tprime 1100
	kViolet-4, kRed-4, kGreen-4, kGray+2,  kBlue-4, kMagenta-4, // Tprime 1200
	kViolet-4, kRed-4, kGreen-4, kGray+2,  kBlue-4, kMagenta-4, // Tprime 1300
	kViolet-4, kRed-4, kGreen-4, kGray+2,  kBlue-4, kMagenta-4, // Tprime 1400
	kViolet-4, kRed-4, kGreen-4, kGray+2,  kBlue-4, kMagenta-4  // Tprime 1500
//	kCyan, kBlue, kOrange+4, kViolet+2,  kMagenta, kMagenta+4, // Tprime 400
	};

vector <unsigned int> indices(totalSamples, -1);

//This is for OS:

/*
// const unsigned int NSAMPLES = 10;
// Samples_t samples[NSAMPLES] = {
// Data, TTbar, SingleTop, DrellYan,
// Tprime800_bWbW, Tprime800_bWtH, Tprime800_tZtZ, Tprime800_tHtH, Tprime800_bWtZ, Tprime800_tHtZ
// //ZJets, DY1050
// };
*/
const unsigned int NSAMPLES = 16;
 Samples_t samples[NSAMPLES] = {
 Data,
 TTbar, ZJets, DY1050 , //WJets,
  T_tW,    Tbar_tW, T_t,     Tbar_t, T_s,   Tbar_s,
 Tprime800_bWbW, Tprime800_bWtH, Tprime800_tZtZ, Tprime800_tHtH, Tprime800_bWtZ, Tprime800_tHtZ
// // Tprime500_bWbW, Tprime500_bWtH, Tprime500_tZtZ, Tprime500_tHtH, Tprime500_bWtZ, Tprime500_tHtZ
 };

/*
//This is for SS:
const unsigned int NSAMPLES = 15;
Samples_t samples[NSAMPLES] = {
Data, 
WZ, ZZ,
WWm, WWp,     WWW,
TTW,      TTZ, ChargeMisID, FakeRate, //, TTZ FakeRate
Tprime800_bWtH, Tprime800_tZtZ, Tprime800_tHtH, Tprime800_bWtZ, Tprime800_tHtZ
//Tprime1000_bWtH, Tprime1000_tZtZ, Tprime1000_tHtH, Tprime1000_bWtZ, Tprime1000_tHtZ
};
*/
//This is for multi-lepton:
/*
 const unsigned int NSAMPLES = 20;
 Samples_t samples[NSAMPLES] = {
 Data, 
 WWm,     WWp,     WWW, TTW,      TTWW, ZZ, WZ,
 ZZZNoGs, WWZNoGs, WZZNoGs, TTbar, ZJets, DY1050,  TTZ, //FakeRate,
 Tprime800_bWtH, Tprime800_tZtZ, Tprime800_tHtH, Tprime800_bWtZ, Tprime800_tHtZ
//  Tprime1000_bWtH, Tprime1000_tZtZ, Tprime1000_tHtH, Tprime1000_bWtZ, Tprime1000_tHtZ
  };*/
// const unsigned int NSAMPLES = 11;
// Samples_t samples[NSAMPLES] = {
// Data,
// TTbar, ZJets, DY1050, WJets,
//  T_tW,    Tbar_tW, T_t,     Tbar_t, T_s,   Tbar_s
//  //WW,
// //  WZ,
// //      ZZ,
// // // // WGEl,    WGMu,
// // // For Tri: 16
// // Tprime500_tZtZ, Tprime500_bWtZ, Tprime500_tHtZ
// // WWm,     WWp,     WWW,
// // TTW,     TTZ};
// };
//


float dySF[4] = {1.13,1.,1.04,1.1};


vector <TString> sNames;

vector <double> xSection;
vector <TH1D*>  Histos;
vector <double> sScale;
vector <double> systUnc;


/**
#DrawOneHisto parameter list:
#1: Input directory in root file
#2: Channel (ElEl, ElMu or MuMu
#3: Prefix for histograms
#4: Name of histogram
#5: Log or not log
#6: x-Axis name
#7: Optional: rebin by this number
*/

float getMyEntries(TH1* histo)
{
  float tot=0.;
  for (int i=0;i<histo->GetNbinsX()+1;++i) {
    if (histo->GetBinContent(i)>0.) tot+=histo->GetBinContent(i);
  }
  return tot;
}

float getMinEntries(TH1* histo)
{
  float min=999999.;
  for (int i=0;i<histo->GetNbinsX()+1;++i) {
    if ((histo->GetBinContent(i)>0.) && (histo->GetBinContent(i)<min)) min=histo->GetBinContent(i);
  }
  return min;
}
float getMinEntries(TH1* histo, float min)
{
  for (int i=0;i<histo->GetNbinsX()+1;++i) {
    if ((histo->GetBinContent(i)>0.) && (histo->GetBinContent(i)<min)) min=histo->GetBinContent(i);
  }
  return min;
}

double uncertainty (Samples_t sampleIndex );

int main( int argc, const char* argv[] ){
  if (argc < 9){
    cout<<"Need at least 8 arguments. Only "<<argc<<" found."<<endl;
    return 1;
  }

  setTDRStyle();

  //Root file information
  TString rootFile  = argv[1];
  TString inDir  = argv[2];
  TString inSuff = argv[3];
  TString inPref = argv[4];
  TString inName = argv[5];
cout << "a\n";

  //Drawing information
  int yLog = atoi(argv[6]);
  TString xTitle = argv[7];
  TString yTitle  = argv[8];
cout << "a\n"<<argc<<endl;

  int rebin = 1;
 
  if (argc > 9) rebin = atoi(argv[9]);
  float minX = -999, maxX;
  if (argc > 10) {
    minX = atof(argv[10]);
    maxX = atof(argv[11]);
  }

  TString plotType;
  if(argc > 12) plotType = argv[12];
cout << "a\n";

  cout<<inDir<<" ";
  cout<<inSuff<<" ";
  cout<<inPref<<" ";
  cout<<inName<<" ";
  cout<<yLog<<" ";
  cout<<xTitle<<" ";
  cout<<rebin<<" ";
  cout<<minX<<" ";
  cout<<maxX<<endl;

  const unsigned int NCHAN = 4;

  TFile *file0 = TFile::Open(rootFile);
  if (file0==0) {
    cout << "\n\nERROR: "<<rootFile<< " does not exist\n\n";
    return 1;
  }
  file0->cd();

  TH1F* histArray[NSAMPLES][NCHAN];
  if (inPref != "Top") inPref = "_"+inPref;
  else                 inPref = "";
  if (inDir != "Top") file0->cd(inDir);
  for (unsigned int isam = 0; isam < NSAMPLES; isam++){
//   cout << "ab "<<samples[isam]<< indices.size()<< endl;
    indices[samples[isam]] = isam;
    sNames.push_back(allNames[samples[isam]]);

    TString sampleName;
    if (isam < NSAMPLES) sampleName = sNames[isam];
    else                 sampleName = "FakeRate";

    cout << isam<<" "<<allNames[samples[isam]]<<" "<<color[samples[isam]]<<" "<<inName<<inPref<<"_"<<sampleName<<endl;;

    histArray[isam][0] = (TH1F*)gDirectory->Get(inName+inPref+"_"+sampleName+"_"+"MuMu");
    histArray[isam][1] = (TH1F*)gDirectory->Get(inName+inPref+"_"+sampleName+"_"+"ElMu");
    histArray[isam][2] = (TH1F*)gDirectory->Get(inName+inPref+"_"+sampleName+"_"+"ElEl");
    if ((histArray[isam][0]==0) || (histArray[isam][0]==0) || (histArray[isam][0]==0)) {
      cout << "\n\nERROR: One of the histograms do not exist: "<< inName<<inPref<<"_"<<sampleName<<"\n\n";
      gDirectory->ls();
      return 1;
    }
// cout <<  histArray[isam][0]->GetName()<<endl;
// cout <<  histArray[isam][0]->GetNbinsX()<<endl;
    if (rebin != 1){

      if(plotType=="Thesis"){

       histArray[isam][0]->Rebin(rebin);
       histArray[isam][1]->Rebin(rebin);
       histArray[isam][2]->Rebin(rebin);


       }

     if(plotType=="MET_Fig2"){
       Double_t rebin_array[27] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400, 420, 440, 460, 480, 500, 600};
       histArray[isam][0] = (TH1F*)histArray[isam][0]->Rebin(26, "temp1", rebin_array);  
       histArray[isam][0]->SetBinError(26, histArray[isam][0]->GetBinError(26)/5.);
       histArray[isam][0]->SetBinContent(26, histArray[isam][0]->GetBinContent(26)/5.);
       histArray[isam][1] = (TH1F*)histArray[isam][1]->Rebin(26, "temp2", rebin_array); 
       histArray[isam][1]->SetBinError(26, histArray[isam][1]->GetBinError(26)/5.);
       histArray[isam][1]->SetBinContent(26, histArray[isam][1]->GetBinContent(26)/5.);
       histArray[isam][2] = (TH1F*)histArray[isam][2]->Rebin(26, "temp3", rebin_array); 
       histArray[isam][2]->SetBinError(26, histArray[isam][2]->GetBinError(26)/5.);
       histArray[isam][2]->SetBinContent(26, histArray[isam][2]->GetBinContent(26)/5.);


     }

    if(plotType=="MET_Fig3"){
       Double_t rebin_array[26] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400, 420, 440, 460, 480, 600};
       histArray[isam][0] = (TH1F*)histArray[isam][0]->Rebin(25, "temp1", rebin_array);
       histArray[isam][0]->SetBinError(25, histArray[isam][0]->GetBinError(25)/6.);
       histArray[isam][0]->SetBinContent(25, histArray[isam][0]->GetBinContent(25)/6.);
       histArray[isam][1] = (TH1F*)histArray[isam][1]->Rebin(25, "temp2", rebin_array);
       histArray[isam][1]->SetBinError(25, histArray[isam][1]->GetBinError(25)/6.);
       histArray[isam][1]->SetBinContent(25, histArray[isam][1]->GetBinContent(25)/6.);
       histArray[isam][2] = (TH1F*)histArray[isam][2]->Rebin(25, "temp3", rebin_array);
       histArray[isam][2]->SetBinError(25, histArray[isam][2]->GetBinError(25)/6.);
       histArray[isam][2]->SetBinContent(25, histArray[isam][2]->GetBinContent(25)/6.);


     }


    if(plotType=="Jet1Pt_Fig3"){

       histArray[isam][0]->Rebin(rebin);
       histArray[isam][1]->Rebin(rebin);
       histArray[isam][2]->Rebin(rebin);

     }

    if(plotType=="Jet2Pt_Fig3"){
      Double_t rebin_array[15] = {0.0, 50.0, 100.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0, 550.0, 600.0, 800.0, 1000.0};
      histArray[isam][0] = (TH1F*)histArray[isam][0]->Rebin(14, "temp1", rebin_array);
      histArray[isam][0]->SetBinError(13, histArray[isam][0]->GetBinError(13)/4.);
      histArray[isam][0]->SetBinContent(13, histArray[isam][0]->GetBinContent(13)/4.);
      histArray[isam][0]->SetBinError(14, histArray[isam][0]->GetBinError(14)/4.);
      histArray[isam][0]->SetBinContent(14, histArray[isam][0]->GetBinContent(14)/4.);
      histArray[isam][1] = (TH1F*)histArray[isam][1]->Rebin(14, "temp2", rebin_array);
      histArray[isam][1]->SetBinError(13, histArray[isam][1]->GetBinError(13)/4.);
      histArray[isam][1]->SetBinContent(13, histArray[isam][1]->GetBinContent(13)/4.);
      histArray[isam][1]->SetBinError(14, histArray[isam][1]->GetBinError(14)/4.);
      histArray[isam][1]->SetBinContent(14, histArray[isam][1]->GetBinContent(14)/4.);
      histArray[isam][2] = (TH1F*)histArray[isam][2]->Rebin(14, "temp3", rebin_array);
      histArray[isam][2]->SetBinError(13, histArray[isam][2]->GetBinError(13)/4.);
      histArray[isam][2]->SetBinContent(13, histArray[isam][2]->GetBinContent(13)/4.);
      histArray[isam][2]->SetBinError(14, histArray[isam][2]->GetBinError(14)/4.);
      histArray[isam][2]->SetBinContent(14, histArray[isam][2]->GetBinContent(14)/4.);

       }

    if(plotType=="Lep1Pt_Fig4"){
   /* 
     histArray[isam][0]->Rebin(rebin);
     histArray[isam][1]->Rebin(rebin);
     histArray[isam][2]->Rebin(rebin);
*/
  
     Double_t rebin_array[18] = {0, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 400, 450, 500};

     histArray[isam][0] = (TH1F*)histArray[isam][0]->Rebin(17, "temp1", rebin_array);
     histArray[isam][0]->SetBinError(15, histArray[isam][0]->GetBinError(15)/2.);
     histArray[isam][0]->SetBinContent(15, histArray[isam][0]->GetBinContent(15)/2.);
     histArray[isam][0]->SetBinError(16, histArray[isam][0]->GetBinError(16)/2.);
     histArray[isam][0]->SetBinContent(16, histArray[isam][0]->GetBinContent(16)/2.);
     histArray[isam][0]->SetBinError(17, histArray[isam][0]->GetBinError(17)/2.);
     histArray[isam][0]->SetBinContent(17, histArray[isam][0]->GetBinContent(17)/2.);

     histArray[isam][1] = (TH1F*)histArray[isam][1]->Rebin(17, "temp1", rebin_array);
     histArray[isam][1]->SetBinError(15, histArray[isam][1]->GetBinError(15)/2.);
     histArray[isam][1]->SetBinContent(15, histArray[isam][1]->GetBinContent(15)/2.);
     histArray[isam][1]->SetBinError(16, histArray[isam][1]->GetBinError(16)/2.);
     histArray[isam][1]->SetBinContent(16, histArray[isam][1]->GetBinContent(16)/2.);
     histArray[isam][1]->SetBinError(17, histArray[isam][1]->GetBinError(17)/2.);
     histArray[isam][1]->SetBinContent(17, histArray[isam][1]->GetBinContent(17)/2.);

     histArray[isam][2] = (TH1F*)histArray[isam][2]->Rebin(17, "temp1", rebin_array);
     histArray[isam][2]->SetBinError(15, histArray[isam][2]->GetBinError(15)/2.);
     histArray[isam][2]->SetBinContent(15, histArray[isam][2]->GetBinContent(15)/2.);
     histArray[isam][2]->SetBinError(16, histArray[isam][2]->GetBinError(16)/2.);
     histArray[isam][2]->SetBinContent(16, histArray[isam][2]->GetBinContent(16)/2.);
     histArray[isam][2]->SetBinError(17, histArray[isam][2]->GetBinError(17)/2.);
     histArray[isam][2]->SetBinContent(17, histArray[isam][2]->GetBinContent(17)/2.);



/*
      Double_t rebin_array[24] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400, 440, 480, 520};

      histArray[isam][0] = (TH1F*)histArray[isam][0]->Rebin(23, "temp1", rebin_array);
      histArray[isam][0]->SetBinError(21, histArray[isam][0]->GetBinError(21)/2.);
      histArray[isam][0]->SetBinContent(21, histArray[isam][0]->GetBinContent(21)/2.);
      histArray[isam][0]->SetBinError(22, histArray[isam][0]->GetBinError(22)/2.);
      histArray[isam][0]->SetBinContent(22, histArray[isam][0]->GetBinContent(22)/2.);
      histArray[isam][0]->SetBinError(23, histArray[isam][0]->GetBinError(23)/2.);
      histArray[isam][0]->SetBinContent(23, histArray[isam][0]->GetBinContent(23)/2.);

      histArray[isam][1] = (TH1F*)histArray[isam][1]->Rebin(23, "temp2", rebin_array);
      histArray[isam][1]->SetBinError(21, histArray[isam][1]->GetBinError(21)/2.);
      histArray[isam][1]->SetBinContent(21, histArray[isam][1]->GetBinContent(21)/2.);
      histArray[isam][1]->SetBinError(22, histArray[isam][1]->GetBinError(22)/2.);
      histArray[isam][1]->SetBinContent(22, histArray[isam][1]->GetBinContent(22)/2.);
      histArray[isam][1]->SetBinError(23, histArray[isam][1]->GetBinError(23)/2.);
      histArray[isam][1]->SetBinContent(23, histArray[isam][1]->GetBinContent(23)/2.);

      histArray[isam][2] = (TH1F*)histArray[isam][2]->Rebin(23, "temp3", rebin_array);
      histArray[isam][2]->SetBinError(21, histArray[isam][2]->GetBinError(21)/2.);
      histArray[isam][2]->SetBinContent(21, histArray[isam][2]->GetBinContent(21)/2.);
      histArray[isam][2]->SetBinError(22, histArray[isam][2]->GetBinError(22)/2.);
      histArray[isam][2]->SetBinContent(22, histArray[isam][2]->GetBinContent(22)/2.);
      histArray[isam][2]->SetBinError(23, histArray[isam][2]->GetBinError(23)/2.);
      histArray[isam][2]->SetBinContent(23, histArray[isam][2]->GetBinContent(23)/2.);
*/
      }

   if(plotType=="Lep2Pt_Fig4"){ 
    Double_t rebin_array[19] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 400};
    histArray[isam][0] = (TH1F*)histArray[isam][0]->Rebin(18, "temp1", rebin_array);
    histArray[isam][0]->SetBinError(18, histArray[isam][0]->GetBinError(18)/3.);
    histArray[isam][0]->SetBinContent(18, histArray[isam][0]->GetBinContent(18)/3.);
    histArray[isam][1] = (TH1F*)histArray[isam][1]->Rebin(18, "temp2", rebin_array);
    histArray[isam][1]->SetBinError(18, histArray[isam][1]->GetBinError(18)/3.);
    histArray[isam][1]->SetBinContent(18, histArray[isam][1]->GetBinContent(18)/3.);
    histArray[isam][2] = (TH1F*)histArray[isam][2]->Rebin(18, "temp3", rebin_array);
    histArray[isam][2]->SetBinError(18, histArray[isam][2]->GetBinError(18)/3.);
    histArray[isam][2]->SetBinContent(18, histArray[isam][2]->GetBinContent(18)/3.);

     }

   if(plotType=="minMlb_Fig4"){
     Double_t rebin_array[19] = {0, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375, 400, 600, 900};    //rebinning for sum lepton pt in OS 
     histArray[isam][0] = (TH1F*)histArray[isam][0]->Rebin(18, "temp1", rebin_array);
     histArray[isam][0]->SetBinError(17, histArray[isam][0]->GetBinError(17)/8.);
     histArray[isam][0]->SetBinContent(17, histArray[isam][0]->GetBinContent(17)/8.);
     histArray[isam][0]->SetBinError(18, histArray[isam][0]->GetBinError(18)/12.);
     histArray[isam][0]->SetBinContent(18, histArray[isam][0]->GetBinContent(18)/12.);
     histArray[isam][1] = (TH1F*)histArray[isam][1]->Rebin(18, "temp2", rebin_array);
     histArray[isam][1]->SetBinError(17, histArray[isam][1]->GetBinError(17)/8.);
     histArray[isam][1]->SetBinContent(17, histArray[isam][1]->GetBinContent(17)/8.);
     histArray[isam][1]->SetBinError(18, histArray[isam][1]->GetBinError(18)/12.);
     histArray[isam][1]->SetBinContent(18, histArray[isam][1]->GetBinContent(18)/12.);
     histArray[isam][2] = (TH1F*)histArray[isam][2]->Rebin(18, "temp3", rebin_array);
     histArray[isam][2]->SetBinError(17, histArray[isam][2]->GetBinError(17)/8.);
     histArray[isam][2]->SetBinContent(17, histArray[isam][2]->GetBinContent(17)/8.);
     histArray[isam][2]->SetBinError(18, histArray[isam][2]->GetBinError(18)/12.);
     histArray[isam][2]->SetBinContent(18, histArray[isam][2]->GetBinContent(18)/12.);

    }

   if(plotType=="HT_Fig5"){
     Double_t rebin_array[16] = {0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1200, 1400, 1600, 2000, 2500}; //HT
     histArray[isam][0] = (TH1F*)histArray[isam][0]->Rebin(15, "temp1", rebin_array);
     histArray[isam][0]->SetBinError(12, histArray[isam][0]->GetBinError(12)/2.);
     histArray[isam][0]->SetBinContent(12, histArray[isam][0]->GetBinContent(12)/2.);  
     histArray[isam][0]->SetBinError(13, histArray[isam][0]->GetBinError(13)/2.);
     histArray[isam][0]->SetBinContent(13, histArray[isam][0]->GetBinContent(13)/2.);
     histArray[isam][0]->SetBinError(14, histArray[isam][0]->GetBinError(14)/4.);
     histArray[isam][0]->SetBinContent(14, histArray[isam][0]->GetBinContent(14)/4.);
     histArray[isam][0]->SetBinError(15, histArray[isam][0]->GetBinError(15)/5.);
     histArray[isam][0]->SetBinContent(15, histArray[isam][0]->GetBinContent(15)/5.);
     histArray[isam][1] = (TH1F*)histArray[isam][1]->Rebin(15, "temp2", rebin_array);
     histArray[isam][1]->SetBinError(12, histArray[isam][1]->GetBinError(12)/2.);
     histArray[isam][1]->SetBinContent(12, histArray[isam][1]->GetBinContent(12)/2.);
     histArray[isam][1]->SetBinError(13, histArray[isam][1]->GetBinError(13)/2.);
     histArray[isam][1]->SetBinContent(13, histArray[isam][1]->GetBinContent(13)/2.);
     histArray[isam][1]->SetBinError(14, histArray[isam][1]->GetBinError(14)/4.);
     histArray[isam][1]->SetBinContent(14, histArray[isam][1]->GetBinContent(14)/4.);
     histArray[isam][1]->SetBinError(15, histArray[isam][1]->GetBinError(15)/5.);
     histArray[isam][1]->SetBinContent(15, histArray[isam][1]->GetBinContent(15)/5.);
     histArray[isam][2] = (TH1F*)histArray[isam][2]->Rebin(15, "temp3", rebin_array);
     histArray[isam][2]->SetBinError(12, histArray[isam][2]->GetBinError(12)/2.);
     histArray[isam][2]->SetBinContent(12, histArray[isam][2]->GetBinContent(12)/2.);
     histArray[isam][2]->SetBinError(13, histArray[isam][2]->GetBinError(13)/2.);
     histArray[isam][2]->SetBinContent(13, histArray[isam][2]->GetBinContent(13)/2.);
     histArray[isam][2]->SetBinError(14, histArray[isam][2]->GetBinError(14)/4.);
     histArray[isam][2]->SetBinContent(14, histArray[isam][2]->GetBinContent(14)/4.);
     histArray[isam][2]->SetBinError(15, histArray[isam][2]->GetBinError(15)/5.);
     histArray[isam][2]->SetBinContent(15, histArray[isam][2]->GetBinContent(15)/5.);


     }

   if(plotType=="sumPtL_Fig5"){
     Double_t rebin_array[20] = {0, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 350, 400, 500, 600, 700, 800, 900};    //rebinning for sum lepton pt 
     histArray[isam][0] = (TH1F*)histArray[isam][0]->Rebin(19, "temp1", rebin_array);
     histArray[isam][0]->SetBinError(13, histArray[isam][0]->GetBinError(13)/2.);
     histArray[isam][0]->SetBinContent(13, histArray[isam][0]->GetBinContent(13)/2.);  
     histArray[isam][0]->SetBinError(14, histArray[isam][0]->GetBinError(14)/2.);
     histArray[isam][0]->SetBinContent(14, histArray[isam][0]->GetBinContent(14)/2.);
     histArray[isam][0]->SetBinError(15, histArray[isam][0]->GetBinError(15)/4.);
     histArray[isam][0]->SetBinContent(15, histArray[isam][0]->GetBinContent(15)/4.);
     histArray[isam][0]->SetBinError(16, histArray[isam][0]->GetBinError(16)/4.);
     histArray[isam][0]->SetBinContent(16, histArray[isam][0]->GetBinContent(16)/4.);
     histArray[isam][0]->SetBinError(17, histArray[isam][0]->GetBinError(17)/4.);
     histArray[isam][0]->SetBinContent(17, histArray[isam][0]->GetBinContent(17)/4.);
     histArray[isam][0]->SetBinError(18, histArray[isam][0]->GetBinError(18)/4.);
     histArray[isam][0]->SetBinContent(18, histArray[isam][0]->GetBinContent(18)/4.);
     histArray[isam][0]->SetBinError(19, histArray[isam][0]->GetBinError(19)/4.);
     histArray[isam][0]->SetBinContent(19, histArray[isam][0]->GetBinContent(19)/4.);
     histArray[isam][1] = (TH1F*)histArray[isam][1]->Rebin(19, "temp2", rebin_array);
     histArray[isam][1]->SetBinError(13, histArray[isam][1]->GetBinError(13)/2.);
     histArray[isam][1]->SetBinContent(13, histArray[isam][1]->GetBinContent(13)/2.);
     histArray[isam][1]->SetBinError(14, histArray[isam][1]->GetBinError(14)/2.);
     histArray[isam][1]->SetBinContent(14, histArray[isam][1]->GetBinContent(14)/2.);
     histArray[isam][1]->SetBinError(15, histArray[isam][1]->GetBinError(15)/4.);
     histArray[isam][1]->SetBinContent(15, histArray[isam][1]->GetBinContent(15)/4.);
     histArray[isam][1]->SetBinError(16, histArray[isam][1]->GetBinError(16)/4.);
     histArray[isam][1]->SetBinContent(16, histArray[isam][1]->GetBinContent(16)/4.);
     histArray[isam][1]->SetBinError(17, histArray[isam][1]->GetBinError(17)/4.);
     histArray[isam][1]->SetBinContent(17, histArray[isam][1]->GetBinContent(17)/4.);
     histArray[isam][1]->SetBinError(18, histArray[isam][1]->GetBinError(18)/4.);
     histArray[isam][1]->SetBinContent(18, histArray[isam][1]->GetBinContent(18)/4.);
     histArray[isam][1]->SetBinError(19, histArray[isam][1]->GetBinError(19)/4.);
     histArray[isam][1]->SetBinContent(19, histArray[isam][1]->GetBinContent(19)/4.);
     histArray[isam][2] = (TH1F*)histArray[isam][2]->Rebin(19, "temp3", rebin_array);
     histArray[isam][2]->SetBinError(13, histArray[isam][2]->GetBinError(13)/2.);
     histArray[isam][2]->SetBinContent(13, histArray[isam][2]->GetBinContent(13)/2.);
     histArray[isam][2]->SetBinError(14, histArray[isam][2]->GetBinError(14)/2.);
     histArray[isam][2]->SetBinContent(14, histArray[isam][2]->GetBinContent(14)/2.);
     histArray[isam][2]->SetBinError(15, histArray[isam][2]->GetBinError(15)/4.);
     histArray[isam][2]->SetBinContent(15, histArray[isam][2]->GetBinContent(15)/4.);
     histArray[isam][2]->SetBinError(16, histArray[isam][2]->GetBinError(16)/4.);
     histArray[isam][2]->SetBinContent(16, histArray[isam][2]->GetBinContent(16)/4.);
     histArray[isam][2]->SetBinError(17, histArray[isam][2]->GetBinError(17)/4.);
     histArray[isam][2]->SetBinContent(17, histArray[isam][2]->GetBinContent(17)/4.);
     histArray[isam][2]->SetBinError(18, histArray[isam][2]->GetBinError(18)/4.);
     histArray[isam][2]->SetBinContent(18, histArray[isam][2]->GetBinContent(18)/4.);
     histArray[isam][2]->SetBinError(19, histArray[isam][2]->GetBinError(19)/4.);
     histArray[isam][2]->SetBinContent(19, histArray[isam][2]->GetBinContent(19)/4.);


}

  if(plotType=="ST_Fig5" || plotType=="leptonJetsSum_Fig5"){


     histArray[isam][0]->Rebin(rebin);
     histArray[isam][1]->Rebin(rebin);
     histArray[isam][2]->Rebin(rebin);



     //Double_t rebin_array[38] = {0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200, 1250, 1300, 1350, 1400, 1450, 1500, 1550, 1600, 1650, 1700, 1800, 1900, 2000}; //rebinning ST
     /*Double_t rebin_array[38] = {0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200, 1250, 1300, 1350, 1400, 1450, 1500, 1550, 1600, 1650, 1700, 1800, 1900, 2000};
     histArray[isam][0] = (TH1F*)histArray[isam][0]->Rebin(37, "temp1", rebin_array);
     histArray[isam][0]->SetBinError(35, histArray[isam][0]->GetBinError(35)/2.);
     histArray[isam][0]->SetBinContent(35, histArray[isam][0]->GetBinContent(35)/2.);
     histArray[isam][0]->SetBinError(36, histArray[isam][0]->GetBinError(36)/2.);
     histArray[isam][0]->SetBinContent(36, histArray[isam][0]->GetBinContent(36)/2.);
     histArray[isam][0]->SetBinError(37, histArray[isam][0]->GetBinError(37)/2.);
     histArray[isam][0]->SetBinContent(37, histArray[isam][0]->GetBinContent(37)/2.);
     histArray[isam][1] = (TH1F*)histArray[isam][1]->Rebin(37, "temp2", rebin_array);
     histArray[isam][1]->SetBinError(35, histArray[isam][1]->GetBinError(35)/2.);
     histArray[isam][1]->SetBinContent(35, histArray[isam][1]->GetBinContent(35)/2.);
     histArray[isam][1]->SetBinError(36, histArray[isam][1]->GetBinError(36)/2.);
     histArray[isam][1]->SetBinContent(36, histArray[isam][1]->GetBinContent(36)/2.);
     histArray[isam][1]->SetBinError(37, histArray[isam][1]->GetBinError(37)/2.);
     histArray[isam][1]->SetBinContent(37, histArray[isam][1]->GetBinContent(37)/2.);
     histArray[isam][2] = (TH1F*)histArray[isam][2]->Rebin(37, "temp3", rebin_array);
     histArray[isam][2]->SetBinError(35, histArray[isam][2]->GetBinError(35)/2.);
     histArray[isam][2]->SetBinContent(35, histArray[isam][2]->GetBinContent(35)/2.);
     histArray[isam][2]->SetBinError(36, histArray[isam][2]->GetBinError(36)/2.);
     histArray[isam][2]->SetBinContent(36, histArray[isam][2]->GetBinContent(36)/2.);
     histArray[isam][2]->SetBinError(37, histArray[isam][2]->GetBinError(37)/2.);
     histArray[isam][2]->SetBinContent(37, histArray[isam][2]->GetBinContent(37)/2.);*/

}

  if(plotType=="leptonJetsSum_Fig6"){
 
    Double_t rebin_array[23] = {0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1500, 2000}; //rebinning ST
    histArray[isam][0] = (TH1F*)histArray[isam][0]->Rebin(22, "temp1", rebin_array);
    histArray[isam][0]->SetBinError(21, histArray[isam][0]->GetBinError(21)/10.);
    histArray[isam][0]->SetBinContent(21, histArray[isam][0]->GetBinContent(21)/10.);
    histArray[isam][0]->SetBinError(22, histArray[isam][0]->GetBinError(22)/10.);
    histArray[isam][0]->SetBinContent(22, histArray[isam][0]->GetBinContent(22)/10.);
    histArray[isam][1] = (TH1F*)histArray[isam][1]->Rebin(22, "temp2", rebin_array);
    histArray[isam][1]->SetBinError(21, histArray[isam][1]->GetBinError(21)/10.);
    histArray[isam][1]->SetBinContent(21, histArray[isam][1]->GetBinContent(21)/10.);
    histArray[isam][1]->SetBinError(22, histArray[isam][1]->GetBinError(22)/10.);
    histArray[isam][1]->SetBinContent(22, histArray[isam][1]->GetBinContent(22)/10.);
    histArray[isam][2] = (TH1F*)histArray[isam][2]->Rebin(22, "temp3", rebin_array);
    histArray[isam][2]->SetBinError(21, histArray[isam][2]->GetBinError(21)/10.);
    histArray[isam][2]->SetBinContent(21, histArray[isam][2]->GetBinContent(21)/10.);
    histArray[isam][2]->SetBinError(22, histArray[isam][2]->GetBinError(22)/10.);
    histArray[isam][2]->SetBinContent(22, histArray[isam][2]->GetBinContent(22)/10.);

  }


   if(plotType=="MET_Fig13"){

     Double_t rebin_array[18] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 400, 600};
     histArray[isam][0] = (TH1F*)histArray[isam][0]->Rebin(17, "temp1", rebin_array);
     histArray[isam][0]->SetBinError(16, histArray[isam][0]->GetBinError(16)/5.);
     histArray[isam][0]->SetBinContent(16, histArray[isam][0]->GetBinContent(16)/5.);
     histArray[isam][0]->SetBinError(17, histArray[isam][0]->GetBinError(17)/5.);
     histArray[isam][0]->SetBinContent(17, histArray[isam][0]->GetBinContent(17)/5.);
     histArray[isam][1] = (TH1F*)histArray[isam][1]->Rebin(17, "temp2", rebin_array);
     histArray[isam][1]->SetBinError(16, histArray[isam][1]->GetBinError(16)/5.);
     histArray[isam][1]->SetBinContent(16, histArray[isam][1]->GetBinContent(16)/5.);
     histArray[isam][1]->SetBinError(17, histArray[isam][1]->GetBinError(17)/5.);
     histArray[isam][1]->SetBinContent(17, histArray[isam][1]->GetBinContent(17)/5.);
     histArray[isam][2] = (TH1F*)histArray[isam][2]->Rebin(17, "temp3", rebin_array);
     histArray[isam][2]->SetBinError(16, histArray[isam][2]->GetBinError(16)/5.);
     histArray[isam][2]->SetBinContent(16, histArray[isam][2]->GetBinContent(16)/5.);
     histArray[isam][2]->SetBinError(17, histArray[isam][2]->GetBinError(17)/5.);
     histArray[isam][2]->SetBinContent(17, histArray[isam][2]->GetBinContent(17)/5.);

     }

   if(plotType=="Jet1Pt_Fig14"){
      Double_t rebin_array[14] = {0.0, 50.0, 100.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0, 600.0, 700.0, 1000.0};
      histArray[isam][0] = (TH1F*)histArray[isam][0]->Rebin(13, "temp1", rebin_array);
      histArray[isam][0]->SetBinError(11, histArray[isam][0]->GetBinError(11)/2.);
      histArray[isam][0]->SetBinContent(11, histArray[isam][0]->GetBinContent(11)/2.);
      histArray[isam][0]->SetBinError(12, histArray[isam][0]->GetBinError(12)/2.);
      histArray[isam][0]->SetBinContent(12, histArray[isam][0]->GetBinContent(12)/2.);
      histArray[isam][0]->SetBinError(13, histArray[isam][0]->GetBinError(13)/6.);
      histArray[isam][0]->SetBinContent(13, histArray[isam][0]->GetBinContent(13)/6.);
      histArray[isam][1] = (TH1F*)histArray[isam][1]->Rebin(13, "temp2", rebin_array);
      histArray[isam][1]->SetBinError(11, histArray[isam][1]->GetBinError(11)/2.);
      histArray[isam][1]->SetBinContent(11, histArray[isam][1]->GetBinContent(11)/2.);
      histArray[isam][1]->SetBinError(12, histArray[isam][1]->GetBinError(12)/2.);
      histArray[isam][1]->SetBinContent(12, histArray[isam][1]->GetBinContent(12)/2.);
      histArray[isam][1]->SetBinError(13, histArray[isam][1]->GetBinError(13)/6.);
      histArray[isam][1]->SetBinContent(13, histArray[isam][1]->GetBinContent(13)/6.);
      histArray[isam][2] = (TH1F*)histArray[isam][2]->Rebin(13, "temp3", rebin_array);
      histArray[isam][2]->SetBinError(11, histArray[isam][2]->GetBinError(11)/2.);
      histArray[isam][2]->SetBinContent(11, histArray[isam][2]->GetBinContent(11)/2.);
      histArray[isam][2]->SetBinError(12, histArray[isam][2]->GetBinError(12)/2.);
      histArray[isam][2]->SetBinContent(12, histArray[isam][2]->GetBinContent(12)/2.);
      histArray[isam][2]->SetBinError(13, histArray[isam][2]->GetBinError(13)/6.);
      histArray[isam][2]->SetBinContent(13, histArray[isam][2]->GetBinContent(13)/6.);

       }

   if(plotType=="Jet2Pt_Fig14"){
      Double_t rebin_array[9] = {0.0, 50.0, 100.0, 150.0, 200.0, 250.0, 300.0, 350.0, 500.0};
      histArray[isam][0] = (TH1F*)histArray[isam][0]->Rebin(8, "temp1", rebin_array);
      histArray[isam][0]->SetBinError(8, histArray[isam][0]->GetBinError(8)/3.);
      histArray[isam][0]->SetBinContent(8, histArray[isam][0]->GetBinContent(8)/3.);
      histArray[isam][1] = (TH1F*)histArray[isam][1]->Rebin(8, "temp2", rebin_array);
      histArray[isam][1]->SetBinError(8, histArray[isam][1]->GetBinError(8)/3.);
      histArray[isam][1]->SetBinContent(8, histArray[isam][1]->GetBinContent(8)/3.);
      histArray[isam][2] = (TH1F*)histArray[isam][2]->Rebin(8, "temp3", rebin_array);
      histArray[isam][2]->SetBinError(8, histArray[isam][2]->GetBinError(8)/3.);
      histArray[isam][2]->SetBinContent(8, histArray[isam][2]->GetBinContent(8)/3.);

     }


   if(plotType=="Lep1Pt_Fig14"){

    Double_t rebin_array[16] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 300, 400, 500};
    histArray[isam][0] = (TH1F*)histArray[isam][0]->Rebin(15, "temp1", rebin_array);
    histArray[isam][0]->SetBinError(13, histArray[isam][0]->GetBinError(13)/3.);
    histArray[isam][0]->SetBinContent(13, histArray[isam][0]->GetBinContent(13)/3.);
    histArray[isam][0]->SetBinError(14, histArray[isam][0]->GetBinError(14)/5.);
    histArray[isam][0]->SetBinContent(14, histArray[isam][0]->GetBinContent(14)/5.);
    histArray[isam][0]->SetBinError(15, histArray[isam][0]->GetBinError(15)/5.);
    histArray[isam][0]->SetBinContent(15, histArray[isam][0]->GetBinContent(15)/5.);
    histArray[isam][1] = (TH1F*)histArray[isam][1]->Rebin(15, "temp2", rebin_array);
    histArray[isam][1]->SetBinError(13, histArray[isam][1]->GetBinError(13)/3.);
    histArray[isam][1]->SetBinContent(13, histArray[isam][1]->GetBinContent(13)/3.);
    histArray[isam][1]->SetBinError(14, histArray[isam][1]->GetBinError(14)/5.);
    histArray[isam][1]->SetBinContent(14, histArray[isam][1]->GetBinContent(14)/5.);
    histArray[isam][1]->SetBinError(15, histArray[isam][1]->GetBinError(15)/5.);
    histArray[isam][1]->SetBinContent(15, histArray[isam][1]->GetBinContent(15)/5.);
    histArray[isam][2] = (TH1F*)histArray[isam][2]->Rebin(15, "temp3", rebin_array);
    histArray[isam][2]->SetBinError(13, histArray[isam][2]->GetBinError(13)/3.);
    histArray[isam][2]->SetBinContent(13, histArray[isam][2]->GetBinContent(13)/3.);
    histArray[isam][2]->SetBinError(14, histArray[isam][2]->GetBinError(14)/5.);
    histArray[isam][2]->SetBinContent(14, histArray[isam][2]->GetBinContent(14)/5.);
    histArray[isam][2]->SetBinError(15, histArray[isam][2]->GetBinError(15)/5.);
    histArray[isam][2]->SetBinContent(15, histArray[isam][2]->GetBinContent(15)/5.);


    } 

   if(plotType=="Lep2Pt_Fig14"){

    Double_t rebin_array[10] = {0, 20, 40, 60, 80, 100, 120, 140, 180, 260};
    histArray[isam][0] = (TH1F*)histArray[isam][0]->Rebin(9, "temp1", rebin_array);
    histArray[isam][0]->SetBinError(8, histArray[isam][0]->GetBinError(8)/2.);
    histArray[isam][0]->SetBinContent(8, histArray[isam][0]->GetBinContent(8)/2.);
    histArray[isam][0]->SetBinError(9, histArray[isam][0]->GetBinError(9)/4.);
    histArray[isam][0]->SetBinContent(9, histArray[isam][0]->GetBinContent(9)/4.);
    histArray[isam][1] = (TH1F*)histArray[isam][1]->Rebin(9, "temp2", rebin_array);
    histArray[isam][1]->SetBinError(8, histArray[isam][1]->GetBinError(8)/2.);
    histArray[isam][1]->SetBinContent(8, histArray[isam][1]->GetBinContent(8)/2.);
    histArray[isam][1]->SetBinError(9, histArray[isam][1]->GetBinError(9)/4.);
    histArray[isam][1]->SetBinContent(9, histArray[isam][1]->GetBinContent(9)/4.);
    histArray[isam][2] = (TH1F*)histArray[isam][2]->Rebin(9, "temp3", rebin_array);
    histArray[isam][2]->SetBinError(8, histArray[isam][2]->GetBinError(8)/2.);
    histArray[isam][2]->SetBinContent(8, histArray[isam][2]->GetBinContent(8)/2.);
    histArray[isam][2]->SetBinError(9, histArray[isam][2]->GetBinError(9)/4.);
    histArray[isam][2]->SetBinContent(9, histArray[isam][2]->GetBinContent(9)/4.);

    }

   if(plotType=="MET_Fig15"){

     Double_t rebin_array[12] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 200, 300, 400};
     histArray[isam][0] = (TH1F*)histArray[isam][0]->Rebin(11, "temp1", rebin_array);
     histArray[isam][0]->SetBinError(9, histArray[isam][0]->GetBinError(9)/2.);
     histArray[isam][0]->SetBinContent(9, histArray[isam][0]->GetBinContent(9)/2.);
     histArray[isam][0]->SetBinError(10, histArray[isam][0]->GetBinError(10)/5.);
     histArray[isam][0]->SetBinContent(10, histArray[isam][0]->GetBinContent(10)/5.);
     histArray[isam][0]->SetBinError(11, histArray[isam][0]->GetBinError(11)/5.);
     histArray[isam][0]->SetBinContent(11, histArray[isam][0]->GetBinContent(11)/5.);
     histArray[isam][1] = (TH1F*)histArray[isam][1]->Rebin(11, "temp2", rebin_array);
     histArray[isam][1]->SetBinError(9, histArray[isam][1]->GetBinError(9)/2.);
     histArray[isam][1]->SetBinContent(9, histArray[isam][1]->GetBinContent(9)/2.);
     histArray[isam][1]->SetBinError(10, histArray[isam][1]->GetBinError(10)/5.);
     histArray[isam][1]->SetBinContent(10, histArray[isam][1]->GetBinContent(10)/5.);
     histArray[isam][1]->SetBinError(11, histArray[isam][1]->GetBinError(11)/5.);
     histArray[isam][1]->SetBinContent(11, histArray[isam][1]->GetBinContent(11)/5.);
     histArray[isam][2] = (TH1F*)histArray[isam][2]->Rebin(11, "temp3", rebin_array);
     histArray[isam][2]->SetBinError(9, histArray[isam][2]->GetBinError(9)/2.);
     histArray[isam][2]->SetBinContent(9, histArray[isam][2]->GetBinContent(9)/2.);
     histArray[isam][2]->SetBinError(10, histArray[isam][2]->GetBinError(10)/5.);
     histArray[isam][2]->SetBinContent(10, histArray[isam][2]->GetBinContent(10)/5.);
     histArray[isam][2]->SetBinError(11, histArray[isam][2]->GetBinError(11)/5.);
     histArray[isam][2]->SetBinContent(11, histArray[isam][2]->GetBinContent(11)/5.);

     }

   if(plotType=="minMlb_Fig15"){
     
     Double_t rebin_array[11] = {0, 25, 50, 75, 100, 125, 150, 175, 200, 300, 500};    //rebinning for sum lepton pt in OS 
     histArray[isam][0] = (TH1F*)histArray[isam][0]->Rebin(10, "temp1", rebin_array);
     histArray[isam][0]->SetBinError(9, histArray[isam][0]->GetBinError(9)/4.);
     histArray[isam][0]->SetBinContent(9, histArray[isam][0]->GetBinContent(9)/4.);
     histArray[isam][0]->SetBinError(10, histArray[isam][0]->GetBinError(10)/8.);
     histArray[isam][0]->SetBinContent(10, histArray[isam][0]->GetBinContent(10)/8.);
     histArray[isam][1] = (TH1F*)histArray[isam][1]->Rebin(10, "temp2", rebin_array);
     histArray[isam][1]->SetBinError(9, histArray[isam][1]->GetBinError(9)/4.);
     histArray[isam][1]->SetBinContent(9, histArray[isam][1]->GetBinContent(9)/4.);
     histArray[isam][1]->SetBinError(10, histArray[isam][1]->GetBinError(10)/8.);
     histArray[isam][1]->SetBinContent(10, histArray[isam][1]->GetBinContent(10)/8.);
     histArray[isam][2] = (TH1F*)histArray[isam][2]->Rebin(10, "temp3", rebin_array);
     histArray[isam][2]->SetBinError(9, histArray[isam][2]->GetBinError(9)/4.);
     histArray[isam][2]->SetBinContent(9, histArray[isam][2]->GetBinContent(9)/4.);
     histArray[isam][2]->SetBinError(10, histArray[isam][2]->GetBinError(10)/8.);
     histArray[isam][2]->SetBinContent(10, histArray[isam][2]->GetBinContent(10)/8.);

    }

    if(plotType=="HT_Fig16"){
     Double_t rebin_array[11] = {0, 100, 200, 300, 400, 500, 600, 700, 800, 1000, 1600}; //HT
     histArray[isam][0] = (TH1F*)histArray[isam][0]->Rebin(10, "temp1", rebin_array);
     histArray[isam][0]->SetBinError(9, histArray[isam][0]->GetBinError(9)/2.);
     histArray[isam][0]->SetBinContent(9, histArray[isam][0]->GetBinContent(9)/2.);
     histArray[isam][0]->SetBinError(10, histArray[isam][0]->GetBinError(10)/6.);
     histArray[isam][0]->SetBinContent(10, histArray[isam][0]->GetBinContent(10)/6.);
     
     histArray[isam][1] = (TH1F*)histArray[isam][1]->Rebin(10, "temp2", rebin_array);
     histArray[isam][1]->SetBinError(9, histArray[isam][1]->GetBinError(9)/2.);
     histArray[isam][1]->SetBinContent(9, histArray[isam][1]->GetBinContent(9)/2.);
     histArray[isam][1]->SetBinError(10, histArray[isam][1]->GetBinError(10)/6.);
     histArray[isam][1]->SetBinContent(10, histArray[isam][1]->GetBinContent(10)/6.);
     
     histArray[isam][2] = (TH1F*)histArray[isam][2]->Rebin(10, "temp3", rebin_array);
     histArray[isam][2]->SetBinError(9, histArray[isam][2]->GetBinError(9)/2.);
     histArray[isam][2]->SetBinContent(9, histArray[isam][2]->GetBinContent(9)/2.);
     histArray[isam][2]->SetBinError(10, histArray[isam][2]->GetBinError(10)/6.);
     histArray[isam][2]->SetBinContent(10, histArray[isam][2]->GetBinContent(10)/6.);

     }

   if(plotType=="sumPtL_Fig16"){

     Double_t rebin_array[13] = {0, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 300, 500};    //rebinning for sum lepton pt 
     histArray[isam][0] = (TH1F*)histArray[isam][0]->Rebin(12, "temp1", rebin_array);
     histArray[isam][0]->SetBinError(11, histArray[isam][0]->GetBinError(11)/2.);
     histArray[isam][0]->SetBinContent(11, histArray[isam][0]->GetBinContent(11)/2.);
     histArray[isam][0]->SetBinError(12, histArray[isam][0]->GetBinError(12)/8.);
     histArray[isam][0]->SetBinContent(12, histArray[isam][0]->GetBinContent(12)/8.);
     histArray[isam][1] = (TH1F*)histArray[isam][1]->Rebin(12, "temp2", rebin_array);
     histArray[isam][1]->SetBinError(11, histArray[isam][1]->GetBinError(11)/2.);
     histArray[isam][1]->SetBinContent(11, histArray[isam][1]->GetBinContent(11)/2.);
     histArray[isam][1]->SetBinError(12, histArray[isam][1]->GetBinError(12)/8.);
     histArray[isam][1]->SetBinContent(12, histArray[isam][1]->GetBinContent(12)/8.);
     histArray[isam][2] = (TH1F*)histArray[isam][2]->Rebin(12, "temp3", rebin_array);
     histArray[isam][2]->SetBinError(11, histArray[isam][2]->GetBinError(11)/2.);
     histArray[isam][2]->SetBinContent(11, histArray[isam][2]->GetBinContent(11)/2.);
     histArray[isam][2]->SetBinError(12, histArray[isam][2]->GetBinError(12)/8.);
     histArray[isam][2]->SetBinContent(12, histArray[isam][2]->GetBinContent(12)/8.);

}

   if(plotType=="ST_Fig16"){

    Double_t rebin_array[17] = {200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 900, 1000, 1200, 1500};
    histArray[isam][0] = (TH1F*)histArray[isam][0]->Rebin(16, "temp1", rebin_array);
    histArray[isam][0]->SetBinError(13, histArray[isam][0]->GetBinError(13)/2.);
    histArray[isam][0]->SetBinContent(13, histArray[isam][0]->GetBinContent(13)/2.);
    histArray[isam][0]->SetBinError(14, histArray[isam][0]->GetBinError(14)/2.);
    histArray[isam][0]->SetBinContent(14, histArray[isam][0]->GetBinContent(14)/2.);
    histArray[isam][0]->SetBinError(15, histArray[isam][0]->GetBinError(15)/4.);
    histArray[isam][0]->SetBinContent(15, histArray[isam][0]->GetBinContent(15)/4.);
    histArray[isam][0]->SetBinError(16, histArray[isam][0]->GetBinError(16)/6.);
    histArray[isam][0]->SetBinContent(16, histArray[isam][0]->GetBinContent(16)/6.);
    
    histArray[isam][1] = (TH1F*)histArray[isam][1]->Rebin(16, "temp2", rebin_array);
    histArray[isam][1]->SetBinError(13, histArray[isam][1]->GetBinError(13)/2.);
    histArray[isam][1]->SetBinContent(13, histArray[isam][1]->GetBinContent(13)/2.);
    histArray[isam][1]->SetBinError(14, histArray[isam][1]->GetBinError(14)/2.);
    histArray[isam][1]->SetBinContent(14, histArray[isam][1]->GetBinContent(14)/2.);
    histArray[isam][1]->SetBinError(15, histArray[isam][1]->GetBinError(15)/4.);
    histArray[isam][1]->SetBinContent(15, histArray[isam][1]->GetBinContent(15)/4.);
    histArray[isam][1]->SetBinError(16, histArray[isam][1]->GetBinError(16)/6.);
    histArray[isam][1]->SetBinContent(16, histArray[isam][1]->GetBinContent(16)/6.);

    histArray[isam][2] = (TH1F*)histArray[isam][2]->Rebin(16, "temp3", rebin_array);
    histArray[isam][2]->SetBinError(13, histArray[isam][2]->GetBinError(13)/2.);
    histArray[isam][2]->SetBinContent(13, histArray[isam][2]->GetBinContent(13)/2.);
    histArray[isam][2]->SetBinError(14, histArray[isam][2]->GetBinError(14)/2.);
    histArray[isam][2]->SetBinContent(14, histArray[isam][2]->GetBinContent(14)/2.);
    histArray[isam][2]->SetBinError(15, histArray[isam][2]->GetBinError(15)/4.);
    histArray[isam][2]->SetBinContent(15, histArray[isam][2]->GetBinContent(15)/4.);
    histArray[isam][2]->SetBinError(16, histArray[isam][2]->GetBinError(16)/6.);
    histArray[isam][2]->SetBinContent(16, histArray[isam][2]->GetBinContent(16)/6.);

   }

 
  if(plotType=="leptonJetsSum_Fig16"){

     Double_t rebin_array[19] = {0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 1000, 1500}; //rebinning ST
     histArray[isam][0] = (TH1F*)histArray[isam][0]->Rebin(18, "temp1", rebin_array);
     histArray[isam][0]->SetBinError(17, histArray[isam][0]->GetBinError(17)/4.);
     histArray[isam][0]->SetBinContent(17, histArray[isam][0]->GetBinContent(17)/4.);
     histArray[isam][0]->SetBinError(18, histArray[isam][0]->GetBinError(18)/10.);
     histArray[isam][0]->SetBinContent(18, histArray[isam][0]->GetBinContent(18)/10.);
     histArray[isam][1] = (TH1F*)histArray[isam][1]->Rebin(18, "temp2", rebin_array);
     histArray[isam][1]->SetBinError(17, histArray[isam][1]->GetBinError(17)/4.);
     histArray[isam][1]->SetBinContent(17, histArray[isam][1]->GetBinContent(17)/4.);
     histArray[isam][1]->SetBinError(18, histArray[isam][1]->GetBinError(18)/10.);
     histArray[isam][1]->SetBinContent(18, histArray[isam][1]->GetBinContent(18)/10.);
     histArray[isam][2] = (TH1F*)histArray[isam][2]->Rebin(18, "temp3", rebin_array);
     histArray[isam][2]->SetBinError(17, histArray[isam][2]->GetBinError(17)/4.);
     histArray[isam][2]->SetBinContent(17, histArray[isam][2]->GetBinContent(17)/4.);
     histArray[isam][2]->SetBinError(18, histArray[isam][2]->GetBinError(18)/10.);
     histArray[isam][2]->SetBinContent(18, histArray[isam][2]->GetBinContent(18)/10.);

}

     if(plotType=="MET_Fig18"){

     Double_t rebin_array[18] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 400, 600};
     histArray[isam][0] = (TH1F*)histArray[isam][0]->Rebin(17, "temp1", rebin_array);
     histArray[isam][0]->SetBinError(16, histArray[isam][0]->GetBinError(16)/5.);
     histArray[isam][0]->SetBinContent(16, histArray[isam][0]->GetBinContent(16)/5.);
     histArray[isam][0]->SetBinError(17, histArray[isam][0]->GetBinError(17)/5.);
     histArray[isam][0]->SetBinContent(17, histArray[isam][0]->GetBinContent(17)/5.);
     histArray[isam][1] = (TH1F*)histArray[isam][1]->Rebin(17, "temp2", rebin_array);
     histArray[isam][1]->SetBinError(16, histArray[isam][1]->GetBinError(16)/5.);
     histArray[isam][1]->SetBinContent(16, histArray[isam][1]->GetBinContent(16)/5.);
     histArray[isam][1]->SetBinError(17, histArray[isam][1]->GetBinError(17)/5.);
     histArray[isam][1]->SetBinContent(17, histArray[isam][1]->GetBinContent(17)/5.);
     histArray[isam][2] = (TH1F*)histArray[isam][2]->Rebin(17, "temp3", rebin_array);
     histArray[isam][2]->SetBinError(16, histArray[isam][2]->GetBinError(16)/5.);
     histArray[isam][2]->SetBinContent(16, histArray[isam][2]->GetBinContent(16)/5.);
     histArray[isam][2]->SetBinError(17, histArray[isam][2]->GetBinError(17)/5.);
     histArray[isam][2]->SetBinContent(17, histArray[isam][2]->GetBinContent(17)/5.);

     }

    if(plotType=="Jet1Pt_Fig19"){
    
      Double_t rebin_array[6] = {0.0, 100.0, 200.0, 400.0, 600.0, 1000.0};
      histArray[isam][0] = (TH1F*)histArray[isam][0]->Rebin(5, "temp1", rebin_array);
      histArray[isam][0]->SetBinError(3, histArray[isam][0]->GetBinError(3)/2.);
      histArray[isam][0]->SetBinContent(3, histArray[isam][0]->GetBinContent(3)/2.);
      histArray[isam][0]->SetBinError(4, histArray[isam][0]->GetBinError(4)/2.);
      histArray[isam][0]->SetBinContent(4, histArray[isam][0]->GetBinContent(4)/2.);
      histArray[isam][0]->SetBinError(5, histArray[isam][0]->GetBinError(5)/4.);
      histArray[isam][0]->SetBinContent(5, histArray[isam][0]->GetBinContent(5)/4.);
      histArray[isam][1] = (TH1F*)histArray[isam][1]->Rebin(5, "temp2", rebin_array);
      histArray[isam][1]->SetBinError(3, histArray[isam][1]->GetBinError(3)/2.);
      histArray[isam][1]->SetBinContent(3, histArray[isam][1]->GetBinContent(3)/2.);
      histArray[isam][1]->SetBinError(4, histArray[isam][1]->GetBinError(4)/2.);
      histArray[isam][1]->SetBinContent(4, histArray[isam][1]->GetBinContent(4)/2.);
      histArray[isam][1]->SetBinError(5, histArray[isam][1]->GetBinError(5)/4.);
      histArray[isam][1]->SetBinContent(5, histArray[isam][1]->GetBinContent(5)/4.);
      histArray[isam][2] = (TH1F*)histArray[isam][2]->Rebin(5, "temp3", rebin_array);
      histArray[isam][2]->SetBinError(3, histArray[isam][2]->GetBinError(3)/2.);
      histArray[isam][2]->SetBinContent(3, histArray[isam][2]->GetBinContent(3)/2.);
      histArray[isam][2]->SetBinError(4, histArray[isam][2]->GetBinError(4)/2.);
      histArray[isam][2]->SetBinContent(4, histArray[isam][2]->GetBinContent(4)/2.);
      histArray[isam][2]->SetBinError(5, histArray[isam][2]->GetBinError(5)/4.);
      histArray[isam][2]->SetBinContent(5, histArray[isam][2]->GetBinContent(5)/4.);

       }

    if(plotType=="Jet2Pt_Fig19"){
      Double_t rebin_array[4] = {0.0, 100.0, 200.0, 600.0};
      histArray[isam][0] = (TH1F*)histArray[isam][0]->Rebin(3, "temp1", rebin_array);
      histArray[isam][0]->SetBinError(3, histArray[isam][0]->GetBinError(3)/4.);
      histArray[isam][0]->SetBinContent(3, histArray[isam][0]->GetBinContent(3)/4.);
      histArray[isam][1] = (TH1F*)histArray[isam][1]->Rebin(3, "temp2", rebin_array);
      histArray[isam][1]->SetBinError(3, histArray[isam][1]->GetBinError(3)/4.);
      histArray[isam][1]->SetBinContent(3, histArray[isam][1]->GetBinContent(3)/4.);
      histArray[isam][2] = (TH1F*)histArray[isam][2]->Rebin(3, "temp3", rebin_array);
      histArray[isam][2]->SetBinError(3, histArray[isam][2]->GetBinError(3)/4.);
      histArray[isam][2]->SetBinContent(3, histArray[isam][2]->GetBinContent(3)/4.);

     }

   if(plotType=="Lep1Pt_Fig19"){

      Double_t rebin_array[6] = {0, 50, 100, 150, 300, 500};
      histArray[isam][0] = (TH1F*)histArray[isam][0]->Rebin(5, "temp1", rebin_array);
      histArray[isam][0]->SetBinError(4, histArray[isam][0]->GetBinError(4)/3.0);
      histArray[isam][0]->SetBinContent(4, histArray[isam][0]->GetBinContent(4)/3.0);  
      histArray[isam][0]->SetBinError(5, histArray[isam][0]->GetBinError(5)/4.0);
      histArray[isam][0]->SetBinContent(5, histArray[isam][0]->GetBinContent(5)/4.0);
      histArray[isam][1] = (TH1F*)histArray[isam][1]->Rebin(5, "temp2", rebin_array);
      histArray[isam][1]->SetBinError(4, histArray[isam][1]->GetBinError(4)/3.0);
      histArray[isam][1]->SetBinContent(4, histArray[isam][1]->GetBinContent(4)/3.0);
      histArray[isam][1]->SetBinError(5, histArray[isam][1]->GetBinError(5)/4.0);
      histArray[isam][1]->SetBinContent(5, histArray[isam][1]->GetBinContent(5)/4.0);
      histArray[isam][2] = (TH1F*)histArray[isam][2]->Rebin(5, "temp3", rebin_array);
      histArray[isam][2]->SetBinError(4, histArray[isam][2]->GetBinError(4)/3.0);
      histArray[isam][2]->SetBinContent(4, histArray[isam][2]->GetBinContent(4)/3.0);
      histArray[isam][2]->SetBinError(5, histArray[isam][2]->GetBinError(5)/4.0);
      histArray[isam][2]->SetBinContent(5, histArray[isam][2]->GetBinContent(5)/4.0);



    }


   if(plotType=="Lep2Pt_Fig19"){//check this
    //100, 150, 260
 
    Double_t rebin_array[5] = {0.0, 50.0, 100, 150, 260};

    histArray[isam][0] = (TH1F*)histArray[isam][0]->Rebin(4, "temp1", rebin_array);
    histArray[isam][0]->SetBinError(4, histArray[isam][0]->GetBinError(4)*(5.0/11.0));
    histArray[isam][0]->SetBinContent(4, histArray[isam][0]->GetBinContent(4)*(5.0/11.0));
    histArray[isam][1] = (TH1F*)histArray[isam][1]->Rebin(4, "temp2", rebin_array);
    histArray[isam][1]->SetBinError(4, histArray[isam][1]->GetBinError(4)*(5.0/11.0));
    histArray[isam][1]->SetBinContent(4, histArray[isam][1]->GetBinContent(4)*(5.0/11.0));
    histArray[isam][2] = (TH1F*)histArray[isam][2]->Rebin(4, "temp3", rebin_array);
    histArray[isam][2]->SetBinError(4, histArray[isam][2]->GetBinError(4)*(5.0/11.0));
    histArray[isam][2]->SetBinContent(4, histArray[isam][2]->GetBinContent(4)*(5.0/11.0));



    }

   if(plotType=="MET_Fig20"){


     //180, 220, 280, 360
     Double_t rebin_array[13] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 220, 280, 360};
     histArray[isam][0] = (TH1F*)histArray[isam][0]->Rebin(12, "temp1", rebin_array);
     histArray[isam][0]->SetBinError(10, histArray[isam][0]->GetBinError(10)/2.);
     histArray[isam][0]->SetBinContent(10, histArray[isam][0]->GetBinContent(10)/2.);
     histArray[isam][0]->SetBinError(11, histArray[isam][0]->GetBinError(11)/3.);
     histArray[isam][0]->SetBinContent(11, histArray[isam][0]->GetBinContent(11)/3.);
     histArray[isam][0]->SetBinError(12, histArray[isam][0]->GetBinError(12)/4.);
     histArray[isam][0]->SetBinContent(12, histArray[isam][0]->GetBinContent(12)/4.);
     histArray[isam][1] = (TH1F*)histArray[isam][1]->Rebin(12, "temp2", rebin_array);
     histArray[isam][1]->SetBinError(10, histArray[isam][1]->GetBinError(10)/2.);
     histArray[isam][1]->SetBinContent(10, histArray[isam][1]->GetBinContent(10)/2.);
     histArray[isam][1]->SetBinError(11, histArray[isam][1]->GetBinError(11)/3.);
     histArray[isam][1]->SetBinContent(11, histArray[isam][1]->GetBinContent(11)/3.);
     histArray[isam][1]->SetBinError(12, histArray[isam][1]->GetBinError(12)/4.);
     histArray[isam][1]->SetBinContent(12, histArray[isam][1]->GetBinContent(12)/4.);
     histArray[isam][2] = (TH1F*)histArray[isam][2]->Rebin(12, "temp3", rebin_array);
     histArray[isam][2]->SetBinError(10, histArray[isam][2]->GetBinError(10)/2.);
     histArray[isam][2]->SetBinContent(10, histArray[isam][2]->GetBinContent(10)/2.);
     histArray[isam][2]->SetBinError(11, histArray[isam][2]->GetBinError(11)/3.);
     histArray[isam][2]->SetBinContent(11, histArray[isam][2]->GetBinContent(11)/3.);
     histArray[isam][2]->SetBinError(12, histArray[isam][2]->GetBinError(12)/4.);
     histArray[isam][2]->SetBinContent(12, histArray[isam][2]->GetBinContent(12)/4.);

     }

   if(plotType=="minMlb_Fig20"){

     //100, 250, 500
     Double_t rebin_array[5] = {0, 50, 100, 250, 500};    //rebinning minMlb for TRI
     histArray[isam][0] = (TH1F*)histArray[isam][0]->Rebin(4, "temp1", rebin_array);
     histArray[isam][0]->SetBinError(3, histArray[isam][0]->GetBinError(3)/3.0);
     histArray[isam][0]->SetBinContent(3, histArray[isam][0]->GetBinContent(3)/3.0);
     histArray[isam][0]->SetBinError(4, histArray[isam][0]->GetBinError(4)/5.0);
     histArray[isam][0]->SetBinContent(4, histArray[isam][0]->GetBinContent(4)/5.0);
     histArray[isam][1] = (TH1F*)histArray[isam][1]->Rebin(4, "temp2", rebin_array);
     histArray[isam][1]->SetBinError(3, histArray[isam][1]->GetBinError(3)/3.0);
     histArray[isam][1]->SetBinContent(3, histArray[isam][1]->GetBinContent(3)/3.0);
     histArray[isam][1]->SetBinError(4, histArray[isam][1]->GetBinError(4)/5.0);
     histArray[isam][1]->SetBinContent(4, histArray[isam][1]->GetBinContent(4)/5.0);
     histArray[isam][2] = (TH1F*)histArray[isam][2]->Rebin(4, "temp3", rebin_array);
     histArray[isam][2]->SetBinError(3, histArray[isam][2]->GetBinError(3)/3.0);
     histArray[isam][2]->SetBinContent(3, histArray[isam][2]->GetBinContent(3)/3.0);
     histArray[isam][2]->SetBinError(4, histArray[isam][2]->GetBinError(4)/5.0);
     histArray[isam][2]->SetBinContent(4, histArray[isam][2]->GetBinContent(4)/5.0);
 

    }
 

    if(plotType=="HT_Fig21"){
     //600, 800, 1000, 1200, 1600
     Double_t rebin_array[11] = {0, 100, 200, 300, 400, 500, 600, 800, 1000, 1200, 1600}; //HT
     histArray[isam][0] = (TH1F*)histArray[isam][0]->Rebin(10, "temp1", rebin_array);
     histArray[isam][0]->SetBinError(7, histArray[isam][0]->GetBinError(7)/2.);
     histArray[isam][0]->SetBinContent(7, histArray[isam][0]->GetBinContent(7)/2.);
     histArray[isam][0]->SetBinError(8, histArray[isam][0]->GetBinError(8)/2.);
     histArray[isam][0]->SetBinContent(8, histArray[isam][0]->GetBinContent(8)/2.);
     histArray[isam][0]->SetBinError(9, histArray[isam][0]->GetBinError(9)/2.0);
     histArray[isam][0]->SetBinContent(9, histArray[isam][0]->GetBinContent(9)/2.0);
     histArray[isam][0]->SetBinError(10, histArray[isam][0]->GetBinError(10)/4.);
     histArray[isam][0]->SetBinContent(10, histArray[isam][0]->GetBinContent(10)/4.);
     histArray[isam][1] = (TH1F*)histArray[isam][1]->Rebin(10, "temp2", rebin_array);
     histArray[isam][1]->SetBinError(7, histArray[isam][1]->GetBinError(7)/2.);
     histArray[isam][1]->SetBinContent(7, histArray[isam][1]->GetBinContent(7)/2.);
     histArray[isam][1]->SetBinError(8, histArray[isam][1]->GetBinError(8)/2.);
     histArray[isam][1]->SetBinContent(8, histArray[isam][1]->GetBinContent(8)/2.);
     histArray[isam][1]->SetBinError(9, histArray[isam][1]->GetBinError(9)/2.0);
     histArray[isam][1]->SetBinContent(9, histArray[isam][1]->GetBinContent(9)/2.0);
     histArray[isam][1]->SetBinError(10, histArray[isam][1]->GetBinError(10)/4.);
     histArray[isam][1]->SetBinContent(10, histArray[isam][1]->GetBinContent(10)/4.);
     histArray[isam][2] = (TH1F*)histArray[isam][2]->Rebin(10, "temp3", rebin_array);
     histArray[isam][2]->SetBinError(7, histArray[isam][2]->GetBinError(7)/2.);
     histArray[isam][2]->SetBinContent(7, histArray[isam][2]->GetBinContent(7)/2.);
     histArray[isam][2]->SetBinError(8, histArray[isam][2]->GetBinError(8)/2.);
     histArray[isam][2]->SetBinContent(8, histArray[isam][2]->GetBinContent(8)/2.);
     histArray[isam][2]->SetBinError(9, histArray[isam][2]->GetBinError(9)/2.0);
     histArray[isam][2]->SetBinContent(9, histArray[isam][2]->GetBinContent(9)/2.0);
     histArray[isam][2]->SetBinError(10, histArray[isam][2]->GetBinError(10)/4.);
     histArray[isam][2]->SetBinContent(10, histArray[isam][2]->GetBinContent(10)/4.);

     }

    if(plotType=="sumPtL_Fig21"){

     //250, 300, 400, 500
     Double_t rebin_array[5] = {200, 250, 300, 400, 500};    //rebinning for sum lepton pt 
     histArray[isam][0] = (TH1F*)histArray[isam][0]->Rebin(4, "temp1", rebin_array);
     histArray[isam][0]->SetBinError(3, histArray[isam][0]->GetBinError(3)/2.);
     histArray[isam][0]->SetBinContent(3, histArray[isam][0]->GetBinContent(3)/2.);
     histArray[isam][0]->SetBinError(4, histArray[isam][0]->GetBinError(4)/2.);
     histArray[isam][0]->SetBinContent(4, histArray[isam][0]->GetBinContent(4)/2.);
     
     histArray[isam][1] = (TH1F*)histArray[isam][1]->Rebin(4, "temp2", rebin_array);
     histArray[isam][1]->SetBinError(3, histArray[isam][1]->GetBinError(3)/2.);
     histArray[isam][1]->SetBinContent(3, histArray[isam][1]->GetBinContent(3)/2.);
     histArray[isam][1]->SetBinError(4, histArray[isam][1]->GetBinError(4)/2.);
     histArray[isam][1]->SetBinContent(4, histArray[isam][1]->GetBinContent(4)/2.);

     histArray[isam][2] = (TH1F*)histArray[isam][2]->Rebin(4, "temp3", rebin_array);
     histArray[isam][2]->SetBinError(3, histArray[isam][2]->GetBinError(3)/2.);
     histArray[isam][2]->SetBinContent(3, histArray[isam][2]->GetBinContent(3)/2.);
     histArray[isam][2]->SetBinError(4, histArray[isam][2]->GetBinError(4)/2.);
     histArray[isam][2]->SetBinContent(4, histArray[isam][2]->GetBinContent(4)/2.);

}
 
   if(plotType=="ST_Fig21" || plotType=="leptonJetsSum_Fig21"){

    // 800, 1000, 1500, 2000
    Double_t rebin_array[16] = {200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 1000, 1500, 2000};
    histArray[isam][0] = (TH1F*)histArray[isam][0]->Rebin(15, "temp1", rebin_array);
    histArray[isam][0]->SetBinError(13, histArray[isam][0]->GetBinError(13)/4.);
    histArray[isam][0]->SetBinContent(13, histArray[isam][0]->GetBinContent(13)/4.);
    histArray[isam][0]->SetBinError(14, histArray[isam][0]->GetBinError(14)/10.);
    histArray[isam][0]->SetBinContent(14, histArray[isam][0]->GetBinContent(14)/10.);
    histArray[isam][0]->SetBinError(15, histArray[isam][0]->GetBinError(15)/10.);
    histArray[isam][0]->SetBinContent(15, histArray[isam][0]->GetBinContent(15)/10.);

    histArray[isam][1] = (TH1F*)histArray[isam][1]->Rebin(15, "temp2", rebin_array);
    histArray[isam][1]->SetBinError(13, histArray[isam][1]->GetBinError(13)/4.);
    histArray[isam][1]->SetBinContent(13, histArray[isam][1]->GetBinContent(13)/4.);
    histArray[isam][1]->SetBinError(14, histArray[isam][1]->GetBinError(14)/10.);
    histArray[isam][1]->SetBinContent(14, histArray[isam][1]->GetBinContent(14)/10.);
    histArray[isam][1]->SetBinError(15, histArray[isam][1]->GetBinError(15)/10.);
    histArray[isam][1]->SetBinContent(15, histArray[isam][1]->GetBinContent(15)/10.);

    histArray[isam][2] = (TH1F*)histArray[isam][2]->Rebin(15, "temp3", rebin_array);
    histArray[isam][2]->SetBinError(13, histArray[isam][2]->GetBinError(13)/4.);
    histArray[isam][2]->SetBinContent(13, histArray[isam][2]->GetBinContent(13)/4.);
    histArray[isam][2]->SetBinError(14, histArray[isam][2]->GetBinError(14)/10.);
    histArray[isam][2]->SetBinContent(14, histArray[isam][2]->GetBinContent(14)/10.);
    histArray[isam][2]->SetBinError(15, histArray[isam][2]->GetBinError(15)/10.);
    histArray[isam][2]->SetBinContent(15, histArray[isam][2]->GetBinContent(15)/10.);

   }





      //histArray[isam][0]->Rebin(rebin);
      //histArray[isam][1]->Rebin(rebin);
      //histArray[isam][2]->Rebin(rebin);
    


     //Double_t rebin_array[24] = {0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 1000, 1100, 1200, 1300, 1500}; //rebinning ST
//ST binning for SS
/*
      Double_t rebin_array[17] = {200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 900, 1000, 1200, 1500};
      histArray[isam][0] = (TH1F*)histArray[isam][0]->Rebin(16, "temp1", rebin_array);  
      histArray[isam][0]->SetBinError(13, histArray[isam][0]->GetBinError(13)/2.);
      histArray[isam][0]->SetBinContent(13, histArray[isam][0]->GetBinContent(13)/2.);
      histArray[isam][0]->SetBinError(14, histArray[isam][0]->GetBinError(14)/2.);
      histArray[isam][0]->SetBinContent(14, histArray[isam][0]->GetBinContent(14)/2.);
      histArray[isam][0]->SetBinError(15, histArray[isam][0]->GetBinError(15)/4.);
      histArray[isam][0]->SetBinContent(15, histArray[isam][0]->GetBinContent(15)/4.);
      histArray[isam][0]->SetBinError(16, histArray[isam][0]->GetBinError(16)/6.);
      histArray[isam][0]->SetBinContent(16, histArray[isam][0]->GetBinContent(16)/6.);

      histArray[isam][1] = (TH1F*)histArray[isam][1]->Rebin(16, "temp1", rebin_array);
      histArray[isam][1]->SetBinError(13, histArray[isam][1]->GetBinError(13)/2.);
      histArray[isam][1]->SetBinContent(13, histArray[isam][1]->GetBinContent(13)/2.);
      histArray[isam][1]->SetBinError(14, histArray[isam][1]->GetBinError(14)/2.);
      histArray[isam][1]->SetBinContent(14, histArray[isam][1]->GetBinContent(14)/2.);
      histArray[isam][1]->SetBinError(15, histArray[isam][1]->GetBinError(15)/4.);
      histArray[isam][1]->SetBinContent(15, histArray[isam][1]->GetBinContent(15)/4.);
      histArray[isam][1]->SetBinError(16, histArray[isam][1]->GetBinError(16)/6.);
      histArray[isam][1]->SetBinContent(16, histArray[isam][1]->GetBinContent(16)/6.);

      histArray[isam][2] = (TH1F*)histArray[isam][2]->Rebin(16, "temp1", rebin_array);
      histArray[isam][2]->SetBinError(13, histArray[isam][2]->GetBinError(13)/2.);
      histArray[isam][2]->SetBinContent(13, histArray[isam][2]->GetBinContent(13)/2.);
      histArray[isam][2]->SetBinError(14, histArray[isam][2]->GetBinError(14)/2.);
      histArray[isam][2]->SetBinContent(14, histArray[isam][2]->GetBinContent(14)/2.);
      histArray[isam][2]->SetBinError(15, histArray[isam][2]->GetBinError(15)/4.);
      histArray[isam][2]->SetBinContent(15, histArray[isam][2]->GetBinContent(15)/4.);
      histArray[isam][2]->SetBinError(16, histArray[isam][2]->GetBinError(16)/6.);
      histArray[isam][2]->SetBinContent(16, histArray[isam][2]->GetBinContent(16)/6.);

*/
  
/*
     Double_t rebin_array[14] = {0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1200, 1400, 1600}; //HT
     
     histArray[isam][0] = (TH1F*)histArray[isam][0]->Rebin(13, "temp1", rebin_array);
     histArray[isam][0]->SetBinError(11, histArray[isam][0]->GetBinError(11)/2.);
     histArray[isam][0]->SetBinContent(11, histArray[isam][0]->GetBinContent(11)/2.);  
     histArray[isam][0]->SetBinError(12, histArray[isam][0]->GetBinError(12)/2.);
     histArray[isam][0]->SetBinContent(12, histArray[isam][0]->GetBinContent(12)/2.);
     histArray[isam][0]->SetBinError(13, histArray[isam][0]->GetBinError(13)/2.);
     histArray[isam][0]->SetBinContent(13, histArray[isam][0]->GetBinContent(13)/2.);
     histArray[isam][1] = (TH1F*)histArray[isam][1]->Rebin(13, "temp2", rebin_array);
     histArray[isam][1]->SetBinError(11, histArray[isam][1]->GetBinError(11)/2.);
     histArray[isam][1]->SetBinContent(11, histArray[isam][1]->GetBinContent(11)/2.);
     histArray[isam][1]->SetBinError(12, histArray[isam][1]->GetBinError(12)/2.);
     histArray[isam][1]->SetBinContent(12, histArray[isam][1]->GetBinContent(12)/2.);
     histArray[isam][1]->SetBinError(13, histArray[isam][1]->GetBinError(13)/2.);
     histArray[isam][1]->SetBinContent(13, histArray[isam][1]->GetBinContent(13)/2.);
     histArray[isam][2] = (TH1F*)histArray[isam][2]->Rebin(13, "temp3", rebin_array);
     histArray[isam][2]->SetBinError(11, histArray[isam][2]->GetBinError(11)/2.);
     histArray[isam][2]->SetBinContent(11, histArray[isam][2]->GetBinContent(11)/2.);
     histArray[isam][2]->SetBinError(12, histArray[isam][2]->GetBinError(12)/2.);
     histArray[isam][2]->SetBinContent(12, histArray[isam][2]->GetBinContent(12)/2.);
     histArray[isam][2]->SetBinError(13, histArray[isam][2]->GetBinError(13)/2.);
     histArray[isam][2]->SetBinContent(13, histArray[isam][2]->GetBinContent(13)/2.); 

     Double_t rebin_array[15] = {0, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 400, 500};    //rebinning for sum lepton pt 
//and minMlb

     histArray[isam][0] = (TH1F*)histArray[isam][0]->Rebin(14, "temp1", rebin_array);
     histArray[isam][0]->SetBinError(12, histArray[isam][0]->GetBinError(12)/4.);
     histArray[isam][0]->SetBinContent(12, histArray[isam][0]->GetBinContent(12)/4.);  
     histArray[isam][0]->SetBinError(13, histArray[isam][0]->GetBinError(13)/4.);
     histArray[isam][0]->SetBinContent(13, histArray[isam][0]->GetBinContent(13)/4.);
     histArray[isam][0]->SetBinError(14, histArray[isam][0]->GetBinError(14)/4.);
     histArray[isam][0]->SetBinContent(14, histArray[isam][0]->GetBinContent(14)/4.);
     histArray[isam][1] = (TH1F*)histArray[isam][1]->Rebin(14, "temp2", rebin_array);
     histArray[isam][1]->SetBinError(12, histArray[isam][1]->GetBinError(12)/4.); 
     histArray[isam][1]->SetBinContent(12, histArray[isam][1]->GetBinContent(12)/4.);
     histArray[isam][1]->SetBinError(13, histArray[isam][1]->GetBinError(13)/4.);
     histArray[isam][1]->SetBinContent(13, histArray[isam][1]->GetBinContent(13)/4.);
     histArray[isam][1]->SetBinError(14, histArray[isam][1]->GetBinError(14)/4.);
     histArray[isam][1]->SetBinContent(14, histArray[isam][1]->GetBinContent(14)/4.);
     histArray[isam][2] = (TH1F*)histArray[isam][2]->Rebin(14, "temp3", rebin_array);
     histArray[isam][2]->SetBinError(12, histArray[isam][2]->GetBinError(12)/4.);
     histArray[isam][2]->SetBinContent(12, histArray[isam][2]->GetBinContent(12)/4.);
     histArray[isam][2]->SetBinError(13, histArray[isam][2]->GetBinError(13)/4.);
     histArray[isam][2]->SetBinContent(13, histArray[isam][2]->GetBinContent(13)/4.);
     histArray[isam][2]->SetBinError(14, histArray[isam][2]->GetBinError(14)/4.);
     histArray[isam][2]->SetBinContent(14, histArray[isam][2]->GetBinContent(14)/4.); 


     Double_t rebin_array[29] = {0, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375, 400, 425, 450, 475, 500, 550, 600, 650, 700, 750, 800, 850, 900};    //rebinning for sum lepton pt in OS 
//and minMlb

     histArray[isam][0] = (TH1F*)histArray[isam][0]->Rebin(28, "temp1", rebin_array);
     histArray[isam][0]->SetBinError(20, histArray[isam][0]->GetBinError(20)/2.);
     histArray[isam][0]->SetBinContent(20, histArray[isam][0]->GetBinContent(20)/2.);  
     histArray[isam][0]->SetBinError(21, histArray[isam][0]->GetBinError(21)/2.);
     histArray[isam][0]->SetBinContent(21, histArray[isam][0]->GetBinContent(21)/2.);
     histArray[isam][0]->SetBinError(22, histArray[isam][0]->GetBinError(22)/2.);
     histArray[isam][0]->SetBinContent(22, histArray[isam][0]->GetBinContent(22)/2.);
     histArray[isam][0]->SetBinError(23, histArray[isam][0]->GetBinError(23)/2.);
     histArray[isam][0]->SetBinContent(23, histArray[isam][0]->GetBinContent(23)/2.);
     histArray[isam][0]->SetBinError(24, histArray[isam][0]->GetBinError(24)/2.);
     histArray[isam][0]->SetBinContent(24, histArray[isam][0]->GetBinContent(24)/2.);
     histArray[isam][0]->SetBinError(25, histArray[isam][0]->GetBinError(25)/2.);
     histArray[isam][0]->SetBinContent(25, histArray[isam][0]->GetBinContent(25)/2.);
     histArray[isam][0]->SetBinError(26, histArray[isam][0]->GetBinError(26)/2.);
     histArray[isam][0]->SetBinContent(26, histArray[isam][0]->GetBinContent(26)/2.);
     histArray[isam][0]->SetBinError(27, histArray[isam][0]->GetBinError(27)/2.);
     histArray[isam][0]->SetBinContent(27, histArray[isam][0]->GetBinContent(27)/2.);
     histArray[isam][0]->SetBinError(28, histArray[isam][0]->GetBinError(28)/2.);
     histArray[isam][0]->SetBinContent(28, histArray[isam][0]->GetBinContent(28)/2.);
     histArray[isam][1] = (TH1F*)histArray[isam][1]->Rebin(28, "temp1", rebin_array);
     histArray[isam][1]->SetBinError(20, histArray[isam][1]->GetBinError(20)/2.);
     histArray[isam][1]->SetBinContent(20, histArray[isam][1]->GetBinContent(20)/2.);
     histArray[isam][1]->SetBinError(21, histArray[isam][1]->GetBinError(21)/2.);
     histArray[isam][1]->SetBinContent(21, histArray[isam][1]->GetBinContent(21)/2.);
     histArray[isam][1]->SetBinError(22, histArray[isam][1]->GetBinError(22)/2.);
     histArray[isam][1]->SetBinContent(22, histArray[isam][1]->GetBinContent(22)/2.);
     histArray[isam][1]->SetBinError(23, histArray[isam][1]->GetBinError(23)/2.);
     histArray[isam][1]->SetBinContent(23, histArray[isam][1]->GetBinContent(23)/2.);
     histArray[isam][1]->SetBinError(24, histArray[isam][1]->GetBinError(24)/2.);
     histArray[isam][1]->SetBinContent(24, histArray[isam][1]->GetBinContent(24)/2.);
     histArray[isam][1]->SetBinError(25, histArray[isam][1]->GetBinError(25)/2.);
     histArray[isam][1]->SetBinContent(25, histArray[isam][1]->GetBinContent(25)/2.);
     histArray[isam][1]->SetBinError(26, histArray[isam][1]->GetBinError(26)/2.);
     histArray[isam][1]->SetBinContent(26, histArray[isam][1]->GetBinContent(26)/2.);
     histArray[isam][1]->SetBinError(27, histArray[isam][1]->GetBinError(27)/2.);
     histArray[isam][1]->SetBinContent(27, histArray[isam][1]->GetBinContent(27)/2.);
     histArray[isam][1]->SetBinError(28, histArray[isam][1]->GetBinError(28)/2.);
     histArray[isam][1]->SetBinContent(28, histArray[isam][1]->GetBinContent(28)/2.);
     histArray[isam][2] = (TH1F*)histArray[isam][2]->Rebin(28, "temp1", rebin_array);
     histArray[isam][2]->SetBinError(20, histArray[isam][2]->GetBinError(20)/2.);
     histArray[isam][2]->SetBinContent(20, histArray[isam][2]->GetBinContent(20)/2.);
     histArray[isam][2]->SetBinError(21, histArray[isam][2]->GetBinError(21)/2.);
     histArray[isam][2]->SetBinContent(21, histArray[isam][2]->GetBinContent(21)/2.);
     histArray[isam][2]->SetBinError(22, histArray[isam][2]->GetBinError(22)/2.);
     histArray[isam][2]->SetBinContent(22, histArray[isam][2]->GetBinContent(22)/2.);
     histArray[isam][2]->SetBinError(23, histArray[isam][2]->GetBinError(23)/2.);
     histArray[isam][2]->SetBinContent(23, histArray[isam][2]->GetBinContent(23)/2.);
     histArray[isam][2]->SetBinError(24, histArray[isam][2]->GetBinError(24)/2.);
     histArray[isam][2]->SetBinContent(24, histArray[isam][2]->GetBinContent(24)/2.);
     histArray[isam][2]->SetBinError(25, histArray[isam][2]->GetBinError(25)/2.);
     histArray[isam][2]->SetBinContent(25, histArray[isam][2]->GetBinContent(25)/2.);
     histArray[isam][2]->SetBinError(26, histArray[isam][2]->GetBinError(26)/2.);
     histArray[isam][2]->SetBinContent(26, histArray[isam][2]->GetBinContent(26)/2.);
     histArray[isam][2]->SetBinError(27, histArray[isam][2]->GetBinError(27)/2.);
     histArray[isam][2]->SetBinContent(27, histArray[isam][2]->GetBinContent(27)/2.);
     histArray[isam][2]->SetBinError(28, histArray[isam][2]->GetBinError(28)/2.);
     histArray[isam][2]->SetBinContent(28, histArray[isam][2]->GetBinContent(28)/2.);
*/


      }
    if (minX >-990){
      histArray[isam][0]->SetAxisRange(minX,maxX);
      histArray[isam][1]->SetAxisRange(minX,maxX);
      histArray[isam][2]->SetAxisRange(minX,maxX);
    }

    if ( (samples[isam]!=Data) && (samples[isam]!=FakeRate) && (samples[isam]!=ChargeMisID)){
//       histArray[isam][0]->Scale(19.62/16.52);
//       histArray[isam][1]->Scale(19.62/16.52);
//       histArray[isam][2]->Scale(19.62/16.52);
    }

    histArray[isam][3] = (TH1F*)histArray[isam][0]->Clone("histArray_"+sampleName);
    histArray[isam][3]->Add(histArray[isam][1]);
    histArray[isam][3]->Add(histArray[isam][2]);
  }

  for (int iChan = 0;iChan<4;++iChan) {
    switch( iChan ) {
      case 0: inSuff = "MuMu"; break;
      case 1: inSuff = "ElMu"; break;
      case 2: inSuff = "ElEl"; break;
      case 3: inSuff = "All" ; break;
    }
    vector <TH1F*> histos;
    for (unsigned int isam = 0; isam < NSAMPLES; isam++) {
      histos.push_back(histArray[isam][iChan]);
    }

    //Format data histogram
    histos.at(indices[Data])->SetMarkerStyle(20);
    for (int ibin = 0; ibin < histos.at(indices[Data])->GetNbinsX() + 1; ibin++){
      if (histos.at(indices[Data])->GetBinContent(ibin) == 0) histos.at(indices[Data])->SetBinContent(ibin, -10);
    }

    THStack* mystack = new THStack("mystack","mystack");
    TH1F * histo1D_mc = (TH1F *) histos.at(Data)->Clone("MC");
    histo1D_mc->Reset();
    TH1F * histo1D_multiBoson = (TH1F *) histos.at(Data)->Clone("MB");
    histo1D_multiBoson->Reset();
    TH1F * histo1D_ttBoson = (TH1F *) histos.at(Data)->Clone("ttv");
    histo1D_ttBoson->Reset();
    TH1F * histo1D_signal = (TH1F *) histos.at(Data)->Clone("signal");
    histo1D_signal->Reset();
    //histo1D_signal->SetLineColor(kBlack);
    //histo1D_signal->SetLineStyle(2);
    //histo1D_signal->SetLineWidth(4);

    histo1D_signal->SetMarkerStyle(22);
    histo1D_signal->SetMarkerSize(1.5);
    histo1D_signal->SetMarkerColor(kBlue);
    TH1F * ratio = (TH1F *) histos.at(Data)->Clone("ratio");

// Legend:

  //TLegend * leg = new TLegend(0.65, 0.6, 0.90, 0.94, "","brNDC");
  TLegend * leg = new TLegend(0.65, 0.70, 0.90, 0.94, "","brNDC");
  //TLegend * leg = new TLegend(0.65, 0.74, 0.90, 0.94, "","brNDC");
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  leg->SetLineColor(1);
  leg->SetLineStyle(0);
  leg->SetLineWidth(1);
  leg->SetFillColor(10);
  leg->SetFillStyle(0);

 /* if(plotType=="ST_Fig5" || plotType=="leptonJetsSum_Fig5"){
     TLegend * leg = new TLegend(0.65, 0.82, 0.90, 0.98, "","brNDC");
     leg->SetBorderSize(0);
     leg->SetTextSize(0.03);
     leg->SetLineColor(1);
     leg->SetLineStyle(0);
     leg->SetLineWidth(1);
     leg->SetFillColor(10);
     leg->SetFillStyle(0);
  
   }*/
  leg->AddEntry(histos.at(Data), "Data", "lp");
  leg->AddEntry(histo1D_signal, "t'800 All Modes", "lp");
    //h_lumiBand->Sumw2();

  cout << "Yield : "<< inSuff<<endl;
  cout << "Yield: Data  "<<getMyEntries(histos.at(Data))<<endl;
  float BW=0.5;
  float TH=0.25;
  float TZ=0.25;

  float singeTopYield = 0.;
  float dyYield = 0.;
  float totalBackground = 0.;
  float minimum=999999;
  double signalMax = 0.;
  bool mb=false;
  for (unsigned int isam = 0; isam < NSAMPLES; isam++){
    if (samples[isam]!=Data) {
      histos.at(isam)->SetFillColor(color[samples[isam]]);
      if (samples[isam]<Tprime400_bWbW) {

        minimum=getMinEntries(histos.at(isam),minimum);
        totalBackground+=histos.at(isam)->Integral();

	if (!((samples[isam]>=T_tW)&& (samples[isam]<=Tbar_s))) {
	  if ((samples[isam]!= DY1050) && (samples[isam]!= ZJets))
            cout << "Yield: "<< std::fixed<<std::setprecision(1)<<legendNames[samples[isam]] << " "<< histos.at(isam)->Integral()<<" "<<
	    endl;
	    else {dyYield+=histos.at(isam)->Integral();
	      histos.at(isam)->Scale(dySF[iChan]);
	    }
	} else singeTopYield+=histos.at(isam)->Integral();


	if (samples[isam]==T_tW)
	  leg->AddEntry(histos.at(isam), "Single Top", "f");
	if (samples[isam]==ZJets)
	  leg->AddEntry(histos.at(isam), "Drell Yan", "f");
// 	if (samples[isam]==WWp)
// 	  leg->AddEntry(histos.at(isam), "W^{#pm}W^{#pm}", "f");

	if ((samples[isam]<=WJets) || 
	  ((samples[isam]>=T_tW)&& (samples[isam]<=Tbar_s)))
	  {
            cout << "Adding: "<< std::fixed<<std::setprecision(1)<<legendNames[samples[isam]] << " "<< histos.at(isam)->Integral()<<" "<<
	    histos.at(isam)->GetBinContent(1)<<
	    endl;
	  	mystack->Add(histos.at(isam));
		 if (!((samples[isam]>=T_tW)&& (samples[isam]<=Tbar_s)) &&
		 (samples[isam]!= DY1050) && (samples[isam]!= ZJets) &&
		 (samples[isam]!= WWp) && (samples[isam]!= WWm) && !(samples[isam]==WJets))
		   leg->AddEntry(histos.at(isam), legendNames[samples[isam]], "f");
		}
	else if ((samples[isam]==WW) || (samples[isam]==WZ) || (samples[isam]==ZZ) ||
	(samples[isam]==WWm) || (samples[isam]==WWp) || (samples[isam]==WWW) ||
	(samples[isam]>=ZZZNoGs)){
	  if (!mb) leg->AddEntry(histo1D_multiBoson, "Multi-bosons", "f");
	  mb=true;
	  histo1D_multiBoson->Add(histos.at(isam));
	} else if ((samples[isam]==TTW) || (samples[isam]==TTZ) || (samples[isam]==TTWW)){
	if (histo1D_ttBoson->Integral()==0.) leg->AddEntry(histo1D_ttBoson, "t#bar{t}+bosons", "f");
		histo1D_ttBoson->Add(histos.at(isam));
    
        } else {
	cout << " Unknown : " << legendNames[samples[isam]]<<endl;
	}
	histo1D_mc->Add(histos.at(isam));

      } else {
      // This is for the data:
	//histos.at(isam)->Scale(500.0);
	histos.at(isam)->Scale(1.0);
        histos.at(isam)->SetMarkerStyle(22);
	histos.at(isam)->SetMarkerSize(1.5);
	histos.at(isam)->SetMarkerColor(color[samples[isam]]);
	histos.at(isam)->SetLineStyle(2);
	histos.at(isam)->SetLineWidth(4);
	histos.at(isam)->SetLineColor(color[samples[isam]]);
	histos.at(isam)->SetFillColor(0);
	float br=0.;
	if (allNames[samples[isam]].Contains("BWBW")) br = BW*BW;
	if (allNames[samples[isam]].Contains("BWTH")) br = 2*BW*TH;
	if (allNames[samples[isam]].Contains("BWTZ")) br = 2*BW*TZ;
	if (allNames[samples[isam]].Contains("THTH")) br = TH*TH  ;
	if (allNames[samples[isam]].Contains("THTZ")) br = 2*TH*TZ;
	if (allNames[samples[isam]].Contains("TZTZ")) br = TZ*TZ  ;

	histo1D_signal->Add(histos.at(isam),br);
	cout << histo1D_signal->Integral()<<" "<< histos.at(isam)->Integral()<<endl;
	signalMax = max(signalMax,histos.at(isam)->GetMaximum());
        leg->AddEntry(histos.at(isam), legendNames[samples[isam]], "l");
            cout << "Yield: "<< legendNames[samples[isam]] << " "<< getMyEntries(histos.at(isam))/500<<
	    endl;


      }
    }
  }

  cout << histo1D_multiBoson->Integral() <<" "<< histo1D_multiBoson->GetBinContent(1)<<
  endl;
  if (histo1D_multiBoson->Integral()>0.) {
    histo1D_multiBoson->SetFillColor(color[WW]);
    mystack->Add(histo1D_multiBoson);
  }
  cout << "ttV int "<<histo1D_ttBoson->Integral()<<" "<< histo1D_ttBoson->GetBinContent(1)<<endl;
  if (histo1D_ttBoson->Integral()>0.)    {
    histo1D_ttBoson->SetFillColor(color[TTW]);
    mystack->Add(histo1D_ttBoson);
  }


  cout << "Yield: DY: "<< std::fixed<<std::setprecision(1)<<dyYield<<endl;
  cout << "Yield: Single top: "<< singeTopYield<<endl;
  cout << "Yield: total background: "<< std::fixed<<std::setprecision(1)<<totalBackground<<endl;
  cout << "Yield: Signal Tprime 500: "<< histo1D_signal->Integral()/500.<<endl;

    //Calculate systematics
    vector <TH1F*> selectedHistos;
    vector <double> selectedUnc;

    for (unsigned int isam = 0; isam < NSAMPLES; isam++){
      if (samples[isam]!=Data  && samples[isam]<Tprime400_bWbW) {
	selectedHistos.push_back(histos.at(isam));
	selectedUnc.push_back(uncertainty(samples[isam]));
	//cout <<  "Unc: "<< allNames[samples[isam]]<<" "<<uncertainty(samples[isam])<<endl;
      }
    }

    TH1F* h_lumiBand = (TH1F*) histo1D_mc->Clone("h_lumiBand");
    
    for (int ibin = 1; ibin < h_lumiBand->GetNbinsX()+1; ibin++){

      double uncStat = 0;
      double uncSyst = 0;
      double uncTot  = 0;

      for (unsigned int ih = 0; ih < selectedHistos.size(); ih++){
        //if (selectedHistos.at(ih)->GetBinError(ibin)>0) 
	uncStat += selectedHistos.at(ih)->GetBinError(ibin)*selectedHistos.at(ih)->GetBinError(ibin);
     	uncSyst += selectedHistos.at(ih)->GetBinContent(ibin)*selectedHistos.at(ih)->GetBinContent(ibin)*selectedUnc.at(ih)*selectedUnc.at(ih);
        
       }

      uncTot = sqrt(uncStat + uncSyst);
      h_lumiBand->SetBinError(ibin, uncTot);
      cout << "Bin "<< ibin << " "<< h_lumiBand->GetBinContent(ibin) << " sqrt(uncStat) =  "<< sqrt(uncStat) << " sqrt(uncSyst) "<<sqrt(uncSyst)<< " uncTot "<<uncTot<<endl;
    }

  if (mystack->GetMaximum() > histos.at(indices[Data])->GetMaximum() + sqrt(histos.at(indices[Data])->GetMaximum())){
    if (yLog == 0) histos.at(indices[Data])->SetMaximum((mystack->GetMaximum() + sqrt(histos.at(indices[Data])->GetMaximum())) * 1.05);
    else           histos.at(indices[Data])->SetMaximum((mystack->GetMaximum() + sqrt(histos.at(indices[Data])->GetMaximum())) * 1.40);
  }


  //Larger description of channel
  TString schan = "";
  if (inSuff == "ElEl") schan = " - ee";
  if (inSuff == "ElMu") schan = " - e#mu";
  if (inSuff == "MuMu") schan = " - #mu#mu";
  if (inSuff == "All")  schan = "";

  if (yLog == 0) histos.at(indices[Data])->SetMinimum(0);
  else           histos.at(indices[Data])->SetMinimum(getMinEntries(histo1D_mc)/5);

  double maximum = max(histos.at(indices[Data])->GetMaximum(),signalMax);
  maximum = max(maximum,histo1D_mc->GetMaximum());
  cout << "max: "<<maximum<<endl;
  if (yLog == 0) histos.at(indices[Data])->SetMaximum(maximum * 1.2);
  else           histos.at(indices[Data])->SetMaximum(maximum * 5);

  //histos.at(indices[Data])->GetXaxis()->SetTitle(xTitle+schan);
  histos.at(indices[Data])->GetYaxis()->SetTitle(yTitle);
  //histos.at(indices[Data])->GetYaxis()->SetTitle("Events");
  //histos.at(indices[Data])->GetYaxis()->SetTitle("Events/50 GeV"); //for sum plots
  //histos.at(indices[Data])->GetYaxis()->SetTitle("Events/100 GeV"); //for HT
  //histos.at(indices[Data])->GetYaxis()->SetTitle("Events/25 GeV"); //for sumlepton pt and minMlb
  histos.at(indices[Data])->GetYaxis()->SetTitleOffset(1.3);
  //histos.at(indices[Data])->GetXaxis()->SetTitleOffset(1.3);

  TCanvas* canv = new TCanvas("canv", "canv", 600, 800);
  canv->SetLeftMargin(0.1);
  canv->SetRightMargin(0.03);
  canv->SetTopMargin(0.05);
  canv->SetBottomMargin(0.3);
  canv->cd();

  histos.at(indices[Data])->SetMarkerSize(1.2);
  histos.at(indices[Data])->GetXaxis()->SetLabelSize(0);
  histos.at(indices[Data])->GetXaxis()->SetTitleSize(0);
  //histos.at(indices[Data])->GetYaxis()->SetRangeUser(0.001, 80000.0);
  //histos.at(indices[Data])->GetYaxis()->SetRangeUser(0.1, 300000.0);
  //histos.at(indices[Data])->GetYaxis()->SetRangeUser(0.01, 2000000.0);
  //histos.at(indices[Data])->GetYaxis()->SetRangeUser(1.0, 2000000.0);
  //histos.at(indices[Data])->GetYaxis()->SetRangeUser(1.0, 3000000.0);
  //histos.at(indices[Data])->GetYaxis()->SetRangeUser(0.001, 80000.0);
  histos.at(indices[Data])->Draw("P E x0");
  mystack->Draw("hist same");
//   mystack->GetHistogram()->Draw("same");
  histos.at(indices[Data])->Draw("P E same");

     for (unsigned int isam = 0; isam < NSAMPLES; isam++){
       if (samples[isam]>=Tprime400_bWbW) {
// // 	histos.at(isam)->Draw("P E same x0");
  	histos.at(isam)->Draw("hist same");
       }
     }

  //histo1D_signal->Draw("hist same x0");
  //histo1D_signal->Draw("PE same x0");

  histo1D_signal->Draw("hist p same");
  gStyle->SetHatchesLineWidth(4);
  gStyle->SetErrorX(0.5);

  h_lumiBand->SetFillStyle(3005);
  h_lumiBand->SetFillColor(1);
  h_lumiBand->SetMarkerStyle(1);
  h_lumiBand->Draw("same e2");

  leg->Draw();

  //Mandatory text boxes
  TLatex* text1 = new TLatex(1.570061,23.08044,"CMS Preliminary, 19.6 fb^{-1} at #sqrt{s} = 8 TeV");
  text1->SetNDC();
  
  text1->SetTextAlign(13);
  text1->SetX(0.12);
  text1->SetY(1.00);
  //text1->SetY(0.982);
  text1->SetTextFont(42);
  text1->SetTextSizePixels(20);
  text1->Draw();


  //Double_t AXIS[24] = {0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 1000, 1100, 1200, 1300, 1500};
  //Double_t AXIS[17] = {200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 900, 1000, 1200, 1500};
  //histos.at(indices[Data])->Draw("AXISsame");


  if (getMyEntries(histos.at(Data))>0. && (yLog == 1) ) gPad->SetLogy(yLog);

    TPad *pad = new TPad("pad", "pad", 0.0, 0.0, 1.0, 1.0);
    pad->SetTopMargin(0.70);
    pad->SetLeftMargin(0.1);
    pad->SetRightMargin(0.03);
    pad->SetFillColor(0);
    pad->SetFillStyle(0);
//     pad->SetTickx(1);
//     pad->SetTicky(1);
    pad->Draw();
    pad->cd(0);
    pad->SetGridy(true);
    ratio->Divide(h_lumiBand);
//    ratio->Divide(histo1D_mc);
    ratio->GetXaxis()->SetTitle(xTitle+schan);
    ratio->GetYaxis()->SetTitleOffset(1.3);
    ratio->GetYaxis()->SetTitle("Data/MC     ");
    ratio->GetYaxis()->SetNdivisions(505);
    //ratio->GetYaxis()->SetTitle("Data/MC");
    ratio->Draw("p e x0");
    ratio->SetMaximum( 2.5);
    ratio->SetMinimum( -0.5);

//    histos.at(Data)->Draw("same");


  TString outName = inName+inPref;
  if (inSuff != "All") outName += "_"+inSuff;


  canv->SaveAs(outName+".png");
//   canv->SaveAs(outName+".eps");
  canv->SaveAs(outName+".pdf");
  canv->SaveAs(outName+".C");
  }
}
// int main( int argc, const char* argv[] ){
//
// TString channel[4] = {"ElEl","ElMu","MuMu", "All"}
//
//   TString rootFile  = argv[1];
//   TString inDir  = argv[2];
//   TString inSuff = argv[3];
//   TString inPref = argv[4];
//   TString inName = argv[5];
//
//   //Drawing information
//   int yLog = atoi(argv[6]);
//   TString xTitle = argv[7];
//
//   vector<TString> legend(NHISTOS);
//   legend.at(Lep1Pt)    = TString("Leading lepton p_{T} [GeV/c]");
// legend.at(Lep2Pt)    = TString("Second lepton p_{T} [GeV/c]");
// legend.at(Lep3Pt)    = TString("Third lepton p_{T} [GeV/c]");
//
// legend.at(ElPt)	  = TString("Electron p_{T} [GeV/c]");
// legend.at(MuPt)	  = TString("Muon p_{T} [GeV/c]");
//
// legend.at(ElEta)	  = TString("");
// legend.at(MuEta)	  = TString("");
//
// legend.at(Jet1Pt)    = TString("Leading jet p_{T} [GeV/c]");
// legend.at(Jet2Pt)    = TString("Second jet p_{T} [GeV/c]");
//
// legend.at(nJets)	  = TString("Number of Jets"	    );
// legend.at(nBJets)    = TString("Number of b-tagged jets");
// legend.at(nElec)	  = TString("Number of Electrons");
// legend.at(nMuon)	  = TString("Number of Muons");
// legend.at(nLept)	  = TString("Number of leptons");
// legend.at(sumPtL)    = TString("Sum of lepton p_{T} [GeV/c]");
// legend.at(sumPtJ)    = TString("Sum of jet p_{T} [GeV/c]");
// legend.at(MET)	  = TString("MET [GeV]");
//
// legend.at(HT)	  = TString("HT [GeV]");
//
// legend.at(LepInvM)   = TString("Invariant Mass of Lepton Pair [GeV/c^{2}]");
// legend.at(LepJInvM)  = TString("");
// legend.at(dRlept)      = TString("#delta R (l1,l2)");
// legend.at(mindR)	    = TString("Min (#delta R (l,b-jet))");
// legend.at(dR)	    = TString("#delta R (l,b-jet)");
// legend.at(minMlb)      = TString("Min(m_{lb}})");
// legend.at(Mlb)	    = TString("m_{lb}}");
//
// };
//
//
//   for (int=0,i<4;++i)
//   {
//     for (unsigned int ih = 0; ih < NHISTOS; ih++){
//       drawSinglePlot(rootFile, inDir, channel[i] , "Lep1Pt"  1 legend[ih], 10);
//
// }
//
// }

double uncertainty (Samples_t sampleIndex )
{
  double mcUnc = sqrt (2 * 0.006 * 0.006 + //Trigger
		       2 * 0.03  * 0.03);  //Lepton Efficiency
  double systUnc;

      switch( sampleIndex )
      {
        case Data:	   systUnc = 0.  ; break;
	case ChargeMisID:  systUnc = 0.20; break;
	case FakeRate:	   systUnc = 0.50; break;

	case TTbar:	   systUnc = 0.08; break;
	case ZJets:	   systUnc = 0.30; break;
	case DY1050:	   systUnc = 0.30; break;
	case WJets:	   systUnc = 0.30; break;

	case T_tW:	   systUnc = 0.20;  break;
	case Tbar_tW:	   systUnc = 0.20;  break;
	case T_t:	   systUnc = 0.20;  break;
	case Tbar_t:	   systUnc = 0.20;  break;
	case T_s:	   systUnc = 0.20;  break;
	case Tbar_s:	   systUnc = 0.20;  break;

	case WZ:	   systUnc = 0.17; break;
	case ZZ:	   systUnc = 0.07; break;

	case WWm:	   systUnc = 0.50; break;
	case WWp:	   systUnc = 0.50; break;
	case WWW:	   systUnc = 0.50; break;

	case TTW:	   systUnc = 0.32; break;
	case TTWW:	   systUnc = 0.50; break;
	case TTZ:	   systUnc = 0.50; break;

        case Tprime400_bWbW: systUnc = mcUnc; break;
	case Tprime400_bWtH: systUnc = mcUnc; break;
	case Tprime400_bWtZ: systUnc = mcUnc; break;
	case Tprime400_tZtZ: systUnc = mcUnc; break;
	case Tprime400_tHtZ: systUnc = mcUnc; break;
        case Tprime400_tHtH: systUnc = mcUnc; break;

        case Tprime500_bWbW: systUnc = mcUnc; break;
	case Tprime500_bWtH: systUnc = mcUnc; break;
	case Tprime500_bWtZ: systUnc = mcUnc; break;
	case Tprime500_tZtZ: systUnc = mcUnc; break;
	case Tprime500_tHtZ: systUnc = mcUnc; break;
        case Tprime500_tHtH: systUnc = mcUnc; break;
        case Tprime600_bWbW: systUnc = mcUnc; break;
	case Tprime600_bWtH: systUnc = mcUnc; break;
	case Tprime600_bWtZ: systUnc = mcUnc; break;
	case Tprime600_tZtZ: systUnc = mcUnc; break;
	case Tprime600_tHtZ: systUnc = mcUnc; break;
        case Tprime600_tHtH: systUnc = mcUnc; break;
        case Tprime700_bWbW: systUnc = mcUnc; break;
	case Tprime700_bWtH: systUnc = mcUnc; break;
	case Tprime700_bWtZ: systUnc = mcUnc; break;
	case Tprime700_tZtZ: systUnc = mcUnc; break;
	case Tprime700_tHtZ: systUnc = mcUnc; break;
        case Tprime700_tHtH: systUnc = mcUnc; break;
        case Tprime800_bWbW: systUnc = mcUnc; break;
	case Tprime800_bWtH: systUnc = mcUnc; break;
	case Tprime800_bWtZ: systUnc = mcUnc; break;
	case Tprime800_tZtZ: systUnc = mcUnc; break;
	case Tprime800_tHtZ: systUnc = mcUnc; break;
        case Tprime800_tHtH: systUnc = mcUnc; break;
        case Tprime900_bWbW: systUnc = mcUnc; break;
	case Tprime900_bWtH: systUnc = mcUnc; break;
	case Tprime900_bWtZ: systUnc = mcUnc; break;
	case Tprime900_tZtZ: systUnc = mcUnc; break;
	case Tprime900_tHtZ: systUnc = mcUnc; break;
        case Tprime900_tHtH: systUnc = mcUnc; break;
        case Tprime1000_bWbW:  systUnc = mcUnc; break;
	case Tprime1000_bWtH:  systUnc = mcUnc; break;
	case Tprime1000_bWtZ:  systUnc = mcUnc; break;
	case Tprime1000_tZtZ:  systUnc = mcUnc; break;
	case Tprime1000_tHtZ:  systUnc = mcUnc; break;
        case Tprime1000_tHtH:  systUnc = mcUnc; break;
        case Tprime1100_bWbW:  systUnc = mcUnc; break;
	case Tprime1100_bWtH:  systUnc = mcUnc; break;
	case Tprime1100_bWtZ:  systUnc = mcUnc; break;
	case Tprime1100_tZtZ:  systUnc = mcUnc; break;
	case Tprime1100_tHtZ:  systUnc = mcUnc; break;
        case Tprime1100_tHtH:  systUnc = mcUnc; break;
        case Tprime1200_bWbW:  systUnc = mcUnc; break;
	case Tprime1200_bWtH:  systUnc = mcUnc; break;
	case Tprime1200_bWtZ:  systUnc = mcUnc; break;
	case Tprime1200_tZtZ:  systUnc = mcUnc; break;
	case Tprime1200_tHtZ:  systUnc = mcUnc; break;
        case Tprime1200_tHtH:  systUnc = mcUnc; break;
        case Tprime1300_bWbW:  systUnc = mcUnc; break;
	case Tprime1300_bWtH:  systUnc = mcUnc; break;
	case Tprime1300_bWtZ:  systUnc = mcUnc; break;
	case Tprime1300_tZtZ:  systUnc = mcUnc; break;
	case Tprime1300_tHtZ:  systUnc = mcUnc; break;
        case Tprime1300_tHtH:  systUnc = mcUnc; break;
        case Tprime1400_bWbW:  systUnc = mcUnc; break;
	case Tprime1400_bWtH:  systUnc = mcUnc; break;
	case Tprime1400_bWtZ:  systUnc = mcUnc; break;
	case Tprime1400_tZtZ:  systUnc = mcUnc; break;
	case Tprime1400_tHtZ:  systUnc = mcUnc; break;
        case Tprime1400_tHtH:  systUnc = mcUnc; break;
        case Tprime1500_bWbW:  systUnc = mcUnc; break;
	case Tprime1500_bWtH:  systUnc = mcUnc; break;
	case Tprime1500_bWtZ:  systUnc = mcUnc; break;
	case Tprime1500_tZtZ:  systUnc = mcUnc; break;
	case Tprime1500_tHtZ:  systUnc = mcUnc; break;
        case Tprime1500_tHtH:  systUnc = mcUnc; break;
	default:systUnc = mcUnc;
      }
  return systUnc;
}
#endif
