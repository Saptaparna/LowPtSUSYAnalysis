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
#include <algorithm>
#include <TGraphAsymmErrors.h>
#include <TVector2.h>
#include <TF1.h>


TLorentzVector fillTLorentzVector(double pT, double eta, double phi, double E)
{
  TLorentzVector object_p4;
  object_p4.SetPtEtaPhiE(pT, eta, phi, E);
  return object_p4;
}

typedef struct
{
  float pT;
  float eta;
  float phi;
  float energy;
  int charge;
  bool isTight;
  float isolation;
} LeptonInfo;

typedef struct
{
  float pT;
  float eta;
  float phi;
  float energy;
  int charge;
  bool isTight;
  float chIsolation;
  float nuIsolation;
  float phIsolation;
  bool  phIsoTight;
  bool  phIsoMedium;
  bool  phIsoLoose;
} PhotonInfo;

typedef struct
{
  float pT;
  float eta;
  float phi;
  float energy;
  int   PU_mva_loose;
  int   PU_mva_tight;
  int   PU_mva_medium;
  int   PU_cut_loose;
  int   PU_cut_tight;
  int   PU_cut_medium;
} JetInfo;


typedef struct
{
  float Px;
  float Py;
  float Pz;
  float E;
} TriggerInfo;

typedef struct
{
  int matched;
  float pT;
  float eta;
  float phi;
  float energy;
} MatchedLeptonInfo;

double photonSF(double phPt, double phEta)
{
  double sf = 1.0;
  if(phPt > 15.00 and phPt < 20.00)
    {
    if(fabs(phEta) > 0.0 and fabs(phEta) < 0.8) sf = 0.9496;
    else if(fabs(phEta) > 0.8 and fabs(phEta) < 1.4442) sf = 0.9803;
    else if(fabs(phEta) > 1.566 and fabs(phEta) < 2.0) sf = 1.0005;
    else if(fabs(phEta) > 2.0 and fabs(phEta) < 2.5) sf = 1.0171;
    }
  else if(phPt > 20.00 and phPt < 30.00)
    {
    if(fabs(phEta) > 0.0 and fabs(phEta) < 0.8) sf = 0.9672;
    else if(fabs(phEta) > 0.8 and fabs(phEta) < 1.4442) sf = 0.9724;
    else if(fabs(phEta) > 1.566 and fabs(phEta) < 2.0) sf = 0.9867;
    else if(fabs(phEta) > 2.0 and fabs(phEta) < 2.5) sf = 1.0130;
    }
  else if(phPt > 30.00 and phPt < 40.00)
    {
    if(fabs(phEta) > 0.0 and fabs(phEta) < 0.8) sf = 0.9711;
    else if(fabs(phEta) > 0.8 and fabs(phEta) < 1.4442) sf = 0.9688;
    else if(fabs(phEta) > 1.566 and fabs(phEta) < 2.0) sf = 0.9971;
    else if(fabs(phEta) > 2.0 and fabs(phEta) < 2.5) sf = 1.0143;
    }
  else if(phPt > 40.00 and phPt < 50.00)
    {
    if(fabs(phEta) > 0.0 and fabs(phEta) < 0.8) sf = 0.9766;
    else if(fabs(phEta) > 0.8 and fabs(phEta) < 1.4442) sf = 0.9805;
    else if(fabs(phEta) > 1.566 and fabs(phEta) < 2.0) sf = 0.9996;
    else if(fabs(phEta) > 2.0 and fabs(phEta) < 2.5) sf = 1.0129;
    }
  else if(phPt > 50.00)
    {
    if(fabs(phEta) > 0.0 and fabs(phEta) < 0.8) sf = 0.9815;
    else if(fabs(phEta) > 0.8 and fabs(phEta) < 1.4442) sf = 0.9837;
    else if(fabs(phEta) > 1.566 and fabs(phEta) < 2.0) sf = 1.0034;
    else if(fabs(phEta) > 2.0 and fabs(phEta) < 2.5) sf = 1.0128;
    }
  return sf;
}

double electronSF(double elecPt, double elecEta)
{
  double sf = 1.0;
  if(elecPt > 10.00 and elecPt < 15.00)
    {
    if(fabs(elecEta) > 0.0 and fabs(elecEta) < 0.8) sf = 0.838;
    else if(fabs(elecEta) > 0.8 and fabs(elecEta) < 1.44) sf = 0.861;
    else if(fabs(elecEta) > 1.44 and fabs(elecEta) < 1.56) sf = 1.021;
    else if(fabs(elecEta) > 1.56 and fabs(elecEta) < 2.00) sf = 0.951;
    else if(fabs(elecEta) > 2.00 and fabs(elecEta) < 2.50) sf = 1.055;
    }
  else if(elecPt > 15.00 and elecPt < 20.00)
    {
    if(fabs(elecEta) > 0.0 and fabs(elecEta) < 0.8) sf = 0.942;
    else if(fabs(elecEta) > 0.8 and fabs(elecEta) < 1.44) sf = 0.925;
    else if(fabs(elecEta) > 1.44 and fabs(elecEta) < 1.56) sf = 0.889;
    else if(fabs(elecEta) > 1.56 and fabs(elecEta) < 2.00) sf = 0.932;
    else if(fabs(elecEta) > 2.00 and fabs(elecEta) < 2.50) sf = 0.984;
    }
  else if(elecPt > 20.00 and elecPt < 30.00)
    {
    if(fabs(elecEta) > 0.0 and fabs(elecEta) < 0.8) sf = 0.980;
    else if(fabs(elecEta) > 0.8 and fabs(elecEta) < 1.44) sf = 0.955;
    else if(fabs(elecEta) > 1.44 and fabs(elecEta) < 1.56) sf = 0.996;
    else if(fabs(elecEta) > 1.56 and fabs(elecEta) < 2.00) sf = 0.970;
    else if(fabs(elecEta) > 2.00 and fabs(elecEta) < 2.50) sf = 1.027;
    }
  else if(elecPt > 30.00 and elecPt < 40.00)
    {
    if(fabs(elecEta) > 0.0 and fabs(elecEta) < 0.8) sf = 0.982;
    else if(fabs(elecEta) > 0.8 and fabs(elecEta) < 1.44) sf = 0.965;
    else if(fabs(elecEta) > 1.44 and fabs(elecEta) < 1.56) sf = 0.989;
    else if(fabs(elecEta) > 1.56 and fabs(elecEta) < 2.00) sf = 0.968;
    else if(fabs(elecEta) > 2.00 and fabs(elecEta) < 2.50) sf = 1.018;
    }
  else if(elecPt > 40.00 and elecPt < 50.00)
    {
    if(fabs(elecEta) > 0.0 and fabs(elecEta) < 0.8) sf = 0.985;
    else if(fabs(elecEta) > 0.8 and fabs(elecEta) < 1.44) sf = 0.974;
    else if(fabs(elecEta) > 1.44 and fabs(elecEta) < 1.56) sf = 0.961;
    else if(fabs(elecEta) > 1.56 and fabs(elecEta) < 2.00) sf = 0.989;
    else if(fabs(elecEta) > 2.00 and fabs(elecEta) < 2.50) sf = 1.012;
    }
  else if(elecPt > 50.00 and elecPt < 200.00)
    {
    if(fabs(elecEta) > 0.0 and fabs(elecEta) < 0.8) sf = 0.984;
    else if(fabs(elecEta) > 0.8 and fabs(elecEta) < 1.44) sf = 0.978;
    else if(fabs(elecEta) > 1.44 and fabs(elecEta) < 1.56) sf = 0.982;
    else if(fabs(elecEta) > 1.56 and fabs(elecEta) < 2.00) sf = 0.991;
    else if(fabs(elecEta) > 2.00 and fabs(elecEta) < 2.50) sf = 1.008;
    }
  return sf;
}

double muonSF(double muPt, double muEta)
{
  double sf = 1.0;
  if(muPt > 10.00 and muPt < 1000.00)
  {
    if(fabs(muEta) > 0.0 and fabs(muEta) < 0.90) sf = 0.9943;
    else if(fabs(muEta) > 0.9 and fabs(muEta) < 1.20) sf = 0.9933;
    else if(fabs(muEta) > 1.20 and fabs(muEta) < 2.50) sf = 1.0020;
  }
  return sf;
}

bool sortLeptonsInDescendingpT(LeptonInfo lep1, LeptonInfo lep2)
{
  return (lep1.pT > lep2.pT);
}

bool sortJetsInDescendingpT(JetInfo jet1, JetInfo jet2)
{
  return (jet1.pT > jet2.pT);
}

bool sortJetVectorsInDescendingpT(TLorentzVector jet1, TLorentzVector jet2)
{
  return (jet1.Pt() > jet2.Pt());
}

bool sortPhotonsInDescendingpT(PhotonInfo pho1, PhotonInfo pho2)
{
  return (pho1.pT > pho2.pT);
}

bool sortMatchedLeptonsInDescendingpT(MatchedLeptonInfo mlep1, MatchedLeptonInfo mlep2)
{
  return (mlep1.pT > mlep2.pT);
}
