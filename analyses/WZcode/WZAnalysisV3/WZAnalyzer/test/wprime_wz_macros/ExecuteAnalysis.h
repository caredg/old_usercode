//---------------------------------------------
// Author: Edgar Carrera
// 2010-01-13
// This macro will be used for analyzing
// Wprime -> WZ -> lllnu events
// It works on an root-uple of events with
// variables created by the "official" CMS WZ code
//---------------------------------------------
#ifndef _ExecuteAnalysis_h_
#define _ExecuteAnalysis_h_

#include <stdio.h>
#include <math.h>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TRandom.h"
#include "TMath.h"
#include "TEntryList.h"
#include <vector>
#include <algorithm>
#include <limits>

//gROOT->Reset();

using namespace std;




// ++++++++++++++data structure to accumulate info about the file
struct InputFile{
  float x_sect; // cross-section in pb
  float Nprod_evt; // # of events produced
  float weight; // cross-section * integrated luminosity / (# of events produced)
  string pathname;
  string description; // sample description
  TTree * tree; // pointer to ROOT file
  //
  InputFile(float cross_sect)
  {
    x_sect = cross_sect; Nprod_evt = 0; weight = 0; tree = 0;
  }
};

// +++++++++++++++++++location of data files and samples info
//string top_level_dir = "/localdata/ecarrera/analyses/wprime_to_wz/root_uples/";
string top_level_dir = "/media/sda2/fantasia/work/Wprime/data/";

// +++++++++++++++++++Variables to store Branch Addresses:
Int_t           W_flavor;
Int_t           Z_flavor;
Int_t           numberOfZs;
Float_t         WZ_invMassMinPz;
vector<float>   *jet_energy;
vector<float>   *jet_et;
vector<float>   *jet_eta;
vector<float>   *jet_mass;
vector<float>   *jet_phi;
vector<float>   *jet_pt;
vector<float>   *jet_px;
vector<float>   *jet_py;
vector<float>   *jet_pz;
vector<float>   *muon_eta;
vector<float>   *muon_pt;
vector<float>   *electron_eta;
vector<float>   *electron_pt;
Float_t         met_et;
Float_t         Z_mass;
Float_t         Z_pt;
Float_t         W_pt;
Int_t           W_leptonIndex;
Int_t           Z_leptonIndex1;
Int_t           Z_leptonIndex2;

vector<float>   *muon_innerD0;
vector<float>   *muon_innerD0Error;
vector<float>   *muon_caloIso;
vector<float>   *muon_trackIso;

vector<float>   *electron_caloIso;
vector<float>   *electron_deltaEtaIn;
vector<float>   *electron_deltaPhiIn;
vector<float>   *electron_eOverP;
vector<float>   *electron_hOverE;
vector<float>   *electron_sigmaEtaEta;
vector<float>   *electron_trackIso;

// +++++++++++++++++++useful constants
const double PI    = 2.0 * acos(0.);
const double TWOPI = 2.0 * PI;
const bool debugme = false; //print stuff if active
const int PDGMUON = 13;
const int PDGELEC = 11;
const int PDGW = 24;
const int PDGZ = 23;

const float PDGZMASS = 91.1876; //GeV


// +++++++++++++++++++ Histogram Definitions
const int Num_histo_sets = 8; //matches the number of cuts
TH1F * hWZInvMass[Num_histo_sets];
TH1F * hJetMult[Num_histo_sets];
TH1F * hJetsDeltaEta[Num_histo_sets];
TH1F * hHt[Num_histo_sets];
TH1F * hWpt[Num_histo_sets];
TH1F * hZpt[Num_histo_sets];

TEntryList *cutlist[Num_histo_sets];

// +++++++++++++++++++General Cut values
const int cutMaxNumZs = 1;
const float cutMinJetE = 300;
const float cutMinJetPt = 10;
const float cutMinJetsDeltaEta = 4;
const int minNumLeptons = 3;

// +++++++++++++++++++Ht Cuts
const float minHt = 100;
const float minHtMet = 0;

// +++++++++++++++++++ Muon ID Cuts

const float cutZMuon_pt = 10.;
const float cutWMuon_pt = 20.;
const float cutWmuD0 = 8.;
const float cutWmuCombRelIso = 0.1;

const float maxMuonEta = 2.5;
const float minMuonPt = 10; //GeV

// +++++++++++++++++++ Electron ID Cuts
const float minElecPt = 10; //GeV

const float maxEtaBarrel = 1.479;
const float minEtaEndcap = 1.55;
const float maxElecEta = 2.5;

const float cutZElectron_pt = 15;
const float cutWElectron_pt = 20;

const float cutWElCombRelIso = 0.1;

//const float NOCUT = numeric_limits<float>::max();
//just make a big number for now
const float NOCUT = 99999999;


// Cut arrays are for {barrel, endcap}
const float cutDeltaEtaIn[]  = {0.005, 0.007};
const float cutDeltaPhiIn[]  = {0.040, 0.040};
const float cutSigmaEtaEta[] = {0.011, NOCUT};
const float cutEOverP[]      = {0.760, 0.680};
const float cutHOverE[]      = {0.016, 0.025};
const float cutCalRelIso[]   = {0.100, 0.160};
const float cutTrackRelIso[] = {0.100, 0.100};

// ++++++++++++++++++++ Other Cuts
const float minSeparation = 0.1;
const float minJetEta = 2.;
const float ZMASSRES = NOCUT;

// +++++++++++++++++++ Declare the methods that we use:
void RecruitOrderedFiles(vector<InputFile> & files, const int& Nfiles,
                         const int& filenum_low, const int& filenum_step,
                         const string& mask1,
                         const string& mask2, const string& file_desc);
string convertIntToStr(int number);
void getEff(float & eff, float & deff, float Num, float Denom);
double deltaEta(double eta1, double eta2);
double deltaPhi(double phi1, double phi2);
bool inBarrel(float eta);
bool inEndCap(float eta);

void Declare_Histos();
void Declare_Lists();
double deltaR(double eta1, double phi1, double eta2, double phi2);
int Check_Files(unsigned Nfiles, vector<InputFile> & files);
void Load_Input_Files(string file_desc,vector<InputFile> & files,float lumiPb);
int Load_Cross_Sections(vector<InputFile> & wzjj_files,
                        vector<InputFile> & wprime400_files);
void Set_Branch_Addresses(TTree* WZtree);
void Fill_Histos(int index, float weight);
void saveHistos(TFile * fout, string dir);
void printSummary(ofstream & out, const string& dir,
                  const float& Nexp_evt, float Nexp_evt_cut[]);
void Tabulate_Me(int Num_surv_cut[], int& cut_index,const float& weight,const int& evtnum);
void Get_Distributions(vector<InputFile>& files,TFile *fout, 
                       string dir, ofstream & out);
void ExecuteAnalysis();

//methods for the cuts
bool PassTriggers_Cut();
bool HasValidWandZ_Cut();
bool ExeedMaxNumberOfZs_Cut(const int& max_num_Zs);
bool HasTwoEnergeticForwardJets_Cut(const float& cutMinJetE,
                                    const float& cutMinJetPt,
                                    const float& cutMinJetsDeltaEta);
bool PassMuonCut();
bool PassMuonSip(int index);
bool PassMuonRelIso(int index);
bool PassElecCut();
bool PassElecEtaDepCut(int index, int parent);
bool PassElecPtCut(float pt,float eta,float parent);
bool PassElecdEtaCut(float dEta, float eta);
bool PassElecdPhiCut(float dPhi, float eta);
bool PassElecSigmannCut(float sigmann, float eta);
bool PassElecEPCut(float EP,float eta);
bool PassElecHECut(float HE, float eta);
bool PassElecCalIsoCut(float calIso, float eta);
bool PassElecTrkIsoCut(float trkIso,float eta);
bool PassZCut();
bool PassHtCut();
bool PassHtMetCut();

float Calc_Ht();
void Compare_Histos(vector<InputFile>& files,TFile* fout);
#endif//#define _ExecuteAnalysis_h_
