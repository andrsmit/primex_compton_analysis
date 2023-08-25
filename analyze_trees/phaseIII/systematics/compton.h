using namespace std;

#include <stdio.h>
#include <fstream>
#include <math.h>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>

#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TSystem.h>
#include <TDirectory.h>
#include <TMath.h>

#define MAX_CAND 2000
#define N_TAGH_COUNTERS 274
#define N_TAGM_COUNTERS 102
#define N_SIGS 16

//----------   Numerical Constants   ----------//

const double c   = 29.9792; // [cm/ns]
const double pi  = 3.1415926535;
const double m_e = 0.510998928e-3;

const double m_fcalX =    0.408;
const double m_fcalY =    0.027;
const double m_fcalZ =  624.906;

const double m_ccalX =    0.135;
const double m_ccalY =    0.135;
const double m_ccalZ = 1279.376;

double m_beamX, m_beamY, m_beamZ;

//----------   Function Declarations   ----------//

int main(int argc, char **argv);

void reset_event();
void read_event(int evt);
void dump_event() {};

void init_histograms();
void reset_histograms();
void write_histograms(char *name);

void set_cuts();
void load_constants(int group, bool is_first);

void compton_analysis(int run);

//----------   Histograms   ----------//

TH2F *h_deltaK_tagh_fcalE[20],      *h_deltaK_tagm_fcalE[20];
TH2F *h_deltaK_tagh_cut_fcalE[20],  *h_deltaK_tagm_cut_fcalE[20];

TH2F *h_deltaK_tagh_ccalE[13],      *h_deltaK_tagm_ccalE[13];
TH2F *h_deltaK_tagh_cut_ccalE[13],  *h_deltaK_tagm_cut_ccalE[13];

double cut_sigmas[16] = {
1.0, 
2.0, 
3.0, 
3.5, 
4.0, 
4.5, 
5.0, 
5.5, 
6.0, 
6.5, 
7.0, 
8.0, 
9.0, 
10.0, 
15.0, 
20.0
};

TH2F *h_deltaK_tagh_sigE[16],       *h_deltaK_tagm_sigE[16];
TH2F *h_deltaK_tagh_cut_sigE[16],   *h_deltaK_tagm_cut_sigE[16];

TH2F *h_deltaK_tagh_sigPhi[16],     *h_deltaK_tagm_sigPhi[16];
TH2F *h_deltaK_tagh_cut_sigPhi[16], *h_deltaK_tagm_cut_sigPhi[16];

TH2F *h_deltaK_tagh_sigK[16],       *h_deltaK_tagm_sigK[16];

TH2F *h_deltaK_tagh_fcalT[10],      *h_deltaK_tagm_fcalT[10];
TH2F *h_deltaK_tagh_cut_fcalT[10],  *h_deltaK_tagm_cut_fcalT[10];

TH2F *h_deltaK_tagh_ccalT[10],      *h_deltaK_tagm_ccalT[10];
TH2F *h_deltaK_tagh_cut_ccalT[10],  *h_deltaK_tagm_cut_ccalT[10];

static const int N_FID_CUTS = 20;
TH2F *h_deltaK_tagh_fcalfid[N_FID_CUTS];
TH2F *h_deltaK_tagm_fcalfid[N_FID_CUTS];
TH2F *h_deltaK_tagh_ccalfid[N_FID_CUTS];
TH2F *h_deltaK_tagm_ccalfid[N_FID_CUTS];

TH2F *h_deltaK_tagh_fcal_phi[8],   *h_deltaK_tagm_fcal_phi[8];
TH2F *h_deltaK_tagh_ccal_phi[8],   *h_deltaK_tagm_ccal_phi[8];
TH2F *h_deltaK_tagh_fcal_layer[8], *h_deltaK_tagm_fcal_layer[8];
TH2F *h_deltaK_tagh_ccal_layer[5], *h_deltaK_tagm_ccal_layer[5];

TH2F *h_deltaK_tagh_cut_fcal_phi[8],   *h_deltaK_tagm_cut_fcal_phi[8];
TH2F *h_deltaK_tagh_cut_ccal_phi[8],   *h_deltaK_tagm_cut_ccal_phi[8];
TH2F *h_deltaK_tagh_cut_fcal_layer[8], *h_deltaK_tagm_cut_fcal_layer[8];
TH2F *h_deltaK_tagh_cut_ccal_layer[5], *h_deltaK_tagm_cut_ccal_layer[5];

TH2F *h_xy_fcal_phi[8],   *h_xy_ccal_phi[8];
TH2F *h_xy_fcal_layer[8], *h_xy_ccal_layer[8];

TH2F *h_deltaK_tagh_single,    *h_deltaK_tagm_single;

TH2F *h_deltaK_tagh_main,       *h_deltaK_tagm_main;
TH2F *h_deltaK_tagh_main_extra, *h_deltaK_tagm_main_extra;
TH2F *h_deltaK_tagh_acc,        *h_deltaK_tagm_acc;
TH2F *h_deltaK_tagh_acc_single, *h_deltaK_tagm_acc_single;
TH2F *h_deltaK_tagh_acc_extra,  *h_deltaK_tagm_acc_extra;

//----------   Data Objects   ----------//

char rootTree_pathName[256], rootFile_pathName[256];

// DeltaE mu function is 3rd order polynomial
// DeltaE sig/E function is [0] + [1]/sqrt(x) + [2]/x

TF1 *f_deltaE_mu,   *f_deltaE_sig;

double deltaE_mu_p0,  deltaE_mu_p1,  deltaE_mu_p2,  deltaE_mu_p3;
double deltaE_sig_p0, deltaE_sig_p1, deltaE_sig_p2;

// DeltaPhi mu function is 3rd order polynomial
// DeltaPhi sig function is 3rd order polynomial

TF1 *f_deltaPhi_mu, *f_deltaPhi_sig;

double deltaPhi_mu_p0,  deltaPhi_mu_p1,  deltaPhi_mu_p2,  deltaPhi_mu_p3;
double deltaPhi_sig_p0, deltaPhi_sig_p1, deltaPhi_sig_p2, deltaPhi_sig_p3;

// DeltaK mu function is 3rd order polynomial
// DeltaK sig function is 3rd order polynomial

TF1 *f_deltaK_mu,   *f_deltaK_sig;

double deltaK_mu_p0,  deltaK_mu_p1,  deltaK_mu_p2,  deltaK_mu_p3;
double deltaK_sig_p0, deltaK_sig_p1, deltaK_sig_p2, deltaK_sig_p3;

double FCAL_ENERGY_CUT,  CCAL_ENERGY_CUT;
double FCAL_RF_TIME_CUT, CCAL_RF_TIME_CUT;
double deltaE_cut_sig, deltaPhi_cut_sig, deltaK_cut_sig;

double fcal_inner_layer_cut, ccal_inner_layer_cut;

int first;

TTree *tree;

int eventnum, n_candidates;
double rfTime;

int bunch_val[MAX_CAND];
int tag_counter[MAX_CAND], tag_sys[MAX_CAND];

double eb[MAX_CAND], tb[MAX_CAND];

double fcal_e[MAX_CAND];
double fcal_x[MAX_CAND], fcal_y[MAX_CAND], fcal_z[MAX_CAND];
double fcal_t[MAX_CAND];

double ccal_e[MAX_CAND];
double ccal_x[MAX_CAND], ccal_y[MAX_CAND], ccal_z[MAX_CAND];
double ccal_t[MAX_CAND];

double deltaE[MAX_CAND], deltaK[MAX_CAND], deltaPhi[MAX_CAND];
double deltaT[MAX_CAND], deltaR[MAX_CAND];
double deltaK2[MAX_CAND];
