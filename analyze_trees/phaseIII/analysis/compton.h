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

#define MAX_CAND 2000
#define N_TAGH_COUNTERS 274
#define N_TAGM_COUNTERS 102

//----------   Numerical Constants   ----------//

const double c   = 29.9792; // [cm/ns]
const double pi  = 3.1415926535;
const double m_e = 0.510998928e-3;

const double m_fcalX =    0.408;
const double m_fcalY =    0.027;
const double m_fcalZ =  624.906;

const double m_ccalX =    0.184;
const double m_ccalY =    0.110;
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

TH1F *h_fcal_rf_dt, *h_ccal_rf_dt, *h_beam_rf_dt, *h_beam_rf_dt_cut;
TH1F *h_n_cands;
TH2F *h_deltaE_tagh,      *h_deltaE_tagm;
TH2F *h_deltaPhi_tagh,    *h_deltaPhi_tagm;
TH2F *h_deltaK_tagh,      *h_deltaK_tagm;
TH2F *h_deltaK_tagh_cut,  *h_deltaK_tagm_cut;
TH2F *h_deltaK2_tagh,     *h_deltaK2_tagm;
TH2F *h_deltaK2_tagh_cut, *h_deltaK2_tagm_cut;
TH2F *h_fcal_xy,          *h_ccal_xy;
TH2F *h_deltaE_vs_deltaK, *h_deltaE_vs_deltaPhi, *h_deltaPhi_vs_deltaK;
TH2F *h_deltaE_vs_deltaK2;
TH2F *h_deltaH_tagh[4],   *h_deltaH_tagm[4];
TH2F *h_deltaKE_tagh[4],  *h_deltaKE_tagm[4];
TH2F *h_esum_vs_ecomp[4];

TH1F *h_tagh_background, *h_tagh_background_lo, *h_tagh_background_hi;
TH1F *h_tagm_background, *h_tagm_background_lo, *h_tagm_background_hi;

TH2F *h_deltaE_vs_deltaK2_main, *h_deltaE_vs_deltaK2_acc;

TH1F *h_tagged_bkgd, *h_tagged_bkgd_1GeV_cut;
TH2F *h_deltaE_vs_esum;

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

// DeltaK2 mu function is 3rd order polynomial
// DeltaK2 sig function is 3rd order polynomial

TF1 *f_deltaK2_mu,   *f_deltaK2_sig;

double deltaK2_mu_p0,  deltaK2_mu_p1,  deltaK2_mu_p2,  deltaK2_mu_p3;
double deltaK2_sig_p0, deltaK2_sig_p1, deltaK2_sig_p2, deltaK2_sig_p3;

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
