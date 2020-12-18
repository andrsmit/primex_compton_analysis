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

#define MAX_CAND 100
#define N_TAGH_COUNTERS 274
#define N_TAGM_COUNTERS 102
#define N_SIGS 14

//----------   Numerical Constants   ----------//

const double c  = 29.9792; // [cm/ns]
const double pi = 3.1415926535;
const double me = 0.510998928e-3;


//----------   Function Declarations   ----------//

int main(int argc, char **argv);

void reset_event();
void read_event(int evt);
void dump_event() {};

void init_histograms();
void reset_histograms();
void write_histograms(char *name);

void set_cuts();
void load_constants(int group);

void compton_analysis(int run);



//----------   Histograms   ----------//

TH1F *h_n_ccal_total,        *h_n_fcal_total,        *h_n_show_total;
TH1F *h_n_ccal[N_SIGS],      *h_n_fcal[N_SIGS],      *h_n_show[N_SIGS];
TH1F *h_ccal_extraE[N_SIGS], *h_fcal_extraE[N_SIGS], *h_show_extraE[N_SIGS];
TH1F *h_tagh[N_SIGS][N_TAGH_COUNTERS], *h_tagm[N_SIGS][N_TAGM_COUNTERS];


//----------   Data Objects   ----------//


char      cut_pathName[256];
char rootTree_pathName[256], rootFile_pathName[256];

double cut_sigmas[N_SIGS];

// DeltaE mu function is 3rd order polynomial
// DeltaE sig/E function is [0] + [1]/sqrt(x) + [2]/x

TF1 *f_deltaE_mu,   *f_deltaE_sig;

double deltaE_mu_p0  =  1.64677e-01;
double deltaE_mu_p1  = -5.68233e-02;
double deltaE_mu_p2  =  4.79290e-03;
double deltaE_mu_p3  = -1.59059e-04;

double deltaE_sig_p0 =  1.87138e-02;
double deltaE_sig_p1 = -2.37318e-02;
double deltaE_sig_p2 =  7.30177e-02;

// DeltaPhi mu function is 3rd order polynomial
// DeltaPhi sig function is 3rd order polynomial

TF1 *f_deltaPhi_mu, *f_deltaPhi_sig;

double deltaPhi_mu_p0  =  1.79064e+02;
double deltaPhi_mu_p1  =  2.56662e-01;
double deltaPhi_mu_p2  = -3.43324e-02;
double deltaPhi_mu_p3  =  1.39042e-03;

double deltaPhi_sig_p0 =  8.12410e+00;
double deltaPhi_sig_p1 = -5.27691e-01;
double deltaPhi_sig_p2 =  2.97087e-02;
double deltaPhi_sig_p3 = -3.43114e-04;

// DeltaK mu function is 3rd order polynomial
// DeltaK sig function is 3rd order polynomial

TF1 *f_deltaK_mu,   *f_deltaK_sig;

double deltaK_mu_p0  =  5.35336e-02;
double deltaK_mu_p1  =  4.75062e-03;
double deltaK_mu_p2  = -3.51716e-03;
double deltaK_mu_p3  =  1.65015e-04;

double deltaK_sig_p0 =  1.98161e-01;
double deltaK_sig_p1 = -3.59916e-02;
double deltaK_sig_p2 =  7.11062e-03;
double deltaK_sig_p3 = -3.10044e-04;

/*
double deltaE_mu_tagh[N_TAGH_COUNTERS],   deltaE_sig_tagh[N_TAGH_COUNTERS];
ouble deltaE_mu_tagm[N_TAGM_COUNTERS],   deltaE_sig_tagm[N_TAGM_COUNTERS];
double deltaPhi_mu_tagh[N_TAGH_COUNTERS], deltaPhi_sig_tagh[N_TAGH_COUNTERS];
double deltaPhi_mu_tagm[N_TAGM_COUNTERS], deltaPhi_sig_tagm[N_TAGM_COUNTERS];

double deltaK_mu_tagh[N_TAGH_COUNTERS],   deltaK_sig_tagh[N_TAGH_COUNTERS];
double deltaK_mu_tagm[N_TAGM_COUNTERS],   deltaK_sig_tagm[N_TAGM_COUNTERS];
*/


double FCAL_ENERGY_CUT,  CCAL_ENERGY_CUT;
double FCAL_RF_TIME_CUT, CCAL_RF_TIME_CUT;
double deltaE_cut_sig, deltaPhi_cut_sig, deltaK_cut_sig;




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
