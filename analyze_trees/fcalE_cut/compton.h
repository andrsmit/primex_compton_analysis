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
#define N_CUTS 20

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

TH1F *h_tagh[N_CUTS][N_TAGH_COUNTERS], *h_tagm[N_CUTS][N_TAGM_COUNTERS];


//----------   Data Objects   ----------//


char      cut_pathName[256];
char rootTree_pathName[256], rootFile_pathName[256];

double cut_energies[N_CUTS];

double deltaE_mu_tagh[N_TAGH_COUNTERS],   deltaE_sig_tagh[N_TAGH_COUNTERS];
double deltaE_mu_tagm[N_TAGM_COUNTERS],   deltaE_sig_tagm[N_TAGM_COUNTERS];

double deltaPhi_mu_tagh[N_TAGH_COUNTERS], deltaPhi_sig_tagh[N_TAGH_COUNTERS];
double deltaPhi_mu_tagm[N_TAGM_COUNTERS], deltaPhi_sig_tagm[N_TAGM_COUNTERS];

double deltaK_mu_tagh[N_TAGH_COUNTERS],   deltaK_sig_tagh[N_TAGH_COUNTERS];
double deltaK_mu_tagm[N_TAGM_COUNTERS],   deltaK_sig_tagm[N_TAGM_COUNTERS];


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
