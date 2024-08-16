#ifndef _COMPTONCS_
#define _COMPTONCS_

#define N_TAGH_COUNTERS 274
#define N_TAGM_COUNTERS 102

double tagh_flux_unc[274], tagm_flux_unc[102];


// Useful variable switches:

int  PHASE_VAL    = 1;
int  BEAM_CURRENT = 200;
bool IS_BE_TARGET = true;
bool IS_FIELD_OFF = true;

TString TARGET_STR = "Be", FIELD_STR = "FIELDOFF";

/*************************************************************************/
// Data:

TString root_fname, empty_target_root_fname;
TString hname_tagh, hname_tagm;

double tagh_yield[N_TAGH_COUNTERS], tagh_yieldE[N_TAGH_COUNTERS];
double tagm_yield[N_TAGM_COUNTERS], tagm_yieldE[N_TAGM_COUNTERS];

double tagh_cs[N_TAGH_COUNTERS], tagh_csE[N_TAGH_COUNTERS];
double tagm_cs[N_TAGM_COUNTERS], tagm_csE[N_TAGM_COUNTERS];

void get_compton_yield(vector<int> &tagh_counter_vec, vector<int> &tagm_counter_vec);

/*********************************************/
// fitting:

bool FIT_EMPTY   = false;
bool FIT_TRIPLET = false;

bool FIT_USING_BE_EMPTY = true;

bool DEBUG_FITS     = false;
bool DRAW_FITS_TAGH = false;
bool DRAW_FITS_TAGM = false;
bool SAVE_FITS_TAGH = false;
bool SAVE_FITS_TAGM = false;

double DELTA_K_FIT_SIGMA = 0.;

TCanvas *canvas_fit;
TPad *canvas_fit_lin, *canvas_fit_log;

TCanvas *canvas_draw;
TPad *canvas_draw_top, *canvas_draw_bot;

TCanvas *c_debug, *canvas_pull, *canvas_pull_dist;

vector<int> bad_counters_tagh, bad_counters_tagm;
int n_bad_counters_tagh = 0, n_bad_counters_tagm = 0;

vector<TString> mc_templates = {"compton", "pair", "triplet", "empty"};
map<TString, unsigned int> hist_color_map;

map<TString, vector<pair<double,double>>> tagh_fit_fraction, tagh_fit_fraction_exp;
map<TString, vector<pair<double,double>>> tagm_fit_fraction, tagm_fit_fraction_exp;

vector<double> tagh_yieldfit_chi2, tagm_yieldfit_chi2;
vector<pair<double,double>> tagh_yieldfit_pull, tagm_yieldfit_pull;

double bin_size = 16.0/2000.0;
int rebins, n_mev;

/*********************************************/
// simulation:

TString comp_root_dir_name, pair_root_dir_name;

TString hname_tagh_comp, hname_tagm_comp;
TString hname_tagh_pair, hname_tagm_pair;
TString hname_tagh_trip, hname_tagm_trip;

TString comp_mc_dir, pair_mc_dir, triplet_mc_dir;

bool USE_F_ACC = true;

TF1 *f_acc;
TCanvas *canvas_acc;

bool DRAW_ACCEPTANCE   = false;
bool CALC_ACC_FROM_FIT = false; // use a di-Gaussian fit to DeltaK Distribution to get acceptance

double tagh_acc[N_TAGH_COUNTERS], tagh_accE[N_TAGH_COUNTERS];
double tagm_acc[N_TAGM_COUNTERS], tagm_accE[N_TAGM_COUNTERS];

void get_compton_acc();

/*************************************************************************/
// Get energies for each tagger counter:

// run-dependent electron beam energy:
double endpoint_energy, endpoint_energy_calib;

// pathname of files containing fraction energy values for each tagger counter:
TString tagh_xscale_fname, tagm_xscale_fname;

// arrays storing the mid-point energy values of each counter:
double tagh_en[274], tagm_en[102];

// For phase 1 specifically (needed for proper scaling of e+e- simulation):

double endpoint_energy_phase1 = 11.6061, endpoint_energy_calib_phase1 = 11.6061;
double tagh_en_phase1[N_TAGH_COUNTERS], tagm_en_phase1[N_TAGM_COUNTERS];
TString tagh_xscale_fname_phase1 = 
	"/work/halld/home/andrsmit/primex_compton_analysis/photon_flux/phase1/primex_tagh.txt";
TString tagm_xscale_fname_phase1 = 
	"/work/halld/home/andrsmit/primex_compton_analysis/photon_flux/phase1/primex_tagm.txt";

int get_counter_energies();

/*************************************************************************/
// Read in photon flux:

TString              tagh_flux_fname,              tagm_flux_fname;
TString empty_target_tagh_flux_fname, empty_target_tagm_flux_fname;

double tagh_flux[N_TAGH_COUNTERS], tagh_fluxE[N_TAGH_COUNTERS];
double tagm_flux[N_TAGM_COUNTERS], tagm_fluxE[N_TAGM_COUNTERS];
double tagh_flux_empty[N_TAGH_COUNTERS], tagh_fluxE_empty[N_TAGH_COUNTERS];
double tagm_flux_empty[N_TAGM_COUNTERS], tagm_fluxE_empty[N_TAGM_COUNTERS];

int get_flux();

/*************************************************************************/
// NLO Theory Calculation from NIST and from primex_compton generator:

TString nist_theory_fname = "/home/andrsmit/root_macros/compton/nist_cs.dat";

TString theory_cs_pathName;

TF1 *f_theory, *f_nist;

TGraph *g_theory;
TGraph *g_theory_min, *g_theory_max, *g_theory_band;
TGraph *g_dev_min,    *g_dev_max,    *g_dev_band;
TGraph *g_nist;

TCanvas *c_theory;

bool USE_NIST_CALCULATION = false;
bool DRAW_THEORY          = true;

int get_theory_calc();

/*************************************************************************/
// Get e+e- pair production cross section (from NIST):

TString pair_cs_fname;

TF1 *f_pair_cs, *f_triplet_cs;

TCanvas *c_pair_cs;

void get_pair_cs();

/*************************************************************************/
// Get target parameters:

double n_e = 0.0;
double n_Z = 0.0;
double n_A = 0.0;
double target_thickness = 0.0;

double mb = 1.e-27;

double f_abs = 1.0;

void get_target_parameters();

/*************************************************************************/
// Results:

TGraphErrors *gCS;
TCanvas *cCS;

void calc_cs(int tag_sys, int counter, double &loc_cs, double &loc_csE);
void plot_cs(vector<int> tagh_counter_vec, vector<int> tagm_counter_vec);
void plot_fit_fractions(vector<int> tagh_counter_vec, vector<int> tagm_counter_vec);
void plot_chi2(vector<int> tagh_counter_vec, vector<int> tagm_counter_vec);

bool isNaN(double x);

#endif
