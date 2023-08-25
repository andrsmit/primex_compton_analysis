
int PHASE_VAL;

bool IS_BE_TARGET, IS_FIELD_OFF;
int  BEAM_CURRENT;

bool DRAW_THEORY,    DRAW_ACCEPTANCE;
bool DRAW_FITS_TAGH, DRAW_FITS_TAGM;
bool USE_F_ACC;

char TARGET_STR[256], FIELD_STR[256];

/*********************************************/
// photon flux:

double        tagh_flux[274],        tagm_flux[102];
double       tagh_fluxE[274],       tagm_fluxE[102];
double  tagh_flux_empty[274],  tagm_flux_empty[102];
double tagh_fluxE_empty[274], tagm_fluxE_empty[102];

/*********************************************/
// counter energies:

double endpoint_energy, endpoint_energy_calib;
double tagh_en[274], tagm_en[102];

/*********************************************/
// acceptance:

double  tagh_acc[274],  tagm_acc[102];
double tagh_accE[274], tagm_accE[102];

/*********************************************/
// theory:

TF1 *f_theory;
TF1 *f_pair_cs, *f_trip_cs;

/*********************************************/
// target:

double n_e, n_Z;
double f_abs;

double mb = 1.e-27;

/*********************************************/
// data:

TString root_fname, empty_target_root_fname;
TString hname_tagh, hname_tagm;

double  tagh_yield[274],  tagm_yield[102];
double tagh_yieldE[274], tagm_yieldE[102];

/*********************************************/
// fitting:

TCanvas *canvas1;
TPad *pad_lin, *pad_log;

TCanvas *canvas_draw;

vector<int> bad_counters_tagh, bad_counters_tagm;

double tagh_chi2[274], tagm_chi2[102];

double      tagh_comp_frac[274],      tagh_comp_fracE[274];
double      tagm_comp_frac[102],      tagm_comp_fracE[102];

double      tagh_pair_frac[274],     tagh_pair_fracE[274];
double      tagm_pair_frac[102],     tagm_pair_fracE[102];
double  tagh_pair_frac_exp[274],  tagm_pair_frac_exp[102];

double     tagh_empty_frac[274],    tagh_empty_fracE[274];
double     tagm_empty_frac[102],    tagm_empty_fracE[102];
double tagh_empty_frac_exp[274], tagm_empty_frac_exp[102];

double bin_size = 16./2000.;
int rebins, n_mev;

/*********************************************/
// simulation:

TString comp_root_dir_name, pair_root_dir_name, trip_root_dir_name;

TString hname_tagh_comp, hname_tagm_comp;
TString hname_tagh_pair, hname_tagm_pair;
TString hname_tagh_trip, hname_tagm_trip;

TString comp_mc_dir, pair_mc_dir, trip_mc_dir;
