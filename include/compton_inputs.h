#ifndef _COMPTONINPUTS_
#define _COMPTONINPUTS_

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

double tagh_yield[274],  tagm_yield[102];
double tagh_yieldE[274], tagm_yieldE[102];

bool USE_F_ACC = true;

double tagh_cs[274],  tagm_cs[102];
double tagh_csE[274], tagm_csE[102];

/*********************************************/
// fitting:

bool DEBUG_FITS = false;
bool DRAW_FITS_TAGH = false;
bool DRAW_FITS_TAGM = false;
bool SAVE_FITS_TAGH = false;
bool SAVE_FITS_TAGM = false;

TCanvas *canvas_fit;
TPad *canvas_fit_lin, *canvas_fit_log;

TCanvas *canvas_draw;
TPad *canvas_draw_top, *canvas_draw_bot;

TCanvas *c_debug, *canvas_pull, *canvas_pull_dist;

vector<int> bad_counters_tagh, bad_counters_tagm;
int n_bad_counters_tagh = 0, n_bad_counters_tagm = 0;

double tagh_yieldfit_chi2[274],       tagm_yieldfit_chi2[102];
double tagh_yieldfit_pull_mean[274],  tagm_yieldfit_pull_mean[102];
double tagh_yieldfit_pull_stdev[274], tagm_yieldfit_pull_stdev[102];

double tagh_compton_fraction[274],    tagm_compton_fraction[102];
double tagh_compton_fractionE[274],   tagm_compton_fractionE[102];

double tagh_pair_fraction[274],       tagm_pair_fraction[102];
double tagh_pair_fractionE[274],      tagm_pair_fractionE[102];
double tagh_pair_fraction_exp[274],   tagm_pair_fraction_exp[102];

double tagh_empty_fraction[274],      tagm_empty_fraction[102];
double tagh_empty_fractionE[274],     tagm_empty_fractionE[102];
double tagh_empty_fraction_exp[274],  tagm_empty_fraction_exp[102];

double bin_size = 16.0/2000.0;
int rebins, n_mev;

/*********************************************/
// simulation:

TString comp_root_dir_name, pair_root_dir_name;

TString hname_tagh_comp, hname_tagm_comp;
TString hname_tagh_pair, hname_tagm_pair;
TString hname_tagh_trip, hname_tagm_trip;

TString comp_mc_dir, pair_mc_dir, trip_mc_dir;

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
double tagh_en_phase1[274], tagm_en_phase1[102];
TString tagh_xscale_fname_phase1 = "/work/halld/home/andrsmit/primex_compton_analysis/photon_flux/phase1/primex_tagh.txt";
TString tagm_xscale_fname_phase1 = "/work/halld/home/andrsmit/primex_compton_analysis/photon_flux/phase1/primex_tagm.txt";

// routine to read input files and get tagger counter energies:
int get_counter_energies() {
	
	// first check that input filenames exist:
	ifstream inf_tagh(tagh_xscale_fname.Data());
	ifstream inf_tagm(tagm_xscale_fname.Data());
	
	if(!inf_tagh.good() || !inf_tagm.good()) {
		return 1;
	}
	
	int a; double b, c;
	
	for(int i=0; i<274; i++) {
		inf_tagh >> a >> b >> c;
		double deltaE = endpoint_energy - endpoint_energy_calib;
		double emin   = b * endpoint_energy_calib  +  deltaE;
		double emax   = c * endpoint_energy_calib  +  deltaE;
		tagh_en[i] = 0.5 * (emin + emax);
	}
	inf_tagh.close();
	
	for(int i=0; i<102; i++) {
		inf_tagm >> a >> b >> c;
		double deltaE = endpoint_energy - endpoint_energy_calib;
		double emin   = b * endpoint_energy_calib  +  deltaE;
		double emax   = c * endpoint_energy_calib  +  deltaE;
		tagm_en[i] = 0.5 * (emin + emax);
	}
	inf_tagm.close();
	
	// Get Phase 1 Energy Bins:
	
	inf_tagh.open(tagh_xscale_fname_phase1.Data());
	for(int i=0; i<274; i++) {
		inf_tagh >> a >> b >> c;
		double deltaE = endpoint_energy_phase1 - endpoint_energy_calib_phase1;
		double emin   = b * endpoint_energy_calib_phase1  +  deltaE;
		double emax   = c * endpoint_energy_calib_phase1  +  deltaE;
		tagh_en_phase1[i] = 0.5 * (emin + emax);
	}
	inf_tagh.close();
	
	inf_tagm.open(tagm_xscale_fname_phase1.Data());
	for(int i=0; i<102; i++) {
		inf_tagm >> a >> b >> c;
		double deltaE = endpoint_energy_phase1 - endpoint_energy_calib_phase1;
		double emin   = b * endpoint_energy_calib_phase1  +  deltaE;
		double emax   = c * endpoint_energy_calib_phase1  +  deltaE;
		tagm_en_phase1[i] = 0.5 * (emin + emax);
	}
	inf_tagm.close();
	
	return 0;
}

/*************************************************************************/
// Read in photon flux:

TString              tagh_flux_fname,              tagm_flux_fname;
TString empty_target_tagh_flux_fname, empty_target_tagm_flux_fname;

double tagh_flux[274], tagh_fluxE[274];
double tagm_flux[102], tagm_fluxE[102];
double tagh_flux_empty[274], tagh_fluxE_empty[274];
double tagm_flux_empty[274], tagm_fluxE_empty[274];

int get_flux() 
{
	// first check that input filenames exist:
	ifstream inf_tagh(tagh_flux_fname.Data());
	ifstream inf_tagm(tagm_flux_fname.Data());
	
	if(!inf_tagh.good() || !inf_tagm.good()) {
		return 1;
	}
	
	int a; double b, c;
	
	for(int tagh_counter=1; tagh_counter<=274; tagh_counter++) {
		inf_tagh >> a >> b >> c;
		tagh_flux[tagh_counter-1]  = b;
		tagh_fluxE[tagh_counter-1] = c;
	}
	for(int tagm_counter=1; tagm_counter<=102; tagm_counter++) {
		inf_tagm >> a >> b >> c;
		tagm_flux[tagm_counter-1]  = b;
		tagm_fluxE[tagm_counter-1] = c;
	}
	
	inf_tagh.close();
	inf_tagm.close();
	
	// check if empty target filenames exist:
	
	ifstream inf_tagh_empty(empty_target_tagh_flux_fname.Data());
	ifstream inf_tagm_empty(empty_target_tagm_flux_fname.Data());
	
	if(!inf_tagh_empty.good() || !inf_tagm_empty.good()) {
		for(int tagh_counter=1; tagh_counter<=274; tagh_counter++) {
			tagh_flux_empty[tagh_counter-1] = 0.;
			tagh_fluxE_empty[tagh_counter-1] = 0.;
		}
		for(int tagm_counter=1; tagm_counter<=102; tagm_counter++) {
			tagm_flux_empty[tagm_counter-1] = 0.;
			tagm_fluxE_empty[tagm_counter-1] = 0.;
		}
		return 0;
	}
	
	for(int tagh_counter=1; tagh_counter<=274; tagh_counter++) {
		inf_tagh_empty >> a >> b >> c;
		tagh_flux_empty[tagh_counter-1]  = b;
		tagh_fluxE_empty[tagh_counter-1] = c;
	}
	for(int tagm_counter=1; tagm_counter<=102; tagm_counter++) {
		inf_tagm_empty >> a >> b >> c;
		tagm_flux_empty[tagm_counter-1]  = b;
		tagm_fluxE_empty[tagm_counter-1] = c;
	}
	
	return 0;
}

/*************************************************************************/
// Get NLO Theory Calculation from NIST and from primex_compton generator:

TString nist_theory_fname = "/home/andrsmit/root_macros/compton/nist_cs.dat";

TString theory_cs_pathName;

TF1 *f_theory, *f_nist;

TGraph  *g_theory;
TGraph  *g_theory_min, *g_theory_max, *g_theory_band;
TGraph  *g_dev_min,    *g_dev_max,    *g_dev_band;
TGraph  *g_nist;

TCanvas *c_theory;

bool USE_NIST_CALCULATION = false;
bool DRAW_THEORY          = true;

int get_theory_calc() {
	
	// check if filename containing NIST data exists:
	ifstream inf_nist(nist_theory_fname);
	if(!inf_nist.good()) {
		return 1;
	}
	
	// read in cross section from NIST:
	
	double nist_en[22], nist_cs[22];
	for(int ien=0; ien<22; ien++) {
		double aa, bb;
		inf_nist >> aa >> bb;
		nist_en[ien] = 1.e-3*aa;
		nist_cs[ien] = 1.e3*bb/4.;
	}
	g_nist = new TGraph(22, nist_en, nist_cs);
	g_nist->SetLineColor(kRed);
	g_nist->SetMarkerColor(kRed);
	g_nist->SetMarkerStyle(23);
	g_nist->SetMarkerSize(0.7);
	
	int n_cs_files = 0;
	vector<double> gen_enVec, gen_csVec;
	
	char fname[256];
	for(int tagh_counter=1; tagh_counter<=274; tagh_counter++) {
		
		sprintf(fname, "%s/tagh_%03d.dat", theory_cs_pathName.Data(), tagh_counter);
		if(gSystem->AccessPathName(fname)) {
			continue;
		} else {
			
			double aa, bb, cc, dd, ee;
			double loc_cs = 0.;
			
			ifstream locInf(fname);
			locInf >> aa >> bb >> cc >> dd >> ee;
			locInf >> aa >> bb >> cc >> dd >> ee;
			loc_cs += cc;
			locInf >> aa >> bb >> cc >> dd >> ee;
			loc_cs += cc;
			locInf.close();
			
			gen_enVec.push_back(0.5*(aa+bb));
			gen_csVec.push_back(loc_cs);
			n_cs_files++;
		}
	}
	for(int tagm_counter=1; tagm_counter<=102; tagm_counter++) {
		
		sprintf(fname, "%s/tagm_%03d.dat", theory_cs_pathName.Data(), tagm_counter);
		if(gSystem->AccessPathName(fname)) {
			continue;
		} else {
			
			double aa, bb, cc, dd, ee;
			double loc_cs = 0.;
			
			ifstream locInf(fname);
			locInf >> aa >> bb >> cc >> dd >> ee;
			locInf >> aa >> bb >> cc >> dd >> ee;
			loc_cs += cc;
			locInf >> aa >> bb >> cc >> dd >> ee;
			loc_cs += cc;
			locInf.close();
			
			gen_enVec.push_back(0.5*(aa+bb));
			gen_csVec.push_back(loc_cs);
			n_cs_files++;
		}
	}
	
	double *gen_energy = new double[n_cs_files];
	double *gen_cs     = new double[n_cs_files];
	
	for(int i=0; i < n_cs_files; i++) {
		gen_energy[i] = gen_enVec[i];
		gen_cs[i]     = gen_csVec[i];
	}
	
	TGraph *g_theory = new TGraph(n_cs_files, gen_energy, gen_cs);
	g_theory->SetLineColor(kRed); g_theory->SetLineWidth(2);
	g_theory->SetTitle("Generated Compton Cross Section (Born+SV+DH)");
	g_theory->GetXaxis()->SetTitle("Photon Beam Energy [GeV]");
	g_theory->GetXaxis()->SetTitleSize(0.045);
	g_theory->GetXaxis()->SetTitleOffset(1.0);
	g_theory->GetXaxis()->SetLabelSize(0.04);
	g_theory->GetYaxis()->SetTitle("#sigma [mb/electron]");
	g_theory->GetYaxis()->SetTitleSize(0.055);
	g_theory->GetYaxis()->SetTitleOffset(0.8);
	g_theory->GetYaxis()->SetLabelSize(0.05);
	g_theory->SetMarkerStyle(8);
	g_theory->SetMarkerSize(0.7);
	g_theory->SetMarkerColor(kBlack);
	
	f_theory  = new TF1("f_theory", "pol5", 3., 12.);
	g_theory->Fit(f_theory, "R0Q");
	f_theory->SetLineColor(kRed);
	f_theory->SetLineStyle(9);
	
	f_nist = new TF1("f_nist", "pol5", 3., 12.);
	g_nist->Fit(f_nist, "R0Q");
	f_nist->SetLineColor(kBlack);
	f_nist->SetLineStyle(1);
	
	/*
	if(USE_NIST_CALCULATION) {
		g_nist->Fit(f_theory, "R0Q");
	}
	*/
	
	double *gen_cs_min = new double[n_cs_files];
	double *gen_cs_max = new double[n_cs_files];
	
	for(int i=0; i < n_cs_files; i++) {
		gen_cs_min[i] = 0.95*f_theory->Eval(gen_energy[i]);
		gen_cs_max[i] = 1.05*f_theory->Eval(gen_energy[i]);
	}
	
	g_theory_min  = new TGraph(n_cs_files, gen_energy, gen_cs_min);
	g_theory_max  = new TGraph(n_cs_files, gen_energy, gen_cs_max);
	g_theory_band = new TGraph(2*n_cs_files);
	for(int i=0; i<n_cs_files; i++) {
		g_theory_band->SetPoint(i,gen_energy[i],gen_cs_max[i]);
		g_theory_band->SetPoint(n_cs_files+i,
			gen_energy[n_cs_files-i-1], gen_cs_min[n_cs_files-i-1]);
	}
	g_theory_band->SetFillStyle(3013);
	g_theory_band->SetFillColor(kCyan);
	g_theory_band->SetLineWidth(2);
	g_theory_band->SetLineColor(kCyan);
	g_theory_min->SetLineWidth(2);
	g_theory_min->SetLineColor(kCyan);
	g_theory_max->SetLineWidth(2);
	g_theory_max->SetLineColor(kCyan);
	
	double *gen_dev_min = new double[n_cs_files];
	double *gen_dev_max = new double[n_cs_files];
	
	for(int i=0; i < n_cs_files; i++) {
		gen_dev_min[i] = -5.0;
		gen_dev_max[i] =  5.0;
	}
	
	g_dev_min  = new TGraph(n_cs_files, gen_energy, gen_dev_min);
	g_dev_max  = new TGraph(n_cs_files, gen_energy, gen_dev_max);
	g_dev_band = new TGraph(2*n_cs_files);
	for(int i=0; i<n_cs_files; i++) {
		g_dev_band->SetPoint(i,gen_energy[i],gen_dev_max[i]);
		g_dev_band->SetPoint(n_cs_files+i,
			gen_energy[n_cs_files-i-1], gen_dev_min[n_cs_files-i-1]);
	}
	g_dev_band->SetFillStyle(3013);
	g_dev_band->SetFillColor(kCyan);
	g_dev_band->SetLineColor(kCyan);
	g_dev_band->SetLineWidth(2);
	g_dev_min->SetLineColor(kCyan);
	g_dev_min->SetLineWidth(2);
	g_dev_max->SetLineColor(kCyan);
	g_dev_max->SetLineWidth(2);
	
	if(DRAW_THEORY) {
		
		c_theory = new TCanvas("c_theory", "c_theory", 1000., 500.);
		c_theory->SetTickx(); c_theory->SetTicky();
		g_theory->Draw("AP");
		f_theory->Draw("same");
		g_nist->Draw("PC same");
		c_theory->Update();
	}
	
	return 0;
}

/*************************************************************************/
// Get e+e- pair production cross section (from NIST):

TString pair_cs_fname;

TF1 *f_pair_cs, *f_triplet_cs;

TCanvas *c_pair_cs;

void get_pair_cs()
{
	// read in cross section from NIST:
	
	vector<double> photon_energy_vec, pair_cs_vec, triplet_cs_vec;
	ifstream inf(pair_cs_fname.Data());
	double a, b, c;
	while(inf >> a >> b >> c) {
		photon_energy_vec.push_back(a);
		pair_cs_vec.push_back(b);
		triplet_cs_vec.push_back(c);
	}
	inf.close();
	
	int n_ens = (int)photon_energy_vec.size();
	double *photon_energy = new double[n_ens];
	double *pair_cs       = new double[n_ens];
	double *triplet_cs    = new double[n_ens];
	for( int i=0; i<n_ens; i++ ) {
		photon_energy[i] = photon_energy_vec[i] / 1.e3;
		pair_cs[i]       = 1.e3 * pair_cs_vec[i];
		triplet_cs[i]    = 1.e3 * triplet_cs_vec[i];
	}
	
	TGraph *g_pair_cs = new TGraph(n_ens, photon_energy, pair_cs);
	g_pair_cs->SetTitle("e^{+}e^{-} Pair Production CS");
	g_pair_cs->GetXaxis()->SetTitle("E_{#gamma} [GeV]");
	g_pair_cs->GetYaxis()->SetTitle("#sigma [mb / atom]");
	g_pair_cs->SetMarkerColor(kBlue);
	g_pair_cs->SetMarkerStyle(8);
	
	f_pair_cs = new TF1("f_pair_cs", "pol4", 3.0, 12.0);
	f_pair_cs->SetLineColor(kBlue);
	f_pair_cs->SetLineStyle(2);
	g_pair_cs->Fit(f_pair_cs, "R0Q");
	
	TGraph *g_triplet_cs = new TGraph(n_ens, photon_energy, triplet_cs);
	g_triplet_cs->SetTitle("e^{+}e^{-} Triplet Production CS");
	g_triplet_cs->GetXaxis()->SetTitle("E_{#gamma} [GeV]");
	g_triplet_cs->GetYaxis()->SetTitle("#sigma [mb / atom]");
	g_triplet_cs->SetMarkerColor(kBlue);
	g_triplet_cs->SetMarkerStyle(8);
	
	f_triplet_cs = new TF1("f_triplet_cs", "pol4", 3.0, 12.0);
	f_triplet_cs->SetLineColor(kBlue);
	f_triplet_cs->SetLineStyle(2);
	g_triplet_cs->Fit(f_triplet_cs, "R0Q");
	
	if(DRAW_THEORY) {
	
		c_pair_cs = new TCanvas("c_pair_cs", "c_pair_cs", 800, 800);
		c_pair_cs->Divide(1,2);
		
		TPad *p_pair_cs = (TPad*)c_pair_cs->cd(1);
		p_pair_cs->cd();
		p_pair_cs->SetTickx(); p_pair_cs->SetTicky();
		g_pair_cs->Draw("AP");
		f_pair_cs->Draw("same");
		
		TPad *p_triplet_cs = (TPad*)c_pair_cs->cd(2);
		p_triplet_cs->cd();
		p_triplet_cs->SetTickx(); p_triplet_cs->SetTicky();
		g_triplet_cs->Draw("AP");
		f_triplet_cs->Draw("same");
		
		c_pair_cs->Update();
	}
	
	return;
}

/*************************************************************************/
// Get target parameters:

double n_e = 0.0;
double n_Z = 0.0;
double n_A = 0.0;

double mb = 1.e-27;

double f_abs = 1.0;

void get_target_parameters()
{
	if(IS_BE_TARGET) {
		
		// PrimEx-eta Be-9 Target:
		
		n_Z   = 4.0;
		n_A   = 9.0;
		f_abs = 0.981;
		
		double target_thickness = 1.77546;
		double avo_num          = 6.02214076e+23;
		
		double frac_be2c = 1.e-2 * 0.021;
		double frac_fe   = 1.e-2 * 0.105;
		double frac_al   = 1.e-2 * 0.0068;
		double frac_si   = 1.e-2 * 0.0058;
		double frac_mg   = 1.e-2 * 0.0021;
		double frac_mn   = 1.e-2 * 0.0015;
		double frac_be   = 1.0 - (frac_be2c + frac_fe + frac_al + frac_si + frac_mg + frac_mn);
		
		double ne_be2c   = 14.;
		double ne_fe     = 26.;
		double ne_al     = 13.;
		double ne_si     = 14.;
		double ne_mg     = 12.;
		double ne_mn     = 25.;
		double ne_be     =  4.;
		
		double density_be2c = 1.90;
		double density_fe   = 7.874;
		double density_al   = 2.70;
		double density_si   = 2.329;
		double density_mg   = 1.738;
		double density_mn   = 7.26;
		double density_be   = 1.848;
		
		double weight_be2c  = 30.035;
		double weight_fe    = 55.84;
		double weight_al    = 26.981538;
		double weight_si    = 28.085;
		double weight_mg    = 24.305;
		double weight_mn    = 54.938;
		double weight_be    = 9.012183;
		
		ne_be2c *= (density_be2c * target_thickness * (1./weight_be2c) * avo_num);
		ne_fe   *= (density_fe   * target_thickness * (1./weight_fe  ) * avo_num);
		ne_al   *= (density_al   * target_thickness * (1./weight_al  ) * avo_num);
		ne_si   *= (density_si   * target_thickness * (1./weight_si  ) * avo_num);
		ne_mg   *= (density_mg   * target_thickness * (1./weight_mg  ) * avo_num);
		ne_mn   *= (density_mn   * target_thickness * (1./weight_mn  ) * avo_num);
		ne_be   *= (density_be   * target_thickness * (1./weight_be  ) * avo_num);
		
		n_e   = frac_be2c * ne_be2c
			+ frac_fe   * ne_fe
			+ frac_al   * ne_al
			+ frac_si   * ne_si
			+ frac_mg   * ne_mg
			+ frac_mn   * ne_mn
			+ frac_be   * ne_be;
		
	} else {
		
		// PrimEx-eta He-4 Target:
		
		n_Z   = 2.0;
		n_A   = 4.0;
		//f_abs = 0.9908;
		f_abs = 0.9856;
		n_e   = 1.0816e+24; // p = 0.1217 g/cm3, t = 29.535cm
		
		// correction for cold residual gas:
		
		//n_e *= 0.9817;
		
	}
	
	return;
}

#endif
