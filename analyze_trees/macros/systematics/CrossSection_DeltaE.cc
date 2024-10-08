#include "/work/halld/home/andrsmit/primex_compton_analysis/include/compton_cs.h"
#include "/work/halld/home/andrsmit/primex_compton_analysis/include/compton_inputs.cc"
#include "/work/halld/home/andrsmit/primex_compton_analysis/include/compton_yield.cc"
#include "/work/halld/home/andrsmit/primex_compton_analysis/include/compton_acc.cc"
#include "/work/halld/home/andrsmit/primex_compton_analysis/include/plot_results.cc"

//----------   Function Declarations   ----------//

void print_status();

//-----------------------------------------------//

void CrossSection_DeltaE(double loc_cut=5.0)
{
	//==============================================================================//
	// TWO ADJUSTABLE SWITCHES TO CHOOSE DATA SET:
	
	IS_BE_TARGET = true;
	BEAM_CURRENT = 200;
	
	IS_FIELD_OFF = true;
	
	FIT_EMPTY    = false;
	FIT_TRIPLET  = false;
	
	DELTA_K_FIT_SIGMA = 6.0;
	
	for(int i=0; i<274; i++) tagh_flux_unc[i] = 0.;
	for(int i=0; i<102; i++) tagm_flux_unc[i] = 0.;
	
	ifstream inf_flux_unc("tagh_flux_syst.txt");
	int aa; double bb, cc, dd;
	for(int i=0; i<175; i++) {
		inf_flux_unc >> aa >> bb >> cc >> dd;
		tagh_flux_unc[aa-1] = cc;
	}
	inf_flux_unc.close();
	inf_flux_unc.open("tagm_flux_syst.txt");
	for(int i=0; i<102; i++) {
		inf_flux_unc >> aa >> bb >> cc >> dd;
		tagm_flux_unc[aa-1] = cc;
	}
	inf_flux_unc.close();
	
	//==============================================================================//
	
	gStyle->SetOptStat(0); gStyle->SetOptFit(1);
	gStyle->SetEndErrorSize(4);
	
	bad_counters_tagh.clear();
	bad_counters_tagm.clear();
	//bad_counters_tagh = {};//59,  73, 124, 125, 126, 127, 179, 180, 181, 186, 187, 220, 221};
	//bad_counters_tagm = {};// 21,  63,  67,  75,  88, 100, 101, 102};
	
	if(IS_BE_TARGET) {
		bad_counters_tagh = { 54, 73, 125, 186};
		bad_counters_tagm = { 63, 67, 100, 101, 102};
	} else {
		bad_counters_tagh = { 54, 73, 125, 185, 186, 187, 203, 206, 208};
		bad_counters_tagm = { 29, 30, 63, 67, 100, 101, 102};
	}
	const char loc_pathname[256] = "/work/halld/home/andrsmit/primex_compton_analysis";
	
	//----------   Initialize   ----------//
	
	// Adjustable switches for monitoring:
	
	DEBUG_FITS       = false;
	
	DRAW_THEORY      = false;
	DRAW_ACCEPTANCE  = true;
	DRAW_FITS_TAGH   = true;
	DRAW_FITS_TAGM   = true;
	USE_F_ACC        = false;
	
	// adjust bin size. Default: 2000 bins in DeltaK histogram (8MeV per bin)
	
	rebins   = 2;
	bin_size = (16.0/1000.0)*(double)rebins;
	n_mev    = bin_size * 1.e3;
	
	bool APPLY_CORRECTION = false;
	
	//==============================================================================//
	
	if(IS_BE_TARGET) TARGET_STR = "Be";
	else             TARGET_STR = "He";
	
	if(IS_FIELD_OFF) FIELD_STR = "FIELDOFF";
	else             FIELD_STR = "FIELDON";
	
	//----------------------------------------//
	// Photon beam energy:
	
	if(IS_BE_TARGET) {
		endpoint_energy_calib = 11.6061;
		endpoint_energy       = 11.608;
	} else {
		endpoint_energy_calib = 11.1689;
		if(BEAM_CURRENT==200) {
			endpoint_energy = 11.1666;
		} else if(BEAM_CURRENT==100) {
			endpoint_energy = 11.1664;
		} else {
			endpoint_energy = 11.1671;
		}
	}
	
	tagh_xscale_fname = Form("%s/photon_flux/phase1/primex_tagh.txt", loc_pathname);
	tagm_xscale_fname = Form("%s/photon_flux/phase1/primex_tagm.txt", loc_pathname);
	
	//----------------------------------------//
	// Photon flux:
	
	tagh_flux_fname = Form("%s/photon_flux/phase1/%s_%03dnA_%s_flux_tagh.txt", 
		loc_pathname, TARGET_STR.Data(), BEAM_CURRENT, FIELD_STR.Data());
	tagm_flux_fname = Form("%s/photon_flux/phase1/%s_%03dnA_%s_flux_tagm.txt", 
		loc_pathname, TARGET_STR.Data(), BEAM_CURRENT, FIELD_STR.Data());
	/*
	empty_target_tagh_flux_fname = Form("%s/photon_flux/phase1/%s_empty_%s_flux_tagh.txt", 
		loc_pathname, TARGET_STR.Data(), FIELD_STR.Data());
	empty_target_tagm_flux_fname = Form("%s/photon_flux/phase1/%s_empty_%s_flux_tagm.txt", 
		loc_pathname, TARGET_STR.Data(), FIELD_STR.Data());
	*/
	empty_target_tagh_flux_fname = Form("%s/photon_flux/phase1/Be_empty_%s_flux_tagh.txt", 
		loc_pathname, FIELD_STR.Data());
	empty_target_tagm_flux_fname = Form("%s/photon_flux/phase1/Be_empty_%s_flux_tagm.txt", 
		loc_pathname,  FIELD_STR.Data());
	
	//----------------------------------------//
	// Compton & e+e- Pair Production Theory:
	
	theory_cs_pathName = Form("%s/compton_mc/genDir/Run061321/genCS", loc_pathname);
	pair_cs_fname      = Form("%s/photon_absorption/%s_pair_cs.dat",  loc_pathname, TARGET_STR.Data());
	
	//----------------------------------------//
	// Data:
	
	root_fname              = Form("%s/analyze_trees/phase1/analyze_data/rootFiles/systematics/%s_%03dnA_%s.root", 
		loc_pathname, TARGET_STR.Data(), BEAM_CURRENT, FIELD_STR.Data());
	
	//empty_target_root_fname = Form("%s/analyze_trees/phase1/analyze_data/rootFiles/systematics/%s_empty_%s.root", 
	//	loc_pathname, TARGET_STR.Data(), FIELD_STR.Data());
	empty_target_root_fname = Form("%s/analyze_trees/phase1/analyze_data/rootFiles/systematics/Be_empty_%s.root", 
		loc_pathname, FIELD_STR.Data());
	
	int loc_cut_int = 10.0 * loc_cut;
	
	hname_tagh = Form("DeltaE/deltaK_tagh_%03dsigE", loc_cut_int);
	hname_tagm = Form("DeltaE/deltaK_tagm_%03dsigE", loc_cut_int);
	
	//----------------------------------------//
	// Simulation:
	
	int mc_runNumber = 61321;
	if(!IS_BE_TARGET) mc_runNumber = 61866;
	
	// Compton MC:
	comp_mc_dir     = Form("%s/analyze_trees/phase1/analyze_mc/rootFiles/Run%06d/compton/%03dnA_systematics", 
		loc_pathname, mc_runNumber, BEAM_CURRENT);
	hname_tagh_comp = hname_tagh;
	hname_tagm_comp = hname_tagm;
	
	// e+e- Pair Production MC:
	pair_mc_dir     = Form("%s/analyze_trees/phase1/analyze_mc/rootFiles/Run061321/pair/systematics", loc_pathname);
	hname_tagh_pair = hname_tagh;
	hname_tagm_pair = hname_tagm;
	
	//==============================================================================//
	
	// Check that all files exist:
	
	if(gSystem->AccessPathName(root_fname.Data())) {
		cout << "Specified ROOT file does not exist, check file name.\n\n";
		return;
	}
	if(gSystem->AccessPathName(empty_target_root_fname.Data())) {
		cout << "Specified empty target ROOT file does not exist, check file name.\n\n" << endl;
		return;
	}
	if(gSystem->AccessPathName(tagh_flux_fname.Data()) || gSystem->AccessPathName(tagm_flux_fname.Data())) {
		cout << "Specified flux files do not exist, check file name.\n\n" << endl;
		return;
	}
	if(gSystem->AccessPathName(empty_target_tagh_flux_fname.Data()) 
		|| gSystem->AccessPathName(empty_target_tagh_flux_fname.Data())) {
		cout << "Specified empty target flux files do not exist, check file name.\n\n" << endl;
		return;
	}
	
	//==============================================================================//
	
	print_status();
	
	//==============================================================================//
	
	get_target_parameters();
	get_flux();
	get_counter_energies();
	get_theory_calc();
	get_pair_cs();
	get_compton_acc();
	//return;
	//==============================================================================//
	
	vector<int> tagh_vec, tagm_vec;
	tagh_vec.clear();
	tagm_vec.clear();
	
	get_compton_yield(tagh_vec, tagm_vec);
	
	/*
	TF1 *f_correction = new TF1("f_correction", "pol7", 1., 20.);
	f_correction->SetParameters(
		9.196611e-01, -5.843140e-02, 7.641279e-02, -2.962729e-02, 5.772550e-03, -6.164195e-04, 3.448740e-05, -7.920409e-07
	);
	for(int tagh_counter=1; tagh_counter<=274; tagh_counter++) {
		double loc_yield = tagh_yield[tagh_counter-1];
		if(loc_yield <= 0. || loc_cut>10.) continue;
		if(APPLY_CORRECTION) tagh_yield[tagh_counter-1] = loc_yield / f_correction->Eval(loc_cut);
	}
	for(int tagm_counter=1; tagm_counter<=102; tagm_counter++) {
		double loc_yield = tagm_yield[tagm_counter-1];
		if(loc_yield <= 0. || loc_cut>10.) continue;
		if(APPLY_CORRECTION) tagm_yield[tagm_counter-1] = loc_yield / f_correction->Eval(loc_cut);
	}
	*/
	
	TF1 *f_correction = new TF1("f_correction", "expo", 1.0, 12.0);
	f_correction->SetParameters(-1.68865, 0.164773);
	for(int tagh_counter=1; tagh_counter<=274; tagh_counter++) {
		double loc_yield = tagh_yield[tagh_counter-1];
		if(loc_yield <= 0.) continue;
		if(APPLY_CORRECTION) tagh_yield[tagh_counter-1] = loc_yield / (1.0-1.e-2*f_correction->Eval(tagh_en[tagh_counter-1]));
	}
	for(int tagm_counter=1; tagm_counter<=102; tagm_counter++) {
		double loc_yield = tagm_yield[tagm_counter-1];
		if(loc_yield <= 0.) continue;
		if(APPLY_CORRECTION) tagm_yield[tagm_counter-1] = loc_yield / (1.0-1.e-2*f_correction->Eval(tagm_en[tagm_counter-1]));
	}
	//==============================================================================//
	
	// plot results:
	
	plot_cs(tagh_vec, tagm_vec);
	plot_fit_fractions(tagh_vec, tagm_vec);
	plot_chi2(tagh_vec, tagm_vec);
	
	//==============================================================================//
	
	// Write out Data:
	
	char buf[256];
	
	ofstream outf_tagh(Form("%s_%03dnA_tagh_cross_section.txt", TARGET_STR.Data(), BEAM_CURRENT));
	for(int tagh_counter=1; tagh_counter<=274; tagh_counter++) {
		
		double loc_acc, loc_accE, loc_cs, loc_csE;
		if(tagh_yield[tagh_counter-1] <= 0.) {
			loc_acc  = 0.;
			loc_accE = 0.;
			loc_cs   = 0.;
			loc_csE  = 0.;
		} else {
			if(tagh_acc[tagh_counter-1]<=0.) {
				loc_acc  = f_acc->Eval(tagh_en[tagh_counter-1]);
				loc_accE = sqrt(9.e5*loc_acc*(1.0-loc_acc))/9.e5;
			} else {
				loc_acc  = tagh_acc[tagh_counter-1];
				loc_accE = tagh_accE[tagh_counter-1];
			}
			loc_cs  = tagh_cs[tagh_counter-1];
			loc_csE = tagh_csE[tagh_counter-1];
		}
		sprintf(buf, "%e   %e   %e   %e\n", loc_cs, loc_csE, loc_acc, loc_accE);
		outf_tagh << buf;
	}
	outf_tagh.close();
	
	ofstream outf_tagm(Form("%s_%03dnA_tagm_cross_section.txt", TARGET_STR.Data(), BEAM_CURRENT));
	for(int tagm_counter=1; tagm_counter<=102; tagm_counter++) {
		
		double loc_acc, loc_accE, loc_cs, loc_csE;
		if(tagm_yield[tagm_counter-1] <= 0.) {
			loc_acc  = 0.;
			loc_accE = 0.;
			loc_cs   = 0.;
			loc_csE  = 0.;
		} else {
			if(tagm_acc[tagm_counter-1]<=0.) {
				loc_acc  = f_acc->Eval(tagm_en[tagm_counter-1]);
				loc_accE = sqrt(9.e5*loc_acc*(1.0-loc_acc))/9.e5;
			} else {
				loc_acc  = tagm_acc[tagm_counter-1];
				loc_accE = tagm_accE[tagm_counter-1];
			}
			loc_cs  = tagm_cs[tagm_counter-1];
			loc_csE = tagm_csE[tagm_counter-1];
		}
		sprintf(buf, "%e   %e   %e   %e\n", loc_cs, loc_csE, loc_acc, loc_accE);
		outf_tagm << buf;
	}
	outf_tagm.close();
	
	//---------------------------------------------------------------------------------//
	// To write out each component of the cross section calculation for a latex table:
	/*
	ofstream outf_tagh(Form("%s_%03dnA_tagh_cross_section.txt", TARGET_STR.Data(), BEAM_CURRENT));
	for(int tagh_counter=1; tagh_counter<=274; tagh_counter++) {
		if(tagh_yield[tagh_counter-1] <= 0.) continue;
		double loc_acc = tagh_acc[tagh_counter-1];
		if(loc_acc <= 0.) {
			loc_acc = f_acc->Eval(tagh_en[tagh_counter-1]);
		}
		sprintf(buf, "TAGH %03d & %.3f & %.3f & %.4f & %.3f & %.4f & %.2f & \\textbf{%.4f} \\\\\n", 
			tagh_counter, tagh_en[tagh_counter-1], n_e*1.e-23, tagh_flux[tagh_counter-1]*1.e-10, f_abs, loc_acc, 
			tagh_yield[tagh_counter-1], tagh_cs[tagh_counter-1]);
		outf_tagh << buf;
	}
	outf_tagh.close();
	
	ofstream outf_tagm(Form("%s_%03dnA_tagm_cross_section.txt", TARGET_STR.Data(), BEAM_CURRENT));
	for(int tagm_counter=1; tagm_counter<=102; tagm_counter++) {
		if(tagm_yield[tagm_counter-1] <= 0.) continue;
		double loc_acc = tagm_acc[tagm_counter-1];
		if(loc_acc <= 0.) {
			loc_acc = f_acc->Eval(tagm_en[tagm_counter-1]);
		}
		sprintf(buf, "TAGM %03d & %.3f & %.3f & %.4f & %.3f & %.4f & %.2f & \\textbf{%.4f} \\\\\n", 
			tagm_counter, tagm_en[tagm_counter-1], n_e*1.e-23, tagm_flux[tagm_counter-1]*1.e-10, f_abs, loc_acc, 
			tagm_yield[tagm_counter-1], tagm_cs[tagm_counter-1]);
		outf_tagm << buf;
	}
	outf_tagm.close();
	*/
	//---------------------------------------------------------------------------------//
	// To write out the systematic uncertainty from empty target background:
	/*
	ofstream outf(Form("%s_%03dnA_empty_unc_tagh.txt", TARGET_STR.Data(), BEAM_CURRENT));
	for(int tagh_counter=1; tagh_counter<=274; tagh_counter++) {
		double loc_unc = 0.;
		
		if(tagh_yield[tagh_counter-1] > 0.) {
			double loc_frac  = tagh_fit_fraction["empty"][tagh_counter-1].first;
			loc_unc = 1.e2 * (0.05*loc_frac) / (1.0-loc_frac);
		}
		sprintf(buf, "%e\n", loc_unc);
		outf << buf;
	}
	outf.close();
	
	outf.open(Form("%s_%03dnA_empty_unc_tagm.txt", TARGET_STR.Data(), BEAM_CURRENT));
	for(int tagm_counter=1; tagm_counter<=102; tagm_counter++) {
		double loc_unc = 0.;
		
		if(tagm_yield[tagm_counter-1] > 0.) {
			double loc_frac  = tagm_fit_fraction["empty"][tagm_counter-1].first;
			loc_unc = 1.e2 * (0.05*loc_frac) / (1.0-loc_frac);
		}
		sprintf(buf, "%e\n", loc_unc);
		outf << buf;
	}
	outf.close();
	*/
	//---------------------------------------------------------------------------------//
	// To write out the systematic uncertainty from e+e- background:
	/*
	ofstream outf(Form("%s_%03dnA_pair_unc_tagh.txt", TARGET_STR.Data(), BEAM_CURRENT));
	for(int tagh_counter=1; tagh_counter<=274; tagh_counter++) {
		double loc_unc = 0.;
		
		if(tagh_yield[tagh_counter-1] > 0.) {
			double loc_frac  = tagh_fit_fraction["pair"][tagh_counter-1].first;
			double loc_fracE = tagh_fit_fraction["pair"][tagh_counter-1].second;
			
			
			double loc_unc1 = 1.e2 * (loc_fracE) / (1.0-loc_frac);
			double loc_unc2 = 0.4;
			loc_unc = sqrt(pow(loc_unc1,2.0) + pow(loc_unc2,2.0));
		}
		sprintf(buf, "%e\n", loc_unc);
		outf << buf;
	}
	outf.close();
	
	outf.open(Form("%s_%03dnA_pair_unc_tagm.txt", TARGET_STR.Data(), BEAM_CURRENT));
	for(int tagm_counter=1; tagm_counter<=102; tagm_counter++) {
		double loc_unc = 0.;
		
		if(tagm_yield[tagm_counter-1] > 0.) {
			double loc_frac  = tagm_fit_fraction["pair"][tagm_counter-1].first;
			double loc_fracE = tagm_fit_fraction["pair"][tagm_counter-1].second;
			
			
			double loc_unc1 = 1.e2 * (loc_fracE) / (1.0-loc_frac);
			double loc_unc2 = 0.4;
			loc_unc = sqrt(pow(loc_unc1,2.0) + pow(loc_unc2,2.0));
		}
		sprintf(buf, "%e\n", loc_unc);
		outf << buf;
	}
	outf.close();
	*/
	//---------------------------------------------------------------------------------//
	// To write out the cross section for systematic uncertainty analysis:
	/*
	char out_fname[256];
	sprintf(out_fname, "DeltaE/%s_%03dnA/tagh_cross_section_%03dsigE.txt", TARGET_STR.Data(), BEAM_CURRENT, loc_cut_int);
	ofstream outf_tagh(out_fname);
	for(int tagh_counter=1; tagh_counter<=274; tagh_counter++) {
		if(tagh_yield[tagh_counter-1] <= 0.) {
			sprintf(buf, "%03d   %.4f	%.5f   %.5f   %f   %f	%f   %f",
				tagh_counter, tagh_en[tagh_counter-1], 0., 0., 0., 0., 0., 0.);
		} else {
			sprintf(buf, "%03d   %.4f   %.5f   %.5f   %f   %f   %f   %f", 
				tagh_counter, tagh_en[tagh_counter-1], 
				tagh_cs[tagh_counter-1], tagh_csE[tagh_counter-1], 
				tagh_yield[tagh_counter-1], tagh_yieldE[tagh_counter-1], 
				tagh_acc[tagh_counter-1], tagh_accE[tagh_counter-1]);
		}
		outf_tagh << buf << "\n";
	}
	outf_tagh.close();
	
	sprintf(out_fname, "DeltaE/%s_%03dnA/tagm_cross_section_%03dsigE.txt", TARGET_STR.Data(), BEAM_CURRENT, loc_cut_int);
	ofstream outf_tagm(out_fname);
	for(int tagm_counter=1; tagm_counter<=102; tagm_counter++) {
		if(tagm_yield[tagm_counter-1] <= 0.) {
			sprintf(buf, "%03d   %.4f	%.5f   %.5f   %f   %f	%f   %f",
				tagm_counter, tagm_en[tagm_counter-1], 0., 0., 0., 0., 0., 0.);
		} else {
			sprintf(buf, "%03d   %.4f   %.5f   %.5f   %f   %f   %f   %f", 
				tagm_counter, tagm_en[tagm_counter-1], 
				tagm_cs[tagm_counter-1], tagm_csE[tagm_counter-1], 
				tagm_yield[tagm_counter-1], tagm_yieldE[tagm_counter-1], 
				tagm_acc[tagm_counter-1], tagm_accE[tagm_counter-1]);
		}
		outf_tagm << buf << "\n";
	}
	outf_tagm.close();
	*/
	return;
}

void print_status() {
	
	cout << "\n\n";
	cout << "***************************************************************************";
	cout << "*******************" << endl;
	printf("Calculating total Compton Scattering Cross Section for %03dnA %s target (%s):\n\n", 
		BEAM_CURRENT, TARGET_STR.Data(), FIELD_STR.Data());
	printf(" Filled target ROOT File: \n");
	printf("  %s\n", root_fname.Data());
	printf(" Empty target ROOT File: \n");
	printf("  %s\n\n", empty_target_root_fname.Data());
	printf(" Filled target flux accessed from: \n");
	printf("  %s\n", tagh_flux_fname.Data());
	printf("  %s\n", tagm_flux_fname.Data());
	printf(" Empty target flux accessed from: \n");
	printf("  %s\n", empty_target_tagh_flux_fname.Data());
	printf("  %s\n\n", empty_target_tagm_flux_fname.Data());
	printf(" Compton MC directory: \n");
	printf("  %s\n", comp_mc_dir.Data());
	printf(" Theory calculation (primex): \n");
	printf("  %s\n", theory_cs_pathName.Data());
	printf(" Theory calculation (NIST): \n");
	printf("  %s\n\n", nist_theory_fname.Data());
	printf(" e+e- Pair MC directory: \n");
	printf("  %s\n", pair_mc_dir.Data());
	cout << "***************************************************************************";
	cout << "*******************" << endl;
	cout << "\n\n";
	
	return;
}
