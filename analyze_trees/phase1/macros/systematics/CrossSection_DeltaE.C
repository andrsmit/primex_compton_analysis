"TEST"

#include "/work/halld/home/andrsmit/primex_compton_analysis/include/compton_inputs.h"
#include "/work/halld/home/andrsmit/primex_compton_analysis/include/compton_acc.h"
#include "/work/halld/home/andrsmit/primex_compton_analysis/include/compton_yield.h"
#include "/work/halld/home/andrsmit/primex_compton_analysis/include/plot_results.h"

//----------   Function Declarations   ----------//

void print_status();

//-----------------------------------------------//

void CrossSection_DeltaE(double loc_sig=5.0)
{
	//==============================================================================//
	// TWO ADJUSTABLE SWITCHES TO CHOOSE DATA SET:
	
	IS_BE_TARGET = false;
	BEAM_CURRENT = 200;
	
	IS_FIELD_OFF = true;
	
	//==============================================================================//
	
	gStyle->SetOptStat(0); gStyle->SetOptFit(1);
	gStyle->SetEndErrorSize(4);
	
	bad_counters_tagh.clear();
	bad_counters_tagm.clear();
	bad_counters_tagh = {73, 125, 126, 127, 179, 180, 181, 186, 187, 220, 221};
	bad_counters_tagm = {0, 1, 2, 3, 4, 21, 63, 67, 98, 99, 100, 101, 102};
	
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
	
	rebins   =  5;
	bin_size = (16.0/1000.0)*(double)rebins;
	n_mev    = bin_size * 1.e3;
	
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
	
	empty_target_tagh_flux_fname = Form("%s/photon_flux/phase1/%s_empty_%s_flux_tagh.txt", 
			loc_pathname, TARGET_STR.Data(), FIELD_STR.Data());
	empty_target_tagm_flux_fname = Form("%s/photon_flux/phase1/%s_empty_%s_flux_tagm.txt", 
			loc_pathname, TARGET_STR.Data(), FIELD_STR.Data());
	
	//----------------------------------------//
	// Compton & e+e- Pair Production Theory:
	
	theory_cs_pathName = Form("%s/compton_mc/genDir/Run061321/genCS", loc_pathname);
	pair_cs_fname      = Form("%s/photon_absorption/%s_pair_cs.dat",  loc_pathname, TARGET_STR.Data());
	
	//----------------------------------------//
	// Data:
	
	root_fname              = Form("%s/analyze_trees/phase1/analyze_data/rootFiles/systematics/%s_%03dnA_%s.root", 
		loc_pathname, TARGET_STR.Data(), BEAM_CURRENT, FIELD_STR.Data());
	
	empty_target_root_fname = Form("%s/analyze_trees/phase1/analyze_data/rootFiles/systematics/%s_empty_%s.root", 
		loc_pathname, TARGET_STR.Data(), FIELD_STR.Data());
	
	int loc_sig_int = (int)(10.*loc_sig);
	
	hname_tagh = Form("DeltaE/deltaK_tagh_%03dsigE", loc_sig_int);
	hname_tagm = Form("DeltaE/deltaK_tagm_%03dsigE", loc_sig_int);
	
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
	
	//==============================================================================//
	
	vector<int> tagh_vec, tagm_vec;
	tagh_vec.clear();
	tagm_vec.clear();
	
	get_compton_yield(tagh_vec, tagm_vec);
	
	TF1 *f_correction = new TF1("f_correction", "pol7", 1., 20.);
	f_correction->SetParameters(
		9.196611e-01, -5.843140e-02, 7.641279e-02, -2.962729e-02, 5.772550e-03, -6.164195e-04, 3.448740e-05, -7.920409e-07
	);
	/*
	f_correction->SetParameters(
		0.929162, -0.0857315, 0.0953944, -0.0359486, 0.00694032, -0.000738923, 4.13216e-05, -9.49538e-07
	);
	f_correction->SetParameters(
		0.901548, -0.0831837, 0.0925593, -0.0348802, 0.00673406, -0.000716963, 4.00936e-05, -9.21319e-07
	);
	f_correction->SetParameters(
		0.924853, -0.0853339, 0.0949519, -0.0357819, 0.00690813, -0.000735496, 4.113e-05, -9.45134e-07
	);
	*/
	for(int tagh_counter=1; tagh_counter<=274; tagh_counter++) {
		double loc_yield = tagh_yield[tagh_counter-1];
		if(loc_yield <= 0. || loc_sig>10.) continue;
		//tagh_yield[tagh_counter-1] = loc_yield / f_correction->Eval(loc_sig);
	}
	for(int tagm_counter=1; tagm_counter<=102; tagm_counter++) {
		double loc_yield = tagm_yield[tagm_counter-1];
		if(loc_yield <= 0. || loc_sig>10.) continue;
		//tagm_yield[tagm_counter-1] = loc_yield / f_correction->Eval(loc_sig);
	}
	
	//==============================================================================//
	
	// plot results:
	
	plot_cs(tagh_vec, tagm_vec);
	plot_pair_fraction(tagh_vec, tagm_vec);
	plot_empty_fraction(tagh_vec, tagm_vec);
	plot_chi2(tagh_vec, tagm_vec);
	
	//==============================================================================//
	
	// Write out Data:
	
	char out_fname[256];
	sprintf(out_fname, "DeltaE/%s_%03dnA/tagh_cross_section_%03dsigE.txt", TARGET_STR.Data(), BEAM_CURRENT, loc_sig_int);
	ofstream outf_tagh(out_fname);
	for(int tagh_counter=1; tagh_counter<=274; tagh_counter++) {
		char buf[256];
		sprintf(buf, "%03d   %.4f   %.5f   %.5f   %f   %f   %f   %f   %f   %f", 
			tagh_counter, tagh_en[tagh_counter-1], 
			tagh_cs[tagh_counter-1], tagh_csE[tagh_counter-1], 
			tagh_yield[tagh_counter-1], tagh_yieldE[tagh_counter-1], 
			tagh_acc[tagh_counter-1], tagh_accE[tagh_counter-1], 
			tagh_pair_fraction[tagh_counter-1], tagh_pair_fractionE[tagh_counter-1]);
		outf_tagh << buf << "\n";
	}
	outf_tagh.close();
	
	sprintf(out_fname, "DeltaE/%s_%03dnA/tagm_cross_section_%03dsigE.txt", TARGET_STR.Data(), BEAM_CURRENT, loc_sig_int);
	ofstream outf_tagm(out_fname);
	for(int tagm_counter=1; tagm_counter<=102; tagm_counter++) {
		char buf[256];
		sprintf(buf, "%03d   %.4f   %.5f   %.5f   %f   %f   %f   %f   %f   %f", 
			tagm_counter, tagm_en[tagm_counter-1], 
			tagm_cs[tagm_counter-1], tagm_csE[tagm_counter-1], 
			tagm_yield[tagm_counter-1], tagm_yieldE[tagm_counter-1], 
			tagm_acc[tagm_counter-1], tagm_accE[tagm_counter-1], 
			tagm_pair_fraction[tagm_counter-1], tagm_pair_fractionE[tagm_counter-1]);
		outf_tagm << buf << "\n";
	}
	outf_tagm.close();
	
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
