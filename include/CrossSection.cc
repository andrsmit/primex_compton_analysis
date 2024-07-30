#include "/work/halld/home/andrsmit/primex_compton_analysis/include/compton_cs.h"
#include "/work/halld/home/andrsmit/primex_compton_analysis/include/compton_inputs.cc"
#include "/work/halld/home/andrsmit/primex_compton_analysis/include/compton_yield.cc"
#include "/work/halld/home/andrsmit/primex_compton_analysis/include/compton_acc.cc"
#include "/work/halld/home/andrsmit/primex_compton_analysis/include/plot_results.cc"

//----------   Function Declarations   ----------//

void print_status();

//-----------------------------------------------//

void CrossSection()
{
	//==============================================================================//
	// TWO ADJUSTABLE SWITCHES TO CHOOSE DATA SET:
	
	IS_BE_TARGET = true;
	BEAM_CURRENT = 200;
	
	IS_FIELD_OFF = true;
	
	FIT_EMPTY    = false;
	FIT_TRIPLET  = false;
	
	//==============================================================================//
	
	gStyle->SetOptStat(0); gStyle->SetOptFit(1);
	gStyle->SetEndErrorSize(4);
	
	bad_counters_tagh.clear();
	bad_counters_tagm.clear();
	bad_counters_tagh = {73, 186, 187, 220, 221};//1, 73, 124, 125, 126, 127, 179, 180, 181, 186, 187, 220, 221};
	bad_counters_tagm = {21, 63, 67, 100, 101, 102};//21, 63, 67, 100, 101, 102, 90};
	
	const char loc_pathname[256] = "/work/halld/home/andrsmit/primex_compton_analysis";
	
	//----------   Initialize   ----------//
	
	// Adjustable switches for monitoring:
	
	DEBUG_FITS       = false;
	
	DRAW_THEORY      = false;
	DRAW_ACCEPTANCE  = true;
	DRAW_FITS_TAGH   = true;
	DRAW_FITS_TAGM   = true;
	USE_F_ACC        = true;
	
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
	
	hname_tagh = Form("fcalE/deltaK_tagh_35fcalE");
	hname_tagm = Form("fcalE/deltaK_tagm_35fcalE");
	
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
	
	//==============================================================================//
	
	// plot results:
	
	plot_cs(tagh_vec, tagm_vec);
	plot_fit_fractions(tagh_vec, tagm_vec);
	plot_chi2(tagh_vec, tagm_vec);
	
	//==============================================================================//
	
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
