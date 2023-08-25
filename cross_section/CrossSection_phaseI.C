
#include "compton_inc.h"
#include "get_compton_inputs.C"
#include "get_compton_acc.C"
#include "get_compton_yield.C"
#include "calc_cs.C"
#include "plot_cs.C"
#include "plot_pair_fraction.C"
#include "plot_chi2.C"

//----------   Function Declarations   ----------//

void print_status();

//-----------------------------------------------//

void CrossSection_phaseI()
{
	gStyle->SetOptStat(0); gStyle->SetOptFit(1);
	gStyle->SetEndErrorSize(4);
	
	bad_counters_tagh.clear();
	bad_counters_tagm.clear();
	bad_counters_tagh = {1, 73, 124, 125, 126, 127, 186, 187};
	bad_counters_tagm = {21, 63, 67, 100, 101, 102, 90};
	
	const char pathName[256] = "/work/halld/home/andrsmit/primex_compton_analysis";
	
	//----------   Initialize   ----------//
	
	// Set these parameters to select data set:
	
	PHASE_VAL    = 1; // RunPeriod-2019-01
	
	IS_BE_TARGET = true;
	IS_FIELD_OFF = true;
	BEAM_CURRENT = 200;
	
	// Adjustable switches for monitoring:
	
	DRAW_THEORY       = false;
	DRAW_ACCEPTANCE   = true;
	DRAW_FITS_TAGH    = false;
	DRAW_FITS_TAGM    = false;
	USE_F_ACC         = true;
	
	// adjust bin size. Default: 2000 bins in DeltaK2 histogram (8MeV per bin)
	
	rebins    =  5;
	bin_size *=  (double)rebins;
	n_mev     =  4 * rebins;
	
	//==============================================================================//
	
	if(IS_BE_TARGET) sprintf(TARGET_STR, "Be");
	else             sprintf(TARGET_STR, "He");
	
	if(IS_FIELD_OFF) sprintf(FIELD_STR, "FIELDOFF");
	else             sprintf(FIELD_STR, "FIELDON" );
	
	//----------------------------------------//
	// Photon beam energy:
	
	endpoint_energy_calib = 11.6061;
	endpoint_energy       = 11.608;
	
	sprintf(tagh_xscale_fname, "%s/photon_flux/phaseI/primex_tagh.txt", pathName);
	sprintf(tagm_xscale_fname, "%s/photon_flux/phaseI/primex_tagm.txt", pathName);
	
	//----------------------------------------//
	// Photon flux:
	
	sprintf(tagh_flux_fname, "%s/photon_flux/phaseI/%s_%03dnA_%s_flux_tagh.txt", 
		pathName, TARGET_STR, BEAM_CURRENT, FIELD_STR);
	sprintf(tagm_flux_fname, "%s/photon_flux/phaseI/%s_%03dnA_%s_flux_tagm.txt", 
		pathName, TARGET_STR, BEAM_CURRENT, FIELD_STR);
	
	sprintf(empty_target_tagh_flux_fname, 
		"%s/photon_flux/phaseI/%s_empty_%s_flux_tagh.txt", 
			pathName, TARGET_STR, FIELD_STR);
	sprintf(empty_target_tagm_flux_fname, 
		"%s/photon_flux/phaseI/%s_empty_%s_flux_tagm.txt", 
			pathName, TARGET_STR, FIELD_STR);
	
	//----------------------------------------//
	// Compton & e+e- Pair Production Theory:
	
	sprintf(theory_cs_pathName, "%s/compton_mc/phaseI/genDir/Be/genCS", pathName);
	sprintf(pair_cs_fname,      "%s/photon_absorption/%s_pair_cs.dat",  pathName, TARGET_STR);
	
	//----------------------------------------//
	// Data:
	
	char loc_buf[256];
	
	sprintf(loc_buf, 
		"%s/analyze_trees/phaseI/systematics/rootFiles/%s_%03dnA_%s.root", 
		pathName, TARGET_STR, BEAM_CURRENT, FIELD_STR);
	root_fname = loc_buf;
	
	sprintf(loc_buf, 
		"%s/analyze_trees/phaseI/systematics/rootFiles/%s_empty_%s.root", 
		pathName, TARGET_STR, FIELD_STR);
	empty_target_root_fname = loc_buf;
	
	hname_tagh = "fcalE/deltaK_tagh_cut_35fcalE";
	hname_tagm = "fcalE/deltaK_tagm_cut_35fcalE";
	
	//----------------------------------------//
	// Simulation:
	
	// Compton MC:
	
	sprintf(loc_buf, "%s/compton_mc/phaseI/default_geometry/Be/recRootFiles_sys/fcalE_cut_500MeV", 
		pathName);
	comp_mc_dir        = loc_buf;
	comp_root_dir_name = "compton_simulation_systematics_deltaK2";
	hname_tagh_comp    = "fcalE/deltaK_tagh_35fcalE";
	hname_tagm_comp    = "fcalE/deltaK_tagm_35fcalE";
	
	// e+e- Pair Production MC:
	
	sprintf(loc_buf, "%s/pair_sim/phaseI/Be/recRootFiles_pair_sys_deltaK2", pathName);
	pair_mc_dir        = loc_buf;
	pair_root_dir_name = "compton_simulation_systematics_pair_deltaK2";
	hname_tagh_pair    = hname_tagh_comp;
	hname_tagm_pair    = hname_tagm_comp;
	
	// e+e- Triplet Production MC:
	
	sprintf(loc_buf, "%s/pair_sim/phaseI/Be/recRootFiles_trip_sys_deltaK2", pathName);
	trip_mc_dir        = loc_buf;
	trip_root_dir_name = "compton_simulation_systematics_pair_deltaK2";
	hname_tagh_trip    = hname_tagh_comp;
	hname_tagm_trip    = hname_tagm_comp;
	
	//==============================================================================//
	
	// Check that all files exist:
	
	if(gSystem->AccessPathName(root_fname)) {
		cout << "Specified ROOT file does not exist, check file name.\n\n";
		return;
	}
	if(gSystem->AccessPathName(empty_target_root_fname)) {
		cout << "Specified empty target ROOT file does not exist, check file name.\n\n" << endl;
		return;
	}
	if(gSystem->AccessPathName(tagh_flux_fname) || gSystem->AccessPathName(tagm_flux_fname)) {
		cout << "Specified flux files do not exist, check file name.\n\n" << endl;
		return;
	}
	if(gSystem->AccessPathName(empty_target_tagh_flux_fname) 
		|| gSystem->AccessPathName(empty_target_tagh_flux_fname)) {
		cout << "Specified empty target flux files do not exist, check file name.\n\n" << endl;
		return;
	}
	
	//==============================================================================//
	
	print_status();
	
	//==============================================================================//
	
	get_target_parameters(IS_BE_TARGET);
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
	
	//==============================================================================//
	
	// plot results:
	
	plot_cs(tagh_vec, tagm_vec);
	plot_pair_fraction(tagh_vec, tagm_vec);
	plot_chi2(tagh_vec, tagm_vec);
	
	//==============================================================================//
	
	return;
}

void print_status() {
	
	cout << "\n\n";
	cout << "***************************************************************************";
	cout << "*****************************" << endl;
	printf("Calculating total Compton Scattering Cross Section for %03dnA %s target (%s):\n\n", 
		BEAM_CURRENT, TARGET_STR, FIELD_STR);
	printf(" Filled target ROOT File: \n");
	printf("  %s\n", root_fname.Data());
	printf(" Empty target ROOT File: \n");
	printf("  %s\n\n", empty_target_root_fname.Data());
	printf(" Filled target flux accessed from: \n");
	printf("  %s\n", tagh_flux_fname);
	printf("  %s\n", tagm_flux_fname);
	printf(" Empty target flux accessed from: \n");
	printf("  %s\n", empty_target_tagh_flux_fname);
	printf("  %s\n\n", empty_target_tagm_flux_fname);
	printf(" Compton MC directory: \n");
	printf("  %s\n", comp_mc_dir.Data());
	printf(" Theory calculation (primex): \n");
	printf("  %s\n", theory_cs_pathName);
	printf(" Theory calculation (NIST): \n");
	printf("  %s\n\n", nist_theory_fname);
	printf(" e+e- Pair MC directory: \n");
	printf("  %s\n", pair_mc_dir.Data());
	printf(" e+e- Trip MC directory: \n");
	printf("  %s\n", trip_mc_dir.Data());
	cout << "***************************************************************************";
	cout << "*****************************" << endl;
	cout << "\n\n";
	
	return;
}
