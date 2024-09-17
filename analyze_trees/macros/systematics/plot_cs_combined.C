#include "/work/halld/home/andrsmit/primex_compton_analysis/include/compton_inputs.cc"

void plot_cs_combined()
{
	vector<double> eb_vec, cs_vec, cs_diff_vec, stat_vec, sys_p2p_vec, sys_norm_vec;
	
	ifstream inf_tagh_unc("Be_200nA_unc_tagh.txt");
	ifstream inf_tagm_unc("Be_200nA_unc_tagm.txt");
	
	ifstream inf_tagh_cs("Be_200nA_tagh_cross_section.txt");
	ifstream inf_tagm_cs("Be_200nA_tagm_cross_section.txt");
	
	//------------------------------------------------------------------------------------//
	
	const char loc_pathname[256] = "/work/halld/home/andrsmit/primex_compton_analysis";
	
	endpoint_energy_calib = 11.6061;
	endpoint_energy       = 11.608;
	
	tagh_xscale_fname = Form("%s/photon_flux/phase1/primex_tagh.txt", loc_pathname);
	tagm_xscale_fname = Form("%s/photon_flux/phase1/primex_tagm.txt", loc_pathname);
	
	get_counter_energies();
	
	//------------------------//
	
	theory_cs_pathName = Form("%s/compton_mc/genDir/Run061321/genCS", loc_pathname);
	
	DRAW_THEORY = false;
	get_theory_calc();
	
	f_nist->SetLineColor(kAzure-8);
	f_nist->SetLineWidth(2);
	f_nist->SetLineStyle(2);
	
	//------------------------------------------------------------------------------------//
	
	double MIN_BEAM_ENERGY, MAX_BEAM_ENERGY;
	
	// low-energy:
	MIN_BEAM_ENERGY =  6.4;
	MAX_BEAM_ENERGY =  7.9;
	
	// medium-energy:
	MIN_BEAM_ENERGY =  7.9;
	MAX_BEAM_ENERGY =  8.9;
	
	// high-energy:
	MIN_BEAM_ENERGY =  6.4;
	MAX_BEAM_ENERGY = 11.2;
	
	int bin_counter = 0;
	double eb_save = 0., cs_save = 0., weight_save = 0., statE_save = 0., sysE_save = 0., normE_save = 0;
	
	int n_bins_combine_lowE  = 2;
	int n_bins_combine_midE  = 4;
	int n_bins_combine_highE = 2;
	
	int tagh_counter = 1;
	while(tagh_counter<=127) {
		int loc_counter;
		double loc_cs, loc_csE, loc_acc, loc_accE;
		double loc_stat_unc, loc_sys_unc_p2p, loc_sys_unc_norm;
		
		inf_tagh_unc >> loc_counter >> loc_stat_unc >> loc_sys_unc_p2p >> loc_sys_unc_norm;
		inf_tagh_cs  >> loc_cs      >> loc_csE      >> loc_acc         >> loc_accE;
		if(loc_cs <= 0. || loc_stat_unc <= 0.) {
			tagh_counter++;
			continue;
		}
		
		eb_save     += tagh_en[tagh_counter-1]*pow(1.0/loc_stat_unc,2.0);
		statE_save  += pow(loc_stat_unc/2.0 ,2.0);
		cs_save     += loc_cs*pow(1.0/loc_stat_unc,2.0);
		weight_save += pow(1.0/loc_stat_unc,2.0);
		sysE_save   += loc_sys_unc_p2p*pow(1.0/loc_stat_unc,2.0);
		normE_save  += loc_sys_unc_norm*pow(1.0/loc_stat_unc,2.0);
		
		bin_counter++;
		
		// check if stat unc is greater than 1%. If so, combine multiple bins:
		
		if(sqrt(statE_save)/(cs_save/weight_save) < 0.01 || bin_counter==n_bins_combine_highE) {
			
			double loc_eb = eb_save / weight_save;
			double loc_cs = cs_save / weight_save;
			
			double loc_cs_theory = f_theory->Eval(loc_eb);
			
			eb_vec.push_back(loc_eb);
			cs_vec.push_back(loc_cs);
			cs_diff_vec.push_back((loc_cs-loc_cs_theory)/loc_cs_theory);
			stat_vec.push_back(1.0/sqrt(weight_save));
			sys_p2p_vec.push_back(sysE_save/weight_save);
			sys_norm_vec.push_back(normE_save/weight_save);
			
			eb_save     = 0.;
			statE_save  = 0.;
			cs_save     = 0.;
			weight_save = 0.;
			sysE_save   = 0.;
			normE_save  = 0.;
			
			bin_counter = 0;
		}
		tagh_counter++;
	}
	if(eb_save > 0.) {
		
		double loc_eb = eb_save / weight_save;
		double loc_cs = cs_save / weight_save;
		
		double loc_cs_theory = f_theory->Eval(loc_eb);
		
		eb_vec.push_back(loc_eb);
		cs_vec.push_back(loc_cs);
		cs_diff_vec.push_back((loc_cs-loc_cs_theory)/loc_cs_theory);
		stat_vec.push_back(1.0/sqrt(weight_save));
		sys_p2p_vec.push_back(sysE_save/weight_save);
		sys_norm_vec.push_back(normE_save/weight_save);
		
		eb_save     = 0.;
		statE_save  = 0.;
		cs_save     = 0.;
		weight_save = 0.;
		sysE_save   = 0.;
		normE_save  = 0.;
		
		bin_counter = 0;
	}
	
	int tagm_counter = 1;
	while(tagm_counter<=102) {
		int loc_counter;
		double loc_cs, loc_csE, loc_acc, loc_accE;
		double loc_stat_unc, loc_sys_unc_p2p, loc_sys_unc_norm;
		
		inf_tagm_unc >> loc_counter >> loc_stat_unc >> loc_sys_unc_p2p >> loc_sys_unc_norm;
		inf_tagm_cs  >> loc_cs      >> loc_csE      >> loc_acc         >> loc_accE;
		if(loc_cs <= 0. || loc_stat_unc <= 0.) {
			tagm_counter++;
			continue;
		}
		
		eb_save     += tagm_en[tagm_counter-1]*pow(1.0/loc_stat_unc,2.0);
		statE_save  += pow(loc_stat_unc/2.0 ,2.0);
		cs_save     += loc_cs*pow(1.0/loc_stat_unc,2.0);
		weight_save += pow(1.0/loc_stat_unc,2.0);
		sysE_save   += loc_sys_unc_p2p*pow(1.0/loc_stat_unc,2.0);
		normE_save  += loc_sys_unc_norm*pow(1.0/loc_stat_unc,2.0);
		
		bin_counter++;
		
		// check if stat unc is greater than 1%. If so, combine multiple bins:
		
		if(sqrt(statE_save)/(cs_save/weight_save) < 0.01 || bin_counter==n_bins_combine_midE) {
			
			double loc_eb = eb_save / weight_save;
			double loc_cs = cs_save / weight_save;
			
			double loc_cs_theory = f_theory->Eval(loc_eb);
			
			eb_vec.push_back(loc_eb);
			cs_vec.push_back(loc_cs);
			cs_diff_vec.push_back((loc_cs-loc_cs_theory)/loc_cs_theory);
			stat_vec.push_back(1.0/sqrt(weight_save));
			sys_p2p_vec.push_back(sysE_save/weight_save);
			sys_norm_vec.push_back(normE_save/weight_save);
			
			eb_save     = 0.;
			statE_save  = 0.;
			cs_save     = 0.;
			weight_save = 0.;
			sysE_save   = 0.;
			normE_save  = 0.;
			
			bin_counter = 0;
		}
		tagm_counter++;
	}
	if(eb_save > 0.) {
		double loc_eb = eb_save / weight_save;
		double loc_cs = cs_save / weight_save;
		
		double loc_cs_theory = f_theory->Eval(loc_eb);
		
		eb_vec.push_back(loc_eb);
		cs_vec.push_back(loc_cs);
		cs_diff_vec.push_back((loc_cs-loc_cs_theory)/loc_cs_theory);
		stat_vec.push_back(1.0/sqrt(weight_save));
		sys_p2p_vec.push_back(sysE_save/weight_save);
		sys_norm_vec.push_back(normE_save/weight_save);
		
		eb_save     = 0.;
		statE_save  = 0.;
		cs_save     = 0.;
		weight_save = 0.;
		sysE_save   = 0.;
		normE_save  = 0.;
		
		bin_counter = 0;
	}
	
	while(tagh_counter<=274) {
		int loc_counter;
		double loc_cs, loc_csE, loc_acc, loc_accE;
		double loc_stat_unc, loc_sys_unc_p2p, loc_sys_unc_norm;
		
		inf_tagh_unc >> loc_counter >> loc_stat_unc >> loc_sys_unc_p2p >> loc_sys_unc_norm;
		inf_tagh_cs  >> loc_cs      >> loc_csE      >> loc_acc         >> loc_accE;
		if(loc_cs <= 0. || loc_stat_unc <= 0.) {
			tagh_counter++;
			continue;
		}
		
		eb_save     += tagh_en[tagh_counter-1]*pow(1.0/loc_stat_unc,2.0);
		statE_save  += pow(loc_stat_unc/2.0 ,2.0);
		cs_save     += loc_cs*pow(1.0/loc_stat_unc,2.0);
		weight_save += pow(1.0/loc_stat_unc,2.0);
		sysE_save   += loc_sys_unc_p2p*pow(1.0/loc_stat_unc,2.0);
		normE_save  += loc_sys_unc_norm*pow(1.0/loc_stat_unc,2.0);
		
		bin_counter++;
		
		// check if stat unc is greater than 1%. If so, combine multiple bins:
		
		if(sqrt(statE_save)/(cs_save/weight_save) < 0.01 || bin_counter==n_bins_combine_lowE) {
			
			double loc_eb = eb_save / weight_save;
			double loc_cs = cs_save / weight_save;
			
			double loc_cs_theory = f_theory->Eval(loc_eb);
			
			eb_vec.push_back(loc_eb);
			cs_vec.push_back(loc_cs);
			cs_diff_vec.push_back((loc_cs-loc_cs_theory)/loc_cs_theory);
			stat_vec.push_back(1.0/sqrt(weight_save));
			sys_p2p_vec.push_back(sysE_save/weight_save);
			sys_norm_vec.push_back(normE_save/weight_save);
			
			eb_save     = 0.;
			statE_save  = 0.;
			cs_save     = 0.;
			weight_save = 0.;
			sysE_save   = 0.;
			normE_save  = 0.;
			
			bin_counter = 0;
			
			if(loc_eb<7.4) n_bins_combine_lowE = 1;
		}
		tagh_counter++;
	}
	if(eb_save > 0.) {
		double loc_eb = eb_save / weight_save;
		double loc_cs = cs_save / weight_save;
		
		double loc_cs_theory = f_theory->Eval(loc_eb);
		
		eb_vec.push_back(loc_eb);
		cs_vec.push_back(loc_cs);
		cs_diff_vec.push_back((loc_cs-loc_cs_theory)/loc_cs_theory);
		stat_vec.push_back(1.0/sqrt(weight_save));
		sys_p2p_vec.push_back(sysE_save/weight_save);
		sys_norm_vec.push_back(normE_save/weight_save);
		
		eb_save     = 0.;
		statE_save  = 0.;
		cs_save     = 0.;
		weight_save = 0.;
		sysE_save   = 0.;
		normE_save  = 0.;
		
		bin_counter = 0;
	}
	inf_tagh_unc.close();
	inf_tagh_cs.close();
	
	int n_bins = (int)eb_vec.size();
	
	double *energy    = new double[n_bins];
	double *zeros     = new double[n_bins];
	
	double *cs        = new double[n_bins];
	double *statE     = new double[n_bins];
	double *sysE_p2p  = new double[n_bins];
	double *sysE_norm = new double[n_bins];
	double *band_min  = new double[n_bins];
	double *band_max  = new double[n_bins];
	double *norm      = new double[n_bins];
	
	double *cs_diff        = new double[n_bins];
	double *statE_diff     = new double[n_bins];
	double *sysE_p2p_diff  = new double[n_bins];
	double *sysE_norm_diff = new double[n_bins];
	double *band_min_diff  = new double[n_bins];
	double *band_max_diff  = new double[n_bins];
	double *norm_diff      = new double[n_bins];
	
	for(int ib=0; ib<n_bins; ib++) {
		energy[ib]  = eb_vec[ib];
		zeros[ib]   = 0.;
		
		double loc_sysE_p2p  = sys_p2p_vec[ib];
		double loc_sysE_norm = sys_norm_vec[ib];
		
		cs[ib]        = cs_vec[ib];
		statE[ib]     = (1.e-2*stat_vec[ib]) * cs_vec[ib];
		sysE_p2p[ib]  = 1.e-2 * sqrt(pow(loc_sysE_p2p,2.0) + pow(loc_sysE_norm,2.0)) * cs_vec[ib];
		sysE_norm[ib] = 1.e-2 * loc_sysE_norm * cs_vec[ib];
		
		cs_diff[ib]         = 1.e2*cs_diff_vec[ib];
		statE_diff[ib]      = stat_vec[ib];
		sysE_p2p_diff[ib]   = sqrt(pow(loc_sysE_p2p,2.0) + pow(loc_sysE_norm,2.0));
		sysE_norm_diff[ib]  = loc_sysE_norm;
		
		if(1) {
			band_min[ib] = f_theory->Eval(eb_vec[ib]) 
				- (1.e-2 * sqrt(pow(loc_sysE_p2p,2.0) + pow(loc_sysE_norm,2.0)) * cs_vec[ib]);
			band_max[ib] = f_theory->Eval(eb_vec[ib]) 
				+ (1.e-2 * sqrt(pow(loc_sysE_p2p,2.0) + pow(loc_sysE_norm,2.0)) * cs_vec[ib]);
			norm[ib]     = f_theory->Eval(eb_vec[ib]);
		} else {
			band_min[ib] = 0.12;// - sqrt(pow(loc_sysE_p2p,2.0) + pow(loc_sysE_norm,2.0));
			band_max[ib] = 0.12 + (1.e-2 * sqrt(pow(loc_sysE_p2p,2.0) + pow(loc_sysE_norm,2.0)) * cs_vec[ib]);
			norm[ib]     = 0.12;
		}
		
		if(1) {
			band_min_diff[ib] = -0.0 - sqrt(pow(loc_sysE_p2p,2.0) + pow(loc_sysE_norm,2.0));
			band_max_diff[ib] = -0.0 + sqrt(pow(loc_sysE_p2p,2.0) + pow(loc_sysE_norm,2.0));
			norm_diff[ib]     = -0.0;
		} else {
			band_min_diff[ib] = -8.0;// - sqrt(pow(loc_sysE_p2p,2.0) + pow(loc_sysE_norm,2.0));
			band_max_diff[ib] = -8.0 + sqrt(pow(loc_sysE_p2p,2.0) + pow(loc_sysE_norm,2.0));
			norm_diff[ib]     = -8.0;
		}
	}
	
	//-----------------------------------------------------------------//
	// Systematic uncertainty band for Cross Section plot:
	
	TGraph *g_band_min = new TGraph(n_bins, energy, band_min);
	TGraph *g_band_max = new TGraph(n_bins, energy, band_max);
	TGraph *g_band     = new TGraph(2*n_bins);
	for(int ib=0; ib<n_bins; ib++) {
		g_band->SetPoint(ib, energy[ib], band_max[ib]);
		g_band->SetPoint(n_bins+ib, energy[n_bins-ib-1], band_min[n_bins-ib-1]);
	}
	g_band->SetFillColor(kCyan);
	g_band->SetLineWidth(2);
	g_band->SetLineColor(kCyan);
	g_band_min->SetLineWidth(2);
	g_band_min->SetLineColor(kCyan);
	g_band_max->SetLineWidth(2);
	g_band_max->SetLineColor(kCyan);
	
	TGraphErrors *g_sysE = new TGraphErrors(n_bins, energy, norm, zeros, sysE_norm);
	g_sysE->SetFillColor(kCyan);
	
	//-----------------------------------------------------------------//
	// Systematic uncertainty band for Deviation plot:
	
	TGraph *g_band_min_diff = new TGraph(n_bins, energy, band_min_diff);
	TGraph *g_band_max_diff = new TGraph(n_bins, energy, band_max_diff);
	TGraph *g_band_diff     = new TGraph(2*n_bins);
	for(int ib=0; ib<n_bins; ib++) {
		g_band_diff->SetPoint(ib, energy[ib], band_max_diff[ib]);
		g_band_diff->SetPoint(n_bins+ib, energy[n_bins-ib-1], band_min_diff[n_bins-ib-1]);
	}
	g_band_diff->SetFillColor(kCyan);
	g_band_diff->SetLineWidth(2);
	g_band_diff->SetLineColor(kCyan);
	g_band_min_diff->SetLineWidth(2);
	g_band_min_diff->SetLineColor(kCyan);
	g_band_max_diff->SetLineWidth(2);
	g_band_max_diff->SetLineColor(kCyan);
	
	TGraphErrors *g_sysE_diff = new TGraphErrors(n_bins, energy, norm_diff, zeros, sysE_norm_diff);
	g_sysE_diff->SetFillColor(kCyan);
	
	
	//------------------------------------------------------------------//
	// Cross section plot with statistical error bars:
	
	TGraphErrors *gCS = new TGraphErrors(n_bins, energy, cs, zeros, statE);
	gCS->SetTitle("");
	gCS->SetMarkerStyle(24);
	gCS->SetMarkerSize(0.5);
	gCS->SetMarkerColor(kBlue);
	gCS->SetLineColor(kBlue);
	gCS->GetXaxis()->SetTitle("Photon Beam Energy [GeV]");
	gCS->GetXaxis()->SetTitleSize(0.05);
	gCS->GetXaxis()->SetTitleOffset(0.9);
	gCS->GetXaxis()->CenterTitle(true);
	gCS->GetYaxis()->SetTitle("#sigma [mb / electron]");
	gCS->GetYaxis()->SetTitleSize(0.05);
	gCS->GetYaxis()->SetTitleOffset(0.85);
	gCS->GetYaxis()->CenterTitle(true);
	gCS->GetXaxis()->SetRangeUser(6.2,11.4);
	gCS->GetYaxis()->SetRangeUser(0.11,0.24);
	
	gCS->SetLineWidth(2);
	gCS->GetYaxis()->SetTitleSize(0.07);
	gCS->GetYaxis()->SetTitleOffset(0.65);
	gCS->GetXaxis()->SetLabelOffset(1.0);
	
	gCS->GetYaxis()->SetTitleFont(12);
	
	// Systematic error bars (currently unused):
	
	TGraphErrors *gCS_sysE = new TGraphErrors(n_bins, energy, cs, zeros, sysE_p2p);
	gCS_sysE->SetTitle("");
	gCS_sysE->SetMarkerStyle(8);
	gCS_sysE->SetMarkerSize(0.5);
	gCS_sysE->SetMarkerColor(kMagenta);
	gCS_sysE->SetLineColor(kMagenta);
	//gCS_sysE->SetLineWidth(2);
	gCS_sysE->GetXaxis()->SetTitle("Photon Beam Energy [GeV]");
	gCS_sysE->GetXaxis()->SetTitleSize(0.05);
	gCS_sysE->GetXaxis()->SetTitleOffset(0.9);
	gCS_sysE->GetXaxis()->CenterTitle(true);
	gCS_sysE->GetYaxis()->SetTitle("#sigma [mb / electron]");
	gCS_sysE->GetYaxis()->SetTitleSize(0.05);
	gCS_sysE->GetYaxis()->SetTitleOffset(0.85);
	gCS_sysE->GetYaxis()->CenterTitle(true);
	gCS_sysE->GetXaxis()->SetRangeUser(6.2,11.4);
	gCS_sysE->GetYaxis()->SetRangeUser(0.11,0.24);
	
	TGraphErrors *gCS_clone = (TGraphErrors*)gCS->Clone("gCS_clone");
	gCS_clone->SetLineWidth(2);
	TGraphErrors *gCS_sysE_clone = (TGraphErrors*)gCS_sysE->Clone("gCS_sysE_clone");
	gCS_sysE_clone->SetLineWidth(2);
	
	
	//------------------------------------------------------------------//
	// Deviation plot with statistical error bars:
	
	TGraphErrors *gCS_diff = new TGraphErrors(n_bins, energy, cs_diff, zeros, statE_diff);
	gCS_diff->SetTitle("");
	gCS_diff->SetMarkerStyle(24);
	gCS_diff->SetMarkerSize(0.5);
	gCS_diff->SetMarkerColor(kBlue);
	gCS_diff->SetLineColor(kBlue);
	gCS_diff->GetXaxis()->SetTitle("Photon Beam Energy [GeV]");
	gCS_diff->GetXaxis()->SetTitleSize(0.05);
	gCS_diff->GetXaxis()->SetTitleOffset(0.9);
	gCS_diff->GetXaxis()->CenterTitle(true);
	gCS_diff->GetYaxis()->SetTitle("Deviation from theory [%]");
	gCS_diff->GetYaxis()->SetTitleSize(0.05);
	gCS_diff->GetYaxis()->SetTitleOffset(0.85);
	gCS_diff->GetYaxis()->CenterTitle(true);
	gCS_diff->GetXaxis()->SetRangeUser(6.2,11.4);
	gCS_diff->GetYaxis()->SetRangeUser(-7.0,5.0);
	
	gCS_diff->SetLineWidth(2);
	gCS_diff->GetYaxis()->SetTitleSize(0.1);
	gCS_diff->GetYaxis()->SetTitleOffset(0.35);
	gCS_diff->GetXaxis()->SetLabelSize(0.08);
	gCS_diff->GetYaxis()->SetLabelSize(0.075);
	gCS_diff->GetYaxis()->SetLabelOffset(0.01);
	gCS_diff->GetXaxis()->SetTitleSize(0.125);
	gCS_diff->GetXaxis()->SetTitleOffset(0.85);
	
	gCS_diff->GetXaxis()->SetTitleFont(12);
	gCS_diff->GetYaxis()->SetTitleFont(12);
	
	// Systematic error bars (currently unused):
	
	TGraphErrors *gCS_sysE_diff = new TGraphErrors(n_bins, energy, cs_diff, zeros, sysE_p2p_diff);
	gCS_sysE_diff->SetTitle("");
	gCS_sysE_diff->SetMarkerStyle(8);
	gCS_sysE_diff->SetMarkerSize(0.5);
	gCS_sysE_diff->SetMarkerColor(kMagenta);
	gCS_sysE_diff->SetLineColor(kMagenta);
	//gCS_sysE_diff->SetLineWidth(2);
	gCS_sysE_diff->GetXaxis()->SetTitle("Photon Beam Energy [GeV]");
	gCS_sysE_diff->GetXaxis()->SetTitleSize(0.05);
	gCS_sysE_diff->GetXaxis()->SetTitleOffset(0.9);
	gCS_sysE_diff->GetXaxis()->CenterTitle(true);
	gCS_sysE_diff->GetYaxis()->SetTitle("#sigma [mb / electron]");
	gCS_sysE_diff->GetYaxis()->SetTitleSize(0.05);
	gCS_sysE_diff->GetYaxis()->SetTitleOffset(0.85);
	gCS_sysE_diff->GetYaxis()->CenterTitle(true);
	gCS_sysE_diff->GetXaxis()->SetRangeUser(3.8,11.6);
	gCS_sysE_diff->GetYaxis()->SetRangeUser(0.11,0.24);
	
	TGraphErrors *gCS_clone_diff = (TGraphErrors*)gCS_diff->Clone("gCS_clone_diff");
	gCS_clone_diff->SetLineWidth(2);
	TGraphErrors *gCS_sysE_clone_diff = (TGraphErrors*)gCS_sysE_diff->Clone("gCS_sysE_clone_diff");
	gCS_sysE_clone_diff->SetLineWidth(2);
	
	TF1 *f_theory_norm = new TF1("f_theory_norm", "pol0", 5.0, 12.0);
	f_theory_norm->SetParameter(0, 0.0);
	f_theory_norm->SetNpx(1000);
	f_theory_norm->SetLineStyle(9);
	f_theory_norm->SetLineColor(kRed);
	f_theory_norm->SetLineWidth(2);
	
	//------------------------------//
	// Lee, PRL , 2021 comparison:
	
	vector<pair<double,double>> lee_cs_vec;
	double aa, bb, cc;
	ifstream inf("/work/halld/home/andrsmit/primex_compton_analysis/theory/lee_nlo_cs.txt");
	while(inf >> aa >> bb >> cc) {
		if(aa<6.0) continue;
		if(aa>11.2) continue;
		lee_cs_vec.push_back({aa,cc});
	}
	inf.close();
	
	int n_ens = (int)lee_cs_vec.size();
	double *ebeam       = new double[n_ens];
	double *lee_cs      = new double[n_ens];
	double *lee_cs_diff = new double[n_ens];
	
	for(int ib=0; ib<n_ens; ib++) {
		ebeam[ib]  = lee_cs_vec[ib].first;
		lee_cs[ib] = lee_cs_vec[ib].second;
		
		double loc_percent_diff = 1.e2*(lee_cs_vec[ib].second - f_theory->Eval(lee_cs_vec[ib].first))
			/ f_theory->Eval(lee_cs_vec[ib].first);
		lee_cs_diff[ib] = loc_percent_diff;
		//printf("E = %.2f GeV: Percent Difference: %.2f%%\n", lee_cs_vec[ib].first, loc_percent_diff);
	}
	
	TGraph *gLeeCS = new TGraph(n_ens, ebeam, lee_cs);
	gLeeCS->SetMarkerStyle(4);
	gLeeCS->SetMarkerSize(0.1);
	gLeeCS->SetMarkerColor(kBlack);
	gLeeCS->SetLineColor(kBlack);
	gLeeCS->SetLineWidth(2);
	gLeeCS->SetLineStyle(8);
	
	TGraph *gLeeCS_diff = new TGraph(n_ens, ebeam, lee_cs_diff);
	gLeeCS_diff->SetMarkerStyle(4);
	gLeeCS_diff->SetMarkerSize(0.1);
	gLeeCS_diff->SetMarkerColor(kBlack);
	gLeeCS_diff->SetLineColor(kBlack);
	gLeeCS_diff->SetLineWidth(2);
	gLeeCS_diff->SetLineStyle(8);
	
	//------------------------------//
	// NIST Comparison:
	
	int n_ens_nist = g_nist->GetN();
	double *nist_ens     = new double[n_ens_nist];
	double *nist_cs_diff = new double[n_ens_nist];
	for(int ib=0; ib<n_ens_nist; ib++) {
		double loc_eb = g_nist->GetX()[ib];
		double loc_cs = g_nist->GetY()[ib];
		double loc_percent_diff = 1.e2*(loc_cs - f_theory->Eval(loc_eb))/loc_cs;
		nist_ens[ib]     = loc_eb;
		nist_cs_diff[ib] = loc_percent_diff;
	}
	
	TGraph *gNISTCS_diff = new TGraph(n_ens_nist, nist_ens, nist_cs_diff);
	gNISTCS_diff->SetMarkerStyle(4);
	gNISTCS_diff->SetMarkerSize(0.1);
	gNISTCS_diff->SetMarkerColor(kAzure-8);
	gNISTCS_diff->SetLineColor(kAzure-8);
	gNISTCS_diff->SetLineWidth(2);
	gNISTCS_diff->SetLineStyle(2);
	
	
	
	//--------------------------------------------------------------//
	// Cross Section Plot:
	
	TCanvas *cCS = new TCanvas("cCS", "cCS", 1000, 800);
	
	TPad *top_pad = new TPad("top_pad", "top_pad", 0.005, 0.3025, 0.995, 0.995);
	TPad *bot_pad = new TPad("bot_pad", "bot_pad", 0.005, 0.005,  0.995, 0.3025);
	
	top_pad->SetLeftMargin(0.10);
	top_pad->SetRightMargin(0.02);
	top_pad->SetTopMargin(0.075);
	top_pad->SetBottomMargin(0.0);//15);
	top_pad->SetTickx(); top_pad->SetTicky();
	top_pad->SetFrameLineWidth(2);
	
	bot_pad->SetLeftMargin(0.10);
	bot_pad->SetRightMargin(0.02);
	bot_pad->SetTopMargin(0.0);//10);
	bot_pad->SetBottomMargin(0.225);
	bot_pad->SetTickx(); bot_pad->SetTicky();
	bot_pad->SetFrameLineWidth(2);
	
	cCS->cd();
	top_pad->Draw();
	bot_pad->Draw();
	
	top_pad->cd();
	gCS->Draw("AP");
	g_band->Draw("F same");
	f_theory->Draw("same");
	f_nist->Draw("same");
	gLeeCS->Draw("C same");
	gCS->Draw("PE1 same");
	
	
	TLegend *leg = new TLegend(0.495, 0.575, 0.882, 0.853);
	leg->SetBorderSize(0);
	leg->AddEntry(f_theory, "NLO Calculation (PrimEx) [12]", "l");
	leg->AddEntry(gLeeCS_diff, "Lee, #font[52]{et al} (2021) [16]", "l");
	leg->AddEntry(f_nist, "NIST XCOM [17]", "l");
	leg->AddEntry(gCS_clone, "Statistical Uncertainty", "PE");
	leg->AddEntry(g_band, "Systematic Uncertainty", "f");
	leg->Draw();
	
	//---------------------------------------------------------------//
	// Deviation Plot:
	
	bot_pad->cd();
	gCS_diff->Draw("AP");
	g_band_diff->Draw("F same");
	f_theory_norm->Draw("same");
	gLeeCS_diff->Draw("C same");
	gNISTCS_diff->Draw("C same");
	gCS_diff->Draw("PE1 same");
	
	return;
}
