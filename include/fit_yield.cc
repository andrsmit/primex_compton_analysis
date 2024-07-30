#include "compton_cs.h"
#include "get_phase1_energy_bin.cc"

int get_neighbor(int tag_sys, int counter, int sim_case);

int fit_yield(int tag_sys, int counter, TH1F *h1, TH1F *h1_empty, double &yield, double &yieldE, double &chi2) {
	
	map<TString, pair<double,TH1F*>> mc_hists;
	
	int fit_val = 1;
	
	char tag_sys_char[256];
	double eb, loc_flux;
	if(tag_sys==0) { 
		sprintf(tag_sys_char, "tagh");
		eb        = tagh_en[counter-1];
		loc_flux  = tagh_flux[counter-1];
	} else {
		sprintf(tag_sys_char, "tagm");
		eb        = tagm_en[counter-1];
		loc_flux  = tagm_flux[counter-1];
	}
	double loc_acc = f_acc->Eval(eb);
	
	double min_fit_x = h1->GetBinCenter(h1->FindFirstBinAbove()+5);
	double max_fit_x = h1->GetBinCenter(h1->FindLastBinAbove()-5);
	
	double loc_deltaK_mu  = -1.05827e-01 + 8.60222e-02*eb - 1.84818e-02*pow(eb,2.0) + 8.04656e-04*pow(eb,3.0);
	double loc_deltaK_sig = -1.38713e-01 + 1.98660e-01*eb - 1.88848e-02*pow(eb,2.0) + 7.97340e-04*pow(eb,3.0);
	
	h1->GetXaxis()->SetRangeUser(min_fit_x,max_fit_x);
	h1->SetLineColor(kBlack);
	h1->SetMarkerColor(kBlack);
	h1->SetMarkerStyle(8);
	h1->SetMarkerSize(0.6);
	h1->GetXaxis()->SetTitleSize(0.045);
	h1->GetXaxis()->SetTitleOffset(0.95);
	h1->GetXaxis()->SetLabelSize(0.0325);
	
	// Find the energy bin from phase I closest to the new energy bin:
	
	int phase1_tag_counter = 0;
	int phase1_tag_sys     = 0;
	get_phase1_energy_bin(eb, phase1_tag_counter, phase1_tag_sys);
	
	char loc_tag_name[256];
	if(!tag_sys) sprintf(loc_tag_name, "TAGH %d", counter);
	else         sprintf(loc_tag_name, "TAGM %d", counter);
	
	char phase1_tag_name[256];
	if(!phase1_tag_sys) sprintf(phase1_tag_name, "TAGH %d", phase1_tag_counter);
	else                sprintf(phase1_tag_name, "TAGM %d", phase1_tag_counter);
	
	cout << "current bin: " << loc_tag_name << "; phase 1 bin: " << phase1_tag_name << endl;
	
	// How many bins to combine for the e+e- simulation:
	
	int n_bins_combine      = 10; 
	int min_counter_pair_mc = phase1_tag_counter-n_bins_combine;
	int max_counter_pair_mc = phase1_tag_counter+n_bins_combine;
	
	if(phase1_tag_sys==0) {
		if(min_counter_pair_mc<1)   min_counter_pair_mc = 1;
		else if(min_counter_pair_mc>127 && min_counter_pair_mc<179) min_counter_pair_mc = 179;
		
		if(max_counter_pair_mc>221) max_counter_pair_mc = 221;
		else if(max_counter_pair_mc>127 && max_counter_pair_mc<179) max_counter_pair_mc = 127;
	} else {
		if(max_counter_pair_mc<1)   max_counter_pair_mc = 1;
		if(max_counter_pair_mc>102) max_counter_pair_mc = 102;
	}
	
	//=================================================================================//
	// Get the Compton MC Histogram
	
	//cout << "getting compton mc" << endl;
	
	char fname[256];
	if(tag_sys==0) 
		sprintf(fname, "%s/tagh_%03d.root", comp_mc_dir.Data(), counter);
	else 
		sprintf(fname, "%s/tagm_%03d.root", comp_mc_dir.Data(), counter);
	
	
	int neighbor_counter = counter;
	if(gSystem->AccessPathName(fname)) {
		
		//return 0;
		neighbor_counter = get_neighbor(tag_sys, counter, 0);
		if(neighbor_counter<=0) return 0;
		else {
			if(tag_sys==0) 
				sprintf(fname, "%s/tagh_%03d.root", comp_mc_dir.Data(), 
					neighbor_counter);
			else 
				sprintf(fname, "%s/tagm_%03d.root", comp_mc_dir.Data(), 
					neighbor_counter);
		}
	}
	
	TFile *fSim = new TFile(fname, "READ");
	
	TH1F *h_vertex = (TH1F*)fSim->Get("vertex_accepted");
	
	double n_comp_thrown = h_vertex->Integral();
	
	double comp_cs       = f_theory->Eval(eb);
	double comp_flux_sim = n_comp_thrown / (comp_cs * n_e * mb * f_abs);
	
	if(n_comp_thrown < 1.e3) return 0;
	/*
	TH2F *h2_compton_tagh = (TH2F*)fSim->Get(Form("%s", hname_tagh_comp.Data()));
	TH2F *h2_compton_tagm = (TH2F*)fSim->Get(Form("%s", hname_tagm_comp.Data()));
	TH1F *h1_compton      = (TH1F*)h2_compton_tagh->ProjectionY(Form("h1_comp_%d_%d", tag_sys, counter));
	h1_compton->Add((TH1F*)h2_compton_tagm->ProjectionY(Form("h1_comp_tagm_%d_%d", tag_sys, counter)));
	*/
	
	TH2F *h2_compton;
	if(tag_sys==0) 
		h2_compton = (TH2F*)fSim->Get(Form("%s", hname_tagh_comp.Data()));
	else 
		h2_compton = (TH2F*)fSim->Get(Form("%s", hname_tagm_comp.Data()));
	TH1F *h1_compton = (TH1F*)h2_compton->ProjectionY(Form("h1_comp_%d_%d", tag_sys, counter));
	
	h1_compton->SetDirectory(0);
	h1_compton->Rebin(rebins);
	h1_compton->Scale(loc_flux/comp_flux_sim);
	
	fSim->Close();
	
	if(h1_compton->Integral() < 1.e1) return 0;
	
	TH1F *h1_compton_fit = (TH1F*)h1_compton->Clone(Form("h1_compton_fit_%d_%d", tag_sys, counter));
	
	mc_hists["compton"] = {0.0, h1_compton_fit};
	
	//=================================================================================//
	// Get the Pair MC Histogram
	
	TH1F *h1_pair;
	
	// generated flux file:
	sprintf(fname, "/work/halld/home/andrsmit/primex_compton_analysis/bhgen_test/recRootTrees/Run061321/sum.root");
	if(gSystem->AccessPathName(fname)) {
		cout << "Can't access e+e- pair mc generated flux file." << endl;
		return 0;
	}
	TFile *fPairFlux = new TFile(fname, "READ");
	TH1F *h_pair_gen_flux = (TH1F*)fPairFlux->Get("gen_flux");
	
	// histogram file:
	if(FIT_TRIPLET) {
		sprintf(fname, "%s/pair_rec.root",          pair_mc_dir.Data());
	} else {
		sprintf(fname, "%s/pair_rec_combined.root", pair_mc_dir.Data());
	}
	if(gSystem->AccessPathName(fname)) {
		cout << "Can't access e+e- pair mc histogram file." << endl;
		return 0;
	}
	TFile *fPair = new TFile(fname, "READ");
	
	double loc_pair_gen_flux = 0.;
	
	if(phase1_tag_sys==0) {
		
		TH2F *h2_pair = (TH2F*)fPair->Get(hname_tagh_pair.Data());
		h1_pair = (TH1F*)h2_pair->ProjectionY(Form("h1_pair_%d_%d", tag_sys, counter), 
			min_counter_pair_mc, max_counter_pair_mc);
		
		for(int loc_counter = min_counter_pair_mc; loc_counter <= max_counter_pair_mc; loc_counter++) {
			double loc_counter_eb   = tagh_en_phase1[loc_counter-1];
			double loc_counter_flux = h_pair_gen_flux->GetBinContent(h_pair_gen_flux->FindBin(loc_counter_eb));
			loc_pair_gen_flux += loc_counter_flux;
		}
	} else {
		
		TH2F *h2_pair = (TH2F*)fPair->Get(hname_tagm_pair.Data());
		h1_pair = (TH1F*)h2_pair->ProjectionY(Form("h1_pair_%d_%d", tag_sys, counter), 
			min_counter_pair_mc, max_counter_pair_mc);
		
		for(int loc_counter = min_counter_pair_mc; loc_counter <= max_counter_pair_mc; loc_counter++) {
			double loc_counter_eb   = tagm_en_phase1[loc_counter-1];
			double loc_counter_flux = h_pair_gen_flux->GetBinContent(h_pair_gen_flux->FindBin(loc_counter_eb));
			loc_pair_gen_flux += loc_counter_flux;
		}
	}
	h1_pair->SetDirectory(0);
	
	fPair->Close();
	fPairFlux->Close();
	
	double loc_pair_cs   = f_pair_cs->Eval(eb) + f_triplet_cs->Eval(eb);
	double loc_pair_flux = loc_pair_gen_flux / ((n_e/n_Z) * mb * loc_pair_cs);
	
	h1_pair->Scale(loc_flux/loc_pair_flux);
	h1_pair->Rebin(rebins);
	
	mc_hists["pair"] = {0.0, h1_pair};
	
	//=================================================================================//
	// (optionally) Get the Triplet MC Histogram:
	
	if(FIT_TRIPLET) {
		
		TH1F *h1_triplet;
		
		sprintf(fname, "%s/triplet_rec.root", pair_mc_dir.Data());
		if(gSystem->AccessPathName(fname)) {
			cout << "Can't access e+e- triplet mc histogram file." << endl;
			return 0;
		}
		TFile *fTriplet = new TFile(fname, "READ");
		
		if(phase1_tag_sys==0) {
			TH2F *h2_triplet = (TH2F*)fTriplet->Get(hname_tagh_pair.Data());
			h1_triplet = (TH1F*)h2_triplet->ProjectionY(Form("h1_triplet_%d_%d", tag_sys, counter), 
				min_counter_pair_mc, max_counter_pair_mc);
		} else {
			TH2F *h2_triplet = (TH2F*)fTriplet->Get(hname_tagm_pair.Data());
			h1_triplet = (TH1F*)h2_triplet->ProjectionY(Form("h1_triplet_%d_%d", tag_sys, counter), 
				min_counter_pair_mc, max_counter_pair_mc);
		}
		h1_triplet->SetDirectory(0);
		
		h1_triplet->Scale(loc_flux/loc_pair_flux);
		h1_triplet->Rebin(rebins);
		
		mc_hists["triplet"] = {0.0, h1_triplet};
	}
	
	//=================================================================================//
	// (optionally) add the empty target distribution to the fit templates:
	
	if(FIT_EMPTY) mc_hists["empty"] = {0.0, h1_empty};
	else h1->Add(h1_empty, -1.0);
	
	//=================================================================================//
	
	TH1F *h1_fit = (TH1F*)h1->Clone(Form("h1_fit_%d_%d",tag_sys,counter));
	h1_fit->GetXaxis()->SetRangeUser(min_fit_x, max_fit_x);
	double n_total_exp = h1->Integral(h1->FindBin(min_fit_x), h1->FindBin(max_fit_x));
	for(int ib=1; ib<=h1_fit->GetXaxis()->GetNbins(); ib++) {
		if(h1_fit->GetBinContent(ib) < 0.) h1_fit->SetBinContent(ib, 0.);
	}
	
	int n_mc_hists = (int)mc_hists.size();
	map<TString, pair<double, TH1F*>>::iterator mc_hist_it;
	for(mc_hist_it = mc_hists.begin(); mc_hist_it != mc_hists.end(); mc_hist_it++) {
		mc_hist_it->second.second->GetXaxis()->SetTitle("E_{Comp} - E_{#gamma} [GeV]");
		mc_hist_it->second.second->GetXaxis()->SetTitle("E_{Comp} - E_{#gamma} [GeV]");
		mc_hist_it->second.second->GetXaxis()->SetTitleOffset(1.1);
		mc_hist_it->second.second->GetYaxis()->SetTitle( Form("counts / %d MeV", n_mev));
		mc_hist_it->second.second->GetYaxis()->SetTitleOffset(1.3);
		
		unsigned int loc_color = 0;
		map<TString, unsigned int>::iterator color_it = hist_color_map.find(mc_hist_it->first);
		if(color_it != hist_color_map.end()) loc_color = color_it->second;
		
		mc_hist_it->second.second->SetLineColor(loc_color);
		mc_hist_it->second.second->SetFillColor(loc_color);
		mc_hist_it->second.second->SetFillStyle(3004);
		mc_hist_it->second.second->SetLineWidth(2);
		
		for(int ib=1; ib<=mc_hist_it->second.second->GetXaxis()->GetNbins(); ib++) {
			if(mc_hist_it->second.second->GetBinContent(ib)<0.) mc_hist_it->second.second->SetBinContent(ib, 0.);
		}
		
		mc_hist_it->second.second->GetXaxis()->SetRangeUser(min_fit_x, max_fit_x);
		double n_mc_exp = mc_hist_it->second.second->Integral(mc_hist_it->second.second->FindBin(min_fit_x), 
			mc_hist_it->second.second->FindBin(max_fit_x));
		mc_hist_it->second.first = n_mc_exp / n_total_exp;
	}
	
	//-----------------------------------------------------------------------//
	
	if(DEBUG_FITS) {
		if(c_debug==NULL) {
			c_debug = new TCanvas("c_debug", "Fit Debug", 800, 800);
			c_debug->SetTickx(); c_debug->SetTicky();
		}
		
		c_debug->cd();
		h1->Draw("PE");
		
		TH1F *h1_debug = (TH1F*)h1->Clone("h1_debug");
		h1_debug->SetLineColor(kRed);
		h1_debug->SetLineWidth(2);
		
		for(mc_hist_it = mc_hists.begin(); mc_hist_it != mc_hists.end(); mc_hist_it++) {
			mc_hist_it->second.second->Draw("same hist");
			h1_debug->Add(mc_hist_it->second.second);
		}
		h1_debug->Add(h1,-1.0);
		h1_debug->Draw("same hist");
		
		c_debug->Update();
	}
	
	double loc_yield1 = h1_compton->Integral(h1_compton->FindBin(min_fit_x), h1_compton->FindBin(max_fit_x));
	double loc_yield2 = h1_compton->Integral(h1_compton->FindBin(-8.0),      h1_compton->FindBin(8.0));
	
	// If the fraction of pairs is low enough, don't bother with fit:
	
	if(mc_hists["pair"].first < 0.01) {
		
		double fitted_compton_frac = 1.0;
		for(mc_hist_it = mc_hists.begin(); mc_hist_it != mc_hists.end(); mc_hist_it++) {
			fitted_compton_frac -= mc_hist_it->second.first;
		}
		yield  = n_total_exp * fitted_compton_frac;
		yield *= (loc_yield1/loc_yield2);
		yieldE = sqrt(yield);
		
		chi2 = 1.0;
		return fit_val;
	}
	
	//-----   Perform Binned Maximum Likelihood Fit (TFractionFitter)   -----//
	
	TObjArray *mc = new TObjArray(n_mc_hists);
	
	// store the index of each mc template added to the TObjArray in a map with the same keys:
	map<TString, int> mc_index;
	int loc_index = 0;
	for(mc_hist_it = mc_hists.begin(); mc_hist_it != mc_hists.end(); mc_hist_it++) {
		mc->Add(mc_hist_it->second.second);
		mc_index[mc_hist_it->first] = loc_index;
		loc_index++;
	}
	
	TFractionFitter *fit = new TFractionFitter(h1_fit, mc);
	fit->SetRangeX(h1_fit->FindBin(min_fit_x), h1_fit->FindBin(max_fit_x));
	
	for(int ihist=0; ihist<n_mc_hists; ihist++) fit->Constrain(ihist, 0.0, 1.0);
	
	if(FIT_TRIPLET) {
		int triplet_fit_index = mc_index["triplet"];
		double loc_frac = mc_hists["triplet"].first;
		fit->Constrain(triplet_fit_index, 0.90*loc_frac, 1.10*loc_frac);
	}
	
	Int_t fit_status = fit->Fit();
	TH1F *fit_result = (TH1F*)fit->GetPlot();
	
	fit_result->SetLineColor(kRed+1);
	fit_result->GetXaxis()->SetRangeUser(min_fit_x, max_fit_x);
	fit_result->SetMaximum(fit_result->GetMaximum()*1.4);
	
	double Nmc   = fit_result->Integral(fit_result->FindBin(min_fit_x), fit_result->FindBin(max_fit_x));
	double Ndata = h1->Integral(h1->FindBin(min_fit_x), h1->FindBin(max_fit_x));
	
	vector<pair<TString,TH1F*>> FitResults; FitResults.clear();
	
	TH1F *FitResultAll = (TH1F*)h1->Clone(Form("FitResultAll_%d_%d", tag_sys, counter));
	FitResultAll->Reset();
	
	double fitted_compton_frac = 1.0, fitted_compton_fracE = 0.0;
	for(mc_hist_it = mc_hists.begin(); mc_hist_it != mc_hists.end(); mc_hist_it++) {
		
		TString loc_template = mc_hist_it->first;
		map<TString, int>::iterator index_it = mc_index.find(loc_template);
		if(index_it == mc_index.end()) {
			cout << "Unable to find histogram index for MC template: " << loc_template << endl;
			return 0;
		}
		
		double loc_frac, loc_fracE;
		fit->GetResult(index_it->second, loc_frac, loc_fracE);
		
		if(loc_template != "compton") {
			fitted_compton_frac  -= loc_frac;
			fitted_compton_fracE += pow(loc_fracE,2.0);
		}
		
		double loc_frac_exp = mc_hist_it->second.first;
		
		// update the tagh(m)_fit_fraction maps:
		if(tag_sys==0) {
			tagh_fit_fraction[loc_template][counter-1]     = {loc_frac, loc_fracE};
			tagh_fit_fraction_exp[loc_template][counter-1] = {loc_frac_exp, 0.0};
		} else {
			tagm_fit_fraction[loc_template][counter-1]     = {loc_frac, loc_fracE};
			tagm_fit_fraction_exp[loc_template][counter-1] = {loc_frac_exp, 0.0};
		}
		
		unsigned int loc_color = 0;
		map<TString, unsigned int>::iterator color_it = hist_color_map.find(mc_hist_it->first);
		if(color_it != hist_color_map.end()) loc_color = color_it->second;
		
		// Get individual fit result for this MC template:
		/*
		TH1F *loc_mc_fit = (TH1F*)->GetMCPrediction(index_it->second);
		loc_mc_fit->GetXaxiS()->SetRangeUser(min_fit_x, max_fit_x);
		loc_mc_fit->Scale(loc_frac*Ndata/loc_mc_fit->Integral(loc_mc_fit->FindBin(min_fit_x),
			loc_mc-fit->FindBin(max_fit_x)));
		loc_mc_fit->SetLineColor(loc_color);
		loc_mc_fit->SetMarkerColor(loc_color);
		loc_mc_fit->SetFillColor(loc_color);
		loc_mc_fit->SetFillStyle(3004);
		loc_mc_fit->SetLineWidth(2);
		loc_mc_fit->SetMarkerStyle(8);
		loc_mc_fit->SetMarkerSize(0.6);
		*/
		
		TH1F *loc_FitResult = (TH1F*)mc_hist_it->second.second->Clone(
			Form("FitResult_%s_%d_%d", loc_template.Data(), tag_sys, counter));
		loc_FitResult->Scale(loc_frac*Ndata/loc_FitResult->Integral());
		loc_FitResult->SetLineColor(loc_color+2);
		loc_FitResult->SetFillColor(loc_color-9);
		loc_FitResult->SetFillStyle(3001);
		FitResultAll->Add(loc_FitResult);
		
		FitResults.push_back({loc_template,loc_FitResult});
	}
	fitted_compton_fracE = sqrt(fitted_compton_fracE);
	
	h1->GetXaxis()->SetRangeUser(-8.0,8.0);
	
	yield   = n_total_exp * fitted_compton_frac;
	yieldE  = sqrt(pow(sqrt(n_total_exp)*fitted_compton_frac,2.0) + pow(n_total_exp*fitted_compton_fracE,2.0));
	
	yield  /= (loc_yield1/loc_yield2);
	yieldE /= (loc_yield1/loc_yield2);
	
	chi2 = fit->GetChisquare() / fit->GetNDF();
	
	//-----------------------------------------------------//
	// Calculate chi2 manually:
	
	chi2 = 0.;
	double n_bins = 0.;
	for(int ib=h1->FindBin(min_fit_x); ib<=h1_fit->FindBin(max_fit_x); ib++) {
		double loc_data_counts = h1_fit->GetBinContent(ib);
		double loc_fit_counts  = fit_result->GetBinContent(fit_result->FindBin(h1_fit->GetBinCenter(ib)));
		if(loc_fit_counts <= 0.) continue;
		chi2 += pow(loc_data_counts-loc_fit_counts, 2.0) / loc_fit_counts;
		n_bins += 1.0;
	}
	chi2 /= (n_bins/2);
	
	//-----------------------------------------------------//
	// Get the pull distribution:
	
	TH1F *h1_pull = (TH1F*)h1_fit->Clone(Form("h1_pull_%d_%d", tag_sys, counter));
	
	double  mean_pull = 0., n_pull = 0.;
	double stdev_pull = 0.;
	
	for(int ibin=1; ibin<=h1_pull->GetXaxis()->GetNbins(); ibin++) {
		double loc_x     = h1_fit->GetBinCenter(ibin);
		double loc_count = h1_fit->GetBinContent(ibin);
		double loc_fit   = FitResultAll->GetBinContent(FitResultAll->FindBin(loc_x));
		//double loc_fit   = fit_result->GetBinContent(fit_result->FindBin(loc_x));
		
		double loc_stat_unc  = pow(h1->GetBinError(ibin), 2.0);
		for(mc_hist_it = mc_hists.begin(); mc_hist_it != mc_hists.end(); mc_hist_it++) {
			loc_stat_unc += pow(
				mc_hist_it->second.second->GetBinError(mc_hist_it->second.second->FindBin(loc_x)), 2.0);
		}
		loc_stat_unc  = sqrt(loc_stat_unc);
		if(loc_stat_unc <= 1.) loc_stat_unc = 1.0;
		
		double loc_pull = (loc_count-loc_fit) / loc_stat_unc;
		
		mean_pull += loc_pull;
		n_pull    += 1.0;
		
		h1_pull->SetBinContent(ibin, loc_pull);
		h1_pull->SetBinError(ibin, 1.0);
	}
	mean_pull /= n_pull;
	
	for(int ibin=1; ibin<=h1_pull->GetXaxis()->GetNbins(); ibin++) {
		double loc_x     = h1_fit->GetBinCenter(ibin);
		double loc_count = h1_fit->GetBinContent(ibin);
		double loc_fit   = FitResultAll->GetBinContent(FitResultAll->FindBin(loc_x));
		//double loc_fit   = fit_result->GetBinContent(fit_result->FindBin(loc_x));
		
		double loc_stat_unc  = pow(h1->GetBinError(ibin), 2.0);
		for(mc_hist_it = mc_hists.begin(); mc_hist_it != mc_hists.end(); mc_hist_it++) {
			loc_stat_unc += pow(
				mc_hist_it->second.second->GetBinError(mc_hist_it->second.second->FindBin(loc_x)), 2.0);
		}
		loc_stat_unc  = sqrt(loc_stat_unc);
		if(loc_stat_unc <= 1.) loc_stat_unc = 1.0;
		
		double loc_pull = (loc_count-loc_fit) / loc_stat_unc;
		
		stdev_pull += pow(loc_pull-mean_pull, 2.0);
	}
	stdev_pull = sqrt(stdev_pull/(n_pull-1.0));
	
	if(tag_sys==0) {
		tagh_yieldfit_pull[counter-1] = {mean_pull, stdev_pull};
	} else {
		tagm_yieldfit_pull[counter-1] = {mean_pull, stdev_pull};
	}
	
	//----------------------------------------------------------------------//
	// Plot results:
	
	bool loc_DRAW_FITS = false;
	if(!tag_sys) {
		if(DRAW_FITS_TAGH) loc_DRAW_FITS = true;
	} else {
		if(DRAW_FITS_TAGM) loc_DRAW_FITS = true;
	}
	
	fit_result->GetYaxis()->SetTitle(Form("Counts / %d MeV", (int)(bin_size*1.e3)));
	fit_result->GetYaxis()->SetTitleOffset(1.8);
	fit_result->SetTitle("");
	
	FitResultAll->SetMinimum(0);
	FitResultAll->SetMaximum(FitResultAll->GetMaximum()*1.2);
	FitResultAll->SetLineColor(kRed);
	FitResultAll->SetLineWidth(2);
	FitResultAll->GetXaxis()->SetRangeUser(min_fit_x, max_fit_x);
	FitResultAll->GetXaxis()->SetLabelOffset(0.015);
	FitResultAll->GetYaxis()->SetTitle(Form("counts / %d MeV", (int)n_mev));
	
	TH1F *dataInput = (TH1F*)h1->Clone(Form("dataInput_%d_%d",tag_sys,counter));
	dataInput->SetMarkerStyle(21);
	dataInput->SetMarkerSize(0.7);
	dataInput->SetLineColor(1);
	dataInput->SetLineWidth(2);
	
	h1_pull->GetXaxis()->SetRangeUser(min_fit_x, max_fit_x);
	h1_pull->GetYaxis()->SetRangeUser(-5.5, 5.5);
	h1_pull->GetYaxis()->SetTitle("(Data - Fit) / #sigma_{Stat.}");
	h1_pull->GetYaxis()->SetTitleOffset(0.35);
	h1_pull->GetYaxis()->SetTitleSize(0.085);
	h1_pull->GetYaxis()->SetLabelSize(0.05);
	h1_pull->GetYaxis()->CenterTitle(true);
	h1_pull->GetXaxis()->SetTitleSize(0.1);
	h1_pull->GetXaxis()->SetTitleOffset(1.0);
	h1_pull->GetXaxis()->SetLabelSize(0.065);
	h1_pull->GetXaxis()->CenterTitle(true);
	h1_pull->GetXaxis()->SetTitle("#DeltaK = E_{Comp}#left(#theta_{1},#theta_{2}#right) - E_{#gamma} [GeV]");
	h1_pull->SetTitle("");
	
	if(loc_DRAW_FITS) {
		
		/*
		TH1F *h1_lin = (TH1F*)h1->Clone(Form("h1_lin_%d_%d", tag_sys, counter));
		h1_lin->SetMinimum(0.);
		
		canvas_fit_lin->cd();
		h1_lin->Draw("PE");
		fit_result->Draw("same");
		
		canvas_fit_log->cd();
		h1->Draw("PE");
		fit_result->Draw("same");
		
		TLegend *leg = new TLegend( 0.140, 0.545, 0.400, 0.855 );
		leg->AddEntry(h1_lin,     "Data",                 "PE");
		leg->AddEntry(fit_result, "Full Fit",              "l");
		
		for(mc_hist_it = mc_hists.begin(); mc_hist_it != mc_hists.end(); mc_hist_it++) {
			leg->AddEntry(
		}
		
		leg->AddEntry(h_comp_fit, "Compton Signal (fit)",  "l");
		leg->AddEntry(h_pair_fit, "e+e- Background (fit)", "l");
		leg->SetBorderSize(0);
		
		h_comp_fit->Draw("same hist");
		h_pair_fit->Draw("same hist");
		h1_lin->Draw("PE same");
		leg->Draw();
		
		h_comp_fit->Draw("same hist");
		h_pair_fit->Draw("same hist");
		h1->Draw("PE same");
		*/
		
		canvas_draw_top->cd();
		FitResultAll->Draw("hist");
		dataInput->Draw("PE1same");
		
		TLegend *leg1 = new TLegend(0.132,0.650,0.501,0.883);
		char text[200];
		leg1->SetFillColor(0);
		leg1->SetShadowColor(0);
		leg1->SetFillStyle(0);
		leg1->SetBorderSize(0);
		leg1->SetLineColor(0);
		sprintf(text,"Data");
		leg1->AddEntry(dataInput, text, "pl");
		sprintf(text,"Combined Fit");
		leg1->AddEntry(fit_result, text, "l");
		
		for(int imc=0; imc<(int)FitResults.size(); imc++) {
			FitResults[imc].second->Draw("same hist");
			sprintf(text, "%s = %.2f(%%)", 
				FitResults[imc].first.Data(), 1.e2*FitResults[imc].second->Integral()/fit_result->Integral());
			leg1->AddEntry(FitResults[imc].second, text, "f");
		}
		
		leg1->Draw();
		
		TLatex lat_new;
		lat_new.SetTextColor(kRed-6);
		lat_new.SetTextFont(52);
		lat_new.DrawLatexNDC(0.73, 0.85, Form("#scale[1.0]{E_{#gamma} = %.2f GeV}", eb));
		
		TLine *l0 = new TLine(min_fit_x,  0.0, max_fit_x,  0.0);
		TLine *l1 = new TLine(min_fit_x, -2.0, max_fit_x, -2.0);
		TLine *l2 = new TLine(min_fit_x,  2.0, max_fit_x,  2.0);
		l0->SetLineColor(kBlack);
		l1->SetLineColor(kRed);
		l2->SetLineColor(kRed);
		l0->SetLineWidth(2);
		l1->SetLineWidth(2);
		l2->SetLineWidth(2);
		
		canvas_draw_bot->cd();
		h1_pull->Draw("PE");
		l0->Draw("same");
		l1->Draw("same");
		l2->Draw("same");
		
		lat_new.DrawLatexNDC(0.115, 0.375, Form("#scale[1.5]{Mean of Pull Distribution: %.2f}", mean_pull));
		lat_new.DrawLatexNDC(0.115, 0.275, Form("#scale[1.5]{Std Dev. of Pull Distribution: %.2f}", stdev_pull));
	}
	
	//h_vertex->Delete();
	//h2_compton->Delete();
	//h2_pair->Delete();
	//h_pair_gen_flux->Delete();
	
	return fit_val;
}

int get_neighbor(int tag_sys, int counter, int sim_case) {
	
	char mc_dir[256];
	
	if(sim_case==0)      sprintf(mc_dir, "%s", comp_mc_dir.Data());
	else if(sim_case==1) sprintf(mc_dir, "%s", pair_mc_dir.Data());
	else if(sim_case==2) sprintf(mc_dir, "%s", triplet_mc_dir.Data());
	else return 0;
	
	int trial_counter = counter;
	
	if(tag_sys==0) {
		
		//if(counter>219) return 0;
		if(counter>127&&counter<179) return 0;
		else {
			
			int keep_going = 1;
			int loc_incr   = 1;
			
			while(keep_going) {
				
				int pos_incr;
				
				if(fabs(tagh_en[counter+loc_incr]-tagh_en[counter])
					< fabs(tagh_en[counter-loc_incr]-tagh_en[counter])) {
					
					pos_incr = 1;
					trial_counter = counter+loc_incr;
					
					if((trial_counter>127&&trial_counter<179) 
						|| trial_counter>219) {
						trial_counter = counter-loc_incr;
						if((trial_counter>127&&trial_counter<179) 
							|| trial_counter<1) {
							return 0;
						}
					}
					
				} else {
					
					pos_incr = 0;
					trial_counter = counter-loc_incr;
					
					if((trial_counter>127&&trial_counter<179) 
						|| trial_counter<1) {
						trial_counter = counter+loc_incr;
						if((trial_counter>127&&trial_counter<179) 
							|| trial_counter>219) {
							return 0;
						}
					}
				}
				
				char temp_fname[256];
				sprintf(temp_fname, "%s/tagh_%03d.root", mc_dir, 
					trial_counter);
				if(gSystem->AccessPathName(temp_fname)) {
					if(pos_incr) {
						trial_counter = counter-loc_incr;
					} else {
						trial_counter = counter+loc_incr;
					}
				}
				sprintf(temp_fname, "%s/tagh_%03d.root", mc_dir, 
					trial_counter);
				if(!gSystem->AccessPathName(temp_fname)) keep_going = 0;
				
				loc_incr++;
				if(loc_incr>5) keep_going = 0;
			}
			
			char temp_fname[256];
			sprintf(temp_fname, "%s/tagh_%03d.root", mc_dir, trial_counter);
			if(gSystem->AccessPathName(temp_fname)) return 0;
			
		}
	} else {
		if(counter>101) return 0;
		else if(counter<1) return 0;
		else {
			
			int keep_going = 1;
			int loc_incr = 1;
			while(keep_going) {
				
				int pos_incr;
				
				if(fabs(tagm_en[counter+loc_incr]-tagm_en[counter])
					< fabs(tagm_en[counter-loc_incr]-tagm_en[counter])) {
					
					pos_incr = 1;
					trial_counter = counter+loc_incr;
					
					if(trial_counter>101) {
						trial_counter = counter-loc_incr;
						if(trial_counter<1) {
							return 0;
						}
					}
					
				} else {
					
					pos_incr = 0;
					trial_counter = counter-loc_incr;
					
					if(trial_counter<1) {
						trial_counter = counter+loc_incr;
						if(trial_counter>101) {
							return 0;
						}
					}
				}
				
				char temp_fname[256];
				sprintf(temp_fname, "%s/tagm_%03d.root", mc_dir, 
					trial_counter);
				if(gSystem->AccessPathName(temp_fname)) {
					if(pos_incr) {
						trial_counter = counter-loc_incr;
					} else {
						trial_counter = counter+loc_incr;
					}
				}
				sprintf(temp_fname, "%s/tagm_%03d.root", comp_mc_dir.Data(), 
					trial_counter);
				if(!gSystem->AccessPathName(temp_fname)) keep_going = 0;
				
				loc_incr++;
				if(loc_incr>5) keep_going = 0;
			}
			
			char temp_fname[256];
			sprintf(temp_fname, "%s/tagm_%03d.root", mc_dir, trial_counter);
			if(gSystem->AccessPathName(temp_fname)) return 0;
		}
	}
	
	return trial_counter;
}
