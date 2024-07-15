#include "/work/halld/home/andrsmit/primex_compton_analysis/include/compton_inputs.h"
#include "/work/halld/home/andrsmit/primex_compton_analysis/include/get_phase1_energy_bin.h"

int get_neighbor(int tag_sys, int counter, int sim_case);

int fit_yield(int tag_sys, int counter, TH1F *h1, TH1F *h1_empty, double &yield, double &yieldE, double &chi2) {
	
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
	
	//for(int ibin=h1->GetXaxis()->FindBin(0.); ibin<=h1->GetXaxis()->GetNbins(); ibin++) {
	//	if(h1->GetBinContent(ibin)<=0.) {
	//		max_fit_x = h1->GetXaxis()->GetBinCenter(ibin-1);
	//		break;
	//	}
	//}
	
	double loc_deltaK_mu  = -1.05827e-01 + 8.60222e-02*eb - 1.84818e-02*pow(eb,2.0) + 8.04656e-04*pow(eb,3.0);
	double loc_deltaK_sig = -1.38713e-01 + 1.98660e-01*eb - 1.88848e-02*pow(eb,2.0) + 7.97340e-04*pow(eb,3.0);
	//min_fit_x = loc_deltaK_mu  - 5.0*loc_deltaK_sig;
	//max_fit_x = loc_deltaK_sig + 5.0*loc_deltaK_sig;
	
	if(max_fit_x > 4.25) max_fit_x = 4.25;
	if(eb<9.0 && min_fit_x < -6.0) min_fit_x = -6.0;
	
	h1->GetXaxis()->SetRangeUser(min_fit_x,max_fit_x);
	h1->SetLineColor(kBlack);
	h1->SetMarkerColor(kBlack);
	h1->SetMarkerStyle(8);
	h1->SetMarkerSize(0.6);
	h1->GetXaxis()->SetTitleSize(0.045);
	h1->GetXaxis()->SetTitleOffset(0.95);
	h1->GetXaxis()->SetLabelSize(0.0325);
	
	// Find the energy bin from phase I closest to the new energy bin:
	
	int old_tag_counter = 0;
	int old_tag_sys = 0;
	if(PHASE_VAL==1 && IS_BE_TARGET) {
		old_tag_counter = counter;
		old_tag_sys     = tag_sys;
	} else {
		get_phaseI_energy_bin(eb, old_tag_counter, old_tag_sys);
	}
	
	char loc_tag_name[256];
	if(!tag_sys) sprintf(loc_tag_name, "TAGH %d", counter);
	else         sprintf(loc_tag_name, "TAGM %d", counter);
	
	char old_tag_name[256];
	if(!old_tag_sys) sprintf(old_tag_name, "TAGH %d", old_tag_counter);
	else             sprintf(old_tag_name, "TAGM %d", old_tag_counter);
	
	cout << "current bin: " << loc_tag_name << "; phase I bin: " << old_tag_name << endl;
	
	if(old_tag_sys<0) return 0;
	
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
	double comp_flux_sim = n_comp_thrown / (comp_cs * n_e * mb);
	
	if(n_comp_thrown < 1.e3) return 0;
	TH2F *h2_compton;
	
	if(tag_sys==0) 
		h2_compton = (TH2F*)fSim->Get(Form("%s", hname_tagh_comp.Data()));
	else 
		h2_compton = (TH2F*)fSim->Get(Form("%s", hname_tagm_comp.Data()));
	
	TH1F *h1_compton = (TH1F*)h2_compton->ProjectionY(Form("h1_comp_%d_%d", tag_sys, counter));
	h1_compton->SetDirectory(0);
	
	fSim->Close();
	
	if(h1_compton->Integral() < 1.e1) return 0;
	
	//=================================================================================//
	// Get the Pair MC Histogram
	
	//cout << "getting pair mc" << endl;
	
	// histogram file:
	sprintf(fname, "%s/pair_rec.root", pair_mc_dir.Data());
	if(gSystem->AccessPathName(fname)) {
		cout << "Can't access e+e- pair mc histogram file." << endl;
		return 0;
	}
	TFile *fPair = new TFile(fname, "READ");
	
	TH2F *h2_pair;
	if(tag_sys==0) h2_pair = (TH2F*)fPair->Get(hname_tagh_pair.Data());
	else           h2_pair = (TH2F*)fPair->Get(hname_tagm_pair.Data());
	
	// generated flux file:
	sprintf(fname, "/work/halld/home/andrsmit/primex_compton_analysis/bhgen_test/recRootTrees/Run061321/sum.root");
	if(gSystem->AccessPathName(fname)) {
		cout << "Can't access e+e- pair mc generated flux file." << endl;
		return 0;
	}
	TFile *fPairFlux = new TFile(fname, "READ");
	
	TH1F *h_pair_gen_flux = (TH1F*)fPairFlux->Get("gen_flux");
	
	int n_bins_combine = 1;
	int min_counter = counter-n_bins_combine;
	int max_counter = counter+n_bins_combine;
	
	if(tag_sys==0) {
		if(min_counter<1)   min_counter = 1;
		else if(min_counter>127 && min_counter<179) min_counter = 179;
		
		if(max_counter>221) max_counter = 221;
		else if(max_counter>127 && max_counter<179) max_counter = 127;
	} else {
		if(min_counter<1)   min_counter = 1;
		if(max_counter>102) max_counter = 102;
	}
	
	TH1F *h1_pair = (TH1F*)h2_pair->ProjectionY(Form("h1_pair_%d_%d", tag_sys, counter), min_counter, max_counter);
	h1_pair->SetDirectory(0);
	
	double loc_pair_gen_flux = 0.;
	for(int loc_counter = min_counter; loc_counter <= max_counter; loc_counter++) {
		double loc_counter_eb;
		if(tag_sys==0) loc_counter_eb = tagh_en[loc_counter-1];
		else           loc_counter_eb = tagm_en[loc_counter-1];
		double loc_counter_flux = h_pair_gen_flux->GetBinContent(h_pair_gen_flux->FindBin(loc_counter_eb));
		loc_pair_gen_flux += loc_counter_flux;
	}
	
	double loc_pair_cs   = f_pair_cs->Eval(eb) + f_triplet_cs->Eval(eb);
	double loc_pair_flux = loc_pair_gen_flux / ((n_e/n_Z) * mb * loc_pair_cs);
	
	h1_pair->Scale(0.70);
	
	fPair->Close();
	fPairFlux->Close();
	
	//=================================================================================//
	
	h1_empty->Rebin(rebins);
	h1_empty->GetXaxis()->SetTitle("E_{Comp} - E_{#gamma} [GeV]");
	h1_empty->GetXaxis()->SetTitleOffset(1.1);
	h1_empty->GetYaxis()->SetTitle(Form("counts / %d MeV", n_mev));
	h1_empty->GetYaxis()->SetTitleOffset(1.3);
	h1_empty->SetLineColor(kMagenta);
	h1_empty->SetFillColor(kMagenta);
	h1_empty->SetFillStyle(3004);
	h1_empty->SetLineWidth(2);
	
	h1_compton->Rebin(rebins);
	h1_compton->Scale(loc_flux / comp_flux_sim);
	h1_compton->GetXaxis()->SetTitle("E_{Comp} - E_{#gamma} [GeV]");
	h1_compton->GetXaxis()->SetTitleOffset(1.1);
	h1_compton->GetYaxis()->SetTitle( Form("counts / %d MeV", n_mev));
	h1_compton->GetYaxis()->SetTitleOffset(1.3);
	h1_compton->SetLineColor(kCyan);
	h1_compton->SetFillColor(kCyan);
	h1_compton->SetFillStyle(3004);
	h1_compton->SetLineWidth(2);
	
	h1_pair->Rebin(rebins);
	h1_pair->Scale(loc_flux / loc_pair_flux);
	h1_pair->GetXaxis()->SetTitle("E_{Comp} - E_{#gamma} [GeV]");
	h1_pair->GetXaxis()->SetTitleOffset(1.1);
	h1_pair->GetYaxis()->SetTitle(Form("counts / %d MeV", n_mev));
	h1_pair->GetYaxis()->SetTitleOffset(1.3);
	h1_pair->SetLineColor(kGreen);
	h1_pair->SetFillColor(kGreen);
	h1_pair->SetFillStyle(3004);
	h1_pair->SetLineWidth(2);
	
	for(int ib=1; ib<=h1->GetXaxis()->GetNbins(); ib++) {
		if(h1->GetBinContent(ib)<0.) {
			h1->SetBinContent(ib,0.1);
		}
	}
	for(int ib=1; ib<=h1_empty->GetXaxis()->GetNbins(); ib++) {
		if(h1_empty->GetBinContent(ib)<0.) {
			h1_empty->SetBinContent(ib,0.1);
		}
	}
	for(int ib=1; ib<=h1_compton->GetXaxis()->GetNbins(); ib++) {
		if(h1_compton->GetBinContent(ib)<0.) {
			h1_compton->SetBinContent(ib,0.1);
		}
	}
	for(int ib=1; ib<=h1_pair->GetXaxis()->GetNbins(); ib++) {
		if(h1_pair->GetBinContent(ib)<0.) {
			h1_pair->SetBinContent(ib,0.1);
		}
	}
	
	double min_counts = 5.;
	if(eb<7.0) min_counts=10.;
	// loop through bins and if there are any where all 3 distributions have zero counts, set that to the max range:
	for(int ib=h1->FindBin(0.); ib<h1->FindBin(max_fit_x); ib++) {
		if(h1_empty->GetBinContent(ib)<=min_counts && h1_compton->GetBinContent(ib)<=min_counts && h1_pair->GetBinContent(ib)<=min_counts) {
			max_fit_x = h1->GetXaxis()->GetBinCenter(ib);
			break;
		}
	}
	for(int ib=h1->FindBin(0.); ib>h1->FindBin(min_fit_x); ib--) {
		if(h1_empty->GetBinContent(ib)<=min_counts && h1_compton->GetBinContent(ib)<=min_counts && h1_pair->GetBinContent(ib)<=min_counts) {
			min_fit_x = h1->GetXaxis()->GetBinCenter(ib);
			break;
		}
	}
	
		
	h1->GetXaxis()->SetRangeUser(min_fit_x, max_fit_x);
	h1_empty->GetXaxis()->SetRangeUser(min_fit_x, max_fit_x);
	h1_compton->GetXaxis()->SetRangeUser(min_fit_x, max_fit_x);
	h1_pair->GetXaxis()->SetRangeUser(min_fit_x, max_fit_x);
	
	double n_total_exp   = h1->Integral(h1->FindBin(min_fit_x), h1->FindBin(max_fit_x));
	double n_empty_exp   = h1_empty->Integral(h1_empty->FindBin(min_fit_x), h1_empty->FindBin(max_fit_x));
	double n_compton_exp = h1_compton->Integral(h1_compton->FindBin(min_fit_x), h1_compton->FindBin(max_fit_x));
	double n_pair_exp    = h1_pair->Integral(h1_pair->FindBin(min_fit_x), h1_pair->FindBin(max_fit_x));
	
	double empty_frac_exp = n_empty_exp / n_total_exp;
	double comp_frac_exp  = n_compton_exp / n_total_exp;
	double pair_frac_exp  = n_pair_exp / n_total_exp;
	
	//-----   Perform Binned Maximum Likelihood Fit (TFractionFitter)   -----//
	
	if(DEBUG_FITS) {
		if(c_debug==NULL) {
			c_debug = new TCanvas("c_debug", "Fit Debug", 800, 800);
			c_debug->SetTickx(); c_debug->SetTicky();
		}
		TH1F *h1_debug = (TH1F*)h1->Clone("h1_debug");
		h1_debug->SetLineColor(kRed);
		h1_debug->Add(h1_empty);
		h1_debug->Add(h1_compton);
		h1_debug->Add(h1_pair);
		h1_debug->Add(h1,-1.0);
		h1_debug->SetLineColor(kRed);
		h1_debug->SetLineWidth(2);
		
		c_debug->cd();
		h1->Draw("PE");
		h1_empty->Draw("same hist");
		h1_compton->Draw("same hist");
		h1_pair->Draw("same hist");
		h1_debug->Draw("same hist");
		c_debug->Update();
	}
	
	TObjArray *mc = new TObjArray(3);
	mc->Add(h1_compton);
	mc->Add(h1_pair);
	mc->Add(h1_empty);
	
	TFractionFitter *fit = new TFractionFitter(h1, mc);
	
	fit->Constrain(0, 0.500,               1.000);
	fit->Constrain(1, 0.90*pair_frac_exp,  1.10*pair_frac_exp);
	//fit->Constrain(1, 0.0, 1.0);
	fit->Constrain(2, 0.95*empty_frac_exp, 1.05*empty_frac_exp);
	
	//ROOT::Fit::Fitter *vfit = (ROOT::Fit::Fitter*)fit->GetFitter();
	//vfit->FitConfig.SetParamsSettings(comp_frac_exp, pair_frac_exp, empty_frac_exp);
	
	fit->SetRangeX(h1->FindBin(min_fit_x),h1->FindBin(max_fit_x));
	
	Int_t fit_status = fit->Fit();
	TH1F *fit_result = (TH1F*)fit->GetPlot();
	
	double comp_frac,  comp_fracE;
	double pair_frac,  pair_fracE;
	double empty_frac, empty_fracE;
	fit->GetResult(0, comp_frac, comp_fracE);
	fit->GetResult(1, pair_frac, pair_fracE);
	fit->GetResult(2, empty_frac, empty_fracE);
	
	//yield  = Ncomp;
	//yield  = n_total_exp * (1.0 - pair_frac - empty_frac);
	
	h1->GetXaxis()->SetRangeUser(-8.0,8.0);
	//TH1F *h1_comp_test = (TH1F*)fit->GetMCPrediction(0);
	//double loc_yield = h1_comp_test->Integral();
	
	yield = n_total_exp * (1.0 - pair_frac - empty_frac);
	
	double loc_yield1 = h1_compton->Integral(h1_compton->FindBin(min_fit_x), h1_compton->FindBin(max_fit_x));
	double loc_yield2 = h1_compton->Integral(h1_compton->FindBin(-8.0),      h1_compton->FindBin(8.0));
	
	yield *= (loc_yield1/loc_yield2);
	yieldE = sqrt(yield);
	
	chi2 = fit->GetChisquare() / fit->GetNDF();
	/*
	chi2 = 0.;
	double n_bins = 0.;
	for(int ib=h1->FindBin(min_fit_x); ib<=h1->FindBin(max_fit_x); ib++) {
		double loc_data_counts = h1->GetBinContent(ib);
		double loc_fit_counts  = fit_result->GetBinContent(
			fit_result->FindBin(h1->GetBinCenter(ib)));
		chi2 += pow(loc_data_counts-loc_fit_counts, 2.0) / loc_fit_counts;
		n_bins += 1.0;
	}
	chi2 /= n_bins;
	*/
	if(tag_sys==0) {
		tagh_compton_fraction[counter-1]   = comp_frac;
		tagh_compton_fractionE[counter-1]  = comp_fracE;
		tagh_pair_fraction[counter-1]      = pair_frac;
		tagh_pair_fractionE[counter-1]     = pair_fracE;
		tagh_pair_fraction_exp[counter-1]  = pair_frac_exp;
		tagh_empty_fraction[counter-1]     = empty_frac;
		tagh_empty_fractionE[counter-1]    = empty_fracE;
		tagh_empty_fraction_exp[counter-1] = empty_frac_exp;
	} else {
		tagm_compton_fraction[counter-1]   = comp_frac;
		tagm_compton_fractionE[counter-1]  = comp_fracE;
		tagm_pair_fraction[counter-1]      = pair_frac;
		tagm_pair_fractionE[counter-1]     = pair_fracE;
		tagm_pair_fraction_exp[counter-1]  = pair_frac_exp;
		tagm_empty_fraction[counter-1]     = empty_frac;
		tagm_empty_fractionE[counter-1]    = empty_fracE;
		tagm_empty_fraction_exp[counter-1] = empty_frac_exp;
	}
	
	//----------------------------------------------------------------------//
	
	bool loc_DRAW_FITS = false;
	if(!tag_sys) {
		if(DRAW_FITS_TAGH) loc_DRAW_FITS = true;
	} else {
		if(DRAW_FITS_TAGM) loc_DRAW_FITS = true;
	}
	
	if(loc_DRAW_FITS) {
		
		fit_result->SetLineColor(kRed+1);
		fit_result->GetXaxis()->SetRangeUser(min_fit_x,max_fit_x);
		
		double Nmc   = fit_result->Integral();
		double Ndata = h1->Integral();
		
		TH1F *h_comp_fit, *h_pair_fit, *h_empty_fit;
		
		h_comp_fit = (TH1F*)fit->GetMCPrediction(0);
		h_comp_fit->GetXaxis()->SetRangeUser(min_fit_x, max_fit_x);
		h_pair_fit = (TH1F*)fit->GetMCPrediction(1);
		h_pair_fit->GetXaxis()->SetRangeUser(min_fit_x, max_fit_x);
		h_empty_fit = (TH1F*)fit->GetMCPrediction(2);
		h_empty_fit->GetXaxis()->SetRangeUser(min_fit_x, max_fit_x);
		
		h_comp_fit->Scale(comp_frac*Ndata/h_comp_fit->Integral());
		double Ncomp = h_comp_fit->Integral();
		h_comp_fit->SetLineColor(kCyan);
		h_comp_fit->SetMarkerColor(kCyan);
		h_comp_fit->SetFillColor(kCyan);
		h_comp_fit->SetFillStyle(3004);
		h_comp_fit->SetLineWidth(2);
		h_comp_fit->SetMarkerStyle(8);
		h_comp_fit->SetMarkerSize(0.6);
		
		h_pair_fit->Scale(pair_frac*Ndata/h_pair_fit->Integral());
		double Npair = h_pair_fit->Integral();
		h_pair_fit->SetLineColor(kGreen);
		h_pair_fit->SetMarkerColor(kGreen);
		h_pair_fit->SetFillColor(kGreen);
		h_pair_fit->SetFillStyle(3005);
		h_pair_fit->SetLineWidth(2);
		h_pair_fit->SetMarkerStyle(8);
		h_pair_fit->SetMarkerSize(0.6);
		
		h_empty_fit->Scale(empty_frac*Ndata/h_empty_fit->Integral());
		double Nempty = h_empty_fit->Integral();
		h_empty_fit->SetLineColor(kMagenta);
		h_empty_fit->SetMarkerColor(kMagenta);
		h_empty_fit->SetFillColor(kGreen);
		h_empty_fit->SetFillStyle(3005);
		h_empty_fit->SetLineWidth(2);
		h_empty_fit->SetMarkerStyle(8);
		h_empty_fit->SetMarkerSize(0.6);
		
		TH1F *h1_lin = (TH1F*)h1->Clone(Form("h1_lin_%d_%d", tag_sys, counter));
		h1_lin->SetMinimum(0.);
		
		TLegend *leg = new TLegend( 0.140, 0.545, 0.400, 0.855 );
		leg->AddEntry(h1_lin,     "Data",                 "PE");
		leg->AddEntry(fit_result, "Full Fit",              "l");
		leg->AddEntry(h_comp_fit, "Compton Signal (fit)",  "l");
		leg->AddEntry(h_pair_fit, "e+e- Background (fit)", "l");
		leg->AddEntry(h_empty_fit, "Empty-Target Background (fit)", "l");
		leg->SetBorderSize(0);
		
		canvas_fit_lin->cd();
		h1_lin->Draw("PE");
		fit_result->Draw("same");
		h_comp_fit->Draw("same hist");
		h_pair_fit->Draw("same hist");
		h_empty_fit->Draw("same hist");
		h1_lin->Draw("PE same");
		leg->Draw();
		
		canvas_fit_log->cd();
		h1->Draw("PE");
		fit_result->Draw("same");
		h_comp_fit->Draw("same hist");
		h_pair_fit->Draw("same hist");
		h_empty_fit->Draw("same hist");
		h1->Draw("PE same");
		
		fit_result->SetMaximum(fit_result->GetMaximum()*1.4);
		
		fit_result->GetYaxis()->SetTitle(Form("Counts / %d MeV", (int)(bin_size*1.e3)));
		fit_result->GetYaxis()->SetTitleOffset(1.8);
		fit_result->SetTitle("");
		
		canvas_draw->cd();
		fit_result->Draw("hist");
		//FitResultAll->Draw("hist same");
		
		TH1F *FitResultSignal = (TH1F*)h1_compton->Clone(Form("FitResultSignal_%d_%d", tag_sys, counter));
		FitResultSignal->Scale(1./FitResultSignal->Integral()*comp_frac*Ndata);
		TH1F *FitResultBkgd   = (TH1F*)h1_pair->Clone(Form("FitResultBkgd_%d_%d", tag_sys, counter));
		FitResultBkgd->Scale(1./FitResultBkgd->Integral()*pair_frac*Ndata);
		TH1F *FitResultEmpty  = (TH1F*)h1_empty->Clone(Form("FitResultEmpty_%d_%d", tag_sys, counter));
		FitResultEmpty->Scale(1./FitResultEmpty->Integral()*empty_frac*Ndata);
		
		TH1F *FitResultAll = (TH1F*)h1->Clone(Form("FitResultAll_%d_%d", tag_sys, counter));
		FitResultAll->Reset();
		FitResultAll->Add(FitResultSignal);
		FitResultAll->Add(FitResultBkgd);
		FitResultAll->SetMinimum(0);
		FitResultAll->SetMaximum(FitResultAll->GetMaximum()*1.4);
		FitResultAll->SetLineColor(kRed);
		FitResultAll->SetLineWidth(2);
		
		TH1F *dataInput = (TH1F*)h1->Clone(Form("dataInput_%d_%d",tag_sys,counter));
		dataInput->SetMarkerStyle(21);
		dataInput->SetMarkerSize(0.7);
		dataInput->SetLineColor(1);
		dataInput->SetLineWidth(2);
		dataInput->Draw("PE1same");
		FitResultSignal->SetLineColor(kBlue);
		FitResultSignal->SetFillColor(kBlue);
		FitResultSignal->SetFillStyle(3002);
		FitResultSignal->Draw("same hist");
		FitResultBkgd->SetLineColor(kGreen);
		FitResultBkgd->SetFillColor(kGreen);
		FitResultBkgd->SetFillStyle(3004);
		FitResultBkgd->Draw("same hist");
		FitResultEmpty->SetLineColor(kMagenta);
		FitResultEmpty->SetFillColor(kMagenta);
		FitResultEmpty->SetFillStyle(3004);
		FitResultEmpty->Draw("same hist");
   	
		TLegend *leg1 = new TLegend(0.185,0.673,0.554,0.873);
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
		/*
		sprintf(text,"Data: %5.1f events", dataInput->Integral());
		leg1->AddEntry(dataInput, text, "pl");
		sprintf(text,"Fitted: %5.1f events", fit_result->Integral());
		leg1->AddEntry(fit_result, text, "l");
		*/
		sprintf(text, "Signal = %.2f(%%)", 100.*FitResultSignal->Integral()/fit_result->Integral());
		leg1->AddEntry(FitResultSignal, text, "f");
		sprintf(text, "Bkgd = %.2f(%%)", 100.*FitResultBkgd->Integral()/fit_result->Integral());
		leg1->AddEntry(FitResultBkgd, text, "f");
		sprintf(text, "Empty = %.2f(%%)", 100.*FitResultEmpty->Integral()/fit_result->Integral());
		leg1->AddEntry(FitResultEmpty, text, "f");
		leg1->Draw();
		
		TLatex lat_new;
		lat_new.SetTextColor(kRed-6);
		lat_new.SetTextFont(52);
		lat_new.DrawLatexNDC(0.63, 0.83, Form("#scale[1.0]{E_{#gamma} = %.2f GeV}", eb));
		
		canvas_fit->Update();
		canvas_draw->Update();
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
	else if(sim_case==2) sprintf(mc_dir, "%s", trip_mc_dir.Data());
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
