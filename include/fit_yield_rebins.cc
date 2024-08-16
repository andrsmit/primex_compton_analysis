#include "compton_cs.h"
#include "get_phase1_energy_bin.cc"
#include "fit_yield.cc"

int fit_yield_rebins(int tag_sys, int counter, vector<int> bin_vec, double eb, double loc_flux, TH1F *h1, 
	double &yield, double &yieldE, double &chi2) {
	
	int fit_val = 1;
	
	char tag_sys_char[256];
	if(tag_sys==0) sprintf(tag_sys_char, "tagh");
	else           sprintf(tag_sys_char, "tagm");
	
	double loc_acc = f_acc->Eval(eb);
	
	double min_fit_x = h1->GetBinCenter(h1->FindFirstBinAbove()+2);
	double max_fit_x = h1->GetBinCenter(h1->FindLastBinAbove()-2);
	
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
	
	int n_bins = (int)bin_vec.size();
	char fname[256];
	
	//=================================================================================//
	// Get the Compton MC Histogram
	
	//cout << "getting compton mc" << endl;
	
	TH1F *h1_compton = new TH1F(Form("deltaK_comp_mc_%d", counter), "#DeltaK; [GeV]", 1000, -8.0, 8.0);
	
	double comp_flux_sim = 0.;
	
	for(int ibin=0; ibin<n_bins; ibin++) {
		
		int loc_counter = bin_vec[ibin];
		
		if(tag_sys==0) 
			sprintf(fname, "%s/tagh_%03d.root", comp_mc_dir.Data(), loc_counter);
		else 
			sprintf(fname, "%s/tagm_%03d.root", comp_mc_dir.Data(), loc_counter);
		
		
		int neighbor_counter = counter;
		if(gSystem->AccessPathName(fname)) {
			
			//return 0;
			neighbor_counter = get_neighbor(tag_sys, loc_counter, 0);
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
		
		double comp_cs = f_theory->Eval(eb);
		comp_flux_sim += n_comp_thrown / (comp_cs * n_e * mb);
		
		if(n_comp_thrown < 1.e3) return 0;
		
		TH2F *h2_compton;
		
		if(tag_sys==0) 
			h2_compton = (TH2F*)fSim->Get(Form("%s", hname_tagh_comp.Data()));
		else 
			h2_compton = (TH2F*)fSim->Get(Form("%s", hname_tagm_comp.Data()));
		
		h1_compton->Add((TH1F*)h2_compton->ProjectionY(Form("h1_comp_%d_%d", tag_sys, loc_counter)));
		h1_compton->SetDirectory(0);
		
		fSim->Close();
	}
	
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
	
	// generated flux file:
	sprintf(fname, "/work/halld/home/andrsmit/primex_compton_analysis/bhgen_test/recRootTrees/Run061321/sum.root");
	if(gSystem->AccessPathName(fname)) {
		cout << "Can't access e+e- pair mc generated flux file." << endl;
		return 0;
	}
	TFile *fPairFlux = new TFile(fname, "READ");
	TH1F *h_pair_gen_flux = (TH1F*)fPairFlux->Get("gen_flux");
	
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
	
	
	TH1F *h1_pair;
	double loc_pair_gen_flux = 0.;
	
	int n_bins_combine = 10;
	int min_counter = phase1_tag_counter-n_bins_combine;
	int max_counter = phase1_tag_counter+n_bins_combine;
	
	if(phase1_tag_sys==0) {
		if(min_counter<1)   min_counter = 1;
		else if(min_counter>127 && min_counter<179) min_counter = 179;
		
		if(max_counter>221) max_counter = 221;
		else if(max_counter>127 && max_counter<179) max_counter = 127;
		
		TH2F *h2_pair = (TH2F*)fPair->Get(hname_tagh_pair.Data());
		h1_pair = (TH1F*)h2_pair->ProjectionY(Form("h1_pair_%d_%d", tag_sys, counter), min_counter, max_counter);
		
		for(int loc_counter = min_counter; loc_counter <= max_counter; loc_counter++) {
			double loc_counter_eb   = tagh_en_phase1[loc_counter-1];
			double loc_counter_flux = h_pair_gen_flux->GetBinContent(h_pair_gen_flux->FindBin(loc_counter_eb));
			loc_pair_gen_flux += loc_counter_flux;
		}
		
	} else {
		if(min_counter<1)   min_counter = 1;
		if(max_counter>102) max_counter = 102;
		
		TH2F *h2_pair = (TH2F*)fPair->Get(hname_tagm_pair.Data());
		h1_pair = (TH1F*)h2_pair->ProjectionY(Form("h1_pair_%d_%d", tag_sys, counter), min_counter, max_counter);
		
		for(int loc_counter = min_counter; loc_counter <= max_counter; loc_counter++) {
			double loc_counter_eb   = tagm_en_phase1[loc_counter-1];
			double loc_counter_flux = h_pair_gen_flux->GetBinContent(h_pair_gen_flux->FindBin(loc_counter_eb));
			loc_pair_gen_flux += loc_counter_flux;
		}
	}
	h1_pair->SetDirectory(0);
	
	double loc_pair_cs   = f_pair_cs->Eval(eb) + f_triplet_cs->Eval(eb);
	double loc_pair_flux = loc_pair_gen_flux / ((n_e/n_Z) * mb * loc_pair_cs);
	
	fPair->Close();
	fPairFlux->Close();
	
	//=================================================================================//
	
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
			h1->SetBinContent(ib,0.);
		}
	}
	for(int ib=1; ib<=h1_compton->GetXaxis()->GetNbins(); ib++) {
		if(h1_compton->GetBinContent(ib)<0.) {
			h1_compton->SetBinContent(ib,0.);
		}
	}
	for(int ib=1; ib<=h1_pair->GetXaxis()->GetNbins(); ib++) {
		if(h1_pair->GetBinContent(ib)<0.) {
			h1_pair->SetBinContent(ib,0.);
		}
	}
	
	double min_counts = 0.01*h1->GetMaximum();
	//if(eb<9.0) min_counts=10.;
	// loop through bins and if there are any where all mc distributions have zero counts, set that to the max range:
	for(int ib=h1->FindBin(0.); ib<h1->FindBin(max_fit_x); ib++) {
		if(h1_compton->GetBinContent(ib)<=min_counts && h1_pair->GetBinContent(ib)<=min_counts) {
			max_fit_x = h1->GetXaxis()->GetBinCenter(ib-1);
			break;
		}
	}
	for(int ib=h1->FindBin(0.); ib>h1->FindBin(min_fit_x); ib--) {
		if(h1_compton->GetBinContent(ib)<=min_counts && h1_pair->GetBinContent(ib)<=min_counts) {
			min_fit_x = h1->GetXaxis()->GetBinCenter(ib+1);
			break;
		}
	}
	
	min_fit_x = loc_deltaK_mu - 5.0*loc_deltaK_sig;
	double temp_max_fit_x = loc_deltaK_mu + 5.0*loc_deltaK_sig;
	if(temp_max_fit_x < max_fit_x) max_fit_x = temp_max_fit_x;
	
	min_fit_x = -7.0;
	
	//if(max_fit_x < 2.0) max_fit_x = 2.0;
	
	h1->GetXaxis()->SetRangeUser(min_fit_x, max_fit_x);
	h1_compton->GetXaxis()->SetRangeUser(min_fit_x, max_fit_x);
	h1_pair->GetXaxis()->SetRangeUser(min_fit_x, max_fit_x);
	
	double n_total_exp   = h1->Integral(h1->FindBin(min_fit_x), h1->FindBin(max_fit_x));
	double n_compton_exp = h1_compton->Integral(h1_compton->FindBin(min_fit_x), h1_compton->FindBin(max_fit_x));
	double n_pair_exp    = h1_pair->Integral(h1_pair->FindBin(min_fit_x), h1_pair->FindBin(max_fit_x));
	
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
		h1_debug->Add(h1_compton);
		h1_debug->Add(h1_pair);
		h1_debug->Add(h1,-1.0);
		h1_debug->SetLineColor(kRed);
		h1_debug->SetLineWidth(2);
		
		c_debug->cd();
		h1->Draw("PE");
		h1_compton->Draw("same hist");
		h1_pair->Draw("same hist");
		h1_debug->Draw("same hist");
		c_debug->Update();
	}
	
	double loc_yield1 = h1_compton->Integral(h1_compton->FindBin(min_fit_x), h1_compton->FindBin(max_fit_x));
	double loc_yield2 = h1_compton->Integral(h1_compton->FindBin(-8.0),      h1_compton->FindBin(8.0));
	
	// if fraction of pairs is low enough, don't bother with fit:
	if(pair_frac_exp<0.015) {
		yield  = n_total_exp * (1.0 - pair_frac_exp);
		yield *= (loc_yield1/loc_yield2);
		yieldE = sqrt(yield);
		
		chi2 = 1.0;
		return fit_val;
	}
	
	
	TObjArray *mc = new TObjArray(2);
	mc->Add(h1_compton);
	mc->Add(h1_pair);
	
	TFractionFitter *fit = new TFractionFitter(h1, mc);
	
	fit->Constrain(0, 0.000, 1.000);
	fit->Constrain(1, 0.000, 1.000);
	//fit->Constrain(1, 0.5*pair_frac_exp, 1.5*pair_frac_exp);
	
	fit->SetRangeX(h1->FindBin(min_fit_x),h1->FindBin(max_fit_x));
	
	Int_t fit_status = fit->Fit();
	TH1F *fit_result = (TH1F*)fit->GetPlot();
	
	double comp_frac,  comp_fracE;
	double pair_frac,  pair_fracE;
	fit->GetResult(0, comp_frac, comp_fracE);
	fit->GetResult(1, pair_frac, pair_fracE);
	
	h1->GetXaxis()->SetRangeUser(-8.0,8.0);
	
	yield  = n_total_exp * (1.0 - pair_frac);
	yield /= (loc_yield1/loc_yield2);
	yieldE = sqrt(yield);
	
	yieldE  = sqrt(pow(sqrt(n_total_exp)*(1.0 - pair_frac),2.0) + pow(n_total_exp*pair_fracE,2.0));
	yieldE *= (loc_yield1/loc_yield2);
	
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
		
		double Nmc   = fit_result->Integral(fit_result->FindBin(min_fit_x), fit_result->FindBin(max_fit_x));
		double Ndata = h1->Integral(h1->FindBin(min_fit_x), h1->FindBin(max_fit_x));
		
		TH1F *h_comp_fit, *h_pair_fit;
		
		h_comp_fit = (TH1F*)fit->GetMCPrediction(0);
		h_comp_fit->GetXaxis()->SetRangeUser(min_fit_x, max_fit_x);
		h_pair_fit = (TH1F*)fit->GetMCPrediction(1);
		h_pair_fit->GetXaxis()->SetRangeUser(min_fit_x, max_fit_x);
		
		h_comp_fit->Scale(comp_frac*Ndata/h_comp_fit->Integral(h_comp_fit->FindBin(min_fit_x), h_comp_fit->FindBin(max_fit_x)));
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
		
		TH1F *h1_lin = (TH1F*)h1->Clone(Form("h1_lin_%d_%d", tag_sys, counter));
		h1_lin->SetMinimum(0.);
		
		TLegend *leg = new TLegend( 0.140, 0.545, 0.400, 0.855 );
		leg->AddEntry(h1_lin,     "Data",                 "PE");
		leg->AddEntry(fit_result, "Full Fit",              "l");
		leg->AddEntry(h_comp_fit, "Compton Signal (fit)",  "l");
		leg->AddEntry(h_pair_fit, "e+e- Background (fit)", "l");
		leg->SetBorderSize(0);
		
		canvas_fit_lin->cd();
		h1_lin->Draw("PE");
		fit_result->Draw("same");
		h_comp_fit->Draw("same hist");
		h_pair_fit->Draw("same hist");
		h1_lin->Draw("PE same");
		leg->Draw();
		
		canvas_fit_log->cd();
		h1->Draw("PE");
		fit_result->Draw("same");
		h_comp_fit->Draw("same hist");
		h_pair_fit->Draw("same hist");
		h1->Draw("PE same");
		
		fit_result->SetMaximum(fit_result->GetMaximum()*1.4);
		
		fit_result->GetYaxis()->SetTitle(Form("Counts / %d MeV", (int)(bin_size*1.e3)));
		fit_result->GetYaxis()->SetTitleOffset(1.8);
		fit_result->SetTitle("");
		
		canvas_draw_top->cd();
		fit_result->Draw("hist");
		
		TH1F *FitResultSignal = (TH1F*)h1_compton->Clone(Form("FitResultSignal_%d_%d", tag_sys, counter));
		FitResultSignal->Scale(1./FitResultSignal->Integral()*comp_frac*Ndata);
		TH1F *FitResultBkgd   = (TH1F*)h1_pair->Clone(Form("FitResultBkgd_%d_%d", tag_sys, counter));
		FitResultBkgd->Scale(1./FitResultBkgd->Integral()*pair_frac*Ndata);
		
		TH1F *FitResultAll = (TH1F*)h1->Clone(Form("FitResultAll_%d_%d", tag_sys, counter));
		FitResultAll->Reset();
		FitResultAll->Add(FitResultSignal);
		FitResultAll->Add(FitResultBkgd);
		FitResultAll->SetMinimum(0);
		FitResultAll->SetMaximum(FitResultAll->GetMaximum()*1.2);
		FitResultAll->SetLineColor(kRed);
		FitResultAll->SetLineWidth(2);
		FitResultAll->GetXaxis()->SetRangeUser(min_fit_x, max_fit_x);
		FitResultAll->GetXaxis()->SetLabelOffset(0.015);
		FitResultAll->GetYaxis()->SetTitle(Form("counts / %d MeV", (int)n_mev));
		FitResultAll->Draw("hist");
		
		
		TH1F *dataInput = (TH1F*)h1->Clone(Form("dataInput_%d_%d",tag_sys,counter));
		dataInput->SetMarkerStyle(21);
		dataInput->SetMarkerSize(0.7);
		dataInput->SetLineColor(1);
		dataInput->SetLineWidth(2);
		dataInput->Draw("PE1same");
		FitResultSignal->SetLineColor(kCyan+1);
		FitResultSignal->SetFillColor(kCyan-9);
		FitResultSignal->SetFillStyle(3001);
		FitResultSignal->Draw("same hist");
		FitResultBkgd->SetLineColor(kGreen+2);
		FitResultBkgd->SetFillColor(kGreen-10);
		FitResultBkgd->SetFillStyle(3001);
		FitResultBkgd->Draw("same hist");
   	
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
		/*
		sprintf(text,"Data: %5.1f events", dataInput->Integral());
		leg1->AddEntry(dataInput, text, "pl");
		sprintf(text,"Fitted: %5.1f events", fit_result->Integral());
		leg1->AddEntry(fit_result, text, "l");
		*/
		sprintf(text, "Signal = %.2f(%%)", 100.*FitResultSignal->Integral()/fit_result->Integral());
		leg1->AddEntry(FitResultSignal, text, "f");
		sprintf(text, "e^{+}e^{-} Bkgd = %.2f(%%)", 100.*FitResultBkgd->Integral()/fit_result->Integral());
		leg1->AddEntry(FitResultBkgd, text, "f");
		leg1->Draw();
		
		TLatex lat_new;
		lat_new.SetTextColor(kRed-6);
		lat_new.SetTextFont(52);
		lat_new.DrawLatexNDC(0.73, 0.85, Form("#scale[1.0]{E_{#gamma} = %.2f GeV}", eb));
		
		// Plot the pull distribution:
		
		TH1F *h1_pull      = (TH1F*)h1_lin->Clone("h1_pull");
		TH1F *h1_pull_dist = new TH1F("h1_pull_dist", "Pull Distrtibution; [#sigma]", 50, -10., 10.);
		
		for(int ibin=1; ibin<=h1_pull->GetXaxis()->GetNbins(); ibin++) {
			double loc_x = h1_pull->GetBinCenter(ibin);
			double loc_count = h1_pull->GetBinContent(ibin);
			double loc_fit   = FitResultSignal->GetBinContent(FitResultSignal->FindBin(loc_x)) 
				+ FitResultBkgd->GetBinContent(FitResultBkgd->FindBin(loc_x));
			double loc_stat_unc  = pow(h1_pull->GetBinError(ibin), 2.0);
			loc_stat_unc += pow(h1_compton->GetBinError(h1_compton->FindBin(loc_x)), 2.0);
			loc_stat_unc += pow(h1_pair->GetBinError(h1_pair->FindBin(loc_x)), 2.0);
			loc_stat_unc  = sqrt(loc_stat_unc);
			if(loc_stat_unc <= 1.) loc_stat_unc = 1.0;
			h1_pull->SetBinContent(ibin, (loc_count-loc_fit)/loc_stat_unc);
			h1_pull->SetBinError(ibin, 1.0);
			if(loc_x > min_fit_x && loc_x < max_fit_x) {
				h1_pull_dist->Fill((loc_count-loc_fit)/loc_stat_unc);
			}
		}
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
		
		TF1 *f_pull = new TF1("f_pull", "gaus", -4., 4.);
		h1_pull_dist->Fit(f_pull, "R0QL");
		
		lat_new.DrawLatexNDC(0.115, 0.375, Form("#scale[1.5]{Mean of Pull Distribution: %.2f}", f_pull->GetParameter(1)));
		lat_new.DrawLatexNDC(0.115, 0.275, Form("#scale[1.5]{Std Dev. of Pull Distribution: %.2f}", f_pull->GetParameter(2)));
		
		//canvas_fit->Update();
		canvas_draw->Update();
		
		canvas_pull_dist->cd();
		h1_pull_dist->Draw("hist");
		f_pull->Draw("same");
		//canvas_pull_dist->Update();
	}
	
	//h_vertex->Delete();
	//h2_compton->Delete();
	//h2_pair->Delete();
	//h_pair_gen_flux->Delete();
	
	return fit_val;
}
