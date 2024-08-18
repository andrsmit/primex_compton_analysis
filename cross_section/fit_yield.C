
#include "get_phaseI_energy_bin.C"

int get_neighbor(int tag_sys, int counter, int sim_case);

int fit_yield(int tag_sys, int counter, TH1F *h1, double &yield, double &yieldE, double &chi2) {
	
	int fit_val = 1;
	
	char tag_sys_char[256];
	if(tag_sys==0) sprintf(tag_sys_char, "tagh");
	else           sprintf(tag_sys_char, "tagm");
	
	double min_fit_x = h1->GetBinCenter(h1->FindFirstBinAbove()+2);
	double max_fit_x = h1->GetBinCenter(h1->FindLastBinAbove()-2);
	
	//if(min_fit_x < -3.0) min_fit_x = -3.0;
	
	double eb, loc_flux;
	if(tag_sys==0) { 
		eb        = tagh_en[counter-1];
		loc_flux  = tagh_flux[counter-1];
	} else {
		eb        = tagm_en[counter-1];
		loc_flux  = tagm_flux[counter-1];
	}
	double loc_acc = f_acc->Eval(eb);
	
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
	
	TH1F *hbeam    = (TH1F*)fSim->Get(Form("%s/beam", comp_root_dir_name.Data()))->Clone( 
		Form("h_beam_%d_%d", tag_sys, counter));
	TH1F *h_vertex = (TH1F*)fSim->Get(Form("%s/vertex_accepted", comp_root_dir_name.Data()))->Clone(
		Form("h_vertex_accepted_fit_%d_%d", tag_sys, counter));
	
	double n_comp_thrown = h_vertex->Integral() * (double)hbeam->GetMean();
	
	double comp_cs       = f_theory->Eval(eb);
	double comp_flux_sim = n_comp_thrown / (comp_cs * n_e * mb);
	
	if(n_comp_thrown < 1.e3) return 0;
	TH2F *h2;
	
	if(tag_sys==0) 
		h2 = (TH2F*)fSim->Get(Form("%s/%s", comp_root_dir_name.Data(), 
			hname_tagh_comp.Data()))->Clone(Form("h2_tagh_%d", counter));
	else 
		h2 = (TH2F*)fSim->Get(Form("%s/%s", comp_root_dir_name.Data(), 
			hname_tagm_comp.Data()))->Clone(Form("h2_tagh_%d", counter));
	
	h2->SetDirectory(0);
	
	fSim->Close();
	
	TH1F *h1_comp = (TH1F*)h2->ProjectionY(Form("h1_comp_%d_%d", tag_sys, counter));
	
	if(h1_comp->Integral() < 1.e1) return 0;
	
	//=================================================================================//
	// Get the Pair MC Histogram
	
	//cout << "getting pair mc" << endl;
	
	int num_pair_files_used = 0;
	
	TH1F *hbeam_pair             = new TH1F(Form("h_beam_pair_%d_%d", 
		tag_sys, counter), "", 2, -0.5, 1.5);
	TH1F *h_vertex_accepted_pair = new TH1F(Form("h_vertex_accepted_pair_fit_%d_%d", 
		tag_sys, counter), "", 1000, 0., 100.);
	TH1F *h_weight_pair          = new TH1F(Form("h_weight_pair_%d_%d", 
		tag_sys, counter), "", 1000, 0., 1.e7);
	
	TH2F *h2p;
	if(old_tag_sys==0) {
		h2p = new TH2F(Form("h2p_%d_%d", tag_sys, counter), 
			"#DeltaK; TAGH Counter; [GeV]", 274, 0.5, 274.5, 2000, -8.0, 8.0);
	} else {
		h2p = new TH2F(Form("h2p_%d_%d", tag_sys, counter), 
			"#DeltaK; TAGM Counter; [GeV]", 102, 0.5, 102.5, 2000, -8.0, 8.0);
	}
	
	hbeam_pair->SetDirectory(0);
	h_vertex_accepted_pair->SetDirectory(0);
	h_weight_pair->SetDirectory(0);
	h2p->SetDirectory(0);
	
	for(int ii=0; ii<10; ii++) {
		
		if(old_tag_sys==0) sprintf(fname, "%s/tagh_%03d_%02d.root", pair_mc_dir.Data(), 
			old_tag_counter, ii);
		else               sprintf(fname, "%s/tagm_%03d_%02d.root", pair_mc_dir.Data(), 
			old_tag_counter, ii);
		if(gSystem->AccessPathName(fname)) continue;
		
		num_pair_files_used++;
		
		TFile *fPair = new TFile(fname, "READ");
		
		TH1F *h_pair_vertex = (TH1F*)fPair->Get(Form("%s/vertex", pair_root_dir_name.Data()))->Clone(
			Form("h_pair_vertex_thrown_%d_%d_%d", tag_sys, counter, ii));
		if(h_pair_vertex->Integral() != 1.e6) continue;
		
		hbeam_pair->Add((TH1F*)fPair->Get(Form("%s/beam", pair_root_dir_name.Data())));
		h_vertex_accepted_pair->Add((TH1F*)fPair->Get(Form("%s/vertex_accepted", 
			pair_root_dir_name.Data())));
		h_weight_pair->Add((TH1F*)fPair->Get(Form("%s/pair_event_weight", pair_root_dir_name.Data())));
		
		if(old_tag_sys==0) h2p->Add((TH2F*)fPair->Get(Form("%s/%s", pair_root_dir_name.Data(), 
			hname_tagh_pair.Data())));
		else               h2p->Add((TH2F*)fPair->Get(Form("%s/%s", pair_root_dir_name.Data(), 
			hname_tagm_pair.Data())));
		
		fPair->Close();
	}
	
	if(num_pair_files_used < 1) return 0;
	
	double n_pair_thrown = h_vertex_accepted_pair->Integral() * (double)hbeam_pair->GetMean();
	//double pair_cs     = (double)h_weight_pair->GetMean() / 1.e3;
	double pair_cs       = f_pair_cs->Eval(eb);
	double pair_flux_sim = n_pair_thrown / ((n_e/n_Z) * mb * pair_cs);
	
	TH1F *h1_pair = (TH1F*)h2p->ProjectionY(Form("h1_pair_%d_%d", tag_sys, counter), 
		old_tag_counter, old_tag_counter);
	h1_pair->SetDirectory(0);
	h1_pair->Scale(1./(double)h_weight_pair->GetMean());
	
	//=================================================================================//
	// Get the Triplet MC Histogram
	
	//cout << "getting triplet mc" << endl;
	
	int num_trip_files_used = 0;
	
	TH1F *hbeam_trip             = new TH1F(Form("h_beam_trip_%d_%d", 
		tag_sys, counter), "", 2, -0.5, 1.5);
	TH1F *h_vertex_accepted_trip = new TH1F(Form("h_vertex_accepted_trip_fit_%d_%d", 
		tag_sys, counter), "", 1000, 0., 100.);
	TH1F *h_weight_trip          = new TH1F(Form("h_weight_trip_%d_%d", 
		tag_sys, counter), "", 1000, 0., 1.e7);
	
	TH2F *h2t;
	if(old_tag_sys==0) {
		h2t = new TH2F(Form("h2t_%d_%d", tag_sys, counter), 
			"#DeltaK; TAGH Counter; [GeV]", 274, 0.5, 274.5, 2000, -8.0, 8.0);
	} else {
		h2t = new TH2F(Form("h2t_%d_%d", tag_sys, counter), 
			"#DeltaK; TAGM Counter; [GeV]", 102, 0.5, 102.5, 2000, -8.0, 8.0);
	}
	
	hbeam_trip->SetDirectory(0);
	h_vertex_accepted_trip->SetDirectory(0);
	h_weight_trip->SetDirectory(0);
	h2t->SetDirectory(0);
	
	for(int ii=0; ii<10; ii++) {
		
		if(old_tag_sys==0) sprintf(fname, "%s/tagh_%03d_%02d.root", trip_mc_dir.Data(), 
			old_tag_counter, ii);
		else               sprintf(fname, "%s/tagm_%03d_%02d.root", trip_mc_dir.Data(), 
			old_tag_counter, ii);
		if(gSystem->AccessPathName(fname)) continue;
		
		num_trip_files_used++;
		
		TFile *fTrip = new TFile(fname, "READ");
		
		hbeam_trip->Add((TH1F*)fTrip->Get(Form("%s/beam", trip_root_dir_name.Data())));
		h_vertex_accepted_trip->Add((TH1F*)fTrip->Get(Form("%s/vertex_accepted", 
			trip_root_dir_name.Data())));
		h_weight_trip->Add((TH1F*)fTrip->Get(Form("%s/pair_event_weight", trip_root_dir_name.Data())));
		
		if(old_tag_sys==0) h2t->Add((TH2F*)fTrip->Get(Form("%s/%s", trip_root_dir_name.Data(), 
			hname_tagh_pair.Data())));
		else               h2t->Add((TH2F*)fTrip->Get(Form("%s/%s", trip_root_dir_name.Data(), 
			hname_tagm_pair.Data())));
		
		fTrip->Close();
	}
	
	if(num_trip_files_used < 1) return 0;
	
	double n_trip_thrown = h_vertex_accepted_trip->Integral() * (double)hbeam_trip->GetMean();
	//double trip_cs     = (double)h_weight_trip->GetMean() / 1.e3;
	double trip_cs       = f_trip_cs->Eval(eb);
	double trip_flux_sim = n_trip_thrown / ((n_e/n_Z) * mb * trip_cs);
	
	TH1F *h1_trip = (TH1F*)h2t->ProjectionY(Form("h1_trip_%d_%d", tag_sys, counter), 
		old_tag_counter, old_tag_counter);
	h1_trip->SetDirectory(0);
	h1_trip->Scale(1./(double)h_weight_trip->GetMean());
	
	//=================================================================================//
	
	double pair_cs_ratio = 1.e3*pair_cs / (double)(h_weight_pair->GetMean());
	double trip_cs_ratio = 1.e3*trip_cs / (double)(h_weight_trip->GetMean());
	
	h1_comp->Rebin(rebins);
	h1_comp->Scale(loc_flux / comp_flux_sim);
	h1_comp->GetXaxis()->SetTitle("E_{Comp} - E_{#gamma} [GeV]");
	h1_comp->GetXaxis()->SetTitleOffset(1.1);
	h1_comp->GetYaxis()->SetTitle( Form("counts / %d MeV", n_mev));
	h1_comp->GetYaxis()->SetTitleOffset(1.3);
	h1_comp->SetLineColor(kCyan);
	h1_comp->SetFillColor(kCyan);
	h1_comp->SetFillStyle(3004);
	h1_comp->SetLineWidth(2);
	
	h1_pair->Rebin(rebins);
	h1_pair->Scale(loc_flux / pair_flux_sim);
	h1_pair->GetXaxis()->SetTitle("E_{Comp} - E_{#gamma} [GeV]");
	h1_pair->GetXaxis()->SetTitleOffset(1.1);
	h1_pair->GetYaxis()->SetTitle(Form("counts / %d MeV", n_mev));
	h1_pair->GetYaxis()->SetTitleOffset(1.3);
	h1_pair->SetLineColor(kGreen);
	h1_pair->SetFillColor(kGreen);
	h1_pair->SetFillStyle(3004);
	h1_pair->SetLineWidth(2);
	
	h1_trip->Rebin(rebins);
	h1_trip->Scale(loc_flux / trip_flux_sim);
	h1_trip->GetXaxis()->SetTitle("E_{Comp} - E_{#gamma} [GeV]");
	h1_trip->GetXaxis()->SetTitleOffset(1.1);
	h1_trip->GetYaxis()->SetTitle(Form("counts / %d MeV", n_mev));
	h1_trip->GetYaxis()->SetTitleOffset(1.3);
	h1_trip->SetLineColor(kMagenta);
	h1_trip->SetFillColor(kMagenta);
	h1_trip->SetFillStyle(3004);
	h1_trip->SetLineWidth(2);
	
	for(int ib=1; ib<=h1_pair->GetXaxis()->GetNbins(); ib++) {
		if(h1_comp->GetBinContent(ib)<0.) {
			h1_comp->SetBinContent(ib,0.);
		}
	}
	for(int ib=1; ib<=h1->GetXaxis()->GetNbins(); ib++) {
		if(h1->GetBinContent(ib)<0.) {
			h1->SetBinContent(ib,0.);
		}
	}
	
	h1->GetXaxis()->SetRangeUser(min_fit_x, max_fit_x);
	h1_comp->GetXaxis()->SetRangeUser(min_fit_x, max_fit_x);
	h1_pair->GetXaxis()->SetRangeUser(min_fit_x, max_fit_x);
	h1_trip->GetXaxis()->SetRangeUser(min_fit_x, max_fit_x);
	
	double n_total_exp = h1->Integral(h1->FindBin(min_fit_x), h1->FindBin(max_fit_x));
	double n_comp_exp  = h1_comp->Integral(h1_comp->FindBin(min_fit_x), h1_comp->FindBin(max_fit_x));
	double n_pair_mc   = h1_pair->Integral(h1_pair->FindBin(min_fit_x), h1_pair->FindBin(max_fit_x));
	double n_trip_mc   = h1_trip->Integral(h1_trip->FindBin(min_fit_x), h1_trip->FindBin(max_fit_x));
	
	//n_pair_mc *= (1.e3 * pair_cs / (double)h_weight_pair->GetMean());
	//n_trip_mc *= (1.e3 * trip_cs / (double)h_weight_trip->GetMean());
	
	double comp_frac_exp = n_comp_exp / n_total_exp;
	double pair_frac_exp = n_pair_mc / n_total_exp;
	double trip_frac_exp = n_trip_mc / n_total_exp;
	
	pair_frac_exp += trip_frac_exp;
	h1_pair->Add(h1_trip);
	
	//-----   Perform Binned Maximum Likelihood Fit (TFractionFitter)   -----//
	
	if(DEBUG_FITS) {
		c_debug->cd();
		h1->Draw("PE");
		h1_comp->Draw("same hist");
		h1_pair->Draw("same hist");
		c_debug->Update();
	}
	
	TObjArray *mc = new TObjArray(2);
	mc->Add(h1_comp);
	mc->Add(h1_pair);
	
	TFractionFitter *fit = new TFractionFitter(h1, mc);
	
	fit->Constrain(0, 0.700,               1.000);
	fit->Constrain(1, 0.1*pair_frac_exp,  10.0*pair_frac_exp);
	
	fit->SetRangeX(h1->FindBin(min_fit_x),h1->FindBin(max_fit_x));
	
	Int_t fit_status = fit->Fit();
	TH1F *fit_result = (TH1F*)fit->GetPlot();
	
	double comp_frac, comp_fracE;
	double pair_frac, pair_fracE;
	fit->GetResult(0, comp_frac, comp_fracE);
	fit->GetResult(1, pair_frac, pair_fracE);
	
	//yield  = Ncomp;
	yield  = n_total_exp * (1.0 - pair_frac);
	yieldE = sqrt(yield);
	
	chi2 = fit->GetChisquare() / fit->GetNDF();
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
	
	if(tag_sys==0) {
		tagh_comp_frac[counter-1]     = comp_frac;
		tagh_pair_frac[counter-1]     = pair_frac;
		tagh_comp_fracE[counter-1]    = comp_fracE;
		tagh_pair_fracE[counter-1]    = pair_fracE;
		tagh_pair_frac_exp[counter-1] = pair_frac_exp;
	} else {
		tagm_comp_frac[counter-1]  = comp_frac;
		tagm_pair_frac[counter-1]  = pair_frac;
		tagm_comp_fracE[counter-1] = comp_fracE;
		tagm_pair_fracE[counter-1] = pair_fracE;
		tagm_pair_frac_exp[counter-1] = pair_frac_exp;
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
		
		TH1F *h_comp_fit, *h_pair_fit;
		
		h_comp_fit = (TH1F*)fit->GetMCPrediction(0);
		h_comp_fit->GetXaxis()->SetRangeUser(min_fit_x, max_fit_x);
		h_pair_fit = (TH1F*)fit->GetMCPrediction(1);
		h_pair_fit->GetXaxis()->SetRangeUser(min_fit_x, max_fit_x);
		
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
		
		TH1F *h1_lin = (TH1F*)h1->Clone(Form("h1_lin_%d_%d", tag_sys, counter));
		h1_lin->SetMinimum(0.);
		
		TLegend *leg = new TLegend( 0.140, 0.545, 0.400, 0.855 );
		leg->AddEntry(h1_lin,     "Data",                 "PE");
		leg->AddEntry(fit_result, "Full Fit",              "l");
		leg->AddEntry(h_comp_fit, "Compton Signal (fit)",  "l");
		leg->AddEntry(h_pair_fit, "e+e- Background (fit)", "l");
		leg->SetBorderSize(0);
		
		pad_lin->cd();
		h1_lin->Draw("PE");
		fit_result->Draw("same");
		h_comp_fit->Draw("same hist");
		h_pair_fit->Draw("same hist");
		h1_lin->Draw("PE same");
		leg->Draw();
		
		pad_log->cd();
		h1->Draw("PE");
		fit_result->Draw("same");
		h_comp_fit->Draw("same hist");
		h_pair_fit->Draw("same hist");
		h1->Draw("PE same");
		
		fit_result->SetMaximum(fit_result->GetMaximum()*1.4);
		
		fit_result->GetYaxis()->SetTitle(Form("Counts / %d MeV", (int)(bin_size*1.e3)));
		fit_result->GetYaxis()->SetTitleOffset(1.8);
		fit_result->SetTitle("");
		
		canvas_draw->cd();
		fit_result->Draw("hist");
		//FitResultAll->Draw("hist same");
		
		TH1F *FitResultSignal = (TH1F*)h1_comp->Clone(Form("FitResultSignal_%d_%d", tag_sys, counter));
		FitResultSignal->Scale(1./FitResultSignal->Integral()*comp_frac*Ndata);
		TH1F *FitResultBkgd   = (TH1F*)h1_pair->Clone(Form("FitResultBkgd_%d_%d", tag_sys, counter));
		FitResultBkgd->Scale(1./FitResultBkgd->Integral()*pair_frac*Ndata);
		
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
		leg1->Draw();
		
		TLatex lat_new;
		lat_new.SetTextColor(kRed-6);
		lat_new.SetTextFont(52);
		lat_new.DrawLatexNDC(0.63, 0.83, Form("#scale[1.0]{E_{#gamma} = %.2f GeV}", eb));
		
		canvas1->Update();
		canvas_draw->Update();
	}
	
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


