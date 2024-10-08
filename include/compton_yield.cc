#include "compton_cs.h"
#include "fit_yield.cc"
#include "fit_yield_rebins.cc"

void get_compton_yield(vector<int> &tagh_counter_vec, vector<int> &tagm_counter_vec) {
	
	tagh_counter_vec.clear();
	tagm_counter_vec.clear();
	
	double loc_MIN_BEAM_ENERGY =  6.5;
	double loc_MAX_BEAM_ENERGY = 11.1;
	
	//double loc_MIN_BEAM_ENERGY =  6.5;
	//double loc_MAX_BEAM_ENERGY =  7.8;
	
	for(int tagh_counter=1; tagh_counter<=N_TAGH_COUNTERS; tagh_counter++) {
		tagh_yield[tagh_counter-1]  = 0.;
		tagh_yieldE[tagh_counter-1] = 0.;
		
		tagh_yieldfit_chi2.push_back(0.);
		tagh_yieldfit_pull.push_back({0.,0.});
	}
	for(int tagm_counter=1; tagm_counter<=N_TAGM_COUNTERS; tagm_counter++) {
		tagm_yield[tagm_counter-1]  = 0.;
		tagm_yieldE[tagm_counter-1] = 0.;
		
		tagm_yieldfit_chi2.push_back(0.);
		tagm_yieldfit_pull.push_back({0.,0.});
	}
	
	for(int imc=0; imc<(int)mc_templates.size(); imc++) {
		tagh_fit_fraction[mc_templates[imc]]     = tagh_yieldfit_pull;
		tagh_fit_fraction_exp[mc_templates[imc]] = tagh_yieldfit_pull;
		tagm_fit_fraction[mc_templates[imc]]     = tagm_yieldfit_pull;
		tagm_fit_fraction_exp[mc_templates[imc]] = tagm_yieldfit_pull;
	}
	
	hist_color_map["compton"] = kCyan;
	hist_color_map["pair"]    = kGreen;
	hist_color_map["triplet"] = kMagenta;
	hist_color_map["empty"]   = kBlue;
	
	TFile *fFull  = new TFile(root_fname.Data(),              "READ");
	TFile *fEmpty = new TFile(empty_target_root_fname.Data(), "READ");
	
	TH2F *h2_tagh  = (TH2F*)fFull->Get(hname_tagh.Data())->Clone("h2_tagh");
	TH2F *h2_tagm  = (TH2F*)fFull->Get(hname_tagm.Data())->Clone("h2_tagm");
	TH2F *h2e_tagh = (TH2F*)fEmpty->Get(hname_tagh.Data())->Clone("h2e_tagh");
	TH2F *h2e_tagm = (TH2F*)fEmpty->Get(hname_tagm.Data())->Clone("h2e_tagm");
	
	TH1F *h1_tagh[N_TAGH_COUNTERS], *h1e_tagh[N_TAGH_COUNTERS];
	TH1F *h1_tagm[N_TAGM_COUNTERS], *h1e_tagm[N_TAGM_COUNTERS];
	
	for(int tagh_counter=1; tagh_counter<=274; tagh_counter++) {
		
		double eb             = tagh_en[tagh_counter-1];
		double loc_flux       = tagh_flux[tagh_counter-1];
		double loc_flux_empty = tagh_flux_empty[tagh_counter-1];
		
		if(eb<loc_MIN_BEAM_ENERGY || eb>loc_MAX_BEAM_ENERGY) continue;
		
		int phase1_tag_sys     = 0;
		int phase1_tag_counter = tagh_counter;
		
		if(!IS_BE_TARGET && FIT_USING_BE_EMPTY==true) {
			get_phase1_energy_bin(eb, phase1_tag_counter, phase1_tag_sys);
		}
		
		cout << "\n\n";
		cout << "TAGH counter: " << tagh_counter << endl;
		cout << "Phase 1 counter: " << phase1_tag_counter << " (tag_sys=" << phase1_tag_sys << ")" << endl;
		cout << "\n\n";
		
		if(phase1_tag_sys==0) {
			loc_flux_empty = tagh_flux_empty[phase1_tag_counter-1];
		} else if(phase1_tag_sys==1) {
			loc_flux_empty = tagm_flux_empty[phase1_tag_counter-1];
		} else {
			continue;
		}
		
		if(loc_flux <= 0. || loc_flux_empty <= 0.) continue;
		if(phase1_tag_sys==1 && phase1_tag_counter==21) continue;
		
		h1_tagh[tagh_counter-1]  = (TH1F*)h2_tagh->ProjectionY(Form("h1_tagh_%03d",
			tagh_counter),tagh_counter,tagh_counter);
		h1_tagh[tagh_counter-1]->SetDirectory(0);
		
		if(phase1_tag_sys==0) {
			h1e_tagh[tagh_counter-1] = (TH1F*)h2e_tagh->ProjectionY(Form("h1e_tagh_%03d", 
				tagh_counter), phase1_tag_counter, phase1_tag_counter);
		} else {
			h1e_tagh[tagh_counter-1] = (TH1F*)h2e_tagm->ProjectionY(Form("h1e_tagh_%03d", 
				tagh_counter), phase1_tag_counter, phase1_tag_counter);
		}
		h1e_tagh[tagh_counter-1]->Scale(loc_flux/loc_flux_empty);
		h1e_tagh[tagh_counter-1]->SetDirectory(0);
	}
	
	for(int tagm_counter=1; tagm_counter<=102; tagm_counter++) {
		
		double eb             = tagm_en[tagm_counter-1];
		double loc_flux       = tagm_flux[tagm_counter-1];
		double loc_flux_empty = tagm_flux_empty[tagm_counter-1];
		
		int phase1_tag_sys     = 1;
		int phase1_tag_counter = tagm_counter;
		
		if(!IS_BE_TARGET && FIT_USING_BE_EMPTY==true) {
			get_phase1_energy_bin(eb, phase1_tag_counter, phase1_tag_sys);
		}
		
		if(phase1_tag_sys==0) {
			loc_flux_empty = tagh_flux_empty[phase1_tag_counter-1];
		} else if(phase1_tag_sys==1) {
			loc_flux_empty = tagm_flux_empty[phase1_tag_counter-1];
		} else {
			continue;
		}
		
		if(loc_flux <= 0. || loc_flux_empty <= 0.) continue;
		
		h1_tagm[tagm_counter-1]  = (TH1F*)h2_tagm->ProjectionY(Form("h1_tagm_%03d",
			tagm_counter),tagm_counter,tagm_counter);
		h1_tagm[tagm_counter-1]->SetDirectory(0);
		
		if(phase1_tag_sys==0) {
			h1e_tagm[tagm_counter-1] = (TH1F*)h2e_tagh->ProjectionY(Form("h1e_tagm_%03d", 
				tagm_counter), phase1_tag_counter, phase1_tag_counter);
		} else {
			h1e_tagm[tagm_counter-1] = (TH1F*)h2e_tagm->ProjectionY(Form("h1e_tagm_%03d", 
				tagm_counter), phase1_tag_counter, phase1_tag_counter);
		}
		h1e_tagm[tagm_counter-1]->Scale(loc_flux/loc_flux_empty);
		h1e_tagm[tagm_counter-1]->SetDirectory(0);
	}
	
	fFull->Close();
	fEmpty->Close();
	
	if(DRAW_FITS_TAGH || DRAW_FITS_TAGM) {
		/*
		canvas_fit = new TCanvas("canvas_fit", "Fits", 1200, 600);
		canvas_fit->Divide(2,1);
		
		canvas_fit_lin = (TPad*)canvas_fit->cd(1);
		canvas_fit_lin->SetTickx(); canvas_fit_lin->SetTicky();
		
		canvas_fit_log = (TPad*)canvas_fit->cd(2);
		canvas_fit_log->SetTickx(); canvas_fit_log->SetTicky();
		canvas_fit_log->SetGrid();  canvas_fit_log->SetLogy();
		*/
		canvas_draw = new TCanvas("canvas_draw","canvas_draw",800,800);
		canvas_draw->SetTickx(); canvas_draw->SetTicky();
		canvas_draw->SetTopMargin(0.10);
		canvas_draw->SetBottomMargin(0.10);
		canvas_draw->SetRightMargin(0.05);
		canvas_draw->SetLeftMargin(0.15);
		
		canvas_draw_top = new TPad("canvas_draw_top", "top_pad", 0.005, 0.30, 0.995, 0.995);
		canvas_draw_bot = new TPad("canvas_draw_bot", "bot_pad", 0.005, 0.005,  0.995, 0.2975);
		
		canvas_draw_top->SetLeftMargin(0.10);
		canvas_draw_top->SetRightMargin(0.02);
		canvas_draw_top->SetTopMargin(0.075);
		canvas_draw_top->SetBottomMargin(0.015);
		canvas_draw_top->SetTickx(); canvas_draw_top->SetTicky();
		canvas_draw_top->SetFrameLineWidth(2);
		
		canvas_draw_bot->SetLeftMargin(0.10);
		canvas_draw_bot->SetRightMargin(0.02);
		canvas_draw_bot->SetTopMargin(0.005);
		canvas_draw_bot->SetBottomMargin(0.225);
		canvas_draw_bot->SetTickx(); canvas_draw_bot->SetTicky();
		canvas_draw_bot->SetFrameLineWidth(2);
		
		canvas_draw->cd();
		canvas_draw_bot->Draw();
		canvas_draw_top->Draw();
		
		canvas_pull = new TCanvas("canvas_pull", "Pull of Fit", 1200, 400);
		canvas_pull->SetTickx(); canvas_pull->SetTicky();
		
		canvas_pull_dist = new TCanvas("canvas_pull_dist", "Pull of Fit", 700, 500);
		canvas_pull_dist->SetTickx(); canvas_pull_dist->SetTicky();
	}
	
	//------------------------------------------------------------------------------------------//
	
	for(int tagh_counter = 1; tagh_counter <= 221; tagh_counter++) {
		
		int bad_val = 0;
		for(int ic = 0; ic < bad_counters_tagh.size(); ic++) {
			if(tagh_counter==bad_counters_tagh[ic]) bad_val = 1;
		}
		if(bad_val) continue;
		
		double eb              = tagh_en[tagh_counter-1];
		double loc_flux        = tagh_flux[tagh_counter-1];
		double loc_flux_empty  = tagh_flux_empty[tagh_counter-1];
		
		if(eb<loc_MIN_BEAM_ENERGY || eb>loc_MAX_BEAM_ENERGY) continue;
		
		int phase1_tag_sys     = 0;
		int phase1_tag_counter = tagh_counter;
		
		if(!IS_BE_TARGET) {
			get_phase1_energy_bin(eb, phase1_tag_counter, phase1_tag_sys);
		}
		
		if(phase1_tag_sys==0) {
			loc_flux_empty = tagh_flux_empty[phase1_tag_counter-1];
		} else if(phase1_tag_sys==1) {
			loc_flux_empty = tagm_flux_empty[phase1_tag_counter-1];
		} else {
			continue;
		}
		
		// Skip bins that have no flux:
		
		if(loc_flux <= 0. || loc_flux_empty <= 0.) {
			//cout << "Skipping TAGH counter " << tagh_counter << " (no flux)" << endl;
			continue;
		}
		
		// Skip bins that have no yield:
		
		if(h1_tagh[tagh_counter-1]->Integral() < 1.e1) {
			//cout << "Skipping TAGH counter " << tagh_counter << " (no yield)" << endl;
			continue;
		}
		
		h1_tagh[tagh_counter-1]->SetTitle( 
			Form("TAGH Counter %d (E_{#gamma} = %.3f GeV)", tagh_counter, eb));
		h1_tagh[tagh_counter-1]->Rebin(rebins);
		h1_tagh[tagh_counter-1]->SetLineColor(kBlack);
		h1_tagh[tagh_counter-1]->GetXaxis()->SetTitle("E_{Comp}#left(#theta_{1},#theta_{2}#right) - E_{#gamma} [GeV]");
		h1_tagh[tagh_counter-1]->SetLineWidth(2);
		
		h1e_tagh[tagh_counter-1]->Rebin(rebins);
		
		cout << endl;
		cout << "==================================" << endl;
		cout << "Fitting TAGH Counter " << tagh_counter << endl;
		
		double loc_yield = 0., loc_yieldE = 0., loc_chi2 = 0.;
		int fit_val = fit_yield(0, tagh_counter, h1_tagh[tagh_counter-1], h1e_tagh[tagh_counter-1], 
			loc_yield, loc_yieldE, loc_chi2);
		
		cout << "Done." << endl;
		cout << "==================================" << endl;
		
		if(fit_val <= 0) continue;
		
		tagh_yield[tagh_counter-1]         = loc_yield;
		tagh_yieldE[tagh_counter-1]        = loc_yieldE;
		tagh_yieldfit_chi2[tagh_counter-1] = loc_chi2;
		
		if(SAVE_FITS_TAGH && DRAW_FITS_TAGH) {
			canvas_draw->SaveAs(Form("yield_fit_tagh_%03d.pdf", tagh_counter), "pdf");
		}
		
		tagh_counter_vec.push_back(tagh_counter);
	}
	
	for(int tagm_counter = 5; tagm_counter <= 96; tagm_counter++) {
		
		int bad_val = 0;
		for(int ic = 0; ic < bad_counters_tagm.size(); ic++) {
			if(tagm_counter==bad_counters_tagm[ic]) bad_val = 1;
		}
		if(bad_val) continue;
		
		double eb              = tagm_en[tagm_counter-1];
		double loc_flux        = tagm_flux[tagm_counter-1];
		double loc_flux_empty  = tagm_flux_empty[tagm_counter-1];
		
		int phase1_tag_sys     = 1;
		int phase1_tag_counter = tagm_counter;
		
		if(!IS_BE_TARGET) {
			get_phase1_energy_bin(eb, phase1_tag_counter, phase1_tag_sys);
		}
		
		if(phase1_tag_sys==0) {
			loc_flux_empty = tagh_flux_empty[phase1_tag_counter-1];
		} else if(phase1_tag_sys==1) {
			loc_flux_empty = tagm_flux_empty[phase1_tag_counter-1];
		} else {
			continue;
		}
		
		// Skip bins that have no flux:
		
		if(loc_flux <= 0. || loc_flux_empty <= 0.) {
			//cout << "Skipping TAGM counter " << tagm_counter << " (no flux)" << endl;
			continue;
		}
		
		// Skip bins that have no yield:
		
		if(h1_tagm[tagm_counter-1]->Integral() < 1.e1) {
			//cout << "Skipping TAGM counter " << tagm_counter << " (no yield)" << endl;
			continue;
		}
		
		h1_tagm[tagm_counter-1]->SetTitle( 
			Form("TAGM Counter %d (E_{#gamma} = %.3f GeV)", tagm_counter, eb));
		h1_tagm[tagm_counter-1]->Rebin(rebins);
		h1_tagm[tagm_counter-1]->SetLineColor(kBlack);
		h1_tagm[tagm_counter-1]->GetXaxis()->SetTitle("E_{Comp}#left(#theta_{1},#theta_{2}#right) - E_{#gamma} [GeV]");
		h1_tagm[tagm_counter-1]->SetLineWidth(2);
		
		h1e_tagm[tagm_counter-1]->Rebin(rebins);
		
		cout << endl;
		cout << "==================================" << endl;
		cout << "Fitting TAGM Counter " << tagm_counter << endl;
		
		double loc_yield = 0., loc_yieldE = 0., loc_chi2 = 0.;
		int fit_val = fit_yield(1, tagm_counter, h1_tagm[tagm_counter-1], h1e_tagm[tagm_counter-1], 
			loc_yield, loc_yieldE, loc_chi2);
		
		if(fit_val <= 0) continue;
		
		tagm_yield[tagm_counter-1]         = loc_yield;
		tagm_yieldE[tagm_counter-1]        = loc_yieldE;
		tagm_yieldfit_chi2[tagm_counter-1] = loc_chi2;
		
		if(SAVE_FITS_TAGM && DRAW_FITS_TAGM) {
			canvas_draw->SaveAs(Form("yield_fit_tagm_%03d.pdf", tagm_counter), "pdf");
		}
		
		tagm_counter_vec.push_back(tagm_counter);
	}
	
	return;
}
