#include "/work/halld/home/andrsmit/primex_compton_analysis/include/compton_inputs.h"

TGraphErrors *gCS;
TCanvas *cCS;

void calc_cs(int tag_sys, int counter, double &loc_cs, double &loc_csE) {
	
	loc_cs = 0., loc_csE = 0.;
	
	double loc_eb;
	double loc_yield,  loc_flux,  loc_acc,  loc_fabs;
	double loc_yieldE, loc_fluxE, loc_accE, loc_fabsE;
	
	if(tag_sys==0) {
		
		loc_eb     = tagh_en[counter-1];
		loc_yield  = tagh_yield[counter-1];
		loc_yieldE = tagh_yieldE[counter-1];
		loc_flux   = tagh_flux[counter-1];
		loc_fluxE  = tagh_fluxE[counter-1];
		if(USE_F_ACC) loc_acc = f_acc->Eval(loc_eb);
		else {
			double test_acc = tagh_acc[counter-1];
			if(test_acc==0.) loc_acc = f_acc->Eval(tagh_en[counter-1]);
			else loc_acc = test_acc;
		}
	} else {
		
		loc_eb     = tagm_en[counter-1];
		loc_yield  = tagm_yield[counter-1];
		loc_yieldE = tagm_yieldE[counter-1];
		loc_flux   = tagm_flux[counter-1];
		loc_fluxE  = tagm_fluxE[counter-1];
		if(USE_F_ACC) loc_acc = f_acc->Eval(loc_eb);
		else {
			double test_acc = tagm_acc[counter-1];
			if(test_acc==0.) loc_acc = f_acc->Eval(tagm_en[counter-1]);
			else loc_acc = test_acc;
		}
	}
	
	if(loc_flux<=0. || loc_acc<=0. || n_e<=0. || mb<=0. || loc_yield<=0.) {
		loc_cs  = 0.;
		loc_csE = 0.;
		return;
	}
	
	loc_accE = 0.;
	
	loc_cs  = loc_yield / (loc_flux * loc_acc * n_e * mb);
	loc_csE = sqrt(
		pow(loc_yieldE / (loc_flux * loc_acc * n_e * mb), 2.0) + 
		pow(loc_fluxE * loc_yield / (loc_flux * loc_flux * loc_acc * n_e * mb), 2.0)
	);
	
	loc_cs  /= f_abs;
	loc_csE /= f_abs;
	
	return;
}

void plot_cs(vector<int> tagh_counter_vec, vector<int> tagm_counter_vec) {
	
	int n_bins1 = (int)tagh_counter_vec.size();
	int n_bins2 = (int)tagm_counter_vec.size();
	int n_bins  = n_bins1+n_bins2;
	
	double *beam_energy       = new double[n_bins];
	double *zeros             = new double[n_bins];
	double *cross_section     = new double[n_bins];
	double *cross_section_err = new double[n_bins];
	double *cs_deviation      = new double[n_bins];
	double *cs_deviation_err  = new double[n_bins];
	
	for(int tagh_counter=1; tagh_counter<=274; tagh_counter++) {
		tagh_cs[tagh_counter-1]  = 0.;
		tagh_csE[tagh_counter-1] = 0.;
	}
	for(int tagm_counter=1; tagm_counter<=102; tagm_counter++) {
		tagm_cs[tagm_counter-1]  = 0.;
		tagm_csE[tagm_counter-1] = 0.;
	}
	
	for(int ib=0; ib<n_bins; ib++) {
		double loc_eb = 0.;
		double loc_cs = 0., loc_cs_err = 0.;
		if(ib<n_bins1) {
			int counter  = tagh_counter_vec[ib];
			loc_eb = tagh_en[counter-1];
			calc_cs(0, counter, loc_cs, loc_cs_err);
			tagh_cs[counter-1]  = loc_cs;
			tagh_csE[counter-1] = loc_cs_err;
		} else {
			int counter  = tagm_counter_vec[ib-n_bins1];
			loc_eb = tagm_en[counter-1];
			calc_cs(1, counter, loc_cs, loc_cs_err);
			tagm_cs[counter-1]  = loc_cs;
			tagm_csE[counter-1] = loc_cs_err;
		}
		beam_energy[ib]       = loc_eb;
		zeros[ib]             = 0.;
		cross_section[ib]     = loc_cs;
		cross_section_err[ib] = loc_cs_err;
		
		double loc_theory_cs = USE_NIST_CALCULATION ? f_nist->Eval(loc_eb) : f_theory->Eval(loc_eb);
		cs_deviation[ib]      = 100. * (loc_cs - loc_theory_cs) / loc_theory_cs;
		cs_deviation_err[ib]  = 100. * loc_cs_err / loc_theory_cs;
	}
	
	gCS = new TGraphErrors(n_bins, beam_energy, cross_section, zeros, 
		cross_section_err);
	gCS->SetMarkerStyle(8);
	gCS->SetMarkerSize(0.7);
	gCS->SetMarkerColor(kBlue);
	gCS->SetLineColor(kBlue);
	gCS->SetTitle("");
	
	gCS->GetXaxis()->SetTitle("Photon Beam Energy [GeV]");
	gCS->GetXaxis()->SetTitleSize(0.05);
	gCS->GetXaxis()->SetTitleOffset(0.9);
	gCS->GetXaxis()->CenterTitle(true);
	
	gCS->GetYaxis()->SetTitle("#sigma [mb / electron]");
	gCS->GetYaxis()->SetTitleSize(0.05);
	gCS->GetYaxis()->SetTitleOffset(0.9);
	gCS->GetYaxis()->CenterTitle(true);
	
	gCS->GetXaxis()->SetRangeUser(5.8, 11.6);
	gCS->GetYaxis()->SetRangeUser(0.11, 0.25);
	
	cCS = new TCanvas("cCS", "cCS", 1200, 500);
	cCS->SetTickx(); cCS->SetTicky();
	gCS->Draw("AP");
	f_theory->Draw("same");
	
	TGraphErrors *gDev = new TGraphErrors(n_bins, beam_energy, cs_deviation, zeros, 
		cs_deviation_err);
	gDev->SetMarkerStyle(8);
	gDev->SetMarkerSize(0.7);
	gDev->SetMarkerColor(kBlue);
	gDev->SetLineColor(kBlue);
	gDev->SetTitle("");
	
	gDev->GetXaxis()->SetTitle("Photon Beam Energy [GeV]");
	gDev->GetXaxis()->SetTitleSize(0.05);
	gDev->GetXaxis()->SetTitleOffset(0.9);
	gDev->GetXaxis()->CenterTitle(true);
	
	gDev->GetYaxis()->SetTitle("100 x (#sigma_{exp} - #sigma_{theory}) / #sigma_{theory} [%]");
	gDev->GetYaxis()->SetTitleSize(0.05);
	gDev->GetYaxis()->SetTitleOffset(0.9);
	gDev->GetYaxis()->CenterTitle(true);
	
	gDev->GetXaxis()->SetRangeUser(5.8, 11.6);
	gDev->GetYaxis()->SetRangeUser(-11., 11.);
	
	TCanvas *cDev = new TCanvas("cDev", "cDev", 1200, 500);
	cDev->SetTickx(); cDev->SetTicky();
	gDev->Draw("AP");
	cDev->Update();
	TLine *l0 = new TLine(gPad->GetUxmin(), 0., gPad->GetUxmax(), 0.);
	l0->SetLineColor(kRed);
	l0->SetLineWidth(2);
	l0->Draw("same");
	
	//---------------------------------------------------------------------------------------//
	
	TCanvas *canvas_cs = new TCanvas("canvas3", "Compton CS", 1200, 900);
	
	TPad *top_pad = new TPad("top_pad", "top_pad", 0.005, 0.3525, 0.995, 0.995);
	TPad *bot_pad = new TPad("bot_pad", "bot_pad", 0.005, 0.005,  0.995, 0.3475);
	
	top_pad->SetLeftMargin(0.10);
	top_pad->SetRightMargin(0.02);
	top_pad->SetTopMargin(0.075);
	top_pad->SetBottomMargin(0.015);
	top_pad->SetTickx(); top_pad->SetTicky();
	top_pad->SetFrameLineWidth(2);
	
	bot_pad->SetLeftMargin(0.10);
	bot_pad->SetRightMargin(0.02);
	bot_pad->SetTopMargin(0.010);
	bot_pad->SetBottomMargin(0.225);
	bot_pad->SetTickx(); bot_pad->SetTicky();
	bot_pad->SetFrameLineWidth(2);
	
	canvas_cs->cd();
	top_pad->Draw();
	bot_pad->Draw();
	
	gDev->GetXaxis()->SetTitleSize(0.1);
	gDev->GetXaxis()->SetTitleOffset(0.9);
	gDev->GetYaxis()->SetTitleSize(0.085);
	gDev->GetYaxis()->SetTitleOffset(0.35);
	gDev->GetYaxis()->SetTitle("(#sigma_{exp}-#sigma_{theory})/#sigma_{theory} [%]");
	gDev->SetTitle(" ");
	
	gCS->GetXaxis()->SetTitleSize(0.065);
	gCS->GetXaxis()->SetLabelOffset(0.015);
	gCS->GetYaxis()->SetTitleOffset(0.7);
	
	TF1 *f_theory_clone = (TF1*)f_theory->Clone("f_theory_clone");
	f_theory_clone->SetLineWidth(4);
	TF1 *f_nist_clone = (TF1*)f_nist->Clone("f_nist_clone");
	f_nist_clone->SetLineWidth(4);
	
	top_pad->cd();
	gCS->Draw("AP");
	f_theory->Draw("same");
	//f_nist->Draw("same");
	
	TLatex lat;
	lat.SetTextColor(kBlue);
	//lat.SetTextSize(44);
	lat.SetTextFont(52);
	lat.DrawLatexNDC(0.138,0.122,"#scale[1.0]{^{9}Be Target}");
	
	TGraphErrors *gCSclone = (TGraphErrors*)gCS->Clone("gCSclone");
	gCSclone->SetMarkerSize(1.0);
	
	TLegend *leg = new TLegend(0.595, 0.645, 0.935, 0.835);
	leg->AddEntry(gCSclone, "Data" , "PE");
	//leg->AddEntry(f_nist_clone, "Theory (NIST)" , "l");
	leg->AddEntry(f_theory_clone, "Theory (PrimEx)" , "l");
	leg->SetBorderSize(0);
	leg->Draw();
	
	bot_pad->cd();
	gDev->Draw("AP");
	
	
	TF1 *f_dev_theory = new TF1("f_dev_theory", 
		"1.e2*([0]-[6] + ([1]-[7])*x + ([2]-[8])*pow(x,2.) + ([3]-[9])*pow(x,3.) + ([4]-[10])*pow(x,4.) + ([5]-[11])*pow(x,5.)) / ([0] + [1]*x + [2]*pow(x,2.) + [3]*pow(x,3.) + [4]*pow(x,4.) + [5]*pow(x,5.))", 5.8, 11.6);
	TF1 *f_dev_nist   = new TF1("f_dev_nist", 
		"1.e2*([0]-[6] + ([1]-[7])*x + ([2]-[8])*pow(x,2.) + ([3]-[9])*pow(x,3.) + ([4]-[10])*pow(x,4.) + ([5]-[11])*pow(x,5.)) / ([0] + [1]*x + [2]*pow(x,2.) + [3]*pow(x,3.) + [4]*pow(x,4.) + [5]*pow(x,5.))", 5.8, 11.6);
	for(int ipar=0; ipar<6; ipar++) {
		f_dev_theory->SetParameter(ipar, f_theory->GetParameter(ipar));
		f_dev_nist->SetParameter(ipar, f_nist->GetParameter(ipar));
	}
	for(int ipar=6; ipar<12; ipar++) {
		if(USE_NIST_CALCULATION) f_dev_theory->SetParameter(ipar, f_nist->GetParameter(ipar-6));
		else f_dev_theory->SetParameter(ipar, f_theory->GetParameter(ipar-6));
		f_dev_nist->SetParameter(ipar, f_nist->GetParameter(ipar-6));
	}
	f_dev_theory->SetRange(5.8, 11.2);
	f_dev_theory->SetLineColor(kRed);
	f_dev_theory->SetLineStyle(9);
	f_dev_theory->SetLineWidth(2);
	f_dev_nist->SetLineColor(kBlack);
	f_dev_nist->SetLineStyle(1);
	f_dev_nist->SetLineWidth(2);
	
	f_dev_theory->Draw("same");
	f_dev_nist->Draw("same");
	
	
	/*
	TF1 *f_dev = new TF1("f_dev", "pol0", 8.0, 11.6);
	gDev->Fit("f_dev", "R0");
	
	f_dev->SetRange(5.8,11.6);
	
	f_dev->SetLineColor(kGreen+1);
	f_dev->SetLineStyle(2);
	f_dev->SetLineWidth(2);
	f_dev->Draw("same");
	
	double f_dev_chi2 = f_dev->GetChisquare() / f_dev->GetNDF();
	
	TLatex lat2;
	lat2.SetTextColor(kGreen+1);
	lat2.DrawLatexNDC(0.13,0.85,Form("#scale[1.1]{Deviation: p_{0} = %.3f%% #pm %.3f%%}", 
		f_dev->GetParameter(0), f_dev->GetParError(0)));
	lat2.DrawLatexNDC(0.13,0.75,Form("#scale[1.1]{#chi^{2}/n.d.f = %.2f}", f_dev_chi2));
	
	bot_pad->Update();
	TLine *lDev = new TLine(gPad->GetUxmin(), 0.0, gPad->GetUxmax(), 0.0);
	lDev->SetLineWidth(2);
	lDev->SetLineColor(kRed);
	//lDev->Draw("same");
	*/
	return;
}

void plot_pair_fraction(vector<int> tagh_counter_vec, vector<int> tagm_counter_vec) {
	
	int n_bins1 = (int)tagh_counter_vec.size();
	int n_bins2 = (int)tagm_counter_vec.size();
	int n_bins  = n_bins1+n_bins2;
	
	double *beam_energy       = new double[n_bins];
	double *zeros             = new double[n_bins];
	double *pair_fraction     = new double[n_bins];
	double *pair_fraction_err = new double[n_bins];
	
	double *pair_fraction_exp = new double[n_bins];
	
	for(int ib=0; ib<n_bins; ib++) {
		double loc_eb        = 0.;
		double loc_pair_frac = 0., loc_pair_frac_err = 0.;
		double loc_pair_frac_exp = 0.;
		if(ib<n_bins1) {
			int counter        = tagh_counter_vec[ib];
			loc_eb             = tagh_en[counter-1];
			loc_pair_frac      = tagh_pair_fraction[counter-1];
			loc_pair_frac_err  = tagh_pair_fractionE[counter-1];
			loc_pair_frac_exp  = tagh_pair_fraction_exp[counter-1];
		} else {
			int counter        = tagm_counter_vec[ib-n_bins1];
			loc_eb             = tagm_en[counter-1];
			loc_pair_frac      = tagm_pair_fraction[counter-1];
			loc_pair_frac_err  = tagm_pair_fractionE[counter-1];
			loc_pair_frac_exp  = tagm_pair_fraction_exp[counter-1];
		}
		beam_energy[ib]       = loc_eb;
		zeros[ib]             = 0.;
		pair_fraction[ib]     = loc_pair_frac;
		pair_fraction_err[ib] = loc_pair_frac_err;
		pair_fraction_exp[ib] = loc_pair_frac_exp;
	}
	
	TGraphErrors *gPairFrac = new TGraphErrors(n_bins, beam_energy, pair_fraction, 
		zeros, pair_fraction_err);
	gPairFrac->SetMarkerStyle(8);
	gPairFrac->SetMarkerSize(0.7);
	gPairFrac->SetMarkerColor(kBlue);
	gPairFrac->SetLineColor(kBlue);
	gPairFrac->SetTitle("");
	
	gPairFrac->GetXaxis()->SetTitle("Photon Beam Energy [GeV]");
	gPairFrac->GetXaxis()->SetTitleSize(0.05);
	gPairFrac->GetXaxis()->SetTitleOffset(0.9);
	gPairFrac->GetXaxis()->CenterTitle(true);
	
	gPairFrac->GetYaxis()->SetTitle("p_{Pair} [a.u.]");
	gPairFrac->GetYaxis()->SetTitleSize(0.05);
	gPairFrac->GetYaxis()->SetTitleOffset(0.9);
	gPairFrac->GetYaxis()->CenterTitle(true);
	
	gPairFrac->GetXaxis()->SetRangeUser(5.6, 11.8);
	
	TGraph *gPairFracExp = new TGraph(n_bins, beam_energy, pair_fraction_exp);
	gPairFracExp->SetMarkerStyle(8);
	gPairFracExp->SetMarkerSize(0.7);
	gPairFracExp->SetMarkerColor(kGreen);
	gPairFracExp->SetLineColor(kGreen);
	gPairFracExp->SetTitle("");
	
	TCanvas *cPairFrac = new TCanvas("cPairFrac", "Pair Fraction (fit)", 1200, 500);
	cPairFrac->SetTickx(); cPairFrac->SetTicky();
	gPairFrac->Draw("AP");
	gPairFracExp->Draw("P same");
	
	return;
}

void plot_empty_fraction(vector<int> tagh_counter_vec, vector<int> tagm_counter_vec) {
	
	int n_bins1 = (int)tagh_counter_vec.size();
	int n_bins2 = (int)tagm_counter_vec.size();
	int n_bins  = n_bins1+n_bins2;
	
	double *beam_energy        = new double[n_bins];
	double *zeros              = new double[n_bins];
	double *empty_fraction     = new double[n_bins];
	double *empty_fraction_err = new double[n_bins];
	
	double *empty_fraction_exp = new double[n_bins];
	
	for(int ib=0; ib<n_bins; ib++) {
		double loc_eb         = 0.;
		double loc_empty_frac = 0., loc_empty_frac_err = 0.;
		double loc_empty_frac_exp = 0.;
		if(ib<n_bins1) {
			int counter         = tagh_counter_vec[ib];
			loc_eb              = tagh_en[counter-1];
			loc_empty_frac      = tagh_empty_fraction[counter-1];
			loc_empty_frac_err  = tagh_empty_fractionE[counter-1];
			loc_empty_frac_exp  = tagh_empty_fraction_exp[counter-1];
		} else {
			int counter         = tagm_counter_vec[ib-n_bins1];
			loc_eb              = tagm_en[counter-1];
			loc_empty_frac      = tagm_empty_fraction[counter-1];
			loc_empty_frac_err  = tagm_empty_fractionE[counter-1];
			loc_empty_frac_exp  = tagm_empty_fraction_exp[counter-1];
		}
		beam_energy[ib]        = loc_eb;
		zeros[ib]              = 0.;
		empty_fraction[ib]     = loc_empty_frac     / loc_empty_frac_exp;
		empty_fraction_err[ib] = loc_empty_frac_err / loc_empty_frac_exp;
		empty_fraction_exp[ib] = loc_empty_frac_exp / loc_empty_frac_exp;
	}
	
	TGraphErrors *gEmptyFrac = new TGraphErrors(n_bins, beam_energy, empty_fraction, 
		zeros, empty_fraction_err);
	gEmptyFrac->SetMarkerStyle(8);
	gEmptyFrac->SetMarkerSize(0.7);
	gEmptyFrac->SetMarkerColor(kBlue);
	gEmptyFrac->SetLineColor(kBlue);
	gEmptyFrac->SetTitle("");
	
	gEmptyFrac->GetXaxis()->SetTitle("Photon Beam Energy [GeV]");
	gEmptyFrac->GetXaxis()->SetTitleSize(0.05);
	gEmptyFrac->GetXaxis()->SetTitleOffset(0.9);
	gEmptyFrac->GetXaxis()->CenterTitle(true);
	
	gEmptyFrac->GetYaxis()->SetTitle("p_{Empty} [a.u.]");
	gEmptyFrac->GetYaxis()->SetTitleSize(0.05);
	gEmptyFrac->GetYaxis()->SetTitleOffset(0.9);
	gEmptyFrac->GetYaxis()->CenterTitle(true);
	
	gEmptyFrac->GetXaxis()->SetRangeUser(5.6, 11.8);
	
	TGraph *gEmptyFracExp = new TGraph(n_bins, beam_energy, empty_fraction_exp);
	gEmptyFracExp->SetMarkerStyle(8);
	gEmptyFracExp->SetMarkerSize(0.7);
	gEmptyFracExp->SetMarkerColor(kGreen);
	gEmptyFracExp->SetLineColor(kGreen);
	gEmptyFracExp->SetTitle("");
	
	TCanvas *cEmptyFrac = new TCanvas("cEmptyFrac", "Empty Fraction (fit)", 1200, 500);
	cEmptyFrac->SetTickx(); cEmptyFrac->SetTicky();
	gEmptyFrac->Draw("AP");
	gEmptyFracExp->Draw("P same");
	
	return;
}


bool isNaN(double x) {
	return x != x;
}

void plot_chi2(vector<int> tagh_counter_vec, vector<int> tagm_counter_vec) {
	
	int n_bins1 = (int)tagh_counter_vec.size();
	int n_bins2 = (int)tagm_counter_vec.size();
	int n_bins  = n_bins1+n_bins2;
	
	double *beam_energy = new double[n_bins];
	double *chi2        = new double[n_bins];
	double *pull_mean   = new double[n_bins];
	double *pull_stdev  = new double[n_bins];
	
	for(int ib=0; ib<n_bins; ib++) {
		double loc_eb         = 0.;
		double loc_chi2       = 0.;
		double loc_pull_mean  = 0.;
		double loc_pull_stdev = 0.;
		if(ib<n_bins1) {
			int counter    = tagh_counter_vec[ib];
			loc_eb         = tagh_en[counter-1];
			loc_chi2       = tagh_yieldfit_chi2[counter-1];
			loc_pull_mean  = tagh_yieldfit_pull_mean[counter-1];
			loc_pull_stdev = tagh_yieldfit_pull_stdev[counter-1];
		} else {
			int counter    = tagm_counter_vec[ib-n_bins1];
			loc_eb         = tagm_en[counter-1];
			loc_chi2       = tagm_yieldfit_chi2[counter-1];
			loc_pull_mean  = tagm_yieldfit_pull_mean[counter-1];
			loc_pull_stdev = tagm_yieldfit_pull_stdev[counter-1];
		}
		if(isNaN(loc_chi2)) loc_chi2 = 0.;
		beam_energy[ib] = loc_eb;
		chi2[ib]        = loc_chi2;
		pull_mean[ib]   = loc_pull_mean;
		pull_stdev[ib]  = loc_pull_stdev;
	}
	
	TGraph *gChi2 = new TGraph(n_bins, beam_energy, chi2);
	gChi2->SetMarkerStyle(8);
	gChi2->SetMarkerSize(0.7);
	gChi2->SetMarkerColor(kBlue);
	gChi2->SetLineColor(kBlue);
	gChi2->SetTitle("");
	gChi2->GetXaxis()->SetTitle("Photon Beam Energy [GeV]");
	gChi2->GetXaxis()->SetTitleSize(0.05);
	gChi2->GetXaxis()->SetTitleOffset(0.9);
	gChi2->GetXaxis()->CenterTitle(true);
	gChi2->GetYaxis()->SetTitle("#chi^{2} / n.d.f.");
	gChi2->GetYaxis()->SetTitleSize(0.05);
	gChi2->GetYaxis()->SetTitleOffset(0.9);
	gChi2->GetYaxis()->CenterTitle(true);
	gChi2->GetXaxis()->SetRangeUser(5.6, 11.8);
	
	TCanvas *cChi2 = new TCanvas("cChi2", "Chi2 of Fit", 1200, 500);
	cChi2->SetTickx(); cChi2->SetTicky();
	gChi2->Draw("AP");
	
	
	TGraph *gPullMean = new TGraph(n_bins, beam_energy, pull_mean);
	gPullMean->SetMarkerStyle(8);
	gPullMean->SetMarkerSize(0.7);
	gPullMean->SetMarkerColor(kBlue);
	gPullMean->SetLineColor(kBlue);
	gPullMean->SetTitle("");
	gPullMean->GetXaxis()->SetTitle("Photon Beam Energy [GeV]");
	gPullMean->GetXaxis()->SetTitleSize(0.05);
	gPullMean->GetXaxis()->SetTitleOffset(0.9);
	gPullMean->GetXaxis()->CenterTitle(true);
	gPullMean->GetYaxis()->SetTitle("Average Pull of Fit");
	gPullMean->GetYaxis()->SetTitleSize(0.05);
	gPullMean->GetYaxis()->SetTitleOffset(0.9);
	gPullMean->GetYaxis()->CenterTitle(true);
	gPullMean->GetXaxis()->SetRangeUser(5.6, 11.8);
	
	TGraph *gPullStdev = new TGraph(n_bins, beam_energy, pull_stdev);
	gPullStdev->SetMarkerStyle(8);
	gPullStdev->SetMarkerSize(0.7);
	gPullStdev->SetMarkerColor(kBlue);
	gPullStdev->SetLineColor(kBlue);
	gPullStdev->SetTitle("");
	gPullStdev->GetXaxis()->SetTitle("Photon Beam Energy [GeV]");
	gPullStdev->GetXaxis()->SetTitleSize(0.05);
	gPullStdev->GetXaxis()->SetTitleOffset(0.9);
	gPullStdev->GetXaxis()->CenterTitle(true);
	gPullStdev->GetYaxis()->SetTitle("Std Dev of Pull");
	gPullStdev->GetYaxis()->SetTitleSize(0.05);
	gPullStdev->GetYaxis()->SetTitleOffset(0.9);
	gPullStdev->GetYaxis()->CenterTitle(true);
	gPullStdev->GetXaxis()->SetRangeUser(5.6, 11.8);
	
	TCanvas *cPull = new TCanvas("cPull", "Pull of Fit", 1200, 1000);
	cPull->Divide(1,2);
	TPad *pPullMean = (TPad*)cPull->cd(1);
	pPullMean->SetTickx(); pPullMean->SetTicky();
	gPullMean->Draw("AP");
	TPad *pPullStdev = (TPad*)cPull->cd(2);
	pPullStdev->SetTickx(); pPullStdev->SetTicky();
	gPullStdev->Draw("AP");
	
	return;
}

