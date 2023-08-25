
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
		double loc_eb        = 0.;
		double loc_empty_frac = 0., loc_empty_frac_err = 0.;
		double loc_empty_frac_exp = 0.;
		if(ib<n_bins1) {
			int counter         = tagh_counter_vec[ib];
			loc_eb              = tagh_en[counter-1];
			loc_empty_frac      = tagh_empty_frac[counter-1];
			loc_empty_frac_err  = tagh_empty_fracE[counter-1];
			loc_empty_frac_exp  = tagh_empty_frac_exp[counter-1];
		} else {
			int counter         = tagm_counter_vec[ib-n_bins1];
			loc_eb              = tagm_en[counter-1];
			loc_empty_frac      = tagm_empty_frac[counter-1];
			loc_empty_frac_err  = tagm_empty_fracE[counter-1];
			loc_empty_frac_exp  = tagm_empty_frac_exp[counter-1];
		}
		beam_energy[ib]        = loc_eb;
		zeros[ib]              = 0.;
		empty_fraction[ib]     = loc_empty_frac;
		empty_fraction_err[ib] = loc_empty_frac_err;
		empty_fraction_exp[ib] = loc_empty_frac_exp;
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
	
	gEmptyFrac->GetXaxis()->SetRangeUser(5.8, 11.6);
	
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
