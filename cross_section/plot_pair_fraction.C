
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
			loc_pair_frac      = tagh_pair_frac[counter-1];
			loc_pair_frac_err  = tagh_pair_fracE[counter-1];
			loc_pair_frac_exp  = tagh_pair_frac_exp[counter-1];
		} else {
			int counter        = tagm_counter_vec[ib-n_bins1];
			loc_eb             = tagm_en[counter-1];
			loc_pair_frac      = tagm_pair_frac[counter-1];
			loc_pair_frac_err  = tagm_pair_fracE[counter-1];
			loc_pair_frac_exp  = tagm_pair_frac_exp[counter-1];
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
