
bool isNaN(double x) {
	return x != x;
}

void plot_chi2(vector<int> tagh_counter_vec, vector<int> tagm_counter_vec) {
	
	int n_bins1 = (int)tagh_counter_vec.size();
	int n_bins2 = (int)tagm_counter_vec.size();
	int n_bins  = n_bins1+n_bins2;
	
	double *beam_energy = new double[n_bins];
	double *chi2        = new double[n_bins];
	
	for(int ib=0; ib<n_bins; ib++) {
		double loc_eb   = 0.;
		double loc_chi2 = 0.;
		if(ib<n_bins1) {
			int counter = tagh_counter_vec[ib];
			loc_eb      = tagh_en[counter-1];
			loc_chi2    = tagh_chi2[counter-1];
		} else {
			int counter = tagm_counter_vec[ib-n_bins1];
			loc_eb      = tagm_en[counter-1];
			loc_chi2    = tagm_chi2[counter-1];
		}
		if(isNaN(loc_chi2)) loc_chi2 = 0.;
		beam_energy[ib] = loc_eb;
		chi2[ib]        = loc_chi2;
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
	
	return;
}
