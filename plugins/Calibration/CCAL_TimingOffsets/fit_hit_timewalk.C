
double pw_x = 1.8;

Double_t tw_fit(Double_t *x, Double_t *par);

void fit_hit_timewalk(double min_fit_x = 0.05, double max_fit_x = 8.0) {
	
	TFile *loc_fIn = new TFile("hd_root.root", "READ");
	
	TH2F *h2 = (TH2F*)loc_fIn->Get("CCAL_TimingOffsets/ccal_rf_dt_vs_E");
	h2->SetDirectory(0);
	
	loc_fIn->Close();
	
	int min_bin_fit = h2->GetXaxis()->FindBin(min_fit_x);
	int max_bin_fit = h2->GetXaxis()->FindBin(max_fit_x);
	
	vector<double> energy_vec, mean_vec, mean_err_vec;
	for(int ibin=min_bin_fit; ibin<=max_bin_fit; ibin++) {
		
		TH1F *loc_h1 = (TH1F*)h2->ProjectionY("loc_h1",ibin,ibin);
		if(loc_h1->Integral()<20.) {
			loc_h1->Delete();
			continue;
		}
		// fit peak to a Gaussian to estimate peak position:
		loc_h1->GetXaxis()->SetRangeUser(-2.0,2.0);
		double mean_guess = loc_h1->GetBinCenter(loc_h1->GetMaximumBin());
		double sigma_guess = 0.2;
		TF1 *loc_f1 = new TF1("loc_f1", "gaus", -2.0, 2.0);
		loc_f1->SetParameters(loc_h1->GetMaximum(), mean_guess, sigma_guess);
		loc_f1->SetRange(mean_guess-1.0, mean_guess+1.0);
		loc_h1->GetXaxis()->SetRangeUser(-5.0, 5.0);
		loc_h1->Fit(loc_f1, "R0QL");
		loc_f1->SetRange(loc_f1->GetParameter(1)-1.5*loc_f1->GetParameter(2), 
			loc_f1->GetParameter(1)+1.5*loc_f1->GetParameter(2));
		loc_h1->Fit(loc_f1, "R0QL");
		
		energy_vec.push_back(h2->GetXaxis()->GetBinCenter(ibin));
		mean_vec.push_back(loc_f1->GetParameter(1));
		mean_err_vec.push_back(loc_f1->GetParError(1));
		
		loc_h1->Delete();
		loc_f1->Delete();
	}
	
	int n_bins = (int)energy_vec.size();
	double *energy     = new double[n_bins];
	double *energy_err = new double[n_bins];
	double *mean       = new double[n_bins];
	double *mean_err   = new double[n_bins];
	for(int ibin=0; ibin<n_bins; ibin++) {
		energy[ibin]     = energy_vec[ibin];
		energy_err[ibin] = h2->GetXaxis()->GetBinCenter(2)-h2->GetXaxis()->GetBinCenter(1);
		mean[ibin]       = mean_vec[ibin];
		mean_err[ibin]   = mean_err_vec[ibin];
	}
	
	TGraphErrors *gMean = new TGraphErrors(n_bins, energy, mean, energy_err, mean_err);
	gMean->GetXaxis()->SetTitle("CCAL Hit Energy (GeV)");
	gMean->GetXaxis()->SetTitleSize(0.05);
	gMean->GetXaxis()->SetTitleOffset(1.0);
	gMean->GetXaxis()->CenterTitle(true);
	gMean->GetYaxis()->SetTitle("t_{CCAL} - t_{RF} (ns)");
	gMean->GetYaxis()->SetTitleSize(0.05);
	gMean->GetYaxis()->SetTitleOffset(1.0);
	gMean->GetYaxis()->CenterTitle(true);
	gMean->SetTitle("CCAL Hit Timewalk");
	gMean->SetMarkerStyle(8);
	gMean->SetMarkerSize(0.7);
	gMean->SetMarkerColor(kBlue+1);
	gMean->SetLineColor(kBlue);
	
	TCanvas *cMean = new TCanvas("cMean", "cMean", 700, 500);
	cMean->SetLeftMargin(0.12);
	cMean->SetRightMargin(0.08);
	cMean->SetBottomMargin(0.12);
	cMean->SetGrid();
	cMean->SetTickx(); cMean->SetTicky();
	gMean->Draw("AP");
	
	TF1 *fMean = new TF1("fMean", tw_fit, min_fit_x, max_fit_x, 7);
	fMean->SetParameters(0.5, -1.0, -1.0, 0.0);
	for(int ipar=4; ipar<7; ipar++) {
		fMean->FixParameter(ipar, 0.);
	}
	fMean->SetRange(min_fit_x, pw_x);
	gMean->Fit(fMean, "R0");
	
	for(int ipar=0; ipar<4; ipar++) {
		fMean->FixParameter(ipar, fMean->GetParameter(ipar));
	}
	for(int ipar=4; ipar<7; ipar++) {
		fMean->ReleaseParameter(ipar);
	}
	fMean->SetParameter(4,  1.0);
	fMean->SetParameter(5, -1.0);
	fMean->SetParameter(6, -1.0);
	
	fMean->SetRange(pw_x, max_fit_x);
	gMean->Fit(fMean, "R0");
	
	fMean->SetRange(min_fit_x, max_fit_x);
	fMean->Draw("same");
	
	/*
	TF1 *fMean[2];
	fMean[0] = new TF1("fMean0", "[0]*exp((x-[1])*[2]) + [3]", min_fit_x,   1.);
	fMean[1] = new TF1("fMean1", "pol3", 1.,   max_fit_x);
	
	fMean[0]->SetParameters(0.5, 1.0, -1.0, 0.0);
	//fMean[1]->SetParameters(1.0, 0.0, -2.0, -1.0);
	
	gMean->Fit(fMean[0], "R0");
	gMean->Fit(fMean[1], "R0");
	
	fMean[0]->Draw("same");
	fMean[1]->Draw("same");
	*/
	
	cout << fMean->Eval(0.015) << endl;
	
	return;
}


Double_t tw_fit(Double_t *x, Double_t *par) {
	
	Double_t xx = x[0];
	
	Double_t p0 = par[0];
	Double_t p1 = par[1];
	Double_t p2 = par[2];
	Double_t p3 = par[3];
	Double_t p4 = par[4];
	Double_t p5 = par[5];
	Double_t p6 = par[6];
	//Double_t p7 = par[7];
	Double_t p7 = p0*exp(p1+p2*pw_x) - p4*exp(p5+p6*pw_x) + p3;
	
	Double_t f;
	if(xx<pw_x) {
		f = p0*exp(p1 + p2*xx) + p3;
	} else {
		f = p4*exp(p5 + p6*xx) + p7;
	}
	
	return f;
}
