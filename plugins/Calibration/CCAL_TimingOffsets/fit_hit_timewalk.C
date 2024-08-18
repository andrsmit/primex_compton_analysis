
double pw_x[4] = {0.5, 1.0, 2.0, 4.0};

Double_t tw_fit(Double_t *x, Double_t *par);

void fit_hit_timewalk(double min_fit_x = 0.0, double max_fit_x = 10.0) {
	
	TFile *loc_fIn = new TFile("rootFiles/061321_v2.root", "READ");
	
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
		energy_err[ibin] = 0.;//h2->GetXaxis()->GetBinCenter(2)-h2->GetXaxis()->GetBinCenter(1);
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
	
	TF1 *fMean = new TF1("fMean", tw_fit, min_fit_x, max_fit_x, 16);
	fMean->SetParameters(1.0, 1.0, -2.0, 0.0);
	for(int ipar=4; ipar<13; ipar++) {
		fMean->FixParameter(ipar, 0.);
	}
	fMean->SetRange(min_fit_x, pw_x[0]);
	gMean->Fit(fMean, "R0");
	
	//-----------------------------------//
	// pw_x1 to pw_x2:
	
	for(int ipar=0; ipar<4; ipar++) {
		fMean->FixParameter(ipar, fMean->GetParameter(ipar));
	}
	for(int ipar=4; ipar<7; ipar++) {
		fMean->ReleaseParameter(ipar);
	}
	fMean->SetParameter(4,  1.0);
	fMean->SetParameter(5,  1.0);
	fMean->SetParameter(6, -2.0);
	
	fMean->SetRange(pw_x[0], pw_x[1]);
	gMean->Fit(fMean, "R0");
	
	//-----------------------------------//
	// pw_x2 to pw_x3:
	
	for(int ipar=4; ipar<7; ipar++) {
		fMean->FixParameter(ipar, fMean->GetParameter(ipar));
	}
	for(int ipar=7; ipar<10; ipar++) {
		fMean->ReleaseParameter(ipar);
	}
	fMean->SetParameter(7, fMean->GetParameter(4));
	fMean->SetParameter(8, fMean->GetParameter(5));
	fMean->SetParameter(9, fMean->GetParameter(6));
	
	fMean->SetRange(pw_x[1], pw_x[2]);
	gMean->Fit(fMean, "R0");
	
	//-----------------------------------//
	// pw_x3 to pw_x4:
	
	for(int ipar=7; ipar<10; ipar++) {
		fMean->FixParameter(ipar, fMean->GetParameter(ipar));
	}
	for(int ipar=10; ipar<13; ipar++) {
		fMean->ReleaseParameter(ipar);
	}
	fMean->SetParameter(10, fMean->GetParameter(7));
	fMean->SetParameter(11, fMean->GetParameter(8));
	fMean->SetParameter(12, fMean->GetParameter(9));
	
	fMean->SetRange(pw_x[2], pw_x[3]);
	gMean->Fit(fMean, "R0");
	
	//-----------------------------------//
	// pw_x4 to max_fit_x:
	
	for(int ipar=10; ipar<13; ipar++) {
		fMean->FixParameter(ipar, fMean->GetParameter(ipar));
	}
	for(int ipar=13; ipar<16; ipar++) {
		fMean->ReleaseParameter(ipar);
	}
	fMean->SetParameter(13, fMean->GetParameter(10));
	fMean->SetParameter(14, fMean->GetParameter(11));
	fMean->SetParameter(15, fMean->GetParameter(12));
	
	fMean->SetRange(pw_x[3], max_fit_x);
	gMean->Fit(fMean, "R0");
	
	fMean->SetRange(min_fit_x, max_fit_x);
	fMean->Draw("same");
	
	double fit_pars[5][4];
	
	fit_pars[0][0] = fMean->GetParameter(0);
	fit_pars[0][1] = fMean->GetParameter(1);
	fit_pars[0][2] = fMean->GetParameter(2);
	fit_pars[0][3] = fMean->GetParameter(3);
	
	fit_pars[1][0] = fMean->GetParameter(4);
	fit_pars[1][1] = fMean->GetParameter(5);
	fit_pars[1][2] = fMean->GetParameter(6);
	
	fit_pars[2][0] = fMean->GetParameter(7);
	fit_pars[2][1] = fMean->GetParameter(8);
	fit_pars[2][2] = fMean->GetParameter(9);
	
	fit_pars[3][0] = fMean->GetParameter(10);
	fit_pars[3][1] = fMean->GetParameter(11);
	fit_pars[3][2] = fMean->GetParameter(12);
	
	fit_pars[4][0] = fMean->GetParameter(13);
	fit_pars[4][1] = fMean->GetParameter(14);
	fit_pars[4][2] = fMean->GetParameter(15);
	
	for(int ipar=1; ipar<5; ipar++) {
		fit_pars[ipar][3] = 
			(fit_pars[ipar-1][0]*exp(fit_pars[ipar-1][1] + fit_pars[ipar-1][2]*pw_x[ipar-1]) + fit_pars[ipar-1][3]) 
			- (fit_pars[ipar][0]*exp(fit_pars[ipar][1] + fit_pars[ipar][2]*pw_x[ipar-1]));
	}
	
	ofstream outf("tw_correction.txt");
	for(int irow=0; irow<5; irow++) {
		char buf[256];
		sprintf(buf, "%f  %f  %f  %f", fit_pars[irow][0], fit_pars[irow][1], fit_pars[irow][2], fit_pars[irow][3]);
		//sprintf(buf, "%f  %f  %f  %f", 0., 0., 0., 0.);
		outf << buf << "\n";
	}
	outf.close();
	
	// double check we calculated the parameters correctly:
	
	TF1 *f_compare[5];
	for(int i=0; i<5; i++) {
		double   min_x;
		if(i==0) min_x = min_fit_x;
		else     min_x = pw_x[i-1];
		double   max_x;
		if(i==4) max_x = max_fit_x;
		else     max_x = pw_x[i];
		f_compare[i] = new TF1(Form("f_compare_%d",i), "[0]*exp([1] + [2]*x) + [3]", min_x, max_x);
		f_compare[i]->SetParameters(fit_pars[i][0], fit_pars[i][1], fit_pars[i][2], fit_pars[i][3]);
		f_compare[i]->SetLineColor(i+3);
		f_compare[i]->SetLineStyle(2);
		f_compare[i]->Draw("same");
	}
	
	//-----------------------------------//
	// make a plot of the residuals:
	
	double *res     = new double[n_bins];
	double *res_err = new double[n_bins];
	for(int ibin=0; ibin<n_bins; ibin++) {
		res[ibin]     = mean[ibin] - fMean->Eval(energy[ibin]);
		res_err[ibin] = mean_err[ibin];
	}
	
	TGraphErrors *gRes = new TGraphErrors(n_bins, energy, res, energy_err, res_err);
	gRes->GetXaxis()->SetTitle("CCAL Hit Energy (GeV)");
	gRes->GetXaxis()->SetTitleSize(0.05);
	gRes->GetXaxis()->SetTitleOffset(1.0);
	gRes->GetXaxis()->CenterTitle(true);
	gRes->GetYaxis()->SetTitle("Measured - Fit");
	gRes->GetYaxis()->SetTitleSize(0.05);
	gRes->GetYaxis()->SetTitleOffset(1.0);
	gRes->GetYaxis()->CenterTitle(true);
	gRes->SetTitle("CCAL Timewalk Correction Residuals");
	gRes->SetMarkerStyle(8);
	gRes->SetMarkerSize(0.7);
	gRes->SetMarkerColor(kBlue+1);
	gRes->SetLineColor(kBlue);
	
	TCanvas *cRes = new TCanvas("cRes", "cRes", 700, 500);
	cRes->SetLeftMargin(0.12);
	cRes->SetRightMargin(0.08);
	cRes->SetBottomMargin(0.12);
	cRes->SetGrid();
	cRes->SetTickx(); cRes->SetTicky();
	gRes->Draw("AP");
	
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
	
	Double_t p8  = par[7];
	Double_t p9  = par[8];
	Double_t p10 = par[9];
	
	Double_t p12 = par[10];
	Double_t p13 = par[11];
	Double_t p14 = par[12];
	
	Double_t p16 = par[13];
	Double_t p17 = par[14];
	Double_t p18 = par[15];
	
	// four continuous exponentials:
	Double_t p7  =  p0*exp(p1 +p2 *pw_x[0]) + p3  - p4 *exp(p5  + p6 *pw_x[0]);
	Double_t p11 =  p4*exp(p5 +p6 *pw_x[1]) + p7  - p8 *exp(p9  + p10*pw_x[1]);
	Double_t p15 =  p8*exp(p9 +p10*pw_x[2]) + p11 - p12*exp(p13 + p14*pw_x[2]);
	Double_t p19 = p12*exp(p13+p14*pw_x[3]) + p15 - p16*exp(p17 + p18*pw_x[3]);
	
	Double_t f;
	if(xx<pw_x[0]) {
		f = p0*exp(p1 + p2*xx) + p3;
	} else if(xx<pw_x[1]) {
		f = p4*exp(p5 + p6*xx) + p7;
	} else if(xx<pw_x[2]) {
		f = p8*exp(p9 + p10*xx) + p11;
	} else if(xx<pw_x[3]) {
		f = p12*exp(p13 + p14*xx) + p15;
	} else {
		f = p16*exp(p17 + p18*xx) + p19;
	}
	
	return f;
}
