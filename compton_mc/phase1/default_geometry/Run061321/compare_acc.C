char tagh_xscale_fname1[256], tagm_xscale_fname1[256];
char tagh_xscale_fname2[256], tagm_xscale_fname2[256];
char        hname_tagh[256], hname_tagm[256];
char       target_char[256];
char          mc_dir_1[256], mc_dir_2[256];

int target_type, beam_current;
double endpoint_energy1, endpoint_energy_calib1;
double endpoint_energy2, endpoint_energy_calib2;
char name1[256], name2[256];

double   tagh_en1[274],   tagm_en1[102];
double   tagh_en2[274],   tagm_en2[102];
double  tagh_acc1[274],  tagm_acc1[102];
double tagh_acc1E[274], tagm_acc1E[102];
double  tagh_acc2[274],  tagm_acc2[102];
double tagh_acc2E[274], tagm_acc2E[102];

double bin_size = 16. / 4000.;
int rebins, n_mev;

bool DRAW_ACCEPTANCE;
bool CALC_ACC_FROM_FIT;

TF1 *f_acc_1, *f_acc_2;
TCanvas *canvas_acc;

TLegend *leg;

//----------   Function Declarations   ----------//

void get_counter_energies1();
void get_counter_energies2();
void get_acc_1();
void get_acc_2();

Double_t        bkgd_fit(Double_t *x, Double_t *par);
Double_t double_gaus_fit(Double_t *x, Double_t *par);

//-----------------------------------------------//

void compare_acc()
{
	gStyle->SetOptStat(0); gStyle->SetOptFit(0);
	
	char loc_hname_tagh[256], loc_hname_tagm[256];
	sprintf(loc_hname_tagh, "deltaK/deltaK_tagh_0");
	sprintf(loc_hname_tagm, "deltaK/deltaK_tagm_0");
	
	//----------   Initialize   ----------//
	
	if(true) {
		
		// Be Target:
		
		target_type  =   0;
		beam_current = 200;
		sprintf(target_char, "^{9}Be");
	} else {
		
	 	// He Target:
		
		target_type  =   1;
		beam_current = 200;
		sprintf(target_char, "^{4}He");
	}
	
	rebins    =  5;
	bin_size *=  (double)rebins;
	n_mev     =  4 * rebins;
	
	// Adjustable switches:
	
	DRAW_ACCEPTANCE   = true;
	CALC_ACC_FROM_FIT = false;
	
	//-----------------------------------------------------------------------------//
	
	const char pathName[256] = "/work/halld/home/andrsmit/primex_compton_analysis";
	
	if(!target_type) {
		
		//endpoint_energy_calib1 = 11.6061;
		//endpoint_energy1       = 11.608;
		//sprintf(mc_dir_1, "%s/compton_mc/phaseI/small_beam_spot/Be/recRootFiles_sys_new", 
		//	pathName);
		
		endpoint_energy_calib1 = 11.6061;
		endpoint_energy1       = 11.6061;
		sprintf(mc_dir_1, "%s/compton_mc/phaseI/default_geometry/Be/recRootFiles_TOF", 
			pathName);
		sprintf(name1, "Be Target, no TOF Veto");
		
		endpoint_energy_calib2 = 11.6061;
		endpoint_energy2       = 11.6061;
		sprintf(mc_dir_2, "%s/compton_mc/phaseI/default_geometry/Be/recRootFiles_TOF", 
			pathName);
		sprintf(name2, "Be Target, TOF Veto applied");
		
		
	} else {
		
		return;
	}
	
	// file containing the xscales for the tagh and tagm counters:
	
	sprintf(tagh_xscale_fname1, "%s/photon_flux/phase1/primex_tagh.txt", pathName);
	sprintf(tagm_xscale_fname1, "%s/photon_flux/phase1/primex_tagm.txt", pathName);
	sprintf(tagh_xscale_fname2, "%s/photon_flux/phase1/primex_tagh.txt", pathName);
	sprintf(tagm_xscale_fname2, "%s/photon_flux/phase1/primex_tagm.txt", pathName);
	
	//-----------------------------------------------------------------------------//
	
	for( int i = 0; i < 274; i++ ) {
		
		tagh_en1[i]   = 0.;
		tagh_en2[i]   = 0.;
		tagh_acc1[i]  = 0.;
		tagh_acc1E[i] = 0.;
		tagh_acc2[i]  = 0.;
		tagh_acc2E[i] = 0.;
	}
	for( int i = 0; i < 102; i++ ) {
		
		tagm_en1[i]   = 0.;
		tagm_en2[i]   = 0.;
		tagm_acc1[i]  = 0.;
		tagm_acc1E[i] = 0.;
		tagm_acc2[i]  = 0.;
		tagm_acc2E[i] = 0.;
	}
	
	//-----------------------------------------------------------------------------//
	
	leg = new TLegend(0.508,0.652,0.834,0.831);
	leg->SetBorderSize(0);
	
	get_counter_energies1();
	get_counter_energies2();
	sprintf(hname_tagh, "compton_analysis_TOF/deltaK_cut/deltaK_tagh_cut_0");
	sprintf(hname_tagm, "compton_analysis_TOF/deltaK_cut/deltaK_tagm_cut_0");
	get_acc_1();
	sprintf(hname_tagh, "compton_analysis_TOF/deltaK_cut/deltaK_tagh_cut_1");
	sprintf(hname_tagm, "compton_analysis_TOF/deltaK_cut/deltaK_tagm_cut_1");
	get_acc_2();
	
	canvas_acc->cd();
	leg->Draw();
	
	double xiter = 6.0;
	vector<double> xvals, yvals;
	while(xiter<tagh_en1[0]) {
		double yiter = f_acc_2->Eval(xiter) / f_acc_1->Eval(xiter);
		xvals.push_back(xiter);
		yvals.push_back(yiter);
		xiter += 0.01;
	}
	
	int n_bins = (int)xvals.size();
	double *x_val = new double[n_bins];
	double *y_val = new double[n_bins];
	for(int i=0; i<n_bins; i++) {
		x_val[i] = xvals[i];
		y_val[i] = yvals[i];
	}
	
	TGraph *g = new TGraph(n_bins, x_val, y_val);
	g->SetMarkerStyle(8);
	g->SetMarkerSize(0.5);
	g->GetXaxis()->SetTitle("Photon Beam Energy [GeV]");
	g->GetXaxis()->SetTitleSize(0.05);
	g->GetXaxis()->SetTitleOffset(0.8);
	g->GetYaxis()->SetTitle("#epsilon_{Veto} / #epsilon_{No Veto}");
	g->GetYaxis()->SetTitleSize(0.05);
	g->GetYaxis()->SetTitleOffset(0.9);
	g->GetYaxis()->SetRangeUser(0.9,1.1);
	
	g->SetTitle("Compton Acceptance Ratio (With TOF Veto / Without Veto)");
	
	TF1 *ff = new TF1("ff","pol5",5.0,12.0);
	ff->SetLineColor(kRed);
	g->Fit(ff,"R0");
	
	TCanvas *cc = new TCanvas("cc","cc",1000,600);
	g->Draw("AP");
	ff->Draw("same");
	
	return;
}


void get_counter_energies1() 
{
	
	int a; double b, c;
	
	ifstream inf1( tagh_xscale_fname1 );
	for( int i=0; i<274; i++ ) {
		inf1 >> a >> b >> c;
		double deltaE  = endpoint_energy1 - endpoint_energy_calib1;
		double emin    = b * endpoint_energy_calib1  +  deltaE;
		double emax    = c * endpoint_energy_calib1  +  deltaE;
		tagh_en1[i] = 0.5 * (emin + emax);
	}
	inf1.close();
	
	ifstream inf2( tagm_xscale_fname1 );
	for( int i=0; i<102; i++ ) {
		inf2 >> a >> b >> c;
		double deltaE  = endpoint_energy1 - endpoint_energy_calib1;
		double emin    = b * endpoint_energy_calib1  +  deltaE;
		double emax    = c * endpoint_energy_calib1  +  deltaE;
		tagm_en1[i] = 0.5 * (emin + emax);
	}
	inf2.close();
	
	
	return;
}

void get_counter_energies2() 
{
	
	int a; double b, c;
	
	ifstream inf1( tagh_xscale_fname2 );
	for( int i=0; i<274; i++ ) {
		inf1 >> a >> b >> c;
		double deltaE  = endpoint_energy2 - endpoint_energy_calib2;
		double emin    = b * endpoint_energy_calib2  +  deltaE;
		double emax    = c * endpoint_energy_calib2  +  deltaE;
		tagh_en2[i] = 0.5 * (emin + emax);
	}
	inf1.close();
	
	ifstream inf2( tagm_xscale_fname2 );
	for( int i=0; i<102; i++ ) {
		inf2 >> a >> b >> c;
		double deltaE  = endpoint_energy2 - endpoint_energy_calib2;
		double emin    = b * endpoint_energy_calib2  +  deltaE;
		double emax    = c * endpoint_energy_calib2  +  deltaE;
		tagm_en2[i] = 0.5 * (emin + emax);
	}
	inf2.close();
	
	
	return;
}



Double_t bkgd_fit( Double_t *x, Double_t *par ) {
	
	Double_t xx = x[0];
	
	
	if( xx > -1.5 && xx < 1.0 ) {
		TF1::RejectPoint();
		return 0.;
	}
	
	Double_t f = par[0] + 
		     par[1] * x[0] +
		     par[2] * x[0] * x[0] + 
		     par[3] * x[0] * x[0] * x[0];
	
	return f;
}



Double_t double_gaus_fit( Double_t *x, Double_t *par ) 
{
	
	Double_t xx  =   x[0];
	
	Double_t N    = par[0];
	Double_t z    = par[1];
	Double_t mu1  = par[2];
	Double_t dmu  = par[3];
	Double_t sig1 = par[4];
	Double_t sig2 = par[5];
	
	Double_t p0   = par[6];
	Double_t p1   = par[7];
	Double_t p2   = par[8];
	Double_t p3   = par[9];
	
	Double_t mu2  = mu1 - dmu;
	
	Double_t f = N * ( (1.-z)*exp( -0.5*pow((xx-mu1)/sig1,2.0) ) 
		+ z*exp( -0.5*pow((xx-mu2)/sig2,2.0) ) );
	
	f += 1.0*( p0 + p1*xx + p2*xx*xx + p3*xx*xx*xx );
	
	
	return f;
}


void get_acc_1()
{
	vector<double> tagh_enVec, tagh_accVec, tagh_accEVec;
	vector<double> tagm_enVec, tagm_accVec, tagm_accEVec;
	
	for(int tagh_counter = 1; tagh_counter <= 274; tagh_counter++) {
		
		char fname[256];
		sprintf(fname, "%s/tagh_%03d.root", mc_dir_1, tagh_counter);
		
		if(gSystem->AccessPathName(fname)) continue;
		
		//cout << "Processing TAGH Counter " << tagh_counter << endl;
		
		double loc_eb = tagh_en1[tagh_counter-1];
		
		TFile *fSim = new TFile(fname, "READ");
		
		TH1F *h_vertex = (TH1F*)fSim->Get(
			"compton_analysis_TOF/vertex_accepted")->Clone(
			Form("h_vertex_tagh_%d",tagh_counter));
		
		double n_gen = h_vertex->Integral();
		
		//-----------------------------------------------------//
		
		TH2F *h2_tagh = new TH2F(Form("h2_tagh_tagh_sim_%d",tagh_counter), "DeltaK", 
			274, 0.5, 274.5, 4000, -8.0, 8.0);
		TH2F *h2_tagm = new TH2F(Form("h2_tagh_tagm_sim_%d",tagh_counter), "DeltaK", 
			102, 0.5, 102.5, 4000, -8.0, 8.0);
		
		h2_tagh->Add((TH2F*)fSim->Get(Form("%s",hname_tagh)));
		h2_tagm->Add((TH2F*)fSim->Get(Form("%s",hname_tagm)));
		
		h2_tagh->SetDirectory(0);
		h2_tagm->SetDirectory(0);
		
		fSim->Close();
		
		TH1F *h1 = (TH1F*)h2_tagh->ProjectionY(Form("h1_sim_tagh_tagh_%d",tagh_counter));
		h1->Add((TH1F*)h2_tagm->ProjectionY(Form("h1_sim_tagh_tagm_%d",tagh_counter)));
		if(h1->Integral() < 1.e1) continue;
		
		h1->Rebin(rebins);
		h1->GetXaxis()->SetTitle("E_{Comp} - E_{#gamma} [GeV]");
		h1->GetXaxis()->SetTitleOffset( 1.1 );
		h1->GetXaxis()->SetRangeUser(-1.5, 1.5);
		h1->GetYaxis()->SetTitle( Form("counts / %d MeV", n_mev));
		h1->GetYaxis()->SetTitleOffset(1.3);
		h1->SetLineColor(kBlack);
		h1->SetLineWidth(2);
		h1->GetXaxis()->SetRangeUser(-3., 2.);
		
		double yield;
		
		if(CALC_ACC_FROM_FIT) {
			
			TF1 *f_gaus = new TF1(Form("f_gaus_tagh_acc_%d",tagh_counter), 
				"gaus", -2., 2.);
			
			f_gaus->SetParameters(h1->GetMaximum(), 
				h1->GetBinCenter(h1->GetMaximumBin()), 0.2);
			f_gaus->SetParLimits(2, 0., 1.);
			
			h1->Fit(f_gaus, "R0QL");
			f_gaus->SetRange(f_gaus->GetParameter(1)-0.2, f_gaus->GetParameter(1)+0.2);
			h1->Fit(f_gaus, "R0QL");
			
			TF1 *f_fit = new TF1(Form("f_fit_tagh_acc_%d",tagh_counter), 
				double_gaus_fit, -2., 2., 10);	
			
			f_fit->SetParName(0, "A");
			f_fit->SetParName(1, "z");
			f_fit->SetParName(2, "#mu_{1}");
			f_fit->SetParName(3, "#mu_{1}-#mu_{2}");
			f_fit->SetParName(4, "#sigma_{1}");
			f_fit->SetParName(5, "#sigma_{2}");
			
			f_fit->SetParLimits(1, 0.0, 0.5);
			f_fit->SetParLimits(4, 0.0, 1.0);
			f_fit->SetParLimits(5, 0.0, 2.0);
			
			f_fit->SetParameters(f_gaus->GetParameter(0), 0., f_gaus->GetParameter(1), 
				0., f_gaus->GetParameter(2), 2.*f_gaus->GetParameter(2));
			
			f_fit->FixParameter(6, 0.0);
			f_fit->FixParameter(7, 0.0);
			f_fit->FixParameter(8, 0.0);
			f_fit->FixParameter(9, 0.0);
			
			f_fit->SetRange(-2.0, 2.0);
			h1->Fit(f_fit, "R0QL");
			f_fit->SetRange(-2., 1.5);
			h1->Fit(f_fit, "R0QL");
			
			yield  = f_fit->GetParameter(0) * sqrt(2.*TMath::Pi()) 
				* ((1.-f_fit->GetParameter(1))*f_fit->GetParameter(4) 
				+ f_fit->GetParameter(1)*f_fit->GetParameter(5));
			yield /= bin_size;
		} else {
			
			yield = h1->Integral();
		}
		
		double loc_acc  = yield / n_gen;
		double loc_accE = sqrt(n_gen*loc_acc*(1.-loc_acc)) / n_gen;
		
		tagh_acc1[tagh_counter-1]  = loc_acc;
		tagh_acc1E[tagh_counter-1] = loc_accE;
		
		tagh_enVec.push_back(loc_eb);
		tagh_accVec.push_back(loc_acc);
		tagh_accEVec.push_back(loc_accE);
	}
	
	for(int tagm_counter = 1; tagm_counter <= 102; tagm_counter++) {
		
		if(tagm_counter>=71 && tagm_counter<=73) continue;
		
		char fname[256];
		sprintf(fname, "%s/tagm_%03d.root", mc_dir_1, tagm_counter);
		
		if(gSystem->AccessPathName(fname)) continue;
		
		//cout << "Processing TAGM Counter " << tagm_counter << endl;
		
		double loc_eb = tagm_en1[tagm_counter-1];
		
		TFile *fSim = new TFile(fname, "READ");
		
		TH1F *h_vertex = (TH1F*)fSim->Get("compton_analysis_TOF/vertex_accepted")->Clone(
			Form("h_vertex_tagm_%d",tagm_counter));
		
		double n_gen = h_vertex->Integral();
		
		//-----------------------------------------------------//
		
		TH2F *h2_tagh = new TH2F(Form("h2_tagm_tagh_sim_%d",tagm_counter), "DeltaK", 
			274, 0.5, 274.5, 4000, -8.0, 8.0);
		TH2F *h2_tagm = new TH2F(Form("h2_tagm_tagm_sim_%d",tagm_counter), "DeltaK", 
			102, 0.5, 102.5, 4000, -8.0, 8.0);
		
		h2_tagh->Add((TH2F*)fSim->Get(Form("%s",hname_tagh)));
		h2_tagm->Add((TH2F*)fSim->Get(Form("%s",hname_tagm)));
		
		h2_tagh->SetDirectory(0);
		h2_tagm->SetDirectory(0);
		
		fSim->Close();
		
		TH1F *h1 = (TH1F*)h2_tagh->ProjectionY(Form("h1_sim_tagm_tagh_%d",tagm_counter));
		h1->Add((TH1F*)h2_tagm->ProjectionY(Form("h1_sim_tagm_tagm_%d",tagm_counter)));
		if(h1->Integral() < 1.e1) continue;
		
		h1->Rebin(rebins);
		h1->GetXaxis()->SetTitle("E_{Comp} - E_{#gamma} [GeV]");
		h1->GetXaxis()->SetTitleOffset(1.1);
		h1->GetXaxis()->SetRangeUser(-1.5, 1.5);
		h1->GetYaxis()->SetTitle(Form("counts / %d MeV", n_mev));
		h1->GetYaxis()->SetTitleOffset(1.3);
		h1->SetLineColor(kBlack);
		h1->SetLineWidth(2);
		h1->GetXaxis()->SetRangeUser(-3., 2.);
		
		double yield;
		
		if(CALC_ACC_FROM_FIT) {
			
			TF1 *f_gaus = new TF1(Form("f_gaus_tagm_acc_%d",tagm_counter), 
				"gaus", -2., 2.);
			
			f_gaus->SetParameters(h1->GetMaximum(), 
				h1->GetBinCenter(h1->GetMaximumBin()), 0.2);
			f_gaus->SetParLimits(2, 0., 1.);
			
			h1->Fit(f_gaus, "R0QL");
			f_gaus->SetRange(f_gaus->GetParameter(1)-0.2, f_gaus->GetParameter(1)+0.2);
			h1->Fit(f_gaus, "R0QL");
			
			TF1 *f_fit = new TF1(Form("f_fit_tagm_acc_%d",tagm_counter), 
				double_gaus_fit, -2., 2., 10);	
			
			f_fit->SetParName(0, "A");
			f_fit->SetParName(1, "z");
			f_fit->SetParName(2, "#mu_{1}");
			f_fit->SetParName(3, "#mu_{1}-#mu_{2}");
			f_fit->SetParName(4, "#sigma_{1}");
			f_fit->SetParName(5, "#sigma_{2}");
			
			f_fit->SetParLimits(1, 0.0, 0.5);
			f_fit->SetParLimits(4, 0.0, 1.0);
			f_fit->SetParLimits(5, 0.0, 2.0);
			
			f_fit->SetParameters(f_gaus->GetParameter(0), 0., f_gaus->GetParameter(1), 
				0., f_gaus->GetParameter(2), 2.*f_gaus->GetParameter(2));
			
			f_fit->FixParameter(6, 0.0);
			f_fit->FixParameter(7, 0.0);
			f_fit->FixParameter(8, 0.0);
			f_fit->FixParameter(9, 0.0);
			
			f_fit->SetRange(-2.0, 2.0);
			h1->Fit(f_fit, "R0QL");
			f_fit->SetRange(-2., 1.5);
			h1->Fit(f_fit, "R0QL");
			
			yield  = f_fit->GetParameter(0) * sqrt(2.*TMath::Pi()) 
				* ((1.-f_fit->GetParameter(1))*f_fit->GetParameter(4) 
				+ f_fit->GetParameter(1)*f_fit->GetParameter(5));
			yield /= bin_size;
		} else {
			
			yield = h1->Integral();
		}
		
		double loc_acc  = yield / n_gen;
		double loc_accE = sqrt(n_gen*loc_acc*(1.-loc_acc)) / n_gen;
		
		tagm_acc1[tagm_counter-1]  = loc_acc;
		tagm_acc1E[tagm_counter-1] = loc_accE;
		
		tagm_enVec.push_back(loc_eb);
		tagm_accVec.push_back(loc_acc);
		tagm_accEVec.push_back(loc_accE);
	}
	
	
	int n_bins1 = (int)tagh_enVec.size();
	double *energy1 = new double[n_bins1];
	double *zero1   = new double[n_bins1];
	double *acc1    = new double[n_bins1];
	double *acc1E   = new double[n_bins1];
	for( int i=0; i<n_bins1; i++ ) {
		energy1[i] = tagh_enVec[i];
		zero1[i]   = 0.;
		acc1[i]    = tagh_accVec[i];
		acc1E[i]   = tagh_accEVec[i];
	}
	
	int n_bins2 = (int)tagm_enVec.size();
	double *energy2 = new double[n_bins2];
	double *zero2   = new double[n_bins2];
	double *acc2    = new double[n_bins2];
	double *acc2E   = new double[n_bins2];
	for( int i=0; i<n_bins2; i++ ) {
		energy2[i] = tagm_enVec[i];
		zero2[i]   = 0.;
		acc2[i]    = tagm_accVec[i];
		acc2E[i]   = tagm_accEVec[i];
	}
	
	int n_bins = n_bins1 + n_bins2;
	double *energy  = new double[n_bins];
	double *zero    = new double[n_bins];
	double *acc     = new double[n_bins];
	double *accE    = new double[n_bins];
	for( int i=0; i<n_bins; i++ ) {
		if( i<n_bins1 ) {
			energy[i] = tagh_enVec[i];
			zero[i]   = 0.;
			acc[i]    = tagh_accVec[i];
			accE[i]   = tagh_accEVec[i];
		} else {
			energy[i] = tagm_enVec[i-n_bins1];
			zero[i]   = 0.;
			acc[i]    = tagm_accVec[i-n_bins1];
			accE[i]   = tagm_accEVec[i-n_bins1];
		}
	}
	
	TGraphErrors *gAcc = new TGraphErrors(n_bins, energy, acc, zero, accE);
	gAcc->SetTitle("Compton Acceptance");
	gAcc->SetMarkerStyle(8);
	gAcc->SetMarkerSize(0.5);
	gAcc->SetMarkerColor(kBlue);
	gAcc->SetLineColor(kBlue);
	gAcc->GetXaxis()->SetTitle("Photon Beam Energy [GeV]");
	gAcc->GetXaxis()->SetTitleSize(0.05);
	gAcc->GetXaxis()->SetTitleOffset(0.8);
	gAcc->GetYaxis()->SetTitle("N_{rec} / N_{gen}");
	gAcc->GetYaxis()->SetTitleSize(0.05);
	gAcc->GetYaxis()->SetTitleOffset(0.8);
	
	f_acc_1 = new TF1("f_acc_1", "pol5", 5.0, 11.3);
	f_acc_1->SetLineColor(kBlue+1);
	gAcc->Fit(f_acc_1, "R0");
	f_acc_1->SetRange(5.0, 12.0);
	
	if(DRAW_ACCEPTANCE) {
		
		canvas_acc = new TCanvas("canvas_acc", "canvas_acc", 1000, 600);
		canvas_acc->SetTickx(); canvas_acc->SetTicky();
		canvas_acc->cd();
		
		gAcc->Draw("AP");
		f_acc_1->Draw("same");
		
		canvas_acc->Update();
	}
	
	leg->AddEntry(gAcc, name1, "PE");
	
	return;
}


void get_acc_2()
{
	vector<double> tagh_enVec, tagh_accVec, tagh_accEVec;
	vector<double> tagm_enVec, tagm_accVec, tagm_accEVec;
	
	for(int tagh_counter = 1; tagh_counter <= 274; tagh_counter++) {
		
		char fname[256];
		sprintf(fname, "%s/tagh_%03d.root", mc_dir_2, tagh_counter);
		
		if(gSystem->AccessPathName(fname)) continue;
		
		//cout << "Processing TAGH Counter " << tagh_counter << endl;
		
		double loc_eb = tagh_en2[tagh_counter-1];
		
		TFile *fSim = new TFile(fname, "READ");
		
		TH1F *h_vertex = (TH1F*)fSim->Get(
			"compton_analysis_TOF/vertex_accepted")->Clone(
			Form("h_vertex_tagh_%d",tagh_counter));
		
		double n_gen = h_vertex->Integral();
		
		//-----------------------------------------------------//
		
		TH2F *h2_tagh = new TH2F(Form("h2_tagh_tagh_sim_%d",tagh_counter), "DeltaK", 
			274, 0.5, 274.5, 4000, -8.0, 8.0);
		TH2F *h2_tagm = new TH2F(Form("h2_tagh_tagm_sim_%d",tagh_counter), "DeltaK", 
			102, 0.5, 102.5, 4000, -8.0, 8.0);
		
		h2_tagh->Add((TH2F*)fSim->Get(Form("%s",hname_tagh)));
		h2_tagm->Add((TH2F*)fSim->Get(Form("%s",hname_tagm)));
		
		h2_tagh->SetDirectory(0);
		h2_tagm->SetDirectory(0);
		
		fSim->Close();
		
		TH1F *h1 = (TH1F*)h2_tagh->ProjectionY(Form("h1_sim_tagh_tagh_%d",tagh_counter));
		h1->Add((TH1F*)h2_tagm->ProjectionY(Form("h1_sim_tagh_tagm_%d",tagh_counter)));
		if(h1->Integral() < 1.e1) continue;
		
		h1->Rebin(rebins);
		h1->GetXaxis()->SetTitle("E_{Comp} - E_{#gamma} [GeV]");
		h1->GetXaxis()->SetTitleOffset( 1.1 );
		h1->GetXaxis()->SetRangeUser(-1.5, 1.5);
		h1->GetYaxis()->SetTitle( Form("counts / %d MeV", n_mev));
		h1->GetYaxis()->SetTitleOffset(1.3);
		h1->SetLineColor(kBlack);
		h1->SetLineWidth(2);
		h1->GetXaxis()->SetRangeUser(-3., 2.);
		
		double yield;
		
		if(CALC_ACC_FROM_FIT) {
			
			TF1 *f_gaus = new TF1(Form("f_gaus_tagh_acc_%d",tagh_counter), 
				"gaus", -2., 2.);
			
			f_gaus->SetParameters(h1->GetMaximum(), 
				h1->GetBinCenter(h1->GetMaximumBin()), 0.2);
			f_gaus->SetParLimits(2, 0., 1.);
			
			h1->Fit(f_gaus, "R0QL");
			f_gaus->SetRange(f_gaus->GetParameter(1)-0.2, f_gaus->GetParameter(1)+0.2);
			h1->Fit(f_gaus, "R0QL");
			
			TF1 *f_fit = new TF1(Form("f_fit_tagh_acc_%d",tagh_counter), 
				double_gaus_fit, -2., 2., 10);	
			
			f_fit->SetParName(0, "A");
			f_fit->SetParName(1, "z");
			f_fit->SetParName(2, "#mu_{1}");
			f_fit->SetParName(3, "#mu_{1}-#mu_{2}");
			f_fit->SetParName(4, "#sigma_{1}");
			f_fit->SetParName(5, "#sigma_{2}");
			
			f_fit->SetParLimits(1, 0.0, 0.5);
			f_fit->SetParLimits(4, 0.0, 1.0);
			f_fit->SetParLimits(5, 0.0, 2.0);
			
			f_fit->SetParameters(f_gaus->GetParameter(0), 0., f_gaus->GetParameter(1), 
				0., f_gaus->GetParameter(2), 2.*f_gaus->GetParameter(2));
			
			f_fit->FixParameter(6, 0.0);
			f_fit->FixParameter(7, 0.0);
			f_fit->FixParameter(8, 0.0);
			f_fit->FixParameter(9, 0.0);
			
			f_fit->SetRange(-2.0, 2.0);
			h1->Fit(f_fit, "R0QL");
			f_fit->SetRange(-2., 1.5);
			h1->Fit(f_fit, "R0QL");
			
			yield  = f_fit->GetParameter(0) * sqrt(2.*TMath::Pi()) 
				* ((1.-f_fit->GetParameter(1))*f_fit->GetParameter(4) 
				+ f_fit->GetParameter(1)*f_fit->GetParameter(5));
			yield /= bin_size;
		} else {
			
			yield = h1->Integral();
		}
		
		double loc_acc  = yield / n_gen;
		double loc_accE = sqrt(n_gen*loc_acc*(1.-loc_acc)) / n_gen;
		
		tagh_acc2[tagh_counter-1]  = loc_acc;
		tagh_acc2E[tagh_counter-1] = loc_accE;
		
		tagh_enVec.push_back(loc_eb);
		tagh_accVec.push_back(loc_acc);
		tagh_accEVec.push_back(loc_accE);
	}
	
	for(int tagm_counter = 1; tagm_counter <= 102; tagm_counter++) {
		
		if(tagm_counter>=71 && tagm_counter<=73) continue;
		
		char fname[256];
		sprintf(fname, "%s/tagm_%03d.root", mc_dir_2, tagm_counter);
		
		if(gSystem->AccessPathName(fname)) continue;
		
		//cout << "Processing TAGM Counter " << tagm_counter << endl;
		
		double loc_eb = tagm_en2[tagm_counter-1];
		
		TFile *fSim = new TFile(fname, "READ");
		
		TH1F *h_vertex = (TH1F*)fSim->Get(
			"compton_analysis_TOF/vertex_accepted")->Clone(
			Form("h_vertex_tagm_%d",tagm_counter));
		
		double n_gen = h_vertex->Integral();
		
		//-----------------------------------------------------//
		
		TH2F *h2_tagh = new TH2F(Form("h2_tagm_tagh_sim_%d",tagm_counter), "DeltaK", 
			274, 0.5, 274.5, 4000, -8.0, 8.0);
		TH2F *h2_tagm = new TH2F(Form("h2_tagm_tagm_sim_%d",tagm_counter), "DeltaK", 
			102, 0.5, 102.5, 4000, -8.0, 8.0);
		
		h2_tagh->Add((TH2F*)fSim->Get(Form("%s",hname_tagh)));
		h2_tagm->Add((TH2F*)fSim->Get(Form("%s",hname_tagm)));
		
		h2_tagh->SetDirectory(0);
		h2_tagm->SetDirectory(0);
		
		fSim->Close();
		
		TH1F *h1 = (TH1F*)h2_tagh->ProjectionY(Form("h1_sim_tagm_tagh_%d",tagm_counter));
		h1->Add((TH1F*)h2_tagm->ProjectionY(Form("h1_sim_tagm_tagm_%d",tagm_counter)));
		if(h1->Integral() < 1.e1) continue;
		
		h1->Rebin(rebins);
		h1->GetXaxis()->SetTitle("E_{Comp} - E_{#gamma} [GeV]");
		h1->GetXaxis()->SetTitleOffset(1.1);
		h1->GetXaxis()->SetRangeUser(-1.5, 1.5);
		h1->GetYaxis()->SetTitle(Form("counts / %d MeV", n_mev));
		h1->GetYaxis()->SetTitleOffset(1.3);
		h1->SetLineColor(kBlack);
		h1->SetLineWidth(2);
		h1->GetXaxis()->SetRangeUser(-3., 2.);
		
		double yield;
		
		if(CALC_ACC_FROM_FIT) {
			
			TF1 *f_gaus = new TF1(Form("f_gaus_tagm_acc_%d",tagm_counter), 
				"gaus", -2., 2.);
			
			f_gaus->SetParameters(h1->GetMaximum(), 
				h1->GetBinCenter(h1->GetMaximumBin()), 0.2);
			f_gaus->SetParLimits(2, 0., 1.);
			
			h1->Fit(f_gaus, "R0QL");
			f_gaus->SetRange(f_gaus->GetParameter(1)-0.2, f_gaus->GetParameter(1)+0.2);
			h1->Fit(f_gaus, "R0QL");
			
			TF1 *f_fit = new TF1(Form("f_fit_tagm_acc_%d",tagm_counter), 
				double_gaus_fit, -2., 2., 10);	
			
			f_fit->SetParName(0, "A");
			f_fit->SetParName(1, "z");
			f_fit->SetParName(2, "#mu_{1}");
			f_fit->SetParName(3, "#mu_{1}-#mu_{2}");
			f_fit->SetParName(4, "#sigma_{1}");
			f_fit->SetParName(5, "#sigma_{2}");
			
			f_fit->SetParLimits(1, 0.0, 0.5);
			f_fit->SetParLimits(4, 0.0, 1.0);
			f_fit->SetParLimits(5, 0.0, 2.0);
			
			f_fit->SetParameters(f_gaus->GetParameter(0), 0., f_gaus->GetParameter(1), 
				0., f_gaus->GetParameter(2), 2.*f_gaus->GetParameter(2));
			
			f_fit->FixParameter(6, 0.0);
			f_fit->FixParameter(7, 0.0);
			f_fit->FixParameter(8, 0.0);
			f_fit->FixParameter(9, 0.0);
			
			f_fit->SetRange(-2.0, 2.0);
			h1->Fit(f_fit, "R0QL");
			f_fit->SetRange(-2., 1.5);
			h1->Fit(f_fit, "R0QL");
			
			yield  = f_fit->GetParameter(0) * sqrt(2.*TMath::Pi()) 
				* ((1.-f_fit->GetParameter(1))*f_fit->GetParameter(4) 
				+ f_fit->GetParameter(1)*f_fit->GetParameter(5));
			yield /= bin_size;
		} else {
			
			yield = h1->Integral();
		}
		
		double loc_acc  = yield / n_gen;
		double loc_accE = sqrt(n_gen*loc_acc*(1.-loc_acc)) / n_gen;
		
		tagm_acc2[tagm_counter-1]  = loc_acc;
		tagm_acc2E[tagm_counter-1] = loc_accE;
		
		tagm_enVec.push_back(loc_eb);
		tagm_accVec.push_back(loc_acc);
		tagm_accEVec.push_back(loc_accE);
	}
	
	
	int n_bins1 = (int)tagh_enVec.size();
	double *energy1 = new double[n_bins1];
	double *zero1   = new double[n_bins1];
	double *acc1    = new double[n_bins1];
	double *acc1E   = new double[n_bins1];
	for( int i=0; i<n_bins1; i++ ) {
		energy1[i] = tagh_enVec[i];
		zero1[i]   = 0.;
		acc1[i]    = tagh_accVec[i];
		acc1E[i]   = tagh_accEVec[i];
	}
	
	int n_bins2 = (int)tagm_enVec.size();
	double *energy2 = new double[n_bins2];
	double *zero2   = new double[n_bins2];
	double *acc2    = new double[n_bins2];
	double *acc2E   = new double[n_bins2];
	for( int i=0; i<n_bins2; i++ ) {
		energy2[i] = tagm_enVec[i];
		zero2[i]   = 0.;
		acc2[i]    = tagm_accVec[i];
		acc2E[i]   = tagm_accEVec[i];
	}
	
	int n_bins = n_bins1 + n_bins2;
	double *energy  = new double[n_bins];
	double *zero    = new double[n_bins];
	double *acc     = new double[n_bins];
	double *accE    = new double[n_bins];
	for( int i=0; i<n_bins; i++ ) {
		if( i<n_bins1 ) {
			energy[i] = tagh_enVec[i];
			zero[i]   = 0.;
			acc[i]    = tagh_accVec[i];
			accE[i]   = tagh_accEVec[i];
		} else {
			energy[i] = tagm_enVec[i-n_bins1];
			zero[i]   = 0.;
			acc[i]    = tagm_accVec[i-n_bins1];
			accE[i]   = tagm_accEVec[i-n_bins1];
		}
	}
	
	TGraphErrors *gAcc = new TGraphErrors(n_bins, energy, acc, zero, accE);
	gAcc->SetTitle("Compton Acceptance");
	gAcc->SetMarkerStyle(8);
	gAcc->SetMarkerSize(0.5);
	gAcc->SetMarkerColor(kRed);
	gAcc->SetLineColor(kRed);
	gAcc->GetXaxis()->SetTitle("Photon Beam Energy [GeV]");
	gAcc->GetXaxis()->SetTitleSize(0.05);
	gAcc->GetXaxis()->SetTitleOffset(0.8);
	gAcc->GetYaxis()->SetTitle("N_{rec} / N_{gen}");
	gAcc->GetYaxis()->SetTitleSize(0.05);
	gAcc->GetYaxis()->SetTitleOffset(0.8);
	
	f_acc_2 = new TF1("f_acc_2", "pol5", 5.0, 11.3);
	f_acc_2->SetLineColor(kRed+1);
	gAcc->Fit(f_acc_2, "R0");
	f_acc_2->SetRange(5.0, 12.0);
	
	if(DRAW_ACCEPTANCE) {
		
		canvas_acc->cd();
		
		gAcc->Draw("P same");
		f_acc_2->Draw("same");
		
		canvas_acc->Update();
	}
	
	leg->AddEntry(gAcc, name2, "PE");
	
	return;
}
