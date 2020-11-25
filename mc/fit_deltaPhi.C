
char   tagh_flux_fname[256],  empty_target_tagh_flux_fname[256];
char   tagm_flux_fname[256],  empty_target_tagm_flux_fname[256];
char tagh_xscale_fname[256],             tagm_xscale_fname[256];
char        hname_tagh[256],                    hname_tagm[256];
char            mc_dir[256];

double endpoint_energy, endpoint_energy_calib;

bool CS_FROM_FIT;
bool FIX_BACKGROUND;

double          tagh_en[274],          tagm_en[102];
double        tagh_flux[274],        tagm_flux[102];
double       tagh_fluxE[274],       tagm_fluxE[102];
double  tagh_flux_empty[274],  tagm_flux_empty[102];
double tagh_fluxE_empty[274], tagm_fluxE_empty[102];


//----------   Function Declarations   ----------//

void get_flux();
void get_counter_energies();

Double_t bkgd_fit( Double_t *x, Double_t *par );
Double_t crys_ball_fit( Double_t *x, Double_t *par );
Double_t line_shape_fit( Double_t *x, Double_t *par );

Double_t fit_pol3( Double_t *x, Double_t *par );

//-----------------------------------------------//


double bin_size = 360. / 3600.;
vector<double> comp_sim_hist;

double tagh_mu[274], tagh_sig[274];
double tagm_mu[102], tagm_sig[102];



void fit_deltaPhi()
{
	
	
	const char pathName[256] = "/work/halld/home/andrsmit/primex_compton_analysis";
	
	// photon flux file names:
	
	sprintf( tagh_flux_fname,              "%s/photon_flux/Be_tagh_flux.txt",       pathName );
	sprintf( tagm_flux_fname,              "%s/photon_flux/Be_tagm_flux.txt",       pathName );
	sprintf( empty_target_tagh_flux_fname, "%s/photon_flux/Be_empty_tagh_flux.txt", pathName );
	sprintf( empty_target_tagm_flux_fname, "%s/photon_flux/Be_empty_tagm_flux.txt", pathName );
	
	// file containing the xscales for the tagh and tagm counters:
	
	sprintf( tagh_xscale_fname, "%s/photon_flux/primex_tagh.txt", pathName );
	sprintf( tagm_xscale_fname, "%s/photon_flux/primex_tagm.txt", pathName );
	
	endpoint_energy       = 11.6061; // Be 200 nA data
	//endpoint_energy       = 11.1671; // He  50 nA data
	//endpoint_energy       = 11.1664; // He 100 nA data
	endpoint_energy_calib = 11.6061;
	
	// Name of histograms for fitting yield:
	
	sprintf( hname_tagh, "DeltaPhi/deltaPhi_tagh_e" );
	sprintf( hname_tagm, "DeltaPhi/deltaPhi_tagm_e" );
	
	// Directory where mc rootFiles are stored:
	
	sprintf( mc_dir, "%s/tag_sim/recRootFiles_v1", pathName );
	
	
	
	//-------------      Initialize     -------------//
	
	
	gStyle->SetOptStat(0); gStyle->SetOptFit(1);
	
	const int rebins  =  2;
	bin_size         *=  (double)rebins;
	double n_deg      =  0.16 * (double)rebins;
	
	const bool DRAW_FITS_TAGH = false;
	const bool DRAW_FITS_TAGM = false;
	
	
	//------------------------------------------------//
	
	
	
	get_flux();
	get_counter_energies();
	
	
	
	TCanvas *canvas = new TCanvas( "canvas", "canvas", 800, 800 );
	canvas->SetTickx(); canvas->SetTicky();
	canvas->SetGrid();
	//canvas->SetLogy();
	
	
	for( int i=0; i<274; i++ ) { tagh_mu[i] = 0.; tagh_sig[i] = 0.; }
	for( int i=0; i<102; i++ ) { tagm_mu[i] = 0.; tagm_sig[i] = 0.; }
	
	
	
	
	//----------   Fit TAGH Counters   --------//
	
	vector<double>  en1Vec,  mu1Vec,  sig1Vec,  alpha1Vec,  n1Vec, chi21Vec;
	vector<double> en1EVec, mu1EVec, sig1EVec, alpha1EVec, n1EVec;
	
	int loc_counter = 0;
	
	for( int tagh_counter = 1; tagh_counter <= 274; tagh_counter++ ) 
	{
		
		double eb = tagh_en[tagh_counter-1];
		
		if( gSystem->AccessPathName(Form("%s/tagh_%03d.root",mc_dir,tagh_counter)) ) {
			cout << "Skipping TAGH counter " << tagh_counter << endl;
			continue;
		}
		
		TFile *fIn = new TFile( Form("%s/tagh_%03d.root", mc_dir, tagh_counter) );
		
		//TH2F  *h2  = (TH2F*)fIn->Get( hname_tagh );
		TH2F  *h2  = new TH2F( Form("h2_tagh_%d",tagh_counter), "#Delta#phi vs. TAGH Counter", 
			280, -0.5, 279.5, 1000, 100.0, 260.0 );
		
		for( int ir=1; ir<=8; ir++ ) {
			h2->Add( (TH2F*)fIn->Get(Form("%s_%d",hname_tagh,ir)) );
		}
		
		
		TH1F  *h1  = (TH1F*)h2->ProjectionY( Form("h1_tagh_%d",tagh_counter) );
		
		if( h1->GetEntries() < 1.e1 ) {
			cout << "Skipping TAGH counter " << tagh_counter << endl;
			fIn->Close();
			continue;
		}
		
		h1->Rebin(rebins);
		h1->SetTitle( Form("TAGH Counter %d (E_{#gamma} = %.3f)", tagh_counter, eb) );
		h1->GetXaxis()->SetTitle( "|#phi_{1} - #phi_{2}| [deg.]" );
		h1->GetXaxis()->SetTitleOffset( 1.1 );
		h1->GetYaxis()->SetTitle( Form("counts / %d MeV", n_mev) );
		h1->GetYaxis()->SetTitleOffset( 1.3 );
		h1->SetLineColor( kBlack );
		h1->SetLineWidth( 2 );
		
		
		canvas->cd();
		h1->Draw();
		
		
		//----------     Gaussian Fit     ----------//
		
		TF1 *f_fit = new TF1( Form("f_fit_tagh_%d",tagh_counter), "gaus", 170., 190. );
		
		f_fit->SetParameters( h1->GetMaximum(), 180.0, 6.5 );
		f_fit->SetLineColor( kRed );
		
		h1->Fit( Form("f_fit_tagh_%d",tagh_counter), "R0QL" );
		
		f_fit->SetRange( f_fit->GetParameter(1) - 2.*f_fit->GetParameter(2), 
			f_fit->GetParameter(1) + 2.*f_fit->GetParameter(2) );
		
		TFitResultPtr result = h1->Fit( f_fit, "SR0QL" );
		double chi2          = result->Chi2() / result->Ndf();
		
		  en1Vec.push_back(eb);
		 en1EVec.push_back(0.);
		  mu1Vec.push_back(f_fit->GetParameter(1));
		 mu1EVec.push_back( f_fit->GetParError(1));
		 sig1Vec.push_back(f_fit->GetParameter(2));
		sig1EVec.push_back( f_fit->GetParError(2));
		chi21Vec.push_back( chi2 );
		
		
		// Draw histogram and fit function:
		
		canvas->cd();
		h1->GetXaxis()->SetRangeUser( 120.0, 240.0 );
		h1->SetMinimum(0.);
		h1->Draw( "hist" );
		
		f_fit->SetRange( 120.0, 240. );
		f_fit->Draw( "same" );
		
		
		// push fit parameters to vectors:
		
		en1Vec.push_back( eb );
		en1EVec.push_back( 0.0 );
		
		mu1Vec.push_back( f_fit->GetParameter(1) );
		mu1EVec.push_back( f_fit->GetParError(1) );
		
		sig1Vec.push_back( f_fit->GetParameter(2) );
		sig1EVec.push_back( f_fit->GetParError(2) );
		
		chi21Vec.push_back( chi2 );
		
		if( DRAW_FITS_TAGH ) {
		
		canvas->Update();
		if( loc_counter==0 ) canvas->Print( "deltaPhi_tagh.pdf(", "pdf" );
		else                  canvas->Print( "deltaPhi_tagh.pdf",  "pdf" );
		
		}
		
		loc_counter++;
		
		tagh_mu[tagh_counter-1]  = f_fit->GetParameter(1);
		tagh_sig[tagh_counter-1] = f_fit->GetParameter(2);
		
		fIn->Close();
	}
	
	if( DRAW_FITS_TAGH ) canvas->Print( "deltaPhi_tagh.pdf)", "pdf" );
	
	
	
	//----------   Fit TAGM Counters   --------//
	
	vector<double>  en2Vec,  mu2Vec,  sig2Vec,  alpha2Vec,  n2Vec, chi22Vec;
	vector<double> en2EVec, mu2EVec, sig2EVec, alpha2EVec, n2EVec;
	
	loc_counter = 0;
	
	for( int tagm_counter = 1; tagm_counter <= 102; tagm_counter++ ) 
	{
		
		double eb = tagm_en[tagm_counter-1];
		
		if( gSystem->AccessPathName(Form("%s/tagm_%03d.root",mc_dir,tagm_counter)) ) {
			cout << "Skipping TAGM counter " << tagm_counter << endl;
			continue;
		}
		
		TFile *fIn = new TFile( Form("%s/tagm_%03d.root", mc_dir, tagm_counter) );
		
		//TH2F  *h2  = (TH2F*)fIn->Get( hname_tagm );
		TH2F  *h2  = new TH2F( Form("h2_tagh_%d",tagm_counter), "#Delta#phi vs. TAGM Counter", 
			110, -0.5, 109.5, 1000, 100.0, 260.0 );
		
		for( int ir=1; ir<=8; ir++ ) {
			h2->Add( (TH2F*)fIn->Get(Form("%s_%d",hname_tagm,ir)) );
		}
		
		TH1F  *h1  = (TH1F*)h2->ProjectionY( Form("h1_tagm_%d",tagm_counter) );
		
		if( h1->GetEntries() < 1.e1 ) {
			cout << "Skipping TAGM counter " << tagm_counter << endl;
			fIn->Close();
			continue;
		}
		
		h1->Rebin(rebins);
		h1->SetTitle( Form("TAGM Counter %d (E_{#gamma} = %.3f)", tagm_counter, eb) );
		h1->GetXaxis()->SetTitle( "|#phi_{1} - #phi_{2}| [deg.]" );
		h1->GetXaxis()->SetTitleOffset( 1.1 );
		h1->GetYaxis()->SetTitle( Form("counts / %d MeV", n_mev) );
		h1->GetYaxis()->SetTitleOffset( 1.3 );
		h1->SetLineColor( kBlack );
		h1->SetLineWidth( 2 );
		
		
		canvas->cd();
		h1->Draw();
		
		
		//----------     Gaussian Fit     ----------//
		
		TF1 *f_fit = new TF1( Form("f_fit_tagm_%d",tagm_counter), "gaus", 170., 190. );
		
		f_fit->SetParameters( h1->GetMaximum(), 180.0, 6.5 );
		f_fit->SetLineColor( kRed );
		
		h1->Fit( Form("f_fit_tagm_%d",tagm_counter), "R0QL" );
		
		f_fit->SetRange( f_fit->GetParameter(1) - 2.*f_fit->GetParameter(2), 
			f_fit->GetParameter(1) + 2.*f_fit->GetParameter(2) );
		
		TFitResultPtr result = h1->Fit( f_fit, "SR0QL" );
		double chi2          = result->Chi2() / result->Ndf();
		
		  en1Vec.push_back(eb);
		 en1EVec.push_back(0.);
		  mu1Vec.push_back(f_fit->GetParameter(1));
		 mu1EVec.push_back( f_fit->GetParError(1));
		 sig1Vec.push_back(f_fit->GetParameter(2));
		sig1EVec.push_back( f_fit->GetParError(2));
		chi21Vec.push_back( chi2 );
		
		
		// Draw histogram and fit function:
		
		canvas->cd();
		h1->GetXaxis()->SetRangeUser( 120.0, 240.0 );
		h1->SetMinimum(0.);
		h1->Draw( "hist" );
		
		f_fit->SetRange( 120.0, 240.0 );
		f_fit->Draw( "same" );
		
		
		// push fit parameters to vectors:
		
		en2Vec.push_back( eb );
		en2EVec.push_back( 0.0 );
		
		mu2Vec.push_back( f_fit->GetParameter(1) );
		mu2EVec.push_back( f_fit->GetParError(1) );
		
		sig2Vec.push_back( f_fit->GetParameter(2) );
		sig2EVec.push_back( f_fit->GetParError(2) );
		
		chi22Vec.push_back( chi2 );
		
		if( DRAW_FITS_TAGM ) {
		
		canvas->Update();
		if( loc_counter==0 ) canvas->Print( "deltaPhi_tagm.pdf(", "pdf" );
		else                  canvas->Print( "deltaPhi_tagm.pdf",  "pdf" );
		
		}
		
		loc_counter++;
		
		tagm_mu[tagm_counter-1]  = f_fit->GetParameter(1);
		tagm_sig[tagm_counter-1] = f_fit->GetParameter(2);
		
		fIn->Close();
	}
	
	if( DRAW_FITS_TAGM ) canvas->Print( "deltaPhi_tagm.pdf)", "pdf" );
	
	
	
	//----------   Plot fit results   --------//
	
	int n_bins1 = (int)en1Vec.size();
	
	double *energy1  = new double[n_bins1];
	double *energy1E = new double[n_bins1];
	
	double *mu1      = new double[n_bins1];
	double *mu1E     = new double[n_bins1];
	
	double *sig1     = new double[n_bins1];
	double *sig1E    = new double[n_bins1];
	
	double *chi21    = new double[n_bins1];
	
	for( int i = 0; i < n_bins1; i++ ) {
		
		energy1[i]  = en1Vec[i];
		energy1E[i] = en1EVec[i];
		
		mu1[i]      = mu1Vec[i];
		mu1E[i]     = mu1EVec[i];
		
		sig1[i]     = sig1Vec[i];
		sig1E[i]    = sig1EVec[i];
		
		chi21[i]    = chi21Vec[i];
		
	}
	
	
	int n_bins2 = (int)en2Vec.size();
	
	double *energy2  = new double[n_bins2];
	double *energy2E = new double[n_bins2];
	
	double *mu2      = new double[n_bins2];
	double *mu2E     = new double[n_bins2];
	
	double *sig2     = new double[n_bins2];
	double *sig2E    = new double[n_bins2];
	
	double *chi22    = new double[n_bins2];
	
	for( int i = 0; i < n_bins2; i++ ) {
		
		energy2[i]  = en2Vec[i];
		energy2E[i] = en2EVec[i];
		
		mu2[i]      = mu2Vec[i];
		mu2E[i]     = mu2EVec[i];
		
		sig2[i]     = sig2Vec[i];
		sig2E[i]    = sig2EVec[i];
		
		chi22[i]    = chi22Vec[i];
		
	}
	
	
	int n_bins = n_bins1 + n_bins2;;
	
	double *energy_full  = new double[n_bins];
	double *energyE_full = new double[n_bins];
	
	double *mu_full      = new double[n_bins];
	double *muE_full     = new double[n_bins];
	
	double *sig_full     = new double[n_bins];
	double *sigE_full    = new double[n_bins];
	
	for( int i = 0; i < n_bins; i++ ) {
		
		if( i < n_bins1 ) {
			
			energy_full[i]  = en1Vec[i];
			energyE_full[i] = en1EVec[i];
			
			mu_full[i]      = mu1Vec[i];
			muE_full[i]     = mu1EVec[i];
			
			sig_full[i]     = sig1Vec[i];
			sigE_full[i]    = sig1EVec[i];
			
		} else {
			
			energy_full[i]  = en2Vec[i-n_bins1];
			energyE_full[i] = en2EVec[i-n_bins1];
			
			mu_full[i]      = mu2Vec[i-n_bins1];
			muE_full[i]     = mu2EVec[i-n_bins1];
			
			sig_full[i]     = sig2Vec[i-n_bins1];
			sigE_full[i]    = sig2EVec[i-n_bins1];
			
		}
		
	}
	
	
	TGraphErrors *gMu1 = new TGraphErrors( n_bins1, energy1, mu1, energy1E, mu1E );
	gMu1->SetTitle( "#mu_{#DeltaE} vs. E_{#gamma}" );
	gMu1->GetXaxis()->SetTitle( "E_{#gamma} [GeV]" );
	gMu1->GetYaxis()->SetTitle( "#mu from Crystal Ball Fit [GeV]" );
	gMu1->GetYaxis()->SetTitleOffset( 1.3 );
	gMu1->SetMarkerStyle( 8 );
	gMu1->SetMarkerColor( kBlue );
	TGraphErrors *gMu2 = new TGraphErrors( n_bins2, energy2, mu2, energy2E, mu2E );
	gMu2->SetTitle( "#mu_{#DeltaE} vs. E_{#gamma}" );
	gMu2->GetXaxis()->SetTitle( "E_{#gamma} [GeV]" );
	gMu2->GetYaxis()->SetTitle( "#mu from Crystal Ball Fit [GeV]" );
	gMu2->GetYaxis()->SetTitleOffset( 1.3 );
	gMu2->SetMarkerStyle( 8 );
	gMu2->SetMarkerColor( kGreen );
	TGraphErrors *gMu = new TGraphErrors( n_bins, energy_full, mu_full, energyE_full, muE_full );
	gMu->SetTitle( "#mu_{#DeltaE} vs. E_{#gamma}" );
	gMu->GetXaxis()->SetTitle( "E_{#gamma} [GeV]" );
	gMu->GetYaxis()->SetTitle( "#mu from Crystal Ball Fit [GeV]" );
	gMu->GetYaxis()->SetTitleOffset( 1.3 );
	gMu->SetMarkerStyle( 8 );
	gMu->SetMarkerColor( kBlue );
	
	TGraphErrors *gSig1 = new TGraphErrors( n_bins1, energy1, sig1, energy1E, sig1E );
	gSig1->SetTitle( "#sigma_{#DeltaE} / E vs. E_{#gamma}" );
	gSig1->GetXaxis()->SetTitle( "E_{#gamma} [GeV]" );
	gSig1->GetYaxis()->SetTitle( "#sigma / E from Crystal Ball Fit [GeV]" );
	gSig1->GetYaxis()->SetTitleOffset( 1.3 );
	gSig1->SetMarkerStyle( 8 );
	gSig1->SetMarkerColor( kBlue );
	TGraphErrors *gSig2 = new TGraphErrors( n_bins2, energy2, sig2, energy2E, sig2E );
	gSig2->SetTitle( "#sigma_{#DeltaE} vs. E_{#gamma}" );
	gSig2->GetXaxis()->SetTitle( "E_{#gamma} [GeV]" );
	gSig2->GetYaxis()->SetTitle( "#sigma from Crystal Ball Fit [GeV]" );
	gSig2->GetYaxis()->SetTitleOffset( 1.3 );
	gSig2->SetMarkerStyle( 8 );
	gSig2->SetMarkerColor( kGreen );
	TGraphErrors *gSig = new TGraphErrors( n_bins, energy_full, sig_full, energyE_full, sigE_full );
	gSig->SetTitle( "#sigma_{#DeltaE} / E vs. E_{#gamma}" );
	gSig->GetXaxis()->SetTitle( "E_{#gamma} [GeV]" );
	gSig->GetYaxis()->SetTitle( "#sigma / E from Crystal Ball Fit [GeV]" );
	gSig->GetYaxis()->SetTitleOffset( 1.3 );
	gSig->SetMarkerStyle( 8 );
	gSig->SetMarkerColor( kBlack );
	
	TGraph *gChi21 = new TGraph( n_bins1, energy1, chi21 );
	gChi21->SetTitle( "#chi^{2} vs. E_{#gamma}" );
	gChi21->GetXaxis()->SetTitle( "E_{#gamma} [GeV]" );
	gChi21->GetYaxis()->SetTitle( "#chi^{2} from Crystal Ball Fit" );
	gChi21->GetYaxis()->SetTitleOffset( 1.3 );
	gChi21->SetMarkerStyle( 8 );
	gChi21->SetMarkerColor( kBlue );
	TGraph *gChi22 = new TGraph( n_bins2, energy2, chi22 );
	gChi22->SetTitle( "#chi^{2} vs. E_{#gamma}" );
	gChi22->GetXaxis()->SetTitle( "E_{#gamma} [GeV]" );
	gChi22->GetYaxis()->SetTitle( "#chi^{2} from Crystal Ball Fit" );
	gChi22->GetYaxis()->SetTitleOffset( 1.3 );
	gChi22->SetMarkerStyle( 8 );
	gChi22->SetMarkerColor( kGreen );
	
	
	TCanvas *cMu = new TCanvas( "cMu", "cMu", 800, 400 );
	cMu->SetTickx(); cMu->SetTicky();
	gMu->Draw( "AP" );
	gMu1->Draw( "P same" );
	gMu2->Draw( "P same" );
	
	TF1 *fMu = new TF1( "fMu", fit_pol3, 6.0, 11.4, 4 );
	gMu->Fit( "fMu", "R0" );
	fMu->Draw("same");
	
	
	
	TCanvas *cSig = new TCanvas( "cSig", "cSig", 800, 400 );
	cSig->SetTickx(); cSig->SetTicky();
	gSig->Draw( "AP" );
	gSig1->Draw( "P same" );
	gSig2->Draw( "P same" );
	
	TF1 *fSig = new TF1( "fSig", fit_pol3, 6.0, 11.4, 4 );
	gSig->Fit( "fSig", "R0" );
	fSig->Draw("same");
	
	TCanvas *cChi2 = new TCanvas( "cChi2", "cChi2", 800, 400 );
	cChi2->SetTickx(); cChi2->SetTicky();
	gChi21->Draw( "AP" );
	gChi22->Draw( "P same" );
	
	
	
	
	//----------   Write mean and widths for each counter to file   ----------//
	
	/*
	ofstream outfile1( "tagh_esig.dat" );
	for( int i=0; i<274; i++ ) outfile1 << "			" << tagh_sig[i] << "\n";
	outfile1.close();
	
	ofstream outfile2( "tagh_emu.dat" );
	for( int i=0; i<274; i++ ) outfile2 << "			" << tagh_mu[i]  << "\n";
	outfile2.close();
	
	ofstream outfile3( "tagm_esig.dat" );
	for( int i=0; i<102; i++ ) outfile3 << "			" << tagm_sig[i] << "\n";
	outfile3.close();
	
	ofstream outfile4( "tagm_emu.dat" );
	for( int i=0; i<102; i++ ) outfile4 << "			" << tagm_mu[i]  << "\n";
	outfile4.close();
	*/
	
	/*
	char buf[256];
	
	ofstream outfile1( "deltaE_tagh.dat" );
	for( int i=0; i<274; i++ ) {
		sprintf( buf, "%03d   %f   %f", i+1, tagh_mu[i], tagh_sig[i] );
		outfile1 << buf << "\n";
	}
	outfile1.close();
	
	ofstream outfile2( "deltaE_tagm.dat" );
	for( int i=0; i<102; i++ ) {
		sprintf( buf, "%03d   %f   %f", i+1, tagm_mu[i], tagm_sig[i] );
		outfile2 << buf << "\n";
	}
	outfile2.close();
	*/
	
	
	return;
}




void get_flux() 
{
	
	int a; double b, c;
	
	ifstream inf1( tagh_flux_fname );
	for( int i=0; i<274; i++ ) {
		inf1 >> a >> b >> c;
		tagh_flux[i]  = b;
		tagh_fluxE[i] = c;
	}
	inf1.close();
	
	ifstream inf2( tagm_flux_fname );
	for( int i=0; i<102; i++ ) {
		inf2 >> a >> b >> c;
		tagm_flux[i]  = b;
		tagm_fluxE[i] = c;
	}
	inf2.close();
	
	ifstream inf3( empty_target_tagh_flux_fname );
	for( int i=0; i<274; i++ ) {
		inf3 >> a >> b >> c;
		tagh_flux_empty[i]  = b;
		tagh_fluxE_empty[i] = c;
	}
	inf3.close();
	
	ifstream inf4( empty_target_tagm_flux_fname );
	for( int i=0; i<102; i++ ) {
		inf4 >> a >> b >> c;
		tagm_flux_empty[i]  = b;
		tagm_fluxE_empty[i] = c;
	}
	inf4.close();
	
	
	return;
}




void get_counter_energies() 
{
	
	int a; double b, c;
	
	ifstream inf1( tagh_xscale_fname );
	for( int i=0; i<274; i++ ) {
		inf1 >> a >> b >> c;
		double deltaE  = endpoint_energy - endpoint_energy_calib;
		double emin    = b * endpoint_energy_calib  +  deltaE;
		double emax    = c * endpoint_energy_calib  +  deltaE;
		tagh_en[i] = 0.5 * (emin + emax);
	}
	inf1.close();
	
	ifstream inf2( tagm_xscale_fname );
	for( int i=0; i<102; i++ ) {
		inf2 >> a >> b >> c;
		double deltaE  = endpoint_energy - endpoint_energy_calib;
		double emin    = b * endpoint_energy_calib  +  deltaE;
		double emax    = c * endpoint_energy_calib  +  deltaE;
		tagm_en[i] = 0.5 * (emin + emax);
	}
	inf2.close();
	
	
	return;
}



Double_t line_shape_fit( Double_t *x, Double_t *par ) {
	
	Double_t fitLS, fq, fitVal;
	
	Int_t xs = (int)((x[0]-par[1])/bin_size + 250);
	fitLS    = par[0] * comp_sim_hist[xs];
	
	fq       = par[2] + 
		   par[3] * x[0] +
		   par[4] * x[0] * x[0] + 
		   par[5] * x[0] * x[0] * x[0];
	
	fitVal = fitLS + fq;
	
	return fitVal;
}



Double_t fit_pol3( Double_t *x, Double_t *par ) {
	
	Double_t xx = x[0];
	
	
	Double_t f = par[0] + 
		     par[1] * x[0] +
		     par[2] * x[0] * x[0] + 
		     par[3] * x[0] * x[0] * x[0];
	
	return f;
}



Double_t bkgd_fit( Double_t *x, Double_t *par ) {
	
	Double_t xx = x[0];
	
	
	if( xx > -1.0 && xx < 1.0 ) {
		TF1::RejectPoint();
		return 0.;
	}
	
	Double_t f = par[0] + 
		     par[1] * x[0] +
		     par[2] * x[0] * x[0] + 
		     par[3] * x[0] * x[0] * x[0];
	
	return f;
}



Double_t crys_ball_fit( Double_t *x, Double_t *par ) 
{
	
	Double_t xx  =   x[0];
	
	Double_t N   = par[0];
	Double_t mu  = par[1];
	Double_t sig = par[2];
	Double_t a   = par[3];
	Double_t n   = par[4];
	
	Double_t p0  = par[5];
	Double_t p1  = par[6];
	Double_t p2  = par[7];
	Double_t p3  = par[8];
	
	
	Double_t A   = 	pow( n/fabs(a), n ) * exp( -0.5*pow(fabs(a),2.0) );
	Double_t B   =  ( n/fabs(a) ) - fabs(a);
	
	Double_t loc_x = ( xx - mu ) / sig;
	Double_t f;
	
	if( loc_x > -a ) {
		f = N * exp( -0.5*pow(loc_x,2.0) );
	} else {
		f = N * A * pow( B - loc_x, -n );
	}
	
	f += p0 + p1*xx + p2*xx*xx + p3*xx*xx*xx;
	
	
	return f;
}




