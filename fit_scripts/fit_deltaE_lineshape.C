
char        root_fname[256],       empty_target_root_fname[256];
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
Double_t line_shape_fit( Double_t *x, Double_t *par );
Double_t crys_ball_fit( Double_t *x, Double_t *par );

//-----------------------------------------------//


double bin_size = 8. / 2000.;
int rebins, n_mev;

double ne, mb;
TF1 *f_theory;

vector<double> comp_sim_hist;

TCanvas *canvas;
TPad *pad_lin, *pad_log;

TCanvas *canvas2;
TPad *top_pad, *bot_pad;


double tagh_mu[274], tagh_sig[274];
double tagm_mu[102], tagm_sig[102];


void fit_deltaE_lineshape()
{
	
	
	const char pathName[256] = "/work/halld/home/andrsmit/primex_compton_analysis";
	
	char root_fname[256], empty_target_root_fname[256];
	
	
	// file names for the full and empty target root files:
	
	sprintf( root_fname,              "%s/data/rootFiles/Be.root",       pathName );
	sprintf( empty_target_root_fname, "%s/data/rootFiles/Be_empty.root", pathName );
	
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
	
	sprintf( hname_tagh, "DeltaE/deltaE_tagh" );
	sprintf( hname_tagm, "DeltaE/deltaE_tagm" );
	
	// Directory where mc rootFiles are stored:
	
	sprintf( mc_dir, "%s/compton_mc/recRootFiles", pathName );
	
	
	
	//-------------      Initialize     -------------//
	
	
	gStyle->SetOptStat(0); gStyle->SetOptFit(0);
	
	rebins    =  2;
	bin_size *=  (double)rebins;
	n_mev     =  4 * rebins;
	
	int n_bins = (int)(2000 / rebins);
	
	comp_sim_hist.clear();
	for( int ib = 0; ib < n_bins; ib++ ) {
		comp_sim_hist.push_back( 0.0 );
	}
	
	const bool DRAW_FITS_TAGH = false;
	const bool DRAW_FITS_TAGM = false;
	
	
	// Number of electrons in target:
	
	ne = 8.77937e+23; // number of electrons per cm^2
	mb = 1.e-27;
	
	
	//------------------------------------------------//
	
	
	
	get_flux();
	get_counter_energies();
	
	
	
	TFile *fFull  = new TFile( root_fname,               "READ" );
	TFile *fEmpty = new TFile( empty_target_root_fname,  "READ" );
	
	TCanvas *canvas = new TCanvas( "canvas", "canvas", 800, 800 );
	canvas->SetTickx(); canvas->SetTicky();
	canvas->SetGrid();
	canvas->SetLogy();
	
	
	
	TH2F *h2_tagh  = new TH2F(  "h2_tagh", "#DeltaE vs. Counter", 
		274, 0.5, 274.5,  2000, -4.0, 4.0 );
	TH2F *h2_tagm  = new TH2F(  "h2_tagm", "#DeltaE vs. Counter", 
		102, 0.5, 102.5,  2000, -4.0, 4.0 );
	
	TH2F *h2e_tagh = new TH2F( "h2e_tagh", "#DeltaE vs. Counter", 
		274, 0.5, 274.5,  2000, -4.0, 4.0 );
	TH2F *h2e_tagm = new TH2F( "h2e_tagm", "#DeltaE vs. Counter", 
		102, 0.5, 102.5,  2000, -4.0, 4.0 );
	
	h2_tagh->Add( (TH2F*)fFull->Get("DeltaE/deltaE_tagh") );
	h2_tagm->Add( (TH2F*)fFull->Get("DeltaE/deltaE_tagm") );
	h2e_tagh->Add( (TH2F*)fEmpty->Get("DeltaE/deltaE_tagh") );
	h2e_tagm->Add( (TH2F*)fEmpty->Get("DeltaE/deltaE_tagm") );
	
	
	
	for( int i=0; i<274; i++ ) { tagh_mu[i] = 0.; tagh_sig[i] = 0.; }
	for( int i=0; i<102; i++ ) { tagm_mu[i] = 0.; tagm_sig[i] = 0.; }
	
	
	
	//----------   Get XS from event generator   ---------//
	
	
	const char cs_pathName[256] = "/work/halld/home/ijaegle/primex_compton/macro/dat";
	char fname[256];
	
	int n_cs_files = 0;
	vector<double> gen_enVec, gen_csVec;
	
	for( int i=60; i<120; i++ ) {
		sprintf( fname, "%s/sd_%d.dat", cs_pathName, 100*i );
		if( gSystem->AccessPathName(fname) ) {
			continue;
		} else {
			double locE = 0.1 * (double)i;
			double a, b, c, d, e;
			double loc_cs = 0.;
			ifstream locInf( fname );
			locInf >> a >> b >> c >> d >> e;
			locInf >> a >> b >> c >> d >> e;
			loc_cs += c;
			locInf >> a >> b >> c >> d >> e;
			loc_cs += c;
			locInf.close();
			gen_enVec.push_back(0.5*(a+b));
			gen_csVec.push_back(loc_cs);
			n_cs_files++;
		}
	}
	
	
	double *gen_energy = new double[n_cs_files];
	double *gen_cs     = new double[n_cs_files];
	
	for( int i=0; i < n_cs_files; i++ ) {
		
		gen_energy[i] = gen_enVec[i];
		gen_cs[i]     = gen_csVec[i];
	}
	
	TGraph *g_theory = new TGraph( n_cs_files, gen_energy, gen_cs );
	g_theory->SetLineColor(kRed); g_theory->SetLineWidth(2);
	g_theory->Draw("same C");
	
	f_theory  = new TF1( "f_theory", "pol5", 5., 12. );
	g_theory->Fit( "f_theory", "R0Q" );
	f_theory->SetLineColor(kRed);
	f_theory->SetLineStyle(2);
	
	
	//----------------------------------------------------//
	
	
	
	
	
	//----------   Fit TAGH Counters   --------//
	
	vector<double>  en1Vec,  mu1Vec,  sig1Vec,  alpha1Vec,  n1Vec, chi21Vec;
	vector<double> en1EVec, mu1EVec, sig1EVec, alpha1EVec, n1EVec;
	
	int loc_counter = 0;
	
	for( int tagh_counter = 1; tagh_counter <= 274; tagh_counter++ ) 
	{
		
		double eb = tagh_en[tagh_counter-1];
		
		if( tagh_flux[tagh_counter-1] <= 0. || tagh_flux_empty[tagh_counter-1] <= 0. ) {
			cout << "Skipping TAGH counter " << tagh_counter << endl;
			continue;
		}
		if( gSystem->AccessPathName(Form("%s/tagh_%03d.root",mc_dir,tagh_counter)) ) {
			cout << "Skipping TAGH counter " << tagh_counter << endl;
			continue;
		}
		
		double loc_flux_ratio = tagh_flux[tagh_counter-1] / tagh_flux_empty[tagh_counter-1];
		
		TH1F *h1  = (TH1F*)h2_tagh->ProjectionY(  Form("h1_tagh_%d", tagh_counter), 
			tagh_counter, tagh_counter );
		TH1F *h1e = (TH1F*)h2e_tagh->ProjectionY( Form("h1e_tagh_%d",tagh_counter), 
			tagh_counter, tagh_counter );
		
		if( h1->GetEntries() < 1.e1 ) {
			cout << "Skipping TAGH counter " << tagh_counter << endl;
			continue;
		}
		
		h1e->Scale( loc_flux_ratio );
		
		h1->Rebin(rebins);
		h1->SetTitle( Form("TAGH Counter %d (E_{#gamma} = %.3f)", tagh_counter, eb) );
		h1->GetXaxis()->SetTitle( "E_{1} + E_{2} - E_{#gamma} [GeV]" );
		h1->GetXaxis()->SetTitleOffset( 1.1 );
		h1->GetYaxis()->SetTitle( Form("counts / %d MeV", n_mev) );
		h1->GetYaxis()->SetTitleOffset( 1.3 );
		h1->SetLineColor( kBlack );
		h1->SetLineWidth( 2 );
		
		h1e->Rebin(rebins);
		h1e->GetXaxis()->SetTitle( "E_{1} + E_{2} - E_{#gamma} [GeV]" );
		h1e->GetXaxis()->SetTitleOffset( 1.1 );
		h1e->GetYaxis()->SetTitle( Form("counts / %d MeV", n_mev) );
		h1e->GetYaxis()->SetTitleOffset( 1.3 );
		h1e->SetLineColor( kBlue );
		h1e->SetLineWidth( 2 );
		
		h1->Add( h1e,-1.0 );
		
		//---------   Get Simulated histogram   --------//
		
		
		TFile *fSim = new TFile( Form("%s/tagh_%03d.root",mc_dir,tagh_counter), "READ" );
		
		TH1F *hbeam = (TH1F*)fSim->Get( "beam" )->Clone(Form("h_beam_tagh_%d",tagh_counter));
		double n_gen = 2.5e5  * (double)hbeam->GetMean();
		
		TH2F  *h2s  = new TH2F( Form("h2s_tagh_%d",tagh_counter), "#DeltaE vs. TAGH Counter", 
			274, 0.5, 274.5, 2000, -4.0, 4.0 );
		h2s->Add( (TH2F*)fSim->Get(hname_tagh) );
		
		TH1F  *h1s  = (TH1F*)h2s->ProjectionY( Form("h1s_tagh_%d",tagh_counter) );
		
		if( h1s->Integral() < 1.e1 || n_gen < 1.e3 ) {
			cout << "Skipping TAGH counter " << tagh_counter << endl;
			fSim->Close();
			continue;
		}
		
		double sim_int = (double)h1s->Integral();
		
		h1s->Rebin(rebins);
		h1s->SetTitle( Form("TAGH Counter %d (E_{#gamma} = %.3f)", tagh_counter, eb) );
		h1s->GetXaxis()->SetTitle( "E_{1} + E_{2} - E_{#gamma} [GeV]" );
		h1s->GetXaxis()->SetTitleOffset( 1.1 );
		h1s->GetYaxis()->SetTitle( Form("counts / %d MeV", n_mev) );
		h1s->GetYaxis()->SetTitleOffset( 1.3 );
		h1s->SetLineColor( kBlue );
		h1s->SetLineWidth( 2 );
		
		//---------   Fit background to a 3rd order polynomial   ----------//
		
		double bkgd_min = h1->GetBinCenter(h1->FindFirstBinAbove());
		double bkgd_max = h1->GetBinCenter(h1->FindLastBinAbove());
		
		TF1 *f_bkgd = new TF1( Form("f_bkgd_tagh_%d",tagh_counter), bkgd_fit, 
			bkgd_min, bkgd_max, 4 );
		h1->Fit( Form("f_bkgd_tagh_%d",tagh_counter), "R0QL" );
		
		
		//----------     Lineshape Fit     ----------//
		
		
		for( unsigned int ib = 0; ib < comp_sim_hist.size(); ib++ ) {
			comp_sim_hist[ib] = (1./sim_int) * (double)h1s->GetBinContent(ib+1);
		}
		
		
		double h1_max      = h1->GetBinCenter( h1->GetMaximumBin() );
		double h1_sim_max  = h1s->GetBinCenter( h1s->GetMaximumBin() );
		
		
		double loc_acc     = sim_int / n_gen;
		double scale_guess = f_theory->Eval(eb) * tagh_flux[tagh_counter-1]*loc_acc*ne*mb;
		
		
		TF1 *f_fit = new TF1( Form("f_fit_tagh_%d",tagh_counter), line_shape_fit, 
			bkgd_min, bkgd_max, 6 );
		
		f_fit->SetLineColor( kRed );
		f_fit->SetNpx( 500. );
		
		f_fit->SetParameters( scale_guess, h1_max - h1_sim_max );
		f_fit->SetParameter(1, 0.0);
		f_fit->SetParameter(2, f_bkgd->GetParameter(0));
		f_fit->SetParameter(3, f_bkgd->GetParameter(1));
		f_fit->SetParameter(4, f_bkgd->GetParameter(2));
		f_fit->SetParameter(5, f_bkgd->GetParameter(3));
		
		h1->Fit( Form("f_fit_tagh_%d",tagh_counter), "R0Q" );
		
		/*
		if( !FIX_BACKGROUND ) {
			f_fit->ReleaseParameter(2);
			f_fit->ReleaseParameter(3);
			f_fit->ReleaseParameter(4);
			f_fit->ReleaseParameter(5);
			
			h1->Fit( Form("f_fit_tagh_%d",tagh_counter), "R0QL" );
		}
		
		f_fit->FixParameter(2, f_fit->GetParameter(2));
		f_fit->FixParameter(3, f_fit->GetParameter(3));
		f_fit->FixParameter(4, f_fit->GetParameter(4));
		f_fit->FixParameter(5, f_fit->GetParameter(5));
		*/
		
		TFitResultPtr result = h1->Fit( Form("f_fit_tagh_%d",tagh_counter), "SR0QL" );
		//if( !((Int_t)result == 0) ) fit_val = 0;
		
		double chi2 = result->Chi2() / result->Ndf();
		
		canvas->cd();
		h1->Draw();
		f_fit->Draw("same");
		
		TF1 *f_comp_ls = new TF1( Form("f_comp_ls_tagh_%d",tagh_counter), line_shape_fit, 
			-4.0, 4.0, 6 );
		f_comp_ls->SetParameters( f_fit->GetParameter(0), f_fit->GetParameter(1), 
			0.0, 0.0, 0.0, 0.0 );
		f_comp_ls->SetLineColor(kBlue);
		
		
		TF1 *f_bkgd_ls = new TF1( Form("f_bkgd_ls_tagh_%d",tagh_counter), line_shape_fit,
			-4.0, 4.0, 6 );
		f_bkgd_ls->SetParameters( 0.0, 0.0, f_fit->GetParameter(2), f_fit->GetParameter(3), 
			f_fit->GetParameter(4), f_fit->GetParameter(5) );
		f_bkgd_ls->SetLineColor(kGreen);
		f_bkgd_ls->SetLineStyle(2);
		
		f_comp_ls->Draw("same");
		f_bkgd_ls->Draw("same");
		
		canvas->Update();
		
		if( DRAW_FITS_TAGH ) {
			if( loc_counter==0 ) canvas->Print( "deltaE_lineshpae_tagh.pdf(", "pdf" );
			else                 canvas->Print( "deltaE_lineshape_tagh.pdf",  "pdf" );
		}
		
		loc_counter++;
		
		fSim->Close();
	}
	
	if( DRAW_FITS_TAGH ) canvas->Print( "deltaE_lineshape_tagh.pdf)", "pdf" );
	
	
	
	//----------   Fit TAGM Counters   --------//
	
	vector<double>  en2Vec,  mu2Vec,  sig2Vec,  alpha2Vec,  n2Vec, chi22Vec;
	vector<double> en2EVec, mu2EVec, sig2EVec, alpha2EVec, n2EVec;
	
	loc_counter = 0;
	
	for( int tagm_counter = 1; tagm_counter <= 102; tagm_counter++ ) 
	{
		
		double eb = tagm_en[tagm_counter-1];
		
		if( tagm_flux[tagm_counter-1] <= 0. || tagm_flux_empty[tagm_counter-1] <= 0. ) {
			cout << "Skipping TAGM counter " << tagm_counter << endl;
			continue;
		}
		if( gSystem->AccessPathName(Form("%s/tagm_%03d.root",mc_dir,tagm_counter)) ) {
			cout << "Skipping TAGM counter " << tagm_counter << endl;
			continue;
		}
		
		double loc_flux_ratio = tagm_flux[tagm_counter-1] / tagm_flux_empty[tagm_counter-1];
		
		TH1F *h1  = (TH1F*)h2_tagm->ProjectionY(  Form("h1_tagm_%d", tagm_counter), 
			tagm_counter, tagm_counter );
		TH1F *h1e = (TH1F*)h2e_tagm->ProjectionY( Form("h1e_tagm_%d",tagm_counter), 
			tagm_counter, tagm_counter );
		
		if( h1->GetEntries() < 1.e1 ) {
			cout << "Skipping TAGM counter " << tagm_counter << endl;
			continue;
		}
		
		h1e->Scale( loc_flux_ratio );
		
		h1->Rebin(rebins);
		h1->SetTitle( Form("TAGM Counter %d (E_{#gamma} = %.3f)", tagm_counter, eb) );
		h1->GetXaxis()->SetTitle( "E_{1} + E_{2} - E_{#gamma} [GeV]" );
		h1->GetXaxis()->SetTitleOffset( 1.1 );
		h1->GetYaxis()->SetTitle( Form("counts / %d MeV", n_mev) );
		h1->GetYaxis()->SetTitleOffset( 1.3 );
		h1->SetLineColor( kBlack );
		h1->SetLineWidth( 2 );
		
		h1e->Rebin(rebins);
		h1e->GetXaxis()->SetTitle( "E_{1} + E_{2} - E_{#gamma} [GeV]" );
		h1e->GetXaxis()->SetTitleOffset( 1.1 );
		h1e->GetYaxis()->SetTitle( Form("counts / %d MeV", n_mev) );
		h1e->GetYaxis()->SetTitleOffset( 1.3 );
		h1e->SetLineColor( kBlue );
		h1e->SetLineWidth( 2 );
		
		h1->Add( h1e,-1.0 );
		
		
		//---------   Get Simulated histogram   --------//
		
		
		TFile *fSim = new TFile( Form("%s/tagm_%03d.root",mc_dir,tagm_counter), "READ" );
		
		TH1F *hbeam = (TH1F*)fSim->Get( "beam" )->Clone(Form("h_beam_tagm_%d",tagm_counter));
		double n_gen = 2.5e5  * (double)hbeam->GetMean();
		
		TH2F  *h2s  = new TH2F( Form("h2s_tagm_%d",tagm_counter), "#DeltaE vs. TAGM Counter", 
			102, 0.5, 102.5, 2000, -4.0, 4.0 );
		h2s->Add( (TH2F*)fSim->Get(hname_tagm) );
		
		TH1F  *h1s  = (TH1F*)h2s->ProjectionY( Form("h1s_tagm_%d",tagm_counter) );
		
		if( h1s->Integral() < 1.e1 || n_gen < 1.e3 ) {
			cout << "Skipping TAGM counter " << tagm_counter << endl;
			fSim->Close();
			continue;
		}
		
		double sim_int = (double)h1s->Integral();
		
		h1s->Rebin(rebins);
		h1s->SetTitle( Form("TAGM Counter %d (E_{#gamma} = %.3f)", tagm_counter, eb) );
		h1s->GetXaxis()->SetTitle( "E_{1} + E_{2} - E_{#gamma} [GeV]" );
		h1s->GetXaxis()->SetTitleOffset( 1.1 );
		h1s->GetYaxis()->SetTitle( Form("counts / %d MeV", n_mev) );
		h1s->GetYaxis()->SetTitleOffset( 1.3 );
		h1s->SetLineColor( kBlue );
		h1s->SetLineWidth( 2 );
		
		
		//---------   Fit background to a 3rd order polynomial   ----------//
		
		
		double bkgd_min = h1->GetBinCenter(h1->FindFirstBinAbove());
		double bkgd_max = h1->GetBinCenter(h1->FindLastBinAbove());
		
		TF1 *f_bkgd = new TF1( Form("f_bkgd_tagm_%d",tagm_counter), bkgd_fit, 
			bkgd_min, bkgd_max, 4 );
		h1->Fit( Form("f_bkgd_tagm_%d",tagm_counter), "R0QL" );
		
		
		//----------     Lineshape Fit     ----------//
		
		
		for( unsigned int ib = 0; ib < comp_sim_hist.size(); ib++ ) {
			comp_sim_hist[ib] = (1./sim_int) * (double)h1s->GetBinContent(ib+1);
		}
		
		
		double h1_max      = h1->GetBinCenter( h1->GetMaximumBin() );
		double h1_sim_max  = h1s->GetBinCenter( h1s->GetMaximumBin() );
		
		
		double loc_acc     = sim_int / n_gen;
		double scale_guess = f_theory->Eval(eb) * tagm_flux[tagm_counter-1]*loc_acc*ne*mb;
		
		
		TF1 *f_fit = new TF1( Form("f_fit_tagm_%d",tagm_counter), line_shape_fit, 
			bkgd_min, bkgd_max, 6 );
		
		f_fit->SetLineColor( kRed );
		f_fit->SetNpx( 500. );
		
		f_fit->SetParameters( scale_guess, h1_max - h1_sim_max );
		f_fit->SetParameter(1, 0.0);
		f_fit->SetParameter(2, f_bkgd->GetParameter(0));
		f_fit->SetParameter(3, f_bkgd->GetParameter(1));
		f_fit->SetParameter(4, f_bkgd->GetParameter(2));
		f_fit->SetParameter(5, f_bkgd->GetParameter(3));
		
		h1->Fit( Form("f_fit_tagm_%d",tagm_counter), "R0Q" );
		
		/*
		if( !FIX_BACKGROUND ) {
			f_fit->ReleaseParameter(2);
			f_fit->ReleaseParameter(3);
			f_fit->ReleaseParameter(4);
			f_fit->ReleaseParameter(5);
			
			h1->Fit( Form("f_fit_tagm_%d",tagm_counter), "R0QL" );
		}
		
		f_fit->FixParameter(2, f_fit->GetParameter(2));
		f_fit->FixParameter(3, f_fit->GetParameter(3));
		f_fit->FixParameter(4, f_fit->GetParameter(4));
		f_fit->FixParameter(5, f_fit->GetParameter(5));
		*/
		
		TFitResultPtr result = h1->Fit( Form("f_fit_tagm_%d",tagm_counter), "SR0QL" );
		//if( !((Int_t)result == 0) ) fit_val = 0;
		
		double chi2 = result->Chi2() / result->Ndf();
		
		canvas->cd();
		h1->Draw();
		f_fit->Draw("same");
		
		TF1 *f_comp_ls = new TF1( Form("f_comp_ls_tagm_%d",tagm_counter), line_shape_fit, 
			-4.0, 4.0, 6 );
		f_comp_ls->SetParameters( f_fit->GetParameter(0), f_fit->GetParameter(1), 
			0.0, 0.0, 0.0, 0.0 );
		f_comp_ls->SetLineColor(kBlue);
		
		
		TF1 *f_bkgd_ls = new TF1( Form("f_bkgd_ls_tagm_%d",tagm_counter), line_shape_fit,
			-4.0, 4.0, 6 );
		f_bkgd_ls->SetParameters( 0.0, 0.0, f_fit->GetParameter(2), f_fit->GetParameter(3), 
			f_fit->GetParameter(4), f_fit->GetParameter(5) );
		f_bkgd_ls->SetLineColor(kGreen);
		f_bkgd_ls->SetLineStyle(2);
		
		f_comp_ls->Draw("same");
		f_bkgd_ls->Draw("same");
		
		canvas->Update();
		
		if( DRAW_FITS_TAGM ) {
			if( loc_counter==0 ) canvas->Print( "deltaE_lineshpae_tagm.pdf(", "pdf" );
			else                 canvas->Print( "deltaE_lineshape_tagm.pdf",  "pdf" );
		}
		
		loc_counter++;
		
		fSim->Close();
		
	}
	
	if( DRAW_FITS_TAGM ) canvas->Print( "deltaE_lineshape_tagm.pdf)", "pdf" );
	
	
	
	/*
	//----------   Plot fit results   --------//
	
	int n_bins1 = (int)en1Vec.size();
	
	double *energy1  = new double[n_bins1];
	double *energy1E = new double[n_bins1];
	
	double *mu1      = new double[n_bins1];
	double *mu1E     = new double[n_bins1];
	
	double *sig1     = new double[n_bins1];
	double *sig1E    = new double[n_bins1];
	
	double *alpha1   = new double[n_bins1];
	double *alpha1E  = new double[n_bins1];
	
	double *nn1      = new double[n_bins1];
	double *nn1E     = new double[n_bins1];
	
	double *chi21    = new double[n_bins1];
	
	for( int i = 0; i < n_bins1; i++ ) {
		
		energy1[i]  = en1Vec[i];
		energy1E[i] = en1EVec[i];
		
		mu1[i]      = mu1Vec[i];
		mu1E[i]     = mu1EVec[i];
		
		sig1[i]     = sig1Vec[i];//  / en1Vec[i];
		sig1E[i]    = sig1EVec[i];// / en1Vec[i];
		
		alpha1[i]   = alpha1Vec[i];
		alpha1E[i]  = alpha1EVec[i];
		
		nn1[i]      = n1Vec[i];
		nn1E[i]     = n1EVec[i];
		
		chi21[i]    = chi21Vec[i];
		
	}
	
	
	int n_bins2 = (int)en2Vec.size();
	
	double *energy2  = new double[n_bins2];
	double *energy2E = new double[n_bins2];
	
	double *mu2      = new double[n_bins2];
	double *mu2E     = new double[n_bins2];
	
	double *sig2     = new double[n_bins2];
	double *sig2E    = new double[n_bins2];
	
	double *alpha2   = new double[n_bins2];
	double *alpha2E  = new double[n_bins2];
	
	double *nn2      = new double[n_bins2];
	double *nn2E     = new double[n_bins2];
	
	double *chi22    = new double[n_bins2];
	
	for( int i = 0; i < n_bins2; i++ ) {
		
		energy2[i]  = en2Vec[i];
		energy2E[i] = en2EVec[i];
		
		mu2[i]      = mu2Vec[i];
		mu2E[i]     = mu2EVec[i];
		
		sig2[i]     = sig2Vec[i];//  / en2Vec[i];
		sig2E[i]    = sig2EVec[i];// / en2Vec[i];
		
		alpha2[i]   = alpha2Vec[i];
		alpha2E[i]  = alpha2EVec[i];
		
		nn2[i]      = n2Vec[i];
		nn2E[i]     = n2EVec[i];
		
		chi22[i]    = chi22Vec[i];
		
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
	
	TGraphErrors *gAlpha1 = new TGraphErrors( n_bins1, energy1, alpha1, energy1E, alpha1E );
	gAlpha1->SetTitle( "#alpha_{#DeltaE} vs. E_{#gamma}" );
	gAlpha1->GetXaxis()->SetTitle( "E_{#gamma} [GeV]" );
	gAlpha1->GetYaxis()->SetTitle( "#alpha from Crystal Ball Fit [GeV]" );
	gAlpha1->GetYaxis()->SetTitleOffset( 1.3 );
	gAlpha1->SetMarkerStyle( 8 );
	gAlpha1->SetMarkerColor( kBlue );
	TGraphErrors *gAlpha2 = new TGraphErrors( n_bins2, energy2, alpha2, energy2E, alpha2E );
	gAlpha2->SetTitle( "#alpha_{#DeltaE} vs. E_{#gamma}" );
	gAlpha2->GetXaxis()->SetTitle( "E_{#gamma} [GeV]" );
	gAlpha2->GetYaxis()->SetTitle( "#alpha from Crystal Ball Fit [GeV]" );
	gAlpha2->GetYaxis()->SetTitleOffset( 1.3 );
	gAlpha2->SetMarkerStyle( 8 );
	gAlpha2->SetMarkerColor( kGreen );
	
	TGraphErrors *gN1 = new TGraphErrors( n_bins1, energy1, nn1, energy1E, nn1E );
	gN1->SetTitle( "n_{#DeltaE} vs. E_{#gamma}" );
	gN1->GetXaxis()->SetTitle( "E_{#gamma} [GeV]" );
	gN1->GetYaxis()->SetTitle( "n from Crystal Ball Fit [GeV]" );
	gN1->GetYaxis()->SetTitleOffset( 1.3 );
	gN1->SetMarkerStyle( 8 );
	gN1->SetMarkerColor( kBlue );
	TGraphErrors *gN2 = new TGraphErrors( n_bins2, energy2, nn2, energy2E, nn2E );
	gN2->SetTitle( "n_{#DeltaE} vs. E_{#gamma}" );
	gN2->GetXaxis()->SetTitle( "E_{#gamma} [GeV]" );
	gN2->GetYaxis()->SetTitle( "n from Crystal Ball Fit [GeV]" );
	gN2->GetYaxis()->SetTitleOffset( 1.3 );
	gN2->SetMarkerStyle( 8 );
	gN2->SetMarkerColor( kGreen );
	
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
	gMu1->Draw( "AP" );
	gMu2->Draw( "P same" );
	
	TCanvas *cSig = new TCanvas( "cSig", "cSig", 800, 400 );
	cSig->SetTickx(); cSig->SetTicky();
	gSig1->Draw( "AP" );
	gSig2->Draw( "P same" );
	
	TCanvas *cAlpha = new TCanvas( "cAlpha", "cAlpha", 800, 400 );
	cAlpha->SetTickx(); cAlpha->SetTicky();
	gAlpha1->Draw( "AP" );
	gAlpha2->Draw( "P same" );
	
	TF1 *fAlpha = new TF1( "fAlpha", "pol3", 6., 11.4 );
	gAlpha1->Fit( "fAlpha", "R0" );
	fAlpha->Draw("same");
	
	TCanvas *cN = new TCanvas( "cN", "cN", 800, 400 );
	cN->SetTickx(); cN->SetTicky();
	gN1->Draw( "AP" );
	gN2->Draw( "P same" );
	
	TF1 *fN = new TF1( "fN", "pol3", 6., 11.4 );
	gN1->Fit( "fN", "R0" );
	fN->Draw("same");
	
	TCanvas *cChi2 = new TCanvas( "cChi2", "cChi2", 800, 400 );
	cChi2->SetTickx(); cChi2->SetTicky();
	gChi21->Draw( "AP" );
	gChi22->Draw( "P same" );
	*/
	
	
	
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
	
	int nbins = (int)comp_sim_hist.size();
	
	Int_t xs = (int)((x[0]-par[1])/bin_size + (nbins/2));
	if( xs <   0 ) xs = 0;
	if( xs > nbins-1 ) xs = nbins-1;
	
	fitLS = par[0] * comp_sim_hist[xs];
	
	fq    = par[2] + 
		par[3] * x[0] + 
		par[4] * x[0] * x[0] + 
		par[5] * x[0] * x[0] * x[0];
	
	fitVal = fitLS + fq;
	
	return fitVal;
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




