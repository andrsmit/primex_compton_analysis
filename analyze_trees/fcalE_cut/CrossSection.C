
char        root_fname[256],       empty_target_root_fname[256];
char   tagh_flux_fname[256],  empty_target_tagh_flux_fname[256];
char   tagm_flux_fname[256],  empty_target_tagm_flux_fname[256];
char tagh_xscale_fname[256],             tagm_xscale_fname[256];
char   tagh_fabs_fname[256],               tagm_fabs_fname[256];
char    hname_tagh_sim[256],                hname_tagm_sim[256];
char            mc_dir[256];

double endpoint_energy, endpoint_energy_calib;

bool CS_FROM_FIT;
bool FIX_BACKGROUND;

double          tagh_en[274],          tagm_en[102];

double        tagh_flux[274],        tagm_flux[102];
double       tagh_fluxE[274],       tagm_fluxE[102];

double  tagh_flux_empty[274],  tagm_flux_empty[102];
double tagh_fluxE_empty[274], tagm_fluxE_empty[102];

double       tagh_yield[274],       tagm_yield[102];
double      tagh_yieldE[274],      tagm_yieldE[102];

double         tagh_acc[274],         tagm_acc[102];
double        tagh_accE[274],        tagm_accE[102];

double        tagh_fabs[274],        tagm_fabs[102];
double       tagh_fabsE[274],       tagm_fabsE[102];

double          tagh_cs[274],          tagm_cs[102];
double         tagh_csE[274],         tagm_csE[102];


//----------   Function Declarations   ----------//

void get_flux();
void get_counter_energies();
void get_fabs();

void calc_cs( int tag_sys, int counter, double &cs, double &csE );

int fit_yield( int tag_sys, int counter, TH1F *h1, double &yield, double &yieldE, double &chi2 );
int get_acc(   int tag_sys, int counter, double &acc, double &accE );
int fit_yield_ls( int tag_sys, int counter, TH1F *h1, double &yield, double &yieldE, double &chi2 );

Double_t bkgd_fit( Double_t *x, Double_t *par );
Double_t crys_ball_fit( Double_t *x, Double_t *par );
Double_t double_gaus_fit( Double_t *x, Double_t *par );
Double_t line_shape_fit( Double_t *x, Double_t *par );

//-----------------------------------------------//


double bin_size = 8. / 2000.;
int rebins, n_mev;

double ne, mb, neE;
TF1 *f_theory;

vector<double> comp_sim_hist;

TCanvas *canvas;
TPad *pad_lin, *pad_log;

TCanvas *canvas2;
TPad *top_pad, *bot_pad;



void CrossSection()
{
	
	
	double fcal_en = 0.95;
	int fcal_en_int = (int)(fcal_en*100);
	
	
	
	const char pathName[256] = "/work/halld/home/andrsmit/primex_compton_analysis";
	
	// file names for the full and empty target root files:
	
	sprintf( root_fname,              "%s/analyze_trees/fcalE_cut/rootFiles/Be.root", 
		pathName );
	sprintf( empty_target_root_fname, "%s/analyze_trees/fcalE_cut/rootFiles/Be_empty.root", 
		pathName );
	
	// photon flux file names:
	
	sprintf( tagh_flux_fname,              "%s/photon_flux/Be_tagh_flux.txt",       pathName );
	sprintf( tagm_flux_fname,              "%s/photon_flux/Be_tagm_flux.txt",       pathName );
	sprintf( empty_target_tagh_flux_fname, "%s/photon_flux/Be_empty_tagh_flux.txt", pathName );
	sprintf( empty_target_tagm_flux_fname, "%s/photon_flux/Be_empty_tagm_flux.txt", pathName );
	
	// file containing the xscales for the tagh and tagm counters:
	
	sprintf( tagh_xscale_fname, "%s/photon_flux/primex_tagh.txt", pathName );
	sprintf( tagm_xscale_fname, "%s/photon_flux/primex_tagm.txt", pathName );
	
	// files containing the absorption factors for the tagh and tagm counters:
	
	sprintf( tagh_fabs_fname, "%s/photon_absorption/Be_tagh_fabs.dat", pathName );
	sprintf( tagm_fabs_fname, "%s/photon_absorption/Be_tagm_fabs.dat", pathName );
	
	endpoint_energy_calib   = 11.6061;
	endpoint_energy         = 11.6061; // Be 200 nA data
	//endpoint_energy       = 11.1671; // He  50 nA data
	//endpoint_energy       = 11.1664; // He 100 nA data
	
	// Name of histograms for fitting yield:
	
	sprintf( hname_tagh_sim, "compton_systematics_sim/fcalE/deltaK_tagh_%02dfcalE", fcal_en_int );
	sprintf( hname_tagm_sim, "compton_systematics_sim/fcalE/deltaK_tagm_%02dfcalE", fcal_en_int );
	
	// Directory where mc rootFiles are stored:
	
	sprintf( mc_dir, "%s/compton_mc/recRootFiles_systematics", pathName );
	
	
	
	//---------------   Initialize   ---------------//
	
	
	gStyle->SetOptStat(0); gStyle->SetOptFit(0);
	
	rebins    =  5;
	bin_size *=  (double)rebins;
	n_mev     =  4 * rebins;
	
	int n_bins = (int)(2000 / rebins);
	
	comp_sim_hist.clear();
	for( int ib = 0; ib < n_bins; ib++ ) {
		comp_sim_hist.push_back( 0.0 );
	}
	
	
	
	// Adjustable switches:
	
	const bool DRAW_FITS_TAGH = false;
	const bool DRAW_FITS_TAGM = false;
	
	CS_FROM_FIT    = true;
	FIX_BACKGROUND = true;
	
	
	// Number of electrons in target:
	
	ne  = 8.77937e+23; // number of electrons per cm^2
	mb  = 1.e-27;
	
	neE = 0.015 * ne;
	
	
	for( int i = 0; i < 274; i++ ) {
		
		tagh_yield[i]  = 0.;
		tagh_yieldE[i] = 0.;
		tagh_acc[i]    = 0.;
		tagh_accE[i]   = 0.;
		tagh_cs[i]     = 0.;
		tagh_csE[i]    = 0.;
		
	}
	for( int i = 0; i < 102; i++ ) {
		
		tagm_yield[i]  = 0.;
		tagm_yieldE[i] = 0.;
		tagm_acc[i]    = 0.;
		tagm_accE[i]   = 0.;
		tagm_cs[i]     = 0.;
		tagm_csE[i]    = 0.;
		
	}
	
	
	
	//------------------------------------------------//
	
	
	
	get_flux();
	get_counter_energies();
	get_fabs();
	
	
	TFile *fFull  = new TFile( root_fname,               "READ" );
	TFile *fEmpty = new TFile( empty_target_root_fname,  "READ" );
	
	
	// Canvas to show fits on linear and log scales:
	
	canvas = new TCanvas( "canvas", "canvas", 1200, 600 );
	canvas->Divide( 2,1 );
	
	pad_lin = (TPad*)canvas->cd(1);
	pad_lin->SetTickx(); pad_lin->SetTicky();
	
	pad_log = (TPad*)canvas->cd(2);
	pad_log->SetTickx(); pad_log->SetTicky();
	pad_log->SetGrid();  pad_log->SetLogy();
	
	
	// Canvas to show deviation between fit and data:
	
	canvas2 = new TCanvas( "canvas2", "canvas2", 800, 800 );
	top_pad = new TPad( "top_pad", "top_pad", 0.005, 0.2025, 0.995, 0.995  );
	bot_pad = new TPad( "bot_pad", "bot_pad", 0.005, 0.005,  0.995, 0.1975 );
	
	top_pad->SetLeftMargin(0.10);
	top_pad->SetRightMargin(0.02);
	top_pad->SetTopMargin(0.075);
	top_pad->SetBottomMargin(0.015);
	top_pad->SetTickx(); top_pad->SetTicky();
	top_pad->SetFrameLineWidth(2);
	
	bot_pad->SetLeftMargin(0.10);
	bot_pad->SetRightMargin(0.02);
	bot_pad->SetTopMargin(0.005);
	bot_pad->SetBottomMargin(0.325);
	bot_pad->SetTickx(); bot_pad->SetTicky();
	bot_pad->SetFrameLineWidth(2);
	
	canvas2->cd();
	top_pad->Draw();
	bot_pad->Draw();
	
	
	
	
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
	
	
	
	TH1F  *h1_tagh[274],  *h1_tagm[102];
	TH1F *h1e_tagh[274], *h1e_tagm[102];
	
	for( int tagh_counter=1; tagh_counter<274; tagh_counter++ ) {
		h1_tagh[tagh_counter-1]  = 
			(TH1F*)fFull->Get(Form("%02dfcalE/tagh_%03d_%02dfcalE",fcal_en_int,
				tagh_counter,fcal_en_int))->Clone(Form("h1_tagh_%03d",tagh_counter));
		h1e_tagh[tagh_counter-1] = 
			(TH1F*)fEmpty->Get(Form("%02dfcalE/tagh_%03d_%02dfcalE",fcal_en_int,
				tagh_counter,fcal_en_int))->Clone(Form("h1e_tagh_%03d",tagh_counter));
	}
	for( int tagm_counter=1; tagm_counter<102; tagm_counter++ ) {
		h1_tagm[tagm_counter-1]  = 
			(TH1F*)fFull->Get(Form("%02dfcalE/tagm_%03d_%02dfcalE",fcal_en_int,
				tagm_counter,fcal_en_int))->Clone(Form("h1_tagm_%03d",tagm_counter));
		h1e_tagm[tagm_counter-1] = 
			(TH1F*)fEmpty->Get(Form("%02dfcalE/tagm_%03d_%02dfcalE",fcal_en_int,
				tagm_counter,fcal_en_int))->Clone(Form("h1e_tagm_%03d",tagm_counter));
	}
	
	
	//----------------------------------------------------//
	
	
	
	vector<int>  tagh_counterVec, tagm_counterVec;
	vector<double>  tagh_chi2Vec, tagm_chi2Vec;
	
	for( int tagh_counter = 1; tagh_counter <= 274; tagh_counter++ ) {
		
		double eb             = tagh_en[tagh_counter-1];
		double loc_flux       = tagh_flux[tagh_counter-1];
		double loc_flux_empty = tagh_flux_empty[tagh_counter-1];
		double loc_fluxE      = tagh_fluxE[tagh_counter-1];
		
		// Skip counters that have no flux:
		
		if( loc_flux <= 0. || loc_flux_empty <= 0. ) {
			cout << "Skipping TAGH counter " << tagh_counter << endl;
			continue;
		}
		
		// Skip counters that don't have a simulation file:
		
		if( gSystem->AccessPathName(Form("%s/tagh_%03d.root", mc_dir, tagh_counter)) ) {
			cout << "Skipping TAGH counter " << tagh_counter << endl;
			continue;
		}
		
		// Skip counters that have no yield:
		
		if( h1_tagh[tagh_counter-1]->Integral() < 1.e1 ) {
			cout << "Skipping TAGH counter " << tagh_counter << endl;
			continue;
		}
		
		h1_tagh[tagh_counter-1]->Add( 
			h1e_tagh[tagh_counter-1], -1.*loc_flux/loc_flux_empty );
		
		h1_tagh[tagh_counter-1]->SetTitle( 
			Form("TAGH Counter %d (E_{#gamma} = %.3f GeV)", tagh_counter, eb) );
		h1_tagh[tagh_counter-1]->Rebin( rebins );
		h1_tagh[tagh_counter-1]->SetLineColor( kBlack );
		h1_tagh[tagh_counter-1]->GetXaxis()->SetTitle( "E_{Comp} - E_{#gamma} [GeV]" );
		h1_tagh[tagh_counter-1]->SetLineWidth( 2 );
		h1_tagh[tagh_counter-1]->GetXaxis()->SetRangeUser( -3., 2. );
		
		
		double loc_yield = 0., loc_yieldE = 0., loc_chi2 = 0.;
		int fit_val = fit_yield( 0, tagh_counter, h1_tagh[tagh_counter-1], 
			loc_yield, loc_yieldE, loc_chi2 );
		if( fit_val <= 0 ) continue;
		
		double loc_acc = 0., loc_accE = 0.;;
		int acc_val  = get_acc( 0, tagh_counter, loc_acc, loc_accE );
		if( acc_val <= 0 ) continue;
		
		if( loc_acc <= 0. ) {
			cout << "zero acceptance in TAGH counter " << tagh_counter << endl;
			continue;
		}
		
		tagh_yield[tagh_counter-1]  = loc_yield;
		tagh_yieldE[tagh_counter-1] = loc_yieldE;
		
		tagh_acc[tagh_counter-1]    = loc_acc;
		tagh_accE[tagh_counter-1]   = loc_accE;
		
		tagh_counterVec.push_back( tagh_counter );
		tagh_chi2Vec.push_back( loc_chi2 );
		
		if( DRAW_FITS_TAGH ) {
			canvas->Update();
			canvas2->Update();
		}
		
	}
	
	
	for( int tagm_counter = 1; tagm_counter <= 102; tagm_counter++ ) {
		
		double eb             = tagm_en[tagm_counter-1];
		double loc_flux       = tagm_flux[tagm_counter-1];
		double loc_flux_empty = tagm_flux_empty[tagm_counter-1];
		double loc_fluxE      = tagm_fluxE[tagm_counter-1];
		
		// Skip counters that have no flux:
		
		if( loc_flux <= 0. || loc_flux_empty <= 0. ) {
			cout << "Skipping TAGM counter " << tagm_counter << endl;
			continue;
		}
		
		// Skip counters that don't have a simulation file:
		
		if( gSystem->AccessPathName(Form("%s/tagm_%03d.root", mc_dir, tagm_counter)) ) {
			cout << "Skipping TAGM counter " << tagm_counter << endl;
			continue;
		}
		
		// Skip counters that have no yield:
		
		if( h1_tagm[tagm_counter-1]->Integral() < 1.e1 ) {
			cout << "Skipping TAGM counter " << tagm_counter << endl;
			continue;
		}
		
		h1_tagm[tagm_counter-1]->Add( 
			h1e_tagm[tagm_counter-1], -1.*loc_flux/loc_flux_empty );
		
		h1_tagm[tagm_counter-1]->SetTitle( 
			Form("TAGM Counter %d (E_{#gamma} = %.3f GeV)", tagm_counter, eb) );
		h1_tagm[tagm_counter-1]->Rebin( rebins );
		h1_tagm[tagm_counter-1]->SetLineColor( kBlack );
		h1_tagm[tagm_counter-1]->GetXaxis()->SetTitle( "E_{Comp} - E_{#gamma} [GeV]" );
		h1_tagm[tagm_counter-1]->SetLineWidth( 2 );
		h1_tagm[tagm_counter-1]->GetXaxis()->SetRangeUser( -3., 2. );
		
		
		double loc_yield = 0., loc_yieldE = 0., loc_chi2 = 0.;
		int fit_val = fit_yield( 1, tagm_counter, h1_tagm[tagm_counter-1], 
			loc_yield, loc_yieldE, loc_chi2 );
		if( fit_val <= 0 ) continue;
		
		double loc_acc = 0., loc_accE = 0.;;
		int acc_val  = get_acc( 1, tagm_counter, loc_acc, loc_accE );
		if( acc_val <= 0 ) continue;
		
		if( loc_acc <= 0. ) {
			cout << "zero acceptance in TAGM counter " << tagm_counter << endl;
			continue;
		}
		
		tagm_yield[tagm_counter-1]  = loc_yield;
		tagm_yieldE[tagm_counter-1] = loc_yieldE;
		
		tagm_acc[tagm_counter-1]    = loc_acc;
		tagm_accE[tagm_counter-1]   = loc_accE;
		
		tagm_counterVec.push_back( tagm_counter );
		tagm_chi2Vec.push_back( loc_chi2 );
		
		if( DRAW_FITS_TAGM ) {
			canvas->Update();
			canvas2->Update();
		}
		
	}
	
	
	
	
	
	int n_bins1 = (int)tagh_counterVec.size();
	
	double *energy1  = new double[n_bins1];
	double *energy1E = new double[n_bins1];
	double *cs1      = new double[n_bins1];
	double *cs1E     = new double[n_bins1];
	
	double *dev1     = new double[n_bins1];
	double *dev1e    = new double[n_bins1];
	
	double *chi21    = new double[n_bins1];
	
	
	for( int ib = 0; ib < n_bins1; ib++ ) {
		
		int counter     = tagh_counterVec[ib];
		
		energy1[ib]     = tagh_en[tagh_counterVec[ib]-1];
		energy1E[ib]    = 0.;
		
		double locCS    = 0., locCSE = 0.;
		calc_cs( 0, counter, locCS, locCSE );
		
		cs1[ib]         = locCS;
		cs1E[ib]        = locCSE;
		
		tagh_cs[counter-1]  = locCS;
		tagh_csE[counter-1] = locCSE;
		
		double theoryCS = f_theory->Eval( tagh_en[counter-1] );
		
		dev1[ib]        = 100. * (locCS - theoryCS) / theoryCS;
		dev1e[ib]       = sqrt( pow(100.*locCSE/theoryCS,2.0) );
		
		chi21[ib]       = tagh_chi2Vec[ib];
	}
	
	
	int n_bins2 = (int)tagm_counterVec.size();
	
	double *energy2  = new double[n_bins2];
	double *energy2E = new double[n_bins2];
	double *cs2      = new double[n_bins2];
	double *cs2E     = new double[n_bins2];
	
	double *dev2     = new double[n_bins2];
	double *dev2e    = new double[n_bins2];
	
	double *chi22    = new double[n_bins2];
	
	for( int ib = 0; ib < n_bins2; ib++ ) {
		
		int counter     = tagm_counterVec[ib];
		
		energy2[ib]     = tagm_en[tagm_counterVec[ib]-1];
		energy2E[ib]    = 0.;
		
		double locCS    = 0., locCSE = 0.;
		calc_cs( 1, counter, locCS, locCSE );
		
		cs2[ib]         = locCS;
		cs2E[ib]        = locCSE;
		
		tagm_cs[counter-1]  = locCS;
		tagm_csE[counter-1] = locCSE;
		
		double theoryCS = f_theory->Eval( tagm_en[counter-1] );
		
		dev2[ib]        = 100. * (locCS - theoryCS) / theoryCS;
		dev2e[ib]       = sqrt( pow(100.*locCSE/theoryCS,2.0) );
			
		chi22[ib]       = tagm_chi2Vec[ib];
	}
	
	
	
	//-----   Plot Cross Section   -----//
	
	TGraphErrors *gCS1 = new TGraphErrors( n_bins1, energy1, cs1, energy1E, cs1E );
	gCS1->GetXaxis()->SetTitle( "Photon Beam Energy [GeV]" );
	gCS1->GetYaxis()->SetTitle( "#sigma [mb / electron]" );
	gCS1->SetTitle( "Compton Scattering Cross Section" );
	gCS1->SetMarkerStyle( 8 );
	gCS1->SetMarkerSize( 0.6 );
	gCS1->SetMarkerColor( kBlue );
	gCS1->SetLineColor( kBlue );
	
	TGraphErrors *gCS2 = new TGraphErrors( n_bins2, energy2, cs2, energy2E, cs2E );
	gCS2->GetXaxis()->SetTitle( "Photon Beam Energy [GeV]" );
	gCS2->GetYaxis()->SetTitle( "#sigma [mb / electron]" );
	gCS2->SetTitle( "Compton Scattering Cross Section" );
	gCS2->SetMarkerStyle( 8 );
	gCS2->SetMarkerSize( 0.6 );
	gCS2->SetMarkerColor( kGreen );
	gCS2->SetLineColor( kGreen );
	
	gCS1->GetYaxis()->SetRangeUser( 0.1, 0.26 );
	
	TCanvas *cCS = new TCanvas( "cCS", "cCS", 1200, 500 );
	cCS->SetTickx(); cCS->SetTicky();
	gCS1->Draw( "AP" );
	gCS2->Draw( "P same" );
	g_theory->Draw( "same" );
	
	
	//-----   Plot Cross Section Residuals  -----//
	
	TGraphErrors *gDev1 = new TGraphErrors( n_bins1, energy1, dev1, energy1E, dev1e );
	gDev1->SetTitle( "Relative Error of Compton Cross Section" );
	gDev1->GetXaxis()->SetTitle( "Photon Beam Energy [GeV]" );
	gDev1->GetYaxis()->SetTitle( "100 x (#sigma_{exp} - #sigma_{theory}) / #sigma_{theory} [%]" );
	gDev1->SetMarkerStyle( 8 );
	gDev1->SetMarkerSize( 0.8 );
	gDev1->SetMarkerColor( kBlue );
	gDev1->SetLineColor( kBlue );
	
	TGraphErrors *gDev2 = new TGraphErrors( n_bins2, energy2, dev2, energy2E, dev2e );
	gDev2->SetTitle( "Relative Error of Compton Cross Section" );
	gDev2->GetXaxis()->SetTitle( "Photon Beam Energy [GeV]" );
	gDev2->GetYaxis()->SetTitle( "(#sigma_{exp} - #sigma_{theory}) / #sigma_{theory}" );
	gDev2->SetMarkerStyle( 8 );
	gDev2->SetMarkerSize( 0.8 );
	gDev2->SetMarkerColor( kGreen );
	gDev2->SetLineColor( kGreen );
	
	gDev1->GetXaxis()->SetRangeUser(   6.0, 11.1 );
	gDev1->GetYaxis()->SetRangeUser( -20.0, 20.0 );
	
	TCanvas *cDev = new TCanvas( "cDev", "cDev", 1200, 350 );
	cDev->SetTickx(); cDev->SetTicky(); //cDev->SetGrid();
	gDev1->Draw( "AP" );
	gDev2->Draw( "P same" );
	
	cDev->Update();
	TLine *lDev = new TLine( gPad->GetUxmin(), 0.0, gPad->GetUxmax(), 0.0 );
	lDev->SetLineWidth(2);
	lDev->Draw("same");
	
	
	
	
	//-----   Plot Yield Fit Chi-2   -----//
	
	TGraph *gChi21 = new TGraph( n_bins1, energy1, chi21 );
	gChi21->GetXaxis()->SetTitle( "Photon Beam Energy [GeV]" );
	gChi21->GetYaxis()->SetTitle( "#chi^{2} / n.d.f" );
	gChi21->SetTitle( "#chi^{2} From Yield Fit" );
	gChi21->SetMarkerStyle(8);
	gChi21->SetMarkerColor( kBlue );
	
	TGraph *gChi22 = new TGraph( n_bins2, energy2, chi22 );
	gChi22->GetXaxis()->SetTitle( "Photon Beam Energy [GeV]" );
	gChi22->GetYaxis()->SetTitle( "#chi^{2} / n.d.f" );
	gChi22->SetTitle( "#chi^{2} From Yield Fit" );
	gChi22->SetMarkerStyle(8);
	gChi22->SetMarkerColor( kGreen );
	
	TCanvas *cChi2 = new TCanvas( "cChi2", "cChi2", 1000, 600 );
	cChi2->SetTickx(); cChi2->SetTicky();
	gChi21->Draw( "AP" );
	gChi22->Draw( "P same" );
	
	
	
	
	
	
	TCanvas *canvas3 = new TCanvas( "canvas3", "canvas3", 1200, 900 );
	TPad *top_pad3 = new TPad( "top_pad", "top_pad", 0.005, 0.3525, 0.995, 0.995  );
	TPad *bot_pad3 = new TPad( "bot_pad", "bot_pad", 0.005, 0.005,  0.995, 0.3475 );
	
	top_pad3->SetLeftMargin(0.10);
	top_pad3->SetRightMargin(0.02);
	top_pad3->SetTopMargin(0.075);
	top_pad3->SetBottomMargin(0.015);
	top_pad3->SetTickx(); top_pad3->SetTicky();
	top_pad3->SetFrameLineWidth(2);
	
	bot_pad3->SetLeftMargin(0.10);
	bot_pad3->SetRightMargin(0.02);
	bot_pad3->SetTopMargin(0.010);
	bot_pad3->SetBottomMargin(0.225);
	bot_pad3->SetTickx(); bot_pad3->SetTicky();
	bot_pad3->SetFrameLineWidth(2);
	
	canvas3->cd();
	top_pad3->Draw();
	bot_pad3->Draw();
	
	gCS1->GetYaxis()->SetRangeUser(0.,0.35);
	gCS1->GetXaxis()->SetRangeUser(6.,11.);
	gCS1->GetXaxis()->SetTitleOffset(1);
	gCS1->GetXaxis()->SetLabelOffset(1);
	gCS1->GetYaxis()->SetTitleSize(0.045);
	gCS1->GetYaxis()->SetTitleOffset(0.8);
	
	gDev1->GetXaxis()->SetTitleSize(0.065);
	gDev1->GetXaxis()->SetLabelSize(0.06);
	gDev1->GetYaxis()->SetTitleSize(0.07);
	gDev1->GetYaxis()->SetTitleOffset(0.4);
	gDev1->GetYaxis()->SetTitle( "(#sigma_{exp}-#sigma_{theory})/#sigma_{theory} [%]");
	gDev1->SetTitle( " " );
	
	top_pad3->cd();
	gCS1->Draw("AP");
	gCS2->Draw("P same");
	g_theory->Draw("same");
	
	bot_pad3->cd();
	gDev1->Draw("AP");
	gDev2->Draw("P same");
	
	bot_pad3->Update();
	TLine *lDev3 = new TLine( gPad->GetUxmin(), 0.0, gPad->GetUxmax(), 0.0 );
	lDev3->SetLineWidth(2);
	lDev3->Draw("same");
	
	
	
	
	char buf[256];
	
	ofstream outf_tagh( Form("tagh_fcalE_%02d.txt", fcal_en_int) );
	for( int tagh_counter = 1; tagh_counter <= 274; tagh_counter++ ) {
		sprintf( buf, "%03d   %2.4f   %f   %f   %f   %f   %f   %f", tagh_counter, 
			tagh_en[tagh_counter-1], tagh_cs[tagh_counter-1], tagh_csE[tagh_counter-1], 
			tagh_yield[tagh_counter-1], tagh_yieldE[tagh_counter-1], 
			tagh_acc[tagh_counter-1], tagh_accE[tagh_counter-1] );
		outf_tagh << buf << "\n";
	}
	outf_tagh.close();
	
	ofstream outf_tagm( Form("tagm_fcalE_%02d.txt", fcal_en_int) );
	for( int tagm_counter = 1; tagm_counter <= 102; tagm_counter++ ) {
		sprintf( buf, "%03d   %2.4f   %f   %f   %f   %f   %f   %f", tagm_counter, 
			tagm_en[tagm_counter-1], tagm_cs[tagm_counter-1], tagm_csE[tagm_counter-1], 
			tagm_yield[tagm_counter-1], tagm_yieldE[tagm_counter-1], 
			tagm_acc[tagm_counter-1], tagm_accE[tagm_counter-1] );
		outf_tagm << buf << "\n";
	}
	outf_tagm.close();
	
	
	
	
	
	return;
}




void get_fabs() 
{
	
	int a; double b, c, d;
	
	for( int i=0; i<274; i++ ) {
		tagh_fabs[i]  = 0.;
		tagh_fabsE[i] = 0.;
	}
	for( int i=0; i<102; i++ ) {
		tagm_fabs[i]  = 0.;
		tagm_fabsE[i] = 0.;
	}
	
	ifstream inf1( tagh_fabs_fname );
	for( int i=0; i<274; i++ ) {
		inf1 >> a >> b >> c >> d;
		tagh_fabs[i]  = b;
		tagh_fabsE[i] = fabs(0.5*(c-d));
	}
	inf1.close();
	
	ifstream inf2( tagm_fabs_fname );
	for( int i=0; i<102; i++ ) {
		inf2 >> a >> b >> c >> d;
		tagm_fabs[i]  = b;
		tagm_fabsE[i] = fabs(0.5*(c-d));
	}
	inf2.close();
	
	
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



void calc_cs( int tag_sys, int counter, double &loc_cs, double &loc_csE )
{
	
	loc_cs = 0., loc_csE = 0.;
	
	double loc_yield,  loc_flux,  loc_acc,  loc_fabs;
	double loc_yieldE, loc_fluxE, loc_accE, loc_fabsE;
	
	if( tag_sys==0 ) {
		
		loc_yield  = tagh_yield[counter-1];
		loc_yieldE = tagh_yieldE[counter-1];
		
		loc_flux   = tagh_flux[counter-1];
		loc_fluxE  = tagh_fluxE[counter-1];
		
		loc_acc    = tagh_acc[counter-1];
		loc_accE   = tagh_accE[counter-1];
		
		loc_fabs   = 1.0 - tagh_fabs[counter-1];
		loc_fabsE  = tagh_fabsE[counter-1];
		
	} else {
		
		loc_yield  = tagm_yield[counter-1];
		loc_yieldE = tagm_yieldE[counter-1];
		
		loc_flux   = tagm_flux[counter-1];
		loc_fluxE  = tagm_fluxE[counter-1];
		
		loc_acc    = tagm_acc[counter-1];
		loc_accE   = tagm_accE[counter-1];
		
		loc_fabs   = 1.0 - tagm_fabs[counter-1];
		loc_fabsE  = tagm_fabsE[counter-1];
		
	}
	
	loc_cs   = loc_yield / (loc_flux * loc_acc * ne * mb);
	loc_csE  = sqrt( 
			  pow( loc_yieldE / (loc_flux*loc_acc*ne*mb), 2.0)
			+ pow( loc_yield*loc_fluxE/(loc_flux*loc_flux*loc_acc*ne*mb), 2.0 )
			+ pow( loc_yield*loc_accE/(loc_flux*loc_acc*loc_acc*ne*mb), 2.0 )
			+ pow( loc_yield*neE/(loc_flux*loc_acc*ne*ne*mb), 2.0 ) );
	
	loc_cs  /= loc_fabs;
	
	
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



Double_t line_shape_fit( Double_t *x, Double_t *par )
{
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



int fit_yield( int tag_sys, int counter, TH1F *h1, double &yield, double &yieldE, double &chi2 )
{
	
	/*
	Steps to fitting for the yield:
	
	1. Get simulated deltaK histogram and fit lineshape to the data with a polynomial background
	2. Use a double gaussian fit with polynomial background to fit for yield. The background
	    is fixed from the fit in step 1.
	
	*/
	
	
	int fit_val = 1;
	
	char tag_sys_char[256];
	if( tag_sys==0 ) sprintf( tag_sys_char, "tagh" );
	else             sprintf( tag_sys_char, "tagm" );
	
	
	double min_fit_x = h1->GetBinCenter(h1->FindFirstBinAbove());
	double max_fit_x = h1->GetBinCenter(h1->FindLastBinAbove());
	
	
	//-----   Get lineshape from simulation   -----//
	
	char fname[256];
	if( tag_sys==0 ) 
		sprintf( fname, "%s/tagh_%03d.root", mc_dir, counter );
	else 
		sprintf( fname, "%s/tagm_%03d.root", mc_dir, counter );
	
	if( gSystem->AccessPathName(fname) ) {
		return 0;
	}
	
	TFile *fSim = new TFile( fname, "READ" );
	
	TH1F *hbeam  = (TH1F*)fSim->Get( "compton_systematics_sim/beam" )->Clone( 
		Form("h_beam_%d_%d",tag_sys,counter) );
	double n_gen = 2.5e5 * (double)hbeam->GetMean();
	
	if( n_gen < 1.e3 ) return 0;
	
	TH2F *h2;
	
	if( tag_sys == 0 ) {
		
		h2 = new TH2F( Form("h2_tagh_%d",counter), "DeltaK", 
			274, 0.5, 274.5, 2000, -4.0, 4.0 );
		
		h2->Add( (TH2F*)fSim->Get(Form("%s",hname_tagh_sim)) );
	} else {
		
		h2 = new TH2F( Form("h2_tagm_%d",counter), "DeltaK", 
			102, 0.5, 102.5, 2000, -4.0, 4.0 );
		
		h2->Add( (TH2F*)fSim->Get(Form("%s",hname_tagm_sim)) );
	}
	
	h2->SetDirectory(0);
	
	fSim->Close();
	
	TH1F *h1_sim = (TH1F*)h2->ProjectionY( Form("h1_sim_%d_%d",tag_sys,counter) );
	if( h1_sim->Integral() < 1.e1 ) return 0;
	
	double sim_int = (double)h1_sim->Integral();
	
	h1_sim->Rebin(rebins);
	h1_sim->GetXaxis()->SetTitle( "E_{Comp} - E_{#gamma} [GeV]" );
	h1_sim->GetXaxis()->SetTitleOffset( 1.1 );
	h1_sim->GetYaxis()->SetTitle( Form("counts / %d MeV", n_mev) );
	h1_sim->GetYaxis()->SetTitleOffset( 1.3 );
	h1_sim->SetLineColor( kBlack );
	h1_sim->SetLineWidth( 2 );
	
	for( unsigned int ib = 0; ib < comp_sim_hist.size(); ib++ ) {
		comp_sim_hist[ib] = (1./sim_int) * (double)h1_sim->GetBinContent(ib+1);
	}
	
	double h1_max      = h1->GetBinCenter( h1->GetMaximumBin() );
	double h1_sim_max  = h1_sim->GetBinCenter( h1_sim->GetMaximumBin() );
	
	double eb, loc_flux;
	if( tag_sys==0 ) { 
		eb       = tagh_en[counter-1];
		loc_flux = tagh_flux[counter-1];
	} else {
		eb       = tagm_en[counter-1];
		loc_flux = tagm_flux[counter-1];
	}
	
	double loc_acc = sim_int / n_gen;
	
	double scale_guess = f_theory->Eval(eb) * loc_flux * loc_acc * ne * mb;
	
	TF1 *f_ls_fit  = new TF1( Form("f_ls_fit_%s_%d",tag_sys_char,counter), 
		line_shape_fit, min_fit_x, max_fit_x, 6 );
	
	f_ls_fit->SetLineColor( kRed );
	f_ls_fit->SetNpx( 500. );
	
	f_ls_fit->SetParameters( scale_guess, h1_max - h1_sim_max );
	f_ls_fit->SetParameter(1, 0.0);
	f_ls_fit->FixParameter(2, 0.0);
	f_ls_fit->FixParameter(3, 0.0);
	f_ls_fit->FixParameter(4, 0.0);
	f_ls_fit->FixParameter(5, 0.0);
	
	h1->Fit( Form("f_ls_fit_%s_%d",tag_sys_char,counter), "R0Q" );
	
	f_ls_fit->FixParameter(0, f_ls_fit->GetParameter(0));
	f_ls_fit->FixParameter(1, f_ls_fit->GetParameter(1));
	f_ls_fit->ReleaseParameter(2);
	f_ls_fit->ReleaseParameter(3);
	f_ls_fit->ReleaseParameter(4);
	f_ls_fit->ReleaseParameter(5);
	
	h1->Fit( Form("f_ls_fit_%s_%d",tag_sys_char,counter), "R0QL" );
	
	f_ls_fit->ReleaseParameter(0);
	f_ls_fit->ReleaseParameter(1);
	
	h1->Fit( Form("f_ls_fit_%s_%d",tag_sys_char,counter), "R0QL" );
	
	
	/*
	//-----   Local Gaussian fit to get initial parameter guesses   -----//
	
	TF1 *f_gaus = new TF1( Form("f_gaus_%s_%d",tag_sys_char,counter), "gaus", 
		min_fit_x, max_fit_x );
	
	f_gaus->SetParameters( h1->GetMaximum(), h1->GetBinCenter(h1->GetMaximumBin()), 0.2 );
	f_gaus->SetParLimits( 2, 0.0, 1.0 );
	
	h1->Fit( Form("f_gaus_%s_%d",tag_sys_char,counter), "R0QL" );
	
	f_gaus->SetRange( f_gaus->GetParameter(1)-0.2, f_gaus->GetParameter(1)+0.2 );
	
	h1->Fit( Form("f_gaus_%s_%d",tag_sys_char,counter), "R0QL" );
	
	
	//-----   Full Yield Fit   -----//
	
	TF1 *f_fit  = new TF1( Form("f_fit_%s_%d",tag_sys_char,counter), 
		double_gaus_fit, min_fit_x, max_fit_x, 10 );
	
	f_fit->SetParName( 0, "A" );
	f_fit->SetParName( 1, "z" );
	f_fit->SetParName( 2, "#mu_{1}" );
	f_fit->SetParName( 3, "#mu_{1}-#mu_{2}" );
	f_fit->SetParName( 4, "#sigma_{1}" );
	f_fit->SetParName( 5, "#sigma_{2}" );
	
	f_fit->SetParLimits( 1, 0.0, 0.5 );
	f_fit->SetParLimits( 4, 0.0, 1.0 );
	f_fit->SetParLimits( 5, 0.0, 2.0 );
	
	f_fit->SetParameters( f_gaus->GetParameter(0), 0., f_gaus->GetParameter(1), 
		0., f_gaus->GetParameter(2), 2.*f_gaus->GetParameter(2) );
	
	f_fit->FixParameter(6, f_ls_fit->GetParameter(2) );
	f_fit->FixParameter(7, f_ls_fit->GetParameter(3) );
	f_fit->FixParameter(8, f_ls_fit->GetParameter(4) );
	f_fit->FixParameter(9, f_ls_fit->GetParameter(5) );
	
	h1->Fit( Form("f_fit_%s_%d",tag_sys_char,counter), "R0QL" );
	
	if( !FIX_BACKGROUND ) {
		f_fit->ReleaseParameter(6);
		f_fit->ReleaseParameter(7);
		f_fit->ReleaseParameter(8);
		f_fit->ReleaseParameter(9);
		
		h1->Fit( Form("f_fit_%s_%d",tag_sys_char,counter), "R0QL" );
	}
	
	TFitResultPtr result = h1->Fit( Form("f_fit_%s_%d",tag_sys_char,counter), "SR0QL" );
	if( !((Int_t)result == 0) ) return 0;
	
	chi2 = result->Chi2() / result->Ndf();
	
	h1->GetXaxis()->SetRangeUser( min_fit_x, max_fit_x );
	
	
	if( CS_FROM_FIT ) {
	
	yield  = f_fit->GetParameter(0) * sqrt(2.*TMath::Pi()) 
		* ( (1.-f_fit->GetParameter(1))*f_fit->GetParameter(4)
		+ f_fit->GetParameter(1)*f_fit->GetParameter(5) );
	yield /= bin_size;
	
	yieldE = sqrt(yield);
	
	} else {
			
	TF1 *fInt     = new TF1( Form("fInt_tagh_%d",counter), "pol4", -3., 3. );
	fInt->SetParameters( 0., f_fit->GetParameter(6), 0.5*f_fit->GetParameter(7), 
		(1./3.)*f_fit->GetParameter(8), 0.25*f_fit->GetParameter(9) );
	
	double intsub = fInt->Eval(1.5) - fInt->Eval(-2.);
	intsub       /= bin_size;
	
	yield  = h1->Integral();// - intsub;
	yieldE = sqrt(yield);
	
	}
	
	
	TF1 *f_draw = new TF1( Form("f_draw_%d_%d",tag_sys,counter), 
		"pol3", min_fit_x, max_fit_x );
	f_draw->SetParameter(0, f_fit->GetParameter(6));
	f_draw->SetParameter(1, f_fit->GetParameter(7));
	f_draw->SetParameter(2, f_fit->GetParameter(8));
	f_draw->SetParameter(3, f_fit->GetParameter(9));
	f_draw->SetLineColor(kGreen);
	f_draw->SetLineStyle(2);
	
	
	
	pad_lin->cd();
	h1->Draw();
	f_fit->Draw("same");
	f_draw->Draw("same");
	
	TH1F *h1_log = (TH1F*)h1->Clone( Form("h1_log_%d_%d",tag_sys,counter) );
	h1_log->SetMinimum(1.0);
	
	pad_log->cd();
	h1_log->Draw();
	f_fit->Draw("same");
	f_draw->Draw("same");
	
	
	TH1F *h1_dev = new TH1F( Form("h1_dev_%d_%d",tag_sys,counter), "", 
		h1->GetXaxis()->GetNbins(), -4.0, 4.0 );
	
	for( int ib=1; ib<=h1_dev->GetXaxis()->GetNbins(); ib++ ) {
		double loc_x   = h1->GetBinCenter(ib);
		double loc_c   = h1->GetBinContent(ib);
		double loc_err;
		if( loc_c > 0. ) loc_err = sqrt(loc_c);
		else loc_err = 1.0;
		
		double loc_dev = (loc_c - f_fit->Eval(loc_x))/loc_err;
		h1_dev->SetBinContent(ib,loc_dev);
	}
	
	h1->SetTitleFont(22);
	h1->GetXaxis()->SetRangeUser( min_fit_x, max_fit_x );
	h1->GetXaxis()->SetTitle( "" );
	h1->GetXaxis()->SetLabelOffset(1.0);
	h1->GetYaxis()->SetTitle( Form("counts / %d MeV",n_mev) );
	h1->GetYaxis()->SetTitleSize( 0.045 );
	h1->GetYaxis()->SetTitleOffset( 1.0 );
	h1->GetYaxis()->SetTitleFont( 22 );
	h1->SetMinimum(0.);
	
	h1_dev->GetYaxis()->SetTitle( "(data - fit)/err" );
	h1_dev->GetYaxis()->SetTitleSize( 0.125 );
	h1_dev->GetYaxis()->SetTitleOffset( 0.275 );
	h1_dev->GetYaxis()->SetLabelSize( 0.10 );
	h1_dev->GetYaxis()->SetTitleFont( 22 );
	
	h1_dev->GetXaxis()->SetTitle( "E_{Comp} - E_{#gamma} [GeV]" );
	h1_dev->GetXaxis()->SetTitleSize( 0.175 );
	h1_dev->GetXaxis()->SetTitleOffset( 0.90 );
	h1_dev->GetXaxis()->SetLabelSize( 0.13 );
	h1_dev->GetXaxis()->SetTitleFont( 22 );
	
	h1_dev->SetLineWidth(2);
	h1_dev->SetLineColor(kBlack);
	h1_dev->SetMarkerColor(kBlack);
	h1_dev->GetXaxis()->SetRangeUser( min_fit_x, max_fit_x );
	h1_dev->SetTitle( " " );
	h1_dev->GetYaxis()->SetRangeUser( -12., 12. );
	
	top_pad->cd();
	h1->Draw();
	f_fit->Draw("same");
	
	bot_pad->cd();
	h1_dev->Draw("PE");
	*/
	//canvas2->Update();
	
	
	
	pad_lin->cd();
	h1->Draw();
	f_ls_fit->Draw("same");
	
	
	yield  = f_ls_fit->GetParameter(0);
	yieldE = f_ls_fit->GetParError(0);
	
	
	return fit_val;
}




int get_acc(  int tag_sys, int counter, double &acc, double &accE )
{
	
	int acc_val = 1;
	
	char fname[256];
	if( tag_sys==0 ) 
		sprintf( fname, "%s/tagh_%03d.root", mc_dir, counter );
	else 
		sprintf( fname, "%s/tagm_%03d.root", mc_dir, counter );
	
	
	if( gSystem->AccessPathName(fname) ) {
		return 0;
	}
	
	
	TFile *fSim = new TFile( fname, "READ" );
	
	TH1F *hbeam  = (TH1F*)fSim->Get( "compton_systematics_sim/beam" )->Clone(
		Form("h_beam_%d_%d",tag_sys,counter) );
	double n_gen = 2.5e5 * (double)hbeam->GetMean();
	
	TH2F *h2;
	
	if( tag_sys == 0 ) {
		
		h2 = new TH2F( Form("h2_tagh_%d",counter), "DeltaK", 
			274, 0.5, 274.5, 2000, -4.0, 4.0 );
		
		h2->Add( (TH2F*)fSim->Get(Form("%s",hname_tagh_sim)) );
	} else {
		
		h2 = new TH2F( Form("h2_tagm_%d",counter), "DeltaK", 
			102, 0.5, 102.5, 2000, -4.0, 4.0 );
		
		h2->Add( (TH2F*)fSim->Get(Form("%s",hname_tagm_sim)) );
	}
	
	h2->SetDirectory(0);
	
	fSim->Close();
	
	
	TH1F *h1 = (TH1F*)h2->ProjectionY( Form("h1_sim_%d_%d",tag_sys,counter) );
	if( h1->Integral() < 1.e1 ) return 0;
	
	
	h1->Rebin(rebins);
	h1->GetXaxis()->SetTitle( "E_{Comp} - E_{#gamma} [GeV]" );
	h1->GetXaxis()->SetTitleOffset( 1.1 );
	h1->GetXaxis()->SetRangeUser( -1.5, 1.5 );
	h1->GetYaxis()->SetTitle( Form("counts / %d MeV", n_mev) );
	h1->GetYaxis()->SetTitleOffset( 1.3 );
	h1->SetLineColor( kBlack );
	h1->SetLineWidth( 2 );
	h1->GetXaxis()->SetRangeUser( -3., 2. );
	
	
	
	TF1 *f_gaus = new TF1( Form("f_gaus_%d_%d",tag_sys,counter), "gaus", -2., 2. );
	
	f_gaus->SetParameters( h1->GetMaximum(), 
		h1->GetBinCenter(h1->GetMaximumBin()), 0.2 );
	f_gaus->SetParLimits( 2, 0., 1. );
	
	h1->Fit( Form("f_gaus_%d_%d",tag_sys,counter), "R0QL" );
	
	f_gaus->SetRange( f_gaus->GetParameter(1)-0.2, f_gaus->GetParameter(1)+0.2 );
	
	h1->Fit( Form("f_gaus_%d_%d",tag_sys,counter), "R0QL" );
	
	
	TF1 *f_fit = new TF1( Form("f_fit_%d_%d",tag_sys,counter), double_gaus_fit, -2., 2., 10 );
	
	f_fit->SetParName( 0, "A" );
	f_fit->SetParName( 1, "z" );
	f_fit->SetParName( 2, "#mu_{1}" );
	f_fit->SetParName( 3, "#mu_{1}-#mu_{2}" );
	f_fit->SetParName( 4, "#sigma_{1}" );
	f_fit->SetParName( 5, "#sigma_{2}" );
	
	f_fit->SetParLimits( 1, 0.0, 0.5 );
	f_fit->SetParLimits( 4, 0.0, 1.0 );
	f_fit->SetParLimits( 5, 0.0, 2.0 );
	
	f_fit->SetParameters( f_gaus->GetParameter(0), 0., f_gaus->GetParameter(1), 
		0., f_gaus->GetParameter(2), 2.*f_gaus->GetParameter(2) );
	
	f_fit->FixParameter( 6, 0.0 );
	f_fit->FixParameter( 7, 0.0 );
	f_fit->FixParameter( 8, 0.0 );
	f_fit->FixParameter( 9, 0.0 );
	
	f_fit->SetRange( -2.0, 2.0 );
	
	h1->Fit( Form("f_fit_%d_%d",tag_sys,counter), "R0QL" );
	
	f_fit->SetRange( -2., 1.5 );
	
	TFitResultPtr result = h1->Fit( Form("f_fit_%d_%d",tag_sys,counter), "SR0QL" );
	if( !((Int_t)result == 0) ) return 0;
	
	h1->Fit( Form("f_fit_%d_%d",tag_sys,counter), "R0QL" );
	
	
	
	
	
	double yield = f_fit->GetParameter(0) * sqrt(2.*TMath::Pi()) 
		* ( (1.-f_fit->GetParameter(1))*f_fit->GetParameter(4)
		+ f_fit->GetParameter(1)*f_fit->GetParameter(5) );
	yield /= bin_size;
	
	yield = h1->Integral();
	
	acc  = yield / n_gen;
	accE = sqrt( n_gen*acc*(1.-acc) ) / n_gen;
	
	
	
	return acc_val;
}
