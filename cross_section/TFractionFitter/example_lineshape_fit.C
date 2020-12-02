
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

double tagh_acc[274], tagh_accE[274];
double tagm_acc[102], tagm_accE[102];

//----------   Function Declarations   ----------//

void get_flux();
void get_counter_energies();

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

double ne, mb;
TF1 *f_theory;

vector<double> comp_sim_hist;
vector<double> pair_sim_hist;
vector<double> trip_sim_hist;

TCanvas *canvas1;
TPad *pad_lin, *pad_log;

TCanvas *canvas2;
TPad *top_pad, *bot_pad;



void example_lineshape_fit() 
{
	
	
	const char pathName[256] = "/work/halld/home/andrsmit/primex_compton_analysis";
	
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
	
	sprintf( hname_tagh, "DeltaK/deltaK_tagh_ep" );
	sprintf( hname_tagm, "DeltaK/deltaK_tagm_ep" );
	
	// Directory where mc rootFiles are stored:
	
	sprintf( mc_dir, "%s/compton_mc/recRootFiles", pathName );
	
	
	
	//---------------   Initialize   ---------------//
	
	
	gStyle->SetOptStat(0); gStyle->SetOptFit(0);
	
	rebins    =  5;
	bin_size *=  (double)rebins;
	n_mev     =  4 * rebins;
	
	int n_bins = (int)(2000 / rebins);
	
	comp_sim_hist.clear();
	pair_sim_hist.clear();
	trip_sim_hist.clear();
	
	for( int ib = 0; ib < n_bins; ib++ ) {
		comp_sim_hist.push_back( 0.0 );
		pair_sim_hist.push_back( 0.0 );
		trip_sim_hist.push_back( 0.0 );
	}
	
	
	double f_abs = 1.0;
	
	
	// Adjustable switches:
	
	const bool DRAW_FITS_TAGH = false;
	const bool DRAW_FITS_TAGM = false;
	
	CS_FROM_FIT    = true;
	FIX_BACKGROUND = false;
	
	
	// Number of electrons in target:
	
	ne = 8.77937e+23; // number of electrons per cm^2
	mb = 1.e-27;
	
	
	//------------------------------------------------//
	
	
	
	get_flux();
	get_counter_energies();
	
	
	TFile *fFull  = new TFile( root_fname,               "READ" );
	TFile *fEmpty = new TFile( empty_target_root_fname,  "READ" );
	
	
	// Canvas to show fits on linear and log scales:
	
	canvas1 = new TCanvas( "canvas1", "canvas1", 1200, 600 );
	canvas1->Divide( 2,1 );
	
	pad_lin = (TPad*)canvas1->cd(1);
	pad_lin->SetTickx(); pad_lin->SetTicky();
	
	pad_log = (TPad*)canvas1->cd(2);
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
	
	
	// Data:
	
	TH2F *h2_tagh  = new TH2F( "h2_tagh",  "#DeltaK vs. Counter", 
		274, 0.5, 274.5,  2000, -4.0, 4.0 );
	TH2F *h2_tagm  = new TH2F( "h2_tagm",  "#DeltaK vs. Counter", 
		102, 0.5, 102.5,  2000, -4.0, 4.0 );
	
	TH2F *h2e_tagh = new TH2F( "h2e_tagh", "#DeltaK vs. Counter", 
		274, 0.5, 274.5,  2000, -4.0, 4.0 );
	TH2F *h2e_tagm = new TH2F( "h2e_tagm", "#DeltaK vs. Counter", 
		102, 0.5, 102.5,  2000, -4.0, 4.0 );
	
	h2_tagh->Add( (TH2F*)fFull->Get(hname_tagh) );
	h2_tagm->Add( (TH2F*)fFull->Get(hname_tagm) );
	h2e_tagh->Add( (TH2F*)fEmpty->Get(hname_tagh) );
	h2e_tagm->Add( (TH2F*)fEmpty->Get(hname_tagm) );
	
	
	// Pair Production Simulation:
	
	TFile *fPair = new TFile( 
	"/work/halld/home/andrsmit/compton_analysis/pair_sim/rootFiles_pair/Be.root", "READ" );
	TFile *fTrip = new TFile( 
	"/work/halld/home/andrsmit/compton_analysis/pair_sim/rootFiles_trip/Be.root", "READ" );
	
	TH2F *h2p_tagh = new TH2F( "h2p_tagh", "#DeltaK vs. TAGH Counter", 280, -0.5, 279.5, 
		2000, -4.0, 4.0 );
	TH2F *h2p_tagm = new TH2F( "h2p_tagm", "#DeltaK vs. TAGM Counter", 110, -0.5, 109.5, 
		2000, -4.0, 4.0 );
	
	TH2F *h2t_tagh = new TH2F( "h2t_tagh", "#DeltaK vs. TAGH Counter", 280, -0.5, 279.5, 
		2000, -4.0, 4.0 );
	TH2F *h2t_tagm = new TH2F( "h2t_tagm", "#DeltaK vs. TAGM Counter", 110, -0.5, 109.5, 
		2000, -4.0, 4.0 );
	
	for( int ir=1; ir<=10; ir++ ) {
		h2p_tagh->Add( (TH2F*)fPair->Get(Form("DeltaK/deltaK_tagh_pe_weight_%d",ir)) );
		h2p_tagm->Add( (TH2F*)fPair->Get(Form("DeltaK/deltaK_tagm_pe_weight_%d",ir)) );
		h2t_tagh->Add( (TH2F*)fTrip->Get(Form("DeltaK/deltaK_tagh_pe_weight_%d",ir)) );
		h2t_tagm->Add( (TH2F*)fTrip->Get(Form("DeltaK/deltaK_tagm_pe_weight_%d",ir)) );
	}
	
	
	// Compton Simulation:
	
	TFile *fComp = new TFile( 
	"/work/halld/home/andrsmit/primex_compton_analysis/compton_mc/recRootFiles/Be_sim.root", 
	"READ" );
	
	TH2F *h2c_tagh = new TH2F( "h2c_tagh", "#DeltaK vs. TAGH Counter", 274, 0.5, 274.5, 
		2000, -4.0, 4.0 );
	TH2F *h2c_tagm = new TH2F( "h2c_tagm", "#DeltaK vs. TAGM Counter", 102, 0.5, 102.5, 
		2000, -4.0, 4.0 );
	
	h2c_tagh->Add( (TH2F*)fComp->Get("DeltaK/deltaK_tagh_ep") );
	h2c_tagm->Add( (TH2F*)fComp->Get("DeltaK/deltaK_tagm_ep") );
	
	
	//----------------------------------------------------//
	
	
	
	
	
	// TAGM Counter 44:
	
	TH1F *hData = (TH1F*)h2_tagm->ProjectionY( "hData", 44, 44 );
	TH1F *hComp = (TH1F*)h2c_tagm->ProjectionY( "hComp", 44, 44 );
	TH1F *hPair = (TH1F*)h2p_tagm->ProjectionY( "hPair", 45, 45 );
	TH1F *hTrip = (TH1F*)h2t_tagm->ProjectionY( "hTrip", 45, 45 );
	
	hData->Rebin(rebins);
	hData->SetLineColor(kBlack);
	hData->SetMinimum(0.);
	hData->GetXaxis()->SetTitle( "E_{Compton} - E_{#gamma} [GeV]" );
	hData->SetTitle( Form("TAGM Counter 44 (E_{#gamma} = %.4f GeV)", tagm_en[43]) );
	
	hComp->Rebin(rebins);
	hComp->SetLineColor(kBlue);
	
	hPair->Rebin(rebins);
	hPair->SetLineColor(kGreen);
	
	hTrip->Rebin(rebins);
	hTrip->SetLineColor(kMagenta);
	
	TCanvas *canvas = new TCanvas( "canvas", "canvas", 800, 800 );
	canvas->SetTickx(); canvas->SetTicky();
	hData->Draw();
	hComp->Draw("same");
	hPair->Draw("same");
	hTrip->Draw("same");
	
	
	TObjArray *mctot = new TObjArray(3);
	mctot->Add(hComp);
	mctot->Add(hPair);
	//mctot->Add(hTrip);
	TFractionFitter* fracfit = new TFractionFitter(hData, mctot);
	
	fracfit->SetRangeX(101,300);
	Int_t status = fracfit->Fit();
	
	
	
	
	TCanvas *canvas0 = new TCanvas( "canvas0", "canvas0", 800, 800 );
	canvas0->SetTickx(); canvas0->SetTicky();
	
	TH1F *result = (TH1F*)fracfit->GetPlot();
	result->SetLineColor(kRed);
	
	TH1F *mcp0 = (TH1F*)fracfit->GetMCPrediction(0);
	mcp0->SetLineColor(kBlue);
	
	TH1F *mcp1 = (TH1F*)fracfit->GetMCPrediction(1);
	mcp1->SetLineColor(kGreen);
	
	
	hData->Draw("PE");
	result->Draw("same");
	mcp0->Draw("same");
	mcp1->Draw("same");
	
	cout << mcp0->Integral() << endl;
	cout << mcp1->Integral() << endl;
	
	
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
	
	int fit_val = 1;
	
	
	char tag_sys_char[256];
	if( tag_sys==0 ) sprintf( tag_sys_char, "tagh" );
	else             sprintf( tag_sys_char, "tagm" );
	
	
	double bkgd_min = h1->GetBinCenter( h1->FindFirstBinAbove() );
	double bkgd_max = h1->GetBinCenter(  h1->FindLastBinAbove() );
	
	TF1 *f_bkgd = new TF1( Form("f_bkgd_%s_%d",tag_sys_char,counter), 
		bkgd_fit, bkgd_min, 1.5, 4 );
	h1->Fit( Form("f_bkgd_%s_%d",tag_sys_char,counter), "R0QL" );
	
	
	
	TF1 *f_gaus = new TF1( Form("f_gaus_%s_%d",tag_sys_char,counter), "gaus", -2.0, 2.0 );
	
	f_gaus->SetParameters( h1->GetMaximum(), h1->GetBinCenter(h1->GetMaximumBin()), 0.2 );
	f_gaus->SetParLimits( 2, 0.0, 1.0 );
	
	h1->Fit( Form("f_gaus_%s_%d",tag_sys_char,counter), "R0QL" );
	
	f_gaus->SetRange( f_gaus->GetParameter(1)-0.2, f_gaus->GetParameter(1)+0.2 );
	
	h1->Fit( Form("f_gaus_%s_%d",tag_sys_char,counter), "R0QL" );
	
	
	
	TF1 *f_fit  = new TF1( Form("f_fit_%s_%d",tag_sys_char,counter), 
		double_gaus_fit, -2., 2., 10 );
	
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
	
	f_fit->FixParameter(6, f_bkgd->GetParameter(0) );
	f_fit->FixParameter(7, f_bkgd->GetParameter(1) );
	f_fit->FixParameter(8, f_bkgd->GetParameter(2) );
	f_fit->FixParameter(9, f_bkgd->GetParameter(3) );
	
	f_fit->SetRange( -2.0, 2.0 );
	
	h1->Fit( Form("f_fit_%s_%d",tag_sys_char,counter), "R0QL" );
	
	if( !FIX_BACKGROUND ) {
		f_fit->ReleaseParameter(6);
		f_fit->ReleaseParameter(7);
		f_fit->ReleaseParameter(8);
		f_fit->ReleaseParameter(9);
		
		//f_fit->SetRange(-2., 2.);
		h1->Fit( Form("f_fit_%s_%d",tag_sys_char,counter), "R0QL" );
	}
	
	f_fit->SetRange( -2., 1.5 );
	/*
	f_fit->FixParameter(6, f_fit->GetParameter(6));
	f_fit->FixParameter(7, f_fit->GetParameter(7));
	f_fit->FixParameter(8, f_fit->GetParameter(8));
	f_fit->FixParameter(9, f_fit->GetParameter(9));
	*/
	
	TFitResultPtr result = h1->Fit( Form("f_fit_%s_%d",tag_sys_char,counter), "SR0QL" );
	if( !((Int_t)result == 0) ) return 0;
	
	/*
	f_fit->SetParameter(6, f_bkgd->GetParameter(0) );
	f_fit->SetParameter(7, f_bkgd->GetParameter(1) );
	f_fit->SetParameter(8, f_bkgd->GetParameter(2) );
	f_fit->SetParameter(9, f_bkgd->GetParameter(3) );
	
	f_fit->SetRange( bkgd_min, 2.0 );
	
	h1->Fit( Form("f_fit_%s_%d",tag_sys_char,counter), "R0QL" );
	
	f_fit->FixParameter(6, f_fit->GetParameter(6));
	f_fit->FixParameter(7, f_fit->GetParameter(7));
	f_fit->FixParameter(8, f_fit->GetParameter(8));
	f_fit->FixParameter(9, f_fit->GetParameter(9));
	
	
	if( !FIX_BACKGROUND ) {
		f_fit->ReleaseParameter(6);
		f_fit->ReleaseParameter(7);
		f_fit->ReleaseParameter(8);
		f_fit->ReleaseParameter(9);
		
		h1->Fit( Form("f_fit_%s_%d",tag_sys_char,counter), "R0QL" );
	}
	
	f_fit->SetRange( -2., 1.5 );
	
	TFitResultPtr result = h1->Fit( Form("f_fit_%s_%d",tag_sys_char,counter), "SR0QL" );
	if( !((Int_t)result == 0) ) fit_val = 0;
	*/
	
	chi2 = result->Chi2() / result->Ndf();
	
	
	h1->Fit( Form("f_fit_%s_%d",tag_sys_char,counter), "R0QL" );
	h1->GetXaxis()->SetRangeUser( -2.0, 1.5 );
	
	
	
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
	
	
	TF1 *f_draw = new TF1( Form("f_draw_%d_%d",tag_sys,counter), "pol3", -3., 3. );
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
	
	pad_log->cd();
	h1->Draw();
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
	h1->GetXaxis()->SetRangeUser( -2.0, bkgd_max );
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
	h1_dev->GetXaxis()->SetRangeUser( -2.0, bkgd_max );
	h1_dev->SetTitle( " " );
	h1_dev->GetYaxis()->SetRangeUser( -12., 12. );
	
	top_pad->cd();
	h1->Draw();
	f_fit->Draw("same");
	
	bot_pad->cd();
	h1_dev->Draw("PE");
	
	//canvas2->Update();
	
	
	
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
	
	TH1F *hbeam  = (TH1F*)fSim->Get( "beam" )->Clone( Form("h_beam_%d_%d",tag_sys,counter) );
	double n_gen = 2.5e5 * (double)hbeam->GetMean();
	
	TH2F *h2;
	
	if( tag_sys == 0 ) {
		
		h2 = new TH2F( Form("h2_tagh_%d",counter), "DeltaK", 
			274, 0.5, 274.5, 2000, -4.0, 4.0 );
		
		h2->Add( (TH2F*)fSim->Get(Form("%s",hname_tagh)) );
	} else {
		
		h2 = new TH2F( Form("h2_tagm_%d",counter), "DeltaK", 
			102, 0.5, 102.5, 2000, -4.0, 4.0 );
		
		h2->Add( (TH2F*)fSim->Get(Form("%s",hname_tagm)) );
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




int fit_yield_ls( int tag_sys, int counter, TH1F *h1, double &yield, double &yieldE, double &chi2 )
{
	
	int fit_val = 1;
	
	
	char tag_sys_char[256];
	if( tag_sys==0 ) sprintf( tag_sys_char, "tagh" );
	else             sprintf( tag_sys_char, "tagm" );
	
	
	
	
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
	
	TH1F *hbeam  = (TH1F*)fSim->Get( "beam" )->Clone( Form("h_beam_%d_%d",tag_sys,counter) );
	double n_gen = 2.5e5 * (double)hbeam->GetMean();
	
	if( n_gen < 1.e3 ) return 0;
	
	TH2F *h2;
	
	if( tag_sys == 0 ) {
		
		h2 = new TH2F( Form("h2_tagh_%d",counter), "DeltaK", 
			274, 0.5, 274.5, 2000, -4.0, 4.0 );
		
		h2->Add( (TH2F*)fSim->Get(Form("%s",hname_tagh)) );
	} else {
		
		h2 = new TH2F( Form("h2_tagm_%d",counter), "DeltaK", 
			102, 0.5, 102.5, 2000, -4.0, 4.0 );
		
		h2->Add( (TH2F*)fSim->Get(Form("%s",hname_tagm)) );
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
	
	
	
	//---------------------------------------------//
	
	
	
	
	double bkgd_min = h1->GetBinCenter( h1->FindFirstBinAbove() );
	double bkgd_max = h1->GetBinCenter( h1->FindLastBinAbove() );
	
	TF1 *f_bkgd = new TF1( Form("f_bkgd_%s_%d",tag_sys_char,counter), bkgd_fit, 
		bkgd_min, bkgd_max, 4 );
	h1->Fit( Form("f_bkgd_%s_%d",tag_sys_char,counter), "R0QL" );
	
	
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
	
	
	TF1 *f_fit  = new TF1( Form("f_fit_%s_%d",tag_sys_char,counter), 
		line_shape_fit, bkgd_min, 2., 6 );
	
	f_fit->SetLineColor( kRed );
	f_fit->SetNpx( 500. );
	
	f_fit->SetParameters( scale_guess, h1_max - h1_sim_max );
	f_fit->FixParameter(1, 0.0);
	f_fit->FixParameter(2, f_bkgd->GetParameter(0));
	f_fit->FixParameter(3, f_bkgd->GetParameter(1));
	f_fit->FixParameter(4, f_bkgd->GetParameter(2));
	f_fit->FixParameter(5, f_bkgd->GetParameter(3));
	
	h1->Fit( Form("f_fit_%s_%d",tag_sys_char,counter), "R0Q" );
	
	if( !FIX_BACKGROUND ) {
		f_fit->ReleaseParameter(2);
		f_fit->ReleaseParameter(3);
		f_fit->ReleaseParameter(4);
		f_fit->ReleaseParameter(5);
		
		h1->Fit( Form("f_fit_%s_%d",tag_sys_char,counter), "R0QL" );
	}
	
	f_fit->FixParameter(2, f_fit->GetParameter(2));
	f_fit->FixParameter(3, f_fit->GetParameter(3));
	f_fit->FixParameter(4, f_fit->GetParameter(4));
	f_fit->FixParameter(5, f_fit->GetParameter(5));
	
	
	TFitResultPtr result = h1->Fit( Form("f_fit_%s_%d",tag_sys_char,counter), "SR0QL" );
	//if( !((Int_t)result == 0) ) fit_val = 0;
	
	chi2 = result->Chi2() / result->Ndf();
	
	
	
	if( CS_FROM_FIT ) {
	
	h1_sim->Scale( f_fit->GetParameter(0) );
	
	yield  = f_fit->GetParameter(0);
	yieldE = f_fit->GetParError(0);
	
	} else {
	
	TF1 *fInt     = new TF1( Form("fInt_tagh_%d",counter), "pol4", -3., 3. );
	fInt->SetParameters( 0., f_fit->GetParameter(2), 0.5*f_fit->GetParameter(3), 
		(1./3.)*f_fit->GetParameter(4), 0.25*f_fit->GetParameter(5) );
	
	double intsub = fInt->Eval(1.5) - fInt->Eval(-2.);
	intsub       /= bin_size;
	
	yield  = h1->Integral() - intsub;
	yieldE = sqrt(yield);
	
	}
	
	
	
	// Draw fits:
	
	
	TF1 *f_comp_ls = new TF1( Form("f_comp_ls_%d_%d",tag_sys,counter), line_shape_fit, 
		-4.0, 4.0, 6 );
	f_comp_ls->SetParameters( f_fit->GetParameter(0), f_fit->GetParameter(1), 
		0.0, 0.0, 0.0, 0.0 );
	f_comp_ls->SetLineColor(kBlue);
	
	
	TF1 *f_bkgd_ls = new TF1( Form("f_bkgd_ls_%d_%d",tag_sys,counter), line_shape_fit,
		-4.0, 4.0, 6 );
	f_bkgd_ls->SetParameters( 0.0, 0.0, f_fit->GetParameter(2), f_fit->GetParameter(3), 
		f_fit->GetParameter(4), f_fit->GetParameter(5) );
	f_bkgd_ls->SetLineColor(kGreen);
	f_bkgd_ls->SetLineStyle(2);
	
	
	TH1F *h1_lin = (TH1F*)h1->Clone( Form("h1_lin_%d_%d",tag_sys,counter) );
	h1_lin->SetMinimum(0.);
	
	pad_lin->cd();
	h1_lin->Draw();
	f_fit->Draw("same");
	f_comp_ls->Draw("same");
	f_bkgd_ls->Draw("same");
	
	
	pad_log->cd();
	h1->Draw();
	f_fit->Draw("same");
	f_comp_ls->Draw("same");
	f_bkgd_ls->Draw("same");
	
	//canvas1->Update();
	
	
	
	
	
	
	
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
	h1->GetXaxis()->SetRangeUser( -2.0, bkgd_max );
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
	h1_dev->GetXaxis()->SetRangeUser( -2.0, bkgd_max );
	h1_dev->SetTitle( " " );
	h1_dev->GetYaxis()->SetRangeUser( -12., 12. );
	
	top_pad->cd();
	h1->Draw();
	f_fit->Draw("same");
	
	bot_pad->cd();
	h1_dev->Draw("PE");
	
	//canvas2->Update();
	
	
	
	return fit_val;
}


