#include "/work/halld/home/andrsmit/primex_compton_analysis/include/compton_inputs.cc"
#include "/work/halld/home/andrsmit/primex_compton_analysis/include/get_phase1_energy_bin.cc"

bool FIX_BACKGROUND;
bool SUBTRACT_EMPTY;
bool SUBTRACT_PAIR;

vector<int> tagh_fitVec, tagm_fitVec;
vector<double> tagh_mu,      tagm_mu;
vector<double> tagh_muE,     tagm_muE;
vector<double> tagh_sigma,   tagm_sigma;
vector<double> tagh_sigmaE,  tagm_sigmaE;
vector<double> tagh_frac,    tagm_frac;
vector<double> tagh_fracE,   tagm_fracE;
vector<double> tagh_dmu,     tagm_dmu;
vector<double> tagh_dmuE,    tagm_dmuE;
vector<double> tagh_csigma,  tagm_csigma;
vector<double> tagh_csigmaE, tagm_csigmaE;
vector<double> tagh_chi2,    tagm_chi2;

vector<int> tagh_fitVec_mc, tagm_fitVec_mc;
vector<double> tagh_mu_mc,      tagm_mu_mc;
vector<double> tagh_muE_mc,     tagm_muE_mc;
vector<double> tagh_sigma_mc,   tagm_sigma_mc;
vector<double> tagh_sigmaE_mc,  tagm_sigmaE_mc;
vector<double> tagh_frac_mc,    tagm_frac_mc;
vector<double> tagh_fracE_mc,   tagm_fracE_mc;
vector<double> tagh_dmu_mc,     tagm_dmu_mc;
vector<double> tagh_dmuE_mc,    tagm_dmuE_mc;
vector<double> tagh_csigma_mc,  tagm_csigma_mc;
vector<double> tagh_csigmaE_mc, tagm_csigmaE_mc;
vector<double> tagh_chi2_mc,    tagm_chi2_mc;

//----------   Function Declarations   ----------//

TF1 *f_frac_mc, *f_dmu_mc, *f_csigma_mc;

Double_t     bkgd_fit(Double_t *x, Double_t *par);
Double_t dbl_gaus_fit(Double_t *x, Double_t *par);

TH1F *h_bkgd_fit;
Double_t bkgd_fit_lineshape(Double_t *x, Double_t *par);

void fit_mc(bool draw_fits=false);
void fit_data(bool draw_fits=false);

int deltaK_fit(int is_mc, int tag_sys, int tag_counter, TH1F *h1, 
	double &mu,     double &muE, 
	double &sigma,  double &sigmaE,
	double &frac,   double &fracE,
	double &dmu,    double &dmuE,
	double &csigma, double &csigmaE, 
	double &chi2,  bool draw_fits = false);

void get_tgrapherrors_mc(TGraphErrors **g, 
	vector<double> tagh_vec, vector<double> taghE_vec, 
	vector<double> tagm_vec, vector<double> tagmE_vec);
void get_tgrapherrors(TGraphErrors **g, 
	vector<double> tagh_vec, vector<double> taghE_vec, 
	vector<double> tagm_vec, vector<double> tagmE_vec);

void get_tgraph_mc(TGraph **g, vector<double> tagh_vec, vector<double> tagm_vec);
void get_tgraph(TGraph **g, vector<double> tagh_vec, vector<double> tagm_vec);

int get_data_hists();
int get_pair_mc_hists(TString hists_fname, TString flux_fname);

//-----------------------------------------------//

TCanvas *canvas1;
TPad *pad_lin, *pad_log;

TCanvas *canvas2;
TPad *top_pad, *bot_pad;

TCanvas *canvas3;

int cut_index;
TString hist_name;

TH2F *h2_tagh,       *h2_tagm;
TH2F *h2_tagh_empty, *h2_tagm_empty;
TH2F *h2_tagh_pair,  *h2_tagm_pair;
TH1F *h_pair_gen_flux;

void fit_deltaK(int loc_cut_index=0)
{
	gStyle->SetOptStat(0); gStyle->SetOptFit(0);
	
	//----------   Initialize   ----------//
	
	// Adjustable switches:
	
	const bool DRAW_FITS_TAGH = false;
	const bool DRAW_FITS_TAGM = false;
	
	DEBUG_FITS     = true;
	
	FIX_BACKGROUND = false;
	SUBTRACT_EMPTY = true;
	SUBTRACT_PAIR  = true;
	
	//----------   Set up path names   ----------//
	
	const char loc_pathname[256] = "/work/halld/home/andrsmit/primex_compton_analysis";
	
	if(false) {
		
		// Beryllium Target
		
		IS_BE_TARGET = true;
		
		TARGET_STR = "Be";
		
		endpoint_energy_calib = 11.6061;
		endpoint_energy       = 11.608;
		
		BEAM_CURRENT = 200;
		
		// Directory with histograms from Compton MC:
		comp_mc_dir = Form("%s/analyze_trees/analyze_mc/rootFiles/phase1/Run061321/compton/%03dnA", loc_pathname, BEAM_CURRENT);
		
	} else {
		
		// Helium Target
		
		IS_BE_TARGET = false;
		
		endpoint_energy_calib = 11.1689;
		
		/*
		BEAM_CURRENT    = 200;
		endpoint_energy = 11.1666;
		
		BEAM_CURRENT    = 100;
		endpoint_energy = 11.1664;
		
		BEAM_CURRENT    = 50;
		endpoint_energy = 11.1671;
		*/
		
		BEAM_CURRENT    = 200;
		endpoint_energy = 11.1666;
		
		TARGET_STR = "He";
		
		// Directory with histograms from Compton MC:
		
		comp_mc_dir = Form("%s/analyze_trees/analyze_mc/rootFiles/phase1/Run061866/compton/%03dnA", loc_pathname, BEAM_CURRENT);
	}
	
	// ROOT filenames for full and empty target data:
	
	root_fname              = Form("%s/analyze_trees/analyze_data/rootFiles/phase1/%s_%03dnA_FIELDOFF.root", 
		loc_pathname, TARGET_STR.Data(), BEAM_CURRENT);
	empty_target_root_fname = Form("%s/analyze_trees/analyze_data/rootFiles/phase1/%s_empty_FIELDOFF.root", 
		loc_pathname, TARGET_STR.Data());
	
	// Filled target flux filenames:
	
	tagh_flux_fname = Form("%s/photon_flux/phase1/%s_%03dnA_FIELDOFF_flux_tagh.txt", 
		loc_pathname, TARGET_STR.Data(), BEAM_CURRENT);
	tagm_flux_fname = Form("%s/photon_flux/phase1/%s_%03dnA_FIELDOFF_flux_tagm.txt", 
		loc_pathname, TARGET_STR.Data(), BEAM_CURRENT);
	
	// Empty target flux filenames:
	
	empty_target_tagh_flux_fname = Form("%s/photon_flux/phase1/Be_empty_FIELDOFF_flux_tagh.txt", loc_pathname);
	empty_target_tagm_flux_fname = Form("%s/photon_flux/phase1/Be_empty_FIELDOFF_flux_tagm.txt", loc_pathname);
	
	// file containing the xscales for the tagh and tagm counters:
	
	tagh_xscale_fname = Form("%s/photon_flux/phase1/primex_tagh.txt", loc_pathname);
	tagm_xscale_fname = Form("%s/photon_flux/phase1/primex_tagm.txt", loc_pathname);
	
	//------------------------------------------------//
	// Get flux, tagger energies, and e+e- cross section (from NIST):
	
	get_counter_energies();
	get_flux();
	get_target_parameters();
	
	DRAW_THEORY = true;
	pair_cs_fname = Form("%s/photon_absorption/%s_pair_cs.dat", loc_pathname, TARGET_STR.Data());
	get_pair_cs();
	
	cut_index = loc_cut_index;
	hist_name = Form("deltaK/deltaK");
	
	//------------------------------------------------//
	// Get histograms from data:
	
	if(get_data_hists()) return;
	
	//------------------------------------------------//
	// Get e+e- distributions:
	
	TString pair_mc_hist_fname = Form("%s/analyze_trees/analyze_mc/rootFiles/phase1/Run061321/pair/pair_rec_combined.root", 
		loc_pathname);
	TString pair_mc_flux_fname = Form("%s/bhgen_test/recRootTrees/Run061321/sum.root", loc_pathname);
	if(get_pair_mc_hists(pair_mc_hist_fname, pair_mc_flux_fname)) return;
	
	//------------------------------------------------//
	// Do fits to the Compton MC:
	
	rebins = 2;
	
	bool DRAW_FITS_MC = false;
	cout << "Fitting MC..." << flush;
	fit_mc(DRAW_FITS_MC);
	cout << "done." << endl;
	
	TGraphErrors *g_mu_mc, *g_sigma_mc, *g_frac_mc, *g_dmu_mc, *g_csigma_mc;
	get_tgrapherrors_mc(&g_mu_mc,     tagh_mu_mc,     tagh_muE_mc,     tagm_mu_mc,     tagm_muE_mc);
	get_tgrapherrors_mc(&g_sigma_mc,  tagh_sigma_mc,  tagh_sigmaE_mc,  tagm_sigma_mc,  tagm_sigmaE_mc);
	get_tgrapherrors_mc(&g_frac_mc,   tagh_frac_mc,   tagh_fracE_mc,   tagm_frac_mc,   tagm_fracE_mc);
	get_tgrapherrors_mc(&g_dmu_mc,    tagh_dmu_mc,    tagh_dmuE_mc,    tagm_dmu_mc,    tagm_dmuE_mc);
	get_tgrapherrors_mc(&g_csigma_mc, tagh_csigma_mc, tagh_csigmaE_mc, tagm_csigma_mc, tagm_csigmaE_mc);
	
	TGraph *g_chi2_mc;
	get_tgraph_mc(&g_chi2_mc, tagh_chi2_mc, tagm_chi2_mc);
	
	f_frac_mc = new TF1("f_frac_mc", "pol3", 6., 11.5);
	g_frac_mc->Fit(f_frac_mc, "R0Q");
	
	f_dmu_mc = new TF1("f_dmu_mc", "pol3", 6., 11.5);
	g_dmu_mc->Fit(f_dmu_mc, "R0Q");
	
	f_csigma_mc = new TF1("f_csigma_mc", "pol3", 6., 11.5);
	g_csigma_mc->Fit(f_csigma_mc, "R0Q");
	
	//------------------------------------------------//
	// Do fits to the real data:
	
	bool DRAW_FITS_DATA = false;
	cout << "Fitting data..." << flush;
	fit_data(DRAW_FITS_DATA);
	cout << "done." << endl;
	
	TGraphErrors *g_mu, *g_sigma, *g_frac, *g_dmu, *g_csigma;
	get_tgrapherrors(&g_mu,     tagh_mu,     tagh_muE,     tagm_mu,     tagm_muE);
	get_tgrapherrors(&g_sigma,  tagh_sigma,  tagh_sigmaE,  tagm_sigma,  tagm_sigmaE);
	get_tgrapherrors(&g_frac,   tagh_frac,   tagh_fracE,   tagm_frac,   tagm_fracE);
	get_tgrapherrors(&g_dmu,    tagh_dmu,    tagh_dmuE,    tagm_dmu,    tagm_dmuE);
	get_tgrapherrors(&g_csigma, tagh_csigma, tagh_csigmaE, tagm_csigma, tagm_csigmaE);
	
	TGraph *g_chi2;
	get_tgraph(&g_chi2, tagh_chi2, tagm_chi2);
	
	//------------------------------------------------//
	// Plot the resulting fit parameters:
	
	// Mu:
	g_mu->GetYaxis()->SetRangeUser(-1.25,0.10);
	g_mu->SetTitle("#mu from Di-Gaussian Fit");
	g_mu->GetYaxis()->SetTitle("#mu_{#DeltaK} [GeV]");
	
	TCanvas *c_mu = new TCanvas("c_mu", "Mu", 1000, 600);
	c_mu->SetTickx(); c_mu->SetTicky();
	g_mu->Draw("AP");
	g_mu_mc->Draw("P same");
	
	cout << "\nDeltaK mu fit results for data:" << endl;
	TF1 *f_mu = new TF1("f_mu", "[0] + [1]*x + [2]*x*x + [3]*x*x*x", 6.0, 11.1);
	if(loc_cut_index%4==0 || loc_cut_index%4==1) f_mu->FixParameter(3, 0.);
	g_mu->Fit(f_mu, "R0");
	f_mu->SetLineColor(kBlue);
	f_mu->Draw("same");
	
	cout << "\nDeltaK mu fit results for mc:" << endl;
	TF1 *f_mu_mc = new TF1("f_mu_mc", "[0] + [1]*x + [2]*x*x + [3]*x*x*x", 6.0, 11.1);
	if(loc_cut_index%4==0 || loc_cut_index%4==1) f_mu_mc->FixParameter(3, 0.);
	g_mu_mc->Fit(f_mu_mc, "R0");
	f_mu_mc->SetLineColor(kGreen);
	f_mu_mc->Draw("same");
	
	
	// Sigma:
	g_sigma->GetYaxis()->SetRangeUser(0.5,1.1);
	g_sigma->SetTitle("#sigma from Di-Gaussian Fit");
	g_sigma->GetYaxis()->SetTitle("#sigma_{#DeltaK} [GeV]");
	
	TCanvas *c_sigma = new TCanvas("c_sigma", "Sigma", 1000, 600);
	c_sigma->SetTickx(); c_sigma->SetTicky();
	g_sigma->Draw("AP");
	g_sigma_mc->Draw("P same");
	
	cout << "\nDeltaK sigma fit results for data:" << endl;
	TF1 *f_sigma = new TF1("f_sigma", "[0] + [1]*x + [2]*x^2 + [3]*x^3", 6.0, 11.4);
	g_sigma->Fit(f_sigma, "R0");
	f_sigma->SetLineColor(kBlue);
	f_sigma->Draw("same");
	
	cout << "\nDeltaK sigma fit results for mc:" << endl;
	TF1 *f_sigma_mc = new TF1("f_sigma_mc", "[0] + [1]*x + [2]*x^2 + [3]*x^3", 6.0, 11.4);
	g_sigma_mc->Fit(f_sigma_mc, "R0");
	f_sigma_mc->SetLineColor(kGreen);
	f_sigma_mc->Draw("same");
	
	
	// Fraction of second Gaussian:
	g_frac->GetYaxis()->SetRangeUser(0., 1.0);
	g_frac->SetTitle("Fraction of Second Gaussian");
	g_frac->GetYaxis()->SetTitle("frac");
	
	TCanvas *c_frac = new TCanvas("c_frac", "Fraction", 1000, 600);
	c_frac->SetTickx(); c_frac->SetTicky();
	g_frac->Draw("AP");
	g_frac_mc->Draw("P same");
	
	
	// DeltaMu of second Gaussian:
	g_dmu->GetYaxis()->SetRangeUser(0., 1.5);
	g_dmu->SetTitle("#Delta#mu of Second Gaussian");
	g_dmu->GetYaxis()->SetTitle("#Delta#mu [GeV]");
	
	TCanvas *c_dmu = new TCanvas("c_dmu", "DeltaMu", 1000, 600);
	c_dmu->SetTickx(); c_dmu->SetTicky();
	g_dmu->Draw("AP");
	g_dmu_mc->Draw("P same");
	
	
	// cSigma of second Gaussian:
	g_csigma->GetYaxis()->SetRangeUser(0., 3.5);
	g_csigma->SetTitle("#sigma ratio from Di-Gaussian Fit");
	g_csigma->GetYaxis()->SetTitle("#sigma_{2} / #sigma_{1}");
	
	TCanvas *c_csigma = new TCanvas("c_csigma", "cSigma", 1000, 600);
	c_csigma->SetTickx(); c_csigma->SetTicky();
	g_csigma->Draw("AP");
	g_csigma_mc->Draw("P same");
	
	
	// Chi2/Ndf:
	g_chi2->GetYaxis()->SetRangeUser(0., 5.0);
	g_chi2->SetTitle("#chi^{2}/N.d.f from Crystal Ball Fit");
	g_chi2->GetYaxis()->SetTitle("#chi^{2}/N.d.f.");
	
	TCanvas *c_chi2 = new TCanvas("c_chi2", "c_chi2", 1000, 600);
	c_chi2->SetTickx(); c_chi2->SetTicky();
	g_chi2->Draw("AP");
	g_chi2_mc->Draw("P same");
	
	
	//------------------------------------------------------------------------------------------//
	
	TF1 *f_mu_compare = new TF1("f_mu_compare", "pol3", 5.0, 12.0);
	f_mu_compare->SetParameters(2.377e-02, -8.578e-03, 4.108e-04);
	f_mu_compare->SetParameters(f_mu->GetParameters());
	f_mu_compare->SetLineColor(kBlack);
	f_mu_compare->SetLineStyle(2);
	
	TF1 *f_sigma_compare = new TF1("f_sigma_compare", "pol3", 5.0, 12.0);
	f_sigma_compare->SetParameters(9.016e-03, 4.514e-02); // Be_200nA
	f_sigma_compare->SetParameters(f_sigma->GetParameters());
	f_sigma_compare->SetLineColor(kBlack);
	f_sigma_compare->SetLineStyle(2);
	
	TCanvas *c_deltaK = new TCanvas("c_deltaK", "c_deltaK", 1200, 800);
	TPad *pMu  = new TPad("pMu",  "pMu",  0.0, 0.5, 1.0, 1.0);
	TPad *pSig = new TPad("pSig", "pSig", 0.0, 0.0, 1.0, 0.5);
	
	pMu->SetTickx(); pMu->SetTicky();
	pMu->SetTopMargin(0.02);
	pMu->SetBottomMargin(0.11);
	pMu->SetLeftMargin(0.07);
	pMu->SetRightMargin(0.05);
	
	pSig->SetTickx(); pSig->SetTicky();
	pSig->SetBottomMargin(0.11);
	pSig->SetTopMargin(0.02);
	pSig->SetLeftMargin(0.07);
	pSig->SetRightMargin(0.05);
	
	c_deltaK->cd();
	pMu->Draw();
	pSig->Draw();
	
	TLatex lat;
	lat.SetTextFont(52);
	lat.SetTextColor(kRed);
	
	g_mu->SetTitle("");
	g_mu->GetXaxis()->SetRangeUser(5.8, 11.6);
	g_mu->GetYaxis()->SetRangeUser(-0.5, 0.1);
	g_mu->SetMarkerColor(kGreen);
	g_mu->SetLineColor(kGreen);
	g_mu->GetYaxis()->CenterTitle(true);
	g_mu->GetYaxis()->SetTitleOffset(0.4);
	
	g_sigma->SetTitle("");
	g_sigma->GetXaxis()->SetRangeUser(5.8, 11.6);
	g_sigma->GetYaxis()->SetRangeUser(0.4, 1.0);
	g_sigma->SetMarkerColor(kGreen);
	g_sigma->SetLineColor(kGreen);
	g_sigma->GetYaxis()->CenterTitle(true);
	g_sigma->GetYaxis()->SetTitleOffset(0.4);
	
	pMu->cd();
	g_mu->Draw("AP");
	f_mu_compare->Draw("same");
	
	lat.DrawLatexNDC(0.1035, 0.8865, Form("#scale[1.5]{%s Target Data}", TARGET_STR.Data()));
	
	TPaveText *ptMu = new TPaveText(0.424, 0.173, 0.595, 0.418, "NDC");
	ptMu->SetTextAlign(12);
	ptMu->SetBorderSize(0);
	ptMu->SetFillColor(0);
	ptMu->SetTextFont(52);
	ptMu->AddText(Form("p_{0} =  %.3e",f_mu_compare->GetParameter(0)));
	ptMu->AddText(Form("p_{1} = %.3e", f_mu_compare->GetParameter(1)));
	ptMu->AddText(Form("p_{2} =  %.3e",f_mu_compare->GetParameter(2)));
	//ptMu->AddText(Form("p_{3} =  %.3e",f_mu_compare->GetParameter(3)));
	ptMu->Draw();
	
	TPaveText *ptMuEqn = new TPaveText(0.091, 0.259, 0.418, 0.329, "NDC");
	ptMuEqn->SetBorderSize(0);
	ptMuEqn->SetFillColor(0);
	ptMuEqn->SetTextFont(52);
	ptMuEqn->SetTextColor(kRed+2);
	ptMuEqn->AddText("#mu_{#DeltaK}#left(E_{#gamma}#right) = p_{0} + p_{1}E_{#gamma} + p_{2}E_{#gamma}^{2}"); 
	ptMuEqn->Draw();
	
	pSig->cd();
	g_sigma->Draw("AP");
	f_sigma_compare->Draw("same");
	
	lat.DrawLatexNDC(0.1035, 0.8865, Form("#scale[1.5]{%s Target Data}", TARGET_STR.Data()));
	
	TPaveText *ptSig = new TPaveText(0.4282, 0.1688, 0.704, 0.4240, "NDC");
	ptSig->SetTextAlign(12);
	ptSig->SetBorderSize(0);
	ptSig->SetFillColor(0);
	ptSig->SetTextFont(52);
	ptSig->AddText(Form("p_{0} = %.3e", f_sigma_compare->GetParameter(0)));
	ptSig->AddText(Form("p_{1} = %.3e", f_sigma_compare->GetParameter(1)));
	ptSig->AddText(Form("p_{2} = %.3e", f_sigma_compare->GetParameter(2)));
	ptSig->AddText(Form("p_{3} =%.3e", f_sigma_compare->GetParameter(3)));
	ptSig->Draw();
	
	TPaveText *ptSigEqn = new TPaveText(0.091, 0.259, 0.418, 0.329, "NDC");
	ptSigEqn->SetBorderSize(0);
	ptSigEqn->SetFillColor(0);
	ptSigEqn->SetTextFont(52);
	ptSigEqn->SetTextColor(kRed+2);
	ptSigEqn->AddText("#sigma_{#DeltaK}#left(E_{#gamma}#right) = p_{0} + p_{1}E_{#gamma} + p_{2}E_{#gamma}^{2} + p_{3}E_{#gamma}^{3}"); 
	ptSigEqn->Draw();
	
	
	
	return;
}

int get_data_hists() {
	
	if(gSystem->AccessPathName(root_fname.Data())) {
		cout << "Specified ROOT file for filled target data does not exist." << endl;
		return 1;
	}
	else if(gSystem->AccessPathName(empty_target_root_fname.Data())) {
		cout << "Specified ROOT file for empty target data does not exist." << endl;
		return 1;
	}
	
	TFile *fFull  = new TFile(root_fname.Data(),             "READ");
	TFile *fEmpty = new TFile(empty_target_root_fname.Data(),"READ");
	
	h2_tagh       = (TH2F*)fFull->Get(Form("%s_tagh_%d",hist_name.Data(),cut_index))->Clone("h2_tagh_data");
	h2_tagm       = (TH2F*)fFull->Get(Form("%s_tagm_%d",hist_name.Data(),cut_index))->Clone("h2_tagm_data");
	h2_tagh->SetDirectory(0);
	h2_tagm->SetDirectory(0);
	
	h2_tagh_empty = (TH2F*)fEmpty->Get(Form("%s_tagh_%d",hist_name.Data(),cut_index))->Clone("h2_tagh_empty");
	h2_tagm_empty = (TH2F*)fEmpty->Get(Form("%s_tagm_%d",hist_name.Data(),cut_index))->Clone("h2_tagm_empty");
	h2_tagh_empty->SetDirectory(0);
	h2_tagm_empty->SetDirectory(0);
	
	fFull->Close(); fEmpty->Close();
	
	return 0;
}

int get_pair_mc_hists(TString hist_fname, TString flux_fname) {
	
	if(gSystem->AccessPathName(hist_fname.Data())) {
		cout << "Specified ROOT file for e+e- pair simulation does not exist." << endl;
		return 1;
	}
	else if(gSystem->AccessPathName(flux_fname.Data())) {
		cout << "Specified ROOT file for generated e+e- pair flux does not exist." << endl;
		return 1;
	}
	
	TFile *f_pair = new TFile(hist_fname.Data(), "READ");
	h2_tagh_pair = (TH2F*)f_pair->Get(Form("%s_tagh_%d",hist_name.Data(),cut_index))->Clone("h2_tagh_pair");
	h2_tagm_pair = (TH2F*)f_pair->Get(Form("%s_tagm_%d",hist_name.Data(),cut_index))->Clone("h2_tagm_pair");
	h2_tagh_pair->SetDirectory(0);
	h2_tagm_pair->SetDirectory(0);
	f_pair->Close();
	
	TFile *f_pair_gen = new TFile(flux_fname.Data(), "READ");
	h_pair_gen_flux = (TH1F*)f_pair_gen->Get("gen_flux");
	h_pair_gen_flux->SetDirectory(0);
	f_pair_gen->Close();
	
	return 0;
}


void get_tgrapherrors_mc(TGraphErrors **g, vector<double> tagh_vec, vector<double> taghE_vec, 
	vector<double> tagm_vec, vector<double> tagmE_vec) {
	
	int n_bins_tagh = (int)tagh_fitVec_mc.size();
	int n_bins_tagm = (int)tagm_fitVec_mc.size();
	int n_bins      = n_bins_tagh + n_bins_tagm;
	
	double *loc_x     = new double[n_bins];
	double *loc_y     = new double[n_bins];
	double *loc_x_err = new double[n_bins];
	double *loc_y_err = new double[n_bins];
	
	for(int ib=0; ib<n_bins; ib++) {
		if(ib<n_bins_tagh) {
			loc_x[ib]     = tagh_en[tagh_fitVec_mc[ib]-1];
			loc_y[ib]     = tagh_vec[ib];
			loc_x_err[ib] = 0.;
			loc_y_err[ib] = taghE_vec[ib];
		} else {
			loc_x[ib]     = tagm_en[tagm_fitVec_mc[ib-n_bins_tagh]-1];
			loc_y[ib]     = tagm_vec[ib-n_bins_tagh];
			loc_x_err[ib] = 0.;
			loc_y_err[ib] = tagmE_vec[ib-n_bins_tagh];
		}
	}
	
	*g = new TGraphErrors(n_bins, loc_x, loc_y, loc_x_err, loc_y_err);
	(*g)->GetXaxis()->SetTitle("Photon Beam Energy (GeV)");
	(*g)->GetXaxis()->SetTitleSize(0.05);
	(*g)->GetXaxis()->SetTitleOffset(1.0);
	(*g)->GetXaxis()->CenterTitle(true);
	(*g)->GetYaxis()->SetTitleSize(0.05);
	(*g)->GetYaxis()->SetTitleOffset(1.0);
	(*g)->GetYaxis()->SetTitleOffset(1.0);
	(*g)->GetYaxis()->CenterTitle(true);
	(*g)->SetMarkerSize(0.6);
	(*g)->SetMarkerColor(kRed);
	(*g)->SetMarkerStyle(8);
	(*g)->SetLineColor(kRed);
	
	return;
}

void get_tgrapherrors(TGraphErrors **g, vector<double> tagh_vec, vector<double> taghE_vec, 
	vector<double> tagm_vec, vector<double> tagmE_vec) {
	
	int n_bins_tagh = (int)tagh_fitVec.size();
	int n_bins_tagm = (int)tagm_fitVec.size();
	int n_bins      = n_bins_tagh + n_bins_tagm;
	
	double *loc_x     = new double[n_bins];
	double *loc_y     = new double[n_bins];
	double *loc_x_err = new double[n_bins];
	double *loc_y_err = new double[n_bins];
	
	for(int ib=0; ib<n_bins; ib++) {
		if(ib<n_bins_tagh) {
			loc_x[ib]     = tagh_en[tagh_fitVec[ib]-1];
			loc_y[ib]     = tagh_vec[ib];
			loc_x_err[ib] = 0.;
			loc_y_err[ib] = taghE_vec[ib];
		} else {
			loc_x[ib]     = tagm_en[tagm_fitVec[ib-n_bins_tagh]-1];
			loc_y[ib]     = tagm_vec[ib-n_bins_tagh];
			loc_x_err[ib] = 0.;
			loc_y_err[ib] = tagmE_vec[ib-n_bins_tagh];
		}
	}
	
	*g = new TGraphErrors(n_bins, loc_x, loc_y, loc_x_err, loc_y_err);
	(*g)->GetXaxis()->SetTitle("Photon Beam Energy (GeV)");
	(*g)->GetXaxis()->SetTitleSize(0.05);
	(*g)->GetXaxis()->SetTitleOffset(1.0);
	(*g)->GetYaxis()->SetTitleSize(0.065);
	(*g)->GetYaxis()->SetTitleOffset(1.0);
	(*g)->GetYaxis()->SetTitleOffset(1.0);
	(*g)->SetMarkerSize(0.6);
	(*g)->SetMarkerColor(kBlue);
	(*g)->SetMarkerStyle(8);
	(*g)->SetLineColor(kBlue);
	
	return;
}

void get_tgraph_mc(TGraph **g, vector<double> tagh_vec, vector<double> tagm_vec) {
	
	int n_bins_tagh = (int)tagh_fitVec_mc.size();
	int n_bins_tagm = (int)tagm_fitVec_mc.size();
	int n_bins      = n_bins_tagh + n_bins_tagm;
	
	double *loc_x = new double[n_bins];
	double *loc_y = new double[n_bins];
	
	for(int ib=0; ib<n_bins; ib++) {
		if(ib<n_bins_tagh) {
			loc_x[ib]     = tagh_en[tagh_fitVec_mc[ib]-1];
			loc_y[ib]     = tagh_vec[ib];
		} else {
			loc_x[ib]     = tagm_en[tagm_fitVec_mc[ib-n_bins_tagh]-1];
			loc_y[ib]     = tagm_vec[ib-n_bins_tagh];
		}
	}
	
	*g = new TGraph(n_bins, loc_x, loc_y);
	(*g)->GetXaxis()->SetTitle("Photon Beam Energy (GeV)");
	(*g)->GetXaxis()->SetTitleSize(0.05);
	(*g)->GetXaxis()->SetTitleOffset(1.0);
	(*g)->GetXaxis()->CenterTitle(true);
	(*g)->GetYaxis()->SetTitleSize(0.05);
	(*g)->GetYaxis()->SetTitleOffset(1.0);
	(*g)->GetYaxis()->SetTitleOffset(1.0);
	(*g)->GetYaxis()->CenterTitle(true);
	(*g)->SetMarkerSize(0.6);
	(*g)->SetMarkerColor(kRed);
	(*g)->SetMarkerStyle(8);
	(*g)->SetLineColor(kRed);
	
	return;
}

void get_tgraph(TGraph **g, vector<double> tagh_vec, vector<double> tagm_vec) {
	
	int n_bins_tagh = (int)tagh_fitVec.size();
	int n_bins_tagm = (int)tagm_fitVec.size();
	int n_bins      = n_bins_tagh + n_bins_tagm;
	
	double *loc_x = new double[n_bins];
	double *loc_y = new double[n_bins];
	
	for(int ib=0; ib<n_bins; ib++) {
		if(ib<n_bins_tagh) {
			loc_x[ib]     = tagh_en[tagh_fitVec[ib]-1];
			loc_y[ib]     = tagh_vec[ib];
		} else {
			loc_x[ib]     = tagm_en[tagm_fitVec[ib-n_bins_tagh]-1];
			loc_y[ib]     = tagm_vec[ib-n_bins_tagh];
		}
	}
	
	*g = new TGraph(n_bins, loc_x, loc_y);
	(*g)->GetXaxis()->SetTitle("Photon Beam Energy (GeV)");
	(*g)->GetXaxis()->SetTitleSize(0.05);
	(*g)->GetXaxis()->SetTitleOffset(1.0);
	(*g)->GetYaxis()->SetTitleSize(0.05);
	(*g)->GetYaxis()->SetTitleOffset(1.0);
	(*g)->GetYaxis()->SetTitleOffset(1.0);
	(*g)->SetMarkerSize(0.6);
	(*g)->SetMarkerColor(kBlue);
	(*g)->SetMarkerStyle(8);
	(*g)->SetLineColor(kBlue);
	
	return;
}


void fit_mc(bool draw_fits=false) {
	
	char loc_fname[256];
	
	// Loop over all TAGH counters:
	for(int tagh_counter=1; tagh_counter<=221; tagh_counter++) {
		
		sprintf(loc_fname, "%s/tagh_%03d.root", comp_mc_dir.Data(), tagh_counter);
		if(gSystem->AccessPathName(loc_fname)) continue;
		
		TFile *loc_fIn = new TFile(loc_fname, "READ");
		TH2F *loc_h2 = (TH2F*)loc_fIn->Get(Form("%s_tagh_%d",hist_name.Data(),cut_index));
		TH1F *loc_h1 = (TH1F*)loc_h2->ProjectionY("loc_h1",tagh_counter,tagh_counter);
		loc_h1->SetDirectory(0);
		loc_fIn->Close();
		
		if(loc_h1->Integral() < 1.e3) continue;
		
		double loc_eb = tagh_en[tagh_counter-1];
		
		if(loc_eb>8.770 && loc_eb<8.780) continue;
		if(loc_eb>8.903 && loc_eb<8.905) continue;
		if(loc_eb>9.120 && loc_eb<8.125) continue;
		
		loc_h1->SetTitle(Form("TAGH Counter %d (E_{#gamma} = %.3f GeV)", tagh_counter, loc_eb));
		loc_h1->Rebin(rebins);
		loc_h1->SetLineColor(kBlack);
		loc_h1->GetXaxis()->SetTitle("E_{1} + E_{2} - E_{#gamma} [GeV]");
		loc_h1->SetLineWidth(2);
		
		double loc_mu    = 0, loc_sigma  = 0, loc_frac  = 0., loc_dmu  = 0., loc_csigma  = 0.;
		double loc_muE   = 0, loc_sigmaE = 0, loc_fracE = 0., loc_dmuE = 0., loc_csigmaE = 0.;
		double loc_chi2  = 0;
		
		//cout << "\n";
		//cout << "Fitting DeltaE distribution in TAGH counter " << tagh_counter << endl;
		
		int fit_val = deltaK_fit(1, 0, tagh_counter, loc_h1, 
			loc_mu, loc_muE, loc_sigma, loc_sigmaE, loc_frac, loc_fracE, loc_dmu, loc_dmuE, loc_csigma, loc_csigmaE, 
			loc_chi2, draw_fits);
		
		if(fit_val  <= 0 ) continue;
		if(loc_chi2 > 50.) continue;
		
		tagh_mu_mc.push_back(loc_mu);
		tagh_muE_mc.push_back(loc_muE);
		tagh_sigma_mc.push_back(loc_sigma);
		tagh_sigmaE_mc.push_back(loc_sigmaE);
		tagh_frac_mc.push_back(loc_frac);
		tagh_fracE_mc.push_back(loc_fracE);
		tagh_dmu_mc.push_back(loc_dmu);
		tagh_dmuE_mc.push_back(loc_dmuE);
		tagh_csigma_mc.push_back(loc_csigma);
		tagh_csigmaE_mc.push_back(loc_csigmaE);
		tagh_chi2_mc.push_back(loc_chi2);
		
		tagh_fitVec_mc.push_back(tagh_counter);
	}
	
	// Do the same for the TAGM counters:
	for(int tagm_counter=1; tagm_counter<=102; tagm_counter++) {
		
		sprintf(loc_fname, "%s/tagm_%03d.root", comp_mc_dir.Data(), tagm_counter);
		if(gSystem->AccessPathName(loc_fname)) continue;
		
		TFile *loc_fIn = new TFile(loc_fname, "READ");
		TH2F *loc_h2 = (TH2F*)loc_fIn->Get(Form("%s_tagm_%d",hist_name.Data(),cut_index));
		TH1F *loc_h1 = (TH1F*)loc_h2->ProjectionY("loc_h1",tagm_counter,tagm_counter);
		loc_h1->SetDirectory(0);
		loc_fIn->Close();
		
		if(loc_h1->Integral() < 1.e3) continue;
		
		double loc_eb = tagm_en[tagm_counter-1];
		
		if(loc_eb>8.770 && loc_eb<8.780) continue;
		if(loc_eb>8.903 && loc_eb<8.905) continue;
		if(loc_eb>9.120 && loc_eb<8.125) continue;
		
		loc_h1->SetTitle(Form("TAGM Counter %d (E_{#gamma} = %.3f GeV)", tagm_counter, loc_eb));
		loc_h1->Rebin(rebins);
		loc_h1->SetLineColor(kBlack);
		loc_h1->GetXaxis()->SetTitle("E_{1} + E_{2} - E_{#gamma} [GeV]");
		loc_h1->SetLineWidth(2);
		
		double loc_mu    = 0, loc_sigma  = 0, loc_frac  = 0., loc_dmu  = 0., loc_csigma  = 0.;
		double loc_muE   = 0, loc_sigmaE = 0, loc_fracE = 0., loc_dmuE = 0., loc_csigmaE = 0.;
		double loc_chi2  = 0;
		
		//cout << "\n";
		//cout << "Fitting DeltaE distribution in TAGM counter " << tagm_counter << endl;
		
		int fit_val = deltaK_fit(1, 1, tagm_counter, loc_h1, 
			loc_mu, loc_muE, loc_sigma, loc_sigmaE, loc_frac, loc_fracE, loc_dmu, loc_dmuE, loc_csigma, loc_csigmaE, 
			loc_chi2, draw_fits);
		
		if(fit_val  <= 0 ) continue;
		if(loc_chi2 > 50.) continue;
		
		tagm_mu_mc.push_back(loc_mu);
		tagm_muE_mc.push_back(loc_muE);
		tagm_sigma_mc.push_back(loc_sigma);
		tagm_sigmaE_mc.push_back(loc_sigmaE);
		tagm_frac_mc.push_back(loc_frac);
		tagm_fracE_mc.push_back(loc_fracE);
		tagm_dmu_mc.push_back(loc_dmu);
		tagm_dmuE_mc.push_back(loc_dmuE);
		tagm_csigma_mc.push_back(loc_csigma);
		tagm_csigmaE_mc.push_back(loc_csigmaE);
		tagm_chi2_mc.push_back(loc_chi2);
		
		tagm_fitVec_mc.push_back(tagm_counter);
	}
	
	return;
}


void fit_data(bool draw_fits=false) {
	
	char loc_fname[256];
	
	// Loop over all TAGH counters:
	for(int tagh_counter=5; tagh_counter<=221; tagh_counter++) {
		
		TH1F *loc_h1 = (TH1F*)h2_tagh->ProjectionY("loc_h1",tagh_counter,tagh_counter);
		if(loc_h1->Integral() < 1.e3) continue;
		
		if(SUBTRACT_EMPTY) {
			TH1F *loc_h1_empty = (TH1F*)h2_tagh_empty->ProjectionY("loc_h1_empty", tagh_counter, tagh_counter);
			loc_h1->Add(loc_h1_empty, -1.0*tagh_flux[tagh_counter-1]/tagh_flux_empty[tagh_counter-1]);
			loc_h1_empty->Delete();
		}
		if(loc_h1->Integral() < 1.e3) continue;
		
		double loc_eb = tagh_en[tagh_counter-1];
		
		if(loc_eb>8.8 && loc_eb<9.0) continue;
		
		loc_h1->SetTitle(Form("TAGH Counter %d (E_{#gamma} = %.3f GeV)", tagh_counter, loc_eb));
		loc_h1->Rebin(rebins);
		loc_h1->SetLineColor(kBlack);
		loc_h1->GetXaxis()->SetTitle("E_{1} + E_{2} - E_{#gamma} [GeV]");
		loc_h1->SetLineWidth(2);
		
		//--------------------------------------------------------------//
		// Scale the e+e- simulation to match data and subtract:
		
		// phase 1 energy bin:
		int phase1_tag_counter = 0;
		int phase1_tag_sys     = 0;
		get_phase1_energy_bin(loc_eb, phase1_tag_counter, phase1_tag_sys);
		
		int n_bins_combine = 5;
		int min_counter = phase1_tag_counter-n_bins_combine;
		int max_counter = phase1_tag_counter+n_bins_combine;
		
		TH1F *loc_h1_pair;
		double loc_pair_gen_flux = 0.;
		
		if(phase1_tag_sys==0) {
			if(min_counter<1) min_counter = 1;
			else if(min_counter>127 && min_counter<179) min_counter = 179;
			
			if(max_counter>221) max_counter = 221;
			else if(max_counter>127 && max_counter<179) max_counter = 127;
			
			loc_h1_pair = (TH1F*)h2_tagh_pair->ProjectionY("loc_h1_pair", min_counter, max_counter);
			
			for(int loc_counter = min_counter; loc_counter <= max_counter; loc_counter++) {
				double loc_counter_eb   = tagh_en_phase1[loc_counter-1];
				double loc_counter_flux = h_pair_gen_flux->GetBinContent(h_pair_gen_flux->FindBin(loc_counter_eb));
				loc_pair_gen_flux += loc_counter_flux;
			}
			
		} else {
			if(min_counter<1)   min_counter = 1;
			if(max_counter>102) max_counter = 102;
			
			loc_h1_pair = (TH1F*)h2_tagm_pair->ProjectionY("loc_h1_pair", min_counter, max_counter);
			
			for(int loc_counter = min_counter; loc_counter <= max_counter; loc_counter++) {
				double loc_counter_eb   = tagm_en_phase1[loc_counter-1];
				double loc_counter_flux = h_pair_gen_flux->GetBinContent(h_pair_gen_flux->FindBin(loc_counter_eb));
				loc_pair_gen_flux += loc_counter_flux;
			}
		}
		
		double loc_pair_cs   = f_pair_cs->Eval(loc_eb) + f_triplet_cs->Eval(loc_eb);
		double loc_pair_flux = loc_pair_gen_flux / ((n_e/n_Z) * mb * loc_pair_cs);
		
		loc_h1_pair->Scale(tagh_flux[tagh_counter-1]/loc_pair_flux);
		loc_h1_pair->Rebin(rebins);
		loc_h1_pair->SetLineColor(kGreen);
		
		// fit e+e- to tail of deltaK:
		
		TF1 *f_bkgd_fit = new TF1("f_bkgd_fit", bkgd_fit_lineshape, -8.0, -3.0, 1);
		f_bkgd_fit->SetParameter(0, 1);
		h_bkgd_fit = (TH1F*)loc_h1_pair->Clone("loc_bkgd_fit");
		loc_h1->Fit(f_bkgd_fit, "R0QL");
		
		loc_h1_pair->Scale(f_bkgd_fit->GetParameter(0));
		//h_bkgd_fit->Delete();
		
		if(DEBUG_FITS && tagh_counter==45) {
			
			if(canvas3==NULL) {
				canvas3 = new TCanvas("canvas3", "Test", 800, 800);
				canvas3->SetTickx(); canvas3->SetTicky();
				canvas3->SetLeftMargin(0.13); canvas3->SetRightMargin(0.07);
		
			}
			
			TH1F *loc_h1_clone      = (TH1F*)loc_h1->Clone("loc_h1_clone_tagh_045");
			TH1F *loc_h1_pair_clone = (TH1F*)loc_h1_pair->Clone("loc_h1_pair_clone_tagh_045");
			loc_h1_clone->SetMinimum(0.);
			loc_h1_clone->GetXaxis()->SetRangeUser(-7.8,6.1);
			n_mev = 16*rebins;
			loc_h1_clone->GetYaxis()->SetTitle(Form("counts / %d MeV",n_mev));
			loc_h1_clone->GetYaxis()->SetTitleSize(0.05);
			loc_h1_clone->GetYaxis()->SetTitleOffset(1.0);
			loc_h1_clone->GetYaxis()->CenterTitle(false);
			loc_h1_clone->GetXaxis()->SetTitle("#DeltaK [GeV]");
			loc_h1_clone->GetXaxis()->SetTitleSize(0.05);
			loc_h1_clone->GetXaxis()->SetTitleOffset(0.9);
			loc_h1_clone->GetXaxis()->CenterTitle(true);
			
			//canvas3->SetLogy();
			loc_h1_clone->SetLineColor(kBlue);
			loc_h1_clone->SetMarkerColor(kBlue);
			loc_h1_pair_clone->SetLineColor(kGreen);
			loc_h1_pair_clone->SetFillColor(kGreen);
			loc_h1_pair_clone->SetMarkerColor(kGreen);
			loc_h1_pair_clone->SetFillStyle(3004);
			loc_h1_clone->Draw();
			loc_h1_pair_clone->Draw("same hist");
			
			TLatex lat;
			//lat.DrawLatexNDC(0.15, 0.835, Form("Scale Factor: %.3f", f_bkgd_fit->GetParameter(0)));
			
			TLegend *leg = new TLegend(0.15, 0.6, 0.4, 0.8);
			leg->AddEntry(loc_h1_clone,      "Experimental Data", "PE");
			leg->AddEntry(loc_h1_pair_clone, "Simulated Pair Bkgd", "PE");
			leg->Draw();
			
			canvas3->Update();
		}
		
		if(SUBTRACT_PAIR) loc_h1->Add(loc_h1_pair, -1.0);
		
		//--------------------------------------------------------------//
		
		double loc_mu    = 0, loc_sigma  = 0, loc_frac  = 0., loc_dmu  = 0., loc_csigma  = 0.;
		double loc_muE   = 0, loc_sigmaE = 0, loc_fracE = 0., loc_dmuE = 0., loc_csigmaE = 0.;
		double loc_chi2  = 0;
		
		//cout << "\n";
		//cout << "Fitting DeltaE distribution in TAGH counter " << tagh_counter << endl;
		
		int fit_val = deltaK_fit(0, 0, tagh_counter, loc_h1, 
			loc_mu, loc_muE, loc_sigma, loc_sigmaE, loc_frac, loc_fracE, loc_dmu, loc_dmuE, loc_csigma, loc_csigmaE, 
			loc_chi2, draw_fits);
		
		if(fit_val  <= 0 ) continue;
		if(loc_chi2 > 50.) continue;
		
		tagh_mu.push_back(loc_mu);
		tagh_muE.push_back(loc_muE);
		tagh_sigma.push_back(loc_sigma);
		tagh_sigmaE.push_back(loc_sigmaE);
		tagh_frac.push_back(loc_frac);
		tagh_fracE.push_back(loc_fracE);
		tagh_dmu.push_back(loc_dmu);
		tagh_dmuE.push_back(loc_dmuE);
		tagh_csigma.push_back(loc_csigma);
		tagh_csigmaE.push_back(loc_csigmaE);
		tagh_chi2.push_back(loc_chi2);
		
		tagh_fitVec.push_back(tagh_counter);
	}
	
	// Do the same for the TAGM counters:
	for(int tagm_counter=1; tagm_counter<=102; tagm_counter++) {
		
		TH1F *loc_h1 = (TH1F*)h2_tagm->ProjectionY("loc_h1",tagm_counter,tagm_counter);
		if(loc_h1->Integral() < 1.e3) continue;
		
		if(SUBTRACT_EMPTY) {
			TH1F *loc_h1_empty = (TH1F*)h2_tagm_empty->ProjectionY("loc_h1_empty", tagm_counter, tagm_counter);
			loc_h1->Add(loc_h1_empty, -1.0*tagm_flux[tagm_counter-1]/tagm_flux_empty[tagm_counter-1]);
			loc_h1_empty->Delete();
		}
		if(loc_h1->Integral() < 1.e3) continue;
		
		double loc_eb = tagm_en[tagm_counter-1];
		
		if(loc_eb>8.8 && loc_eb<9.0) continue;
		
		loc_h1->SetTitle(Form("TAGM Counter %d (E_{#gamma} = %.3f GeV)", tagm_counter, loc_eb));
		loc_h1->Rebin(rebins);
		loc_h1->SetLineColor(kBlack);
		loc_h1->GetXaxis()->SetTitle("E_{1} + E_{2} - E_{#gamma} [GeV]");
		loc_h1->SetLineWidth(2);
		
		//--------------------------------------------------------------//
		// Scale the e+e- simulation to match data and subtract:
		
		// phase 1 energy bin:
		int phase1_tag_counter = 0;
		int phase1_tag_sys     = 0;
		get_phase1_energy_bin(loc_eb, phase1_tag_counter, phase1_tag_sys);
		
		int n_bins_combine = 5;
		int min_counter = phase1_tag_counter-n_bins_combine;
		int max_counter = phase1_tag_counter+n_bins_combine;
		
		TH1F *loc_h1_pair;
		double loc_pair_gen_flux = 0.;
		
		if(phase1_tag_sys==0) {
			if(min_counter<1) min_counter = 1;
			else if(min_counter>127 && min_counter<179) min_counter = 179;
			
			if(max_counter>221) max_counter = 221;
			else if(max_counter>127 && max_counter<179) max_counter = 127;
			
			loc_h1_pair = (TH1F*)h2_tagh_pair->ProjectionY("loc_h1_pair", min_counter, max_counter);
			
			for(int loc_counter = min_counter; loc_counter <= max_counter; loc_counter++) {
				double loc_counter_eb   = tagh_en_phase1[loc_counter-1];
				double loc_counter_flux = h_pair_gen_flux->GetBinContent(h_pair_gen_flux->FindBin(loc_counter_eb));
				loc_pair_gen_flux += loc_counter_flux;
			}
			
		} else {
			if(min_counter<1)   min_counter = 1;
			if(max_counter>102) max_counter = 102;
			
			loc_h1_pair = (TH1F*)h2_tagm_pair->ProjectionY("loc_h1_pair", min_counter, max_counter);
			
			for(int loc_counter = min_counter; loc_counter <= max_counter; loc_counter++) {
				double loc_counter_eb   = tagm_en_phase1[loc_counter-1];
				double loc_counter_flux = h_pair_gen_flux->GetBinContent(h_pair_gen_flux->FindBin(loc_counter_eb));
				loc_pair_gen_flux += loc_counter_flux;
			}
		}
		
		double loc_pair_cs   = f_pair_cs->Eval(loc_eb) + f_triplet_cs->Eval(loc_eb);
		double loc_pair_flux = loc_pair_gen_flux / ((n_e/n_Z) * mb * loc_pair_cs);
		
		loc_h1_pair->Scale(tagm_flux[tagm_counter-1]/loc_pair_flux);
		loc_h1_pair->Rebin(rebins);
		loc_h1_pair->SetLineColor(kGreen);
		
		// fit e+e- to tail of deltaK:
		
		TF1 *f_bkgd_fit = new TF1("f_bkgd_fit", bkgd_fit_lineshape, -8.0, -3.0, 1);
		f_bkgd_fit->SetParameter(0, 1);
		h_bkgd_fit = (TH1F*)loc_h1_pair->Clone("loc_bkgd_fit");
		loc_h1->Fit(f_bkgd_fit, "R0QL");
		
		loc_h1_pair->Scale(f_bkgd_fit->GetParameter(0));
		//h_bkgd_fit->Delete();
		
		if(0) {
			
			if(canvas3==NULL) {
				canvas3 = new TCanvas("canvas3", "Test", 1000, 800);
			}
			loc_h1->Draw();
			loc_h1_pair->Draw("same hist");
			
			TLatex lat;
			lat.DrawLatexNDC(0.15, 0.835, Form("Scale Factor: %.3f", f_bkgd_fit->GetParameter(0)));
			
			canvas3->Update();
		}
		
		if(SUBTRACT_PAIR) loc_h1->Add(loc_h1_pair, -1.0);
		
		//--------------------------------------------------------------//
		
		double loc_mu    = 0, loc_sigma  = 0, loc_frac  = 0., loc_dmu  = 0., loc_csigma  = 0.;
		double loc_muE   = 0, loc_sigmaE = 0, loc_fracE = 0., loc_dmuE = 0., loc_csigmaE = 0.;
		double loc_chi2  = 0;
		
		//cout << "\n";
		//cout << "Fitting DeltaE distribution in TAGM counter " << tagm_counter << endl;
		
		int fit_val = deltaK_fit(0, 1, tagm_counter, loc_h1, 
			loc_mu, loc_muE, loc_sigma, loc_sigmaE, loc_frac, loc_fracE, loc_dmu, loc_dmuE, loc_csigma, loc_csigmaE, 
			loc_chi2, draw_fits);
		
		if(fit_val  <= 0 ) continue;
		if(loc_chi2 > 50.) continue;
		
		tagm_mu.push_back(loc_mu);
		tagm_muE.push_back(loc_muE);
		tagm_sigma.push_back(loc_sigma);
		tagm_sigmaE.push_back(loc_sigmaE);
		tagm_frac.push_back(loc_frac);
		tagm_fracE.push_back(loc_fracE);
		tagm_dmu.push_back(loc_dmu);
		tagm_dmuE.push_back(loc_dmuE);
		tagm_csigma.push_back(loc_csigma);
		tagm_csigmaE.push_back(loc_csigmaE);
		tagm_chi2.push_back(loc_chi2);
		
		tagm_fitVec.push_back(tagm_counter);
	}
	
	return;
}

int deltaK_fit(int is_mc, int tag_sys, int tag_counter, TH1F *h1, 
	double &mu,     double &muE, 
	double &sigma,  double &sigmaE,
	double &frac,   double &fracE,
	double &dmu,    double &dmuE,
	double &csigma, double &csigmaE, 
	double &chi2,   bool draw_fits=false)
{
	int fit_val = 1;
	
	double min_fit_x = h1->GetBinCenter(h1->FindFirstBinAbove()+1)+1.0;
	double max_fit_x = h1->GetBinCenter(h1->FindLastBinAbove() -1);
	if(max_fit_x<2.0) max_fit_x = 2.0;
	
	double eb;
	if(tag_sys==0) eb = tagh_en[tag_counter-1];
	else           eb = tagm_en[tag_counter-1];
	
	//------------------------------------------------------------------------//
	// First, fit background with polynomial
	
	TF1 *f_bkgd = new TF1("f_bkgd", bkgd_fit, min_fit_x, max_fit_x, 4);
	f_bkgd->FixParameter(2, 0.);
	f_bkgd->FixParameter(3, 0.);
	
	h1->Fit(f_bkgd, "R0QL");
	
	//------------------------------------------------------------------------//
	// Use a Gaussian to get initial guesses for mean and width
	
	double peak_guess = h1->GetBinCenter(h1->GetMaximumBin());
	if(fabs(peak_guess) > 1.0) peak_guess = 0.;
	
	double amp_guess  = h1->GetMaximum();
	if(amp_guess < 1) amp_guess = 1.;
	
	TF1 *f_gaus = new TF1("f_gaus", "gaus", -2., 2.);
	f_gaus->SetParameters(amp_guess, peak_guess, 0.65);
	f_gaus->SetRange(f_gaus->GetParameter(1) - 0.7, f_gaus->GetParameter(1) + 0.7);
	
	f_gaus->SetParLimits(0, 1.e0, 1.e6);
	f_gaus->SetParLimits(1, -1.0,  1.0);
	f_gaus->SetParLimits(2,  0.0,  1.2);
	
	h1->Fit(f_gaus, "R0QL");
	
	f_gaus->SetRange(f_gaus->GetParameter(1) - f_gaus->GetParameter(2), 
		f_gaus->GetParameter(1) + f_gaus->GetParameter(2));
	
	f_gaus->SetParameters(amp_guess, peak_guess, 0.65);
	
	f_gaus->SetParLimits(0, 1.e0, 1.e6);
	f_gaus->SetParLimits(1, -1.0,  1.0);
	f_gaus->SetParLimits(2,  0.0,  1.2);
	
	h1->Fit(f_gaus, "R0QL");
	
	//------------------------------------------------------------------------//
	// Now fit full distribution
	/*
	double loc_max_fit = f_gaus->GetParameter(1) + f_gaus->GetParameter(2);
	TF1 *f_fit = new TF1(Form("f_fit_%d_%d_%d",tag_sys,tag_counter,is_mc), "gaus(0)+pol3(3)", min_fit_x, loc_max_fit);
	
	f_fit->SetParName(0, "Const" );
	f_fit->SetParName(1, "#mu"   );
	f_fit->SetParName(2, "#sigma");
	f_fit->SetParName(3, "p0"    );
	f_fit->SetParName(4, "p1"    );
	f_fit->SetParName(5, "p2"    );
	f_fit->SetParName(6, "p3"    );
	
	f_fit->SetParameters(f_gaus->GetParameter(0), f_gaus->GetParameter(1), f_gaus->GetParameter(2));
	
	f_fit->SetParLimits(0, 1.e0, 1.e6);
	f_fit->SetParLimits(1, -1.0,  1.0);
	f_fit->SetParLimits(2,  0.0,  1.2);
	
	if(FIX_BACKGROUND) {
		f_fit->FixParameter(3, f_bkgd->GetParameter(0));
		f_fit->FixParameter(4, f_bkgd->GetParameter(1));
		f_fit->FixParameter(5, f_bkgd->GetParameter(2));
		f_fit->FixParameter(6, f_bkgd->GetParameter(3));
	} else {
		f_fit->SetParameter(3, f_bkgd->GetParameter(0));
		f_fit->SetParameter(4, f_bkgd->GetParameter(1));
		f_fit->SetParameter(5, f_bkgd->GetParameter(2));
		f_fit->SetParameter(6, f_bkgd->GetParameter(3));
	}
	
	f_fit->FixParameter(5, 0.0);
	f_fit->FixParameter(6, 0.0);
	*/
	
	double loc_max_fit_x = max_fit_x;
	//if(cut_index%4==3 && eb>9.0) {
	//	loc_max_fit_x = f_gaus->GetParameter(1) + 1.5*f_gaus->GetParameter(2);
	//}
	
	TF1 *f_fit = new TF1(Form("f_fit_%d_%d_%d",tag_sys,tag_counter,is_mc), dbl_gaus_fit, min_fit_x, loc_max_fit_x, 10);
	
	f_fit->SetParName(0, "Const" );
	f_fit->SetParName(1, "#mu"   );
	f_fit->SetParName(2, "#sigma");
	f_fit->SetParName(3, "frac");
	f_fit->SetParName(4, "#Delta#mu");
	f_fit->SetParName(5, "c_{#sigma}");
	f_fit->SetParName(6, "p0"    );
	f_fit->SetParName(7, "p1"    );
	f_fit->SetParName(8, "p2"    );
	f_fit->SetParName(9, "p3"    );
	
	f_fit->SetParameters(f_gaus->GetParameter(0), f_gaus->GetParameter(1), f_gaus->GetParameter(2));
	f_fit->SetParameter(3, 0.0);
	f_fit->SetParameter(4, 0.0);
	f_fit->SetParameter(5, 2.0);
	
	f_fit->SetParLimits(0, 1.e0, 1.e6);
	f_fit->SetParLimits(1, -1.0,  1.0);
	f_fit->SetParLimits(2,  0.0,  1.2);
	f_fit->SetParLimits(3,  0.0,  0.4);
	f_fit->SetParLimits(4, -0.5,  1.0);
	f_fit->SetParLimits(5,  1.0,  5.0);
	
	if(FIX_BACKGROUND) {
		f_fit->FixParameter(6, f_bkgd->GetParameter(0));
		f_fit->FixParameter(7, f_bkgd->GetParameter(1));
		f_fit->FixParameter(8, f_bkgd->GetParameter(2));
		f_fit->FixParameter(9, f_bkgd->GetParameter(3));
	} else {
		f_fit->SetParameter(6, f_bkgd->GetParameter(0));
		f_fit->SetParameter(7, f_bkgd->GetParameter(1));
		f_fit->SetParameter(8, f_bkgd->GetParameter(2));
		f_fit->SetParameter(9, f_bkgd->GetParameter(3));
	}
	
	f_fit->FixParameter(8, 0.0);
	f_fit->FixParameter(9, 0.0);
	
	if(!is_mc) {
		f_fit->FixParameter(3, f_frac_mc->Eval(eb));
		f_fit->FixParameter(4, f_dmu_mc->Eval(eb));
		f_fit->FixParameter(5, f_csigma_mc->Eval(eb));
	}
	
	
	//if(cut_index%4==3) {
	if(1) {
		f_fit->FixParameter(3, 0.);
		f_fit->FixParameter(4, 0.);
		f_fit->FixParameter(5, 1.);
		for(int ipar=6; ipar<10; ipar++) f_fit->FixParameter(ipar, 0.);
		
		f_fit->SetRange(f_fit->GetParameter(1) - 1.5*f_fit->GetParameter(2), f_fit->GetParameter(1) 
			+ 1.5*f_fit->GetParameter(2));
	}
	
	
	TFitResultPtr result = h1->Fit(f_fit, "SR0QL");
	chi2                 = result->Chi2() / result->Ndf();
	
	mu      = f_fit->GetParameter(1);
	muE     = f_fit->GetParError(1);
	
	sigma   = f_fit->GetParameter(2);
	sigmaE  = f_fit->GetParError(2);
	
	frac    = f_fit->GetParameter(3);
	fracE   = f_fit->GetParError(3);
	
	dmu     = f_fit->GetParameter(4);
	dmuE    = f_fit->GetParError(4);
	
	csigma  = f_fit->GetParameter(5);
	csigmaE = f_fit->GetParError(5);
	
	//------------------------------------------------------------------------//
	// Draw:
	
	if(!is_mc && tag_sys==0 && tag_counter==45) {
		
		TCanvas *loc_c1 = new TCanvas("loc_c1", "c1", 800, 800);
		loc_c1->SetTickx(); loc_c1->SetTicky();
		loc_c1->SetLeftMargin(0.13); loc_c1->SetRightMargin(0.07);
		
		TH1F *h1_clone = (TH1F*)h1->Clone("h1_tagh_045");
		h1_clone->SetMarkerColor(kBlack);
		h1_clone->SetLineColor(kBlack);
		h1_clone->SetMinimum(0.);
		
		h1_clone->SetMinimum(0.);
		h1_clone->GetXaxis()->SetRangeUser(-7.8,6.1);
		n_mev = 16*rebins;
		h1_clone->GetYaxis()->SetTitle(Form("counts / %d MeV",n_mev));
		h1_clone->GetYaxis()->SetTitleSize(0.045);
		h1_clone->GetYaxis()->SetTitleOffset(1.25);
		h1_clone->GetYaxis()->CenterTitle(false);
		h1_clone->GetXaxis()->SetTitle("#DeltaK [GeV]");
		h1_clone->GetXaxis()->SetTitleSize(0.05);
		h1_clone->GetXaxis()->SetTitleOffset(0.9);
		h1_clone->GetXaxis()->CenterTitle(true);
		
		//TF1 *f1_clone = (TF1*)f_fit->Clone("f1_tagh_045");
		TF1 *f1_clone = new TF1("f1_clone_tagh_045", "gaus", -1., 1.);
		f1_clone->SetRange(f_fit->GetParameter(1) - 1.5*f_fit->GetParameter(2), f_fit->GetParameter(1) 
			+ 1.5*f_fit->GetParameter(2));
		f1_clone->SetParameters(f_fit->GetParameter(0), f_fit->GetParameter(1), f_fit->GetParameter(2));
		h1_clone->Fit(f1_clone, "R0QL");
		f1_clone->SetRange(-3.0,3.0);
		
		h1_clone->Draw();
		f1_clone->Draw("same");
	}
	
	
	if(draw_fits) {
		
		//max_fit_x = 2.0;
		f_fit->SetNpx(1000);
		
		TF1 *f_bkgd_draw = new TF1(Form("f_bkgd_draw_%d_%d_%d",tag_sys,tag_counter,is_mc), 
			dbl_gaus_fit, min_fit_x, max_fit_x, 10);
		f_bkgd_draw->SetParameters(f_fit->GetParameters());
		f_bkgd_draw->SetParameter(0, 0.);
		f_bkgd_draw->SetLineColor(kGreen);
		f_bkgd_draw->SetLineWidth(2);
		f_bkgd_draw->SetLineStyle(2);
		
		TF1 *f_gaus_draw = new TF1(Form("f_gaus_draw_%d_%d_%d",tag_sys,tag_counter,is_mc), "gaus", min_fit_x, max_fit_x);
		f_gaus_draw->SetParameter(0, f_fit->GetParameter(0)*f_fit->GetParameter(3));
		f_gaus_draw->SetParameter(1, f_fit->GetParameter(1)+f_fit->GetParameter(4));
		f_gaus_draw->SetParameter(2, f_fit->GetParameter(2)*f_fit->GetParameter(5));
		f_gaus_draw->SetLineColor(kCyan);
		f_gaus_draw->SetLineWidth(2);
		f_gaus_draw->SetLineStyle(2);
		
		TH1F *h1_dev = new TH1F(Form("h1_dev_%d_%d_%d",tag_sys,tag_counter,is_mc), "", h1->GetXaxis()->GetNbins(), -8.0, 8.0);
		
		for(int ib=1; ib<=h1_dev->GetXaxis()->GetNbins(); ib++) {
			double loc_x   = h1->GetBinCenter(ib);
			double loc_c   = h1->GetBinContent(ib);
			double loc_err = h1->GetBinError(ib);
			/*
			if(loc_c > 0.) loc_err = sqrt(loc_c);
			else loc_err = 1.0;
			*/
			double loc_dev = (loc_c - f_fit->Eval(loc_x))/loc_err;
			h1_dev->SetBinContent(ib,loc_dev);
		}
		
		h1->GetXaxis()->SetRangeUser(min_fit_x, max_fit_x);
		
		TH1F *h1_lin = (TH1F*)h1->Clone(Form("h1_lin_%d_%d_%d",tag_sys,tag_counter,is_mc));
		h1_lin->SetMinimum(0.);
		h1_lin->GetXaxis()->SetRangeUser(min_fit_x, max_fit_x);
		h1_lin->SetTitleFont(22);
		h1_lin->GetXaxis()->SetTitle("");
		h1_lin->GetXaxis()->SetLabelOffset(1.0);
		h1_lin->GetYaxis()->SetTitle(Form("counts / %d MeV",n_mev));
		h1_lin->GetYaxis()->SetTitleSize(0.045);
		h1_lin->GetYaxis()->SetTitleOffset(1.25);
		h1_lin->GetYaxis()->SetTitleFont(22);
		
		h1_dev->GetYaxis()->SetTitle("(data - fit)/err");
		h1_dev->GetYaxis()->SetTitleSize(0.125);
		h1_dev->GetYaxis()->SetTitleOffset(0.275);
		h1_dev->GetYaxis()->SetLabelSize(0.10);
		h1_dev->GetYaxis()->SetTitleFont(22);
		
		h1_dev->GetXaxis()->SetTitle("E_{1} + E_{2} - E_{#gamma} [GeV]");
		h1_dev->GetXaxis()->SetTitleSize(0.175);
		h1_dev->GetXaxis()->SetTitleOffset(0.90);
		h1_dev->GetXaxis()->SetLabelSize(0.13);
		h1_dev->GetXaxis()->SetTitleFont(22);
		
		h1_dev->SetLineWidth(2);
		h1_dev->SetLineColor(kBlack);
		h1_dev->SetMarkerColor(kBlack);
		h1_dev->GetXaxis()->SetRangeUser(min_fit_x, max_fit_x);
		h1_dev->SetTitle(" ");
		h1_dev->GetYaxis()->SetRangeUser(-12., 12.);
		
		// Make sure our canvases have been initialized:
		
		if(canvas1==NULL) {
			canvas1 = new TCanvas("canvas1", "canvas1", 1200, 600);
			canvas1->Divide(2,1);
			
			pad_lin = (TPad*)canvas1->cd(1);
			pad_lin->SetTickx(); pad_lin->SetTicky();
			
			pad_log = (TPad*)canvas1->cd(2);
			pad_log->SetTickx(); pad_log->SetTicky();
			pad_log->SetGrid();  pad_log->SetLogy();
			
			// Canvas to show deviation between fit and data:
			
			canvas2 = new TCanvas("canvas2", "canvas2", 800, 800);
			top_pad = new TPad("top_pad", "top_pad", 0.005, 0.2025, 0.995, 0.995);
			bot_pad = new TPad("bot_pad", "bot_pad", 0.005, 0.005,  0.995, 0.1975);
			
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
		}
		
		pad_lin->cd();
		h1_lin->Draw();
		f_fit->Draw("same");
		f_bkgd_draw->Draw("same");
		f_gaus_draw->Draw("same");
		
		pad_log->cd();
		h1->Draw();
		f_fit->Draw("same");
		f_bkgd_draw->Draw("same");
		f_gaus_draw->Draw("same");
		
		canvas1->Update();
		/*
		top_pad->cd();
		h1_lin->Draw();
		f_fit->Draw("same");
		f_bkgd_draw->Draw("same");
		f_gaus_draw->Draw("same");
		
		bot_pad->cd();
		h1_dev->Draw("PE");
		
		canvas2->Update();
		*/
	}
	
	f_bkgd->Delete();
	f_gaus->Delete();
	
	return fit_val;
}



Double_t bkgd_fit(Double_t *x, Double_t *par) {
	
	Double_t xx = x[0];
	
	if(xx > -4.0 && xx < 4.0) {
		TF1::RejectPoint();
		return 0.;
	}
	
	Double_t f = par[0] + 
		     par[1] * x[0] +
		     par[2] * x[0] * x[0] + 
		     par[3] * x[0] * x[0] * x[0];
	
	//Double_t f = par[0] * exp(par[1]*x[0] + par[2]);
	
	return f;
}

Double_t dbl_gaus_fit(Double_t *x, Double_t *par) {
	
	Double_t xx = x[0];
	
	Double_t A1   = par[0];
	Double_t mu1  = par[1];
	Double_t sig1 = par[2];
	Double_t A2   = par[3] * A1;
	Double_t mu2  = par[4] + mu1;
	Double_t sig2 = par[5] * sig1;
	
	Double_t loc_x1 = (xx-mu1)/sig1;
	Double_t f1 = A1 * exp(-0.5*pow(loc_x1,2.0));
	
	Double_t loc_x2 = (xx-mu2)/sig2;
	Double_t f2 = A2 * exp(-0.5*pow(loc_x2,2.0));
	
	Double_t p0  = par[6];
	Double_t p1  = par[7];
	Double_t p2  = par[8];
	Double_t p3  = par[9];
	
	Double_t f = f1 + f2 + p0 + p1*xx + p2*xx*xx + p3*xx*xx*xx;
	
	return f;
}

Double_t bkgd_fit_lineshape(Double_t *x, Double_t *par) {
	
	Double_t xx = x[0];
	
	return (par[0]*h_bkgd_fit->GetBinContent(h_bkgd_fit->FindBin(xx)));
}
