
#include "/work/halld/home/andrsmit/primex_compton_analysis/include/compton_inputs.h"

//----------   Function Declarations   ----------//

TString hist_name;
int cut_index = 0;

int get_data_hists();
int get_pair_mc_hists(TString hists_fname, TString flux_fname);
int get_compton_mc_hists();

//-----------------------------------------------//

TH2F *h2;
TH2F *h2_empty;
TH2F *h2_pair;
TH2F *h2_compton;

TH1F *h_pair_gen_flux;

void compare_2d(int loc_cut_index=0)
{
	gStyle->SetOptStat(0); gStyle->SetOptFit(0);
	
	//----------   Set up path names   ----------//
	
	const char loc_pathname[256] = "/work/halld/home/andrsmit/primex_compton_analysis";
	
	if(true) {
		
		endpoint_energy_calib = 11.6061;
		endpoint_energy       = 11.6061;
		
		IS_BE_TARGET = true;
		IS_FIELD_OFF = true;
		BEAM_CURRENT = 200;
		
		// ROOT filenames for full and empty target data:
		root_fname = Form("%s/analyze_trees/phase1/analyze_data/rootFiles/Be_%03dnA_FIELDOFF.root", 
			loc_pathname, BEAM_CURRENT);
		empty_target_root_fname = Form("%s/analyze_trees/phase1/analyze_data/rootFiles/Be_empty_FIELDOFF.root", 
			loc_pathname);
		
		// Filled target flux filenames:
		tagh_flux_fname = Form("%s/photon_flux/phase1/Be_%03dnA_FIELDOFF_flux_tagh.txt", loc_pathname, BEAM_CURRENT);
		tagm_flux_fname = Form("%s/photon_flux/phase1/Be_%03dnA_FIELDOFF_flux_tagm.txt", loc_pathname, BEAM_CURRENT);
		
		// Empty target flux filenames:
		empty_target_tagh_flux_fname = Form("%s/photon_flux/phase1/Be_empty_FIELDOFF_flux_tagh.txt", loc_pathname);
		empty_target_tagm_flux_fname = Form("%s/photon_flux/phase1/Be_empty_FIELDOFF_flux_tagm.txt", loc_pathname);
		
		// Directory with histograms from Compton MC:
		comp_mc_dir = Form("%s/analyze_trees/phase1/analyze_mc/rootFiles/Run061321/compton", loc_pathname);
		
	} else {
		
		return;
	}
	
	// file containing the xscales for the tagh and tagm counters:
	
	tagh_xscale_fname = Form("%s/photon_flux/phase1/primex_tagh.txt", loc_pathname);
	tagm_xscale_fname = Form("%s/photon_flux/phase1/primex_tagm.txt", loc_pathname);
	
	//------------------------------------------------//
	// Get flux, tagger energies, and e+e- cross section (from NIST):
	
	get_counter_energies();
	get_flux();
	get_target_parameters();
	
	DRAW_THEORY = false;
	
	pair_cs_fname = Form("%s/photon_absorption/Be_pair_cs.dat", loc_pathname);
	get_pair_cs();
	
	theory_cs_pathName = Form("%s/compton_mc/genDir/Run061321/genCS", loc_pathname);
	get_theory_calc();
	
	cut_index = loc_cut_index;
	hist_name = "deltaK_vs_deltaE/deltaK_vs_deltaE";
	
	//------------------------------------------------//
	// Get histograms from data:
	
	if(get_data_hists()) return;
	
	//------------------------------------------------//
	// Get e+e- distributions (scaled to realistic flux):
	
	TString pair_mc_hist_fname = Form("%s/analyze_trees/phase1/analyze_mc/rootFiles/Run061321/pair/pair_rec.root", 
		loc_pathname);
	TString pair_mc_flux_fname = Form("%s/bhgen_test/recRootTrees/Run061321/sum.root", loc_pathname);
	if(get_pair_mc_hists(pair_mc_hist_fname, pair_mc_flux_fname)) return;
	
	//------------------------------------------------//
	// Get Compton distributions (scaled to realistic flux):
	
	if(get_compton_mc_hists()) return;
	
	//------------------------------------------------//
	
	double deltaE_cut_min = h2->GetYaxis()->FindBin(-2.5);
	double deltaE_cut_max = h2->GetYaxis()->FindBin( 2.5);
	
	TH1F *h_data  = (TH1F*)h2->ProjectionX("h_data", deltaE_cut_min, deltaE_cut_max);
	h_data->SetLineColor(kBlack);
	h_data->SetLineWidth(2);
	
	TH1F *h_empty = (TH1F*)h2_empty->ProjectionX("h_empty", deltaE_cut_min, deltaE_cut_max);
	h_empty->Scale(1.885);
	h_empty->SetLineColor(kBlue);
	h_empty->SetMarkerColor(kBlue);
	
	h_data->Add(h_empty,-1.0);
	
	TH1F *h_pair = (TH1F*)h2_pair->ProjectionX("h_pair", deltaE_cut_min, deltaE_cut_max);
	h_pair->SetLineColor(kGreen);
	h_pair->SetMarkerColor(kGreen);
	h_pair->SetLineWidth(2);
	h_pair->SetLineStyle(2);
	
	TH1F *h_compton = (TH1F*)h2_compton->ProjectionX("h_compton", deltaE_cut_min, deltaE_cut_max);
	h_compton->SetLineColor(kCyan);
	h_compton->SetMarkerColor(kCyan);
	h_compton->SetLineWidth(2);
	h_compton->SetLineStyle(2);
	
	rebins = 1;
	h_data->Rebin(rebins);
	h_pair->Rebin(rebins);
	h_compton->Rebin(rebins);
	
	h_pair->Scale(1.0);
	h_compton->Scale(1.0);
	
	TH1F *h_sum = (TH1F*)h_compton->Clone("h_sum");
	h_sum->Add(h_pair);
	h_sum->SetLineColor(kRed);
	h_sum->SetLineWidth(2);
	h_sum->SetLineStyle(1);
	
	TCanvas *c_compare = new TCanvas("c_compare", "c_compare", 1000, 800);
	c_compare->SetTickx(); c_compare->SetTicky();
	
	h_data->SetMinimum(0.);
	h_data->SetMarkerStyle(8);
	h_data->SetMarkerSize(0.7);
	//h_data->GetXaxis()->SetRangeUser(-3.0,1.25);
	
	h_data->Draw();
	h_pair->Draw("same hist");
	h_compton->Draw("same hist");
	h_sum->Draw("same hist");
	
	TLegend *leg = new TLegend(0.135, 0.655, 0.437, 0.852);
	leg->SetBorderSize(0);
	leg->AddEntry(h_data, "Data", "PE");
	leg->AddEntry(h_pair, "Pair MC", "l");
	leg->AddEntry(h_compton, "Compton MC", "l");
	leg->AddEntry(h_sum, "Compton+Pair MC", "l");
	leg->Draw();
	
	TH1F *h_diff = (TH1F*)h_data->Clone("h_diff");
	h_diff->Add(h_sum,-1.0);
	
	TCanvas *c_diff = new TCanvas("c_diff", "c_diff", 1000, 800);
	c_diff->SetTickx(); c_diff->SetTicky();
	h_diff->Draw("PE");
	
	
	return;
}

int get_data_hists() {
	
	if(gSystem->AccessPathName(root_fname)) {
		cout << "Specified ROOT file for filled target data does not exist." << endl;
		return 1;
	}
	else if(gSystem->AccessPathName(empty_target_root_fname)) {
		cout << "Specified ROOT file for empty target data does not exist." << endl;
		return 1;
	}
	
	TFile *fFull  = new TFile(root_fname,             "READ");
	TFile *fEmpty = new TFile(empty_target_root_fname,"READ");
	
	h2 = (TH2F*)fFull->Get(Form("%s_%d",hist_name.Data(),cut_index))->Clone("h2_data");
	h2->SetDirectory(0);
	
	h2_empty = (TH2F*)fEmpty->Get(Form("%s_%d",hist_name.Data(),cut_index))->Clone("h2_empty");
	h2_empty->SetDirectory(0);
	
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
	
	TFile *f_pair = new TFile(hist_fname, "READ");
	h2_pair = (TH2F*)f_pair->Get(Form("%s_%d",hist_name.Data(),cut_index))->Clone("h2_pair");
	h2_pair->SetDirectory(0);
	f_pair->Close();
	
	TFile *f_pair_gen = new TFile(flux_fname, "READ");
	h_pair_gen_flux = (TH1F*)f_pair_gen->Get("gen_flux");
	h_pair_gen_flux->SetDirectory(0);
	f_pair_gen->Close();
	
	double loc_pair_gen  = h_pair_gen_flux->Integral();
	double loc_pair_cs   = 175.0;
	double loc_flux = 0.;
	for(int tagh_counter=1; tagh_counter<=274; tagh_counter++) {
		loc_flux += tagh_flux[tagh_counter-1];
	}
	for(int tagm_counter=1; tagm_counter<=102; tagm_counter++) {
		loc_flux += tagm_flux[tagm_counter-1];
	}
	double loc_pair_flux = loc_pair_gen / ((n_e/n_Z) * mb * loc_pair_cs);
	h2_pair->Scale(loc_flux/loc_pair_flux);
	
	return 0;
}

int get_compton_mc_hists() {
	
	int n_bins_x      = h2->GetXaxis()->GetNbins();
	double bin_size_x = h2->GetXaxis()->GetBinCenter(2) - h2->GetXaxis()->GetBinCenter(1);
	double min_x      = h2->GetXaxis()->GetBinCenter(1) - bin_size_x/2.0;
	double max_x      = h2->GetXaxis()->GetBinCenter(n_bins_x) + bin_size_x/2.0;
	
	int n_bins_y      = h2->GetYaxis()->GetNbins();
	double bin_size_y = h2->GetYaxis()->GetBinCenter(2) - h2->GetYaxis()->GetBinCenter(1);
	double min_y      = h2->GetYaxis()->GetBinCenter(1) - bin_size_y/2.0;
	double max_y      = h2->GetYaxis()->GetBinCenter(n_bins_y) + bin_size_y/2.0;
	
	h2_compton = new TH2F("h2_compton", "", n_bins_x, min_x, max_x, n_bins_y, min_y, max_y);
	
	char loc_fname[256];
	
	double loc_compton_flux = 0.;
	double loc_total_flux = 0.;
	
	// Loop over all TAGH counters:
	for(int tagh_counter=1; tagh_counter<=230; tagh_counter++) {
		
		double loc_eb   = tagh_en[tagh_counter-1];
		double loc_flux = tagh_flux[tagh_counter-1];
		if(loc_flux <= 0.) continue;
		
		loc_total_flux += tagh_flux[tagh_counter-1];
		
		sprintf(loc_fname, "%s/tagh_%03d.root", comp_mc_dir.Data(), tagh_counter);
		if(gSystem->AccessPathName(loc_fname)) continue;
		
		loc_compton_flux += tagh_flux[tagh_counter-1];
		
		TFile *loc_fIn  = new TFile(loc_fname, "READ");
		TH2F  *loc_h2   = (TH2F*)loc_fIn->Get(Form("%s_%d",hist_name.Data(),cut_index));
		TH1F  *h_vertex = (TH1F*)loc_fIn->Get("vertex_accepted");
		
		double loc_compton_gen  = h_vertex->Integral();
		double loc_compton_cs   = f_theory->Eval(loc_eb);
		double loc_compton_flux = loc_compton_gen / (n_e * mb * loc_compton_cs);
		
		loc_h2->Scale(loc_flux/loc_compton_flux);
		
		h2_compton->Add(loc_h2);
		
		loc_h2->Delete();
		loc_fIn->Close();
	}
	
	// Loop over all TAGM counters:
	for(int tagm_counter=1; tagm_counter<=102; tagm_counter++) {
		
		double loc_eb   = tagm_en[tagm_counter-1];
		double loc_flux = tagm_flux[tagm_counter-1];
		if(loc_flux <= 0.) continue;
		
		loc_total_flux += tagm_flux[tagm_counter-1];
		
		sprintf(loc_fname, "%s/tagm_%03d.root", comp_mc_dir.Data(), tagm_counter);
		if(gSystem->AccessPathName(loc_fname)) continue;
		
		loc_compton_flux += tagm_flux[tagm_counter-1];
		
		TFile *loc_fIn  = new TFile(loc_fname, "READ");
		TH2F  *loc_h2   = (TH2F*)loc_fIn->Get(Form("%s_%d",hist_name.Data(),cut_index));
		TH1F  *h_vertex = (TH1F*)loc_fIn->Get("vertex_accepted");
		
		double loc_compton_gen  = h_vertex->Integral();
		double loc_compton_cs   = f_theory->Eval(loc_eb);
		double loc_compton_flux = loc_compton_gen / (n_e * mb * loc_compton_cs);
		
		loc_h2->Scale(loc_flux/loc_compton_flux);
		
		h2_compton->Add(loc_h2);
		
		loc_h2->Delete();
		loc_fIn->Close();
	}
	
	h2_compton->Scale(loc_total_flux/loc_compton_flux);
	
	return 0;
}
