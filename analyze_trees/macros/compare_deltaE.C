#include "/work/halld/home/andrsmit/primex_compton_analysis/include/compton_inputs.cc"

//----------   Function Declarations   ----------//

TString hist_name;
int cut_index = 0;

int get_data_hists();
int get_pair_mc_hists(TString hists_fname, TString flux_fname);
int get_compton_mc_hists();

//-----------------------------------------------//

TH2F *h2_tagh,         *h2_tagm;
TH2F *h2_tagh_empty,   *h2_tagm_empty;
TH2F *h2_tagh_pair,    *h2_tagm_pair;
TH2F *h2_tagh_compton, *h2_tagm_compton;

TH1F *h_pair_gen_flux;

double m_deltaE_mu_pars[4],     m_deltaE_sigma_pars[4];
double m_deltaE_mu_pars_two[4], m_deltaE_sigma_pars_two[4];

int loadCutParameters();

void compare_deltaE(int loc_cut_index=0)
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
		root_fname = Form("%s/analyze_trees/analyze_data/rootFiles/phase1/Be_%03dnA_FIELDOFF.root", 
			loc_pathname, BEAM_CURRENT);
		empty_target_root_fname = Form("%s/analyze_trees/analyze_data/rootFiles/phase1/Be_empty_FIELDOFF.root", 
			loc_pathname);
		
		// Filled target flux filenames:
		tagh_flux_fname = Form("%s/photon_flux/phase1/Be_%03dnA_FIELDOFF_flux_tagh.txt", loc_pathname, BEAM_CURRENT);
		tagm_flux_fname = Form("%s/photon_flux/phase1/Be_%03dnA_FIELDOFF_flux_tagm.txt", loc_pathname, BEAM_CURRENT);
		
		// Empty target flux filenames:
		empty_target_tagh_flux_fname = Form("%s/photon_flux/phase1/Be_empty_FIELDOFF_flux_tagh.txt", loc_pathname);
		empty_target_tagm_flux_fname = Form("%s/photon_flux/phase1/Be_empty_FIELDOFF_flux_tagm.txt", loc_pathname);
		
		// Directory with histograms from Compton MC:
		comp_mc_dir = Form("%s/analyze_trees/analyze_mc/rootFiles/phase1/Run061321/compton/200nA", loc_pathname);
		
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
	hist_name = "deltaE/deltaE";
	
	//------------------------------------------------//
	// Get histograms from data:
	
	if(get_data_hists()) return;
	
	//------------------------------------------------//
	// Get e+e- distributions (scaled to realistic flux):
	
	TString pair_mc_hist_fname = Form("%s/analyze_trees/analyze_mc/rootFiles/phase1/Run061321/pair/pair_rec_combined.root", 
		loc_pathname);
	TString pair_mc_flux_fname = Form("%s/bhgen_test/recRootTrees/Run061321/sum.root", loc_pathname);
	if(get_pair_mc_hists(pair_mc_hist_fname, pair_mc_flux_fname)) return;
	
	//------------------------------------------------//
	// Get Compton distributions (scaled to realistic flux):
	
	if(get_compton_mc_hists()) return;
	
	//------------------------------------------------//
	
	rebins = 2;
	n_mev  = 8*rebins;
	
	int min_tagh_counter =  54;
	int max_tagh_counter =  67;
	int min_tagm_counter = 103;
	int max_tagm_counter = 102;
	
	double min_eb   = 12.0;
	double max_eb   =  6.0;
	double avg_eb   =  0.0;
	double sum_flux =  0.0;
	
	TH1F *h_data    = new TH1F("h_data",    ";#DeltaE = #left(E_{1}+E_{2}#right) - #left(E_{#gamma}-m_{e}#right) [GeV]", 
		h2_tagh->GetYaxis()->GetNbins(),         -8.0, 8.0);
	TH1F *h_pair    = new TH1F("h_pair",    ";#DeltaE = #left(E_{1}+E_{2}#right) - #left(E_{#gamma}-m_{e}#right) [GeV]", 
		h2_tagh_pair->GetYaxis()->GetNbins(),    -8.0, 8.0);
	TH1F *h_compton = new TH1F("h_compton", ";#DeltaE = #left(E_{1}+E_{2}#right) - #left(E_{#gamma}-m_{e}#right) [GeV]", 
		h2_tagh_compton->GetYaxis()->GetNbins(), -8.0, 8.0);
	
	for(int tagh_counter=min_tagh_counter; tagh_counter<=max_tagh_counter; tagh_counter++) {
		if(gSystem->AccessPathName(Form("%s/tagh_%03d.root",comp_mc_dir.Data(),tagh_counter))) continue;
		
		double loc_eb = tagh_en[tagh_counter-1];
		if(loc_eb<min_eb) min_eb = loc_eb;
		if(loc_eb>max_eb) max_eb = loc_eb;
		
		double loc_flux       = tagh_flux[tagh_counter-1];
		double loc_flux_empty = tagh_flux_empty[tagh_counter-1];
		
		avg_eb               += loc_eb*loc_flux;
		sum_flux             += loc_flux;
		
		TH1F *loc_h_data      = (TH1F*)h2_tagh->ProjectionY(Form("loc_h_data_tagh_%03d",tagh_counter), 
			tagh_counter, tagh_counter);
		TH1F *loc_h_empty     = (TH1F*)h2_tagh_empty->ProjectionY(Form("loc_h_empty_tagh_%03d",tagh_counter), 
			tagh_counter, tagh_counter);
		TH1F *loc_h_pair      = (TH1F*)h2_tagh_pair->ProjectionY(Form("loc_h_pair_tagh_%03d",tagh_counter), 
			tagh_counter, tagh_counter);
		TH1F *loc_h_compton   = (TH1F*)h2_tagh_compton->ProjectionY(Form("loc_h_compton_tagh_%03d",tagh_counter), 
			tagh_counter, tagh_counter);
		
		loc_h_empty->Scale(loc_flux/loc_flux_empty);
		
		h_data->Add(loc_h_data);
		h_data->Add(loc_h_empty, -1.0);
		h_pair->Add(loc_h_pair);
		h_compton->Add(loc_h_compton);
	}
	for(int tagm_counter=min_tagm_counter; tagm_counter<=max_tagm_counter; tagm_counter++) {
		if(gSystem->AccessPathName(Form("%s/tagm_%03d.root",comp_mc_dir.Data(),tagm_counter))) continue;
		
		double loc_eb = tagm_en[tagm_counter-1];
		if(loc_eb<min_eb) min_eb = loc_eb;
		if(loc_eb>max_eb) max_eb = loc_eb;
		
		double loc_flux       = tagm_flux[tagm_counter-1];
		double loc_flux_empty = tagm_flux_empty[tagm_counter-1];
		
		avg_eb               += loc_eb*loc_flux;
		sum_flux             += loc_flux;
		
		TH1F *loc_h_data      = (TH1F*)h2_tagm->ProjectionY(Form("loc_h_data_tagm_%03d",tagm_counter), 
			tagm_counter, tagm_counter);
		TH1F *loc_h_empty     = (TH1F*)h2_tagm_empty->ProjectionY(Form("loc_h_empty_tagm_%03d",tagm_counter), 
			tagm_counter, tagm_counter);
		TH1F *loc_h_pair      = (TH1F*)h2_tagm_pair->ProjectionY(Form("loc_h_pair_tagm_%03d",tagm_counter), 
			tagm_counter, tagm_counter);
		TH1F *loc_h_compton   = (TH1F*)h2_tagm_compton->ProjectionY(Form("loc_h_compton_tagm_%03d",tagm_counter), 
			tagm_counter, tagm_counter);
		
		loc_h_empty->Scale(loc_flux/loc_flux_empty);
		
		h_data->Add(loc_h_data);
		h_data->Add(loc_h_empty, -1.0);
		h_pair->Add(loc_h_pair);
		h_compton->Add(loc_h_compton);
	}
	
	avg_eb /= sum_flux;
	
	h_data->SetLineColor(kBlack);
	h_data->SetLineWidth(2);
	h_data->SetTitle("");//Form("%.2f GeV < E_{#gamma} < %.2f GeV", min_eb, max_eb));
	h_data->GetXaxis()->SetTitle("#DeltaE = #left(E_{1} + E_{2}#right) - #left(E_{#gamma} + m_{e}#right) [GeV]");
	h_data->GetXaxis()->SetTitleSize(0.05);
	h_data->GetXaxis()->CenterTitle(true);
	h_data->GetYaxis()->SetTitle(Form("counts / %d MeV", n_mev));
	h_data->GetYaxis()->SetTitleSize(0.05);
	h_data->GetYaxis()->SetTitleOffset(1.25);
	h_data->GetYaxis()->CenterTitle(false);
	
	//h_data->Add(h_empty,-1.0);
	
	h_pair->SetLineColor(kGreen);
	h_pair->SetMarkerColor(kGreen);
	h_pair->SetLineWidth(2);
	h_pair->SetLineStyle(2);
	
	//TH1F *h_compton = (TH1F*)h2_tagh_compton->ProjectionY("h_compton", min_tag_counter, max_tag_counter);
	h_compton->SetLineColor(kMagenta);
	h_compton->SetMarkerColor(kMagenta);
	h_compton->SetLineWidth(3);
	h_compton->SetLineStyle(7);
	
	h_data->Rebin(rebins);
	h_pair->Rebin(rebins);
	h_compton->Rebin(rebins);
	
	h_compton->Scale(0.98);
	
	TH1F *h_sum = (TH1F*)h_compton->Clone("h_sum");
	h_sum->Add(h_pair);
	h_sum->SetLineColor(kRed);
	h_sum->SetLineWidth(2);
	h_sum->SetLineStyle(1);
	
	TCanvas *c_compare = new TCanvas("c_compare", "c_compare", 1000, 800);
	c_compare->SetTickx(); c_compare->SetTicky();
	c_compare->SetTopMargin(0.08); c_compare->SetBottomMargin(0.12);
	c_compare->SetLeftMargin(0.13); c_compare->SetRightMargin(0.07);
	
	h_data->SetMinimum(0.);
	h_data->SetMarkerStyle(8);
	h_data->SetMarkerSize(0.7);
	//h_data->GetXaxis()->SetRangeUser(-3.0,1.25);
	
	TH1F *h_sub = (TH1F*)h_data->Clone("h_sub");
	h_sub->Add(h_pair,-1.0);
	
	if(loadCutParameters()) {
		cout << "Problem loading cut parameters." << endl;
		return;
	}
	
	double loc_mu = 0., loc_sig = 0.;
	for(int ipar=0; ipar<4; ipar++) {
		loc_mu += m_deltaE_mu_pars[ipar]*pow(avg_eb,(double)ipar);
	}
	loc_sig = sqrt(pow(m_deltaE_sigma_pars[0],2.0) + pow(m_deltaE_sigma_pars[1]/sqrt(avg_eb), 2.0) 
		+ pow(m_deltaE_sigma_pars[2]/avg_eb, 2.0));
	/*
	h_sub->Scale(1.0/h_sub->Integral(h_sub->FindBin(loc_mu - 10.*loc_sig*avg_eb), 
		h_sub->FindBin(loc_mu + 10.*loc_sig*avg_eb)));
	h_compton->Scale(1.0/h_compton->Integral(h_compton->FindBin(loc_mu - 10.*loc_sig*avg_eb), 
		h_compton->FindBin(loc_mu + 10.*loc_sig*avg_eb)));
	
	h_sub->GetXaxis()->SetRangeUser(loc_mu - 10.*loc_sig*avg_eb, loc_mu + 10.*loc_sig*avg_eb);
	h_compton->GetXaxis()->SetRangeUser(loc_mu - 10.*loc_sig*avg_eb, loc_mu + 10.*loc_sig*avg_eb);
	*/
	
	h_data->SetMinimum(0.);
	h_data->GetXaxis()->SetRangeUser(-6.0,2.0);
	h_data->Draw("hist");
	h_compton->Draw("same hist");
	h_pair->Draw("same hist");
	
	TLegend *leg = new TLegend(0.1643, 0.6218, 0.4940, 0.8015);
	leg->SetBorderSize(0);
	leg->AddEntry(h_data, "Data", "PE");
	leg->AddEntry(h_compton, "Compton MC", "l");
	leg->AddEntry(h_pair, "e^{+}e^{-} MC", "l");
	leg->Draw();
	
	TLatex lat_new;
	lat_new.SetTextColor(kRed-6);
	lat_new.SetTextFont(52);
	lat_new.DrawLatexNDC(0.1613, 0.8582, Form("#scale[0.75]{%.2f GeV < E_{#gamma} < %.2f GeV}", min_eb, max_eb));
	
	/*
	h_compton->Draw("same hist");
	h_pair->Draw("same hist");
	h_sum->Draw("same hist");
	
	TFile *fOut = new TFile("out.root", "RECREATE");
	fOut->cd();
	h_data->Write();
	h_compton->Write();
	h_pair->Write();
	fOut->Write();
	fOut->Close();
	
	TLegend *leg = new TLegend(0.135, 0.655, 0.437, 0.852);
	leg->SetBorderSize(1);
	leg->AddEntry(h_data, "Data", "PE");
	leg->AddEntry(h_pair, "Pair MC", "l");
	leg->AddEntry(h_compton, "Compton MC", "l");
	leg->AddEntry(h_sum, "Compton+Pair MC", "l");
	leg->Draw();
	*/
	c_compare->Update();
	
	//double n_compton_sim = h_compton->Integral(1, h_compton->GetXaxis()->GetNbins());
	double n_compton_data = h_sub->Integral(h_sub->FindBin(-3.5), h_sub->FindBin(3.5));
	double n_compton_sim  = h_compton->Integral(h_compton->FindBin(-3.5), h_compton->FindBin(3.5));
	
	vector<double> cut_sigs = {5.0};
	for(int i=0; i<(int)cut_sigs.size(); i++) {
		
		double loc_x1 = loc_mu + loc_sig*avg_eb*cut_sigs[i];
		double loc_x2 = loc_mu - loc_sig*avg_eb*cut_sigs[i];
		
		if(cut_sigs[i]<10.0) {
			TLine *l1 = new TLine(loc_x1, gPad->GetUymin(), loc_x1, gPad->GetUymax());
			l1->SetLineColor(kBlack);
			l1->SetLineStyle(2);
			//l1->Draw("same");
			TLatex lat1;
			lat1.SetTextAngle(90);
			//lat1.DrawLatex(loc_x1-0.025, 0.5*gPad->GetUymax(), Form("#scale[0.5]{%d#sigma}", (int)cut_sigs[i]));
		}
		
		TLine *l2 = new TLine(loc_x2, gPad->GetUymin(), loc_x2, gPad->GetUymax());
		l2->SetLineColor(kRed);
		l2->SetLineStyle(2);
		l2->Draw("same");
		TLatex lat2;
		lat2.SetTextAngle(90);
		lat2.DrawLatex(loc_x2-0.025, 0.5*gPad->GetUymax(), Form("#scale[0.5]{%d#sigma}", -1*(int)cut_sigs[i]));
		
		double loc_fraction_data    = h_sub->Integral(h_sub->FindBin(loc_x2), h_sub->FindBin(8.0)) 
			/ n_compton_data;
		double loc_fraction_compton = h_compton->Integral(h_compton->FindBin(loc_x2), h_compton->FindBin(8.0)) 
			/ n_compton_sim;
		char buf[256];
		cout << "" << endl;
		sprintf(buf, "  Fraction of Data events within +/-%dsigma: %f", (int)cut_sigs[i], loc_fraction_data);
		cout << buf << "\n";
		sprintf(buf, "  Fraction of Compton events within +/-%dsigma: %f", (int)cut_sigs[i], loc_fraction_compton);
		cout << buf << "\n";
	}
	
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
	
	h2_tagh_pair = new TH2F("h2_tagh_pair", "#DeltaE; TAGH Counter; E_{1} + E_{2} - E_{#gamma} [GeV]",
		274, 0.5, 274.5, 2000, -8.0, 8.0);
	h2_tagm_pair = new TH2F("h2_tagm_pair", "#DeltaE; TAGH Counter; E_{1} + E_{2} - E_{#gamma} [GeV]",
		102, 0.5, 102.5, 2000, -8.0, 8.0);
	
	TFile *f_pair = new TFile(hist_fname, "READ");
	TH2F *loc_h2_tagh_pair = (TH2F*)f_pair->Get(Form("%s_tagh_%d",hist_name.Data(),cut_index))->Clone("h2_tagh_pair");
	TH2F *loc_h2_tagm_pair = (TH2F*)f_pair->Get(Form("%s_tagm_%d",hist_name.Data(),cut_index))->Clone("h2_tagm_pair");
	loc_h2_tagh_pair->SetDirectory(0);
	loc_h2_tagm_pair->SetDirectory(0);
	f_pair->Close();
	
	TFile *f_pair_gen = new TFile(flux_fname, "READ");
	h_pair_gen_flux = (TH1F*)f_pair_gen->Get("gen_flux");
	h_pair_gen_flux->SetDirectory(0);
	f_pair_gen->Close();
	
	// scale to match photon flux:
	
	for(int tagh_counter=1; tagh_counter<=274; tagh_counter++) {
		
		double loc_eb   = tagh_en[tagh_counter-1];
		double loc_flux = tagh_flux[tagh_counter-1];
		if(loc_flux <= 0.) continue;
		
		TH1F *loc_h1 = (TH1F*)loc_h2_tagh_pair->ProjectionY("loc_h1", tagh_counter, tagh_counter);
		if(loc_h1->Integral()<1.e1) {
			loc_h1->Delete();
			continue;
		}
		
		double loc_pair_gen  = h_pair_gen_flux->GetBinContent(h_pair_gen_flux->FindBin(loc_eb));
		double loc_pair_cs   = f_pair_cs->Eval(loc_eb) + f_triplet_cs->Eval(loc_eb);
		double loc_pair_flux = loc_pair_gen / ((n_e/n_Z) * mb * loc_pair_cs);
		
		loc_h1->Scale(loc_flux/loc_pair_flux);
		for(int ibin=1; ibin<=loc_h1->GetXaxis()->GetNbins(); ibin++) {
			h2_tagh_pair->SetBinContent(tagh_counter, ibin, loc_h1->GetBinContent(ibin));
			h2_tagh_pair->SetBinError(tagh_counter, ibin, loc_h1->GetBinError(ibin));
		}
		
		loc_h1->Delete();
	}
	for(int tagm_counter=1; tagm_counter<=102; tagm_counter++) {
		
		double loc_eb   = tagm_en[tagm_counter-1];
		double loc_flux = tagm_flux[tagm_counter-1];
		if(loc_flux <= 0.) continue;
		
		TH1F *loc_h1 = (TH1F*)loc_h2_tagm_pair->ProjectionY("loc_h1", tagm_counter, tagm_counter);
		if(loc_h1->Integral()<1.e1) {
			loc_h1->Delete();
			continue;
		}
		
		double loc_pair_gen  = h_pair_gen_flux->GetBinContent(h_pair_gen_flux->FindBin(loc_eb));
		double loc_pair_cs   = f_pair_cs->Eval(loc_eb) + f_triplet_cs->Eval(loc_eb);
		double loc_pair_flux = loc_pair_gen / ((n_e/n_Z) * mb * loc_pair_cs);
		
		loc_h1->Scale(loc_flux/loc_pair_flux);
		for(int ibin=1; ibin<=loc_h1->GetXaxis()->GetNbins(); ibin++) {
			h2_tagm_pair->SetBinContent(tagm_counter, ibin, loc_h1->GetBinContent(ibin));
			h2_tagm_pair->SetBinError(tagm_counter, ibin, loc_h1->GetBinError(ibin));
		}
		
		loc_h1->Delete();
	}
	
	loc_h2_tagh_pair->Delete();
	loc_h2_tagm_pair->Delete();
	
	return 0;
}

int get_compton_mc_hists() {
	
	h2_tagh_compton = new TH2F("h2_tagh_compton", "#DeltaE; TAGH Counter; E_{1} + E_{2} - E_{#gamma} [GeV]",
		274, 0.5, 274.5, 2000, -8.0, 8.0);
	h2_tagm_compton = new TH2F("h2_tagm_compton", "#DeltaE; TAGH Counter; E_{1} + E_{2} - E_{#gamma} [GeV]",
		102, 0.5, 102.5, 2000, -8.0, 8.0);
	
	char loc_fname[256];
	
	// Loop over all TAGH counters:
	for(int tagh_counter=1; tagh_counter<=230; tagh_counter++) {
		
		sprintf(loc_fname, "%s/tagh_%03d.root", comp_mc_dir.Data(), tagh_counter);
		if(gSystem->AccessPathName(loc_fname)) continue;
		
		double loc_eb   = tagh_en[tagh_counter-1];
		double loc_flux = tagh_flux[tagh_counter-1];
		if(loc_flux <= 0.) continue;
		
		TFile *loc_fIn = new TFile(loc_fname, "READ");
		TH2F *loc_h2 = (TH2F*)loc_fIn->Get(Form("%s_tagh_%d",hist_name.Data(),cut_index));
		TH1F *loc_h1 = (TH1F*)loc_h2->ProjectionY("loc_h1",tagh_counter,tagh_counter);
		if(loc_h1->Integral()<1.e2) {
			loc_h2->Delete();
			loc_h1->Delete();
			loc_fIn->Close();
			continue;
		}
		TH1F *h_vertex = (TH1F*)loc_fIn->Get("vertex_accepted");
		
		double loc_compton_gen  = h_vertex->Integral();
		double loc_compton_cs   = f_theory->Eval(loc_eb);
		double loc_compton_flux = loc_compton_gen / (n_e * mb * loc_compton_cs);
		
		loc_h1->Scale(loc_flux/loc_compton_flux);
		for(int ibin=1; ibin<=loc_h1->GetXaxis()->GetNbins(); ibin++) {
			h2_tagh_compton->SetBinContent(tagh_counter, ibin, loc_h1->GetBinContent(ibin));
			h2_tagh_compton->SetBinError(tagh_counter, ibin, loc_h1->GetBinError(ibin));
		}
		
		loc_h2->Delete();
		loc_h1->Delete();
		loc_fIn->Close();
	}
	
	// Loop over all TAGM counters:
	for(int tagm_counter=1; tagm_counter<=102; tagm_counter++) {
		
		sprintf(loc_fname, "%s/tagm_%03d.root", comp_mc_dir.Data(), tagm_counter);
		if(gSystem->AccessPathName(loc_fname)) continue;
		
		double loc_eb   = tagm_en[tagm_counter-1];
		double loc_flux = tagm_flux[tagm_counter-1];
		if(loc_flux <= 0.) continue;
		
		TFile *loc_fIn = new TFile(loc_fname, "READ");
		TH2F *loc_h2 = (TH2F*)loc_fIn->Get(Form("%s_tagm_%d",hist_name.Data(),cut_index));
		TH1F *loc_h1 = (TH1F*)loc_h2->ProjectionY("loc_h1",tagm_counter,tagm_counter);
		if(loc_h1->Integral()<1.e2) {
			loc_h2->Delete();
			loc_h1->Delete();
			loc_fIn->Close();
			continue;
		}
		TH1F *h_vertex = (TH1F*)loc_fIn->Get("vertex_accepted");
		
		double loc_compton_gen  = h_vertex->Integral();
		double loc_compton_cs   = f_theory->Eval(loc_eb);
		double loc_compton_flux = loc_compton_gen / (n_e * mb * loc_compton_cs);
		
		loc_h1->Scale(loc_flux/loc_compton_flux);
		for(int ibin=1; ibin<=loc_h1->GetXaxis()->GetNbins(); ibin++) {
			h2_tagm_compton->SetBinContent(tagm_counter, ibin, loc_h1->GetBinContent(ibin));
			h2_tagm_compton->SetBinError(tagm_counter, ibin, loc_h1->GetBinError(ibin));
		}
		
		loc_h2->Delete();
		loc_h1->Delete();
		loc_fIn->Close();
	}
	
	return 0;
}

int loadCutParameters() {
	
	char cut_dir[256] = "/work/halld/home/andrsmit/primex_compton_analysis/analyze_trees/cuts/phase1/Run061321/200nA";
	
	char buf[256];
	ifstream loc_inf;
	double loc_pars[4];
	
	//---------------------------------------//
	// DeltaE:
	
	// Mu:
	sprintf(buf, "%s/data/deltaE_mu.dat", cut_dir);
	if(gSystem->AccessPathName(buf)) return 1;
	
	loc_inf.open(buf);
	for(int ipar=0; ipar<4; ipar++) {
		loc_inf >> loc_pars[ipar];
		m_deltaE_mu_pars[ipar] = loc_pars[ipar];
	}
	loc_inf.close();
	
	// Sigma:
	sprintf(buf, "%s/data/deltaE_sigma.dat", cut_dir);
	if(gSystem->AccessPathName(buf)) return 2;
	
	loc_inf.open(buf);
	for(int ipar=0; ipar<3; ipar++) {
		loc_inf >> loc_pars[ipar];
		m_deltaE_sigma_pars[ipar] = loc_pars[ipar];
	}
	loc_inf.close();
	
	// Mu (Two layer cut):
	sprintf(buf, "%s/data/deltaE_mu_two.dat", cut_dir);
	if(gSystem->AccessPathName(buf)) return 3;
	
	loc_inf.open(buf);
	for(int ipar=0; ipar<4; ipar++) {
		loc_inf >> loc_pars[ipar];
		m_deltaE_mu_pars_two[ipar] = loc_pars[ipar];
	}
	loc_inf.close();
	
	// Sigma:
	sprintf(buf, "%s/data/deltaE_sigma_two.dat", cut_dir);
	if(gSystem->AccessPathName(buf)) return 4;
	
	loc_inf.open(buf);
	for(int ipar=0; ipar<3; ipar++) {
		loc_inf >> loc_pars[ipar];
		m_deltaE_sigma_pars_two[ipar] = loc_pars[ipar];
	}
	loc_inf.close();
	
	return 0;
}
