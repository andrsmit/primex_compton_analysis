#include "/work/halld/home/andrsmit/primex_compton_analysis/include/compton_inputs.h"

TF1 *f_acc;
TCanvas *canvas_acc;

bool DRAW_ACCEPTANCE   = false;
bool CALC_ACC_FROM_FIT = false; // use a di-Gaussian fit to DeltaK Distribution to get acceptance

double tagh_acc[274], tagh_accE[274];
double tagm_acc[102], tagm_accE[102];

void get_compton_acc() {
	
	vector<double> tagh_enVec, tagh_accVec, tagh_accEVec;
	vector<double> tagm_enVec, tagm_accVec, tagm_accEVec;
	
	for(int tagh_counter = 1; tagh_counter <= 274; tagh_counter++) {
		
		char fname[256];
		sprintf(fname, "%s/tagh_%03d.root", comp_mc_dir.Data(), tagh_counter);
		
		if(gSystem->AccessPathName(fname)) continue;
		
		//cout << "Processing TAGH Counter " << tagh_counter << endl;
		
		double loc_eb = tagh_en[tagh_counter-1];
		if(loc_eb<6.0) continue;
		
		//-----------------------------------------------------//
		
		TFile *fSim = new TFile(fname, "READ");
		
		TH1F *h_vertex = (TH1F*)fSim->Get("vertex_accepted");
		
		// Number of events generated in this energy bin:
		double n_gen = h_vertex->Integral();
		
		//-----------------------------------------------------//
		
		TH2F *h2_tagh = (TH2F*)fSim->Get(hname_tagh_comp.Data());
		TH2F *h2_tagm = (TH2F*)fSim->Get(hname_tagm_comp.Data());
		h2_tagh->SetDirectory(0);
		h2_tagm->SetDirectory(0);
		
		fSim->Close();
		
		//-----------------------------------------------------//
		
		TH1F *h1 = (TH1F*)h2_tagh->ProjectionY("h1_tagh");//, tagh_counter, tagh_counter);
		h1->Add((TH1F*)h2_tagm->ProjectionY("h1_tagm"));
		
		//if(!(PHASE_VAL==1 && !IS_BE_TARGET))
		//	h1->Add((TH1F*)h2_tagm->ProjectionY(Form("h1_sim_tagh_tagm_%d",tagh_counter)));
		
		if(h1->Integral() < 1.e2) continue;
		
		h1->Rebin(20);
		h1->GetXaxis()->SetTitle("E_{Comp}#left(#theta_{1},#theta_{2}#right) - E_{#gamma} [GeV]");
		h1->GetXaxis()->CenterTitle(true);
		h1->GetXaxis()->SetTitleOffset(1.1);
		h1->GetYaxis()->SetTitle(Form("counts / %d MeV", n_mev));
		h1->GetYaxis()->SetTitleOffset(1.3);
		h1->SetLineColor(kBlack);
		h1->SetLineWidth(2);
		
		double n_rec = h1->Integral();
		
		double loc_acc  = n_rec / n_gen;
		double loc_accE = sqrt(n_gen*loc_acc*(1.-loc_acc)) / n_gen;
		
		tagh_acc[tagh_counter-1]  = loc_acc;
		tagh_accE[tagh_counter-1] = loc_accE;
		
		tagh_enVec.push_back(loc_eb);
		tagh_accVec.push_back(loc_acc);
		tagh_accEVec.push_back(loc_accE);
		
		h1->Delete();
		h2_tagh->Delete();
		h2_tagm->Delete();
	}
	
	for(int tagm_counter = 1; tagm_counter <= 102; tagm_counter++) {
		
		char fname[256];
		sprintf(fname, "%s/tagm_%03d.root", comp_mc_dir.Data(), tagm_counter);
		
		if(gSystem->AccessPathName(fname)) continue;
		
		//cout << "Processing TAGM Counter " << tagm_counter << endl;
		
		double loc_eb = tagm_en[tagm_counter-1];
		if(loc_eb<6.0) continue;
		
		//-----------------------------------------------------//
		
		TFile *fSim = new TFile(fname, "READ");
		
		TH1F *h_vertex = (TH1F*)fSim->Get("vertex_accepted");
		
		// Number of events generated in this energy bin:
		double n_gen = h_vertex->Integral();
		
		//-----------------------------------------------------//
		
		TH2F *h2_tagh = (TH2F*)fSim->Get(hname_tagh_comp.Data());
		TH2F *h2_tagm = (TH2F*)fSim->Get(hname_tagm_comp.Data());
		h2_tagh->SetDirectory(0);
		h2_tagm->SetDirectory(0);
		
		fSim->Close();
		
		//-----------------------------------------------------//
		
		TH1F *h1 = (TH1F*)h2_tagm->ProjectionY("h1_tagm");//, tagm_counter, tagm_counter);
		h1->Add((TH1F*)h2_tagh->ProjectionY("h1_tagh"));
		
		//if(!(PHASE_VAL==1 && !IS_BE_TARGET))
		//	h1->Add((TH1F*)h2_tagm->ProjectionY(Form("h1_sim_tagm_tagm_%d",tagm_counter)));
		
		if(h1->Integral() < 1.e2) continue;
		
		h1->Rebin(20);
		h1->GetXaxis()->SetTitle("E_{Comp}#left(#theta_{1},#theta_{2}#right) - E_{#gamma} [GeV]");
		h1->GetXaxis()->CenterTitle(true);
		h1->GetXaxis()->SetTitleOffset(1.1);
		h1->GetYaxis()->SetTitle(Form("counts / %d MeV", n_mev));
		h1->GetYaxis()->SetTitleOffset(1.3);
		h1->SetLineColor(kBlack);
		h1->SetLineWidth(2);
		
		double n_rec = h1->Integral();
		
		double loc_acc  = n_rec / n_gen;
		double loc_accE = sqrt(n_gen*loc_acc*(1.-loc_acc)) / n_gen;
		
		tagm_acc[tagm_counter-1]  = loc_acc;
		tagm_accE[tagm_counter-1] = loc_accE;
		
		tagm_enVec.push_back(loc_eb);
		tagm_accVec.push_back(loc_acc);
		tagm_accEVec.push_back(loc_accE);
		
		h1->Delete();
		h2_tagh->Delete();
		h2_tagm->Delete();
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
	for(int i=0; i<n_bins; i++) {
		if(i<n_bins1) {
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
	gAcc->GetXaxis()->SetTitle("Photon Beam Energy [GeV]");
	gAcc->SetTitle("Compton Scattering Acceptance");
	gAcc->SetMarkerStyle(8);
	gAcc->SetMarkerSize(0.5);
	gAcc->SetMarkerColor(kBlue);
	gAcc->SetLineColor(kBlue);
	gAcc->GetXaxis()->SetTitle("E_{#gamma} [GeV]");
	gAcc->GetXaxis()->SetTitleSize(0.055);
	gAcc->GetXaxis()->SetTitleOffset(0.8);
	gAcc->GetXaxis()->SetLabelSize(0.03);
	gAcc->GetYaxis()->SetTitle("N_{rec} / N_{gen}");
	gAcc->GetYaxis()->SetTitleSize(0.055);
	gAcc->GetYaxis()->SetTitleOffset(0.8);
	gAcc->GetYaxis()->SetLabelSize(0.03);
	
	TGraphErrors *gAcc1 = new TGraphErrors(n_bins1, energy1, acc1, zero1, acc1E);
	gAcc1->GetXaxis()->SetTitle("Photon Beam Energy [GeV]");
	gAcc1->SetTitle("Compton Acceptance");
	gAcc1->SetMarkerStyle(8);
	gAcc1->SetMarkerSize(0.5);
	gAcc1->SetMarkerColor(kBlue);
	gAcc1->SetLineColor(kBlue);
	TGraphErrors *gAcc2 = new TGraphErrors(n_bins2, energy2, acc2, zero2, acc2E);
	gAcc2->GetXaxis()->SetTitle("Photon Beam Energy [GeV]");
	gAcc2->SetTitle("Compton Acceptance");
	gAcc2->SetMarkerStyle(8);
	gAcc2->SetMarkerSize(0.5);
	gAcc2->SetMarkerColor(kBlue);
	gAcc2->SetLineColor(kBlue);
	
	//f_acc = new TF1("f_acc", "[0] + [1]*x + [2]*x^2 + [3]*x^3 + [4]*x^4 + [5]*x^5", 5.0, 11.7);
	f_acc = new TF1("f_acc", "pol5", 5.0, 11.2);
	for(int ipar=0; ipar<6; ipar++) {
		f_acc->SetParName(ipar, Form("p%d", ipar));
	}
	f_acc->SetLineColor(kRed);
	gAcc->Fit("f_acc", "R0ME");
	f_acc->SetRange(5.0, 12.0);
	
	if(DRAW_ACCEPTANCE) {
		
		canvas_acc = new TCanvas("canvas_acc", "canvas_acc", 800, 800);
		canvas_acc->SetTickx(); canvas_acc->SetTicky();
		canvas_acc->SetRightMargin(0.05);
		canvas_acc->SetLeftMargin(0.15);
		canvas_acc->SetTopMargin(0.07);
		canvas_acc->SetBottomMargin(0.13);
		
		canvas_acc->cd();
		
		gAcc->Draw("AP");
		gAcc1->Draw("P same");
		gAcc2->Draw("P same");
		f_acc->Draw("same");
		
		gAcc->SetTitle("");
		//gAcc->GetYaxis()->SetRangeUser(0.02,0.14);
		gAcc->GetYaxis()->SetTitleOffset(1.1);
		gAcc->GetXaxis()->CenterTitle(true);
		gAcc->GetYaxis()->CenterTitle(true);
		TLatex acc_lat;
		acc_lat.SetTextFont(52);
		
		acc_lat.DrawLatexNDC(0.19, 0.20, Form("#color[4]{^{%d}%s Target}", (int)n_A, TARGET_STR.Data()));
		
		canvas_acc->Update();
		
		TLatex lat_form;
		lat_form.SetTextFont(52);
		lat_form.DrawLatexNDC(0.4, 0.8, "f_{acc}#left(E_{#gamma}#right) = #sum_{i=0}^{5}#left(p_{i}E_{#gamma}^{i}#right)");
	}
	
	return;
}
