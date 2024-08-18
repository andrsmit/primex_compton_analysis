
TF1 *f_acc;
TCanvas *canvas_acc;

TF1 *f_acc1, *f_acc2, *f_acc3, *f_acc4, *f_acc5;


bool CALC_ACC_FROM_FIT;

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
		if(7.08<loc_eb && loc_eb<7.10) continue;
		if(loc_eb>11.3) continue;
		
		TFile *fSim = new TFile(fname, "READ");
		
		TH1F *hbeam    = (TH1F*)fSim->Get(Form("%s/beam", comp_root_dir_name.Data()))->Clone(
			Form("h_beam_tagh_%d",   tagh_counter));
		TH1F *h_vertex = (TH1F*)fSim->Get(Form("%s/vertex_accepted", 
			comp_root_dir_name.Data()))->Clone(Form("h_vertex_tagh_%d", tagh_counter));
		
		double n_gen = h_vertex->Integral() * (double)hbeam->GetMean();
		//if(hbeam->GetMean() != 1.0) continue; 
		
		//-----------------------------------------------------//
		
		TH2F *h2_tagh = new TH2F(Form("h2_tagh_tagh_sim_%d",tagh_counter), "DeltaK", 
			274, 0.5, 274.5, 2000, -8.0, 8.0);
		TH2F *h2_tagm = new TH2F(Form("h2_tagh_tagm_sim_%d",tagh_counter), "DeltaK", 
			102, 0.5, 102.5, 2000, -8.0, 8.0);
		
		h2_tagh->Add((TH2F*)fSim->Get(Form("%s/%s", comp_root_dir_name.Data(), 
			hname_tagh_comp.Data())));
		h2_tagm->Add((TH2F*)fSim->Get(Form("%s/%s", comp_root_dir_name.Data(), 
			hname_tagm_comp.Data())));
		
		h2_tagh->SetDirectory(0);
		h2_tagm->SetDirectory(0);
		
		fSim->Close();
		
		TH1F *h1 = (TH1F*)h2_tagh->ProjectionY(Form("h1_sim_tagh_tagh_%d",tagh_counter));
		if(!(PHASE_VAL==1 && !IS_BE_TARGET))
			h1->Add((TH1F*)h2_tagm->ProjectionY(Form("h1_sim_tagh_tagm_%d",tagh_counter)));
		if(h1->Integral() < 1.e1) continue;
		
		h1->Rebin(rebins);
		h1->GetXaxis()->SetTitle("E_{Comp} - E_{#gamma} [GeV]");
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
	}
	
	for(int tagm_counter = 1; tagm_counter <= 102; tagm_counter++) {
		
		if(tagm_counter>=71 && tagm_counter<=73) continue;
		
		char fname[256];
		sprintf(fname, "%s/tagm_%03d.root", comp_mc_dir.Data(), tagm_counter);
		
		if(gSystem->AccessPathName(fname)) continue;
		
		//cout << "Processing TAGM Counter " << tagm_counter << endl;
		
		double loc_eb = tagm_en[tagm_counter-1];
		if(loc_eb<6.0) continue;
		
		TFile *fSim = new TFile(fname, "READ");
		
		TH1F *hbeam    = (TH1F*)fSim->Get(Form("%s/beam", comp_root_dir_name.Data()))->Clone(
			Form("h_beam_tagm_%d",   tagm_counter));
		TH1F *h_vertex = (TH1F*)fSim->Get(Form("%s/vertex_accepted", 
			comp_root_dir_name.Data()))->Clone(Form("h_vertex_tagm_%d", tagm_counter));
		
		double n_gen = h_vertex->Integral() * (double)hbeam->GetMean();
		//if(hbeam->GetMean() != 1.0) continue; 
		
		//-----------------------------------------------------//
		
		TH2F *h2_tagh = new TH2F(Form("h2_tagm_tagh_sim_%d",tagm_counter), "DeltaK", 
			274, 0.5, 274.5, 2000, -8.0, 8.0);
		TH2F *h2_tagm = new TH2F(Form("h2_tagm_tagm_sim_%d",tagm_counter), "DeltaK", 
			102, 0.5, 102.5, 2000, -8.0, 8.0);
		
		h2_tagh->Add((TH2F*)fSim->Get(Form("%s/%s", comp_root_dir_name.Data(), 
			hname_tagh_comp.Data())));
		h2_tagm->Add((TH2F*)fSim->Get(Form("%s/%s", comp_root_dir_name.Data(), 
			hname_tagm_comp.Data())));
		
		h2_tagh->SetDirectory(0);
		h2_tagm->SetDirectory(0);
		
		fSim->Close();
		
		TH1F *h1 = (TH1F*)h2_tagm->ProjectionY(Form("h1_sim_tagm_tagm_%d",tagm_counter));
		if(!(PHASE_VAL==1 && !IS_BE_TARGET))
			h1->Add((TH1F*)h2_tagh->ProjectionY(Form("h1_sim_tagm_tagh_%d",tagm_counter)));
		if(h1->Integral() < 1.e1) continue;
		
		h1->Rebin(rebins);
		h1->GetXaxis()->SetTitle("E_{Comp} - E_{#gamma} [GeV]");
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
	gAcc->GetXaxis()->SetTitle("Photon Beam Energy [GeV]");
	gAcc->SetTitle("Compton Scattering Acceptance");
	gAcc->SetMarkerStyle(8);
	gAcc->SetMarkerSize(0.5);
	gAcc->SetMarkerColor(kBlue);
	gAcc->SetLineColor(kBlue);
	gAcc->GetXaxis()->SetTitle("E_{#gamma} [GeV]");
	gAcc->GetXaxis()->SetTitleSize(0.05);
	gAcc->GetXaxis()->SetTitleOffset(0.8);
	gAcc->GetXaxis()->SetLabelSize(0.03);
	gAcc->GetYaxis()->SetTitle("N_{rec} / N_{gen}");
	gAcc->GetYaxis()->SetTitleSize(0.05);
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
	f_acc = new TF1("f_acc", "pol5", 5.0, 11.7);
	for(int ipar=0; ipar<6; ipar++) {
		f_acc->SetParName(ipar, Form("p_{%d}", ipar));
	}
	f_acc->SetLineColor(kRed);
	gAcc->Fit("f_acc", "R0ME");
	f_acc->SetRange(5.0, 12.0);
	
	if(DRAW_ACCEPTANCE) {
		
		canvas_acc = new TCanvas("canvas_acc", "canvas_acc", 800, 800);
		canvas_acc->SetTickx(); canvas_acc->SetTicky();
		canvas_acc->cd();
		
		gAcc->Draw("AP");
		gAcc1->Draw("P same");
		gAcc2->Draw("P same");
		f_acc->Draw("same");
		
		gAcc->SetTitle("");
		gAcc->GetYaxis()->SetRangeUser(0.02,0.14);
		TLatex acc_lat;
		acc_lat.DrawLatexNDC(0.12,0.15,"^{9}Be Target");
		
		canvas_acc->Update();
		
		TLatex lat_form;
		lat_form.SetTextFont(52);
		lat_form.DrawLatexNDC(0.4, 0.8, "f_{acc}#left(E_{#gamma}#right) = #sum_{i=0}^{5}#left(p_{i}E_{#gamma}^{i}#right)");
		
	}
	
	return;
}
