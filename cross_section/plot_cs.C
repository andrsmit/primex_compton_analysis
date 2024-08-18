
void plot_cs(vector<int> tagh_counter_vec, vector<int> tagm_counter_vec) {
	
	int n_bins1 = (int)tagh_counter_vec.size();
	int n_bins2 = (int)tagm_counter_vec.size();
	int n_bins  = n_bins1+n_bins2;
	
	double *beam_energy       = new double[n_bins];
	double *zeros             = new double[n_bins];
	double *cross_section     = new double[n_bins];
	double *cross_section_err = new double[n_bins];
	double *cs_deviation      = new double[n_bins];
	double *cs_deviation_err  = new double[n_bins];
	
	for(int tagh_counter=1; tagh_counter<=274; tagh_counter++) {
		tagh_cs[tagh_counter-1]  = 0.;
		tagh_csE[tagh_counter-1] = 0.;
	}
	for(int tagm_counter=1; tagm_counter<=102; tagm_counter++) {
		tagm_cs[tagm_counter-1]  = 0.;
		tagm_csE[tagm_counter-1] = 0.;
	}
	
	for(int ib=0; ib<n_bins; ib++) {
		double loc_eb = 0.;
		double loc_cs = 0., loc_cs_err = 0.;
		if(ib<n_bins1) {
			int counter  = tagh_counter_vec[ib];
			loc_eb = tagh_en[counter-1];
			calc_cs(0, counter, loc_cs, loc_cs_err);
			tagh_cs[counter-1]  = loc_cs;
			tagh_csE[counter-1] = loc_cs_err;
		} else {
			int counter  = tagm_counter_vec[ib-n_bins1];
			loc_eb = tagm_en[counter-1];
			calc_cs(1, counter, loc_cs, loc_cs_err);
			tagm_cs[counter-1]  = loc_cs;
			tagm_csE[counter-1] = loc_cs_err;
		}
		beam_energy[ib]       = loc_eb;
		zeros[ib]             = 0.;
		cross_section[ib]     = loc_cs;
		cross_section_err[ib] = loc_cs_err;
		
		double loc_theory_cs  = f_theory->Eval(loc_eb);
		cs_deviation[ib]      = 100. * (loc_cs - loc_theory_cs) / loc_theory_cs;
		cs_deviation_err[ib]  = 100. * loc_cs_err / loc_theory_cs;
	}
	
	TGraphErrors *gCS = new TGraphErrors(n_bins, beam_energy, cross_section, zeros, 
		cross_section_err);
	gCS->SetMarkerStyle(8);
	gCS->SetMarkerSize(0.7);
	gCS->SetMarkerColor(kBlue);
	gCS->SetLineColor(kBlue);
	gCS->SetTitle("");
	
	gCS->GetXaxis()->SetTitle("Photon Beam Energy [GeV]");
	gCS->GetXaxis()->SetTitleSize(0.05);
	gCS->GetXaxis()->SetTitleOffset(0.9);
	gCS->GetXaxis()->CenterTitle(true);
	
	gCS->GetYaxis()->SetTitle("#sigma [mb / electron]");
	gCS->GetYaxis()->SetTitleSize(0.05);
	gCS->GetYaxis()->SetTitleOffset(0.9);
	gCS->GetYaxis()->CenterTitle(true);
	
	gCS->GetXaxis()->SetRangeUser(5.8, 11.6);
	gCS->GetYaxis()->SetRangeUser(0.11, 0.25);
	
	TCanvas *cCS = new TCanvas("cCS", "cCS", 1200, 500);
	cCS->SetTickx(); cCS->SetTicky();
	gCS->Draw("AP");
	f_theory->Draw("same");
	
	TGraphErrors *gDev = new TGraphErrors(n_bins, beam_energy, cs_deviation, zeros, 
		cs_deviation_err);
	gDev->SetMarkerStyle(8);
	gDev->SetMarkerSize(0.7);
	gDev->SetMarkerColor(kBlue);
	gDev->SetLineColor(kBlue);
	gDev->SetTitle("");
	
	gDev->GetXaxis()->SetTitle("Photon Beam Energy [GeV]");
	gDev->GetXaxis()->SetTitleSize(0.05);
	gDev->GetXaxis()->SetTitleOffset(0.9);
	gDev->GetXaxis()->CenterTitle(true);
	
	gDev->GetYaxis()->SetTitle("100 x (#sigma_{exp} - #sigma_{theory}) / #sigma_{theory} [%]");
	gDev->GetYaxis()->SetTitleSize(0.05);
	gDev->GetYaxis()->SetTitleOffset(0.9);
	gDev->GetYaxis()->CenterTitle(true);
	
	gDev->GetXaxis()->SetRangeUser(5.8, 11.6);
	gDev->GetYaxis()->SetRangeUser(-11., 11.);
	
	TCanvas *cDev = new TCanvas("cDev", "cDev", 1200, 500);
	cDev->SetTickx(); cDev->SetTicky();
	gDev->Draw("AP");
	cDev->Update();
	TLine *l0 = new TLine(gPad->GetUxmin(), 0., gPad->GetUxmax(), 0.);
	l0->SetLineColor(kRed);
	l0->SetLineWidth(2);
	l0->Draw("same");
	
	//---------------------------------------------------------------------------------------//
	
	TCanvas *canvas_cs = new TCanvas("canvas3", "Compton CS", 1200, 900);
	
	TPad *top_pad = new TPad("top_pad", "top_pad", 0.005, 0.3525, 0.995, 0.995);
	TPad *bot_pad = new TPad("bot_pad", "bot_pad", 0.005, 0.005,  0.995, 0.3475);
	
	top_pad->SetLeftMargin(0.10);
	top_pad->SetRightMargin(0.02);
	top_pad->SetTopMargin(0.075);
	top_pad->SetBottomMargin(0.015);
	top_pad->SetTickx(); top_pad->SetTicky();
	top_pad->SetFrameLineWidth(2);
	
	bot_pad->SetLeftMargin(0.10);
	bot_pad->SetRightMargin(0.02);
	bot_pad->SetTopMargin(0.010);
	bot_pad->SetBottomMargin(0.225);
	bot_pad->SetTickx(); bot_pad->SetTicky();
	bot_pad->SetFrameLineWidth(2);
	
	canvas_cs->cd();
	top_pad->Draw();
	bot_pad->Draw();
	
	gDev->GetXaxis()->SetTitleSize(0.1);
	gDev->GetXaxis()->SetTitleOffset(0.9);
	gDev->GetYaxis()->SetTitleSize(0.1);
	gDev->GetYaxis()->SetTitleOffset(0.25);
	gDev->GetYaxis()->SetTitle("(#sigma_{exp}-#sigma_{theory})/#sigma_{theory} [%]");
	gDev->SetTitle(" ");
	
	gCS->GetXaxis()->SetTitleSize(0.065);
	gCS->GetYaxis()->SetTitleOffset(0.5);
	
	TF1 *f_theory_clone = (TF1*)f_theory->Clone("f_theory_clone");
	f_theory_clone->SetLineWidth(4);
	
	top_pad->cd();
	gCS->Draw("AP");
	f_theory->Draw("same");
	
	TLatex lat;
	lat.SetTextColor(kBlue);
	//lat.SetTextSize(44);
	lat.SetTextFont(52);
	lat.DrawLatexNDC(0.138,0.122,"#scale[1.0]{^{9}Be Target}");
	
	TGraphErrors *gCSclone = (TGraphErrors*)gCS->Clone("gCSclone");
	gCSclone->SetMarkerSize(1.0);
	
	TLegend *leg = new TLegend(0.595, 0.645, 0.935, 0.835);
	leg->AddEntry(gCSclone, "Data" , "PE");
	leg->AddEntry(f_theory_clone, "NLO Theory" , "l");
	leg->SetBorderSize(0);
	leg->Draw();
	
	bot_pad->cd();
	gDev->Draw("AP");
	
	TF1 *f_dev = new TF1("f_dev", "pol0", 8.0, 11.6);
	gDev->Fit("f_dev", "R0");
	
	f_dev->SetRange(5.8,11.6);
	
	f_dev->SetLineColor(kGreen+1);
	f_dev->SetLineStyle(2);
	f_dev->SetLineWidth(2);
	f_dev->Draw("same");
	
	double f_dev_chi2 = f_dev->GetChisquare() / f_dev->GetNDF();
	
	TLatex lat2;
	lat2.SetTextColor(kGreen+1);
	lat2.DrawLatexNDC(0.13,0.85,Form("#scale[1.1]{Deviation: p_{0} = %.3f%% #pm %.3f%%}", 
		f_dev->GetParameter(0), f_dev->GetParError(0)));
	lat2.DrawLatexNDC(0.13,0.75,Form("#scale[1.1]{#chi^{2}/n.d.f = %.2f}", f_dev_chi2));
	
	bot_pad->Update();
	TLine *lDev = new TLine(gPad->GetUxmin(), 0.0, gPad->GetUxmax(), 0.0);
	lDev->SetLineWidth(2);
	lDev->SetLineColor(kRed);
	//lDev->Draw("same");
	
	
	return;
}
