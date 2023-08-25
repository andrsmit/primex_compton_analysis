
/*************************************************************************/
// Get energies for each tagger counter:

char tagh_xscale_fname[256], tagm_xscale_fname[256];

void get_counter_energies() {
	
	int a; double b, c;
	
	ifstream inf1(tagh_xscale_fname);
	for(int i=0; i<274; i++) {
		inf1 >> a >> b >> c;
		double deltaE  = endpoint_energy - endpoint_energy_calib;
		double emin    = b * endpoint_energy_calib  +  deltaE;
		double emax    = c * endpoint_energy_calib  +  deltaE;
		tagh_en[i] = 0.5 * (emin + emax);
	}
	inf1.close();
	
	ifstream inf2(tagm_xscale_fname);
	for(int i=0; i<102; i++) {
		inf2 >> a >> b >> c;
		double deltaE  = endpoint_energy - endpoint_energy_calib;
		double emin    = b * endpoint_energy_calib  +  deltaE;
		double emax    = c * endpoint_energy_calib  +  deltaE;
		tagm_en[i] = 0.5 * (emin + emax);
	}
	inf2.close();
	
	return;
}

/*************************************************************************/
// Read in photon flux:

char              tagh_flux_fname[256],              tagm_flux_fname[256];
char empty_target_tagh_flux_fname[256], empty_target_tagm_flux_fname[256];

void get_flux() 
{
	int a; double b, c;
	
	ifstream inf1(tagh_flux_fname);
	for(int i=0; i<274; i++) {
		inf1 >> a >> b >> c;
		tagh_flux[i]  = b;
		tagh_fluxE[i] = c;
	}
	inf1.close();
	
	ifstream inf2(tagm_flux_fname);
	for(int i=0; i<102; i++) {
		inf2 >> a >> b >> c;
		tagm_flux[i]  = b;
		tagm_fluxE[i] = c;
	}
	inf2.close();
	
	ifstream inf3(empty_target_tagh_flux_fname);
	for(int i=0; i<274; i++) {
		inf3 >> a >> b >> c;
		tagh_flux_empty[i]  = b;
		tagh_fluxE_empty[i] = c;
	}
	inf3.close();
	
	ifstream inf4(empty_target_tagm_flux_fname);
	for(int i=0; i<102; i++) {
		inf4 >> a >> b >> c;
		tagm_flux_empty[i]  = b;
		tagm_fluxE_empty[i] = c;
	}
	inf4.close();
	
	return;
}

/*************************************************************************/
// Get NLO Theory Calculation from NIST and from primex_compton generator:

char nist_theory_fname[256] = "/home/andrsmit/root_macros/compton/nist_cs.dat";

char theory_cs_pathName[256];

TGraph  *g_theory;
TGraph  *g_theory_min, *g_theory_max, *g_theory_band;
TGraph  *g_dev_min,    *g_dev_max,    *g_dev_band;
TGraph  *g_nist;

TCanvas *c_theory;

void get_theory_calc() {
	
	// read in cross section from NIST:
	
	double nist_en[22], nist_cs[22];
	ifstream inf_nist(nist_theory_fname);
	for(int i=0; i<22; i++) {
		double aa, bb;
		inf_nist >> aa >> bb;
		nist_en[i] = 1.e-3*aa;
		nist_cs[i] = 1.e3*bb/4.;
	}
	g_nist = new TGraph(22, nist_en, nist_cs);
	g_nist->SetLineColor(kRed);
	g_nist->SetMarkerColor(kRed);
	g_nist->SetMarkerStyle(23);
	g_nist->SetMarkerSize(0.7);
	
	int n_cs_files = 0;
	vector<double> gen_enVec, gen_csVec;
	
	char fname[256];
	for(int tagh_counter=1; tagh_counter<=274; tagh_counter++) {
		
		sprintf(fname, "%s/tagh_%03d.dat", theory_cs_pathName, tagh_counter);
		if(gSystem->AccessPathName(fname)) {
			continue;
		} else {
			
			double aa, bb, cc, dd, ee;
			double loc_cs = 0.;
			
			ifstream locInf(fname);
			locInf >> aa >> bb >> cc >> dd >> ee;
			locInf >> aa >> bb >> cc >> dd >> ee;
			loc_cs += cc;
			locInf >> aa >> bb >> cc >> dd >> ee;
			loc_cs += cc;
			locInf.close();
			
			gen_enVec.push_back(0.5*(aa+bb));
			gen_csVec.push_back(loc_cs);
			n_cs_files++;
		}
	}
	for(int tagm_counter=1; tagm_counter<=102; tagm_counter++) {
		
		sprintf(fname, "%s/tagm_%03d.dat", theory_cs_pathName, tagm_counter);
		if(gSystem->AccessPathName(fname)) {
			continue;
		} else {
			
			double aa, bb, cc, dd, ee;
			double loc_cs = 0.;
			
			ifstream locInf(fname);
			locInf >> aa >> bb >> cc >> dd >> ee;
			locInf >> aa >> bb >> cc >> dd >> ee;
			loc_cs += cc;
			locInf >> aa >> bb >> cc >> dd >> ee;
			loc_cs += cc;
			locInf.close();
			
			gen_enVec.push_back(0.5*(aa+bb));
			gen_csVec.push_back(loc_cs);
			n_cs_files++;
		}
	}
	
	double *gen_energy = new double[n_cs_files];
	double *gen_cs     = new double[n_cs_files];
	
	for(int i=0; i < n_cs_files; i++) {
		gen_energy[i] = gen_enVec[i];
		gen_cs[i]     = gen_csVec[i];
	}
	
	TGraph *g_theory = new TGraph(n_cs_files, gen_energy, gen_cs);
	g_theory->SetLineColor(kRed); g_theory->SetLineWidth(2);
	g_theory->SetTitle("Generated Compton Cross Section (Born+SV+DH)");
	g_theory->GetXaxis()->SetTitle("Photon Beam Energy [GeV]");
	g_theory->GetXaxis()->SetTitleSize(0.045);
	g_theory->GetXaxis()->SetTitleOffset(1.0);
	g_theory->GetXaxis()->SetLabelSize(0.04);
	g_theory->GetYaxis()->SetTitle("#sigma [mb/electron]");
	g_theory->GetYaxis()->SetTitleSize(0.055);
	g_theory->GetYaxis()->SetTitleOffset(0.8);
	g_theory->GetYaxis()->SetLabelSize(0.05);
	g_theory->SetMarkerStyle(8);
	g_theory->SetMarkerSize(0.7);
	g_theory->SetMarkerColor(kBlack);
	
	f_theory  = new TF1("f_theory", "pol5", 3., 12.);
	g_theory->Fit("f_theory", "R0Q");
	f_theory->SetLineColor(kRed);
	f_theory->SetLineStyle(2);
	
	g_nist->Fit("f_theory", "R0Q");
	
	double *gen_cs_min = new double[n_cs_files];
	double *gen_cs_max = new double[n_cs_files];
	
	for(int i=0; i < n_cs_files; i++) {
		gen_cs_min[i] = 0.95*f_theory->Eval(gen_energy[i]);
		gen_cs_max[i] = 1.05*f_theory->Eval(gen_energy[i]);
	}
	
	g_theory_min  = new TGraph(n_cs_files, gen_energy, gen_cs_min);
	g_theory_max  = new TGraph(n_cs_files, gen_energy, gen_cs_max);
	g_theory_band = new TGraph(2*n_cs_files);
	for(int i=0; i<n_cs_files; i++) {
		g_theory_band->SetPoint(i,gen_energy[i],gen_cs_max[i]);
		g_theory_band->SetPoint(n_cs_files+i,
			gen_energy[n_cs_files-i-1], gen_cs_min[n_cs_files-i-1]);
	}
	g_theory_band->SetFillStyle(3013);
	g_theory_band->SetFillColor(kCyan);
	g_theory_band->SetLineWidth(2);
	g_theory_band->SetLineColor(kCyan);
	g_theory_min->SetLineWidth(2);
	g_theory_min->SetLineColor(kCyan);
	g_theory_max->SetLineWidth(2);
	g_theory_max->SetLineColor(kCyan);
	
	double *gen_dev_min = new double[n_cs_files];
	double *gen_dev_max = new double[n_cs_files];
	
	for(int i=0; i < n_cs_files; i++) {
		gen_dev_min[i] = -5.0;
		gen_dev_max[i] =  5.0;
	}
	
	g_dev_min  = new TGraph(n_cs_files, gen_energy, gen_dev_min);
	g_dev_max  = new TGraph(n_cs_files, gen_energy, gen_dev_max);
	g_dev_band = new TGraph(2*n_cs_files);
	for(int i=0; i<n_cs_files; i++) {
		g_dev_band->SetPoint(i,gen_energy[i],gen_dev_max[i]);
		g_dev_band->SetPoint(n_cs_files+i,
			gen_energy[n_cs_files-i-1], gen_dev_min[n_cs_files-i-1]);
	}
	g_dev_band->SetFillStyle(3013);
	g_dev_band->SetFillColor(kCyan);
	g_dev_band->SetLineColor(kCyan);
	g_dev_band->SetLineWidth(2);
	g_dev_min->SetLineColor(kCyan);
	g_dev_min->SetLineWidth(2);
	g_dev_max->SetLineColor(kCyan);
	g_dev_max->SetLineWidth(2);
	
	if(DRAW_THEORY) {
		
		c_theory = new TCanvas("c_theory", "c_theory", 1000., 500.);
		c_theory->SetTickx(); c_theory->SetTicky();
		g_theory->Draw("AP");
		f_theory->Draw("same");
		g_nist->Draw("PC same");
		c_theory->Update();
	}
	
	return;
}

/*************************************************************************/
// Get e+e- pair production cross section (from NIST):

char pair_cs_fname[256];

TCanvas *c_pair_cs;

void get_pair_cs()
{
	// read in cross section from NIST:
	
	vector<double> photon_energy_vec, pair_cs_vec, trip_cs_vec;
	ifstream inf(pair_cs_fname);
	double a, b, c;
	while(inf >> a >> b >> c) {
		photon_energy_vec.push_back(a);
		pair_cs_vec.push_back(b);
		trip_cs_vec.push_back(c);
	}
	inf.close();
	
	int n_ens = (int)photon_energy_vec.size();
	double *photon_energy = new double[n_ens];
	double *pair_cs       = new double[n_ens];
	double *trip_cs       = new double[n_ens];
	for( int i=0; i<n_ens; i++ ) {
		photon_energy[i] = photon_energy_vec[i] / 1.e3;
		pair_cs[i]       = 1.e3 * pair_cs_vec[i];
		trip_cs[i]       = 1.e3 * trip_cs_vec[i];
	}
	
	TGraph *g_pair_cs = new TGraph(n_ens, photon_energy, pair_cs);
	g_pair_cs->SetTitle("e^{+}e^{-} Pair Production CS");
	g_pair_cs->GetXaxis()->SetTitle("E_{#gamma} [GeV]");
	g_pair_cs->GetYaxis()->SetTitle("#sigma [mb / atom]");
	g_pair_cs->SetMarkerColor(kBlue);
	g_pair_cs->SetMarkerStyle(8);
	
	f_pair_cs = new TF1("f_pair_cs", "pol4", 3.0, 12.0);
	f_pair_cs->SetLineColor(kBlue);
	f_pair_cs->SetLineStyle(2);
	g_pair_cs->Fit(f_pair_cs, "R0Q");
	
	TGraph *g_trip_cs = new TGraph(n_ens, photon_energy, trip_cs);
	g_trip_cs->SetTitle("e^{+}e^{-} Triplet Production CS");
	g_trip_cs->GetXaxis()->SetTitle("E_{#gamma} [GeV]");
	g_trip_cs->GetYaxis()->SetTitle("#sigma [mb / atom]");
	g_trip_cs->SetMarkerColor(kBlue);
	g_trip_cs->SetMarkerStyle(8);
	
	f_trip_cs = new TF1("f_trip_cs", "pol4", 3.0, 12.0);
	f_trip_cs->SetLineColor(kBlue);
	f_trip_cs->SetLineStyle(2);
	g_trip_cs->Fit(f_trip_cs, "R0Q");
	
	if(DRAW_THEORY) {
	
		c_pair_cs = new TCanvas("c_pair_cs", "c_pair_cs", 800, 800);
		c_pair_cs->Divide(1,2);
		
		TPad *p_pair_cs = (TPad*)c_pair_cs->cd(1);
		p_pair_cs->cd();
		p_pair_cs->SetTickx(); p_pair_cs->SetTicky();
		g_pair_cs->Draw("AP");
		f_pair_cs->Draw("same");
		
		TPad *p_trip_cs = (TPad*)c_pair_cs->cd(2);
		p_trip_cs->cd();
		p_trip_cs->SetTickx(); p_pair_cs->SetTicky();
		g_trip_cs->Draw("AP");
		f_trip_cs->Draw("same");
		
		c_pair_cs->Update();
	}
	
	return;
}

/*************************************************************************/
// Get target parameters:

void get_target_parameters(bool IS_BE_TARGET=true)
{
	if(IS_BE_TARGET) {
		
		// PrimEx-eta Be-9 Target:
		
		n_Z   = 4.0;
		f_abs = 0.987;
		
		double target_thickness = 1.77546;
		double avo_num          = 6.0221e+23;
		
		double frac_be2c = 1.e-2 * 0.021;
		double frac_fe   = 1.e-2 * 0.105;
		double frac_al   = 1.e-2 * 0.0068;
		double frac_si   = 1.e-2 * 0.0058;
		double frac_mg   = 1.e-2 * 0.0021;
		double frac_mn   = 1.e-2 * 0.0015;
		double frac_be   = 1.0 - (frac_be2c + frac_fe + frac_al + frac_si + frac_mg + frac_mn);
		
		double ne_be2c   = 14.;
		double ne_fe     = 26.;
		double ne_al     = 13.;
		double ne_si     = 14.;
		double ne_mg     = 12.;
		double ne_mn     = 25.;
		double ne_be     =  4.;
		
		double density_be2c = 1.90;
		double density_fe   = 7.874;
		double density_al   = 2.70;
		double density_si   = 2.329;
		double density_mg   = 1.738;
		double density_mn   = 7.26;
		double density_be   = 1.85;
		
		double weight_be2c  = 30.035;
		double weight_fe    = 55.84;
		double weight_al    = 26.981538;
		double weight_si    = 28.085;
		double weight_mg    = 24.305;
		double weight_mn    = 54.938;
		double weight_be    = 9.012183;
		
		ne_be2c *= (density_be2c * target_thickness * (1./weight_be2c) * avo_num);
		ne_fe   *= (density_fe   * target_thickness * (1./weight_fe  ) * avo_num);
		ne_al   *= (density_al   * target_thickness * (1./weight_al  ) * avo_num);
		ne_si   *= (density_si   * target_thickness * (1./weight_si  ) * avo_num);
		ne_mg   *= (density_mg   * target_thickness * (1./weight_mg  ) * avo_num);
		ne_mn   *= (density_mn   * target_thickness * (1./weight_mn  ) * avo_num);
		ne_be   *= (density_be   * target_thickness * (1./weight_be  ) * avo_num);
		
		n_e   = frac_be2c * ne_be2c
			+ frac_fe   * ne_fe
			+ frac_al   * ne_al
			+ frac_si   * ne_si
			+ frac_mg   * ne_mg
			+ frac_mn   * ne_mn
			+ frac_be   * ne_be;
		
	} else {
		
		// PrimEx-eta He-4 Target:
		
		n_Z   = 2.0;
		f_abs = 0.989;
		n_e   = 1.0816e+24; // p = 0.1217 g/cm3, t = 29.535cm
		
	}
	
	return;
}
