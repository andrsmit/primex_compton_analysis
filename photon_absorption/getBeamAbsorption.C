
#include "/work/halld/home/andrsmit/primex_compton_analysis/include/compton_inputs.cc"

void getBeamAbsorption()
{
	//-----   Read in Photon Flux and Tagger Counter Energies   -----//
	
	const char pathName[256] = "/work/halld/home/andrsmit/primex_compton_analysis";
	
	// photon flux file names:
	
	tagh_flux_fname = Form("%s/photon_flux/phase1/Be_200nA_FIELDOFF_flux_tagh.txt", pathName);
	tagm_flux_fname = Form("%s/photon_flux/phase1/Be_200nA_FIELDOFF_flux_tagm.txt", pathName);
	
	// file containing the xscales for the tagh and tagm counters:
	
	tagh_xscale_fname = Form("%s/photon_flux/phase1/primex_tagh.txt", pathName);
	tagm_xscale_fname = Form("%s/photon_flux/phase1/primex_tagm.txt", pathName);
	
	endpoint_energy_calib = 11.6061;
	endpoint_energy       = 11.608;
	//endpoint_energy     = 11.6061; // Be 200 nA data
	//endpoint_energy     = 11.1671; // He  50 nA data
	//endpoint_energy     = 11.1664; // He 100 nA data
	//endpoint_energy     = 11.1666; // He 200 nA data
	
	get_flux();
	get_counter_energies();
	
	cout << tagh_en[120] << endl;
	
	double      tagh_fabs[274],      tagm_fabs[102];
	double tagh_fabs_high[274], tagm_fabs_high[102];
	double  tagh_fabs_low[274],  tagm_fabs_low[102];
	
	for(int tagh_counter=1; tagh_counter<=274; tagh_counter++) {
		tagh_fabs[tagh_counter-1]      = 0.;
		tagh_fabs_high[tagh_counter-1] = 0.;
		tagh_fabs_low[tagh_counter-1]  = 0.;
	}
	for(int tagm_counter=1; tagm_counter<=102; tagm_counter++) {
		tagm_fabs[tagm_counter-1]      = 0.;
		tagm_fabs_high[tagm_counter-1] = 0.;
		tagm_fabs_low[tagm_counter-1]  = 0.;
	}
	
	//---------------------------------------------------------------//
	
	const double b2cm2      = 1.e24;         // barn to cm^2
	const double AvogNumber = 6.0221409e+23; // Avogadro's number [atoms/mol]
	
	// Be Target:
	
	TString nist_cs_fname = Form("%s/photon_absorption/Be_pair_cs.dat", pathName);
	
	const double TargetDensity = 1.848;         // g/cm^3
	const double TargetWidth   = 1.77546;       // cm
	const double atomicMass    = 9.012182;      // atomic mass of Beryllium [g/mol]
	
	/*
	// He Target:
	
	char nist_cs_fname[256] = "He_pair_cs.dat";
	
	const double TargetDensity =  0.1217;        // g/cm^3
	const double TargetWidth   = 29.5;           // cm
	const double atomicMass    =  4.0026;        // atomic mass of Beryllium [g/mol]
	*/
	
	//---------------------------------------------------------------//
	
	// Get the pair production cs from file:
	
	ifstream inf_nist(nist_cs_fname.Data());
	if(!inf_nist.good()) {
		cout << "NIST e+e- cross section data does not exist. Check filename." << endl;
		return 1;
	}
	
	vector<double> photon_energy_vec, pair_cs_vec;
	double loc_e, loc_pair_cs, loc_triplet_cs;
	while(inf_nist >> loc_e >> loc_pair_cs >> loc_triplet_cs) {
		photon_energy_vec.push_back(loc_e);
		pair_cs_vec.push_back(loc_pair_cs+loc_triplet_cs);
		//pair_cs_vec.push_back(loc_pair_cs+loc_triplet_cs);
	}
	inf_nist.close();
	
	// Convert the e+e- cross section to an attenuation coefficient:
	
	int n_ens = (int)photon_energy_vec.size();
	double *photon_energy    = new double[n_ens];
	double *atten_coeff      = new double[n_ens];
	double *atten_coeff_high = new double[n_ens];
	double *atten_coeff_low  = new double[n_ens];
	for(int i=0; i<n_ens; i++) {
		photon_energy[i]    = photon_energy_vec[i] / 1.e3;
		double loc_atten    = pair_cs_vec[i] / (atomicMass/AvogNumber) / b2cm2;
		atten_coeff[i]      = loc_atten;
		atten_coeff_high[i] = 1.05*loc_atten;
		atten_coeff_low[i]  = 0.95*loc_atten;
	}
	
	TGraph *gAtten = new TGraph(n_ens, photon_energy, atten_coeff);
	gAtten->SetTitle("Photon Attenuation Coefficient");
	gAtten->SetMarkerColor(kBlue);
	gAtten->SetMarkerStyle(8);
	gAtten->GetYaxis()->SetRangeUser(0.94*gAtten->GetY()[0], 1.06*gAtten->GetY()[n_ens-1]);
	gAtten->GetXaxis()->SetTitle("Photon Beam Energy (GeV)");
	gAtten->GetXaxis()->SetTitleSize(0.05);
	gAtten->GetXaxis()->SetTitleOffset(1.0);
	gAtten->GetXaxis()->CenterTitle(true);
	gAtten->GetYaxis()->SetTitle("#mu [cm^{2}/g]");
	gAtten->GetYaxis()->SetTitleSize(0.05);
	gAtten->GetYaxis()->SetTitleOffset(1.0);
	gAtten->GetYaxis()->CenterTitle(true);
	
	TF1 *fAtten = new TF1("fAtten", "pol4", 3.0, 12.0);
	gAtten->Fit(fAtten, "R0Q");
	
	TGraph *gAtten_high = new TGraph(n_ens, photon_energy, atten_coeff_high);
	gAtten_high->SetTitle("Photon Attenuation Coefficient");
	gAtten_high->GetXaxis()->SetTitle("E_{#gamma} [GeV]");
	gAtten_high->GetYaxis()->SetTitle("#mu [cm^{2}/g]");
	gAtten_high->SetMarkerColor(kGreen);
	gAtten_high->SetMarkerStyle(8);
	
	TF1 *fAtten_high = new TF1("fAtten_high", "pol4", 3.0, 12.0);
	gAtten_high->Fit(fAtten_high, "R0Q");
	
	TGraph *gAtten_low = new TGraph(n_ens, photon_energy, atten_coeff_low);
	gAtten_low->SetTitle("Photon Attenuation Coefficient");
	gAtten_low->GetXaxis()->SetTitle("E_{#gamma} [GeV]");
	gAtten_low->GetYaxis()->SetTitle("#mu [cm^{2}/g]");
	gAtten_low->SetMarkerColor(kGreen);
	gAtten_low->SetMarkerStyle(8);
	
	TF1 *fAtten_low = new TF1("fAtten_low", "pol4", 3.0, 12.0);
	gAtten_low->Fit(fAtten_low, "R0Q");
	
	TCanvas *cAtten = new TCanvas("cAtten", "cAtten", 800, 400);
	cAtten->SetTickx(); cAtten->SetTicky();
	gAtten->Draw("AP"); fAtten->Draw("same");
	gAtten_high->Draw("P same"); fAtten_high->Draw("same");
	gAtten_low->Draw("P same");  fAtten_low->Draw("same");
	
	cAtten->Update();
	
	//-----   Loop over TAGH Counters   -----//
	
	TRandom3 *rand = new TRandom3(0);
	
	for(int tagh_counter = 1; tagh_counter <= 274; tagh_counter++) 
	{
		double loc_eb  = tagh_en[tagh_counter-1];
		
		double rl      = TargetDensity * fAtten->Eval(loc_eb);
		double rl_high = TargetDensity * fAtten_high->Eval(loc_eb);
		double rl_low  = TargetDensity * fAtten_low->Eval(loc_eb);
		
		double abs_counter = 0., total_counter = 0.;
		double abs_counter_high = 0.;
		double abs_counter_low  = 0.;
		
		for(int ievt = 0; ievt < 1.e6; ievt++) {
			
			// Assume Compton interaction vertex is uniformly distributed in z:
			
			double locVertex = TargetWidth * rand->Uniform();
			
			double a_E      = TMath::Exp(-locVertex * rl     );
			double a_E_high = TMath::Exp(-locVertex * rl_high);
			double a_E_low  = TMath::Exp(-locVertex * rl_low );
			
			a_E      = (1./(locVertex*rl)) * (1.0 - a_E);
			a_E_high = (1./(locVertex*rl)) * (1.0 - a_E_high);
			a_E_low  = (1./(locVertex*rl)) * (1.0 - a_E_low);
			
			// a_E is the probability of a pair conversion occuring before this point
			
			double loc_rand = rand->Uniform();
			
			if(loc_rand > a_E)      abs_counter      += 1.0;
			if(loc_rand > a_E_high) abs_counter_high += 1.0;
			if(loc_rand > a_E_low)  abs_counter_low  += 1.0;
			
			total_counter += 1.0;
		}
		
		double absorption_factor      = 1.0 - (abs_counter     /total_counter);
		double absorption_factor_high = 1.0 - (abs_counter_high/total_counter);
		double absorption_factor_low  = 1.0 - (abs_counter_low /total_counter);
		
		tagh_fabs[tagh_counter-1]      = absorption_factor;
		tagh_fabs_high[tagh_counter-1] = absorption_factor_high;
		tagh_fabs_low[tagh_counter-1]  = absorption_factor_low;
		
		cout << "TAGH Counter " << tagh_counter << "; f_abs = " << absorption_factor << endl;
	}
	
	for(int tagm_counter = 1; tagm_counter <= 102; tagm_counter++) 
	{
		double loc_eb  = tagm_en[tagm_counter-1];
		
		double rl      = TargetDensity * fAtten->Eval(loc_eb);
		double rl_high = TargetDensity * fAtten_high->Eval(loc_eb);
		double rl_low  = TargetDensity * fAtten_low->Eval(loc_eb);
		
		double abs_counter = 0., total_counter = 0.;
		double abs_counter_high = 0.;
		double abs_counter_low  = 0.;
		
		for(int ievt = 0; ievt < 1.e6; ievt++) {
			
			// Assume Compton interaction vertex is uniformly distributed in z:
			
			double locVertex = TargetWidth * rand->Uniform();
			
			double a_E      = TMath::Exp(-locVertex * rl     );
			double a_E_high = TMath::Exp(-locVertex * rl_high);
			double a_E_low  = TMath::Exp(-locVertex * rl_low );
			
			a_E      = (1./(locVertex*rl)) * (1.0 - a_E);
			a_E_high = (1./(locVertex*rl)) * (1.0 - a_E_high);
			a_E_low  = (1./(locVertex*rl)) * (1.0 - a_E_low);
			
			// a_E is the probability of a pair conversion occuring before this point
			
			double loc_rand = rand->Uniform();
			
			if(loc_rand > a_E)      abs_counter      += 1.0;
			if(loc_rand > a_E_high) abs_counter_high += 1.0;
			if(loc_rand > a_E_low)  abs_counter_low  += 1.0;
			
			total_counter += 1.0;
		}
		
		double absorption_factor      = 1.0 - (abs_counter     /total_counter);
		double absorption_factor_high = 1.0 - (abs_counter_high/total_counter);
		double absorption_factor_low  = 1.0 - (abs_counter_low /total_counter);
		
		tagm_fabs[tagm_counter-1]      = absorption_factor;
		tagm_fabs_high[tagm_counter-1] = absorption_factor_high;
		tagm_fabs_low[tagm_counter-1]  = absorption_factor_low;
		
		cout << "TAGM Counter " << tagm_counter << "; f_abs = " << absorption_factor << endl;
	}
	
	
	TGraph *gAbsTAGH = new TGraph( 274, tagh_en, tagh_fabs );
	gAbsTAGH->GetXaxis()->SetTitle( "Photon Beam Energy [GeV]" );
	gAbsTAGH->GetYaxis()->SetTitle( "f_{abs}" );
	gAbsTAGH->SetTitle( "Photon Beam Absorption Factor" );
	gAbsTAGH->SetMarkerStyle(8);
	gAbsTAGH->SetMarkerColor(kBlue);
	
	TGraph *gAbsTAGM = new TGraph( 102, tagm_en, tagm_fabs );
	gAbsTAGM->GetXaxis()->SetTitle( "Photon Beam Energy [GeV]" );
	gAbsTAGM->GetYaxis()->SetTitle( "f_{abs}" );
	gAbsTAGM->SetTitle( "Photon Beam Absorption Factor" );
	gAbsTAGM->SetMarkerStyle(8);
	gAbsTAGM->SetMarkerColor(kBlue);
	
	
	gAbsTAGH->GetYaxis()->SetRangeUser( 0.015, 0.025 );
	
	
	TCanvas *cAbs = new TCanvas( "cAbs", "cAbs", 1000, 600 );
	cAbs->SetTickx(); cAbs->SetTicky();
	gAbsTAGH->Draw("AP");
	gAbsTAGM->Draw("P same");
	
	
	TGraph *gAbsTAGH_high = new TGraph( 274, tagh_en, tagh_fabs_high );
	gAbsTAGH_high->GetXaxis()->SetTitle( "Photon Beam Energy [GeV]" );
	gAbsTAGH_high->GetYaxis()->SetTitle( "f_{abs}" );
	gAbsTAGH_high->SetTitle( "Photon Beam Absorption Factor" );
	gAbsTAGH_high->SetMarkerStyle(8);
	gAbsTAGH_high->SetMarkerColor(kRed);
	
	TGraph *gAbsTAGM_high = new TGraph( 102, tagm_en, tagm_fabs_high );
	gAbsTAGM_high->GetXaxis()->SetTitle( "Photon Beam Energy [GeV]" );
	gAbsTAGM_high->GetYaxis()->SetTitle( "f_{abs}" );
	gAbsTAGM_high->SetTitle( "Photon Beam Absorption Factor" );
	gAbsTAGM_high->SetMarkerStyle(8);
	gAbsTAGM_high->SetMarkerColor(kRed);
	
	
	TGraph *gAbsTAGH_low = new TGraph( 274, tagh_en, tagh_fabs_low );
	gAbsTAGH_low->GetXaxis()->SetTitle( "Photon Beam Energy [GeV]" );
	gAbsTAGH_low->GetYaxis()->SetTitle( "f_{abs}" );
	gAbsTAGH_low->SetTitle( "Photon Beam Absorption Factor" );
	gAbsTAGH_low->SetMarkerStyle(8);
	gAbsTAGH_low->SetMarkerColor(kGreen);
	
	TGraph *gAbsTAGM_low = new TGraph( 102, tagm_en, tagm_fabs_low );
	gAbsTAGM_low->GetXaxis()->SetTitle( "Photon Beam Energy [GeV]" );
	gAbsTAGM_low->GetYaxis()->SetTitle( "f_{abs}" );
	gAbsTAGM_low->SetTitle( "Photon Beam Absorption Factor" );
	gAbsTAGM_low->SetMarkerStyle(8);
	gAbsTAGM_low->SetMarkerColor(kGreen);
	
	cAbs->cd();
	gAbsTAGH_high->Draw("P same"); gAbsTAGH_low->Draw("P same");
	gAbsTAGM_high->Draw("P same"); gAbsTAGM_low->Draw("P same");
	
	
	TLegend *leg = new TLegend( 0.15, 0.65, 0.45, 0.85 );
	leg->AddEntry( gAbsTAGH,      "NIST Pair CS",      "P" );
	leg->AddEntry( gAbsTAGH_high, "NIST Pair CS + 5%", "P" );
	leg->AddEntry( gAbsTAGH_low,  "NIST Pair CS - 5%", "P" );
	leg->Draw();
	
	
	char buf[256];
	
	ofstream outf_tagh( "He_tagh_fabs.dat" );
	for( int i=0; i<274; i++ ) {
		sprintf( buf, "%03d   %.7f   %.7f   %.7f", i+1, tagh_fabs[i], tagh_fabs_high[i], 
			tagh_fabs_low[i] );
		outf_tagh << buf << "\n";
	}
	outf_tagh.close();
	
	
	ofstream outf_tagm( "He_tagm_fabs.dat" );
	for( int i=0; i<102; i++ ) {
		sprintf( buf, "%03d   %.7f   %.7f   %.7f", i+1, tagm_fabs[i], tagm_fabs_high[i], 
			tagm_fabs_low[i] );
		outf_tagm << buf << "\n";
	}
	outf_tagm.close();
	
	
	return;
}
