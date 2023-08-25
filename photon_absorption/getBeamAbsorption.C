
char   tagh_flux_fname[256],   tagm_flux_fname[256];
char tagh_xscale_fname[256], tagm_xscale_fname[256];

double    tagh_en[274],    tagm_en[102];
double  tagh_flux[274],  tagm_flux[102];
double tagh_fluxE[274], tagm_fluxE[102];

void get_flux();
void get_counter_energies();

double endpoint_energy, endpoint_energy_calib;


void getBeamAbsorption()
{
	
	//-----   Read in Photon Flux and Tagger Counter Energies   -----//
	
	const char pathName[256] = "/work/halld/home/andrsmit/primex_compton_analysis";
	
	// photon flux file names:
	
	sprintf( tagh_flux_fname, "%s/photon_flux/He_200_tagh_flux.txt", pathName );
	sprintf( tagm_flux_fname, "%s/photon_flux/He_200_tagm_flux.txt", pathName );
	
	// file containing the xscales for the tagh and tagm counters:
	
	sprintf( tagh_xscale_fname, "%s/photon_flux/primex_tagh.txt", pathName );
	sprintf( tagm_xscale_fname, "%s/photon_flux/primex_tagm.txt", pathName );
	
	endpoint_energy_calib = 11.6061;
	//endpoint_energy     = 11.6061; // Be 200 nA data
	//endpoint_energy     = 11.1671; // He  50 nA data
	//endpoint_energy     = 11.1664; // He 100 nA data
	endpoint_energy       = 11.1666; // He 200 nA data
	
	get_flux();
	get_counter_energies();
	
	double      tagh_fabs[274],      tagm_fabs[102];
	double tagh_fabs_high[274], tagm_fabs_high[102];
	double  tagh_fabs_low[274],  tagm_fabs_low[102];
	
	for( int i=0; i<274; i++ ) {
		tagh_fabs[i]      = 0.;
		tagh_fabs_high[i] = 0.;
		tagh_fabs_low[i]  = 0.;
	} for( int i=0; i<102; i++ ) {
		tagm_fabs[i]      = 0.;
		tagm_fabs_high[i] = 0.;
		tagm_fabs_low[i]  = 0.;
	}
	
	//---------------------------------------------------------------//
	
	// Be Target:
	/*
	const double TargetDensity = 1.85;          // g/cm^3
	const double TargetWidth   = 1.77546;       // cm
	const double atomicMass    = 9.012182;      // atomic mass of Beryllium [g/mol]
	const double AvogNumber    = 6.0221409e+23; // Avogadro's number [atoms/mol]
	*/
	
	// He Target:
	
	const double TargetDensity =  0.1217;        // g/cm^3
	const double TargetWidth   = 29.5;           // cm
	const double atomicMass    =  4.0026;        // atomic mass of Beryllium [g/mol]
	const double AvogNumber    =  6.0221409e+23; // Avogadro's number [atoms/mol]
	
	
	
	// Get the pair production cs from file:
	
	vector<double> photon_energy_vec, pair_cs_vec;
	ifstream inf( "He_pair_cs.dat" );
	double a, b, c;
	while( inf >> a >> b >> c ) {
		photon_energy_vec.push_back(a);
		pair_cs_vec.push_back(b+c);
	}
	inf.close();
	
	int n_ens = (int)photon_energy_vec.size();
	double *photon_energy = new double[n_ens];
	double *pair_cs       = new double[n_ens];
	for( int i=0; i<n_ens; i++ ) {
		photon_energy[i] = photon_energy_vec[i] / 1.e3;
		pair_cs[i]       = pair_cs_vec[i];
	}
	
	TGraph *gPairCS = new TGraph( n_ens, photon_energy, pair_cs );
	gPairCS->SetTitle( "He Photon Attenuation Coefficient" );
	gPairCS->GetXaxis()->SetTitle( "E_{#gamma} [GeV]" );
	gPairCS->GetYaxis()->SetTitle( "#mu [cm^{2}/g]" );
	gPairCS->SetMarkerColor(kBlue);
	gPairCS->SetMarkerStyle(8);
	
	TCanvas *canvas = new TCanvas( "canvas", "canvas", 800, 400 );
	canvas->SetTickx(); canvas->SetTicky();
	gPairCS->Draw("AP");
	
	TF1 *fPairCS = new TF1( "fPairCS", "pol4", 3.0, 12.0 );
	gPairCS->Fit( "fPairCS", "R0Q" );
	fPairCS->Draw("same");
	
	double *pair_cs_high  = new double[n_ens];
	double *pair_cs_low   = new double[n_ens];
	for( int i=0; i<n_ens; i++ ) {
		pair_cs_high[i]  = 1.05*pair_cs_vec[i];
		pair_cs_low[i]   = 0.95*pair_cs_vec[i];
	}
	
	TGraph *gPairCS_high = new TGraph( n_ens, photon_energy, pair_cs_high );
	gPairCS_high->SetTitle( "He Photon Attenuation Coefficient" );
	gPairCS_high->GetXaxis()->SetTitle( "E_{#gamma} [GeV]" );
	gPairCS_high->GetYaxis()->SetTitle( "#mu [cm^{2}/g]" );
	gPairCS_high->SetMarkerColor(kGreen);
	gPairCS_high->SetMarkerStyle(8);
	
	TF1 *fPairCS_high = new TF1( "fPairCS_high", "pol4", 3.0, 12.0 );
	gPairCS_high->Fit( "fPairCS_high", "R0Q" );
	
	TGraph *gPairCS_low = new TGraph( n_ens, photon_energy, pair_cs_low );
	gPairCS_low->SetTitle( "He Photon Attenuation Coefficient" );
	gPairCS_low->GetXaxis()->SetTitle( "E_{#gamma} [GeV]" );
	gPairCS_low->GetYaxis()->SetTitle( "#mu [cm^{2}/g]" );
	gPairCS_low->SetMarkerColor(kGreen);
	gPairCS_low->SetMarkerStyle(8);
	
	TF1 *fPairCS_low = new TF1( "fPairCS_low", "pol4", 3.0, 12.0 );
	gPairCS_low->Fit( "fPairCS_low", "R0Q" );
	
	canvas->cd();
	
	gPairCS_high->Draw("P same"); fPairCS_high->Draw("same");
	gPairCS_low->Draw("P same");  fPairCS_low->Draw("same");
	
	canvas->Update();
	
	
	//-----   Loop over TAGH Counters   -----//
	
	TRandom3 *rand = new TRandom3(0);
	
	for( int tagh_counter = 1; tagh_counter <= 274; tagh_counter++ ) 
	{
		
		double loc_eb = tagh_en[tagh_counter-1];
		double loc_cs = fPairCS->Eval(loc_eb);
		
		double n_a         = (AvogNumber/atomicMass)*TargetDensity; // [atoms/cm^3]
		
		double atten_coeff = n_a*loc_cs*(1.e-24) / TargetDensity; // attenuation coefficient
		double rl          = pow( atten_coeff * TargetDensity, -1.0 );
		
		double rl_high = pow(n_a*fPairCS_high->Eval(loc_eb)*(1.e-24), -1.0);
		double rl_low  = pow(n_a* fPairCS_low->Eval(loc_eb)*(1.e-24), -1.0);
		
		double abs_counter = 0., total_counter = 0.;
		double abs_counter_high = 0.;
		double abs_counter_low  = 0.;
		
		for( int ievt = 0; ievt < 1.e6; ievt++ ) {
			
			// for now assume Compton interaction vertex is uniformly distributed in z:
			
			double locVertex = TargetWidth * rand->Uniform();
			
			double a_E      = TMath::Exp(-locVertex/rl     );
			double a_E_high = TMath::Exp(-locVertex/rl_high);
			double a_E_low  = TMath::Exp(-locVertex/rl_low );
			
			// 1-a_E is the probability of pair production occuring before this point
			
			double loc_rand = rand->Uniform();
			
			int abs_val      = 0;
			if( loc_rand < a_E ) abs_val = 1;
			abs_counter      += (double)abs_val;
			
			int abs_val_high = 0;
			if( loc_rand < a_E_high ) abs_val_high = 1;
			abs_counter_high += (double)abs_val_high;
			
			int abs_val_low  = 0;
			if( loc_rand < a_E_low ) abs_val_low = 1;
			abs_counter_low  += (double)abs_val_low;
			
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
	
	for( int tagm_counter = 1; tagm_counter <= 102; tagm_counter++ ) 
	{
		
		double loc_eb = tagm_en[tagm_counter-1];
		double loc_cs = fPairCS->Eval(loc_eb);
		
		double n_a         = (AvogNumber/atomicMass)*TargetDensity; // [atoms/cm^3]
		
		double atten_coeff = n_a*loc_cs*(1.e-24) / TargetDensity; // attenuation coefficient
		double rl          = pow( atten_coeff * TargetDensity, -1.0 );
		
		double rl_high = pow(n_a*fPairCS_high->Eval(loc_eb)*(1.e-24), -1.0);
		double rl_low  = pow(n_a* fPairCS_low->Eval(loc_eb)*(1.e-24), -1.0);
		
		double abs_counter = 0., total_counter = 0.;
		double abs_counter_high = 0.;
		double abs_counter_low  = 0.;
		
		for( int ievt = 0; ievt < 1.e6; ievt++ ) {
			
			// for now assume Compton interaction vertex is uniformly distributed in z:
			
			double locVertex = TargetWidth * rand->Uniform();
			
			double a_E      = TMath::Exp(-locVertex/rl     );
			double a_E_high = TMath::Exp(-locVertex/rl_high);
			double a_E_low  = TMath::Exp(-locVertex/rl_low );
			
			// a_E is the probability of pair production occuring before this point
			
			double loc_rand = rand->Uniform();
			
			int abs_val      = 0;
			if( loc_rand < a_E ) abs_val = 1;
			abs_counter      += (double)abs_val;
			
			int abs_val_high = 0;
			if( loc_rand < a_E_high ) abs_val_high = 1;
			abs_counter_high += (double)abs_val_high;
			
			int abs_val_low  = 0;
			if( loc_rand < a_E_low ) abs_val_low = 1;
			abs_counter_low  += (double)abs_val_low;
			
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
