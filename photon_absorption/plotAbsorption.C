
void plotAbsorption()
{
	
	const double TargetDensity = 1.85; // g/cm^3
	const double TargetWidth   = 1.77546; // cm
	
	
	// Get the photon attenuation coefficient as a function of the beam energy:
	
	vector<double> photon_energy_vec, mu_atten_vec;
	ifstream inf( "Be_atten.dat" );
	double a, b;
	while( inf >> a >> b ) {
		photon_energy_vec.push_back(a);
		mu_atten_vec.push_back(b);
	}
	inf.close();
	
	int n_ens = (int)photon_energy_vec.size();
	double *photon_energy = new double[n_ens];
	double *mu_atten      = new double[n_ens];
	for( int i=0; i<n_ens; i++ ) {
		photon_energy[i] = photon_energy_vec[i] / 1.e3;
		mu_atten[i]      = mu_atten_vec[i];
	}
	
	TGraph *gAtten = new TGraph( n_ens, photon_energy, mu_atten );
	gAtten->SetTitle( "Be Photon Attenuation Coefficient" );
	gAtten->GetXaxis()->SetTitle( "E_{#gamma} [GeV]" );
	gAtten->GetYaxis()->SetTitle( "#mu [cm^{2}/g]" );
	gAtten->SetMarkerColor(kBlue);
	gAtten->SetMarkerStyle(8);
	
	TCanvas *canvas = new TCanvas( "canvas", "canvas", 800, 400 );
	canvas->SetTickx(); canvas->SetTicky();
	gAtten->Draw("AP");
	
	TF1 *fAtten = new TF1( "fAtten", "pol3", 6.0, 12.0 );
	gAtten->Fit( "fAtten", "R0Q" );
	fAtten->Draw("same");
	
	
	/*
	double eGam = 9.0;
	double val  = fAtten->Eval(eGam) * TargetDensity * TargetWidth;
	double a_E  = 1.0 - ((1.0 - TMath::Exp(val))/val);
	
	double absorption_factor = 1.0 - a_E;
	
	cout << "\n\n\n";
	cout << "E_gamma = " << eGam << ", Absorption factor: " << absorption_factor << endl;
	cout << "\n\n\n";
	*/
	
	
	//----------------------------------------------------------------//
	
	TRandom3 *rand = new TRandom3(0);
	
	vector<double> eGamVec, absVec;
	
	for( int ie = 0; ie < 60; ie++ ) {
	
	double eGam = 6.0 + 0.1 * (double)(ie);
	
	double abs_counter = 0.;
	double total_counter = 0.;
	
	for( int ievt = 0; ievt < 1.e7; ievt++ ) {
		
		double locTargetWidth = TargetWidth * rand->Uniform();
		double val            = fAtten->Eval(eGam) * TargetDensity * locTargetWidth;
		//double a_E            = 1.0 - (1.0 - TMath::Exp(-val))/val;
		double a_E            = TMath::Exp(-val);
		
		// a_E is the probability that a photon is absorbed:
		
		int abs_val = 0;
		if( rand->Uniform() > a_E ) {
			abs_val = 1;
		}
		
		abs_counter   += (double)abs_val;
		total_counter += 1.0;
	}
	
	double absorption_factor = 1.0  -  (abs_counter/total_counter);
	
	eGamVec.push_back(eGam);
	absVec.push_back(absorption_factor);
	
	cout << "E_gamma = " << eGam << ", Absorption factor: " << absorption_factor << endl;
	
	}
	
	
	n_ens = (int)eGamVec.size();
	
	double *eGam = new double[n_ens];
	double *abs  = new double[n_ens];
	
	for( int i=0; i<n_ens; i++ ) {
		
		eGam[i] = eGamVec[i];
		abs[i]  = absVec[i];
		
	}
	
	TGraph *graph = new TGraph( n_ens, eGam, abs );
	graph->SetTitle( "Be Target Photon Beam Absorption Correction Factor" );
	graph->GetXaxis()->SetTitle( "E_{#gamma} [GeV]" );
	graph->GetYaxis()->SetTitle( "f_{abs}" );
	graph->SetMarkerStyle(8);
	
	TCanvas *canvas2 = new TCanvas( "canvas2", "canvas2", 800, 500 );
	canvas2->SetTickx(); canvas2->SetTicky();
	graph->Draw("AP");
	
	
	
	
	return;
}
