void get_phaseI_energy_bin(double eb, int &old_tag_counter, int &old_tag_sys) {
	
	double tagh_en_phaseI[274];
	double tagm_en_phaseI[102];
	
	int a; double b, c;
	
	double loc_endpoint_energy       = 11.6061;
	double loc_endpoint_energy_calib = 11.6061;
	
	ifstream inf1(
		"/work/halld/home/andrsmit/primex_compton_analysis/photon_flux/phaseI/primex_tagh.txt");
	for( int i=0; i<274; i++ ) {
		inf1 >> a >> b >> c;
		double deltaE  = loc_endpoint_energy - loc_endpoint_energy_calib;
		double emin    = b * loc_endpoint_energy_calib  +  deltaE;
		double emax    = c * loc_endpoint_energy_calib  +  deltaE;
		tagh_en_phaseI[i] = 0.5 * (emin + emax);
	}
	inf1.close();
	
	ifstream inf2(
		"/work/halld/home/andrsmit/primex_compton_analysis/photon_flux/phaseI/primex_tagm.txt");
	for( int i=0; i<102; i++ ) {
		inf2 >> a >> b >> c;
		double deltaE  = loc_endpoint_energy - loc_endpoint_energy_calib;
		double emin    = b * loc_endpoint_energy_calib  +  deltaE;
		double emax    = c * loc_endpoint_energy_calib  +  deltaE;
		tagm_en_phaseI[i] = 0.5 * (emin + emax);
	}
	inf2.close();
	
	double loc_min_e = 12.0;
	old_tag_sys = -1;
	old_tag_counter = 0;
	
	double old_eb = 0.;
	
	for(int tagh_counter=1; tagh_counter<=274; tagh_counter++) {
		if(tagh_counter>126 && tagh_counter<179) continue;
		double loc_deltaE = fabs(tagh_en_phaseI[tagh_counter-1] - eb);
		if(loc_deltaE < loc_min_e) {
			loc_min_e = loc_deltaE;
			old_tag_sys = 0;
			old_tag_counter = tagh_counter;
			old_eb = tagh_en_phaseI[tagh_counter-1];
		}
	}
	for(int tagm_counter=1; tagm_counter<=102; tagm_counter++) {
		double loc_deltaE = fabs(tagm_en_phaseI[tagm_counter-1] - eb);
		if(loc_deltaE < loc_min_e) {
			loc_min_e = loc_deltaE;
			old_tag_sys = 1;
			old_tag_counter = tagm_counter;
			old_eb = tagm_en_phaseI[tagm_counter-1];
		}
	}
	
	return;
}
