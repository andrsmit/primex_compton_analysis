
void add_flux()
{
	
	/*
	vector<int> run_list = {61321, 61322, 61323, 61325, 61327, 61329, 61330, 61331, 61332, 61333, 
	61334, 61335, 61336, 61337, 61340, 61341, 61343, 61344};
	
	vector<int> run_list = {61345, 61346, 61347, 61348, 61349, 61350, 61352, 61353, 61354};
	
	vector<int> run_list = {61852, 61854, 61855, 61856, 61857, 61858, 61859, 61860, 61861, 61862};
	
	vector<int> run_list = {61914, 61915, 61916, 61918, 61930, 61931, 61933, 61934, 61935, 
	61936, 61937, 61938, 61939};
	
	vector<int> run_list = {61947, 61950, 61951, 61952, 61953, 61954, 61955, 61956};
	*/
	
	
	
	char outfile_tagh[256], outfile_tagm[256];
	
	sprintf( outfile_tagh, "Be_tagh_flux.txt" );
	sprintf( outfile_tagm, "Be_tagm_flux.txt" );
	vector<int> run_list = {61321, 61322, 61323, 61325, 61327, 61329, 61330, 61331, 61332, 61333, 
	61334, 61335, 61336, 61337, 61340, 61341, 61343, 61344};
	
	
	double tagh_flux[274], tagh_fluxE[274];
	for( int i=0; i<274; i++ ) { tagh_flux[i] = 0.; tagh_fluxE[i] = 0.; }
	
	double tagm_flux[102], tagm_fluxE[102];
	for( int i=0; i<102; i++ ) { tagm_flux[i] = 0.; tagm_fluxE[i] = 0.; }
	
	
	for( int ir = 0; ir < (int)run_list.size(); ir++ ) {
		
		int a; double b, c;
		
		if( gSystem->AccessPathName(Form("txtFiles/%d_tagh_ps_acc_cor.txt",run_list[ir])) )
		{
			cout << "No txt files for run " << run_list[ir] << endl;
			continue;
		}
		
		ifstream inf1( Form("txtFiles/%d_tagh_ps_acc_cor.txt",run_list[ir]) );
		for( int i=0; i<274; i++ ) {
			inf1 >> a >> b >> c;
			tagh_flux[i]  += b;
			tagh_fluxE[i] += c*c;
		}
		inf1.close();
		
		ifstream inf2( Form("txtFiles/%d_tagm_ps_acc_cor.txt",run_list[ir]) );
		for( int i=0; i<102; i++ ) {
			inf2 >> a >> b >> c;
			tagm_flux[i]  += b;
			tagm_fluxE[i] += c*c;
		}
		inf2.close();
		
	}
	
	char buf[256];
	
	ofstream outf1( outfile_tagh );
	for( int i=0; i<274; i++ ) {
		sprintf( buf, "%4d  %10.3f  %10.3f \n", i+1, tagh_flux[i], sqrt(tagh_fluxE[i]) );
		outf1 << buf;
	}
	outf1.close();
	
	ofstream outf2( outfile_tagm );
	for( int i=0; i<102; i++ ) {
		sprintf( buf, "%4d  %10.3f  %10.3f \n", i+1, tagm_flux[i], sqrt(tagm_fluxE[i]) );
		outf2 << buf;
	}
	outf2.close();
	
	
	
	
	return;
}
