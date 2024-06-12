/******************************************************************************

Combine the output txt files from 'primex_flux.py' for a group of runs into 
a single file.

*******************************************************************************/

/*
	sprintf(outfile_tagh, "Be_200nA_FIELDOFF_flux_tagh.txt");
	sprintf(outfile_tagm, "Be_200nA_FIELDOFF_flux_tagm.txt");
	vector<int> run_list = {110600, 110601, 110602, 110603, 110605, 110606, 110607, 110608, 
	110609, 110610, 110615};
	
	sprintf(outfile_tagh, "Be_100nA_FIELDOFF_flux_tagh.txt");
	sprintf(outfile_tagm, "Be_100nA_FIELDOFF_flux_tagm.txt");
	vector<int> run_list = {110616, 110617, 110618, 110619, 110621};
	
	sprintf(outfile_tagh, "Be_empty_FIELDOFF_flux_tagh.txt");
	sprintf(outfile_tagm, "Be_empty_FIELDOFF_flux_tagm.txt");
	vector<int> run_list = {110477, 110479, 110481, 110482};
	
	sprintf(outfile_tagh, "Be_200nA_FIELDOFF_good_flux_tagh.txt");
	sprintf(outfile_tagm, "Be_200nA_FIELDOFF_good_flux_tagm.txt");
	vector<int> run_list = {110603, 110605, 110607, 110608, 110610, 110615};
*/

void add_flux()
{
	
	char pathName[256] = "/work/halld/home/andrsmit/primex_compton_analysis/photon_flux/phaseIII";
	
	char outfile_tagh[256], outfile_tagm[256];
	
	sprintf(outfile_tagh, "Be_200nA_FIELDOFF_flux_tagh.txt");
	sprintf(outfile_tagm, "Be_200nA_FIELDOFF_flux_tagm.txt");
	vector<int> run_list = {110600, 110601, 110602, 110603, 110605, 110606, 110607, 110608, 
	110609, 110610, 110615};
	
	double tagh_flux[274], tagh_fluxE[274];
	double tagm_flux[102], tagm_fluxE[102];
	
	for(int i=0; i<274; i++) {
		tagh_flux[i]  = 0.;
		tagh_fluxE[i] = 0.;
	}
	for(int i=0; i<102; i++) {
		tagm_flux[i]  = 0.;
		tagm_fluxE[i] = 0.;
	}
	
	for(int ir = 0; ir < (int)run_list.size(); ir++) {
		
		int a; double b, c;
		
		if(gSystem->AccessPathName(Form("%s/txtFiles/%d_tagh_ps_acc_cor.txt", 
			pathName, run_list[ir])))
		{
			cout << "No txt files for run " << run_list[ir] << endl;
			continue;
		}
		if(gSystem->AccessPathName(Form("%s/txtFiles/%d_tagm_ps_acc_cor.txt", 
			pathName, run_list[ir])))
		{
			cout << "No txt files for run " << run_list[ir] << endl;
			continue;
		}
		
		ifstream inf1(Form("%s/txtFiles/%d_tagh_ps_acc_cor.txt", pathName, run_list[ir]));
		for( int i=0; i<274; i++ ) {
			inf1 >> a >> b >> c;
			tagh_flux[i]  += b;
			tagh_fluxE[i] += c*c;
		}
		inf1.close();
		
		ifstream inf2(Form("%s/txtFiles/%d_tagm_ps_acc_cor.txt", pathName, run_list[ir]));
		for( int i=0; i<102; i++ ) {
			inf2 >> a >> b >> c;
			tagm_flux[i]  += b;
			tagm_fluxE[i] += c*c;
		}
		inf2.close();
	}
	
	char buf[256];
	
	ofstream outf1(outfile_tagh);
	for( int i=0; i<274; i++ ) {
		sprintf(buf, "%4d  %10.3f  %10.3f \n", i+1, tagh_flux[i], sqrt(tagh_fluxE[i]));
		outf1 << buf;
	}
	outf1.close();
	
	ofstream outf2(outfile_tagm);
	for( int i=0; i<102; i++ ) {
		sprintf(buf, "%4d  %10.3f  %10.3f \n", i+1, tagm_flux[i], sqrt(tagm_fluxE[i]));
		outf2 << buf;
	}
	outf2.close();
	
	return;
}
