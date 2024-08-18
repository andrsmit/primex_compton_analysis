/******************************************************************************

Combine the output txt files from 'primex_flux.py' for a group of runs into 
a single file.

*******************************************************************************/

/*
	sprintf(outfile_tagh, "Be_200nA_FIELDOFF_flux_tagh.txt");
	sprintf(outfile_tagm, "Be_200nA_FIELDOFF_flux_tagm.txt");
	vector<int> run_list = {61321, 61322, 61323, 61325, 61327, 61329, 61330, 61331, 61332, 
	61333, 61334, 61335, 61336, 61337, 61340, 61341, 61343, 61344};
	
	sprintf(outfile_tagh, "Be_empty_FIELDOFF_flux_tagh.txt");
	sprintf(outfile_tagm, "Be_empty_FIELDOFF_flux_tagm.txt");
	vector<int> run_list = {61345, 61346, 61347, 61348, 61349, 61350, 61352, 61353, 61354};
	
	sprintf(outfile_tagh, "He_empty_FIELDOFF_flux_tagh.txt");
	sprintf(outfile_tagm, "He_empty_FIELDOFF_flux_tagm.txt");
	vector<int> run_list = {61852, 61854, 61855, 61856, 61857, 61858, 61859, 61860, 61861, 
	61862, 61863};
	vector<int> run_list = {61852, 61854, 61858, 61859, 61860, 61861, 61862, 61863};
	
	i dont see anything wrong with these runs:
	61867, 61874
	why were they excluded previously?
	for run 61879, file 61879_001 fails everytime I try to analyze it, but it may also be the 
	case for the photon flux. maybe it can be included then
	
	sprintf(outfile_tagh, "He_200nA_FIELDOFF_flux_tagh.txt");
	sprintf(outfile_tagm, "He_200nA_FIELDOFF_flux_tagm.txt");
	vector<int> run_list = {61866, 61868, 61875, 61876, 61877, 61878, 61880, 61881, 61882, 
	61883, 61884, 61885, 61887, 61888, 61889, 61890, 61891, 61892, 61893, 61894, 61895, 
	61905, 61906, 61908, 61909, 61910};
	vector<int> run_list = {61866, 61867, 61868, 61874, 61875, 61876, 61877, 61878, 61879, 61880, 
	61881, 61882, 61883, 61884, 61885, 61887, 61888, 61889, 61890, 61891, 61892, 61893, 61894, 
	61895, 61905, 61906, 61908, 61909, 61910};
	
	sprintf(outfile_tagh, "He_100nA_FIELDOFF_flux_tagh.txt");
	sprintf(outfile_tagm, "He_100nA_FIELDOFF_flux_tagm.txt");
	vector<int> run_list = {61950, 61951, 61953, 61954, 61955, 61956};
	vector<int> run_list = {61947, 61950, 61951, 61952, 61953, 61954, 61955, 61956};
	(runs 61947 and 61952 have bad timing alignment - 61947 is REALLY bad)
	
	sprintf(outfile_tagh, "He_50nA_FIELDOFF_flux_tagh.txt");
	sprintf(outfile_tagm, "He_50nA_FIELDOFF_flux_tagm.txt");
	vector<int> run_list = {61914, 61915, 61916, 61917, 61918, 61930, 61931, 61933, 61934, 61935, 
	61936, 61937, 61938, 61939};
*/

void add_flux()
{
	
	char pathName[256] = "/work/halld/home/andrsmit/primex_compton_analysis/photon_flux/phase1";
	
	char outfile_tagh[256], outfile_tagm[256];
	
	sprintf(outfile_tagh, "Be_empty_FIELDOFF_flux_tagh.txt");
	sprintf(outfile_tagm, "Be_empty_FIELDOFF_flux_tagm.txt");
	vector<int> run_list = {61345, 61346, 61347, 61348, 61349, 61350, 61352, 61353, 61354};
	
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
