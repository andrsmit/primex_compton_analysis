
#include "compton.h"


void set_cuts() {
	
	FCAL_ENERGY_CUT  =  0.5;
	CCAL_ENERGY_CUT  =  1.0;
	
	FCAL_RF_TIME_CUT =  3.0;
	CCAL_RF_TIME_CUT =  3.0;
	
	deltaE_cut_sig   =  5.0;
	deltaPhi_cut_sig =  5.0;
	deltaK_cut_sig   =  6.0;
	
	cut_sigmas[0]    =  1.0;
	cut_sigmas[1]    =  2.0;
	cut_sigmas[2]    =  3.0;
	cut_sigmas[3]    =  3.5;
	cut_sigmas[4]    =  4.0;
	cut_sigmas[5]    =  4.5;
	cut_sigmas[6]    =  5.0;
	cut_sigmas[7]    =  5.5;
	cut_sigmas[8]    =  6.0;
	cut_sigmas[9]    =  6.5;
	cut_sigmas[10]   =  7.0;
	cut_sigmas[11]   =  8.0;
	cut_sigmas[12]   =  9.0;
	cut_sigmas[13]   = 10.0;
	
	
	return;
}



int main(int argc, char **argv) {
	
	
	
	if(argc != 3) {
		cout << "Wrong number of arguments. Example format:  ./ana_compton <start_run> ";
		cout << " <end_run>" << endl;
		return 0;
	}
	
	int startRun = atoi(argv[1]);
	int endRun   = atoi(argv[2]);
	
	cout << "Processing runs " << startRun << " - " << endRun << endl;
	
	
	set_cuts();
	
	
	// Directory where ROOT Trees are stored:
	
	sprintf( rootTree_pathName, 
		"/work/halld/home/andrsmit/primex_compton_analysis/data/rootTrees" );
	
	// Directory where output ROOT files will be stored:
	
	sprintf( rootFile_pathName, 
		"/work/halld/home/andrsmit/primex_compton_analysis/analyze_trees/He50/deltaE_cut/rootFiles" );
	
	// Directory where event selection cuts are stored:
	
	sprintf( cut_pathName, 
		"/work/halld/home/andrsmit/primex_compton_analysis/cuts" );
	
	
	
	vector<int> Be_runs       = {61321, 61322, 61323, 61325, 61327, 61329, 61330, 61331, 61332, 
	 61333, 61334, 61335, 61336, 61337, 61340, 61341, 61343, 61344};
	vector<int> Be_empty_runs = {61345, 61346, 61347, 61348, 61349, 61350, 61352, 61353, 61354};
	
	vector<int> He_50_runs    = {61914, 61915, 61916, 61917, 61918, 61930, 61931, 61933, 61934,
	61935, 61936, 61937, 61938, 61939};
	vector<int> He_100_runs   = {61947, 61950, 61951, 61952, 61953, 61954, 61955, 61956};
	vector<int> He_empty_runs = {61852, 61854, 61855, 61856, 61857, 61858, 61859, 61860, 61861, 
	 61862, 61863};
	
	
	
	
	init_histograms();
	
	for(unsigned int irun = startRun; irun <= endRun; irun++) {
		
		
		reset_histograms();
		
		
		//-----   Check which group of runs this run belongs to   -----//
		
		int run_group = 0;
		for( int jr = 0; jr < (int)Be_runs.size(); jr++ ) {
			if( irun == Be_runs[jr] ) {
				run_group = 1;
				break;
			}
		}
		if( !run_group ) {
			for( int jr = 0; jr < (int)Be_empty_runs.size(); jr++ ) {
				if( irun == Be_empty_runs[jr] ) {
					run_group = 2;
					break;
				}
			}
		}
		if( !run_group ) {
			for( int jr = 0; jr < (int)He_empty_runs.size(); jr++ ) {
				if( irun == He_empty_runs[jr] ) {
					run_group = 3;
					break;
				}
			}
		}
		if( !run_group ) {
			for( int jr = 0; jr < (int)He_50_runs.size(); jr++ ) {
				if( irun == He_50_runs[jr] ) {
					run_group = 4;
					break;
				}
			}
		}
		if( !run_group ) {
			for( int jr = 0; jr < (int)He_100_runs.size(); jr++ ) {
				if( irun == He_100_runs[jr] ) {
					run_group = 5;
					break;
				}
			}
		}
		if( !run_group ) continue;
		
		
		cout << "Processing run " << irun << endl;
		
		
		load_constants( run_group );
		
		for( int iext = 0; iext < 100; iext++ ) {
			
			char buf[256];
			sprintf(buf,"%s/%d/%d_%03d.root", rootTree_pathName, irun, irun, iext);
			if( gSystem->AccessPathName(buf) ) continue;
			
			cout << "  ext " << iext << endl;
			
			TFile *f = new TFile(buf,"READ");
			tree = (TTree*)f->Get("primex_compton");
			int maxEvents = tree->GetEntries();
			first = 1;
			
			for(int ievent = 0; ievent < maxEvents; ++ievent) {
				reset_event();
				read_event(ievent);
				compton_analysis(irun);
			}
			f->Close();
		}
		
		char outFileName[256];
		sprintf(outFileName,"%s/%d_compton.root", rootFile_pathName, irun);
		write_histograms(outFileName);
		
	}
	
	cout << "\n\nFinished processing runs.\n\n" << endl;
	
	return 0;
}



//--------------------------------------------
// Analysis
//--------------------------------------------
void compton_analysis(int run) 
{
	
	
	int good_cands = 0; 
	
	vector<double> ccal_showers_total, fcal_showers_total;
	ccal_showers_total.clear();
	fcal_showers_total.clear();
	
	vector<double> ccal_en_vec[N_SIGS], fcal_en_vec[N_SIGS];
	for( int isig=0; isig<N_SIGS; isig++ ) {
		ccal_en_vec[N_SIGS].clear();
		fcal_en_vec[N_SIGS].clear();
	}
	
	
	for( int ic = 0; ic < n_candidates; ic++ ) {
		
		//-----   Minimum Energy Cuts   ----//
		
		int fcal_e_cut = 0, ccal_e_cut = 0;
		if( fcal_e[ic] > FCAL_ENERGY_CUT ) fcal_e_cut = 1;
		if( ccal_e[ic] > CCAL_ENERGY_CUT ) ccal_e_cut = 1;
		
		//-----   Timing Cuts   -----//
		
		int fcal_t_cut = 0, ccal_t_cut = 0;
		if( (fcal_t[ic] - rfTime) < FCAL_RF_TIME_CUT ) fcal_t_cut = 1;
		if( (ccal_t[ic] - rfTime) < CCAL_RF_TIME_CUT ) ccal_t_cut = 1;
		
		if( !fcal_e_cut || !ccal_e_cut ) continue;
		if( !fcal_t_cut || !ccal_t_cut ) continue;
		
		
		int find_val = 0;
		for( unsigned int j = 0; j < ccal_showers_total.size(); j++ ) {
			if( ccal_showers_total[j] == ccal_e[ic] ) find_val++;
		}
		if( !find_val ) ccal_showers_total.push_back( ccal_e[ic] );
		
		find_val = 0;
		for( unsigned int j = 0; j < fcal_showers_total.size(); j++ ) {
			if( fcal_showers_total[j] == fcal_e[ic] ) find_val++;
		}
		if( !find_val ) fcal_showers_total.push_back( fcal_e[ic] );
	}
	
	
	int n_ccal_total = (int)ccal_showers_total.size();
	int n_fcal_total = (int)fcal_showers_total.size();
	
	//if( n_ccal_total != 1 || n_fcal_total != 1 ) return;
	
	
	
	
	for( int ic = 0; ic < n_candidates; ic++ ) {
		
		//-----   Minimum Energy Cuts   ----//
		
		int fcal_e_cut = 0, ccal_e_cut = 0;
		if( fcal_e[ic] > FCAL_ENERGY_CUT ) fcal_e_cut = 1;
		if( ccal_e[ic] > CCAL_ENERGY_CUT ) ccal_e_cut = 1;
		
		
		//-----   Timing Cuts   -----//
		
		int fcal_t_cut = 0, ccal_t_cut = 0;
		if( (fcal_t[ic] - rfTime) < FCAL_RF_TIME_CUT ) fcal_t_cut = 1;
		if( (ccal_t[ic] - rfTime) < CCAL_RF_TIME_CUT ) ccal_t_cut = 1;
		
		
		if( !fcal_e_cut || !ccal_e_cut ) continue;
		if( !fcal_t_cut || !ccal_t_cut ) continue;
		
		
		double fill_weight = 0;
		if( bunch_val[ic] ) {
			if( (tb[ic]-rfTime < 2.004) ) fill_weight = 1.0;
		} else fill_weight = -0.5;
		
		
		
		//-----   Compton Cuts   -----//
		
		double   deltaE_mu,   deltaE_sig;
		double deltaPhi_mu, deltaPhi_sig;
		double   deltaK_mu,   deltaK_sig;
		
		deltaE_mu     = f_deltaE_mu->Eval(eb[ic]);
		deltaE_sig    = eb[ic] * f_deltaE_sig->Eval(eb[ic]);
		
		deltaPhi_mu   = f_deltaPhi_mu->Eval(eb[ic]);
		deltaPhi_sig  = f_deltaPhi_sig->Eval(eb[ic]);
		
		deltaK_mu     = f_deltaK_mu->Eval(eb[ic]);
		deltaK_sig    = f_deltaK_sig->Eval(eb[ic]);
		
		/*
		if( tag_sys[ic]==0 ) {
			
			deltaE_mu    = deltaE_mu_tagh[tag_counter[ic]-1];
			deltaE_sig   = deltaE_sig_tagh[tag_counter[ic]-1];
			
			deltaPhi_mu  = deltaPhi_mu_tagh[tag_counter[ic]-1];
			deltaPhi_sig = deltaPhi_sig_tagh[tag_counter[ic]-1];
			
			deltaK_mu    = deltaK_mu_tagh[tag_counter[ic]-1];
			deltaK_sig   = deltaK_sig_tagh[tag_counter[ic]-1];
			
		} else {
			
			deltaE_mu    = deltaE_mu_tagm[tag_counter[ic]-1];
			deltaE_sig   = deltaE_sig_tagm[tag_counter[ic]-1];
			
			deltaPhi_mu  = deltaPhi_mu_tagm[tag_counter[ic]-1];
			deltaPhi_sig = deltaPhi_sig_tagm[tag_counter[ic]-1];
			
			deltaK_mu    = deltaK_mu_tagm[tag_counter[ic]-1];
			deltaK_sig   = deltaK_sig_tagm[tag_counter[ic]-1];
			
		}
		*/
		
		int p_cut = 0, k_cut = 0;
		if( fabs(deltaPhi[ic] - deltaPhi_mu) < (deltaPhi_cut_sig * deltaPhi_sig) ) p_cut = 1;
		if( fabs(deltaK[ic]   -   deltaK_mu) < (deltaK_cut_sig   *   deltaK_sig) ) k_cut = 1;
		
		
		if( !p_cut || !k_cut ) continue;
		
		
		for( int isig = 0; isig < N_SIGS; isig++ ) {
			
			int e_cut = 0;
			if( fabs(deltaE[ic] - deltaE_mu) < (cut_sigmas[isig]*deltaE_sig) ) e_cut = 1;
			
			if( e_cut ) {
				
				
				int find_val = 0;
				for( unsigned int j = 0; j < ccal_en_vec[isig].size(); j++ ) {
					if( ccal_en_vec[isig][j] == ccal_e[ic] ) find_val++;
				}
				if( !find_val ) ccal_en_vec[isig].push_back( ccal_e[ic] );
				
				find_val = 0;
				for( unsigned int j = 0; j < fcal_en_vec[isig].size(); j++ ) {
					if( fcal_en_vec[isig][j] == fcal_e[ic] ) find_val++;
				}
				if( !find_val ) fcal_en_vec[isig].push_back( fcal_e[ic] );
				
				
				if( tag_sys[ic]==0 ) {
					h_tagh[isig][tag_counter[ic]-1]->Fill( 
						deltaK[ic], fill_weight );
				} else {
					h_tagm[isig][tag_counter[ic]-1]->Fill( 
						deltaK[ic], fill_weight );
				}
				
			}
		}
		
	}
	
	double ccal_en_total = 0.;
	for( int i=0; i<n_ccal_total; i++ ) {
		ccal_en_total += ccal_showers_total[i];
	}
	
	double fcal_en_total = 0.;
	for( int i=0; i<n_fcal_total; i++ ) {
		fcal_en_total += fcal_showers_total[i];
	}
	
	h_n_ccal_total->Fill( n_ccal_total );
	h_n_fcal_total->Fill( n_fcal_total );
	h_n_show_total->Fill( n_ccal_total+n_fcal_total );	
	
	
	for( int isig = 0; isig < N_SIGS; isig++ ) {
		
		int n_ccal_loc = (int)ccal_en_vec[isig].size();
		double ccal_en_loc = 0.;
		for( int i=0; i<n_ccal_loc; i++ ) {
			ccal_en_loc += ccal_en_vec[isig][i];
		}
		
		int n_fcal_loc = (int)fcal_en_vec[isig].size();
		double fcal_en_loc = 0.;
		for( int i=0; i<n_fcal_loc; i++ ) {
			fcal_en_loc += fcal_en_vec[isig][i];
		}
		
		h_n_ccal[isig]->Fill( n_ccal_loc );
		h_n_fcal[isig]->Fill( n_fcal_loc );
		h_n_show[isig]->Fill( n_ccal_loc+n_fcal_loc );
		
		for( int i=0; i<n_ccal_loc; i++ ) 
			h_ccal_extraE[isig]->Fill( ccal_en_loc-ccal_en_vec[isig][i] );
		for( int i=0; i<n_fcal_loc; i++ ) 
			h_fcal_extraE[isig]->Fill( fcal_en_loc-fcal_en_vec[isig][i] );
		
		for( int i=0; i<n_ccal_loc; i++ ) {
			for( int j=0; j<n_fcal_loc; j++ ) {
				double loc_energy = ccal_en_vec[isig][i] + fcal_en_vec[isig][j];
				h_show_extraE[isig]->Fill( ccal_en_loc+fcal_en_loc - loc_energy );
			}
		}
		
	}
	
	
	
	
	
	
	return;
}




//--------------------------------------------
// Reset Event
//--------------------------------------------
void reset_event() 
{
	n_candidates = 0;
	
	return;
}



//--------------------------------------------
// Read Event Data from Tree
//--------------------------------------------
void read_event(int evt) 
{
	
	// Run and Event number:
	
	if( first ) {
		
		first = 0;	  
		
		// Set branch addresses
		
		tree->SetBranchAddress("rfTime", &rfTime);
		tree->SetBranchAddress("nComp", &n_candidates);
		
		tree->SetBranchAddress("tag_counter", &tag_counter);
		tree->SetBranchAddress("tag_sys",     &tag_sys    );
		tree->SetBranchAddress("bunch_val",   &bunch_val  );
		tree->SetBranchAddress("eb", &eb);
		tree->SetBranchAddress("tb", &tb);
		
		tree->SetBranchAddress("e1", &fcal_e);
		tree->SetBranchAddress("x1", &fcal_x);
		tree->SetBranchAddress("y1", &fcal_y);
		tree->SetBranchAddress("z1", &fcal_z);
		tree->SetBranchAddress("t1", &fcal_t);
		
		tree->SetBranchAddress("e2", &ccal_e);
		tree->SetBranchAddress("x2", &ccal_x);
		tree->SetBranchAddress("y2", &ccal_y);
		tree->SetBranchAddress("z2", &ccal_z);
		tree->SetBranchAddress("t2", &ccal_t);
		
		tree->SetBranchAddress("DeltaE",   &deltaE  );
		tree->SetBranchAddress("DeltaPhi", &deltaPhi);
		tree->SetBranchAddress("DeltaK",   &deltaK  );
		tree->SetBranchAddress("DeltaK2",  &deltaK2 );
		tree->SetBranchAddress("DeltaT",   &deltaT  );
		tree->SetBranchAddress("DeltaR",   &deltaR  );
		
	}
	
	tree->GetEvent(evt);
	
	
	return;
}



void reset_histograms() {
	
	h_n_ccal_total->Reset();
	h_n_fcal_total->Reset();
	h_n_show_total->Reset();
	
	for( int isig = 0; isig < N_SIGS; isig++ ) {
		
		h_n_ccal[isig]->Reset();
		h_n_fcal[isig]->Reset();
		h_n_show[isig]->Reset();
		
		h_ccal_extraE[isig]->Reset();
		h_fcal_extraE[isig]->Reset();
		h_show_extraE[isig]->Reset();
		
		for( int tagh_counter=1; tagh_counter<=N_TAGH_COUNTERS; tagh_counter++ )
			h_tagh[isig][tagh_counter-1]->Reset();
		
		for( int tagm_counter=1; tagm_counter<=N_TAGM_COUNTERS; tagm_counter++ )
			h_tagm[isig][tagm_counter-1]->Reset();
		
	}
	
	return;
}


void write_histograms(char *name) {
	
	TFile *fOut = new TFile(name,"RECREATE");
	fOut->cd();
	
	h_n_ccal_total->Write();
	h_n_fcal_total->Write();
	h_n_show_total->Write();
	
	TDirectory *dir_sigs[N_SIGS];
	for( int isig = 0; isig < N_SIGS; isig++ ) {
		
		int sig_int = (int)(cut_sigmas[isig]*10);
		
		dir_sigs[isig] = new TDirectoryFile( Form("%03dsigE",sig_int), 
			Form("%03dsigE",sig_int) );
		dir_sigs[isig]->cd();
		
		h_n_ccal[isig]->Write();
		h_n_fcal[isig]->Write();
		h_n_show[isig]->Write();
		
		h_ccal_extraE[isig]->Write();
		h_fcal_extraE[isig]->Write();
		h_show_extraE[isig]->Write();
		
		for( int tagh_counter=1; tagh_counter<=N_TAGH_COUNTERS; tagh_counter++ )
			h_tagh[isig][tagh_counter-1]->Write();
		
		for( int tagm_counter=1; tagm_counter<=N_TAGM_COUNTERS; tagm_counter++ )
			h_tagm[isig][tagm_counter-1]->Write();
		
		dir_sigs[isig]->cd("../");
		
	}
	fOut->Write();
	
	return;
}


void init_histograms() 
{
	
	h_n_ccal_total = new TH1F( "n_ccal_total", "Number of CCAL Showers", 10, -0.5, 9.5 );
	h_n_fcal_total = new TH1F( "n_fcal_total", "Number of FCAL Showers", 10, -0.5, 9.5 );
	h_n_show_total = new TH1F( "n_show_total", "Number of Showers",      10, -0.5, 9.5 );
	
	for( int isig = 0; isig < N_SIGS; isig++ ) {
		
		int sig_int = (int)(cut_sigmas[isig]*10);
		
		h_n_ccal[isig] = new TH1F( Form("n_ccal_%03dsigE",sig_int), "Number of CCAL Showers", 
			10, -0.5, 9.5 );
		h_n_fcal[isig] = new TH1F( Form("n_fcal_%03dsigE",sig_int), "Number of FCAL Showers", 
			10, -0.5, 9.5 );
		h_n_show[isig] = new TH1F( Form("n_show_%03dsigE",sig_int), "Number of Showers", 
			10, -0.5, 9.5 );
		
		h_ccal_extraE[isig] = new TH1F( Form("ccal_extraE_%03dsigE",sig_int), 
			"Extra CCAL Energy; [GeV]", 1200, 0., 12. );
		h_fcal_extraE[isig] = new TH1F( Form("fcal_extraE_%03dsigE",sig_int), 
			"Extra FCAL Energy; [GeV]", 1200, 0., 12. );
		h_show_extraE[isig] = new TH1F( Form("show_extraE_%03dsigE",sig_int), 
			"Extra Shower Energy; [GeV]", 1200, 0., 12. );
		
		
		for( int tagh_counter=1; tagh_counter<=N_TAGH_COUNTERS; tagh_counter++ ) {
			h_tagh[isig][tagh_counter-1] = new TH1F( Form("tagh_%03d_%03dsigE",
				tagh_counter, sig_int), 
				Form("TAGH Counter %d", tagh_counter), 2000, -4.0, 4.0 );
		}
		for( int tagm_counter=1; tagm_counter<=N_TAGM_COUNTERS; tagm_counter++ ) {
			h_tagm[isig][tagm_counter-1] = new TH1F( Form("tagm_%03d_%03dsigE",
				tagm_counter, sig_int), 
				Form("TAGM Counter %d", tagm_counter), 2000, -4.0, 4.0 );
		}
		
	}
	
	
	return;
}


void load_constants( int group ) 
{
	
	
	f_deltaE_mu = new TF1(   "f_deltaE_mu", "pol3", 5.0, 12.0 );
	f_deltaE_mu->SetParameters( deltaE_mu_p0, deltaE_mu_p1, deltaE_mu_p2, deltaE_mu_p3 );
	
	f_deltaE_sig = new TF1( "f_deltaE_sig", "[0] + [1]/sqrt(x) + [2]/x", 
		5.0, 12.0 );
	f_deltaE_sig->SetParameters( deltaE_sig_p0, deltaE_sig_p1, deltaE_sig_p2 );
	
	
	f_deltaPhi_mu = new TF1(   "f_deltaPhi_mu", "pol3", 5.0, 12.0 );
	f_deltaPhi_mu->SetParameters(   deltaPhi_mu_p0,  deltaPhi_mu_p1,  
		deltaPhi_mu_p2,  deltaPhi_mu_p3 );
	
	f_deltaPhi_sig = new TF1( "f_deltaPhi_sig", "pol3", 5.0, 12.0 );
	f_deltaPhi_sig->SetParameters( deltaPhi_sig_p0, deltaPhi_sig_p1, 
		deltaPhi_sig_p2, deltaPhi_sig_p3 );
	
	
	f_deltaK_mu = new TF1(   "f_deltaK_mu", "pol3", 5.0, 12.0 );
	f_deltaK_mu->SetParameters(   deltaK_mu_p0,  deltaK_mu_p1,  deltaK_mu_p2,  deltaK_mu_p3 );
	
	f_deltaK_sig = new TF1( "f_deltaK_sig", "pol3", 5.0, 12.0 );
	f_deltaK_sig->SetParameters( deltaK_sig_p0, deltaK_sig_p1, deltaK_sig_p2, deltaK_sig_p3 );
	
	
	/*
	char loc_pathName[256];
	
	if( group < 3 )      sprintf( loc_pathName, "%s/Be", cut_pathName );
	else if( group < 5 ) sprintf( loc_pathName, "%s/He50", cut_pathName );
	else                 sprintf( loc_pathName, "%s/He100", cut_pathName );
	
	for( int tagh_counter=1; tagh_counter<N_TAGH_COUNTERS; tagh_counter++ ) {
		deltaE_mu_tagh[tagh_counter-1]    = 0.;
		deltaE_sig_tagh[tagh_counter-1]   = 0.;
		
		deltaPhi_mu_tagh[tagh_counter-1]  = 0.;
		deltaPhi_sig_tagh[tagh_counter-1] = 0.;
		
		deltaK_mu_tagh[tagh_counter-1]    = 0.;
		deltaK_sig_tagh[tagh_counter-1]   = 0.;
	}
	for( int tagm_counter=1; tagm_counter<N_TAGM_COUNTERS; tagm_counter++ ) {
		deltaE_mu_tagm[tagm_counter-1]    = 0.;
		deltaE_sig_tagm[tagm_counter-1]   = 0.;
		
		deltaPhi_mu_tagm[tagm_counter-1]  = 0.;
		deltaPhi_sig_tagm[tagm_counter-1] = 0.;
		
		deltaK_mu_tagm[tagm_counter-1]    = 0.;
		deltaK_sig_tagm[tagm_counter-1]   = 0.;
	}
	
	
	char fname[256];
	
	//-----   DeltaE   -----//
	
	sprintf( fname, "%s/deltaE_tagh.dat", loc_pathName );
	ifstream infile_deltaE_tagh(fname);
	if( infile_deltaE_tagh.good() ) {
		int a; double b, c;
		for( int i=0; i<N_TAGH_COUNTERS; i++ ) {
			infile_deltaE_tagh >> a >> b >> c;
			deltaE_mu_tagh[i]  = b;
			deltaE_sig_tagh[i] = c;
		}
	}
	infile_deltaE_tagh.close();
	
	sprintf( fname, "%s/deltaE_tagm.dat", loc_pathName );
	ifstream infile_deltaE_tagm(fname);
	if( infile_deltaE_tagm.good() ) {
		int a; double b, c;
		for( int i=0; i<N_TAGM_COUNTERS; i++ ) {
			infile_deltaE_tagm >> a >> b >> c;
			deltaE_mu_tagm[i]  = b;
			deltaE_sig_tagm[i] = c;
		}
	}
	infile_deltaE_tagm.close();
	
	
	
	//-----   DeltaPhi   -----//
	
	sprintf( fname, "%s/deltaPhi_tagh.dat", loc_pathName );
	ifstream infile_deltaPhi_tagh(fname);
	if( infile_deltaPhi_tagh.good() ) {
		int a; double b, c;
		for( int i=0; i<N_TAGH_COUNTERS; i++ ) {
			infile_deltaPhi_tagh >> a >> b >> c;
			deltaPhi_mu_tagh[i]  = b;
			deltaPhi_sig_tagh[i] = c;
		}
	}
	infile_deltaPhi_tagh.close();
	
	sprintf( fname, "%s/deltaPhi_tagm.dat", loc_pathName );
	ifstream infile_deltaPhi_tagm(fname);
	if( infile_deltaPhi_tagm.good() ) {
		int a; double b, c;
		for( int i=0; i<N_TAGM_COUNTERS; i++ ) {
			infile_deltaPhi_tagm >> a >> b >> c;
			deltaPhi_mu_tagm[i]  = b;
			deltaPhi_sig_tagm[i] = c;
		}
	}
	infile_deltaPhi_tagm.close();
	
	
	
	//-----   DeltaK   -----//
	
	sprintf( fname, "%s/deltaK_tagh.dat", loc_pathName );
	ifstream infile_deltaK_tagh(fname);
	if( infile_deltaK_tagh.good() ) {
		int a; double b, c;
		for( int i=0; i<N_TAGH_COUNTERS; i++ ) {
			infile_deltaK_tagh >> a >> b >> c;
			deltaK_mu_tagh[i]  = b;
			deltaK_sig_tagh[i] = c;
		}
	}
	infile_deltaK_tagh.close();
	
	sprintf( fname, "%s/deltaK_tagm.dat", loc_pathName );
	ifstream infile_deltaK_tagm(fname);
	if( infile_deltaK_tagm.good() ) {
		int a; double b, c;
		for( int i=0; i<N_TAGM_COUNTERS; i++ ) {
			infile_deltaK_tagm >> a >> b >> c;
			deltaK_mu_tagm[i]  = b;
			deltaK_sig_tagm[i] = c;
		}
	}
	infile_deltaK_tagm.close();
	*/
	
	
	return;
}




