
#include "compton.h"


void set_cuts() {
	
	FCAL_ENERGY_CUT  =  0.5;
	CCAL_ENERGY_CUT  =  1.0;
	
	FCAL_RF_TIME_CUT =  3.0;
	CCAL_RF_TIME_CUT =  3.0;
	
	deltaE_cut_sig   =  5.0;
	deltaPhi_cut_sig =  5.0;
	deltaK_cut_sig   =  5.0;
	
	cut_energies[0]  =  0.0;
	cut_energies[1]  =  0.5;
	cut_energies[2]  =  1.0;
	cut_energies[3]  =  1.5;
	cut_energies[4]  =  2.0;
	cut_energies[5]  =  2.5;
	cut_energies[6]  =  3.0;
	cut_energies[7]  =  3.5;
	cut_energies[8]  =  4.0;
	cut_energies[9]  =  4.5;
	cut_energies[10] =  5.0;
	cut_energies[11] =  5.5;
	cut_energies[12] =  6.0;
	
	
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
		"/work/halld/home/andrsmit/primex_compton_analysis/analyze_trees/ccalE_cut/rootFiles" );
	
	// Directory where event selection cuts are stored:
	
	sprintf( cut_pathName, 
		"/work/halld/home/andrsmit/primex_compton_analysis/cuts" );
	
	
	
	vector<int> Be_runs       = {61321, 61322, 61323, 61325, 61327, 61329, 61330, 61331, 61332, 
	61333, 61334, 61335, 61336, 61337, 61340, 61341, 61343, 61344};
	vector<int> Be_empty_runs = {61345, 61346, 61347, 61348, 61349, 61350, 61352, 61353, 61354};
	
	vector<int> He_50_runs    = {61914, 61915, 61916, 61917, 61918};
	vector<int> He_100_runs   = {61914, 61915, 61916, 61917, 61918};
	vector<int> He_empty_runs = {61800};
	
	
	
	
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
	
	for( int ic = 0; ic < n_candidates; ic++ ) {
		
		//-----   Minimum Energy Cuts   ----//
		
		int fcal_e_cut = 0;
		if( fcal_e[ic] > FCAL_ENERGY_CUT ) fcal_e_cut = 1;
		
		
		//-----   Timing Cuts   -----//
		
		int fcal_t_cut = 0, ccal_t_cut = 0;
		if( (fcal_t[ic] - rfTime) < FCAL_RF_TIME_CUT ) fcal_t_cut = 1;
		if( (ccal_t[ic] - rfTime) < CCAL_RF_TIME_CUT ) ccal_t_cut = 1;
		
		
		if( !fcal_e_cut ) continue;
		if( !fcal_t_cut || !ccal_t_cut ) continue;
		
		
		double fill_weight = 0;
		if( bunch_val[ic] ) {
			if( (tb[ic]-rfTime < 2.004) ) fill_weight = 1.0;
		} else fill_weight = -0.5;
		
		
		
		//-----   Compton Cuts   -----//
		
		double   deltaE_mu,   deltaE_sig;
		double deltaPhi_mu, deltaPhi_sig;
		double   deltaK_mu,   deltaK_sig;
		
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
		
		int e_cut = 0, p_cut = 0, k_cut = 0;
		if( fabs(deltaE[ic]   -   deltaE_mu) < (deltaE_cut_sig   *   deltaE_sig) ) e_cut = 1;
		if( fabs(deltaPhi[ic] - deltaPhi_mu) < (deltaPhi_cut_sig * deltaPhi_sig) ) p_cut = 1;
		if( fabs(deltaK[ic]   -   deltaK_mu) < (deltaK_cut_sig   *   deltaK_sig) ) k_cut = 1;
		
		if( !e_cut || !p_cut || !k_cut ) continue;
		
		
		for( int icut = 0; icut < N_CUTS; icut++ ) {
			
			int ccal_e_cut = 0;
			if( ccal_e[ic] > cut_energies[icut] ) ccal_e_cut = 1;
			
			if( ccal_e_cut ) {
				
				if( tag_sys[ic]==0 ) {
					h_tagh[icut][tag_counter[ic]-1]->Fill( 
						deltaK[ic], fill_weight );
				} else {
					h_tagm[icut][tag_counter[ic]-1]->Fill( 
						deltaK[ic], fill_weight );
				}
				
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
	
	for( int icut = 0; icut < N_CUTS; icut++ ) {
		for( int tagh_counter=1; tagh_counter<=N_TAGH_COUNTERS; tagh_counter++ )
			h_tagh[icut][tagh_counter-1]->Reset();
		
		for( int tagm_counter=1; tagm_counter<=N_TAGM_COUNTERS; tagm_counter++ )
			h_tagm[icut][tagm_counter-1]->Reset();
		
	}
	
	return;
}


void write_histograms(char *name) {
	
	TFile *fOut = new TFile(name,"RECREATE");
	fOut->cd();
	
	TDirectory *dir_sigs[N_CUTS];
	for( int icut = 0; icut < N_CUTS; icut++ ) {
		
		int energy_cut_int = (int)(cut_energies[icut]*100);
		
		dir_sigs[icut] = new TDirectoryFile( Form("%02dfcalE",energy_cut_int), 
			Form("%02dfcalE",energy_cut_int) );
		dir_sigs[icut]->cd();
		
		
		for( int tagh_counter=1; tagh_counter<=N_TAGH_COUNTERS; tagh_counter++ )
			h_tagh[icut][tagh_counter-1]->Write();
		
		for( int tagm_counter=1; tagm_counter<=N_TAGM_COUNTERS; tagm_counter++ )
			h_tagm[icut][tagm_counter-1]->Write();
		
		dir_sigs[icut]->cd("../");
		
	}
	fOut->Write();
	
	return;
}


void init_histograms() 
{
	
	for( int icut = 0; icut < N_CUTS; icut++ ) {
		
		int cut_energy_int = (int)(cut_energies[icut]*100);
		
		for( int tagh_counter=1; tagh_counter<=N_TAGH_COUNTERS; tagh_counter++ ) {
			h_tagh[icut][tagh_counter-1] = new TH1F( Form("tagh_%03d_%02dfcalE",
				tagh_counter, cut_energy_int), 
				Form("TAGH Counter %d", tagh_counter), 2000, -4.0, 4.0 );
		}
		for( int tagm_counter=1; tagm_counter<=N_TAGM_COUNTERS; tagm_counter++ ) {
			h_tagm[icut][tagm_counter-1] = new TH1F( Form("tagm_%03d_%02dfcalE",
				tagm_counter, cut_energy_int), 
				Form("TAGM Counter %d", tagm_counter), 2000, -4.0, 4.0 );
		}
		
	}
	
	
	return;
}


void load_constants( int group ) 
{
	
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
	
	
	
	return;
}




