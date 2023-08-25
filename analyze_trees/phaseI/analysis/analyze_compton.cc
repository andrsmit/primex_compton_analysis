
#include "compton.h"

void set_cuts() {
	
	FCAL_ENERGY_CUT  = 0.50;
	CCAL_ENERGY_CUT  = 3.00;
	
	FCAL_RF_TIME_CUT = 3.0;
	CCAL_RF_TIME_CUT = 3.0;
	
	deltaE_cut_sig   = 5.0;
	deltaPhi_cut_sig = 5.0;
	deltaK_cut_sig   = 5.0;
	
	fcal_inner_layer_cut = 2.5*4.0157;
	ccal_inner_layer_cut = 2.0*2.09;
	
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
	
	sprintf(rootTree_pathName, 
	"/work/halld/home/andrsmit/primex_compton_analysis/data/rootTrees/phaseI");
	
	// Directory where output ROOT files will be stored:
	
	sprintf(rootFile_pathName, 
	"/work/halld/home/andrsmit/primex_compton_analysis/analyze_trees/phaseI/analysis/rootFiles");
	
	//--------------------------------------------------------------------------------//
	
	vector<int> Be_empty_runs  = {
	61345, 61346, 61347, 61348, 61349, 61350, 61352, 61353, 61354};
	vector<int> Be_200nA_runs  = {
	61321, 61322, 61323, 61325, 61327, 61329, 61330, 61331, 61332, 61333, 61334, 61335, 
	61336, 61337, 61340, 61341, 61343, 61344};
	
	vector<int> He_empty_runs = {
	61852, 61854, 61855, 61856, 61857, 61858, 61859, 61860, 61861, 61862, 61863};
	vector<int> He_200nA_runs = {
	61866, 61867, 61868, 61874, 61875, 61876, 61877, 61878, 61889, 61880, 61881, 61882, 
	61883, 61884, 61885, 61887, 61888, 61889, 61890, 61891, 61892, 61893, 61894, 61895, 
	61905, 61906, 61908, 61909, 61910};
	vector<int> He_100nA_runs = {
	61947, 61950, 61951, 61952, 61953, 61954, 61955, 61956};
	vector<int> He_050nA_runs = {
	61914, 61915, 61916, 61917, 61918, 61930, 61931, 61933, 61934, 61935, 61936, 61937, 
	61938, 61939};
	
	//--------------------------------------------------------------------------------//
	
	init_histograms();
	
	bool first_evt = true;
	
	for(unsigned int irun = startRun; irun <= endRun; irun++) {
		
		reset_histograms();
		
		//-----   Check which group of runs this run belongs to   -----//
		
		int run_group = 0;
		for(int jr = 0; jr < (int)Be_empty_runs.size(); jr++) {
			if(irun == Be_empty_runs[jr]) {
				run_group = 1;
				break;
			}
		}
		if(!run_group) {
			for(int jr = 0; jr < (int)Be_200nA_runs.size(); jr++) {
				if(irun == Be_200nA_runs[jr]) {
					run_group = 2;
					break;
				}
			}
		}
		if(!run_group) {
			for(int jr = 0; jr < (int)He_empty_runs.size(); jr++) {
				if(irun == He_empty_runs[jr]) {
					run_group = 3;
					break;
				}
			}
		}
		if(!run_group) {
			for(int jr = 0; jr < (int)He_200nA_runs.size(); jr++) {
				if(irun == He_200nA_runs[jr]) {
					run_group = 4;
					break;
				}
			}
		}
		if(!run_group) {
			for(int jr = 0; jr < (int)He_100nA_runs.size(); jr++) {
				if(irun == He_100nA_runs[jr]) {
					run_group = 5;
					break;
				}
			}
		}
		if(!run_group) {
			for(int jr = 0; jr < (int)He_050nA_runs.size(); jr++) {
				if(irun == He_050nA_runs[jr]) {
					run_group = 6;
					break;
				}
			}
		}
		if(!run_group) continue;
		
		if(irun<61483) {
			m_beamX =  0.190;
			m_beamY = -0.074;
		} else if(irun<61774) {
			m_beamX =  0.165;
			m_beamY = -0.024;
		} else {
			m_beamX =  0.202;
			m_beamY = -0.042;
		}
		
		if(run_group<3) m_beamZ = 64.945;
		else            m_beamZ = 65.;
		
		if(irun<61483) {
			m_ccalX = 0.083;
			m_ccalY = 0.148;
		} else {
			m_ccalX = 0.083;
			m_ccalY = 0.119;
		}
		
		cout << "Processing run " << irun << endl;
		
		load_constants(run_group, first_evt);
		first_evt = false;
		
		for(int iext = 0; iext < 100; iext++ ) {
			
			char buf[256];
			sprintf(buf,"%s/%d/%d_%03d.root", rootTree_pathName, irun, irun, iext);
			if( gSystem->AccessPathName(buf) ) continue;
			
			cout << "  ext " << iext << endl;
			
			TFile *f = new TFile(buf,"READ");
			tree = (TTree*)f->Get("primex_compton")->Clone(Form(
				"primex_compton_%d_%d",irun,iext));
			int maxEvents = tree->GetEntries();
			first = 1;
			
			for(int ievent = 0; ievent < maxEvents; ++ievent) {
				reset_event();
				read_event(ievent);
				compton_analysis(irun);
			}
			f->Close();
			//tree->Delete();
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
	if(n_candidates>100) return;
	
	int n_good_cands = 0; 
	
	for(int ic = 0; ic < n_candidates; ic++) {
		
		//----------     Fiducial Cuts     ----------//
		
		double x1 = fcal_x[ic];
		double y1 = fcal_y[ic];
		double z1 = fcal_z[ic];
		
		double fcal_face_x = m_beamX + (x1 * (m_fcalZ - m_beamZ)/z1) - m_fcalX;
		double fcal_face_y = m_beamY + (y1 * (m_fcalZ - m_beamZ)/z1) - m_fcalY;
		
		int fcal_fid_cut = 0;
		if((-1.*fcal_inner_layer_cut < fcal_face_x && fcal_face_x < fcal_inner_layer_cut)
			&& (-1.*fcal_inner_layer_cut < fcal_face_y 
			&& fcal_face_y < fcal_inner_layer_cut)) fcal_fid_cut = 1;
		
		double x2 = ccal_x[ic];
		double y2 = ccal_y[ic];
		double z2 = ccal_z[ic];
		
		double ccal_face_x = m_beamX + (x2 * (m_ccalZ - m_beamZ)/z2) - m_ccalX;
		double ccal_face_y = m_beamY + (y2 * (m_ccalZ - m_beamZ)/z2) - m_ccalY;
		
		int ccal_fid_cut = 0;
		if((-1.*ccal_inner_layer_cut < ccal_face_x && ccal_face_x < ccal_inner_layer_cut)
			&& (-1.*ccal_inner_layer_cut < ccal_face_y 
			&& ccal_face_y < ccal_inner_layer_cut)) ccal_fid_cut = 1;
		
		if(ccal_face_x < -8.36 || ccal_face_x > 10.45 || 
			ccal_face_y < -10.45 || ccal_face_y > 10.45) ccal_fid_cut = 1;
		
		//----------   Minimum Energy Cuts   ---------//
		
		int fcal_e_cut = 0, ccal_e_cut = 0;
		if(fcal_e[ic] > FCAL_ENERGY_CUT) fcal_e_cut = 1;
		if(ccal_e[ic] > CCAL_ENERGY_CUT) ccal_e_cut = 1;
		
		//----------      Timing Cuts      -----------//
		
		int fcal_t_cut = 0, ccal_t_cut = 0;
		if(fabs(fcal_t[ic] - rfTime) < FCAL_RF_TIME_CUT) fcal_t_cut = 1;
		if(fabs(ccal_t[ic] - rfTime) < CCAL_RF_TIME_CUT) ccal_t_cut = 1;
		
		//--------------------------------------------//
		
		h_fcal_rf_dt->Fill(fcal_t[ic]-rfTime);
		h_ccal_rf_dt->Fill(ccal_t[ic]-rfTime);
		
		if(fcal_fid_cut || ccal_fid_cut) continue;
		if(!fcal_e_cut || !ccal_e_cut)   continue;
		if(!fcal_t_cut || !ccal_t_cut)   continue;
		
		// FCAL square cut to remove dead channels:
		/*
		if((-32. < fcal_face_y && fcal_face_y < -20.) && 
			(-8. < fcal_face_x && fcal_face_x < 4.)) continue;
		*/
		
		//--------------------------------------------//
		
		double fill_weight = 0.;
		if(bunch_val[ic]) {
			
			if(run > 61910 && run < 61947) {
				if(fabs(tb[ic]-rfTime+2.0) < 2.004) fill_weight =  1.0;
			} else if( run==61952 ) {
				if(fabs(tb[ic]-rfTime+2.0) < 2.004) fill_weight =  1.0;
			} else {
				if(fabs(tb[ic]-rfTime+0.0) < 2.004) fill_weight =  1.0;
			}
			
		}  else 
			fill_weight = -0.25;
		
		if(run>61910 && run < 61947) {
			h_beam_rf_dt->Fill(tb[ic]-rfTime+2.0);
		} else if(run==61952) {
			h_beam_rf_dt->Fill(tb[ic]-rfTime+2.0);
		} else {
			h_beam_rf_dt->Fill(tb[ic]-rfTime);
		}
		
		//-----   Compton Cuts   -----//
		
		// re-calculate deltaK with new vertex/positions:
		
		double theta1 = atan2(sqrt(x1*x1 + y1*y1),z1);
		double theta2 = atan2(sqrt(x2*x2 + y2*y2),z2);
		
		double ecomp1 = 1. / ((1./eb[ic]) + (1./m_e)*(1.-cos(theta1)));
		double ecomp2 = 1. / ((1./eb[ic]) + (1./m_e)*(1.-cos(theta2)));
		double loc_deltaK = (ecomp1 + ecomp2) - (eb[ic] + m_e);
		
		double ecomp       = m_e * sin(theta1+theta2) 
			/ (sin(theta1) + sin(theta2) - sin(theta1+theta2));
		double loc_deltaK2 = ecomp - eb[ic];
		
		double   deltaE_mu,   deltaE_sig;
		double deltaPhi_mu, deltaPhi_sig;
		double   deltaK_mu,   deltaK_sig;
		double  deltaK2_mu,  deltaK2_sig;
		
		deltaE_mu    = f_deltaE_mu->Eval(eb[ic]);
		deltaE_sig   = eb[ic] * f_deltaE_sig->Eval(eb[ic]);
		deltaPhi_mu  = f_deltaPhi_mu->Eval(eb[ic]);
		deltaPhi_sig = f_deltaPhi_sig->Eval(eb[ic]);
		deltaK_mu    = f_deltaK_mu->Eval(eb[ic]);
		deltaK_sig   = f_deltaK_sig->Eval(eb[ic]);
		deltaK2_mu   = f_deltaK2_mu->Eval(eb[ic]);
		deltaK2_sig  = f_deltaK2_sig->Eval(eb[ic]);
		
		int e_cut = 0, p_cut = 0, k_cut = 0, k2_cut = 0;
		if(fabs(deltaE[ic]   -   deltaE_mu) < (deltaE_cut_sig   *   deltaE_sig))  e_cut = 1;
		if(fabs(deltaPhi[ic] - deltaPhi_mu) < (deltaPhi_cut_sig * deltaPhi_sig))  p_cut = 1;
		if(fabs(loc_deltaK   -   deltaK_mu) < (deltaK_cut_sig   *   deltaK_sig))  k_cut = 1;
		if(fabs(loc_deltaK2  -  deltaK2_mu) < (deltaK_cut_sig   *  deltaK2_sig)) k2_cut = 1;
		
		if(tag_sys[ic]==0) {
			
			h_deltaE_tagh->Fill(tag_counter[ic], deltaE[ic], fill_weight);
			if(e_cut) {
				h_deltaPhi_tagh->Fill(tag_counter[ic], deltaPhi[ic], fill_weight);
				if(p_cut) {
					h_deltaK_tagh->Fill(tag_counter[ic], loc_deltaK, fill_weight);
					h_deltaK2_tagh->Fill(tag_counter[ic], loc_deltaK2, fill_weight);
					if(k_cut) {
						h_deltaK_tagh_cut->Fill(tag_counter[ic], loc_deltaK, fill_weight);
						h_fcal_xy->Fill(fcal_face_x, fcal_face_y);
						h_ccal_xy->Fill(ccal_face_x, ccal_face_y);
						h_beam_rf_dt_cut->Fill(tb[ic]-rfTime);
					}
					if(k2_cut) {
						h_deltaK2_tagh_cut->Fill(tag_counter[ic], loc_deltaK2, fill_weight);
					}
				}
			}
			
		} else {
			
			h_deltaE_tagm->Fill(tag_counter[ic], deltaE[ic], fill_weight);
			if(e_cut) {
				h_deltaPhi_tagm->Fill(tag_counter[ic], deltaPhi[ic], fill_weight);
				if(p_cut) {
					h_deltaK_tagm->Fill(tag_counter[ic], loc_deltaK, fill_weight);
					h_deltaK2_tagm->Fill(tag_counter[ic], loc_deltaK2, fill_weight);
					if(k_cut) {
						h_deltaK_tagm_cut->Fill(tag_counter[ic], loc_deltaK, fill_weight);
						h_fcal_xy->Fill(fcal_face_x, fcal_face_y);
						h_ccal_xy->Fill(ccal_face_x, ccal_face_y);
						h_beam_rf_dt_cut->Fill(tb[ic]-rfTime);
					}
					if(k2_cut) {
						h_deltaK2_tagm_cut->Fill(tag_counter[ic], loc_deltaK2, fill_weight);
					}
				}
			}
		}
		
		h_deltaE_vs_deltaK->Fill(loc_deltaK, deltaE[ic], fill_weight);
		h_deltaE_vs_deltaPhi->Fill(deltaPhi[ic], deltaE[ic], fill_weight);
		h_deltaPhi_vs_deltaK->Fill(loc_deltaK, deltaPhi[ic], fill_weight);
		h_deltaE_vs_deltaK2->Fill(loc_deltaK2, deltaE[ic], fill_weight);
		
		if(bunch_val[ic]) {
			h_deltaE_vs_deltaK2_main->Fill(loc_deltaK2, deltaE[ic], 1.0);
		} else {
			h_deltaE_vs_deltaK2_acc->Fill(loc_deltaK2, deltaE[ic], 1.0);
		}
		
		if(e_cut && p_cut && k_cut) n_good_cands++;
	}
	
	h_n_cands->Fill(n_good_cands);
	
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
	
	if(first) {
		
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
	
	h_fcal_rf_dt->Reset();
	h_ccal_rf_dt->Reset();
	h_beam_rf_dt->Reset();
	h_beam_rf_dt_cut->Reset();
	
	h_n_cands->Reset();
	
	h_deltaE_tagh->Reset();
	h_deltaPhi_tagh->Reset();
	h_deltaK_tagh->Reset();
	h_deltaK_tagh_cut->Reset();
	h_deltaK2_tagh->Reset();
	h_deltaK2_tagh_cut->Reset();
	
	h_deltaE_tagm->Reset();
	h_deltaPhi_tagm->Reset();
	h_deltaK_tagm->Reset();
	h_deltaK_tagm_cut->Reset();
	h_deltaK2_tagm->Reset();
	h_deltaK2_tagm_cut->Reset();
	
	h_fcal_xy->Reset();
	h_ccal_xy->Reset();
	
	h_deltaE_vs_deltaK->Reset();
	h_deltaE_vs_deltaPhi->Reset();
	h_deltaPhi_vs_deltaK->Reset();
	h_deltaE_vs_deltaK2->Reset();
	h_deltaE_vs_deltaK2_main->Reset();
	h_deltaE_vs_deltaK2_acc->Reset();
	
	return;
}


void write_histograms(char *name) {
	
	TFile *fOut = new TFile(name,"RECREATE");
	fOut->cd();
	
	h_fcal_rf_dt->Write();
	h_ccal_rf_dt->Write();
	h_beam_rf_dt->Write();
	h_beam_rf_dt_cut->Write();
	
	h_n_cands->Write();
	
	h_deltaE_tagh->Write();
	h_deltaPhi_tagh->Write();
	h_deltaK_tagh->Write();
	h_deltaK_tagh_cut->Write();
	h_deltaK2_tagh->Write();
	h_deltaK2_tagh_cut->Write();
	
	h_deltaE_tagm->Write();
	h_deltaPhi_tagm->Write();
	h_deltaK_tagm->Write();
	h_deltaK_tagm_cut->Write();
	h_deltaK2_tagm->Write();
	h_deltaK2_tagm_cut->Write();
	
	h_fcal_xy->Write();
	h_ccal_xy->Write();
	
	h_deltaE_vs_deltaK->Write();
	h_deltaE_vs_deltaPhi->Write();
	h_deltaPhi_vs_deltaK->Write();
	h_deltaE_vs_deltaK2->Write();
	h_deltaE_vs_deltaK2_main->Write();
	h_deltaE_vs_deltaK2_acc->Write();
	
	fOut->Write();
	
	return;
}


void init_histograms() 
{
	
	h_fcal_rf_dt = new TH1F("fcal_rf_dt", "t_{FCAL} - t_{RF}; [ns]", 2000, -100., 100.);
	h_ccal_rf_dt = new TH1F("ccal_rf_dt", "t_{FCAL} - t_{RF}; [ns]", 2000, -100., 100.);
	h_beam_rf_dt = new TH1F("beam_rf_dt", "t_{FCAL} - t_{RF}; [ns]", 2000, -100., 100.);
	h_beam_rf_dt_cut = new TH1F("beam_rf_dt_cut", "t_{FCAL} - t_{RF}; [ns]", 2000, -100., 100.);
	
	h_n_cands = new TH1F("n_cands", "Number of Compton Candidates per Event", 15, -0.5, 14.5);
	
	h_deltaE_tagh = new TH2F("deltaE_tagh",   
		"#DeltaE; TAGH Counter; E_{1} + E_{2} - E_{#gamma} [GeV]", 
		274, 0.5, 274.5, 2000, -4.0, 4.0);
	h_deltaE_tagm = new TH2F("deltaE_tagm",   
		"#DeltaE; TAGM Counter; E_{1} + E_{2} - E_{#gamma} [GeV]", 
		102, 0.5, 102.5, 2000, -4.0, 4.0);
	
	h_deltaPhi_tagh = new TH2F("deltaPhi_tagh", 
		"#Delta#phi; TAGH Counter; |#phi_{1} - #phi_{2}| [deg.]", 
		274, 0.5, 274.5, 3600, 0.0, 360.0);
	h_deltaPhi_tagm = new TH2F("deltaPhi_tagm", 
		"#Delta#phi; TAGM Counter; |#phi_{1} - #phi_{2}| [deg.]", 
		102, 0.5, 102.5, 3600, 0.0, 360.0);
	
	h_deltaK_tagh = new TH2F("deltaK_tagh", 
		"#DeltaK; TAGH Counter; E_{Compton} - E_{#gamma} [GeV]", 
		274, 0.5, 274.5, 2000, -4.0, 4.0);
	h_deltaK_tagm = new TH2F("deltaK_tagm", 
		"#DeltaK; TAGM Counter; E_{Compton} - E_{#gamma} [GeV]", 
		102, 0.5, 102.5, 2000, -4.0, 4.0);
	
	h_deltaK_tagh_cut = new TH2F("deltaK_tagh_cut", 
		"#DeltaK; TAGH Counter; E_{Compton} - E_{#gamma} [GeV]", 
		274, 0.5, 274.5, 2000, -4.0, 4.0);
	h_deltaK_tagm_cut = new TH2F("deltaK_tagm_cut", 
		"#DeltaK; TAGM Counter; E_{Compton} - E_{#gamma} [GeV]", 
		102, 0.5, 102.5, 2000, -4.0, 4.0);
	
	h_deltaK2_tagh = new TH2F("deltaK2_tagh", 
		"#DeltaK; TAGH Counter; E_{Compton} - E_{#gamma} [GeV]", 
		274, 0.5, 274.5, 2000, -8.0, 8.0);
	h_deltaK2_tagm = new TH2F("deltaK2_tagm", 
		"#DeltaK; TAGM Counter; E_{Compton} - E_{#gamma} [GeV]", 
		102, 0.5, 102.5, 2000, -8.0, 8.0);
	
	h_deltaK2_tagh_cut = new TH2F("deltaK2_tagh_cut", 
		"#DeltaK; TAGH Counter; E_{Compton} - E_{#gamma} [GeV]", 
		274, 0.5, 274.5, 2000, -8.0, 8.0);
	h_deltaK2_tagm_cut = new TH2F("deltaK2_tagm_cut", 
		"#DeltaK; TAGM Counter; E_{Compton} - E_{#gamma} [GeV]", 
		102, 0.5, 102.5, 2000, -8.0, 8.0);
	
	
	h_fcal_xy = new TH2F("fcal_xy", "FCAL Y vs. X; x_{FCAL} [cm]; y_{FCAL} [cm]", 
		1000, -100., 100., 1000, -100., 100.);
	h_ccal_xy = new TH2F("ccal_xy", "CCAL Y vs. X; x_{CCAL} [cm]; y_{CCAL} [cm]", 
		1000,  -15.,  15., 1000,  -15.,  15.);
	
	h_deltaE_vs_deltaK = new TH2F("deltaE_vs_deltaK", 
	"#DeltaE vs. #DeltaK; E_{Compton} - E_{#gamma} [GeV]; (E_{1} + E_{2}) - E_{#gamma} [GeV]", 
	500, -4., 4., 500, -4., 4.);
	
	h_deltaE_vs_deltaPhi = new TH2F("deltaE_vs_deltaPhi", 
	"#DeltaE vs. #Delta#phi; |#phi_{1} - #phi_{2}| [#circ]; (E_{1} + E_{2}) - E_{#gamma} [GeV]", 
	500, 0., 360., 500, -4., 4.);
	
	h_deltaPhi_vs_deltaK = new TH2F("deltaPhi_vs_deltaK", 
	"#DeltaE vs. #DeltaK; E_{Compton} - E_{#gamma} [GeV]; |#phi_{1} - #phi_{2}| [#circ]", 
	500, -4., 4., 500, 0., 360.);
	
	h_deltaE_vs_deltaK2 = new TH2F("deltaE_vs_deltaK2", 
	"#DeltaE vs. #DeltaK; E_{Compton} - E_{#gamma} [GeV]; (E_{1} + E_{2}) - E_{#gamma} [GeV]", 
	500, -8., 8., 500, -4., 4.);
	
	h_deltaE_vs_deltaK2_main = new TH2F("deltaE_vs_deltaK2_main", 
	"#DeltaE vs. #DeltaK (Main RF Bunch); E_{Compton} - E_{#gamma} [GeV]; (E_{1} + E_{2}) - E_{#gamma} [GeV]", 
	500, -8., 8., 500, -4., 4.);
	
	h_deltaE_vs_deltaK2_acc = new TH2F("deltaE_vs_deltaK2_acc", 
	"#DeltaE vs. #DeltaK (Accidentals); E_{Compton} - E_{#gamma} [GeV]; (E_{1} + E_{2}) - E_{#gamma} [GeV]", 
	500, -8., 8., 500, -4., 4.);
	
	return;
}


void load_constants(int group, bool is_first) 
{
	/*
	9/20/21:
	  - improved distribution fits and changed cuts for Be and He-200nA runs.
	  - still need to do the same for the low intensity runs
	  - cuts for simulation were also updated in the compton_simulation and 
	  	compton_simulation_systematics plugins
	*/
	
	if(group < 3) { // Be Target
		
		deltaE_mu_p0    =  8.33517e-03;
		deltaE_mu_p1    =  2.09025e-03;
		deltaE_mu_p2    = -1.09342e-04;
		deltaE_mu_p3    =  0.;
		//---------------------------//
		deltaE_sig_p0   =  8.37004e-03;
		deltaE_sig_p1   =  4.56259e-02;
		deltaE_sig_p2   =  0.;
		
		deltaPhi_mu_p0  =  1.79943e+02;
		deltaPhi_mu_p1  = -2.11766e-02;
		deltaPhi_mu_p2  =  0.;
		deltaPhi_mu_p3  =  0.;
		//---------------------------//
		deltaPhi_sig_p0 =  1.20139e+01;
		deltaPhi_sig_p1 = -1.75486e+00;
		deltaPhi_sig_p2 =  1.67515e-01;
		deltaPhi_sig_p3 = -5.48316e-03;
		
		deltaK_mu_p0    = -1.35516e-02;
		deltaK_mu_p1    =  1.83806e-02;
		deltaK_mu_p2    = -4.03904e-03;
		deltaK_mu_p3    =  1.54466e-04;
		//---------------------------// 
		deltaK_sig_p0   = -6.63872e-02;
		deltaK_sig_p1   =  5.36909e-02;
		deltaK_sig_p2   = -3.35098e-03;
		deltaK_sig_p3   =  1.07868e-04;
		
		deltaK2_mu_p0   = -9.36095e-02;
		deltaK2_mu_p1   =  5.48923e-02;
		deltaK2_mu_p2   = -1.19844e-02;
		deltaK2_mu_p3   =  4.38188e-04;
		//---------------------------// 
		deltaK2_sig_p0  =  6.68283e-01;
		deltaK2_sig_p1  = -8.45642e-02;
		deltaK2_sig_p2  =  1.61255e-02;
		deltaK2_sig_p3  = -5.93363e-04;
		
	} else if(group < 5) { // He 200nA
		
		deltaE_mu_p0    =  4.10435e-01;
		deltaE_mu_p1    = -4.20598e-02;
		deltaE_mu_p2    = -6.19015e-04;
		deltaE_mu_p3    =  0.;
		//---------------------------//
		deltaE_sig_p0   =  1.26554e-02;
		deltaE_sig_p1   =  4.28404e-02;
		deltaE_sig_p2   =  0.;
		
		deltaPhi_mu_p0  =  1.79777e+02;
		deltaPhi_mu_p1  = -7.73534e-03;
		deltaPhi_mu_p2  =  0.;
		deltaPhi_mu_p3  =  0.;
		//---------------------------// 
		deltaPhi_sig_p0 =  1.11023e+01;
		deltaPhi_sig_p1 = -1.58498e+00;
		deltaPhi_sig_p2 =  1.53958e-01;
		deltaPhi_sig_p3 = -5.09573e-03;
		
		deltaK_mu_p0    = -4.23755e-02;
		deltaK_mu_p1    =  5.52975e-02;
		deltaK_mu_p2    = -9.11057e-03;
		deltaK_mu_p3    =  3.35119e-04;
		//---------------------------// 
		deltaK_sig_p0   = -3.52834e-02;
		deltaK_sig_p1   =  4.58032e-02;
		deltaK_sig_p2   = -2.87025e-03;
		deltaK_sig_p3   =  1.05750e-04;
		
		deltaK2_mu_p0   = -1.74619e-01;
		deltaK2_mu_p1   =  1.01328e-01;
		deltaK2_mu_p2   = -1.84593e-02;
		deltaK2_mu_p3   =  7.13390e-04;
		//---------------------------// 
		deltaK2_sig_p0  =  3.88777e-01;
		deltaK2_sig_p1  =  2.50588e-02;
		deltaK2_sig_p2  =  2.20573e-03;
		deltaK2_sig_p3  = -2.76677e-05;
		
	} else if(group < 6) { // He 100nA
		
		deltaE_mu_p0    =  4.10435e-01;
		deltaE_mu_p1    = -4.20598e-02;
		deltaE_mu_p2    = -6.19015e-04;
		deltaE_mu_p3    =  0.;
		//---------------------------//
		deltaE_sig_p0   =  1.26554e-02;
		deltaE_sig_p1   =  4.28404e-02;
		deltaE_sig_p2   =  0.;
		
		deltaPhi_mu_p0  =  1.79777e+02;
		deltaPhi_mu_p1  = -7.73534e-03;
		deltaPhi_mu_p2  =  0.;
		deltaPhi_mu_p3  =  0.;
		//---------------------------// 
		deltaPhi_sig_p0 =  1.11023e+01;
		deltaPhi_sig_p1 = -1.58498e+00;
		deltaPhi_sig_p2 =  1.53958e-01;
		deltaPhi_sig_p3 = -5.09573e-03;
		
		deltaK_mu_p0    = -4.23755e-02;
		deltaK_mu_p1    =  5.52975e-02;
		deltaK_mu_p2    = -9.11057e-03;
		deltaK_mu_p3    =  3.35119e-04;
		//---------------------------// 
		deltaK_sig_p0   = -3.52834e-02;
		deltaK_sig_p1   =  4.58032e-02;
		deltaK_sig_p2   = -2.87025e-03;
		deltaK_sig_p3   =  1.05750e-04;
		
		deltaK2_mu_p0   = -1.74619e-01;
		deltaK2_mu_p1   =  1.01328e-01;
		deltaK2_mu_p2   = -1.84593e-02;
		deltaK2_mu_p3   =  7.13390e-04;
		//---------------------------// 
		deltaK2_sig_p0  =  3.88777e-01;
		deltaK2_sig_p1  =  2.50588e-02;
		deltaK2_sig_p2  =  2.20573e-03;
		deltaK2_sig_p3  = -2.76677e-05;
		
	} else if(group < 7) { // He 050nA
		
		deltaE_mu_p0    =  4.10435e-01;
		deltaE_mu_p1    = -4.20598e-02;
		deltaE_mu_p2    = -6.19015e-04;
		deltaE_mu_p3    =  0.;
		//---------------------------//
		deltaE_sig_p0   =  1.26554e-02;
		deltaE_sig_p1   =  4.28404e-02;
		deltaE_sig_p2   =  0.;
		
		deltaPhi_mu_p0  =  1.79777e+02;
		deltaPhi_mu_p1  = -7.73534e-03;
		deltaPhi_mu_p2  =  0.;
		deltaPhi_mu_p3  =  0.;
		//---------------------------// 
		deltaPhi_sig_p0 =  1.11023e+01;
		deltaPhi_sig_p1 = -1.58498e+00;
		deltaPhi_sig_p2 =  1.53958e-01;
		deltaPhi_sig_p3 = -5.09573e-03;
		
		deltaK_mu_p0    = -4.23755e-02;
		deltaK_mu_p1    =  5.52975e-02;
		deltaK_mu_p2    = -9.11057e-03;
		deltaK_mu_p3    =  3.35119e-04;
		//---------------------------// 
		deltaK_sig_p0   = -3.52834e-02;
		deltaK_sig_p1   =  4.58032e-02;
		deltaK_sig_p2   = -2.87025e-03;
		deltaK_sig_p3   =  1.05750e-04;
		
		deltaK2_mu_p0   = -1.74619e-01;
		deltaK2_mu_p1   =  1.01328e-01;
		deltaK2_mu_p2   = -1.84593e-02;
		deltaK2_mu_p3   =  7.13390e-04;
		//---------------------------// 
		deltaK2_sig_p0  =  3.88777e-01;
		deltaK2_sig_p1  =  2.50588e-02;
		deltaK2_sig_p2  =  2.20573e-03;
		deltaK2_sig_p3  = -2.76677e-05;
		
	} else {
		
		cout << "\n\n\nInvalid run group specified.\n\n" << endl;
		exit(-1);
	}
	
	if(is_first) {
		f_deltaE_mu    = new TF1("f_deltaE_mu",    "pol3", 5.0, 12.0);
		f_deltaE_sig   = new TF1("f_deltaE_sig",   
			"sqrt([0]*[0] + ([1]/sqrt(x))*([1]/sqrt(x)) + ([2]/x)*([2]/x))", 5.0, 12.0);
		f_deltaPhi_mu  = new TF1("f_deltaPhi_mu",  "pol3", 5.0, 12.0);
		f_deltaPhi_sig = new TF1("f_deltaPhi_sig", "pol3", 5.0, 12.0);
		f_deltaK_mu    = new TF1("f_deltaK_mu",    "pol3", 5.0, 12.0);
		f_deltaK_sig   = new TF1("f_deltaK_sig",   "pol3", 5.0, 12.0);
		f_deltaK2_mu   = new TF1("f_deltaK2_mu",   "pol3", 5.0, 12.0);
		f_deltaK2_sig  = new TF1("f_deltaK2_sig",  "pol3", 5.0, 12.0);
	}
	
	f_deltaE_mu->SetParameters(deltaE_mu_p0, deltaE_mu_p1, deltaE_mu_p2, deltaE_mu_p3);
	f_deltaE_sig->SetParameters(deltaE_sig_p0, deltaE_sig_p1, deltaE_sig_p2);
	
	f_deltaPhi_mu->SetParameters(deltaPhi_mu_p0,  deltaPhi_mu_p1,  
		deltaPhi_mu_p2,  deltaPhi_mu_p3);
	f_deltaPhi_sig->SetParameters(deltaPhi_sig_p0, deltaPhi_sig_p1, 
		deltaPhi_sig_p2, deltaPhi_sig_p3);
	
	f_deltaK_mu->SetParameters(deltaK_mu_p0,  deltaK_mu_p1,  deltaK_mu_p2,  deltaK_mu_p3);
	f_deltaK_sig->SetParameters(deltaK_sig_p0, deltaK_sig_p1, deltaK_sig_p2, deltaK_sig_p3);
	
	f_deltaK2_mu->SetParameters(deltaK2_mu_p0,  deltaK2_mu_p1,  deltaK2_mu_p2,  deltaK2_mu_p3);
	f_deltaK2_sig->SetParameters(deltaK2_sig_p0, deltaK2_sig_p1, deltaK2_sig_p2, deltaK2_sig_p3);
	
	return;
}
