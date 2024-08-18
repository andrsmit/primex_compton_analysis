
#include "compton.h"

void set_cuts() {
	
	FCAL_ENERGY_CUT  = 0.35;
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
	"/work/halld/home/andrsmit/primex_compton_analysis/data/rootTrees/phaseIII");
	
	// Directory where output ROOT files will be stored:
	
	sprintf(rootFile_pathName, 
	"/work/halld/home/andrsmit/primex_compton_analysis/analyze_trees/phaseIII/analysis/rootFiles");
	
	//--------------------------------------------------------------------------------//
	
	vector<int> Be_empty_FIELDON_runs  = {
	110447, 110448, 110449, 110450, 110451, 110453, 110454, 110455};
	vector<int> Be_200nA_FIELDON_runs  = {
	110497, 110498, 110499, 110500, 110501, 110503, 110504, 110505, 110515, 110516, 110517, 
	110532, 110533, 110534, 110535, 110536, 110537, 110538, 110539, 110540, 110541, 110544, 
	110545, 110546, 110547, 110548, 110549, 110551, 110552, 110553};
	vector<int> Be_100nA_FIELDON_runs  = {
	110554, 110555, 110556, 110557, 110558, 110559, 110560, 110561, 110562, 110563, 110564};
	vector<int> Be_050nA_FIELDON_runs  = {
	110565, 110567, 110576, 110577, 110578, 110579, 110580, 110581, 110582, 110583};
	
	vector<int> Be_empty_FIELDOFF_runs = {
	110476, 110477, 110478, 110479, 110481, 110482};
	vector<int> Be_200nA_FIELDOFF_runs = {
	110600, 110601, 110602, 110603, 110605, 110606, 110607, 110608, 110609, 110610, 110615};
	vector<int> Be_100nA_FIELDOFF_runs = {
	110616, 110617, 110618, 110619, 110620, 110621};
	
	//--------------------------------------------------------------------------------//
	
	vector<int> He_empty_FIELDON_runs  = {
	110737, 110738, 110739, 110740, 110741, 110743, 110744, 110745, 110746, 110747, 110749, 
	110750, 110751, 110752, 110753, 110754, 110755, 
	110828, 110829, 110830, 110831, 110833, 110834, 110835, 110836, 110837, 
	110907, 110909, 110910, 110911, 110912, 110913, 110914, 110915, 110916, 110917, 110918, 
	110919, 110920, 110921, 110922, 110923, 110924, 110925, 110926, 110927, 
	111013, 111014, 111015, 111016, 111017, 111018, 111019, 111020, 111021, 111022, 111023, 
	111024, 111025, 111026, 111027, 111028, 111029, 111030, 111031, 111032, 111033, 111034, 
	111035, 111036, 111037, 
	111096, 111097, 111098, 111099, 111100, 111101, 111102, 111103, 111106, 111107, 111108, 
	111112, 111113, 111114, 111115, 111116, 111117, 111118, 111119, 
	111180, 111182, 111183, 111184, 111185, 111186, 111187, 111188, 111189, 111190, 111191, 
	111192, 111193, 111194, 111195, 111196, 111197, 111201, 111202, 111203, 111204, 111205, 
	111206, 111207, 111208, 111209, 111210, 111211, 111212, 111215
	};
	vector<int> He_200nA_FIELDON_runs  = {
	110641, 110642, 110643, 110644, 110657, 110672, 110675, 110676, 110677, 110678, 110679, 
	110680, 110682, 110683, 110684, 110685, 110686, 110687, 110688, 110689, 110690, 110691, 
	110692, 110693, 110694, 110695, 110696, 110697, 110698, 110699, 110700, 110702, 110703, 
	110704, 110705, 110706, 110707, 110708, 110709, 110710, 110711, 110713, 110714, 110715, 
	110716, 110718, 110719, 110724, 110725, 
	110771, 110772, 110773, 110774, 110775, 110776, 
	110782, 110783, 110784, 110785, 110786, 110787, 110788, 110789, 110790, 110791, 110792, 
	110793, 110794, 110795, 110796, 110798, 110799, 110800, 110801, 110802, 110803, 110804, 
	110805, 110806, 110807, 110808, 110809, 110810, 110811, 110812, 110814, 110815, 110816, 
	110817, 110818, 110819, 110820, 110821, 110822, 110824, 110825, 
	110841, 110842, 110843, 110844, 110845, 110846, 110847, 110848, 110849, 110850, 110851, 
	110852, 110853, 110854, 110855, 110856, 110857, 110858, 110859, 110860, 110861, 110862, 
	110863, 110864, 110865, 110866, 110867, 110868, 110869, 110870, 110871, 110872, 110873, 
	110874, 110875, 110876, 110877, 110878, 110879, 110881, 110882, 110896, 110897, 110899, 
	110900, 110901, 110902, 110903, 110904, 110905, 
	110929, 110930, 110931, 110932, 110933, 110934, 110935, 110936, 110937, 110942, 110943, 
	110944, 110945, 110946, 110949, 110950, 110951, 110953, 110954, 110955, 110956, 110957, 
	110958, 110959, 110960, 110961, 110962, 110963, 110964, 110965, 110966, 110967, 110968, 
	110969, 110970, 110971, 110972, 110973, 
	110987, 110988, 110989, 110990, 110991, 
	110993, 110994, 110995, 110996, 110997, 110998, 110999, 111000, 111001, 111002, 111003, 
	111004, 111005, 111006, 111007, 111008, 111009, 111010, 
	111041, 111042, 111043, 111044, 111045, 111046, 111047, 111048, 111050, 111052, 111053, 
	111054, 111055, 111056, 111057, 
	111063, 111064, 111065, 111067, 111068, 111069, 111070, 111071, 111072, 111073, 111074, 
	111075, 111077, 111078, 111079, 111080, 111081, 111082, 111083, 111084, 
	111085, 111087, 111088, 111089, 111090, 111092, 111093, 111094, 111095, 
	111120, 111121, 111122, 111123, 111124, 111126, 111127, 111128, 111129, 111131, 111132, 
	111133, 111134, 111135, 111136, 
	111137, 111140, 111141, 111142, 111143, 111144, 111145, 
	111146, 111147, 111148, 111149, 111150, 111151, 111152, 111153, 111154, 111155, 111156, 
	111157, 111158, 111159, 111160, 111161, 111162, 111163, 111164, 111165, 111166, 111167, 
	111168, 111171, 111172, 111173, 111174, 111175, 111176, 111177, 111178, 111179
	};
	vector<int> He_100nA_FIELDON_runs  = {
	};
	vector<int> He_050nA_FIELDON_runs  = {
	};
	vector<int> He_250nA_FIELDON_runs  = {
	111061, 111062};
	
	vector<int> He_empty_FIELDOFF_runs = {
	};
	vector<int> He_200nA_FIELDOFF_runs = {
	111969, 111970, 111971, 111972, 111973, 111974, 111975, 111976, 111977, 111978, 111979, 
	111984, 111985, 111986, 111987, 111988, 111989, 111990, 111991, 111992, 111993, 111994, 
	111995, 111996, 111997, 111998, 111999, 112000, 112001};
	vector<int> He_100nA_FIELDOFF_runs = {
	};
	vector<int> He_050nA_FIELDOFF_runs = {
	};
	vector<int> He_250nA_FIELDOFF_runs = {
	};
	
	//--------------------------------------------------------------------------------//
	
	init_histograms();
	
	bool first_evt = true;
	
	for(unsigned int irun = startRun; irun <= endRun; irun++) {
		
		reset_histograms();
		
		//-----   Check which group of runs this run belongs to   -----//
		
		int run_group = 0;
		for(int jr = 0; jr < (int)Be_200nA_FIELDOFF_runs.size(); jr++) {
			if(irun == Be_200nA_FIELDOFF_runs[jr]) {
				run_group = 1;
				break;
			}
		}
		if(!run_group) {
			for(int jr = 0; jr < (int)Be_100nA_FIELDOFF_runs.size(); jr++) {
				if(irun == Be_100nA_FIELDOFF_runs[jr]) {
					run_group = 2;
					break;
				}
			}
		}
		if(!run_group) {
			for(int jr = 0; jr < (int)Be_empty_FIELDOFF_runs.size(); jr++) {
				if(irun == Be_empty_FIELDOFF_runs[jr]) {
					run_group = 3;
					break;
				}
			}
		}
		if(!run_group) {
			for(int jr = 0; jr < (int)Be_200nA_FIELDON_runs.size(); jr++) {
				if(irun == Be_200nA_FIELDON_runs[jr]) {
					run_group = 4;
					break;
				}
			}
		}
		if(!run_group) {
			for(int jr = 0; jr < (int)Be_100nA_FIELDON_runs.size(); jr++) {
				if(irun == Be_100nA_FIELDON_runs[jr]) {
					run_group = 5;
					break;
				}
			}
		}
		if(!run_group) {
			for(int jr = 0; jr < (int)Be_050nA_FIELDON_runs.size(); jr++) {
				if(irun == Be_050nA_FIELDON_runs[jr]) {
					run_group = 6;
					break;
				}
			}
		}
		if(!run_group) {
			for(int jr = 0; jr < (int)Be_empty_FIELDON_runs.size(); jr++) {
				if(irun == Be_empty_FIELDON_runs[jr]) {
					run_group = 7;
					break;
				}
			}
		}
		if(!run_group) {
			for(int jr = 0; jr < (int)He_200nA_FIELDOFF_runs.size(); jr++) {
				if(irun == He_200nA_FIELDOFF_runs[jr]) {
					run_group = 8;
					break;
				}
			}
		}
		if(!run_group) {
			for(int jr = 0; jr < (int)He_200nA_FIELDON_runs.size(); jr++) {
				if(irun == He_200nA_FIELDON_runs[jr]) {
					run_group = 9;
					break;
				}
			}
		}
		if(!run_group) {
			for(int jr = 0; jr < (int)He_250nA_FIELDON_runs.size(); jr++) {
				if(irun == He_250nA_FIELDON_runs[jr]) {
					run_group = 10;
					break;
				}
			}
		}
		if(!run_group) {
			for(int jr = 0; jr < (int)He_empty_FIELDON_runs.size(); jr++) {
				if(irun == He_empty_FIELDON_runs[jr]) {
					run_group = 11;
					break;
				}
			}
		}
		if(!run_group) continue;
		
		//m_beamX = 0.146;
		//m_beamY = 0.017;
		m_beamX = 0.151;
		m_beamY = 0.012;
		if(run_group < 8) m_beamZ = 64.935;
		else 							m_beamZ = 65.0;
		
		cout << "Processing run " << irun << endl;
		
		load_constants(run_group, first_evt);
		first_evt = false;
		
		for(int iext = 0; iext < 200; iext++ ) {
			
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
		
		/*
		Update FCAL position:
		old beam position: (0.146, 0.017)
		new beam position: (0.151, 0.012)
		
		old CCAL position: (0.135, 0.135)
		new CCAL position: (0.184, 0.110)
		*/
		
		//x1 = x1 - (0.151-0.146);
		//y1 = y1 - (0.012-0.017);
		
		double fcal_face_x = m_beamX + (x1 * (m_fcalZ - m_beamZ)/z1) - m_fcalX;
		double fcal_face_y = m_beamY + (y1 * (m_fcalZ - m_beamZ)/z1) - m_fcalY;
		
		int fcal_fid_cut = 0;
		if((-1.*fcal_inner_layer_cut < fcal_face_x && fcal_face_x < fcal_inner_layer_cut)
			&& (-1.*fcal_inner_layer_cut < fcal_face_y 
			&& fcal_face_y < fcal_inner_layer_cut)) fcal_fid_cut = 1;
		
		double x2 = ccal_x[ic];
		double y2 = ccal_y[ic];
		double z2 = ccal_z[ic];
		
		//x2 = x2 - (0.151-0.146) + (0.184-0.135);
		//y2 = y2 - (0.012-0.017) + (0.110-0.135);
		
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
		if(fcal_e[ic] > FCAL_ENERGY_CUT) {
			h_fcal_rf_dt->Fill(fcal_t[ic]-rfTime);
			fcal_e_cut = 1;
		}
		if(ccal_e[ic] > CCAL_ENERGY_CUT) {
			h_ccal_rf_dt->Fill(ccal_t[ic]-rfTime);
			ccal_e_cut = 1;
		}
		
		//----------      Timing Cuts      -----------//
		
		int fcal_t_cut = 0, ccal_t_cut = 0;
		if(fabs(fcal_t[ic] - rfTime) < FCAL_RF_TIME_CUT) fcal_t_cut = 1;
		if(fabs(ccal_t[ic] - rfTime) < CCAL_RF_TIME_CUT) ccal_t_cut = 1;
		
		//--------------------------------------------//
		
		if(fcal_fid_cut || ccal_fid_cut) continue;
		if(!fcal_e_cut || !ccal_e_cut)   continue;
		if(!fcal_t_cut || !ccal_t_cut)   continue;
		
		//--------------------------------------------//
		
		double fill_weight = 0.;
		if(bunch_val[ic]) {
			if(fabs(tb[ic]-rfTime) < 2.004) fill_weight = 1.0;
		} else {
			fill_weight = -0.1;
		}
		
		h_beam_rf_dt->Fill(tb[ic]-rfTime);
		
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
		
		
		double deltaKE = (fcal_e[ic]+ccal_e[ic])-ecomp;
		
		// rotated deltaE vs. deltaK:
		
		double deltaH = loc_deltaK2*cos(3.14159/4.0) - deltaE[ic]*sin(3.14159/4.0);
		
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
		
		// plot tagged photon distribution for background events:
		
		double x_cord = loc_deltaK2;
		double y_cord = deltaE[ic];
		
		// y=mx+b where m=1 and b=0.8 or -0.8:
		
		double cut_lo = x_cord - 0.8;
		double cut_hi = x_cord + 0.8;
		if(cut_lo<y_cord && y_cord<cut_hi && p_cut) {
			if(fabs(deltaE[ic])>1.0) {
				if(tag_sys[ic]==0) {
					h_tagh_background->Fill(tag_counter[ic]);
					if(deltaE[ic]<0.) h_tagh_background_lo->Fill(tag_counter[ic]);
					else              h_tagh_background_hi->Fill(tag_counter[ic]);
				} else {
					h_tagm_background->Fill(tag_counter[ic]);
					if(deltaE[ic]<0.) h_tagm_background_lo->Fill(tag_counter[ic]);
					else              h_tagm_background_hi->Fill(tag_counter[ic]);
				}
			}
			h_deltaE_vs_esum->Fill(fcal_e[ic]+ccal_e[ic], deltaE[ic], fill_weight);
			if((fcal_e[ic]+ccal_e[ic])>9.) {
				if(fabs(deltaE[ic])>0.5 && fabs(deltaE[ic])<1.5) {
					if(deltaE[ic]<0.) h_tagged_bkgd->Fill(0.,fill_weight);
					else              h_tagged_bkgd->Fill(1.,fill_weight);
				}
				if(fabs(deltaE[ic])>1.0 && fabs(deltaE[ic])<2.0) {
					if(deltaE[ic]<0.) h_tagged_bkgd_1GeV_cut->Fill(0.,fill_weight);
					else              h_tagged_bkgd_1GeV_cut->Fill(1.,fill_weight);
				}
			}
		}
		
		
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
			
			if(p_cut) {
				h_deltaH_tagh[0]->Fill(tag_counter[ic], deltaH, fill_weight);
				h_deltaKE_tagh[0]->Fill(tag_counter[ic], deltaKE, fill_weight);
				h_esum_vs_ecomp[0]->Fill(ecomp, fcal_e[ic]+ccal_e[ic], fill_weight);
				if(deltaE[ic]<-1.0) {
					h_deltaH_tagh[1]->Fill(tag_counter[ic], deltaH, fill_weight);
					h_deltaKE_tagh[1]->Fill(tag_counter[ic], deltaKE, fill_weight);
					h_esum_vs_ecomp[1]->Fill(ecomp, fcal_e[ic]+ccal_e[ic], fill_weight);
				}
				if(e_cut) {
					h_deltaH_tagh[2]->Fill(tag_counter[ic], deltaH, fill_weight);
					h_deltaKE_tagh[2]->Fill(tag_counter[ic], deltaKE, fill_weight);
					h_esum_vs_ecomp[2]->Fill(ecomp, fcal_e[ic]+ccal_e[ic], fill_weight);
					if(k2_cut) {
						h_deltaH_tagh[3]->Fill(tag_counter[ic], deltaH, fill_weight);
						h_deltaKE_tagh[3]->Fill(tag_counter[ic], deltaKE, fill_weight);
						h_esum_vs_ecomp[3]->Fill(ecomp, fcal_e[ic]+ccal_e[ic], fill_weight);
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
			
			if(p_cut) {
				h_deltaH_tagm[0]->Fill(tag_counter[ic], deltaH, fill_weight);
				h_deltaKE_tagm[0]->Fill(tag_counter[ic], deltaKE, fill_weight);
				h_esum_vs_ecomp[0]->Fill(ecomp, fcal_e[ic]+ccal_e[ic], fill_weight);
				if(deltaE[ic]<-1.0) {
					h_deltaH_tagm[1]->Fill(tag_counter[ic], deltaH, fill_weight);
					h_deltaKE_tagm[1]->Fill(tag_counter[ic], deltaKE, fill_weight);
					h_esum_vs_ecomp[1]->Fill(ecomp, fcal_e[ic]+ccal_e[ic], fill_weight);
				}
				if(e_cut) {
					h_deltaH_tagm[2]->Fill(tag_counter[ic], deltaH, fill_weight);
					h_deltaKE_tagm[2]->Fill(tag_counter[ic], deltaKE, fill_weight);
					h_esum_vs_ecomp[2]->Fill(ecomp, fcal_e[ic]+ccal_e[ic], fill_weight);
					if(k2_cut) {
						h_deltaH_tagm[3]->Fill(tag_counter[ic], deltaH, fill_weight);
						h_deltaKE_tagm[3]->Fill(tag_counter[ic], deltaKE, fill_weight);
						h_esum_vs_ecomp[3]->Fill(ecomp, fcal_e[ic]+ccal_e[ic], fill_weight);
					}
				}
			}
		}
		
		h_deltaE_vs_deltaK->Fill(loc_deltaK, deltaE[ic], fill_weight);
		h_deltaE_vs_deltaPhi->Fill(deltaPhi[ic], deltaE[ic], fill_weight);
		h_deltaPhi_vs_deltaK->Fill(loc_deltaK, deltaPhi[ic], fill_weight);
		if((fcal_e[ic]+ccal_e[ic])>8.0) h_deltaE_vs_deltaK2->Fill(loc_deltaK2, deltaE[ic], fill_weight);
		
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
	
	for(int i=0; i<4; i++) {
		h_deltaH_tagh[i]->Reset();
		h_deltaH_tagm[i]->Reset();
		h_deltaKE_tagh[i]->Reset();
		h_deltaKE_tagm[i]->Reset();
		h_esum_vs_ecomp[i]->Reset();
	}
	
	h_fcal_xy->Reset();
	h_ccal_xy->Reset();
	
	h_deltaE_vs_deltaK->Reset();
	h_deltaE_vs_deltaPhi->Reset();
	h_deltaPhi_vs_deltaK->Reset();
	h_deltaE_vs_deltaK2->Reset();
	h_deltaE_vs_deltaK2_main->Reset();
	h_deltaE_vs_deltaK2_acc->Reset();
	
	h_tagh_background->Reset();
	h_tagh_background_lo->Reset();
	h_tagh_background_hi->Reset();
	h_tagm_background->Reset();
	h_tagm_background_lo->Reset();
	h_tagm_background_hi->Reset();
	
	h_tagged_bkgd->Reset();
	h_tagged_bkgd_1GeV_cut->Reset();
	h_deltaE_vs_esum->Reset();
	
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
	
	for(int i=0; i<4; i++) {
		h_deltaH_tagh[i]->Write();
		h_deltaH_tagm[i]->Write();
		h_deltaKE_tagh[i]->Write();
		h_deltaKE_tagm[i]->Write();
		h_esum_vs_ecomp[i]->Write();
	}
	
	h_fcal_xy->Write();
	h_ccal_xy->Write();
	
	h_deltaE_vs_deltaK->Write();
	h_deltaE_vs_deltaPhi->Write();
	h_deltaPhi_vs_deltaK->Write();
	h_deltaE_vs_deltaK2->Write();
	h_deltaE_vs_deltaK2_main->Write();
	h_deltaE_vs_deltaK2_acc->Write();
	
	h_tagh_background->Write();
	h_tagh_background_lo->Write();
	h_tagh_background_hi->Write();
	h_tagm_background->Write();
	h_tagm_background_lo->Write();
	h_tagm_background_hi->Write();
	
	h_tagged_bkgd->Write();
	h_tagged_bkgd_1GeV_cut->Write();
	h_deltaE_vs_esum->Write();
	
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
	
	
	for(int ih=0; ih<4; ih++) {
		h_deltaH_tagh[ih] = new TH2F(Form("deltaH_tagh_%d",ih), 
			"#DeltaK cos(#alpha) - #DeltaE sin(#alpha)", 
			274, 0.5, 274.5, 2000, -8.0, 8.0);
		h_deltaH_tagm[ih] = new TH2F(Form("deltaH_tagm_%d",ih), 
			"#DeltaK cos(#alpha) - #DeltaE sin(#alpha)", 
			102, 0.5, 102.5, 2000, -8.0, 8.0);
		
		h_deltaKE_tagh[ih] = new TH2F(Form("deltaKE_tagh_%d",ih), 
			"(E_{1}+E_{2}) - E_{Compton}", 
			274, 0.5, 274.5, 2000, -8.0, 8.0);
		h_deltaKE_tagm[ih] = new TH2F(Form("deltaKE_tagm_%d",ih), 
			"(E_{1}+E_{2}) - E_{Compton}", 
			102, 0.5, 102.5, 2000, -8.0, 8.0);
		
		h_esum_vs_ecomp[ih] = new TH2F(Form("esum_vs_ecomp_%d",ih), 
			"E_{1}+E_{2} vs. E_{Compton}", 
			2000, 0., 12., 2000, 0., 12.);
	}
	
	
	
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
	
	h_tagh_background = new TH1F("tagh_background", 
		"Occupancy of TAGH Counters for Compton Background", 274, 0.5, 274.5);
	h_tagh_background_lo = new TH1F("tagh_background_lo", 
		"Occupancy of TAGH Counters for Compton Background (#DeltaE < 0)", 274, 0.5, 274.5);
	h_tagh_background_hi = new TH1F("tagh_background_hi", 
		"Occupancy of TAGH Counters for Compton Background (#DeltaE > 0)", 274, 0.5, 274.5);
	
	h_tagm_background = new TH1F("tagm_background", 
		"Occupancy of TAGM Counters for Compton Background", 102, 0.5, 102.5);
	h_tagm_background_lo = new TH1F("tagm_background_lo", 
		"Occupancy of TAGM Counters for Compton Background (#DeltaE < 0)", 102, 0.5, 102.5);
	h_tagm_background_hi = new TH1F("tagm_background_hi", 
		"Occupancy of TAGM Counters for Compton Background (#DeltaE > 0)", 102, 0.5, 102.5);
	
	h_tagged_bkgd          = new TH1F("tagged_bkgd", 
		"Is the tagged photon energy too low?", 2, -0.5, 1.5);
	h_tagged_bkgd_1GeV_cut = new TH1F("tagged_bkgd_1GeV_cut", 
		"Is the tagged photon energy too low?", 2, -0.5, 1.5);
	
	h_deltaE_vs_esum = new TH2F("deltaE_vs_esum", 
		"#DeltaE vs. E_{1}+E_{2}; E_{1}+E_{2} [GeV]; #DeltaE [GeV]", 
		500, 0., 12., 2000, -4., 4.);
	
	return;
}


void load_constants(int group, bool is_first) 
{
	if(group < 4) { // Be FIELD OFF
		
		deltaE_mu_p0    = -1.16581e-01;
		deltaE_mu_p1    =  3.40885e-02;
		deltaE_mu_p2    = -2.77050e-03;
		deltaE_mu_p3    =  2.09654e-05;
		//---------------------------//
		deltaE_sig_p0   =  1.65687e-02;
		deltaE_sig_p1   =  2.65247e-02;
		deltaE_sig_p2   =  4.38584e-02;
		/*
		deltaE_mu_p0    = -2.85756e-01;
		deltaE_mu_p1    =  8.91005e-02;
		deltaE_mu_p2    = -8.73318e-03;
		deltaE_mu_p3    =  2.45143e-04;
		//---------------------------//
		deltaE_sig_p0   =  1.81080e-02;
		deltaE_sig_p1   =  1.38460e-02;
		deltaE_sig_p2   =  5.76009e-02;
		*/
		
		deltaPhi_mu_p0  =  1.79825e+02;
		deltaPhi_mu_p1  = -1.10610e-02;
		deltaPhi_mu_p2  =  0.;
		deltaPhi_mu_p3  =  0.;
		//---------------------------// 
		deltaPhi_sig_p0 =  1.33099e+01;
		deltaPhi_sig_p1 = -2.14366e+00;
		deltaPhi_sig_p2 =  2.06214e-01;
		deltaPhi_sig_p3 = -6.71967e-03;
		
		deltaK_mu_p0    =  9.54559e-03;
		deltaK_mu_p1    =  6.67955e-03;
		deltaK_mu_p2    = -2.25633e-03;
		deltaK_mu_p3    =  7.71613e-05;
		//---------------------------// 
		deltaK_sig_p0   = -1.17831e-01;
		deltaK_sig_p1   =  6.81072e-02;
		deltaK_sig_p2   = -4.77039e-03;
		deltaK_sig_p3   =  1.57138e-04;
		
		deltaK2_mu_p0   =  2.09726e-01;
		deltaK2_mu_p1   = -4.71020e-02;
		deltaK2_mu_p2   =  3.81198e-04;
		deltaK2_mu_p3   = -4.82342e-05;
		//---------------------------// 
		deltaK2_sig_p0  =  5.19636e-01;
		deltaK2_sig_p1  = -3.35925e-02;
		deltaK2_sig_p2  =  1.03144e-02;
		deltaK2_sig_p3  = -3.63616e-04;
		
	} 
	/*else if(group < 4) { // Be-Empty FIELD OFF
		
		deltaE_mu_p0    = -6.60845e-01;
		deltaE_mu_p1    =  2.56921e-01;
		deltaE_mu_p2    = -2.55507e-02;
		deltaE_mu_p3    =  8.94288e-04;
		//---------------------------//
		deltaE_sig_p0   =  1.67980e-02;
		deltaE_sig_p1   =  3.37356e-02;
		deltaE_sig_p2   =  8.58059e-02;
		
		deltaPhi_mu_p0  =  1.80011e+02;
		deltaPhi_mu_p1  = -8.18407e-03;
		deltaPhi_mu_p2  =  0.;
		deltaPhi_mu_p3  =  0.;
		//---------------------------// 
		deltaPhi_sig_p0 =  1.22089e+01;
		deltaPhi_sig_p1 = -1.81467e+00;
		deltaPhi_sig_p2 =  1.72870e-01;
		deltaPhi_sig_p3 = -5.55077e-03;
		
		deltaK_mu_p0    =  9.54559e-03;
		deltaK_mu_p1    =  6.67955e-03;
		deltaK_mu_p2    = -2.25633e-03;
		deltaK_mu_p3    =  7.71613e-05;
		//---------------------------// 
		deltaK_sig_p0   = -1.17831e-01;
		deltaK_sig_p1   =  6.81072e-02;
		deltaK_sig_p2   = -4.77039e-03;
		deltaK_sig_p3   =  1.57138e-04;
		
		deltaK2_mu_p0   =  3.22415e-01;
		deltaK2_mu_p1   = -8.24914e-02;
		deltaK2_mu_p2   =  4.25586e-03;
		deltaK2_mu_p3   = -1.88916e-04;
		//---------------------------// 
		deltaK2_sig_p0  =  5.70095e-01;
		deltaK2_sig_p1  = -4.93059e-02;
		deltaK2_sig_p2  =  1.19542e-02;
		deltaK2_sig_p3  = -4.17058e-04;
		
	} */
	else if(group < 8) { // Be FIELD ON
		
		deltaE_mu_p0    = -0.05;
		deltaE_mu_p1    =  0.;
		deltaE_mu_p2    =  0.;
		deltaE_mu_p3    =  0.;
		//---------------------------//
		deltaE_sig_p0   =  1.40025e-02;
		deltaE_sig_p1   =  3.55640e-02;
		deltaE_sig_p2   =  5.15193e-04;
		
		deltaPhi_mu_p0  =  1.80000e+02;
		deltaPhi_mu_p1  =  0.;
		deltaPhi_mu_p2  =  0.;
		deltaPhi_mu_p3  =  0.;
		//---------------------------// 
		deltaPhi_sig_p0 =  1.50000e+01;
		deltaPhi_sig_p1 =  0.0;
		deltaPhi_sig_p2 =  0.0;
		deltaPhi_sig_p3 =  0.0;
		
		deltaK_mu_p0    =  6.89191e-04;
		deltaK_mu_p1    =  1.65777e-02;
		deltaK_mu_p2    = -4.57181e-03;
		deltaK_mu_p3    =  2.15276e-04;
		//---------------------------// 
		deltaK_sig_p0   = -4.12018e-02;
		deltaK_sig_p1   =  4.22817e-02;
		deltaK_sig_p2   = -2.05753e-03;
		deltaK_sig_p3   =  6.12380e-05;
		
		deltaK2_mu_p0   =  6.89191e-04;
		deltaK2_mu_p1   =  1.65777e-02;
		deltaK2_mu_p2   = -4.57181e-03;
		deltaK2_mu_p3   =  2.15276e-04;
		//---------------------------// 
		deltaK2_sig_p0  = -4.12018e-02;
		deltaK2_sig_p1  =  4.22817e-02;
		deltaK2_sig_p2  = -2.05753e-03;
		deltaK2_sig_p3  =  6.12380e-05;
		
	} else if(group < 9) { // He FIELD OFF
		
		deltaE_mu_p0    = -2.22458e-01;
		deltaE_mu_p1    =  8.24677e-02;
		deltaE_mu_p2    = -8.61126e-03;
		deltaE_mu_p3    =  2.32165e-04;
		//---------------------------// 
		deltaE_sig_p0   =  1.39151e-02;
		deltaE_sig_p1   =  1.71639e-02;
		deltaE_sig_p2   =  7.21903e-02;
		
		deltaPhi_mu_p0  =  1.80814e+02;
		deltaPhi_mu_p1  = -4.90680e-01;
		deltaPhi_mu_p2  =  6.14097e-02;
		deltaPhi_mu_p3  = -2.70582e-03;
		//---------------------------// 
		deltaPhi_sig_p0 =  6.26580e+00;
		deltaPhi_sig_p1 =  5.24376e-01;
		deltaPhi_sig_p2 = -1.21566e-01;
		deltaPhi_sig_p3 =  6.88625e-03;
		
		deltaK_mu_p0    =  5.04883e-01;
		deltaK_mu_p1    = -2.02498e-01;
		deltaK_mu_p2    =  2.62968e-02;
		deltaK_mu_p3    = -1.28288e-03;
		//---------------------------// 
		deltaK_sig_p0   =  3.15228e-02;
		deltaK_sig_p1   =  8.24456e-03;
		deltaK_sig_p2   =  2.83940e-03;
		deltaK_sig_p3   = -2.01125e-04;
		
		deltaK2_mu_p0   =  6.89191e-04;
		deltaK2_mu_p1   =  1.65777e-02;
		deltaK2_mu_p2   = -4.57181e-03;
		deltaK2_mu_p3   =  2.15276e-04;
		//---------------------------// 
		deltaK2_sig_p0  = -4.12018e-02;
		deltaK2_sig_p1  =  4.22817e-02;
		deltaK2_sig_p2  = -2.05753e-03;
		deltaK2_sig_p3  =  6.12380e-05;
		
	} else if(group < 12) { // He FIELD ON
		
		deltaE_mu_p0    = -1.94558e-01;
		deltaE_mu_p1    =  6.55092e-02;
		deltaE_mu_p2    = -6.18180e-03;
		deltaE_mu_p3    =  1.12281e-04;
		//---------------------------//
		deltaE_sig_p0   =  9.03329e-03;
		deltaE_sig_p1   =  4.26364e-02;
		deltaE_sig_p2   =  1.65108e-05;
		
		deltaPhi_mu_p0  =  1.80000e+02;
		deltaPhi_mu_p1  =  0.;
		deltaPhi_mu_p2  =  0.;
		deltaPhi_mu_p3  =  0.;
		//---------------------------// 
		deltaPhi_sig_p0 =  1.50000e+01;
		deltaPhi_sig_p1 =  0.0;
		deltaPhi_sig_p2 =  0.0;
		deltaPhi_sig_p3 =  0.0;
		
		deltaK_mu_p0    =  5.34489e-01;
		deltaK_mu_p1    = -2.12861e-01;
		deltaK_mu_p2    =  2.75625e-02;
		deltaK_mu_p3    = -1.32998e-03;
		//---------------------------//
		deltaK_sig_p0   =  1.11873e-02;
		deltaK_sig_p1   =  1.64423e-02;
		deltaK_sig_p2   =  1.67522e-03;
		deltaK_sig_p3   = -1.50906e-04;
		
		deltaK2_mu_p0   =  6.89191e-04;
		deltaK2_mu_p1   =  1.65777e-02;
		deltaK2_mu_p2   = -4.57181e-03;
		deltaK2_mu_p3   =  2.15276e-04;
		//---------------------------// 
		deltaK2_sig_p0  = -4.12018e-02;
		deltaK2_sig_p1  =  4.22817e-02;
		deltaK2_sig_p2  = -2.05753e-03;
		deltaK2_sig_p3  =  6.12380e-05;
		
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




