
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
	
	char loc_pathname[256] = "/work/halld/home/andrsmit/primex_compton_analysis";
	
	// Directory where ROOT Trees are stored:
	
	sprintf(rootTree_pathName, "%s/data/rootTrees/phaseIII", loc_pathname);
	
	// Directory where output ROOT files will be stored:
	
	sprintf(rootFile_pathName, "%s/analyze_trees/phaseIII/systematics/rootFiles", loc_pathname);
	
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
		
		//----------------------------------------------------------//
		
		//m_beamX = 0.146;
		//m_beamY = 0.017;
		m_beamX = 0.151;
		m_beamY = 0.012;
		if(run_group < 8) m_beamZ = 64.935;
		else 							m_beamZ = 65.0;
		
		cout << "Processing run " << irun << endl;
		
		load_constants(run_group, first_evt);
		first_evt = false;
		
		for(int iext = 0; iext < 200; iext++) {
			
			char buf[256];
			sprintf(buf,"%s/%d/%d_%03d.root", rootTree_pathName, irun, irun, iext);
			if(gSystem->AccessPathName(buf)) continue;
			
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
	if(n_candidates > 100) return;
	
	int    n_good_cands    =  0;
	int    single_tag_sys  = -1, single_tag_counter = -1;
	double single_deltaK   =  0.;
	
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
		
		
		double phi1 = atan2(y1,x1) * (180./TMath::Pi());
		double phi2 = atan2(y2,x2) * (180./TMath::Pi());
		
		int fcal_phi_sect;
		if(     -180. < phi1 && phi1 <= -135.) fcal_phi_sect = 0;
		else if(-135. < phi1 && phi1 <=  -90.) fcal_phi_sect = 1;
		else if( -90. < phi1 && phi1 <=  -45.) fcal_phi_sect = 2;
		else if( -45. < phi1 && phi1 <=    0.) fcal_phi_sect = 3;
		else if(   0. < phi1 && phi1 <=   45.) fcal_phi_sect = 4;
		else if(  45. < phi1 && phi1 <=   90.) fcal_phi_sect = 5;
		else if(  90. < phi1 && phi1 <=  135.) fcal_phi_sect = 6;
		else if( 135. < phi1 && phi1 <=  180.) fcal_phi_sect = 7;
		else fcal_phi_sect = 8;
		
		int fcal_layer;
		if((-4.0157*(1.5+1.) < fcal_face_x && fcal_face_x < 4.0157*(1.5+1.)) 
			&& (-4.0157*(1.5+1.) < fcal_face_y && fcal_face_y < 4.0157*(1.5+1.))) 
			fcal_layer = 1;
		else if((-4.0157*(1.5+2.) < fcal_face_x && fcal_face_x < 4.0157*(1.5+2.)) 
			&& (-4.0157*(1.5+2.) < fcal_face_y && fcal_face_y < 4.0157*(1.5+2.))) 
			fcal_layer = 2;
		else if((-4.0157*(1.5+3.) < fcal_face_x && fcal_face_x < 4.0157*(1.5+3.)) 
			&& (-4.0157*(1.5+3.) < fcal_face_y && fcal_face_y < 4.0157*(1.5+3.))) 
			fcal_layer = 3;
		else if((-4.0157*(1.5+4.) < fcal_face_x && fcal_face_x < 4.0157*(1.5+4.)) 
			&& (-4.0157*(1.5+4.) < fcal_face_y && fcal_face_y < 4.0157*(1.5+4.))) 
			fcal_layer = 4;
		else if((-4.0157*(1.5+5.) < fcal_face_x && fcal_face_x < 4.0157*(1.5+5.)) 
			&& (-4.0157*(1.5+5.) < fcal_face_y && fcal_face_y < 4.0157*(1.5+5.))) 
			fcal_layer = 5;
		else if((-4.0157*(1.5+6.) < fcal_face_x && fcal_face_x < 4.0157*(1.5+6.)) 
			&& (-4.0157*(1.5+6.) < fcal_face_y && fcal_face_y < 4.0157*(1.5+6.))) 
			fcal_layer = 6;
		else if((-4.0157*(1.5+7.) < fcal_face_x && fcal_face_x < 4.0157*(1.5+7.)) 
			&& (-4.0157*(1.5+7.) < fcal_face_y && fcal_face_y < 4.0157*(1.5+7.))) 
			fcal_layer = 7;
		else fcal_layer = 8;
		
		int ccal_phi_sect;
		if(     -180. < phi2 && phi2 <= -135.) ccal_phi_sect = 0;
		else if(-135. < phi2 && phi2 <=  -90.) ccal_phi_sect = 1;
		else if( -90. < phi2 && phi2 <=  -45.) ccal_phi_sect = 2;
		else if( -45. < phi2 && phi2 <=    0.) ccal_phi_sect = 3;
		else if(   0. < phi2 && phi2 <=   45.) ccal_phi_sect = 4;
		else if(  45. < phi2 && phi2 <=   90.) ccal_phi_sect = 5;
		else if(  90. < phi2 && phi2 <=  135.) ccal_phi_sect = 6;
		else if( 135. < phi2 && phi2 <=  180.) ccal_phi_sect = 7;
		else ccal_phi_sect = 8;
		
		int ccal_layer;
		if((-2.09*(1.+1.) < ccal_face_x && ccal_face_x < 2.09*(1.+1.)) 
			&& (-2.09*(1.+1.) < ccal_face_y && ccal_face_y < 2.09*(1.+1.))) 
			ccal_layer = 1;
		else if((-2.09*(1.+2.) < ccal_face_x && ccal_face_x < 2.09*(1.+2.)) 
			&& (-2.09*(1.+2.) < ccal_face_y && ccal_face_y < 2.09*(1.+2.))) 
			ccal_layer = 2;
		else if((-2.09*(1.+3.) < ccal_face_x && ccal_face_x < 2.09*(1.+3.)) 
			&& (-2.09*(1.+3.) < ccal_face_y && ccal_face_y < 2.09*(1.+3.))) 
			ccal_layer = 3;
		else if((-2.09*(1.+4.) < ccal_face_x && ccal_face_x < 2.09*(1.+4.)) 
			&& (-2.09*(1.+4.) < ccal_face_y && ccal_face_y < 2.09*(1.+4.))) 
			ccal_layer = 4;
		else ccal_layer = 5;
		
		if(fcal_phi_sect == 8 || ccal_phi_sect == 8) {
			cout << "fcal_phi = " << fcal_phi_sect << endl;
			cout << "ccal_phi = " << ccal_phi_sect << endl;
		}
		
		//--------------------------------------------//
		/*
		if((-32. < fcal_face_y && fcal_face_y < -20.) && 
			(-8. < fcal_face_x && fcal_face_x < 4.)) continue;
		*/
		//--------------------------------------------//
		
		double fill_weight = 0.;
		if(bunch_val[ic]) {
			if(fabs(tb[ic]-rfTime) < 2.004) fill_weight =  1.0;
		} else {
			fill_weight = -0.1;
		}
		
		//-----   Compton Cuts   -----//
		
		// re-calculate deltaK with new vertex/positions:
		
		double theta1 = atan2(sqrt(x1*x1 + y1*y1),z1);
		double theta2 = atan2(sqrt(x2*x2 + y2*y2),z2);
		
		//double ecomp1 = 1. / ((1./eb[ic]) + (1./m_e)*(1.-cos(theta1)));
		//double ecomp2 = 1. / ((1./eb[ic]) + (1./m_e)*(1.-cos(theta2)));
		//double loc_deltaK = (ecomp1 + ecomp2) - (eb[ic] + m_e);
		
		double ecomp      = m_e * sin(theta1+theta2) 
			/ (sin(theta1) + sin(theta2) - sin(theta1+theta2));
		double loc_deltaK = ecomp - eb[ic];
		
		double   deltaE_mu,   deltaE_sig;
		double deltaPhi_mu, deltaPhi_sig;
		double   deltaK_mu,   deltaK_sig;
		
		deltaE_mu    = f_deltaE_mu->Eval(eb[ic]);
		deltaE_sig   = eb[ic] * f_deltaE_sig->Eval(eb[ic]);
		deltaPhi_mu  = f_deltaPhi_mu->Eval(eb[ic]);
		deltaPhi_sig = f_deltaPhi_sig->Eval(eb[ic]);
		deltaK_mu    = f_deltaK_mu->Eval(eb[ic]);
		deltaK_sig   = f_deltaK_sig->Eval(eb[ic]);
		
		int phi_cut[16], e_cut[16], k_cut[16];
		for(int icut=0; icut<16; icut++) {
			double loc_cut = cut_sigmas[icut];
			
			if(fabs(deltaPhi[ic]-deltaPhi_mu) < loc_cut*deltaPhi_sig) {
				phi_cut[icut] = 1;
			} else {
				phi_cut[icut] = 0;
			}
			
			if(fabs(deltaE[ic]-deltaE_mu) < loc_cut*deltaE_sig) {
				e_cut[icut] = 1;
			} else {
				e_cut[icut] = 0;
			}
			
			if(fabs(loc_deltaK-deltaK_mu) < loc_cut*deltaK_sig) {
				k_cut[icut] = 1;
			} else {
				k_cut[icut] = 0;
			}
		}
		
		int e5_cut =   e_cut[6];
		int k5_cut =   k_cut[6];
		int p5_cut = phi_cut[6];
		
		//---------------------------------------------------------------------------//
		// RF Timing Cuts:
		
		int fcal_t_cut = 0, ccal_t_cut = 0;
		if(fabs(fcal_t[ic] - rfTime) < FCAL_RF_TIME_CUT) fcal_t_cut = 1;
		if(fabs(ccal_t[ic] - rfTime) < CCAL_RF_TIME_CUT) ccal_t_cut = 1;
		
		int fcal_t_cut_vec[10], ccal_t_cut_vec[10];
		for(int icut=0; icut<10; icut++) {
			double loc_cut = 0.5 * (double)(icut+1);
			
			if(fabs(fcal_t[ic]-rfTime) < loc_cut) {
				fcal_t_cut_vec[icut] = 1;
			} else {
				fcal_t_cut_vec[icut] = 0;
			}
			
			if(fabs(ccal_t[ic]-rfTime) < loc_cut) {
				ccal_t_cut_vec[icut] = 1;
			} else {
				ccal_t_cut_vec[icut] = 0;
			}
		}
		
		//---------------------------------------------------------------------------//
		// Minimum Energy Cuts:
		
		int fcal_e_cut = 0, ccal_e_cut = 0;
		if(fcal_e[ic] > FCAL_ENERGY_CUT) fcal_e_cut = 1;
		if(ccal_e[ic] > CCAL_ENERGY_CUT) ccal_e_cut = 1;
		
		//---------------------------------------------------------------------------//
		
		if(e5_cut && p5_cut && fcal_e_cut && ccal_e_cut && !fcal_fid_cut && !ccal_fid_cut
			&& fcal_t_cut && ccal_t_cut) 
		{
			n_good_cands++;
			single_deltaK      = loc_deltaK;
			single_tag_sys     = tag_sys[ic];
			single_tag_counter = tag_counter[ic];
			
			if(tag_sys[ic]==0) {
				h_deltaK_tagh_fcal_phi[fcal_phi_sect]->Fill(tag_counter[ic], 
					loc_deltaK, fill_weight);
				h_deltaK_tagh_ccal_phi[ccal_phi_sect]->Fill(tag_counter[ic], 
					loc_deltaK, fill_weight);
				h_deltaK_tagh_fcal_layer[fcal_layer-1]->Fill(tag_counter[ic], 
					loc_deltaK, fill_weight);
				h_deltaK_tagh_ccal_layer[ccal_layer-1]->Fill(tag_counter[ic], 
					loc_deltaK, fill_weight);
			} else {
				h_deltaK_tagm_fcal_phi[fcal_phi_sect]->Fill(tag_counter[ic], 
					loc_deltaK, fill_weight);
				h_deltaK_tagm_ccal_phi[ccal_phi_sect]->Fill(tag_counter[ic], 
					loc_deltaK, fill_weight);
				h_deltaK_tagm_fcal_layer[fcal_layer-1]->Fill(tag_counter[ic], 
					loc_deltaK, fill_weight);
				h_deltaK_tagm_ccal_layer[ccal_layer-1]->Fill(tag_counter[ic], 
					loc_deltaK, fill_weight);
			}
			
			if(k5_cut) {
			
				if(tag_sys[ic]==0) {
					h_deltaK_tagh_cut_fcal_phi[fcal_phi_sect]->Fill(
						tag_counter[ic], loc_deltaK, fill_weight);
					h_deltaK_tagh_cut_ccal_phi[ccal_phi_sect]->Fill(
						tag_counter[ic], loc_deltaK, fill_weight);
					h_deltaK_tagh_cut_fcal_layer[fcal_layer-1]->Fill(
						tag_counter[ic], loc_deltaK, fill_weight);
					h_deltaK_tagh_cut_ccal_layer[ccal_layer-1]->Fill(
						tag_counter[ic], loc_deltaK, fill_weight);
				} else {
					h_deltaK_tagm_cut_fcal_phi[fcal_phi_sect]->Fill(
						tag_counter[ic], loc_deltaK, fill_weight);
					h_deltaK_tagm_cut_ccal_phi[ccal_phi_sect]->Fill(
						tag_counter[ic], loc_deltaK, fill_weight);
					h_deltaK_tagm_cut_fcal_layer[fcal_layer-1]->Fill(
						tag_counter[ic], loc_deltaK, fill_weight);
					h_deltaK_tagm_cut_ccal_layer[ccal_layer-1]->Fill(
						tag_counter[ic], loc_deltaK, fill_weight);
				}
				
				h_xy_fcal_phi[fcal_phi_sect]->Fill(fcal_face_x, fcal_face_y, 
					fill_weight);
				h_xy_ccal_phi[ccal_phi_sect]->Fill(ccal_face_x, ccal_face_y, 
					fill_weight);
				h_xy_fcal_layer[fcal_layer-1]->Fill(fcal_face_x, fcal_face_y, 
					fill_weight);
				h_xy_ccal_layer[ccal_layer-1]->Fill(ccal_face_x, ccal_face_y, 
					fill_weight);
				
				double loc_brfdt = tb[ic]-rfTime;
				if(run>61910 && run<61947) loc_brfdt += 2.0;
				else if(run==61952)        loc_brfdt += 2.0;
				
				if(tag_sys[ic]==0) {
					if(bunch_val[ic]) {
						h_deltaK_tagh_main_extra->Fill(tag_counter[ic], 
							loc_deltaK, 1.00);
						if(fill_weight==1.0) {
							h_deltaK_tagh_main->Fill(tag_counter[ic], 
								loc_deltaK, 1.00);
						}
					}
					else {
						if((-26.052 < loc_brfdt < -22.044) || 
							(22.044 < loc_brfdt < 26.052)) {
							h_deltaK_tagh_acc_single->Fill(tag_counter[ic], 								loc_deltaK, 0.5);
						}
						if((-26.052 < loc_brfdt < -18.036) || 
							(18.036 < loc_brfdt < 26.052)) {
							h_deltaK_tagh_acc->Fill(tag_counter[ic], 
								loc_deltaK, 0.25);
						}
						if((-26.052 < loc_brfdt < -18.036) || 
							(18.036 < loc_brfdt < 26.052)) {
							h_deltaK_tagh_acc_extra->Fill(tag_counter[ic], 
								loc_deltaK, 0.75);
						}
					}
				} else {
					if(bunch_val[ic]) {
						h_deltaK_tagm_main_extra->Fill(tag_counter[ic], 
							loc_deltaK, 1.00);
						if(fill_weight==1.0) {
							h_deltaK_tagm_main->Fill(tag_counter[ic], 
								loc_deltaK, 1.00);
						}
					}
					else {
						if((-26.052 < loc_brfdt < -22.044) || 
							(22.044 < loc_brfdt < 26.052)) {
							h_deltaK_tagm_acc_single->Fill(tag_counter[ic], 								loc_deltaK, 0.5);
						}
						if((-26.052 < loc_brfdt < -18.036) || 
							(18.036 < loc_brfdt < 26.052)) {
							h_deltaK_tagm_acc->Fill(tag_counter[ic], 
								loc_deltaK, 0.25);
						}
						if((-26.052 < loc_brfdt < -18.036) || 
							(18.036 < loc_brfdt < 26.052)) {
							h_deltaK_tagm_acc_extra->Fill(tag_counter[ic], 
								loc_deltaK, 0.75);
						}
					}
				}
			}
		}
		
		
		
		if(tag_sys[ic]==0) { // TAGH
			
			if(!fcal_fid_cut && !ccal_fid_cut) {
				
				if(e5_cut && p5_cut && ccal_e_cut && fcal_t_cut && ccal_t_cut) {
					
					for(int icut=0; icut<20; icut++) {
						double loc_cut = 0.05 * (double)(icut);
						if(fcal_e[ic] > loc_cut) {
							h_deltaK_tagh_fcalE[icut]->Fill(tag_counter[ic], loc_deltaK, fill_weight);
							if(k5_cut) {
								h_deltaK_tagh_cut_fcalE[icut]->Fill(tag_counter[ic], loc_deltaK, fill_weight);
							}
						}
					}
				}
				
				if(e5_cut && p5_cut && fcal_e_cut && fcal_t_cut && ccal_t_cut) {
					
					for(int icut=0; icut<13; icut++) {
						double loc_cut = 0.05 * (double)(icut);
						if(ccal_e[ic] > loc_cut) {
							h_deltaK_tagh_ccalE[icut]->Fill(tag_counter[ic], loc_deltaK, fill_weight);
							if(k5_cut) {
								h_deltaK_tagh_cut_ccalE[icut]->Fill(tag_counter[ic], loc_deltaK, fill_weight);
							}
						}
					}
				}
				
				if(e5_cut && fcal_e_cut && ccal_e_cut && fcal_t_cut && ccal_t_cut) {
					
					for(int icut=0; icut<16; icut++) {
						if(phi_cut[icut]) {
							h_deltaK_tagh_sigPhi[icut]->Fill(tag_counter[ic], loc_deltaK, fill_weight);
							if(k5_cut) {
								h_deltaK_tagh_cut_sigPhi[icut]->Fill(tag_counter[ic], loc_deltaK, fill_weight);
							}
						}
					}
				}
				
				if(p5_cut && fcal_e_cut && ccal_e_cut && fcal_t_cut && ccal_t_cut) {
					
					for(int icut=0; icut<16; icut++) {
						if(e_cut[icut]) {
							h_deltaK_tagh_sigE[icut]->Fill(tag_counter[ic], loc_deltaK, fill_weight);
							if(k5_cut) {
								h_deltaK_tagh_cut_sigE[icut]->Fill(tag_counter[ic], loc_deltaK, fill_weight);
							}
						}
					}
				}
				
				if(e5_cut && p5_cut && fcal_e_cut && ccal_e_cut && ccal_t_cut) {
					
					for(int icut=0; icut<10; icut++) {
						if(fcal_t_cut_vec[icut]) {
							h_deltaK_tagh_fcalT[icut]->Fill(tag_counter[ic], loc_deltaK, fill_weight);
							if(k5_cut) {
								h_deltaK_tagh_cut_fcalT[icut]->Fill(tag_counter[ic], loc_deltaK, fill_weight);
							}
						}
					}
				}
				
				if(e5_cut && p5_cut && fcal_e_cut && ccal_e_cut && fcal_t_cut) {
					
					for(int icut=0; icut<10; icut++) {
						if(ccal_t_cut_vec[icut]) {
							h_deltaK_tagh_ccalT[icut]->Fill(tag_counter[ic], loc_deltaK, fill_weight);
							if(k5_cut) {
								h_deltaK_tagh_cut_ccalT[icut]->Fill(tag_counter[ic], loc_deltaK, fill_weight);
							}
						}
					}
				}
				
				if(p5_cut && e5_cut && fcal_e_cut && ccal_e_cut && fcal_t_cut && ccal_t_cut) {
					
					for(int icut=0; icut<16; icut++) {
						if(k_cut[icut]) {
							h_deltaK_tagh_sigK[icut]->Fill(tag_counter[ic], loc_deltaK, fill_weight);
						}
					}
				}
			} // if(!fcal_fid_cut && !ccal_fid_cut)
			
		} else { // TAGM
			
			if(!fcal_fid_cut && !ccal_fid_cut) {
				
				if(e5_cut && p5_cut && ccal_e_cut && fcal_t_cut && ccal_t_cut) {
					
					for(int icut=0; icut<20; icut++) {
						double loc_cut = 0.05 * (double)(icut);
						if(fcal_e[ic] > loc_cut) {
							h_deltaK_tagm_fcalE[icut]->Fill(tag_counter[ic], loc_deltaK, fill_weight);
							if(k5_cut) {
								h_deltaK_tagm_cut_fcalE[icut]->Fill(tag_counter[ic], loc_deltaK, fill_weight);
							}
						}
					}
				}
				
				if(e5_cut && p5_cut && fcal_e_cut && fcal_t_cut && ccal_t_cut) {
					
					for(int icut=0; icut<13; icut++) {
						double loc_cut = 0.05 * (double)(icut);
						if(ccal_e[ic] > loc_cut) {
							h_deltaK_tagm_ccalE[icut]->Fill(tag_counter[ic], loc_deltaK, fill_weight);
							if(k5_cut) {
								h_deltaK_tagm_cut_ccalE[icut]->Fill(tag_counter[ic], loc_deltaK, fill_weight);
							}
						}
					}
				}
				
				if(e5_cut && fcal_e_cut && ccal_e_cut && fcal_t_cut && ccal_t_cut) {
					
					for(int icut=0; icut<16; icut++) {
						if(phi_cut[icut]) {
							h_deltaK_tagm_sigPhi[icut]->Fill(tag_counter[ic], loc_deltaK, fill_weight);
							if(k5_cut) {
								h_deltaK_tagm_cut_sigPhi[icut]->Fill(tag_counter[ic], loc_deltaK, fill_weight);
							}
						}
					}
				}
				
				if(p5_cut && fcal_e_cut && ccal_e_cut && fcal_t_cut && ccal_t_cut) {
					
					for(int icut=0; icut<16; icut++) {
						if(e_cut[icut]) {
							h_deltaK_tagm_sigE[icut]->Fill(tag_counter[ic], loc_deltaK, fill_weight);
							if(k5_cut) {
								h_deltaK_tagm_cut_sigE[icut]->Fill(tag_counter[ic], loc_deltaK, fill_weight);
							}
						}
					}
				}
				
				if(e5_cut && p5_cut && fcal_e_cut && ccal_e_cut && ccal_t_cut) {
					
					for(int icut=0; icut<10; icut++) {
						if(fcal_t_cut_vec[icut]) {
							h_deltaK_tagm_fcalT[icut]->Fill(tag_counter[ic], loc_deltaK, fill_weight);
							if(k5_cut) {
								h_deltaK_tagm_cut_fcalT[icut]->Fill(tag_counter[ic], loc_deltaK, fill_weight);
							}
						}
					}
				}
				
				if(e5_cut && p5_cut && fcal_e_cut && ccal_e_cut && fcal_t_cut) {
					
					for(int icut=0; icut<10; icut++) {
						if(ccal_t_cut_vec[icut]) {
							h_deltaK_tagm_ccalT[icut]->Fill(tag_counter[ic], loc_deltaK, fill_weight);
							if(k5_cut) {
								h_deltaK_tagm_cut_ccalT[icut]->Fill(tag_counter[ic], loc_deltaK, fill_weight);
							}
						}
					}
				}
				
				if(p5_cut && e5_cut && fcal_e_cut && ccal_e_cut && fcal_t_cut && ccal_t_cut) {
					
					for(int icut=0; icut<16; icut++) {
						if(k_cut[icut]) {
							h_deltaK_tagm_sigK[icut]->Fill(tag_counter[ic], loc_deltaK, fill_weight);
						}
					}
				}
			} // if(!fcal_fid_cut && !ccal_fid_cut)
			
		}
		
		
		if(e5_cut && p5_cut && k5_cut && fcal_e_cut && ccal_e_cut 
			&& fcal_t_cut && ccal_t_cut) {
			
			for(int icut = 0; icut < N_FID_CUTS; icut++) {
				
				int loc_fcal_fid_cut = 0;
				double loc_fcal_inner_layer_cut = 2.0*4.0157 + 
					(4.0157/10.)*(double)icut;
				if((-1.*loc_fcal_inner_layer_cut < fcal_face_x
					&& fcal_face_x < loc_fcal_inner_layer_cut)
					&& (-1.*loc_fcal_inner_layer_cut < fcal_face_y 
					&& fcal_face_y < loc_fcal_inner_layer_cut)) 
						loc_fcal_fid_cut = 1;
				
				if(!loc_fcal_fid_cut && !ccal_fid_cut) {
					
					if(tag_sys[ic]==0) 
						h_deltaK_tagh_fcalfid[icut]->Fill(tag_counter[ic], 
							loc_deltaK, fill_weight);
					else 
						h_deltaK_tagm_fcalfid[icut]->Fill(tag_counter[ic], 
							loc_deltaK, fill_weight);
				}
			}
			
			for(int icut = 0; icut < N_FID_CUTS; icut++) {
				
				int loc_ccal_fid_cut = 0;
				double loc_ccal_inner_layer_cut = 1.0*2.09 + (2.09/10.)*(double)icut;
				if((-1.*loc_ccal_inner_layer_cut < ccal_face_x
					&& ccal_face_x < loc_ccal_inner_layer_cut)
					&& (-1.*loc_ccal_inner_layer_cut < ccal_face_y 
					&& ccal_face_y < loc_ccal_inner_layer_cut)) 
						loc_ccal_fid_cut = 1;
				
				if(!loc_ccal_fid_cut && !fcal_fid_cut) {
					
					if(tag_sys[ic]==0) 
						h_deltaK_tagh_ccalfid[icut]->Fill(tag_counter[ic], 
							loc_deltaK, fill_weight);
					else 
						h_deltaK_tagm_ccalfid[icut]->Fill(tag_counter[ic], 
							loc_deltaK, fill_weight);
				}
			}
		}
	} // for( int ic = 0; ic < n_candidates; ic++ )
	
	if(n_good_cands==1) {
		if(single_tag_sys==0) 
			h_deltaK_tagh_single->Fill(single_tag_counter, single_deltaK, 1.0);
		else
			h_deltaK_tagm_single->Fill(single_tag_counter, single_deltaK, 1.0);
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
	
	for(int ihist=0; ihist<20; ihist++) {
		h_deltaK_tagh_fcalE[ihist]->Reset();
		h_deltaK_tagm_fcalE[ihist]->Reset();
		h_deltaK_tagh_cut_fcalE[ihist]->Reset();
		h_deltaK_tagm_cut_fcalE[ihist]->Reset();
	}
	for(int ihist=0; ihist<13; ihist++) {
		h_deltaK_tagh_ccalE[ihist]->Reset();
		h_deltaK_tagm_ccalE[ihist]->Reset();
		h_deltaK_tagh_cut_ccalE[ihist]->Reset();
		h_deltaK_tagm_cut_ccalE[ihist]->Reset();
	}
	for(int ihist=0; ihist<16; ihist++) {
		h_deltaK_tagh_sigE[ihist]->Reset();
		h_deltaK_tagm_sigE[ihist]->Reset();
		h_deltaK_tagh_cut_sigE[ihist]->Reset();
		h_deltaK_tagm_cut_sigE[ihist]->Reset();
		
		h_deltaK_tagh_sigPhi[ihist]->Reset();
		h_deltaK_tagm_sigPhi[ihist]->Reset();
		h_deltaK_tagh_cut_sigPhi[ihist]->Reset();
		h_deltaK_tagm_cut_sigPhi[ihist]->Reset();
		
		h_deltaK_tagh_sigK[ihist]->Reset();
		h_deltaK_tagm_sigK[ihist]->Reset();
	}
	for(int ihist=0; ihist<10; ihist++) {
		h_deltaK_tagh_fcalT[ihist]->Reset();
		h_deltaK_tagm_fcalT[ihist]->Reset();
		h_deltaK_tagh_cut_fcalT[ihist]->Reset();
		h_deltaK_tagm_cut_fcalT[ihist]->Reset();
		
		h_deltaK_tagh_ccalT[ihist]->Reset();
		h_deltaK_tagm_ccalT[ihist]->Reset();
		h_deltaK_tagh_cut_ccalT[ihist]->Reset();
		h_deltaK_tagm_cut_ccalT[ihist]->Reset();
	}
	
	for( int icut = 0; icut < N_FID_CUTS; icut++ ) {
		h_deltaK_tagh_fcalfid[icut]->Reset();
		h_deltaK_tagm_fcalfid[icut]->Reset();
		h_deltaK_tagh_ccalfid[icut]->Reset();
		h_deltaK_tagm_ccalfid[icut]->Reset();
	}
	
	for( int ip = 0; ip < 8; ip++ ) {
		h_deltaK_tagh_fcal_phi[ip]->Reset();
		h_deltaK_tagm_fcal_phi[ip]->Reset();
		h_deltaK_tagh_ccal_phi[ip]->Reset();
		h_deltaK_tagm_ccal_phi[ip]->Reset();
		
		h_deltaK_tagh_cut_fcal_phi[ip]->Reset();
		h_deltaK_tagm_cut_fcal_phi[ip]->Reset();
		h_deltaK_tagh_cut_ccal_phi[ip]->Reset();
		h_deltaK_tagm_cut_ccal_phi[ip]->Reset();
		
		h_xy_fcal_phi[ip]->Reset();
		h_xy_ccal_phi[ip]->Reset();
	}
	
	for( int il = 0; il < 8; il++ ) {
		h_deltaK_tagh_fcal_layer[il]->Reset();
		h_deltaK_tagm_fcal_layer[il]->Reset();
		
		h_deltaK_tagh_cut_fcal_layer[il]->Reset();
		h_deltaK_tagm_cut_fcal_layer[il]->Reset();
		
		h_xy_fcal_layer[il]->Reset();
	}
	for( int il = 0; il < 5; il++ ) {
		h_deltaK_tagh_ccal_layer[il]->Reset();
		h_deltaK_tagm_ccal_layer[il]->Reset();
		
		h_deltaK_tagh_cut_ccal_layer[il]->Reset();
		h_deltaK_tagm_cut_ccal_layer[il]->Reset();
		
		h_xy_ccal_layer[il]->Reset();
	}
	
	h_deltaK_tagh_single->Reset();
	h_deltaK_tagm_single->Reset();
	
	h_deltaK_tagh_main->Reset();
	h_deltaK_tagm_main->Reset();
	h_deltaK_tagh_main_extra->Reset();
	h_deltaK_tagm_main_extra->Reset();
	h_deltaK_tagh_acc->Reset();
	h_deltaK_tagm_acc->Reset();
	h_deltaK_tagh_acc_single->Reset();
	h_deltaK_tagm_acc_single->Reset();
	h_deltaK_tagh_acc_extra->Reset();
	h_deltaK_tagm_acc_extra->Reset();
	
	return;
}


void write_histograms(char *name) {
	
	TFile *fOut = new TFile(name,"RECREATE");
	fOut->cd();
	
	TDirectory *dir_fcalE = new TDirectoryFile("fcalE", "fcalE");
	dir_fcalE->cd();
	for(int ihist=0; ihist<20; ihist++) {
		h_deltaK_tagh_fcalE[ihist]->Write();
		h_deltaK_tagm_fcalE[ihist]->Write();
		h_deltaK_tagh_cut_fcalE[ihist]->Write();
		h_deltaK_tagm_cut_fcalE[ihist]->Write();
	}
	dir_fcalE->cd("../");
	
	TDirectory *dir_ccalE = new TDirectoryFile("ccalE", "ccalE");
	dir_ccalE->cd();
	for(int ihist=0; ihist<13; ihist++) {
		h_deltaK_tagh_ccalE[ihist]->Write();
		h_deltaK_tagm_ccalE[ihist]->Write();
		h_deltaK_tagh_cut_ccalE[ihist]->Write();
		h_deltaK_tagm_cut_ccalE[ihist]->Write();
	}
	dir_ccalE->cd("../");
	
	TDirectory *dir_fcalT = new TDirectoryFile("fcalT", "fcalT");
	dir_fcalT->cd();
	for(int ihist=0; ihist<10; ihist++) {
		h_deltaK_tagh_fcalT[ihist]->Write();
		h_deltaK_tagm_fcalT[ihist]->Write();
		h_deltaK_tagh_cut_fcalT[ihist]->Write();
		h_deltaK_tagm_cut_fcalT[ihist]->Write();
	}
	dir_fcalT->cd("../");
	
	TDirectory *dir_ccalT = new TDirectoryFile("ccalT", "ccalT");
	dir_ccalT->cd();
	for(int ihist=0; ihist<10; ihist++) {
		h_deltaK_tagh_ccalT[ihist]->Write();
		h_deltaK_tagm_ccalT[ihist]->Write();
		h_deltaK_tagh_cut_ccalT[ihist]->Write();
		h_deltaK_tagm_cut_ccalT[ihist]->Write();
	}
	dir_ccalT->cd("../");
	
	TDirectory *dir_sigE = new TDirectoryFile("DeltaE", "DeltaE");
	dir_sigE->cd();
	for(int ihist=0; ihist<16; ihist++) {
		h_deltaK_tagh_sigE[ihist]->Write();
		h_deltaK_tagm_sigE[ihist]->Write();
		h_deltaK_tagh_cut_sigE[ihist]->Write();
		h_deltaK_tagm_cut_sigE[ihist]->Write();
	}
	dir_sigE->cd("../");
	
	TDirectory *dir_sigPhi = new TDirectoryFile("DeltaPhi", "DeltaPhi");
	dir_sigPhi->cd();
	for(int ihist=0; ihist<16; ihist++) {
		h_deltaK_tagh_sigPhi[ihist]->Write();
		h_deltaK_tagm_sigPhi[ihist]->Write();
		h_deltaK_tagh_cut_sigPhi[ihist]->Write();
		h_deltaK_tagm_cut_sigPhi[ihist]->Write();
	}
	dir_sigPhi->cd("../");
	
	TDirectory *dir_sigK = new TDirectoryFile("DeltaK", "DeltaK");
	dir_sigK->cd();
	for(int ihist=0; ihist<16; ihist++) {
		h_deltaK_tagh_sigK[ihist]->Write();
		h_deltaK_tagm_sigK[ihist]->Write();
	}
	dir_sigK->cd("../");
	
	TDirectory *dir_fcal_fid = new TDirectoryFile("fcal_fid", "fcal_fid");
	dir_fcal_fid->cd();
	for(int icut = 0; icut < N_FID_CUTS; icut++) {
		h_deltaK_tagh_fcalfid[icut]->Write();
		h_deltaK_tagm_fcalfid[icut]->Write();
	}
	dir_fcal_fid->cd("../");
	
	TDirectory *dir_ccal_fid = new TDirectoryFile("ccal_fid", "ccal_fid");
	dir_ccal_fid->cd();
	for(int icut = 0; icut < N_FID_CUTS; icut++) {
		h_deltaK_tagh_ccalfid[icut]->Write();
		h_deltaK_tagm_ccalfid[icut]->Write();
	}
	dir_ccal_fid->cd("../");
	
	TDirectory *dir_fcal_phi = new TDirectoryFile("fcal_phi", "fcal_phi");
	dir_fcal_phi->cd();
	for(int ip = 0; ip < 8; ip++) {
		h_deltaK_tagh_fcal_phi[ip]->Write();
		h_deltaK_tagm_fcal_phi[ip]->Write();
		
		h_deltaK_tagh_cut_fcal_phi[ip]->Write();
		h_deltaK_tagm_cut_fcal_phi[ip]->Write();
		
		h_xy_fcal_phi[ip]->Write();
	}
	dir_fcal_phi->cd("../");
	
	TDirectory *dir_ccal_phi = new TDirectoryFile("ccal_phi", "ccal_phi");
	dir_ccal_phi->cd();
	for(int ip = 0; ip < 8; ip++) {
		h_deltaK_tagh_ccal_phi[ip]->Write();
		h_deltaK_tagm_ccal_phi[ip]->Write();
		
		h_deltaK_tagh_cut_ccal_phi[ip]->Write();
		h_deltaK_tagm_cut_ccal_phi[ip]->Write();
		
		h_xy_ccal_phi[ip]->Write();
	}
	dir_ccal_phi->cd("../");
	
	TDirectory *dir_fcal_layer = new TDirectoryFile("fcal_layer", "fcal_layer");
	dir_fcal_layer->cd();
	for( int il = 0; il < 8; il++ ) {
		h_deltaK_tagh_fcal_layer[il]->Write();
		h_deltaK_tagm_fcal_layer[il]->Write();
		
		h_deltaK_tagh_cut_fcal_layer[il]->Write();
		h_deltaK_tagm_cut_fcal_layer[il]->Write();
		
		h_xy_fcal_layer[il]->Write();
	}
	dir_fcal_layer->cd("../");
	
	TDirectory *dir_ccal_layer = new TDirectoryFile("ccal_layer", "ccal_layer");
	dir_ccal_layer->cd();
	for(int il = 0; il < 5; il++) {
		h_deltaK_tagh_ccal_layer[il]->Write();
		h_deltaK_tagm_ccal_layer[il]->Write();
		
		h_deltaK_tagh_cut_ccal_layer[il]->Write();
		h_deltaK_tagm_cut_ccal_layer[il]->Write();
		
		h_xy_ccal_layer[il]->Write();
	}
	dir_ccal_layer->cd("../");
	
	h_deltaK_tagh_single->Write();
	h_deltaK_tagm_single->Write();
	
	h_deltaK_tagh_main->Write();
	h_deltaK_tagm_main->Write();
	h_deltaK_tagh_main_extra->Write();
	h_deltaK_tagm_main_extra->Write();
	h_deltaK_tagh_acc->Write();
	h_deltaK_tagm_acc->Write();
	h_deltaK_tagh_acc_single->Write();
	h_deltaK_tagm_acc_single->Write();
	h_deltaK_tagh_acc_extra->Write();
	h_deltaK_tagm_acc_extra->Write();
	
	fOut->Write();
	
	return;
}


void init_histograms() 
{
	//-------------------------------------------------------------------------------//
	// FCAL Energy Cut
	
	for(int ihist=0; ihist<20; ihist++) {
		
		int loc_cut_int = 5*ihist;
		double loc_cut  = 0.05 * (double)(ihist);
		
		h_deltaK_tagh_fcalE[ihist] = new TH2F(Form("deltaK_tagh_%02dfcalE", loc_cut_int), 
			Form("#DeltaK (E_{FCAL} > %.2f GeV)", loc_cut), 274, 0.5, 274.5, 2000, -8., 8.);
		h_deltaK_tagm_fcalE[ihist] = new TH2F(Form("deltaK_tagm_%02dfcalE", loc_cut_int), 
			Form("#DeltaK (E_{FCAL} > %.2f GeV)", loc_cut), 102, 0.5, 102.5, 2000, -8., 8.);
		
		h_deltaK_tagh_cut_fcalE[ihist] = new TH2F(Form("deltaK_tagh_cut_%02dfcalE", loc_cut_int), 
			Form("#DeltaK (E_{FCAL} > %.2f GeV)", loc_cut), 274, 0.5, 274.5, 2000, -8., 8.);
		h_deltaK_tagm_cut_fcalE[ihist] = new TH2F(Form("deltaK_tagm_cut_%02dfcalE", loc_cut_int), 
			Form("#DeltaK (E_{FCAL} > %.2f GeV)", loc_cut), 102, 0.5, 102.5, 2000, -8., 8.);
	}
	
	//-------------------------------------------------------------------------------//
	// CCAL Energy Cut
	
	for(int ihist=0; ihist<13; ihist++) {
		
		int loc_cut_int = 5*ihist;
		double loc_cut  = 0.5 * (double)(ihist);
		
		h_deltaK_tagh_ccalE[ihist] = new TH2F(Form("deltaK_tagh_%02dccalE", loc_cut_int), 
			Form("#DeltaK (E_{CCAL} > %.2f GeV)", loc_cut), 274, 0.5, 274.5, 2000, -8., 8.);
		h_deltaK_tagm_ccalE[ihist] = new TH2F(Form("deltaK_tagm_%02dccalE", loc_cut_int), 
			Form("#DeltaK (E_{CCAL} > %.2f GeV)", loc_cut), 102, 0.5, 102.5, 2000, -8., 8.);
		
		h_deltaK_tagh_cut_ccalE[ihist] = new TH2F(Form("deltaK_tagh_cut_%02dccalE", loc_cut_int), 
			Form("#DeltaK (E_{CCAL} > %.2f GeV)", loc_cut), 274, 0.5, 274.5, 2000, -8., 8.);
		h_deltaK_tagm_cut_ccalE[ihist] = new TH2F(Form("deltaK_tagm_cut_%02dccalE", loc_cut_int), 
			Form("#DeltaK (E_{CCAL} > %.2f GeV)", loc_cut), 102, 0.5, 102.5, 2000, -8., 8.);
	}
	
	//-------------------------------------------------------------------------------//
	// FCAL Timing Cut
	
	for(int ihist=0; ihist<10; ihist++) {
		
		int loc_cut_int = 5*(ihist+1);
		double loc_cut  = 0.5 * (double)(ihist+1);
		
		h_deltaK_tagh_fcalT[ihist] = new TH2F(Form("deltaK_tagh_%02dfcalT", loc_cut_int), 
			Form("#DeltaK (|t_{FCAL}-t_{RF}| < %.1f ns)", loc_cut), 274, 0.5, 274.5, 2000, -8., 8.);
		h_deltaK_tagm_fcalT[ihist] = new TH2F(Form("deltaK_tagm_%02dfcalT", loc_cut_int), 
			Form("#DeltaK (|t_{FCAL}-t_{RF}| < %.1f ns)", loc_cut), 102, 0.5, 102.5, 2000, -8., 8.);
		
		h_deltaK_tagh_cut_fcalT[ihist] = new TH2F(Form("deltaK_tagh_cut_%02dfcalT", loc_cut_int), 
			Form("#DeltaK (|t_{FCAL}-t_{RF}| < %.1f ns)", loc_cut), 274, 0.5, 274.5, 2000, -8., 8.);
		h_deltaK_tagm_cut_fcalT[ihist] = new TH2F(Form("deltaK_tagm_cut_%02dfcalT", loc_cut_int), 
			Form("#DeltaK (|t_{FCAL}-t_{RF}| < %.1f ns)", loc_cut), 102, 0.5, 102.5, 2000, -8., 8.);
	}
	
	//-------------------------------------------------------------------------------//
	// CCAL Timing Cut
	
	for(int ihist=0; ihist<10; ihist++) {
		
		int loc_cut_int = 5*ihist;
		double loc_cut  = 0.5 * (double)(ihist);
		
		h_deltaK_tagh_ccalT[ihist] = new TH2F(Form("deltaK_tagh_%02dccalT", loc_cut_int), 
			Form("#DeltaK (|t_{CCAL}-t_{RF}| < %.1f ns)", loc_cut), 274, 0.5, 274.5, 2000, -8., 8.);
		h_deltaK_tagm_ccalT[ihist] = new TH2F(Form("deltaK_tagm_%02dccalT", loc_cut_int), 
			Form("#DeltaK (|t_{CCAL}-t_{RF}| < %.1f ns)", loc_cut), 102, 0.5, 102.5, 2000, -8., 8.);
		
		h_deltaK_tagh_cut_ccalT[ihist] = new TH2F(Form("deltaK_tagh_cut_%02dccalT", loc_cut_int), 
			Form("#DeltaK (|t_{CCAL}-t_{RF}| < %.1f ns)", loc_cut), 274, 0.5, 274.5, 2000, -8., 8.);
		h_deltaK_tagm_cut_ccalT[ihist] = new TH2F(Form("deltaK_tagm_cut_%02dccalT", loc_cut_int), 
			Form("#DeltaK (|t_{CCAL}-t_{RF}| < %.1f ns)", loc_cut), 102, 0.5, 102.5, 2000, -8., 8.);
	}
	
	//-------------------------------------------------------------------------------//
	// DeltaE Cut
	
	for(int ihist=0; ihist<16; ihist++) {
		
		double loc_cut  = cut_sigmas[ihist];
		int loc_cut_int = (int)(loc_cut*10.);
		
		h_deltaK_tagh_sigE[ihist] = new TH2F(Form("deltaK_tagh_%03dsigE", loc_cut_int), 
			Form("#DeltaK (%.1f#sigma #DeltaE Cut)", loc_cut), 274, 0.5, 274.5, 2000, -8., 8.);
		h_deltaK_tagm_sigE[ihist] = new TH2F(Form("deltaK_tagm_%03dsigE", loc_cut_int), 
			Form("#DeltaK (%.1f#sigma #DeltaE Cut)", loc_cut), 102, 0.5, 102.5, 2000, -8., 8.);
		
		h_deltaK_tagh_cut_sigE[ihist] = new TH2F(Form("deltaK_tagh_cut_%03dsigE", loc_cut_int), 
			Form("#DeltaK (%.1f#sigma #DeltaE Cut)", loc_cut), 274, 0.5, 274.5, 2000, -8., 8.);
		h_deltaK_tagm_cut_sigE[ihist] = new TH2F(Form("deltaK_tagm_cut_%03dsigE", loc_cut_int), 
			Form("#DeltaK (%.1f#sigma #DeltaE Cut)", loc_cut), 102, 0.5, 102.5, 2000, -8., 8.);
	}
	
	//-------------------------------------------------------------------------------//
	// DeltaPhi Cut
	
	for(int ihist=0; ihist<16; ihist++) {
		
		double loc_cut  = cut_sigmas[ihist];
		int loc_cut_int = (int)(loc_cut*10.);
		
		h_deltaK_tagh_sigPhi[ihist] = new TH2F(Form("deltaK_tagh_%03dsigPhi", loc_cut_int), 
			Form("#DeltaK (%.1f#sigma #Delta#phi Cut)", loc_cut), 274, 0.5, 274.5, 2000, -8., 8.);
		h_deltaK_tagm_sigPhi[ihist] = new TH2F(Form("deltaK_tagm_%03dsigPhi", loc_cut_int), 
			Form("#DeltaK (%.1f#sigma #Delta#phi Cut)", loc_cut), 102, 0.5, 102.5, 2000, -8., 8.);
		
		h_deltaK_tagh_cut_sigPhi[ihist] = new TH2F(Form("deltaK_tagh_cut_%03dsigPhi", loc_cut_int), 
			Form("#DeltaK (%.1f#sigma #Delta#phi Cut)", loc_cut), 274, 0.5, 274.5, 2000, -8., 8.);
		h_deltaK_tagm_cut_sigPhi[ihist] = new TH2F(Form("deltaK_tagm_cut_%03dsigPhi", loc_cut_int), 
			Form("#DeltaK (%.1f#sigma #Delta#phi Cut)", loc_cut), 102, 0.5, 102.5, 2000, -8., 8.);
	}
	
	//-------------------------------------------------------------------------------//
	// DeltaK Cut
	
	for(int ihist=0; ihist<16; ihist++) {
		
		double loc_cut  = cut_sigmas[ihist];
		int loc_cut_int = (int)(loc_cut*10.);
		
		h_deltaK_tagh_sigK[ihist] = new TH2F(Form("deltaK_tagh_%03dsigK", loc_cut_int), 
			Form("#DeltaK (%.1f#sigma #DeltaK Cut)", loc_cut), 274, 0.5, 274.5, 2000, -8., 8.);
		h_deltaK_tagm_sigK[ihist] = new TH2F(Form("deltaK_tagm_%03dsigK", loc_cut_int), 
			Form("#DeltaK (%.1f#sigma #DeltaK Cut)", loc_cut), 102, 0.5, 102.5, 2000, -8., 8.);
	}
	
	//-------------------------------------------------------------------------------//
	// Fiducial Cuts
	
	for(int icut = 0; icut < N_FID_CUTS; icut++) {
		
		h_deltaK_tagh_fcalfid[icut] = new TH2F(Form("deltaK_tagh_fcalfid_%02d",icut), 
			"#DeltaK", 274, 0.5, 274.5, 2000, -8., 8.);
		h_deltaK_tagm_fcalfid[icut] = new TH2F(Form("deltaK_tagm_fcalfid_%02d",icut), 
			"#DeltaK", 102, 0.5, 102.5, 2000, -8., 8.);
	}
	
	for(int icut = 0; icut < N_FID_CUTS; icut++) {
		
		h_deltaK_tagh_ccalfid[icut] = new TH2F(Form("deltaK_tagh_ccalfid_%02d",icut), 
			"#DeltaK", 274, 0.5, 274.5, 2000, -8., 8.);
		h_deltaK_tagm_ccalfid[icut] = new TH2F(Form("deltaK_tagm_ccalfid_%02d",icut), 
			"#DeltaK", 102, 0.5, 102.5, 2000, -8., 8.);
	}
	
	for(int ip = 0; ip < 8; ip++) {
		
		h_deltaK_tagh_fcal_phi[ip]     = new TH2F(Form("h_deltaK_tagh_fcal_phi_%d",ip), 
			Form("#DeltaK (FCAL Phi Sect. %d)",ip), 274, 0.5, 274.5, 2000, -8., 8.);
		h_deltaK_tagm_fcal_phi[ip]     = new TH2F(Form("h_deltaK_tagm_fcal_phi_%d",ip), 
			Form("#DeltaK (FCAL Phi Sect. %d)",ip), 102, 0.5, 102.5, 2000, -8., 8.);
		h_deltaK_tagh_cut_fcal_phi[ip] = new TH2F(Form("h_deltaK_tagh_cut_fcal_phi_%d",ip), 
			Form("#DeltaK (FCAL Phi Sect. %d)",ip), 274, 0.5, 274.5, 2000, -8., 8.);
		h_deltaK_tagm_cut_fcal_phi[ip] = new TH2F(Form("h_deltaK_tagm_cut_fcal_phi_%d",ip), 
			Form("#DeltaK (FCAL Phi Sect. %d)",ip), 102, 0.5, 102.5, 2000, -8., 8.);
		
		h_deltaK_tagh_ccal_phi[ip]     = new TH2F(Form("h_deltaK_tagh_ccal_phi_%d",ip), 
			Form("#DeltaK (CCAL Phi Sect. %d)",ip), 274, 0.5, 274.5, 2000, -8., 8.);
		h_deltaK_tagm_ccal_phi[ip]     = new TH2F(Form("h_deltaK_tagm_ccal_phi_%d",ip), 
			Form("#DeltaK (CCAL Phi Sect. %d)",ip), 102, 0.5, 102.5, 2000, -8., 8.);
		h_deltaK_tagh_cut_ccal_phi[ip] = new TH2F(Form("h_deltaK_tagh_cut_ccal_phi_%d",ip), 
			Form("#DeltaK (CCAL Phi Sect. %d)",ip), 274, 0.5, 274.5, 2000, -8., 8.);
		h_deltaK_tagm_cut_ccal_phi[ip] = new TH2F(Form("h_deltaK_tagm_cut_ccal_phi_%d",ip), 
			Form("#DeltaK (CCAL Phi Sect. %d)",ip), 102, 0.5, 102.5, 2000, -8., 8.);
		
		h_xy_fcal_phi[ip] = new TH2F(Form("h_xy_fcal_phi_%d",ip), 
			Form("FCAL Y vs. X (Phi Sect. %d)",ip), 500, -100., 100., 500, -100., 100.);
		h_xy_ccal_phi[ip] = new TH2F(Form("h_xy_ccal_phi_%d",ip), 
			Form("CCAL Y vs. X (Phi Sect. %d)",ip), 500,  -13.,  13., 500,  -13.,  13.);
	}
	
	for(int ip = 0; ip < 8; ip++) {
		
		h_deltaK_tagh_fcal_layer[ip] = new TH2F(Form("h_deltaK_tagh_fcal_layer_%d",ip+1), 
			Form("#DeltaK (FCAL Layer %d)",ip+1), 274, 0.5, 274.5, 2000, -8., 8.);
		h_deltaK_tagm_fcal_layer[ip] = new TH2F(Form("h_deltaK_tagm_fcal_layer_%d",ip+1), 
			Form("#DeltaK (FCAL Layer %d)",ip+1), 102, 0.5, 102.5, 2000, -8., 8.);
		
		h_deltaK_tagh_cut_fcal_layer[ip] = new TH2F(
			Form("h_deltaK_tagh_cut_fcal_layer_%d",ip+1), 
			Form("#DeltaK (FCAL Layer %d)",ip+1), 274, 0.5, 274.5, 2000, -8., 8.);
		h_deltaK_tagm_cut_fcal_layer[ip] = new TH2F(
			Form("h_deltaK_tagm_cut_fcal_layer_%d",ip+1), 
			Form("#DeltaK (FCAL Layer %d)",ip+1), 102, 0.5, 102.5, 2000, -8., 8.);
		
		h_xy_fcal_layer[ip] = new TH2F(Form("h_xy_fcal_layer_%d",ip+1), 
			Form("FCAL Y vs. X (Layer %d)",ip+1), 500, -100., 100., 500, -100., 100.);
	}
	
	for( int ip = 0; ip < 5; ip++ ) {
		
		h_deltaK_tagh_ccal_layer[ip] = new TH2F(Form("h_deltaK_tagh_ccal_layer_%d",ip+1), 
			Form("#DeltaK (CCAL Layer %d)",ip+1), 274, 0.5, 274.5, 2000, -8., 8.);
		h_deltaK_tagm_ccal_layer[ip] = new TH2F(Form("h_deltaK_tagm_ccal_layer_%d",ip+1), 
			Form("#DeltaK (CCAL Layer %d)",ip+1), 102, 0.5, 102.5, 2000, -8., 8.);
		
		h_deltaK_tagh_cut_ccal_layer[ip] = new TH2F(
			Form("h_deltaK_tagh_cut_ccal_layer_%d",ip+1), 
			Form("#DeltaK (CCAL Layer %d)",ip+1), 274, 0.5, 274.5, 2000, -8., 8.);
		h_deltaK_tagm_cut_ccal_layer[ip] = new TH2F(
			Form("h_deltaK_tagm_cut_ccal_layer_%d",ip+1), 
			Form("#DeltaK (CCAL Layer %d)",ip+1), 102, 0.5, 102.5, 2000, -8., 8.);
		
		h_xy_ccal_layer[ip] = new TH2F( Form("h_xy_ccal_layer_%d",ip+1), 
			Form("CCAL Y vs. X (Layer %d)",ip+1), 500, -13., 13., 500, -13., 13.);
	}
	
	//-------------------------------------------------------------------------------//
	
	h_deltaK_tagh_single = new TH2F( "deltaK_tagh_single", 
		"#DeltaK (only one CCAL Shower, FCAL Shower, and DBeamPhoton); TAGH Counter; [GeV]", 
		274, 0.5, 274.5, 2000, -8.0, 8.0 );
	h_deltaK_tagm_single = new TH2F( "deltaK_tagm_single", 
		"#DeltaK (only one CCAL Shower, FCAL Shower, and DBeamPhoton); TAGM Counter; [GeV]", 
		102, 0.5, 102.5, 2000, -8.0, 8.0 );
	
	h_deltaK_tagh_main       = new TH2F("deltaK_tagh_main", 
		"#DeltaK (Main RF Bunch); TAGH Counter; E_{Comp} - E_{#gamma} [GeV]", 
		274, 0.5, 274.5, 2000, -8.0, 8.0);
	h_deltaK_tagm_main       = new TH2F("deltaK_tagm_main", 
		"#DeltaK (Main RF Bunch); TAGM Counter; E_{Comp} - E_{#gamma} [GeV]", 
		102, 0.5, 102.5, 2000, -8.0, 8.0);
	h_deltaK_tagh_main_extra = new TH2F("deltaK_tagh_main_extra", 
		"#DeltaK (Main RF Bunches); TAGH Counter; E_{Comp} - E_{#gamma} [GeV]", 
		274, 0.5, 274.5, 2000, -8.0, 8.0);
	h_deltaK_tagm_main_extra = new TH2F("deltaK_tagm_main_extra", 
		"#DeltaK (Main RF Bunches); TAGM Counter; E_{Comp} - E_{#gamma} [GeV]", 
		102, 0.5, 102.5, 2000, -8.0, 8.0);
	
	h_deltaK_tagh_acc        = new TH2F("deltaK_tagh_acc", 
		"#DeltaK (Side RF Bunches); TAGH Counter; E_{Comp} - E_{#gamma} [GeV]", 
		274, 0.5, 274.5, 2000, -8.0, 8.0);
	h_deltaK_tagm_acc        = new TH2F("deltaK_tagm_acc", 
		"#DeltaK (Side RF Bunches); TAGM Counter; E_{Comp} - E_{#gamma} [GeV]", 
		102, 0.5, 102.5, 2000, -8.0, 8.0);
	h_deltaK_tagh_acc_single = new TH2F("deltaK_tagh_acc_single", 
		"#DeltaK (Side RF Bunch); TAGH Counter; E_{Comp} - E_{#gamma} [GeV]", 
		274, 0.5, 274.5, 2000, -8.0, 8.0);
	h_deltaK_tagm_acc_single = new TH2F("deltaK_tagm_acc_single", 
		"#DeltaK (Side RF Bunch); TAGM Counter; E_{Comp} - E_{#gamma} [GeV]", 
		102, 0.5, 102.5, 2000, -8.0, 8.0);
	h_deltaK_tagh_acc_extra  = new TH2F("deltaK_tagh_acc_extra", 
		"#DeltaK (Side RF Bunches); TAGH Counter; E_{Comp} - E_{#gamma} [GeV]", 
		274, 0.5, 274.5, 2000, -8.0, 8.0);
	h_deltaK_tagm_acc_extra  = new TH2F("deltaK_tagm_acc_extra", 
		"#DeltaK (Side RF Bunches); TAGM Counter; E_{Comp} - E_{#gamma} [GeV]", 
		102, 0.5, 102.5, 2000, -8.0, 8.0);
	
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
		
		deltaPhi_mu_p0  =  1.79825e+02;
		deltaPhi_mu_p1  = -1.10610e-02;
		deltaPhi_mu_p2  =  0.;
		deltaPhi_mu_p3  =  0.;
		//---------------------------// 
		deltaPhi_sig_p0 =  1.33099e+01;
		deltaPhi_sig_p1 = -2.14366e+00;
		deltaPhi_sig_p2 =  2.06214e-01;
		deltaPhi_sig_p3 = -6.71967e-03;
		
		deltaK_mu_p0    =  2.09726e-01;
		deltaK_mu_p1    = -4.71020e-02;
		deltaK_mu_p2    =  3.81198e-04;
		deltaK_mu_p3    = -4.82342e-05;
		//---------------------------// 
		deltaK_sig_p0   =  5.19636e-01;
		deltaK_sig_p1   = -3.35925e-02;
		deltaK_sig_p2   =  1.03144e-02;
		deltaK_sig_p3   = -3.63616e-04;
		
	} else if(group < 8) { // Be FIELD ON (stand-in values, need to update)
		
		deltaE_mu_p0    = -1.15974e-01;
		deltaE_mu_p1    =  3.09559e-02;
		deltaE_mu_p2    = -2.18315e-03;
		deltaE_mu_p3    =  0.;
		//---------------------------//
		deltaE_sig_p0   =  1.59557e-02;
		deltaE_sig_p1   =  3.36354e-02;
		deltaE_sig_p2   =  0.;
		
		deltaPhi_mu_p0  =  1.80032e+02;
		deltaPhi_mu_p1  = -9.61299e-03;
		deltaPhi_mu_p2  =  0.;
		deltaPhi_mu_p3  =  0.;
		//---------------------------// 
		deltaPhi_sig_p0 =  1.26943e+01;
		deltaPhi_sig_p1 = -1.95900e+00;
		deltaPhi_sig_p2 =  1.87668e-01;
		deltaPhi_sig_p3 = -6.08281e-03;
		
		deltaK_mu_p0    =  3.52591e-01;
		deltaK_mu_p1    = -9.12544e-02;
		deltaK_mu_p2    =  5.17767e-03;
		deltaK_mu_p3    = -2.20872e-04;
		//---------------------------// 
		deltaK_sig_p0   =  5.97017e-01;
		deltaK_sig_p1   = -5.77323e-02;
		deltaK_sig_p2   =  1.29216e-02;
		deltaK_sig_p3   = -4.54356e-04;
		
	} else if(group < 9) { // He FIELD OFF (stand-in values, need to update)
		
		deltaE_mu_p0    = -1.15974e-01;
		deltaE_mu_p1    =  3.09559e-02;
		deltaE_mu_p2    = -2.18315e-03;
		deltaE_mu_p3    =  0.;
		//---------------------------//
		deltaE_sig_p0   =  1.59557e-02;
		deltaE_sig_p1   =  3.36354e-02;
		deltaE_sig_p2   =  0.;
		
		deltaPhi_mu_p0  =  1.80032e+02;
		deltaPhi_mu_p1  = -9.61299e-03;
		deltaPhi_mu_p2  =  0.;
		deltaPhi_mu_p3  =  0.;
		//---------------------------// 
		deltaPhi_sig_p0 =  1.26943e+01;
		deltaPhi_sig_p1 = -1.95900e+00;
		deltaPhi_sig_p2 =  1.87668e-01;
		deltaPhi_sig_p3 = -6.08281e-03;
		
		deltaK_mu_p0    =  3.52591e-01;
		deltaK_mu_p1    = -9.12544e-02;
		deltaK_mu_p2    =  5.17767e-03;
		deltaK_mu_p3    = -2.20872e-04;
		//---------------------------// 
		deltaK_sig_p0   =  5.97017e-01;
		deltaK_sig_p1   = -5.77323e-02;
		deltaK_sig_p2   =  1.29216e-02;
		deltaK_sig_p3   = -4.54356e-04;
		
	} else if(group < 12) { // He FIELD ON (stand-in values, need to update)
		
		deltaE_mu_p0    = -1.15974e-01;
		deltaE_mu_p1    =  3.09559e-02;
		deltaE_mu_p2    = -2.18315e-03;
		deltaE_mu_p3    =  0.;
		//---------------------------//
		deltaE_sig_p0   =  1.59557e-02;
		deltaE_sig_p1   =  3.36354e-02;
		deltaE_sig_p2   =  0.;
		
		deltaPhi_mu_p0  =  1.80032e+02;
		deltaPhi_mu_p1  = -9.61299e-03;
		deltaPhi_mu_p2  =  0.;
		deltaPhi_mu_p3  =  0.;
		//---------------------------// 
		deltaPhi_sig_p0 =  1.26943e+01;
		deltaPhi_sig_p1 = -1.95900e+00;
		deltaPhi_sig_p2 =  1.87668e-01;
		deltaPhi_sig_p3 = -6.08281e-03;
		
		deltaK_mu_p0    =  3.52591e-01;
		deltaK_mu_p1    = -9.12544e-02;
		deltaK_mu_p2    =  5.17767e-03;
		deltaK_mu_p3    = -2.20872e-04;
		//---------------------------// 
		deltaK_sig_p0   =  5.97017e-01;
		deltaK_sig_p1   = -5.77323e-02;
		deltaK_sig_p2   =  1.29216e-02;
		deltaK_sig_p3   = -4.54356e-04;
		
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
	}
	
	f_deltaE_mu->SetParameters(deltaE_mu_p0, deltaE_mu_p1, deltaE_mu_p2, deltaE_mu_p3);
	f_deltaE_sig->SetParameters(deltaE_sig_p0, deltaE_sig_p1, deltaE_sig_p2);
	
	f_deltaPhi_mu->SetParameters(deltaPhi_mu_p0,  deltaPhi_mu_p1,  
		deltaPhi_mu_p2,  deltaPhi_mu_p3);
	f_deltaPhi_sig->SetParameters(deltaPhi_sig_p0, deltaPhi_sig_p1, 
		deltaPhi_sig_p2, deltaPhi_sig_p3);
	
	f_deltaK_mu->SetParameters(deltaK_mu_p0,  deltaK_mu_p1,  deltaK_mu_p2,  deltaK_mu_p3);
	f_deltaK_sig->SetParameters(deltaK_sig_p0, deltaK_sig_p1, deltaK_sig_p2, deltaK_sig_p3);
	
	return;
}
