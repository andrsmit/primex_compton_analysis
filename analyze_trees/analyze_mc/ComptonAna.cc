#include "ComptonAna.h"

// Default Constructor:

ComptonAna::ComptonAna() {
	
	m_SHIFT_DISTRIBUTIONS = 0;
	m_SMEAR_DISTRIBUTIONS = 0; // note: if m_SHIFT_DISTRIBUTIONS is 0, then m_SMEAR_DISTRIBUTIONS is not checked
	
	// set default run number to 61321:
	
	m_runNumber = 61321;
	
	m_random = new TRandom3(0);
	
	// set defaults for cut values:
	
	m_cut_fcalE    = 0.35;
	m_cut_ccalE    = 3.0;
	m_cut_fcalrfdt = 2.004;
	m_cut_ccalrfdt = 2.004;
	m_cut_beamrfdt = 2.004;
	m_cut_deltaE   = 5.0;
	m_cut_deltaK   = 5.0;
	m_cut_deltaPhi = 5.0;
	
	// Defaults for Geometry from CCDB:
	
	m_phase_val      = 1;
	m_target_length  = 1.7755;
	m_target_density = 1.848;
	m_target_atten   = 0.01172;
	
	m_fcal_face.SetXYZ(0.529, -0.002, 624.906);
	m_fcal_correction.SetXYZ(0., 0., 0);
	
	m_ccal_face.SetXYZ(-0.0225, 0.0073, 1279.376);
	m_ccal_correction.SetXYZ(0., 0., 0.);
	
	m_vertex.SetXYZ(0.1914, -0.0769, 65.);
	
	// Set event number to zero on initialization:
	
	m_event = 0;
}

void ComptonAna::runAnalysis(TString infname) {
	
	m_infile = new TFile(infname.Data(), "READ");
	int n_events_total = loadTree();
	
	while(m_event < n_events_total) {
		readEvent();
		comptonAnalysis();
		m_event++;
	}
	
	m_infile->Close();
	
	return;
}

void ComptonAna::comptonAnalysis() {
	
	if(m_mc_e[0] <= 0.) h_reaction_weight->Fill(m_reaction_weight);
	else h_reaction_weight_double->Fill(m_reaction_weight);
	
	if(m_nmc>0) {
		h_vertex->Fill(m_mc_z[0]);
		if(acceptRejectEvent()) return;
		h_vertex_accepted->Fill(m_mc_z[0]);
	}
	
	if(m_nfcal > MAX_FCAL || m_nccal > MAX_CCAL || m_nbeam > MAX_BEAM) {
		cout << "Skipping event " << m_event << endl;
		return;
	}
	
	int is_triplet = false;
	if(m_nmc=3) {
		if(m_mc_pdgtype[2]==11.0) is_triplet = true;
	}
	//if(!is_triplet) return;
	
	if(m_event > 5.e5) return;
	
	// get the angle of the scattered photon:
	
	double loc_thrown_photon_angle = m_mc_theta[1];
	h_thrown_photon_angle->Fill(loc_thrown_photon_angle);
	
	// calculate invariant mass of thrown pair:
	int part1 = 0;
	int part2 = 1;
	//if(is_triplet && m_mc_e[2]>m_mc_e[1]) part2 = 2;
	
	double  e_thrown = m_mc_e[part1] + m_mc_e[part2];
	
	double px_thrown = m_mc_p[part1]*sin(m_mc_theta[part1]*TMath::Pi()/180.0)*cos(m_mc_phi[part1]*TMath::Pi()/180.0)
		+ m_mc_p[part2]*sin(m_mc_theta[part2]*TMath::Pi()/180.0)*cos(m_mc_phi[part2]*TMath::Pi()/180.0);
	
	double py_thrown = m_mc_p[part1]*sin(m_mc_theta[part1]*TMath::Pi()/180.0)*sin(m_mc_phi[part1]*TMath::Pi()/180.0)
		+ m_mc_p[part2]*sin(m_mc_theta[part2]*TMath::Pi()/180.0)*sin(m_mc_phi[part2]*TMath::Pi()/180.0);
	
	double pz_thrown = m_mc_p[part1]*cos(m_mc_theta[part1]*TMath::Pi()/180.0) 
		+ m_mc_p[part2]*cos(m_mc_theta[part2]*TMath::Pi()/180.0);
	
	double invmass_thrown = sqrt(pow(e_thrown,2.0) - pow(px_thrown,2.0) - pow(py_thrown,2.0) - pow(pz_thrown,2.0));
	h_thrown_mass->Fill(invmass_thrown);
	
	//
	// Make a list of good FCAL showers to use in analysis:
	//
	vector<int> locGoodFCALShowers;
	locGoodFCALShowers.clear();
	int locNFCALShowers = 0, locNGoodFCALShowers = 0;
	for(int ishow=0; ishow<m_nfcal; ishow++) {
		
		TVector3 loc_pos = getFCALPosition(ishow);
		double loc_t = m_fcal_t[ishow] - (loc_pos.Mag()/m_c) - m_rfTime;
		
		h_fcal_rf_dt->Fill(loc_t);
		//locGoodFCALShowers.push_back(ishow);
		
		int loc_fid_cut = fcal_fiducial_cut(loc_pos, 1.0);
		if(fabs(loc_t) < m_cut_fcalrfdt) {
			locNFCALShowers++;
			if(m_fcal_e[ishow] > m_cut_fcalE) {
				locNGoodFCALShowers++;
				if(!loc_fid_cut) {
					locGoodFCALShowers.push_back(ishow);
				}
			}
		}
	}
	
	//
	// Make a list of good CCAL showers to use in analysis:
	//
	vector<int> locGoodCCALShowers;
	locGoodCCALShowers.clear();
	int locNCCALShowers = 0, locNGoodCCALShowers = 0;
	for(int ishow=0; ishow<m_nccal; ishow++) {
		
		TVector3 loc_pos = getCCALPosition(ishow);
		double loc_t = m_ccal_t[ishow] - (loc_pos.Mag()/m_c) - m_rfTime;
		
		h_ccal_rf_dt->Fill(loc_t);
		//locGoodCCALShowers.push_back(ishow);
		
		int loc_fid_cut = ccal_fiducial_cut(loc_pos);
		if(fabs(loc_t) < m_cut_ccalrfdt) {
			locNCCALShowers++;
			if(m_ccal_e[ishow] > m_cut_ccalE) {
				locNGoodCCALShowers++;
				if(!loc_fid_cut) {
					locGoodCCALShowers.push_back(ishow);
				}
			}
		}
	}
	
	//---------------------------------------------------------------------------//
	// Make a list of prompt and selected-sideband beam photons to use in analysis:
	
	vector<int> locGoodBeamPhotons;
	vector<double> locGoodBeamPhotons_weight;
	locGoodBeamPhotons.clear();
	locGoodBeamPhotons_weight.clear();
	
	for(int igam=0; igam<m_nbeam; igam++) {
		
		double loc_dt     = m_beam_t[igam] - m_rfTime;
		double loc_weight = 0.0;
		
		h_beam_rf_dt->Fill(loc_dt);
		
		double loc_beam_cut = m_beam_bunches_main*m_cut_beamrfdt;
		
		if(fabs(loc_dt) < loc_beam_cut) loc_weight = 1.0;
		else if(
			((m_beam_bunches_main+5.5)*4.008 < fabs(loc_dt)) && 
			(fabs(loc_dt) < (m_beam_bunches_main+5.5+m_beam_bunches_acc)*4.008)
		) loc_weight = -1.0/(2.0*m_beam_bunches_acc);
		else continue;
		
		//if(loc_weight < 0.0) loc_weight *= m_acc_scale_factor[igam];
		if(m_beam_e[igam] > 6.0) {
			locGoodBeamPhotons.push_back(igam);
			locGoodBeamPhotons_weight.push_back(loc_weight);
		}
	}
	
	for(auto &ifcal : locGoodFCALShowers) {
		
		double e1 = m_fcal_e[ifcal];
		
		TVector3 pos1 = getFCALPosition(ifcal);
		double t1 = m_fcal_t[ifcal] - (pos1.Mag()/m_c) - m_rfTime;
		
		double theta1 = pos1.Theta();
		
		//--------------------------------------------------------------------------------//
		// check for matches with TOF:
		
		double tof_dx, tof_dy, tof_dt;
		check_TOF_match(pos1, tof_dx, tof_dy, tof_dt, 1.0);
		double tof_dr = sqrt(pow(tof_dx,2.0) + pow(tof_dy,2.0));
		
		int tof_match = 0;
		if(tof_dr < 8.0) tof_match = 1;
		
		//--------------------------------------------------------------------------------//
		
		int   fid_cut_1 = fcal_fiducial_cut(pos1, 1.0);
		int min_e_cut_1 = e1 > 0.35 ? 1 : 0;
		int  rf_t_cut_1 = fabs(t1) < 2.004 ? 1 : 0;
		
		for(auto &iccal : locGoodCCALShowers) {
			
			double e2 = m_ccal_e[iccal];
			
			TVector3 pos2 = getCCALPosition(iccal);
			double t2 = m_ccal_t[iccal] - (pos2.Mag()/m_c) - m_rfTime;
			
			double theta2 = pos2.Theta();
			
			// calculate deltaPhi:
			
			double deltaPhi = fabs(pos2.Phi() - pos1.Phi()) * (180./TMath::Pi());
			
			// calculate invariant mass (assuming massless particles):
			
			double cos12   = (pos1.X()*pos2.X() + pos1.Y()*pos2.Y() + pos1.Z()*pos2.Z()) 
				/ (pos1.Mag()*pos2.Mag());
			double invmass = sqrt(2.0*e1*e2*(1.-cos12));
			
			double opangle = acos(cos12) * (180./TMath::Pi());
			
			int   fid_cut_2 = ccal_fiducial_cut(pos2);
			int min_e_cut_2 = e2 > 3.0 ? 1 : 0;
			int  rf_t_cut_2 = fabs(t2) < 2.004 ? 1 : 0;
			
			// loop over beam photons:
			
			for(int igam=0; igam<(int)locGoodBeamPhotons.size(); igam++) {
				
				int loc_beam_index  = locGoodBeamPhotons[igam];
				double eb           = m_beam_e[loc_beam_index];
				double fill_weight  = locGoodBeamPhotons_weight[igam];
				
				int loc_tag_sys     = m_tag_sys[loc_beam_index];
				int loc_tag_counter = m_tag_counter[loc_beam_index];
				
				double deltaE = (e1+e2) - (eb+m_e);
				
				double deltaK = m_e*(sin(theta1) / (2.*pow(sin(theta1/2.),2.)*tan(theta2))) - m_e;
				deltaK       -= eb;
				
				//-------------------------------------------------------------//
				// Cuts:
				
				double deltaE_smeared   = smear_deltaE(deltaE, eb);
				double deltaPhi_smeared = smear_deltaPhi(deltaPhi, eb);
				double deltaK_smeared   = smear_deltaK(deltaK, eb);
				
				double deltaE_smeared_two   = smear_deltaE_two(deltaE, eb);
				double deltaPhi_smeared_two = smear_deltaPhi_two(deltaPhi, eb);
				double deltaK_smeared_two   = smear_deltaK_two(deltaK, eb);
				
				int   e_cut = cut_deltaE(  deltaE_smeared,   eb, m_cut_deltaE,   1.e2);
				int phi_cut = cut_deltaPhi(deltaPhi_smeared, eb, m_cut_deltaPhi, m_cut_deltaPhi);
				int   k_cut = cut_deltaK(  deltaK_smeared,   eb, m_cut_deltaK,   m_cut_deltaK);
				
				int   e_cut_two = cut_deltaE_two(  deltaE_smeared_two,   eb, m_cut_deltaE,   1.e2);
				int phi_cut_two = cut_deltaPhi_two(deltaPhi_smeared_two, eb, m_cut_deltaPhi, m_cut_deltaPhi);
				int   k_cut_two = cut_deltaK_two(  deltaK_smeared_two,   eb, m_cut_deltaK,   m_cut_deltaK);
				
				//-------------------------------------------------------------//
				//
				// We want to look at how the following set of cuts change the data:
				//
				//    1. Multiplicity cut (3 options)
				//        a. No extra showers
				//        b. Allow extra showers if they are below some energy threshold
				//        c. Allow any ammount of extra showers
				//    2. TOF Veto
				//        a. No match betweeen FCAL and TOF
				//        b. Maybe a match, maybe not (no information)
				//    3. FCAL fiducial cut
				//        a. One layer cut
				//        b. Two layer cut
				//
				
				int cut_vals[m_n_cuts];
				for(int icut=0; icut<m_n_cuts; icut++) cut_vals[icut] = 0;
				
				int loc_fid_cut = fcal_fiducial_cut(pos1, 2.0);
				
				cut_vals[0] = 1;
				if(!tof_match) {
					cut_vals[1] = 1;
					if(!loc_fid_cut) {
						cut_vals[2] = 1;
					}
				}
				if(!loc_fid_cut) cut_vals[3] = 1;
				
				if(locNGoodFCALShowers==1 && locNGoodCCALShowers==1) {
					cut_vals[4] = 1;
					if(!tof_match) {
						cut_vals[5] = 1;
						if(!loc_fid_cut) {
							cut_vals[6] = 1;
						}
					}
					if(!loc_fid_cut) cut_vals[7] = 1;
				}
				
				if(locNFCALShowers==1 && locNCCALShowers==1) {
					cut_vals[8] = 1;
					if(!tof_match) {
						cut_vals[9] = 1;
						if(!loc_fid_cut) {
							cut_vals[10] = 1;
						}
					}
					if(!loc_fid_cut) cut_vals[11] = 1;
				}
				
				for(int icut = 0; icut < m_n_cuts; icut++) {
					
					int loc_e_cut   = e_cut;
					int loc_phi_cut = phi_cut;
					int loc_k_cut   = k_cut;
					
					double loc_deltaE   = deltaE_smeared;
					double loc_deltaPhi = deltaPhi_smeared;
					double loc_deltaK   = deltaK_smeared;
					
					if(icut%4==2 || icut%4==3) {
						loc_e_cut = e_cut_two;
						loc_phi_cut = phi_cut_two;
						loc_k_cut = k_cut_two;
						
						loc_deltaE   = deltaE_smeared_two;
						loc_deltaPhi = deltaPhi_smeared_two;
						loc_deltaK   = deltaK_smeared_two;
					}
					
					if(loc_tag_sys==0) {
						
						h_opangle_tagh[icut]->Fill(loc_tag_counter, opangle, fill_weight);
						if(loc_e_cut) {
							h_opangle_tagh_ecut[icut]->Fill(loc_tag_counter, opangle, fill_weight);
							if(loc_k_cut) {
								h_opangle_tagh_ekcut[icut]->Fill(loc_tag_counter, opangle, fill_weight);
							}
						}
						
						if(cut_vals[icut]) {
							h_deltaE_tagh[icut]->Fill(loc_tag_counter, loc_deltaE, fill_weight);
							if(loc_e_cut) {
								h_deltaPhi_tagh[icut]->Fill(loc_tag_counter, loc_deltaPhi, 
									fill_weight);
								if(loc_phi_cut) {
									h_deltaK_tagh[icut]->Fill(loc_tag_counter, loc_deltaK, 
										fill_weight);
								}
							}
						}
					} else if(loc_tag_sys==1) {
						
						h_opangle_tagm[icut]->Fill(loc_tag_counter, opangle, fill_weight);
						if(loc_e_cut) {
							h_opangle_tagm_ecut[icut]->Fill(loc_tag_counter, opangle, fill_weight);
							if(loc_k_cut) {
								h_opangle_tagm_ekcut[icut]->Fill(loc_tag_counter, opangle, fill_weight);
							}
						}
						
						if(cut_vals[icut]) {
							h_deltaE_tagm[icut]->Fill(loc_tag_counter, loc_deltaE, fill_weight);
							if(loc_e_cut) {
								h_deltaPhi_tagm[icut]->Fill(loc_tag_counter, loc_deltaPhi,
									fill_weight);
								if(loc_phi_cut) {
									h_deltaK_tagm[icut]->Fill(loc_tag_counter, loc_deltaK, 
										fill_weight);
								}
							}
						}
					}
					
					double deltaKCCAL = e2 - (1. / (1./eb + (1.-cos(pos2.Theta()))/m_e));
					double deltaKFCAL = e1 - (1. / (1./eb + (1.-cos(pos1.Theta()))/m_e));
					
					if(cut_vals[icut]) {
						h_elas_vs_deltaE[icut]->Fill(loc_deltaE, 
							(loc_deltaE-loc_deltaK), fill_weight);
						
						//h_mgg_vs_deltaE[icut]->Fill(deltaE_smeared, invmass, fill_weight);
						//h_mgg_vs_deltaK[icut]->Fill(deltaK_smeared, invmass, fill_weight);
						
						if(eb > 9.0) {
							h_deltaK_vs_deltaE[icut]->Fill(loc_deltaE, loc_deltaK, fill_weight);
							
							h_deltaCCAL_vs_deltaE[icut]->Fill(loc_deltaE, deltaKCCAL, fill_weight);
							h_deltaCCAL_vs_deltaK[icut]->Fill(loc_deltaK, deltaKCCAL, fill_weight);
							h_deltaFCAL_vs_deltaE[icut]->Fill(loc_deltaE, deltaKFCAL, fill_weight);
							h_deltaFCAL_vs_deltaK[icut]->Fill(loc_deltaK, deltaKFCAL, fill_weight);
							
							h_deltaFCAL_vs_deltaCCAL[icut]->Fill(deltaKCCAL, deltaKFCAL, fill_weight);
						}
						
						if(loc_e_cut && loc_k_cut && loc_phi_cut) {
							h_ccal_xy[icut]->Fill(pos2.X(), pos2.Y(), fill_weight);
							h_fcal_xy[icut]->Fill(pos1.X(), pos1.Y(), fill_weight);
						}
					}
				}
				
				
				h_deltaE_vs_thrown_mass->Fill(invmass_thrown, deltaE_smeared, fill_weight);
				h_deltaE_vs_rec_mass->Fill(invmass, deltaE_smeared, fill_weight);
				h_deltaK_vs_thrown_mass->Fill(invmass_thrown, deltaK_smeared, fill_weight);
				h_deltaK_vs_rec_mass->Fill(invmass, deltaK_smeared, fill_weight);
				if(invmass_thrown<0.03) {
					h_elas_low_mass->Fill(deltaE_smeared, deltaK_smeared, fill_weight);
				} else {
					h_elas_high_mass->Fill(deltaE_smeared, deltaK_smeared, fill_weight);
				}
				h_rec_vs_thrown_mass->Fill(invmass_thrown, invmass, fill_weight);
				
				if(cut_vals[4]) {
					
					if(e_cut && phi_cut && k_cut) {
						h_thrown_photon_angle_accepted->Fill(loc_thrown_photon_angle, fill_weight);
					}
					
					if(e_cut) {
						double sumPhi = (pos2.Phi()+ pos1.Phi()) * (180./TMath::Pi());
						h_sumPhi_vs_deltaPhi->Fill(deltaPhi/2.0, sumPhi/2.0, fill_weight);
					}
					
					h_ccal_nblocks->Fill(m_ccal_nblocks[iccal], fill_weight);
					h_fcal_nblocks->Fill(m_fcal_nblocks[ifcal], fill_weight);
					if(deltaE > 0.65) {
						h_ccal_xy_highdeltaE->Fill(pos2.X(), pos2.Y(), fill_weight);
						h_fcal_xy_highdeltaE->Fill(pos1.X(), pos1.Y(), fill_weight);
						h_ccal_nblocks_cut->Fill(m_ccal_nblocks[iccal], fill_weight);
						h_fcal_nblocks_cut->Fill(m_fcal_nblocks[ifcal], fill_weight);
					}
					if(deltaE < -3.5 && eb>8.5) {
						h_ccal_xy_lowdeltaE->Fill(pos2.X(), pos2.Y(), fill_weight);
						h_fcal_xy_lowdeltaE->Fill(pos1.X(), pos1.Y(), fill_weight);
					}
				}
				
				if(is_triplet) {
					h_deltaE_triplet[0]->Fill(eb, deltaE_smeared, fill_weight);
					h_deltaK_triplet[0]->Fill(eb, deltaK_smeared, fill_weight);
					h_deltaPhi_triplet[0]->Fill(eb, deltaPhi_smeared, fill_weight);
				} else {
					h_deltaE_pair[0]->Fill(eb, deltaE_smeared, fill_weight);
					h_deltaK_pair[0]->Fill(eb, deltaK_smeared, fill_weight);
					h_deltaPhi_pair[0]->Fill(eb, deltaPhi_smeared, fill_weight);
				}
				if(!fid_cut_1 && !fid_cut_2) {
					if(is_triplet) {
						h_deltaE_triplet[1]->Fill(eb, deltaE_smeared, fill_weight);
						h_deltaK_triplet[1]->Fill(eb, deltaK_smeared, fill_weight);
						h_deltaPhi_triplet[1]->Fill(eb, deltaPhi_smeared, fill_weight);
					} else {
						h_deltaE_pair[1]->Fill(eb, deltaE_smeared, fill_weight);
						h_deltaK_pair[1]->Fill(eb, deltaK_smeared, fill_weight);
						h_deltaPhi_pair[1]->Fill(eb, deltaPhi_smeared, fill_weight);
					}
					if(min_e_cut_1 && min_e_cut_2) {
						if(is_triplet) {
							h_deltaE_triplet[2]->Fill(eb, deltaE_smeared, fill_weight);
							h_deltaK_triplet[2]->Fill(eb, deltaK_smeared, fill_weight);
							h_deltaPhi_triplet[2]->Fill(eb, deltaPhi_smeared, fill_weight);
						} else {
							h_deltaE_pair[2]->Fill(eb, deltaE_smeared, fill_weight);
							h_deltaK_pair[2]->Fill(eb, deltaK_smeared, fill_weight);
							h_deltaPhi_pair[2]->Fill(eb, deltaPhi_smeared, fill_weight);
						}
						if(rf_t_cut_1 && rf_t_cut_2) {
							if(is_triplet) {
								h_deltaE_triplet[3]->Fill(eb, deltaE_smeared, fill_weight);
								h_deltaK_triplet[3]->Fill(eb, deltaK_smeared, fill_weight);
								h_deltaPhi_triplet[3]->Fill(eb, deltaPhi_smeared, fill_weight);
							} else {
								h_deltaE_pair[3]->Fill(eb, deltaE_smeared, fill_weight);
								h_deltaK_pair[3]->Fill(eb, deltaK_smeared, fill_weight);
								h_deltaPhi_pair[3]->Fill(eb, deltaPhi_smeared, fill_weight);
							}
							if(e_cut) {
								if(is_triplet) {
									h_deltaE_triplet[4]->Fill(eb, deltaE_smeared, fill_weight);
									h_deltaK_triplet[4]->Fill(eb, deltaK_smeared, fill_weight);
									h_deltaPhi_triplet[4]->Fill(eb, deltaPhi_smeared, fill_weight);
								} else {
									h_deltaE_pair[4]->Fill(eb, deltaE_smeared, fill_weight);
									h_deltaK_pair[4]->Fill(eb, deltaK_smeared, fill_weight);
									h_deltaPhi_pair[4]->Fill(eb, deltaPhi_smeared, fill_weight);
								}
								if(phi_cut) {
									if(is_triplet) {
										h_deltaE_triplet[5]->Fill(eb, deltaE_smeared, fill_weight);
										h_deltaK_triplet[5]->Fill(eb, deltaK_smeared, fill_weight);
										h_deltaPhi_triplet[5]->Fill(eb, deltaPhi_smeared, fill_weight);
									} else {
										h_deltaE_pair[5]->Fill(eb, deltaE_smeared, fill_weight);
										h_deltaK_pair[5]->Fill(eb, deltaK_smeared, fill_weight);
										h_deltaPhi_pair[5]->Fill(eb, deltaPhi_smeared, fill_weight);
									}
								} // deltaPhi cut
							} // deltaE cut
						} // rf timing cut
					} // min energy cut
				} // fiducial cut
				
			} // end loop over beam photons
		} // end loop over CCAL shower
	} // end loop over FCAL showers
	
	return;
}

int ComptonAna::acceptRejectEvent() {
	//
	// Accept-reject filter to create a realistic z-vertex distribution for simulated events
	//
	if(m_nmc==0) return 0;
	
	int reject = 0;
	
	double vertex_z = m_mc_z[0];
	
	// shift coordinate system so that upstream entrance of target is at z=0:
	
	double loc_z = vertex_z - m_vertex.Z() + (m_target_length/2.0);
	
	// use attenuation length from XCOM database to calculate probability of photon 
	// absorption.
	
	double vertex_weight;
	if(loc_z<0.) {
		vertex_weight = 1.0;
	} else if(loc_z>m_target_length) {
		vertex_weight = TMath::Exp(-m_target_atten * m_target_density * m_target_length);
	} else {
		vertex_weight = TMath::Exp(-m_target_atten * m_target_density * loc_z);
	}
	
	if(vertex_weight < m_random->Uniform()) {
		reject = 1;
	}
	
	return reject;
}

TVector3 ComptonAna::getFCALPosition(int index) {
	
	TVector3 pos(m_fcal_x[index], m_fcal_y[index], m_fcal_z[index]);
	pos = pos + m_fcal_correction - m_vertex;
	
	return pos;
}
TVector3 ComptonAna::getCCALPosition(int index) {
	
	TVector3 pos(m_ccal_x[index], m_ccal_y[index], m_ccal_z[index]);
	pos = pos + m_ccal_correction - m_vertex;
	
	return pos;
}

int ComptonAna::fcal_fiducial_cut(TVector3 pos, double cut_layer) {
	
	int fid_cut = 0;
	
	double fcal_inner_layer_cut = (1.5 + cut_layer) * m_fcal_block_size;
	
	double fcal_face_x = m_vertex.X() + (pos.X() * (m_fcal_face.Z() - m_vertex.Z())/pos.Z()) - m_fcal_face.X();
	double fcal_face_y = m_vertex.Y() + (pos.Y() * (m_fcal_face.Z() - m_vertex.Z())/pos.Z()) - m_fcal_face.Y();
	
	if((fabs(fcal_face_x) < fcal_inner_layer_cut) && (fabs(fcal_face_y) < fcal_inner_layer_cut)) fid_cut = 1;
	
	// only apply the next fiducial cut for runs from phase-I:
	
	if(m_phase_val < 2) {
		if((-32.<fcal_face_y && fcal_face_y<-20.) && (-8.<fcal_face_x && fcal_face_x<4.)) {
			fid_cut = 1;
		}
	}
	
	return fid_cut;
}

int ComptonAna::ccal_fiducial_cut(TVector3 pos, double cut_layer) {
	
	int fid_cut = 0;
	
	double ccal_inner_layer_cut = (1.0 + cut_layer) * m_ccal_block_size;
	
	double ccal_face_x = m_vertex.X() + (pos.X() * (m_ccal_face.Z() - m_vertex.Z())/pos.Z()) - m_ccal_face.X();
	double ccal_face_y = m_vertex.Y() + (pos.Y() * (m_ccal_face.Z() - m_vertex.Z())/pos.Z()) - m_ccal_face.Y();
	
	if((fabs(ccal_face_x) < ccal_inner_layer_cut) && (fabs(ccal_face_y) < ccal_inner_layer_cut)) fid_cut = 1;
	
	if(ccal_face_x<-8.36 || ccal_face_x>10.45 
		|| ccal_face_y<-10.45 || ccal_face_y>10.45) fid_cut = 1;
	
	return fid_cut;
}

void ComptonAna::check_TOF_match(TVector3 pos1, double &dx_min, double &dy_min, double &dt_min, double rf_time_cut) {
	
	dx_min = 1000.;
	dy_min = 1000.;
	dt_min = 1000.;
	
	for(int itof=0; itof<m_ntof; itof++) {
		
		double xt = m_tof_x[itof] - m_vertex.X();
		double yt = m_tof_y[itof] - m_vertex.Y();
		double zt = m_tof_z[itof] - m_vertex.Z();
		double rt = sqrt(xt*xt + yt*yt + zt*zt);
		double tt = m_tof_t[itof] - (rt/m_c);
		double dt = tt - m_rfTime;
		xt *= pos1.Z() / zt;
		yt *= pos1.Z() / zt;
		double dx = pos1.X() - xt;
		double dy = pos1.Y() - yt;
		
		if(fabs(dt) < rf_time_cut) {
			if((dx*dx + dy*dy) < (dx_min*dx_min + dy_min*dy_min)) {
				dx_min = dx;
				dy_min = dy;
				dt_min = dt;
			}
		}
	}
	
	return;
}

int ComptonAna::cut_deltaE(double deltaE, double eb, double n_sigma_left, double n_sigma_right) {
	
	int cut_val = 0;
	
	double loc_mu = 0.;
	for(int ipar=0; ipar<4; ipar++) loc_mu += (m_deltaE_mu_pars[0][ipar]*pow(eb,(double)ipar));
	
	double loc_sigma = sqrt(pow(m_deltaE_sigma_pars[0][0],2.0) + pow(m_deltaE_sigma_pars[0][1]/sqrt(eb),2.0) 
		+ pow(m_deltaE_sigma_pars[0][2]/eb,2.0));
	loc_sigma *= eb;
	
	double loc_diff = deltaE - loc_mu;
	if((-1.*n_sigma_left*loc_sigma < loc_diff) && (loc_diff < n_sigma_right*loc_sigma)) return 1;
	else return 0;
}

int ComptonAna::cut_deltaPhi(double deltaPhi, double eb, double n_sigma_left, double n_sigma_right) {
	
	int cut_val = 0;
	
	double loc_mu = 0.;
	for(int ipar=0; ipar<4; ipar++) loc_mu += (m_deltaPhi_mu_pars[0][ipar]*pow(eb,(double)ipar));
	
	double loc_sigma = 0.;
	for(int ipar=0; ipar<4; ipar++) loc_sigma += (m_deltaPhi_sigma_pars[0][ipar]*pow(eb,(double)ipar));
	
	double loc_diff = deltaPhi - loc_mu;
	if((-1.*n_sigma_left*loc_sigma < loc_diff) && (loc_diff < n_sigma_right*loc_sigma)) return 1;
	else return 0;
}

int ComptonAna::cut_deltaK(double deltaK, double eb, double n_sigma_left, double n_sigma_right) {
	
	int cut_val = 0;
	
	double loc_mu = 0.;
	for(int ipar=0; ipar<4; ipar++) loc_mu += (m_deltaK_mu_pars[0][ipar]*pow(eb,(double)ipar));
	
	double loc_sigma = 0.;
	for(int ipar=0; ipar<4; ipar++) loc_sigma += (m_deltaK_sigma_pars[0][ipar]*pow(eb,(double)ipar));
	
	double loc_diff = deltaK - loc_mu;
	if((-1.*n_sigma_left*loc_sigma < loc_diff) && (loc_diff < n_sigma_right*loc_sigma)) return 1;
	else return 0;
}

int ComptonAna::cut_deltaE_two(double deltaE, double eb, double n_sigma_left, double n_sigma_right) {
	
	int cut_val = 0;
	
	double loc_mu = 0.;
	for(int ipar=0; ipar<4; ipar++) loc_mu += (m_deltaE_mu_pars_two[0][ipar]*pow(eb,(double)ipar));
	
	double loc_sigma = sqrt(pow(m_deltaE_sigma_pars_two[0][0],2.0) + pow(m_deltaE_sigma_pars_two[0][1]/sqrt(eb),2.0) 
		+ pow(m_deltaE_sigma_pars_two[0][2]/eb,2.0));
	loc_sigma *= eb;
	
	double loc_diff = deltaE - loc_mu;
	if((-1.*n_sigma_left*loc_sigma < loc_diff) && (loc_diff < n_sigma_right*loc_sigma)) return 1;
	else return 0;
}

int ComptonAna::cut_deltaPhi_two(double deltaPhi, double eb, double n_sigma_left, double n_sigma_right) {
	
	int cut_val = 0;
	
	double loc_mu = 0.;
	for(int ipar=0; ipar<4; ipar++) loc_mu += (m_deltaPhi_mu_pars_two[0][ipar]*pow(eb,(double)ipar));
	
	double loc_sigma = 0.;
	for(int ipar=0; ipar<4; ipar++) loc_sigma += (m_deltaPhi_sigma_pars_two[0][ipar]*pow(eb,(double)ipar));
	
	double loc_diff = deltaPhi - loc_mu;
	if((-1.*n_sigma_left*loc_sigma < loc_diff) && (loc_diff < n_sigma_right*loc_sigma)) return 1;
	else return 0;
}

int ComptonAna::cut_deltaK_two(double deltaK, double eb, double n_sigma_left, double n_sigma_right) {
	
	int cut_val = 0;
	
	double loc_mu = 0.;
	for(int ipar=0; ipar<4; ipar++) loc_mu += (m_deltaK_mu_pars_two[0][ipar]*pow(eb,(double)ipar));
	
	double loc_sigma = 0.;
	for(int ipar=0; ipar<4; ipar++) loc_sigma += (m_deltaK_sigma_pars_two[0][ipar]*pow(eb,(double)ipar));
	
	double loc_diff = deltaK - loc_mu;
	if((-1.*n_sigma_left*loc_sigma < loc_diff) && (loc_diff < n_sigma_right*loc_sigma)) return 1;
	else return 0;
}

double ComptonAna::smear_deltaE(double deltaE, double eb) {
	
	if(!m_SHIFT_DISTRIBUTIONS) return deltaE;
	
	double loc_mu_data = 0., loc_mu_mc = 0.;
	for(int ipar=0; ipar<4; ipar++) {
		loc_mu_data += (m_deltaE_mu_pars[0][ipar]*pow(eb,(double)ipar));
		loc_mu_mc   += (m_deltaE_mu_pars[1][ipar]*pow(eb,(double)ipar));
	}
	double loc_sigma_data = sqrt(pow(m_deltaE_sigma_pars[0][0],2.0) + pow(m_deltaE_sigma_pars[0][1]/sqrt(eb),2.0) 
		+ pow(m_deltaE_sigma_pars[0][2]/eb,2.0));
	double loc_sigma_mc   = sqrt(pow(m_deltaE_sigma_pars[1][0],2.0) + pow(m_deltaE_sigma_pars[1][1]/sqrt(eb),2.0) 
		+ pow(m_deltaE_sigma_pars[1][2]/eb,2.0));
	
	loc_sigma_data *= eb;
	loc_sigma_mc   *= eb;
	
	double deltaE_smeared = deltaE + (loc_mu_data - loc_mu_mc);
	
	if(!m_SMEAR_DISTRIBUTIONS) return deltaE_smeared;
	
	if(loc_sigma_data > loc_sigma_mc) {
		double loc_smear = sqrt(pow(loc_sigma_data,2.0) - pow(loc_sigma_mc,2.0));
		deltaE_smeared += m_random->Gaus(0.,loc_smear);
	}
	
	return deltaE_smeared;
}

double ComptonAna::smear_deltaPhi(double deltaPhi, double eb) {
	
	if(!m_SHIFT_DISTRIBUTIONS) return deltaPhi;
	
	double loc_mu_data = 0., loc_mu_mc = 0.;
	for(int ipar=0; ipar<4; ipar++) {
		loc_mu_data += (m_deltaPhi_mu_pars[0][ipar]*pow(eb,(double)ipar));
		loc_mu_mc   += (m_deltaPhi_mu_pars[1][ipar]*pow(eb,(double)ipar));
	}
	
	double loc_sigma_data = 0., loc_sigma_mc = 0.;
	for(int ipar=0; ipar<4; ipar++) {
		loc_sigma_data += (m_deltaPhi_sigma_pars[0][ipar]*pow(eb,(double)ipar));
		loc_sigma_mc   += (m_deltaPhi_sigma_pars[1][ipar]*pow(eb,(double)ipar));
	}
	
	double deltaPhi_smeared = deltaPhi + (loc_mu_data - loc_mu_mc);
	
	if(!m_SMEAR_DISTRIBUTIONS) return deltaPhi_smeared;
	
	if(loc_sigma_data > loc_sigma_mc) {
		double loc_smear = sqrt(pow(loc_sigma_data,2.0) - pow(loc_sigma_mc,2.0));
		deltaPhi_smeared += m_random->Gaus(0.,loc_smear);
	}
	
	return deltaPhi_smeared;
}

double ComptonAna::smear_deltaK(double deltaK, double eb) {
	
	if(!m_SHIFT_DISTRIBUTIONS) return deltaK;
	
	double loc_mu_data = 0., loc_mu_mc = 0.;
	for(int ipar=0; ipar<4; ipar++) {
		loc_mu_data += (m_deltaK_mu_pars[0][ipar]*pow(eb,(double)ipar));
		loc_mu_mc   += (m_deltaK_mu_pars[1][ipar]*pow(eb,(double)ipar));
	}
	
	double loc_sigma_data = 0., loc_sigma_mc = 0.;
	for(int ipar=0; ipar<4; ipar++) {
		loc_sigma_data += (m_deltaK_sigma_pars[0][ipar]*pow(eb,(double)ipar));
		loc_sigma_mc   += (m_deltaK_sigma_pars[1][ipar]*pow(eb,(double)ipar));
	}
	
	double deltaK_smeared = deltaK + (loc_mu_data - loc_mu_mc);
	
	if(!m_SMEAR_DISTRIBUTIONS) return deltaK_smeared;
	
	if(loc_sigma_data > loc_sigma_mc) {
		double loc_smear = sqrt(pow(loc_sigma_data,2.0) - pow(loc_sigma_mc,2.0));
		deltaK_smeared += m_random->Gaus(0.,loc_smear);
	}
	
	return deltaK_smeared;
}

double ComptonAna::smear_deltaE_two(double deltaE, double eb) {
	
	if(!m_SHIFT_DISTRIBUTIONS) return deltaE;
	
	double loc_mu_data = 0., loc_mu_mc = 0.;
	for(int ipar=0; ipar<4; ipar++) {
		loc_mu_data += (m_deltaE_mu_pars_two[0][ipar]*pow(eb,(double)ipar));
		loc_mu_mc   += (m_deltaE_mu_pars_two[1][ipar]*pow(eb,(double)ipar));
	}
	double loc_sigma_data = sqrt(pow(m_deltaE_sigma_pars_two[0][0],2.0) + pow(m_deltaE_sigma_pars_two[0][1]/sqrt(eb),2.0) 
		+ pow(m_deltaE_sigma_pars_two[0][2]/eb,2.0));
	double loc_sigma_mc   = sqrt(pow(m_deltaE_sigma_pars_two[1][0],2.0) + pow(m_deltaE_sigma_pars_two[1][1]/sqrt(eb),2.0) 
		+ pow(m_deltaE_sigma_pars_two[1][2]/eb,2.0));
	
	loc_sigma_data *= eb;
	loc_sigma_mc   *= eb;
	
	double deltaE_smeared = deltaE + (loc_mu_data - loc_mu_mc);
	
	if(!m_SMEAR_DISTRIBUTIONS) return deltaE_smeared;
	
	if(loc_sigma_data > loc_sigma_mc) {
		double loc_smear = sqrt(pow(loc_sigma_data,2.0) - pow(loc_sigma_mc,2.0));
		deltaE_smeared += m_random->Gaus(0.,loc_smear);
	}
	
	return deltaE_smeared;
}

double ComptonAna::smear_deltaPhi_two(double deltaPhi, double eb) {
	
	if(!m_SHIFT_DISTRIBUTIONS) return deltaPhi;
	
	double loc_mu_data = 0., loc_mu_mc = 0.;
	for(int ipar=0; ipar<4; ipar++) {
		loc_mu_data += (m_deltaPhi_mu_pars_two[0][ipar]*pow(eb,(double)ipar));
		loc_mu_mc   += (m_deltaPhi_mu_pars_two[1][ipar]*pow(eb,(double)ipar));
	}
	
	double loc_sigma_data = 0., loc_sigma_mc = 0.;
	for(int ipar=0; ipar<4; ipar++) {
		loc_sigma_data += (m_deltaPhi_sigma_pars_two[0][ipar]*pow(eb,(double)ipar));
		loc_sigma_mc   += (m_deltaPhi_sigma_pars_two[1][ipar]*pow(eb,(double)ipar));
	}
	
	double deltaPhi_smeared = deltaPhi + (loc_mu_data - loc_mu_mc);
	
	if(!m_SMEAR_DISTRIBUTIONS) return deltaPhi_smeared;
	
	if(loc_sigma_data > loc_sigma_mc) {
		double loc_smear = sqrt(pow(loc_sigma_data,2.0) - pow(loc_sigma_mc,2.0));
		deltaPhi_smeared += m_random->Gaus(0.,loc_smear);
	}
	
	return deltaPhi_smeared;
}

double ComptonAna::smear_deltaK_two(double deltaK, double eb) {
	
	if(!m_SHIFT_DISTRIBUTIONS) return deltaK;
	
	double loc_mu_data = 0., loc_mu_mc = 0.;
	for(int ipar=0; ipar<4; ipar++) {
		loc_mu_data += (m_deltaK_mu_pars_two[0][ipar]*pow(eb,(double)ipar));
		loc_mu_mc   += (m_deltaK_mu_pars_two[1][ipar]*pow(eb,(double)ipar));
	}
	
	double loc_sigma_data = 0., loc_sigma_mc = 0.;
	for(int ipar=0; ipar<4; ipar++) {
		loc_sigma_data += (m_deltaK_sigma_pars_two[0][ipar]*pow(eb,(double)ipar));
		loc_sigma_mc   += (m_deltaK_sigma_pars_two[1][ipar]*pow(eb,(double)ipar));
	}
	
	double deltaK_smeared = deltaK + (loc_mu_data - loc_mu_mc);
	
	if(!m_SMEAR_DISTRIBUTIONS) return deltaK_smeared;
	
	if(loc_sigma_data > loc_sigma_mc) {
		double loc_smear = sqrt(pow(loc_sigma_data,2.0) - pow(loc_sigma_mc,2.0));
		deltaK_smeared += m_random->Gaus(0.,loc_smear);
	}
	
	return deltaK_smeared;
}

int ComptonAna::loadCutParameters() {
	
	char cut_dir[256];
	sprintf(cut_dir, "/work/halld/home/andrsmit/primex_compton_analysis/analyze_trees/phase1/cuts/Run%06d/200nA", m_runNumber);
	
	char buf[256];
	ifstream loc_inf;
	double loc_pars[4];
	
	vector<string> type_str = {"data","mc"};
	
	int ii = 0;
	for(vector<string>::iterator it = type_str.begin(); it != type_str.end(); it++,ii++) {
		
		//---------------------------------------//
		// DeltaE:
		
		// Mu:
		sprintf(buf, "%s/%s/deltaE_mu.dat", cut_dir, (*it).c_str());
		if(gSystem->AccessPathName(buf)) return 1;
		
		loc_inf.open(buf);
		for(int ipar=0; ipar<4; ipar++) {
			loc_inf >> loc_pars[ipar];
			m_deltaE_mu_pars[ii][ipar] = loc_pars[ipar];
		}
		loc_inf.close();
		
		// Sigma:
		sprintf(buf, "%s/%s/deltaE_sigma.dat", cut_dir, (*it).c_str());
		if(gSystem->AccessPathName(buf)) return 2;
		
		loc_inf.open(buf);
		for(int ipar=0; ipar<3; ipar++) {
			loc_inf >> loc_pars[ipar];
			m_deltaE_sigma_pars[ii][ipar] = loc_pars[ipar];
		}
		loc_inf.close();
		
		// Mu (Two layer cut):
		sprintf(buf, "%s/%s/deltaE_mu_two.dat", cut_dir, (*it).c_str());
		if(gSystem->AccessPathName(buf)) return 3;
		
		loc_inf.open(buf);
		for(int ipar=0; ipar<4; ipar++) {
			loc_inf >> loc_pars[ipar];
			m_deltaE_mu_pars_two[ii][ipar] = loc_pars[ipar];
		}
		loc_inf.close();
		
		// Sigma:
		sprintf(buf, "%s/%s/deltaE_sigma_two.dat", cut_dir, (*it).c_str());
		if(gSystem->AccessPathName(buf)) return 4;
		
		loc_inf.open(buf);
		for(int ipar=0; ipar<3; ipar++) {
			loc_inf >> loc_pars[ipar];
			m_deltaE_sigma_pars_two[ii][ipar] = loc_pars[ipar];
		}
		loc_inf.close();
		
		//---------------------------------------//
		// DeltaPhi:
		
		// Mu:
		sprintf(buf, "%s/%s/deltaPhi_mu.dat", cut_dir, (*it).c_str());
		if(gSystem->AccessPathName(buf)) return 5;
		
		loc_inf.open(buf);
		for(int ipar=0; ipar<4; ipar++) {
			loc_inf >> loc_pars[ipar];
			m_deltaPhi_mu_pars[ii][ipar] = loc_pars[ipar];
		}
		loc_inf.close();
		
		// Sigma:
		sprintf(buf, "%s/%s/deltaPhi_sigma.dat", cut_dir, (*it).c_str());
		if(gSystem->AccessPathName(buf)) return 6;
		
		loc_inf.open(buf);
		for(int ipar=0; ipar<4; ipar++) {
			loc_inf >> loc_pars[ipar];
			m_deltaPhi_sigma_pars[ii][ipar] = loc_pars[ipar];
		}
		loc_inf.close();
		
		// Mu (Two layer cut):
		sprintf(buf, "%s/%s/deltaPhi_mu_two.dat", cut_dir, (*it).c_str());
		if(gSystem->AccessPathName(buf)) return 7;
		
		loc_inf.open(buf);
		for(int ipar=0; ipar<4; ipar++) {
			loc_inf >> loc_pars[ipar];
			m_deltaPhi_mu_pars_two[ii][ipar] = loc_pars[ipar];
		}
		loc_inf.close();
		
		// Sigma:
		sprintf(buf, "%s/%s/deltaPhi_sigma_two.dat", cut_dir, (*it).c_str());
		if(gSystem->AccessPathName(buf)) return 8;
		
		loc_inf.open(buf);
		for(int ipar=0; ipar<4; ipar++) {
			loc_inf >> loc_pars[ipar];
			m_deltaPhi_sigma_pars_two[ii][ipar] = loc_pars[ipar];
		}
		loc_inf.close();
		
		//---------------------------------------//
		// DeltaPhi:
		
		// Mu:
		sprintf(buf, "%s/%s/deltaK_mu.dat", cut_dir, (*it).c_str());
		if(gSystem->AccessPathName(buf)) return 9;
		
		loc_inf.open(buf);
		for(int ipar=0; ipar<4; ipar++) {
			loc_inf >> loc_pars[ipar];
			m_deltaK_mu_pars[ii][ipar] = loc_pars[ipar];
		}
		loc_inf.close();
		
		// Sigma:
		sprintf(buf, "%s/%s/deltaK_sigma.dat", cut_dir, (*it).c_str());
		if(gSystem->AccessPathName(buf)) return 10;
		
		loc_inf.open(buf);
		for(int ipar=0; ipar<4; ipar++) {
			loc_inf >> loc_pars[ipar];
			m_deltaK_sigma_pars[ii][ipar] = loc_pars[ipar];
		}
		loc_inf.close();
		
		// Mu (Two layer cut):
		sprintf(buf, "%s/%s/deltaK_mu_two.dat", cut_dir, (*it).c_str());
		if(gSystem->AccessPathName(buf)) return 11;
		
		loc_inf.open(buf);
		for(int ipar=0; ipar<4; ipar++) {
			loc_inf >> loc_pars[ipar];
			m_deltaK_mu_pars_two[ii][ipar] = loc_pars[ipar];
		}
		loc_inf.close();
		
		// Sigma:
		sprintf(buf, "%s/%s/deltaK_sigma_two.dat", cut_dir, (*it).c_str());
		if(gSystem->AccessPathName(buf)) return 12;
		
		loc_inf.open(buf);
		for(int ipar=0; ipar<4; ipar++) {
			loc_inf >> loc_pars[ipar];
			m_deltaK_sigma_pars_two[ii][ipar] = loc_pars[ipar];
		}
		loc_inf.close();
		
		//---------------------------------------//
	}
	
	return 0;
}

void ComptonAna::setRunNumber(int runNum) { 
	
	m_runNumber = runNum;
	setGeometry();
	int load_val = loadCutParameters();
	if(load_val != 0) {
		cout << "\n\n\n" << endl;
		cout << "Problem loading cut parameters for Run " << m_runNumber << "!!! code: " << load_val << endl;
		cout << "\n\n\n" << endl;
	}
	
	return;
}

void ComptonAna::setGeometry() {
	
	if(m_runNumber>60000 && m_runNumber<69999) {
		
		m_phase_val = 1;
		
		if(m_runNumber<61355) {
			m_target_length  = 1.7755;
			m_target_density = 1.848;
			m_target_atten   = 0.01172;
			m_vertex.SetZ(64.935);
		} else {
			m_target_length  = 29.535;
			m_target_density = 0.1217;
			m_target_atten   = 0.00821;
			m_vertex.SetZ(65.0);
		}
		
		// Correction to alignment after Simon updated beam spot:
		
		m_fcal_face.SetXYZ(0.455, -0.032, 624.906);
		m_fcal_correction.SetXYZ(0.455-0.529, -0.032+0.002, 0.0);
		
		m_ccal_face.SetX(-0.082);
		m_ccal_correction.SetX(-0.082+0.0225);
		
		if(m_runNumber<61483) {
			m_ccal_face.SetY(0.061);
			m_ccal_correction.SetY(0.061-0.0073);
		} else {
			m_ccal_face.SetY(0.051);
			m_ccal_correction.SetY(0.051-0.0073);
		}
		
		if(m_runNumber<61483) {
			m_vertex.SetX( 0.027);
			m_vertex.SetY(-0.128);
		} else if(m_runNumber<61774) {
			m_vertex.SetX( 0.001);
			m_vertex.SetY(-0.077);
		} else {
			m_vertex.SetX( 0.038);
			m_vertex.SetY(-0.095);
		}
	}
	
	return;
}

void ComptonAna::setOutputFileName(string name) { m_output_fname = name; return; }

int ComptonAna::loadTree() {
	
	m_tree = (TTree*)m_infile->Get("primex_compton");
	if(m_tree==NULL) return 0;
	
	// Reset event count to zero when laoding a new Tree:
	m_event = 0;
	
	return m_tree->GetEntries();
}

void ComptonAna::readEvent() {
	
	if(m_event == 0) {
		
		// Set Branch Address on first Event:
		
		m_tree->SetBranchAddress("rfTime",             &m_rfTime);
		m_tree->SetBranchAddress("nbeam",              &m_nbeam);
		m_tree->SetBranchAddress("tag_counter",        &m_tag_counter);
		m_tree->SetBranchAddress("tag_sys",            &m_tag_sys);
		m_tree->SetBranchAddress("beam_e",             &m_beam_e);
		m_tree->SetBranchAddress("beam_t",             &m_beam_t);
		//m_tree->SetBranchAddress("acc_scale_factor",   &m_acc_scale_factor);
		m_tree->SetBranchAddress("nfcal",              &m_nfcal);
		m_tree->SetBranchAddress("fcal_e",             &m_fcal_e);
		m_tree->SetBranchAddress("fcal_x",             &m_fcal_x);
		m_tree->SetBranchAddress("fcal_y",             &m_fcal_y);
		m_tree->SetBranchAddress("fcal_z",             &m_fcal_z);
		m_tree->SetBranchAddress("fcal_t",             &m_fcal_t);
		m_tree->SetBranchAddress("fcal_nblocks",       &m_fcal_nblocks);
		m_tree->SetBranchAddress("nccal",              &m_nccal);
		m_tree->SetBranchAddress("ccal_e",             &m_ccal_e);
		m_tree->SetBranchAddress("ccal_x1",            &m_ccal_x);
		m_tree->SetBranchAddress("ccal_y1",            &m_ccal_y);
		m_tree->SetBranchAddress("ccal_z",             &m_ccal_z);
		m_tree->SetBranchAddress("ccal_t",             &m_ccal_t);
		m_tree->SetBranchAddress("ccal_idmax",         &m_ccal_idmax);
		m_tree->SetBranchAddress("ccal_nblocks",       &m_ccal_nblocks);
		m_tree->SetBranchAddress("ntof",               &m_ntof);
		m_tree->SetBranchAddress("tof_x",              &m_tof_x);
		m_tree->SetBranchAddress("tof_y",              &m_tof_y);
		m_tree->SetBranchAddress("tof_z",              &m_tof_z);
		m_tree->SetBranchAddress("tof_t",              &m_tof_t);
		m_tree->SetBranchAddress("nmc",                &m_nmc);
		m_tree->SetBranchAddress("mc_pdgtype",         &m_mc_pdgtype);
		m_tree->SetBranchAddress("mc_x",               &m_mc_x);
		m_tree->SetBranchAddress("mc_y",               &m_mc_y);
		m_tree->SetBranchAddress("mc_z",               &m_mc_z);
		m_tree->SetBranchAddress("mc_t",               &m_mc_t);
		m_tree->SetBranchAddress("mc_e",               &m_mc_e);
		m_tree->SetBranchAddress("mc_p",               &m_mc_p);
		m_tree->SetBranchAddress("mc_theta",           &m_mc_theta);
		m_tree->SetBranchAddress("mc_phi",             &m_mc_phi);
		m_tree->SetBranchAddress("mc_reaction_type",   &m_reaction_type);
		m_tree->SetBranchAddress("mc_reaction_weight", &m_reaction_weight);
		m_tree->SetBranchAddress("mc_reaction_energy", &m_reaction_energy);
	}
	
	m_tree->GetEvent(m_event);
	
	return;
}

void ComptonAna::initHistograms() {
	
	h_thrown_photon_angle          = new TH1F("thrown_photon_angle", 
		"Thrown Photon Angle; #theta_{#gamma} [#circ]", 1000, 0., 10.);
	h_thrown_photon_angle_accepted = new TH1F("thrown_photon_angle_accepted", 
		"Thrown Photon Angle; #theta_{#gamma} [#circ]", 1000, 0., 10.);
	
	h_reaction_weight        = new TH1F("reaction_weight",        "Reaction Weight",                  10000, 0., 1.e4);
	h_reaction_weight_double = new TH1F("reaction_weight_double", "Reaction Weight (Double Compton)", 10000, 0., 1.e4);
	
	h_vertex          = new TH1F("vertex",          "Vertex Z Position (unweighted)", 1000, 0., 100.);
	h_vertex_accepted = new TH1F("vertex_accepted", "Vertex Z Position (weighted)",   1000, 0., 100.);
	
	h_fcal_rf_dt = new TH1F("fcal_rf_dt", "t_{FCAL} - t_{RF}; [ns]", 2000, -100., 100.);
	h_ccal_rf_dt = new TH1F("ccal_rf_dt", "t_{FCAL} - t_{RF}; [ns]", 2000, -100., 100.);
	h_beam_rf_dt = new TH1F("beam_rf_dt", "t_{FCAL} - t_{RF}; [ns]", 2000, -100., 100.);
	
	h_thrown_mass           = new TH1F("thrown_mass", "; M_{e^{+}e^{-}}(thrown) [GeV/c^{2}]", 1000,  0., 1.);
	
	h_deltaE_vs_thrown_mass = new TH2F("deltaE_vs_thrown_mass", "; M_{e^{+}e^{-}}(thrown) [GeV/c^{2}]; #DeltaE [GeV]",
		1000, 0.0, 1.0, 1000, -8.0, 8.0);
	h_deltaE_vs_rec_mass = new TH2F("deltaE_vs_rec_mass", "; M_{e^{+}e^{-}}(reconstructed) [GeV/c^{2}]; #DeltaE [GeV]",
		1000, 0.0, 1.0, 1000, -8.0, 8.0);
	h_deltaK_vs_thrown_mass = new TH2F("deltaK_vs_thrown_mass", "; M_{e^{+}e^{-}}(thrown) [GeV/c^{2}]; #DeltaK [GeV]",
		1000, 0.0, 1.0, 1000, -8.0, 8.0);
	h_deltaK_vs_rec_mass = new TH2F("deltaK_vs_rec_mass", "; M_{e^{+}e^{-}}(reconstructed) [GeV/c^{2}]; #DeltaK [GeV]",
		1000, 0.0, 1.0, 1000, -8.0, 8.0);
	h_elas_low_mass = new TH2F("elas_low_mass", "; #DeltaE [GeV]; #DeltaK [GeV]",
		1000, -8.0, 8.0, 1000, -8.0, 8.0);
	h_elas_high_mass = new TH2F("elas_high_mass", "; #DeltaE [GeV]; #DeltaK [GeV]",
		1000, -8.0, 8.0, 1000, -8.0, 8.0);
	
	h_rec_vs_thrown_mass = new TH2F("rec_vs_thrown_mass", 
		"; M_{e^{+}e^{-}}(thrown) [GeV/c^{2}]; M_{e^{+}e^{-}}(reconstructed) [GeV/c^{2}]",
		1000, 0.0, 1.0, 1000, 0.0, 1.0);
	
	for(int ihist=0; ihist<6; ihist++) {
		h_deltaE_pair[ihist]   = new TH2F(Form("deltaE_pair_%d",ihist),   "#DeltaE vs. E_{#gamma}; [GeV]; [GeV]", 
			120, 0., 12., 2000, -8.0, 8.0);
		h_deltaK_pair[ihist]   = new TH2F(Form("deltaK_pair_%d",ihist),   "#DeltaK vs. E_{#gamma}; [GeV]; [GeV]", 
			120, 0., 12., 1000, -8.0, 8.0);
		h_deltaPhi_pair[ihist] = new TH2F(Form("deltaPhi_pair_%d",ihist), "#Delta#phi vs. E_{#gamma}; [GeV]; [#circ]", 
			120, 0., 12., 2000, 0.0, 360.0);
		
		h_deltaE_triplet[ihist]   = new TH2F(Form("deltaE_triplet_%d",ihist),   "#DeltaE vs. E_{#gamma}; [GeV]; [GeV]", 
			120, 0., 12., 2000, -8.0, 8.0);
		h_deltaK_triplet[ihist]   = new TH2F(Form("deltaK_triplet_%d",ihist),   "#DeltaK vs. E_{#gamma}; [GeV]; [GeV]", 
			120, 0., 12., 1000, -8.0, 8.0);
		h_deltaPhi_triplet[ihist] = new TH2F(Form("deltaPhi_triplet_%d",ihist), "#Delta#phi vs. E_{#gamma}; [GeV]; [#circ]", 
			120, 0., 12., 2000, 0.0, 360.0);
	}
	
	//---------------------------------------------//
	for(int ihist=0; ihist<m_n_cuts; ihist++) {
		
		h_deltaE_tagh[ihist] = new TH2F(Form("deltaE_tagh_%d",ihist),
			"#DeltaE; TAGH Counter; E_{1} + E_{2} - E_{#gamma} [GeV]",
			274, 0.5, 274.5, 2000, -8.0, 8.0);
		h_deltaE_tagm[ihist] = new TH2F(Form("deltaE_tagm_%d",ihist),
			"#DeltaE; TAGM Counter; E_{1} + E_{2} - E_{#gamma} [GeV]",
			102, 0.5, 102.5, 2000, -8.0, 8.0);
		h_deltaE_tagh[ihist]->Sumw2();
		h_deltaE_tagm[ihist]->Sumw2();
		
		h_deltaPhi_tagh[ihist] = new TH2F(Form("deltaPhi_tagh_%d",ihist),
			"#Delta#phi; TAGH Counter; |#phi_{2} - #phi_{1}| (#circ)",
			274, 0.5, 274.5, 1800, 0.0, 360.0);
		h_deltaPhi_tagm[ihist] = new TH2F(Form("deltaPhi_tagm_%d",ihist),
			"#Delta#phi; TAGM Counter; |#phi_{2} - #phi_{1}| (#circ)",
			102, 0.5, 102.5, 1800, 0.0, 360.0);
		h_deltaPhi_tagh[ihist]->Sumw2();
		h_deltaPhi_tagm[ihist]->Sumw2();
		
		h_deltaK_tagh[ihist] = new TH2F(Form("deltaK_tagh_%d",ihist),
			"#DeltaK; TAGH Counter; E_{Compton}#left(#theta_{1},#theta_{2}#right) - E_{#gamma} [GeV]",
			274, 0.5, 274.5, 1000, -8.0, 8.0);
		h_deltaK_tagm[ihist] = new TH2F(Form("deltaK_tagm_%d",ihist),
			"#DeltaK; TAGM Counter; E_{Compton}#left(#theta_{1},#theta_{2}#right) - E_{#gamma} [GeV]",
			102, 0.5, 102.5, 1000, -8.0, 8.0);
		h_deltaK_tagh[ihist]->Sumw2();
		h_deltaK_tagm[ihist]->Sumw2();
		
		h_elas_vs_deltaE[ihist] = new TH2F(Form("elas_vs_deltaE_%d",ihist),
			"2-D Elasticity; E_{1}+E_{2} - E_{#gamma} [GeV]; E_{1}+E_{2} - E_{Compton} [GeV]",
			500, -8.0, 8.0, 500, -8.0, 8.0);
		
		h_deltaK_vs_deltaE[ihist] = new TH2F(Form("deltaK_vs_deltaE_%d",ihist),
			"2-D Elasticity; E_{1}+E_{2} - E_{Compton} [GeV]; E_{Comp} - E_{#gamma} [GeV]",
			500, -8.0, 8.0, 500, -8.0, 8.0);
		
		h_deltaCCAL_vs_deltaE[ihist] = new TH2F(Form("deltaCCAL_vs_deltaE_%d",ihist),
			"2-D Elasticity; E_{1}+E_{2} - E_{#gamma} [GeV]; E_{CCAL} - E_{Comp} [GeV]",
			500, -8.0, 8.0, 500, -8.0, 8.0);
		h_deltaCCAL_vs_deltaK[ihist] = new TH2F(Form("deltaCCAL_vs_deltaK_%d",ihist),
			"2-D Elasticity; E_{1}+E_{2} - E_{Comp} [GeV]; E_{CCAL} - E_{Comp} [GeV]",
			500, -8.0, 8.0, 500, -8.0, 8.0);
		
		h_deltaFCAL_vs_deltaE[ihist] = new TH2F(Form("deltaFCAL_vs_deltaE_%d",ihist),
			"2-D Elasticity; E_{1}+E_{2} - E_{#gamma} [GeV]; E_{FCAL} - E_{Comp} [GeV]",
			500, -8.0, 8.0, 500, -8.0, 8.0);
		h_deltaFCAL_vs_deltaK[ihist] = new TH2F(Form("deltaFCAL_vs_deltaK_%d",ihist),
			"2-D Elasticity; E_{1}+E_{2} - E_{Comp} [GeV]; E_{FCAL} - E_{Comp} [GeV]",
			500, -8.0, 8.0, 500, -8.0, 8.0);
		
		h_deltaFCAL_vs_deltaCCAL[ihist] = new TH2F(Form("deltaFCAL_vs_deltaCCAL_%d",ihist),
			"E_{FCAL} - E_{Comp} vs. E_{CCAL} - E_{Comp}; E_{CCAL} - E_{Comp} [GeV]; E_{FCAL} - E_{Comp} [GeV]",
			500, -8.0, 8.0, 500, -8.0, 8.0);
		
		h_ccal_xy[ihist] = new TH2F(Form("ccal_xy_%d",ihist),
			"CCAL Shower Occupancy; x_{CCAL} [cm]; y_{CCAL} [cm]",
			500,  -13.,  13.,  500,  -13.,  13.);
		h_fcal_xy[ihist] = new TH2F(Form("fcal_xy_%d",ihist),
			"FCAL Shower Occupancy; x_{FCAL} [cm]; y_{FCAL} [cm]",
			1000, -100., 100., 1000, -100., 100.);
		
		h_opangle_tagh[ihist]       = new TH2F(Form("opangle_tagh_%d",ihist), 
			"Opening Angle; TAGH Counter; #theta_{12} [#circ]", 274, 0.5, 274.5, 1000, 0., 10.);
		h_opangle_tagh_ecut[ihist]  = new TH2F(Form("opangle_tagh_ecut_%d",ihist), 
			"Opening Angle (#DeltaE cut); TAGH Counter; #theta_{12} [#circ]", 274, 0.5, 274.5, 1000, 0., 10.);
		h_opangle_tagh_ekcut[ihist] = new TH2F(Form("opangle_tagh_ekcut_%d",ihist), 
			"Opening Angle (#DeltaE and #DeltaK cuts); TAGH Counter; #theta_{12} [#circ]", 274, 0.5, 274.5, 1000, 0., 10.);
		
		h_opangle_tagm[ihist]       = new TH2F(Form("opangle_tagm_%d",ihist), 
			"Opening Angle; TAGM Counter; #theta_{12} [#circ]", 102, 0.5, 102.5, 1000, 0., 10.);
		h_opangle_tagm_ecut[ihist]  = new TH2F(Form("opangle_tagm_ecut_%d",ihist), 
			"Opening Angle (#DeltaE cut); TAGM Counter; #theta_{12} [#circ]", 102, 0.5, 102.5, 1000, 0., 10.);
		h_opangle_tagm_ekcut[ihist] = new TH2F(Form("opangle_tagm_ekcut_%d",ihist), 
			"Opening Angle (#DeltaE and #DeltaK cuts); TAGM Counter; #theta_{12} [#circ]", 102, 0.5, 102.5, 1000, 0., 10.);
	}
	
	h_ccal_nblocks     = new TH1F("ccal_nblocks",     "Number of hits in CCAL Cluster", 50, -0.5, 50.5);
	h_fcal_nblocks     = new TH1F("fcal_nblocks",     "Number of hits in FCAL Cluster", 50, -0.5, 50.5);
	h_ccal_nblocks_cut = new TH1F("ccal_nblocks_cut", "Number of hits in CCAL Cluster (#DeltaE > 0.65 GeV)", 
		50, -0.5, 50.5);
	h_fcal_nblocks_cut = new TH1F("fcal_nblocks_cut", "Number of hits in FCAL Cluster (#DeltaE > 0.65 GeV)", 
		50, -0.5, 50.5);
	
	h_ccal_xy_highdeltaE = new TH2F("ccal_xy_highdeltaE", 
		"CCAL Shower Occupancy (#DeltaE > 4#sigma); x_{CCAL} [cm]; y_{CCAL} [cm]",  500,   -13.,   13.,  500,   -13.,   13.);
	h_fcal_xy_highdeltaE = new TH2F("fcal_xy_highdeltaE", 
		"FCAL Shower Occupancy (#DeltaE > 4#sigma); x_{CCAL} [cm]; y_{CCAL} [cm]",  500,  -100.,  100.,  500,  -100.,  100.);
	h_ccal_xy_lowdeltaE = new TH2F("ccal_xy_lowdeltaE", 
		"CCAL Shower Occupancy (#DeltaE < -3.5 GeV); x_{CCAL} [cm]; y_{CCAL} [cm]", 500,   -13.,   13.,  500,   -13.,   13.);
	h_fcal_xy_lowdeltaE = new TH2F("fcal_xy_lowdeltaE", 
		"FCAL Shower Occupancy (#DeltaE > -3.5 GeV); x_{CCAL} [cm]; y_{CCAL} [cm]", 500,  -100.,  100.,  500,  -100.,  100.);
	
	h_sumPhi_vs_deltaPhi = new TH2F("sumPhi_vs_deltaPhi", 
		";|#phi_{1}-#phi_{2}| / 2 [#circ]; (#phi_{1}+#phi_{2}) / 2 [#circ]", 500, 50., 130., 1000, -90., 90.);
	h_sumPhi_vs_deltaPhi->GetXaxis()->CenterTitle(true);
	h_sumPhi_vs_deltaPhi->GetYaxis()->CenterTitle(true);
	
	return;
}

void ComptonAna::resetHistograms() {
	
	h_thrown_photon_angle->Reset();
	h_thrown_photon_angle_accepted->Reset();
	
	h_reaction_weight->Reset();
	h_reaction_weight_double->Reset();
	
	h_vertex->Reset();
	h_vertex_accepted->Reset();
	
	h_fcal_rf_dt->Reset();
	h_ccal_rf_dt->Reset();
	h_beam_rf_dt->Reset();
	
	h_thrown_mass->Reset();
	h_deltaE_vs_thrown_mass->Reset();
	h_deltaE_vs_rec_mass->Reset();
	h_deltaK_vs_thrown_mass->Reset();
	h_deltaK_vs_rec_mass->Reset();
	h_elas_low_mass->Reset();
	h_elas_high_mass->Reset();
	h_rec_vs_thrown_mass->Reset();
	
	for(int ihist=0; ihist<6; ihist++) {
		h_deltaE_pair[ihist]->Reset();
		h_deltaK_pair[ihist]->Reset();
		h_deltaPhi_pair[ihist]->Reset();
		h_deltaE_triplet[ihist]->Reset();
		h_deltaK_triplet[ihist]->Reset();
		h_deltaPhi_triplet[ihist]->Reset();
	}
	
	for(int ihist=0; ihist<m_n_cuts; ihist++) {
		h_deltaE_tagh[ihist]->Reset();
		h_deltaE_tagm[ihist]->Reset();
	}
	for(int ihist=0; ihist<m_n_cuts; ihist++) {
		h_deltaPhi_tagh[ihist]->Reset();
		h_deltaPhi_tagm[ihist]->Reset();
	}
	for(int ihist=0; ihist<m_n_cuts; ihist++) {
		h_deltaK_tagh[ihist]->Reset();
		h_deltaK_tagm[ihist]->Reset();
	}
	for(int ihist=0; ihist<m_n_cuts; ihist++) {
		h_elas_vs_deltaE[ihist]->Reset();
	}
	for(int ihist=0; ihist<m_n_cuts; ihist++) {
		h_deltaK_vs_deltaE[ihist]->Reset();
	}
	for(int ihist=0; ihist<m_n_cuts; ihist++) {
		h_deltaCCAL_vs_deltaE[ihist]->Reset();
		h_deltaCCAL_vs_deltaK[ihist]->Reset();
	}
	for(int ihist=0; ihist<m_n_cuts; ihist++) {
		h_deltaFCAL_vs_deltaE[ihist]->Reset();
		h_deltaFCAL_vs_deltaK[ihist]->Reset();
	}
	for(int ihist=0; ihist<m_n_cuts; ihist++) {
		h_deltaFCAL_vs_deltaCCAL[ihist]->Reset();
	}
	for(int ihist=0; ihist<m_n_cuts; ihist++) {
		h_ccal_xy[ihist]->Reset();
		h_fcal_xy[ihist]->Reset();
	}
	for(int ihist=0; ihist<m_n_cuts; ihist++) {
		h_opangle_tagh[ihist]->Reset();
		h_opangle_tagh_ecut[ihist]->Reset();
		h_opangle_tagh_ekcut[ihist]->Reset();
		h_opangle_tagm[ihist]->Reset();
		h_opangle_tagm_ecut[ihist]->Reset();
		h_opangle_tagm_ekcut[ihist]->Reset();
	}
	
	h_ccal_nblocks->Reset();
	h_fcal_nblocks->Reset();
	h_ccal_nblocks_cut->Reset();
	h_fcal_nblocks_cut->Reset();
	
	h_ccal_xy_highdeltaE->Reset();
	h_fcal_xy_highdeltaE->Reset();
	h_ccal_xy_lowdeltaE->Reset();
	h_fcal_xy_lowdeltaE->Reset();
	
	h_sumPhi_vs_deltaPhi->Reset();
	
	return;
}

void ComptonAna::writeHistograms() {
	
	cout << "writing histograms to: " << m_output_fname << "..." << endl;
	TFile *fOut = new TFile(m_output_fname.c_str(), "RECREATE");
	fOut->cd();
	
	h_thrown_photon_angle->Write();
	h_thrown_photon_angle_accepted->Write();
	
	h_reaction_weight->Write();
	h_reaction_weight_double->Write();
	
	h_vertex->Write();
	h_vertex_accepted->Write();
	
	h_fcal_rf_dt->Write();
	h_ccal_rf_dt->Write();
	h_beam_rf_dt->Write();
	
	h_thrown_mass->Write();
	h_deltaE_vs_thrown_mass->Write();
	h_deltaE_vs_rec_mass->Write();
	h_deltaK_vs_thrown_mass->Write();
	h_deltaK_vs_rec_mass->Write();
	h_elas_low_mass->Write();
	h_elas_high_mass->Write();
	h_rec_vs_thrown_mass->Write();
	
	TDirectory *dir_cuts = new TDirectoryFile("cuts", "cuts");
	dir_cuts->cd();
	for(int ihist=0; ihist<6; ihist++) {
		h_deltaE_pair[ihist]->Write();
		h_deltaK_pair[ihist]->Write();
		h_deltaPhi_pair[ihist]->Write();
		h_deltaE_triplet[ihist]->Write();
		h_deltaK_triplet[ihist]->Write();
		h_deltaPhi_triplet[ihist]->Write();
	}
	dir_cuts->cd("../");
	
	TDirectory *dir_deltaE = new TDirectoryFile("deltaE", "deltaE");
	dir_deltaE->cd();
	for(int ihist=0; ihist<m_n_cuts; ihist++) {
		h_deltaE_tagh[ihist]->Write();
		h_deltaE_tagm[ihist]->Write();
	}
	dir_deltaE->cd("../");
	
	TDirectory *dir_deltaPhi = new TDirectoryFile("deltaPhi", "deltaPhi");
	dir_deltaPhi->cd();
	for(int ihist=0; ihist<m_n_cuts; ihist++) {
		h_deltaPhi_tagh[ihist]->Write();
		h_deltaPhi_tagm[ihist]->Write();
	}
	dir_deltaPhi->cd("../");
	
	TDirectory *dir_deltaK = new TDirectoryFile("deltaK", "deltaK");
	dir_deltaK->cd();
	for(int ihist=0; ihist<m_n_cuts; ihist++) {
		h_deltaK_tagh[ihist]->Write();
		h_deltaK_tagm[ihist]->Write();
	}
	dir_deltaK->cd("../");
	
	TDirectory *dir_elas = new TDirectoryFile("elas_vs_deltaE", "elas_vs_deltaE");
	dir_elas->cd();
	for(int ihist=0; ihist<m_n_cuts; ihist++) {
		h_elas_vs_deltaE[ihist]->Write();
	}
	dir_elas->cd("../");
	
	TDirectory *dir_dKvdE = new TDirectoryFile("deltaK_vs_deltaE", "deltaK_vs_deltaE");
	dir_dKvdE->cd();
	for(int ihist=0; ihist<m_n_cuts; ihist++) {
		h_deltaK_vs_deltaE[ihist]->Write();
	}
	dir_dKvdE->cd("../");
	
	TDirectory *dir_dCdE = new TDirectoryFile("deltaCCAL_vs_deltaE", "deltaCCAL_vs_deltaE");
	dir_dCdE->cd();
	for(int ihist=0; ihist<m_n_cuts; ihist++) {
		h_deltaCCAL_vs_deltaE[ihist]->Write();
		h_deltaCCAL_vs_deltaK[ihist]->Write();
	}
	dir_dCdE->cd("../");
	
	TDirectory *dir_dFdE = new TDirectoryFile("deltaFCAL_vs_deltaE", "deltaFCAL_vs_deltaE");
	dir_dFdE->cd();
	for(int ihist=0; ihist<m_n_cuts; ihist++) {
		h_deltaFCAL_vs_deltaE[ihist]->Write();
		h_deltaFCAL_vs_deltaK[ihist]->Write();
	}
	dir_dFdE->cd("../");
	
	TDirectory *dir_dFdC = new TDirectoryFile("deltaFCAL_vs_deltaCCAL", "deltaFCAL_vs_deltaCCAL");
	dir_dFdC->cd();
	for(int ihist=0; ihist<m_n_cuts; ihist++) {
		h_deltaFCAL_vs_deltaCCAL[ihist]->Write();
	}
	dir_dFdC->cd("../");
	
	TDirectory *dir_xy = new TDirectoryFile("xy", "xy");
	dir_xy->cd();
	for(int ihist=0; ihist<m_n_cuts; ihist++) {
		h_ccal_xy[ihist]->Write();
		h_fcal_xy[ihist]->Write();
	}
	dir_xy->cd("../");
	
	TDirectory *dir_opangle = new TDirectoryFile("opangle", "opangle");
	dir_opangle->cd();
	for(int ihist=0; ihist<m_n_cuts; ihist++) {
		h_opangle_tagh[ihist]->Write();
		h_opangle_tagh_ecut[ihist]->Write();
		h_opangle_tagh_ekcut[ihist]->Write();
		h_opangle_tagm[ihist]->Write();
		h_opangle_tagm_ecut[ihist]->Write();
		h_opangle_tagm_ekcut[ihist]->Write();
	}
	dir_opangle->cd("../");
	
	h_ccal_nblocks->Write();
	h_fcal_nblocks->Write();
	h_ccal_nblocks_cut->Write();
	h_fcal_nblocks_cut->Write();
	
	h_ccal_xy_highdeltaE->Write();
	h_fcal_xy_highdeltaE->Write();
	h_ccal_xy_lowdeltaE->Write();
	h_fcal_xy_lowdeltaE->Write();
	
	h_sumPhi_vs_deltaPhi->Write();
	
	fOut->Write();
	
	return;
}
