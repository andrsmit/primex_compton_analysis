#include "ComptonAna.h"

void ComptonAna::runAnalysis_systematics(TString infname) {
	
	m_infile = new TFile(infname.Data(), "READ");
	int n_events_total = loadTree();
	
	while(m_event < n_events_total) {
		readEvent();
		comptonAnalysis_systematics();
		m_event++;
	}
	
	m_infile->Close();
	
	return;
}

void ComptonAna::comptonAnalysis_systematics() {
	
	if(m_nfcal > MAX_FCAL || m_nccal > MAX_CCAL || m_nbeam > MAX_BEAM) {
		cout << "Skipping event " << m_event << endl;
		return;
	}
	
	//---------------------------------------------------------------------------//
	// Make a list of good FCAL showers to use in analysis:
	
	vector<int> locGoodFCALShowers;
	locGoodFCALShowers.clear();
	
	vector<int> locNGoodFCALShowers_e_cut; locNGoodFCALShowers_e_cut.clear();
	for(int i=0; i<m_n_hists_fcalE; i++) locNGoodFCALShowers_e_cut.push_back(0.);
	
	vector<int> locNGoodFCALShowers_fid_cut; locNGoodFCALShowers_fid_cut.clear();
	for(int i=0; i<m_n_fid_cuts; i++) locNGoodFCALShowers_fid_cut.push_back(0.);
	
	vector<int> locNGoodFCALShowers_t_cut; locNGoodFCALShowers_t_cut.clear();
	for(int i=0; i<m_n_hists_fcalT; i++) locNGoodFCALShowers_t_cut.push_back(0.);
	
	int locNFCALShowers = 0, locNGoodFCALShowers = 0;
	for(int ishow=0; ishow<m_nfcal; ishow++) {
		
		TVector3 loc_pos = getFCALPosition(ishow);
		
		// to be consistent between Data and MC, we shoud remove showers from the dead region as soon as possible: 
		if(m_phase_val==1) {
			
			double loc_fcal_face_x = m_vertex.X() - m_fcal_face.X() 
				+ (loc_pos.X() * (m_fcal_face.Z() - m_vertex.Z())/loc_pos.Z());
			double loc_fcal_face_y = m_vertex.Y() - m_fcal_face.Y() 
				+ (loc_pos.Y() * (m_fcal_face.Z() - m_vertex.Z())/loc_pos.Z());
			
			if((-32. < loc_fcal_face_y && loc_fcal_face_y < -20.) 
				&& (-8. < loc_fcal_face_x && loc_fcal_face_x < 4.)) continue;
		}
		
		double loc_t = m_fcal_t[ishow] - (loc_pos.Mag()/m_c) - m_rfTime;
		
		if(fabs(loc_t) < 6.0) {
			locGoodFCALShowers.push_back(ishow);
		} else {
			continue;
		}
		
		// nominal FCAL Energy cut:
		int   e_cut_norm = m_fcal_e[ishow] > m_cut_fcalE ? 1 : 0;
		
		// nominal FCAL fiducial cut:
		int fid_cut_norm = fcal_fiducial_cut(loc_pos, 1.0);
		
		// nominal FCAL timing cut:
		int   t_cut_norm = fabs(loc_t) < m_cut_fcalrfdt ? 1 : 0;
		
		// loop over all minimum energy cuts and count number of showers 
		if(!fid_cut_norm && t_cut_norm) {
			for(int icut=0; icut<m_n_hists_fcalE; icut++) {
				double loc_e_cut = 0.05 * (double)icut;
				if(m_fcal_e[ishow] > loc_e_cut) locNGoodFCALShowers_e_cut[icut]++;
			}
		}
		
		// loop over all fiducial cuts and count number of showers:
		if(e_cut_norm && t_cut_norm) {
			for(int icut=0; icut<m_n_fid_cuts; icut++) {
				int loc_fid_cut = fcal_fiducial_cut(loc_pos, 0.5 + 0.1*(double)icut);
				if(!loc_fid_cut) locNGoodFCALShowers_fid_cut[icut]++;
			}
		}
		
		// loop over all timing cuts and count number of showers:
		if(e_cut_norm && !fid_cut_norm) {
			for(int icut=0; icut<m_n_hists_fcalT; icut++) {
				double loc_t_cut = 0.5*(double)(icut+1);
				if(fabs(loc_t) < loc_t_cut) locNGoodFCALShowers_t_cut[icut]++;
			}
		}
		
		if(!fid_cut_norm && e_cut_norm && t_cut_norm) locNGoodFCALShowers++;
	}
	
	//---------------------------------------------------------------------------//
	// Make a list of good CCAL showers to use in analysis:
	
	vector<int> locGoodCCALShowers;
	locGoodCCALShowers.clear();
	
	vector<int> locNGoodCCALShowers_e_cut; locNGoodCCALShowers_e_cut.clear();
	for(int i=0; i<m_n_hists_ccalE; i++) locNGoodCCALShowers_e_cut.push_back(0.);
	
	vector<int> locNGoodCCALShowers_fid_cut; locNGoodCCALShowers_fid_cut.clear();
	for(int i=0; i<m_n_fid_cuts; i++) locNGoodCCALShowers_fid_cut.push_back(0.);
	
	vector<int> locNGoodCCALShowers_t_cut; locNGoodCCALShowers_t_cut.clear();
	for(int i=0; i<m_n_hists_ccalT; i++) locNGoodCCALShowers_t_cut.push_back(0.);
	
	int locNCCALShowers = 0, locNGoodCCALShowers = 0;
	for(int ishow=0; ishow<m_nccal; ishow++) {
		
		TVector3 loc_pos = getCCALPosition(ishow);
		double loc_t = m_ccal_t[ishow] - (loc_pos.Mag()/m_c) - m_rfTime;
		
		if(fabs(loc_t) < 6.0) {
			locGoodCCALShowers.push_back(ishow);
		} else {
			continue;
		}
		
		// nominal CCAL Energy cut:
		int   e_cut_norm = m_ccal_e[ishow] > m_cut_ccalE ? 1 : 0;
		
		// nominal CCAL fiducial cut:
		int fid_cut_norm = ccal_fiducial_cut(loc_pos, 1.0);
		
		// nominal CCAL timing cut:
		int   t_cut_norm = fabs(loc_t) < m_cut_ccalrfdt ? 1 : 0;
		
		// loop over all minimum energy cuts and count number of showers 
		if(!fid_cut_norm && t_cut_norm) {
			for(int icut=0; icut<m_n_hists_ccalE; icut++) {
				double loc_e_cut = 0.5 * (double)icut;
				if(m_ccal_e[ishow] > loc_e_cut) locNGoodCCALShowers_e_cut[icut]++;
			}
		}
		
		// loop over all fiducial cuts and count number of showers:
		if(e_cut_norm && t_cut_norm) {
			for(int icut=0; icut<m_n_fid_cuts; icut++) {
				int loc_fid_cut = ccal_fiducial_cut(loc_pos, 0.1*(double)icut);
				if(!loc_fid_cut) locNGoodCCALShowers_fid_cut[icut]++;
			}
		}
		
		// loop over all timing cuts and count number of showers:
		if(e_cut_norm && !fid_cut_norm) {
			for(int icut=0; icut<m_n_hists_ccalT; icut++) {
				double loc_t_cut = 0.5*(double)(icut+1);
				if(fabs(loc_t) < loc_t_cut) locNGoodCCALShowers_t_cut[icut]++;
			}
		}
		
		if(!fid_cut_norm && e_cut_norm && t_cut_norm) locNGoodCCALShowers++;
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
		
		double loc_beam_cut = m_beam_bunches_main*m_cut_beamrfdt;
		
		if(fabs(loc_dt) < loc_beam_cut) loc_weight = 1.0;
		else if(
			((m_beam_bunches_main+5.5)*4.008 < fabs(loc_dt)) && 
			(fabs(loc_dt) < (m_beam_bunches_main+5.5+m_beam_bunches_acc)*4.008)
		) loc_weight = -1.0/(2.0*m_beam_bunches_acc);
		else continue;
		
		if(loc_weight < 0.0) loc_weight *= (m_acc_scale_factor[igam]-0.0);
		if(m_beam_e[igam] > 6.0) {
			locGoodBeamPhotons.push_back(igam);
			locGoodBeamPhotons_weight.push_back(loc_weight);
		}
	}
	
	//---------------------------------------------------------------------------//
	// Do full analysis:
	
	for(auto &ifcal : locGoodFCALShowers) {
		
		double e1 = m_fcal_e[ifcal];
		
		TVector3 pos1 = getFCALPosition(ifcal);
		double t1 = m_fcal_t[ifcal] - (pos1.Mag()/m_c) - m_rfTime;
		
		double theta1 = pos1.Theta();
		
		// split FCAL into 8 slices by azimuthal angle:
		
		double loc_fcal_face_x = m_vertex.X() + (pos1.X() * (m_fcal_face.Z() - m_vertex.Z())/pos1.Z()) - m_fcal_face.X();
		double loc_fcal_face_y = m_vertex.Y() + (pos1.Y() * (m_fcal_face.Z() - m_vertex.Z())/pos1.Z()) - m_fcal_face.Y();
		
		double phi1 = atan2(loc_fcal_face_y, loc_fcal_face_x) * (180./TMath::Pi());
		
		int fcal_phi_sect = (int)(phi1+180.)/45.0;
		if(     fcal_phi_sect<0) fcal_phi_sect = 0;
		else if(fcal_phi_sect>7) fcal_phi_sect = 7;
		
		// split FCAL based on layer:
		
		int fcal_layer = 8;
		for(int ilayer=1; ilayer<=7; ilayer++) {
			double loc_cut = (1.5 + (double)ilayer) * m_fcal_block_size;
			if((fabs(loc_fcal_face_x) < loc_cut) && (fabs(loc_fcal_face_y) < loc_cut)) {
				fcal_layer = ilayer;
				break;
			}
		}
		
		// lets place a 2degree cut on the angle of the FCAL shower:
		//if((theta1*(180./TMath::Pi())) < 2.0) continue;
		
		// check if FCAL shower passes nominal energy, timing, and fiducial cuts:
		
		int fcal_e_cut_norm   = e1 > m_cut_fcalE ? 1 : 0;
		int fcal_fid_cut_norm = fcal_fiducial_cut(pos1, 1.0);
		int fcal_t_cut_norm   = fabs(t1) < m_cut_fcalrfdt ? 1 : 0;
		
		//--------------------------------------------------------------------------------//
		
		for(auto &iccal : locGoodCCALShowers) {
			
			double e2 = m_ccal_e[iccal];
			
			TVector3 pos2 = getCCALPosition(iccal);
			double t2 = m_ccal_t[iccal] - (pos2.Mag()/m_c) - m_rfTime;
			
			double theta2 = pos2.Theta();
			
			// split CCAL into 8 slices by azimuthal angle:
			
			double loc_ccal_face_x = m_vertex.X() + (pos2.X() * (m_ccal_face.Z() - m_vertex.Z())/pos2.Z()) - m_ccal_face.X();
			double loc_ccal_face_y = m_vertex.Y() + (pos2.Y() * (m_ccal_face.Z() - m_vertex.Z())/pos2.Z()) - m_ccal_face.Y();
			
			double phi2 = atan2(loc_ccal_face_y, loc_ccal_face_x) * (180./TMath::Pi());
			
			int ccal_phi_sect = (int)(phi2+180.)/45.0;
			if(     ccal_phi_sect<0) ccal_phi_sect = 0;
			else if(ccal_phi_sect>7) ccal_phi_sect = 7;
			
			// split CCAL based on layer:
			
			int ccal_layer = 5;
			for(int ilayer=1; ilayer<=4; ilayer++) {
				double loc_cut = (1.0 + (double)ilayer) * m_ccal_block_size;
				if((fabs(loc_ccal_face_x) < loc_cut) && (fabs(loc_ccal_face_y) < loc_cut)) {
					ccal_layer = ilayer;
					break;
				}
			}
			
			// check if CCAL shower passes nominal energy, timing, and fiducial cuts:
			
			int ccal_e_cut_norm   = e2 > m_cut_ccalE ? 1 : 0;
			int ccal_fid_cut_norm = ccal_fiducial_cut(pos2, 1.0);
			int ccal_t_cut_norm   = fabs(t2) < m_cut_ccalrfdt ? 1 : 0;
			
			// calculate deltaPhi:
			
			double deltaPhi = fabs(pos2.Phi() - pos1.Phi()) * (180./TMath::Pi());
			
			// calculate invariant mass (assuming massless particles):
			
			double cos12   = (pos1.X()*pos2.X() + pos1.Y()*pos2.Y() + pos1.Z()*pos2.Z()) 
				/ (pos1.Mag()*pos2.Mag());
			double invmass = sqrt(2.0*e1*e2*(1.-cos12));
			
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
				
				int   deltaE_cut_norm = cut_deltaE(  deltaE,   eb, m_cut_deltaE,   1.e2);
				int deltaPhi_cut_norm = cut_deltaPhi(deltaPhi, eb, m_cut_deltaPhi, m_cut_deltaPhi);
				int   deltaK_cut_norm = cut_deltaK(  deltaK,   eb, m_cut_deltaK,   m_cut_deltaK);
				
				if(deltaPhi_cut_norm && 
					fcal_e_cut_norm && fcal_t_cut_norm && !fcal_fid_cut_norm && locNGoodFCALShowers==1 &&
					ccal_e_cut_norm && ccal_t_cut_norm && !ccal_fid_cut_norm && locNGoodCCALShowers==1) {
					
					h_mgg_vs_deltaK->Fill(deltaK, invmass, fill_weight);
					if(deltaE_cut_norm) {
						h_mgg_vs_deltaK_cut->Fill(deltaK, invmass, fill_weight);
					}
				}
				
				if(loc_tag_sys==0) {
					
					//
					// Nominal cuts:
					// Timing cuts on showers (already applied), Minimum energy cuts on FCAL and CCAL shower energies, 
					// fiducial cuts on FCAL and CCAL shower positions, DeltaE cut, DeltaPhi cut, Multiplicity Cut
					//
					
					if(
						deltaE_cut_norm &&  deltaPhi_cut_norm && 
						fcal_e_cut_norm && fcal_t_cut_norm && !fcal_fid_cut_norm && locNGoodFCALShowers==1 &&
						ccal_e_cut_norm && ccal_t_cut_norm && !ccal_fid_cut_norm && locNGoodCCALShowers==1)
					{
						h_deltaK_tagh_fcal_phi[fcal_phi_sect]->Fill(loc_tag_counter, deltaK, fill_weight);
						h_deltaK_tagh_ccal_phi[ccal_phi_sect]->Fill(loc_tag_counter, deltaK, fill_weight);
						h_deltaK_tagh_fcal_layer[fcal_layer-1]->Fill(loc_tag_counter, deltaK, fill_weight);
						h_deltaK_tagh_ccal_layer[ccal_layer-1]->Fill(loc_tag_counter, deltaK, fill_weight);
						
						if(deltaK_cut_norm) {
							h_xy_fcal_phi[fcal_phi_sect]->Fill(loc_fcal_face_x, loc_fcal_face_y, fill_weight);
							h_xy_ccal_phi[ccal_phi_sect]->Fill(loc_ccal_face_x, loc_ccal_face_y, fill_weight);
							h_xy_fcal_layer[fcal_layer-1]->Fill(loc_fcal_face_x, loc_fcal_face_y, fill_weight);
							h_xy_ccal_layer[ccal_layer-1]->Fill(loc_ccal_face_x, loc_ccal_face_y, fill_weight);
						}
					}
					
					
					// Vary minuum FCAL shower energy cut:
					
					if(
						deltaE_cut_norm &&  deltaPhi_cut_norm && 
						                   fcal_t_cut_norm && !fcal_fid_cut_norm && 
						ccal_e_cut_norm && ccal_t_cut_norm && !ccal_fid_cut_norm && locNGoodCCALShowers==1)
					{
						for(int icut=0; icut<m_n_hists_fcalE; icut++) {
							double loc_e_cut = 0.05 * (double)icut;
							if(e1 > loc_e_cut && locNGoodFCALShowers_e_cut[icut]==1) {
								h_deltaK_tagh_fcalE[icut]->Fill(loc_tag_counter, deltaK, fill_weight);
							}
						}
					}
					
					// Vary minuum CCAL shower energy cut:
					
					if(
						deltaE_cut_norm &&  deltaPhi_cut_norm && 
						fcal_e_cut_norm && fcal_t_cut_norm && !fcal_fid_cut_norm && locNGoodFCALShowers==1 &&
						                   ccal_t_cut_norm && !ccal_fid_cut_norm)
					{
						for(int icut=0; icut<m_n_hists_ccalE; icut++) {
							double loc_e_cut = 0.5 * (double)icut;
							if(e2 > loc_e_cut && locNGoodCCALShowers_e_cut[icut]==1) {
								h_deltaK_tagh_ccalE[icut]->Fill(loc_tag_counter, deltaK, fill_weight);
							}
						}
					}
					
					// Vary FCAL fiducial cut:
					
					if(
						deltaE_cut_norm &&  deltaPhi_cut_norm && 
						fcal_e_cut_norm && fcal_t_cut_norm && 
						ccal_e_cut_norm && ccal_t_cut_norm && !ccal_fid_cut_norm && locNGoodCCALShowers==1)
					{
						for(int icut=0; icut<m_n_fid_cuts; icut++) {
							int loc_fid_cut = fcal_fiducial_cut(pos1, 0.5 + 0.1*(double)icut);
							if(!loc_fid_cut && locNGoodFCALShowers_fid_cut[icut]==1) {
								h_deltaK_tagh_fcalfid[icut]->Fill(loc_tag_counter, deltaK, fill_weight);
							}
						}
					}
					
					// Vary CCAL fiducial cut:
					
					if(
						deltaE_cut_norm &&  deltaPhi_cut_norm && 
						fcal_e_cut_norm && fcal_t_cut_norm && !fcal_fid_cut_norm && locNGoodFCALShowers==1 &&
						ccal_e_cut_norm && ccal_t_cut_norm)
					{
						for(int icut=0; icut<m_n_fid_cuts; icut++) {
							int loc_fid_cut = ccal_fiducial_cut(pos2, 0.1*(double)icut);
							if(!loc_fid_cut && locNGoodCCALShowers_fid_cut[icut]==1) {
								h_deltaK_tagh_ccalfid[icut]->Fill(loc_tag_counter, deltaK, fill_weight);
							}
						}
					}
					
					// Vary minuum FCAL shower timing cut:
					
					if(
						deltaE_cut_norm &&  deltaPhi_cut_norm && 
						fcal_e_cut_norm &&                    !fcal_fid_cut_norm && 
						ccal_e_cut_norm && ccal_t_cut_norm && !ccal_fid_cut_norm && locNGoodCCALShowers==1)
					{
						for(int icut=0; icut<m_n_hists_fcalT; icut++) {
							double loc_t_cut = 0.5 * (double)(icut+1);
							if(fabs(t1) < loc_t_cut && locNGoodFCALShowers_t_cut[icut]==1) {
								h_deltaK_tagh_fcalT[icut]->Fill(loc_tag_counter, deltaK, fill_weight);
							}
						}
					}
					
					// Vary minuum CCAL shower timing cut:
					
					if(
						deltaE_cut_norm &&  deltaPhi_cut_norm && 
						fcal_e_cut_norm && fcal_t_cut_norm && !fcal_fid_cut_norm && locNGoodFCALShowers==1 &&
						ccal_t_cut_norm &&                    !ccal_fid_cut_norm)
					{
						for(int icut=0; icut<m_n_hists_ccalT; icut++) {
							double loc_t_cut = 0.5 * (double)(icut+1);
							if(fabs(t2) < loc_t_cut && locNGoodCCALShowers_t_cut[icut]==1) {
								h_deltaK_tagh_ccalT[icut]->Fill(loc_tag_counter, deltaK, fill_weight);
							}
						}
					}
					
					// Vary DeltaE cut:
					
					if(
						                    deltaPhi_cut_norm && 
						fcal_e_cut_norm && fcal_t_cut_norm && !fcal_fid_cut_norm && locNGoodFCALShowers==1 &&
						ccal_e_cut_norm && ccal_t_cut_norm && !ccal_fid_cut_norm && locNGoodCCALShowers==1)
					{
						for(int icut=0; icut<m_n_hists_deltaE; icut++) {
							if(cut_deltaE(deltaE, eb, m_deltaE_cuts[icut], 1.e2)) {
								h_deltaK_tagh_sigE[icut]->Fill(loc_tag_counter, deltaK, fill_weight);
							}
						}
					}
					
					// Vary DeltaPhi cut:
					
					if(
						deltaE_cut_norm && 
						fcal_e_cut_norm && fcal_t_cut_norm && !fcal_fid_cut_norm && locNGoodFCALShowers==1 &&
						ccal_e_cut_norm && ccal_t_cut_norm && !ccal_fid_cut_norm && locNGoodCCALShowers==1)
					{
						for(int icut=0; icut<m_n_hists_deltaPhi; icut++) {
							if(cut_deltaPhi(deltaPhi, eb, m_deltaPhi_cuts[icut], m_deltaPhi_cuts[icut])) {
								h_deltaK_tagh_sigPhi[icut]->Fill(loc_tag_counter, deltaK, fill_weight);
							}
						}
					}
					
					
				} else if(loc_tag_sys==1) {
					
					if(
						deltaE_cut_norm &&  deltaPhi_cut_norm && 
						fcal_e_cut_norm && fcal_t_cut_norm && !fcal_fid_cut_norm && locNGoodFCALShowers==1 &&
						ccal_e_cut_norm && ccal_t_cut_norm && !ccal_fid_cut_norm && locNGoodCCALShowers==1)
					{
						h_deltaK_tagm_fcal_phi[fcal_phi_sect]->Fill(loc_tag_counter, deltaK, fill_weight);
						h_deltaK_tagm_ccal_phi[ccal_phi_sect]->Fill(loc_tag_counter, deltaK, fill_weight);
						h_deltaK_tagm_fcal_layer[fcal_layer-1]->Fill(loc_tag_counter, deltaK, fill_weight);
						h_deltaK_tagm_ccal_layer[ccal_layer-1]->Fill(loc_tag_counter, deltaK, fill_weight);
						
						if(deltaK_cut_norm) {
							h_xy_fcal_phi[fcal_phi_sect]->Fill(loc_fcal_face_x, loc_fcal_face_y, fill_weight);
							h_xy_ccal_phi[ccal_phi_sect]->Fill(loc_ccal_face_x, loc_ccal_face_y, fill_weight);
							h_xy_fcal_layer[fcal_layer-1]->Fill(loc_fcal_face_x, loc_fcal_face_y, fill_weight);
							h_xy_ccal_layer[ccal_layer-1]->Fill(loc_ccal_face_x, loc_ccal_face_y, fill_weight);
						}
					}
					
					// Vary minuum FCAL shower energy cut:
					
					if(
						deltaE_cut_norm &&  deltaPhi_cut_norm && 
						                   fcal_t_cut_norm && !fcal_fid_cut_norm && 
						ccal_e_cut_norm && ccal_t_cut_norm && !ccal_fid_cut_norm && locNGoodCCALShowers==1)
					{
						for(int icut=0; icut<m_n_hists_fcalE; icut++) {
							double loc_e_cut = 0.05 * (double)icut;
							if(e1 > loc_e_cut && locNGoodFCALShowers_e_cut[icut]==1) {
								h_deltaK_tagm_fcalE[icut]->Fill(loc_tag_counter, deltaK, fill_weight);
							}
						}
					}
					
					// Vary minuum CCAL shower energy cut:
					
					if(
						deltaE_cut_norm &&  deltaPhi_cut_norm && 
						fcal_e_cut_norm && fcal_t_cut_norm && !fcal_fid_cut_norm && locNGoodFCALShowers==1 &&
						                   ccal_t_cut_norm && !ccal_fid_cut_norm)
					{
						for(int icut=0; icut<m_n_hists_ccalE; icut++) {
							double loc_e_cut = 0.5 * (double)icut;
							if(e2 > loc_e_cut && locNGoodCCALShowers_e_cut[icut]==1) {
								h_deltaK_tagm_ccalE[icut]->Fill(loc_tag_counter, deltaK, fill_weight);
							}
						}
					}
					
					// Vary FCAL fiducial cut:
					
					if(
						deltaE_cut_norm &&  deltaPhi_cut_norm && 
						fcal_e_cut_norm && fcal_t_cut_norm && 
						ccal_e_cut_norm && fcal_t_cut_norm && !ccal_fid_cut_norm && locNGoodCCALShowers==1)
					{
						for(int icut=0; icut<m_n_fid_cuts; icut++) {
							int loc_fid_cut = fcal_fiducial_cut(pos1, 0.5 + 0.1*(double)icut);
							if(!loc_fid_cut && locNGoodFCALShowers_fid_cut[icut]==1) {
								h_deltaK_tagm_fcalfid[icut]->Fill(loc_tag_counter, deltaK, fill_weight);
							}
						}
					}
					
					// Vary CCAL fiducial cut:
					
					if(
						deltaE_cut_norm &&  deltaPhi_cut_norm && 
						fcal_e_cut_norm && fcal_t_cut_norm && !fcal_fid_cut_norm && locNGoodFCALShowers==1 &&
						ccal_e_cut_norm && ccal_t_cut_norm)
					{
						for(int icut=0; icut<m_n_fid_cuts; icut++) {
							int loc_fid_cut = ccal_fiducial_cut(pos2, 0.1*(double)icut);
							if(!loc_fid_cut && locNGoodCCALShowers_fid_cut[icut]==1) {
								h_deltaK_tagm_ccalfid[icut]->Fill(loc_tag_counter, deltaK, fill_weight);
							}
						}
					}
					
					// Vary minuum FCAL shower timing cut:
					
					if(
						deltaE_cut_norm &&  deltaPhi_cut_norm && 
						fcal_e_cut_norm &&                    !fcal_fid_cut_norm && 
						ccal_e_cut_norm && ccal_t_cut_norm && !ccal_fid_cut_norm && locNGoodCCALShowers==1)
					{
						for(int icut=0; icut<m_n_hists_fcalT; icut++) {
							double loc_t_cut = 0.5 * (double)(icut+1);
							if(fabs(t1) < loc_t_cut && locNGoodFCALShowers_t_cut[icut]==1) {
								h_deltaK_tagm_fcalT[icut]->Fill(loc_tag_counter, deltaK, fill_weight);
							}
						}
					}
					
					// Vary minuum CCAL shower timing cut:
					
					if(
						deltaE_cut_norm &&  deltaPhi_cut_norm && 
						fcal_e_cut_norm && fcal_t_cut_norm && !fcal_fid_cut_norm && locNGoodFCALShowers==1 &&
						ccal_t_cut_norm &&                    !ccal_fid_cut_norm)
					{
						for(int icut=0; icut<m_n_hists_ccalT; icut++) {
							double loc_t_cut = 0.5 * (double)(icut+1);
							if(fabs(t2) < loc_t_cut && locNGoodCCALShowers_t_cut[icut]==1) {
								h_deltaK_tagm_ccalT[icut]->Fill(loc_tag_counter, deltaK, fill_weight);
							}
						}
					}
					
					// Vary DeltaE cut:
					
					if(
						                    deltaPhi_cut_norm && 
						fcal_e_cut_norm && fcal_t_cut_norm && !fcal_fid_cut_norm && locNGoodFCALShowers==1 &&
						ccal_e_cut_norm && ccal_t_cut_norm && !ccal_fid_cut_norm && locNGoodCCALShowers==1)
					{
						for(int icut=0; icut<m_n_hists_deltaE; icut++) {
							if(cut_deltaE(deltaE, eb, m_deltaE_cuts[icut], 1.e2)) {
								h_deltaK_tagm_sigE[icut]->Fill(loc_tag_counter, deltaK, fill_weight);
							}
						}
					}
					
					// Vary DeltaPhi cut:
					
					if(
						deltaE_cut_norm && 
						fcal_e_cut_norm && fcal_t_cut_norm && !fcal_fid_cut_norm && locNGoodFCALShowers==1 &&
						ccal_e_cut_norm && ccal_t_cut_norm && !ccal_fid_cut_norm && locNGoodCCALShowers==1)
					{
						for(int icut=0; icut<m_n_hists_deltaPhi; icut++) {
							if(cut_deltaPhi(deltaPhi, eb, m_deltaPhi_cuts[icut], m_deltaPhi_cuts[icut])) {
								h_deltaK_tagm_sigPhi[icut]->Fill(loc_tag_counter, deltaK, fill_weight);
							}
						}
					}
				}
				
			} // end loop over Beam photons
		} // end loop over CCAL shower
	} // end loop over FCAL showers
	
	return;
}


void ComptonAna::initHistograms_systematics() {
	
	h_mgg_vs_deltaK     = new TH2F("mgg_vs_deltaK",     "Inv. Mass vs. #DeltaK; #DeltaK [GeV]; m_{#gamma#gamma} [GeV/c^{2}]",
		1000, -8.0, 8.0, 1000, 0., 1.0);
	h_mgg_vs_deltaK_cut = new TH2F("mgg_vs_deltaK_cut", "Inv. Mass vs. #DeltaK; #DeltaK [GeV]; m_{#gamma#gamma} [GeV/c^{2}]",
		1000, -8.0, 8.0, 1000, 0., 1.0);
	
	//-----------------------------------------//
	// Vary minimum FCAL shower energy cut:
	
	for(int ihist=0; ihist<m_n_hists_fcalE; ihist++) {
		int loc_hist_index = 5*ihist;
		double loc_cut_val = 0.01 * (double)(loc_hist_index);
		
		h_deltaK_tagh_fcalE[ihist] = new TH2F(Form("deltaK_tagh_%02dfcalE", loc_hist_index), 
			Form("#DeltaK (E_{FCAL} > %.2f GeV); TAGH Counter; E_{Compton} - E_{#gamma} [GeV]", loc_cut_val), 
			274, 0.5, 274.5, 1000, -8.0, 8.0);
		h_deltaK_tagh_fcalE[ihist]->Sumw2();
		
		h_deltaK_tagm_fcalE[ihist] = new TH2F(Form("deltaK_tagm_%02dfcalE", loc_hist_index), 
			Form("#DeltaK (E_{FCAL} > %.2f GeV); TAGM Counter; E_{Compton} - E_{#gamma} [GeV]", loc_cut_val), 
			102, 0.5, 102.5, 1000, -8.0, 8.0);
		h_deltaK_tagm_fcalE[ihist]->Sumw2();
	}
	
	//-----------------------------------------//
	// Vary minimum CCAL shower energy cut:
	
	for(int ihist=0; ihist<m_n_hists_ccalE; ihist++) {
		int loc_hist_index = 5*ihist;
		double loc_cut_val = 0.1 * (double)(loc_hist_index);
		
		h_deltaK_tagh_ccalE[ihist] = new TH2F(Form("deltaK_tagh_%02dccalE", loc_hist_index), 
			Form("#DeltaK (E_{CCAL} > %.2f GeV); TAGH Counter; E_{Compton} - E_{#gamma} [GeV]", loc_cut_val), 
			274, 0.5, 274.5, 1000, -8.0, 8.0);
		h_deltaK_tagh_ccalE[ihist]->Sumw2();
		
		h_deltaK_tagm_ccalE[ihist] = new TH2F(Form("deltaK_tagm_%02dccalE", loc_hist_index), 
			Form("#DeltaK (E_{CCAL} > %.2f GeV); TAGM Counter; E_{Compton} - E_{#gamma} [GeV]", loc_cut_val), 
			102, 0.5, 102.5, 1000, -8.0, 8.0);
		h_deltaK_tagm_ccalE[ihist]->Sumw2();
	}
	
	//-----------------------------------------//
	// Vary minimum FCAL shower timing cut:
	
	for(int ihist=0; ihist<m_n_hists_fcalT; ihist++) {
		double loc_cut_val = 0.5 * (double)(ihist+1);
		int loc_hist_index = (int)(10.*loc_cut_val);
		
		h_deltaK_tagh_fcalT[ihist] = new TH2F(Form("deltaK_tagh_%02dfcalT", loc_hist_index), 
			Form("#DeltaK (|t_{FCAL}-t_{RF}| < %.2f ns); TAGH Counter; E_{Compton} - E_{#gamma} [GeV]", loc_cut_val), 
			274, 0.5, 274.5, 1000, -8.0, 8.0);
		h_deltaK_tagh_fcalT[ihist]->Sumw2();
		
		h_deltaK_tagm_fcalT[ihist] = new TH2F(Form("deltaK_tagm_%02dfcalT", loc_hist_index), 
			Form("#DeltaK (|t_{FCAL}-t_{RF}| < %.2f ns); TAGM Counter; E_{Compton} - E_{#gamma} [GeV]", loc_cut_val), 
			102, 0.5, 102.5, 1000, -8.0, 8.0);
		h_deltaK_tagm_fcalT[ihist]->Sumw2();
	}
	
	//-----------------------------------------//
	// Vary minimum CCAL shower timing cut:
	
	for(int ihist=0; ihist<m_n_hists_ccalT; ihist++) {
		double loc_cut_val = 0.5 * (double)(ihist+1);
		int loc_hist_index = (int)(10.*loc_cut_val);
		
		h_deltaK_tagh_ccalT[ihist] = new TH2F(Form("deltaK_tagh_%02dccalT", loc_hist_index), 
			Form("#DeltaK (|t_{CCAL}-t_{RF}| < %.2f ns); TAGH Counter; E_{Compton} - E_{#gamma} [GeV]", loc_cut_val), 
			274, 0.5, 274.5, 1000, -8.0, 8.0);
		h_deltaK_tagh_ccalT[ihist]->Sumw2();
		
		h_deltaK_tagm_ccalT[ihist] = new TH2F(Form("deltaK_tagm_%02dccalT", loc_hist_index), 
			Form("#DeltaK (|t_{CCAL}-t_{RF}| < %.2f ns); TAGM Counter; E_{Compton} - E_{#gamma} [GeV]", loc_cut_val), 
			102, 0.5, 102.5, 1000, -8.0, 8.0);
		h_deltaK_tagm_ccalT[ihist]->Sumw2();
	}
	
	//-----------------------------------------//
	// Vary size of square fiducial beam-hole cuts:
	
	for(int icut=0; icut<m_n_fid_cuts; icut++) {
		
		h_deltaK_tagh_fcalfid[icut] = new TH2F(Form("deltaK_tagh_fcalfid_%02d",icut), 
			"#DeltaK", 274, 0.5, 274.5, 1000, -8.0, 8.0);
		h_deltaK_tagm_fcalfid[icut] = new TH2F(Form("deltaK_tagm_fcalfid_%02d",icut), 
			"#DeltaK", 102, 0.5, 102.5, 1000, -8.0, 8.0);
		h_deltaK_tagh_fcalfid[icut]->Sumw2();
		h_deltaK_tagm_fcalfid[icut]->Sumw2();
	}
	
	for(int icut=0; icut<m_n_fid_cuts; icut++) {
		
		h_deltaK_tagh_ccalfid[icut] = new TH2F(Form("deltaK_tagh_ccalfid_%02d",icut), 
			"#DeltaK", 274, 0.5, 274.5, 1000, -8.0, 8.0);
		h_deltaK_tagm_ccalfid[icut] = new TH2F(Form("deltaK_tagm_ccalfid_%02d",icut), 
			"#DeltaK", 102, 0.5, 102.5, 1000, -8.0, 8.0);
		h_deltaK_tagh_ccalfid[icut]->Sumw2();
		h_deltaK_tagm_ccalfid[icut]->Sumw2();
	}
	
	//-----------------------------------------//
	// Vary width of DeltaE cut:
	
	for(int ihist=0; ihist<m_n_hists_deltaE; ihist++) {
		int loc_hist_index = (int)(10.*m_deltaE_cuts[ihist]);
		double loc_cut_val = m_deltaE_cuts[ihist];
		
		h_deltaK_tagh_sigE[ihist] = new TH2F(Form("deltaK_tagh_%03dsigE", loc_hist_index), 
			Form("#DeltaK (|#DeltaE| < %.1f#sigma); TAGH Counter; E_{Compton} - E_{#gamma} [GeV]", loc_cut_val), 
			274, 0.5, 274.5, 1000, -8.0, 8.0);
		h_deltaK_tagh_sigE[ihist]->Sumw2();
		
		h_deltaK_tagm_sigE[ihist] = new TH2F(Form("deltaK_tagm_%03dsigE", loc_hist_index), 
			Form("#DeltaK (|#DeltaE| < %.1f#sigma); TAGM Counter; E_{Compton} - E_{#gamma} [GeV]", loc_cut_val), 
			102, 0.5, 102.5, 1000, -8.0, 8.0);
		h_deltaK_tagm_sigE[ihist]->Sumw2();
	}
	
	//-----------------------------------------//
	// Vary width of DeltaPhi cut:
	
	for(int ihist=0; ihist<m_n_hists_deltaPhi; ihist++) {
		int loc_hist_index = (int)(10.*m_deltaPhi_cuts[ihist]);
		double loc_cut_val = m_deltaPhi_cuts[ihist];
		
		h_deltaK_tagh_sigPhi[ihist] = new TH2F(Form("deltaK_tagh_%03dsigPhi", loc_hist_index), 
			Form("#DeltaK (|#Delta#phi| < %.1f#sigma); TAGH Counter; E_{Compton} - E_{#gamma} [GeV]", loc_cut_val), 
			274, 0.5, 274.5, 1000, -8.0, 8.0);
		h_deltaK_tagh_sigPhi[ihist]->Sumw2();
		
		h_deltaK_tagm_sigPhi[ihist] = new TH2F(Form("deltaK_tagm_%03dsigPhi", loc_hist_index), 
			Form("#DeltaK (|#Delta#phi| < %.1f#sigma); TAGM Counter; E_{Compton} - E_{#gamma} [GeV]", loc_cut_val), 
			102, 0.5, 102.5, 1000, -8.0, 8.0);
		h_deltaK_tagm_sigPhi[ihist]->Sumw2();
	}
	
	//-----------------------------------------//
	// Measure cross-section in different regions of FCAL and CCAL:
	
	for(int iphi=0; iphi<8; iphi++) {
		
		h_deltaK_tagh_fcal_phi[iphi] = new TH2F(Form("deltaK_tagh_fcal_phi_%d",iphi), 
			Form("#DeltaK (FCAL Phi Sect. %d)",iphi), 274, 0.5, 274.5, 1000, -8.0, 8.0);
		h_deltaK_tagm_fcal_phi[iphi] = new TH2F(Form("deltaK_tagm_fcal_phi_%d",iphi), 
			Form("#DeltaK (FCAL Phi Sect. %d)",iphi), 102, 0.5, 102.5, 1000, -8.0, 8.0);
		
		h_deltaK_tagh_fcal_phi[iphi]->Sumw2();
		h_deltaK_tagm_fcal_phi[iphi]->Sumw2();
		
		h_xy_fcal_phi[iphi] = new TH2F(Form("xy_fcal_phi_%d",iphi), 
			Form("FCAL Y vs. X (Phi Sect. %d)",iphi), 500, -100., 100., 500, -100., 100.);
	}
	
	for(int iphi=0; iphi<8; iphi++) {
		
		h_deltaK_tagh_ccal_phi[iphi] = new TH2F(Form("deltaK_tagh_ccal_phi_%d",iphi), 
			Form("#DeltaK (CCAL Phi Sect. %d)",iphi), 274, 0.5, 274.5, 1000, -8.0, 8.0);
		h_deltaK_tagm_ccal_phi[iphi] = new TH2F(Form("deltaK_tagm_ccal_phi_%d",iphi), 
			Form("#DeltaK (CCAL Phi Sect. %d)",iphi), 102, 0.5, 102.5, 1000, -8.0, 8.0);
		
		h_deltaK_tagh_ccal_phi[iphi]->Sumw2();
		h_deltaK_tagm_ccal_phi[iphi]->Sumw2();
		
		h_xy_ccal_phi[iphi] = new TH2F(Form("xy_ccal_phi_%d",iphi), 
			Form("CCAL Y vs. X (Phi Sect. %d)",iphi), 500,  -13.,  13., 500,  -13.,  13.);
	}
	
	for(int ilayer=1; ilayer<=8; ilayer++) {
		
		h_deltaK_tagh_fcal_layer[ilayer-1] = new TH2F(Form("deltaK_tagh_fcal_layer_%d",ilayer), 
			Form("#DeltaK (FCAL Layer %d)",ilayer), 274, 0.5, 274.5, 1000, -8.0, 8.0);
		h_deltaK_tagm_fcal_layer[ilayer-1] = new TH2F(Form("deltaK_tagm_fcal_layer_%d",ilayer), 
			Form("#DeltaK (FCAL Layer %d)",ilayer), 102, 0.5, 102.5, 1000, -8.0, 8.0);
		
		h_deltaK_tagh_fcal_layer[ilayer-1]->Sumw2();
		h_deltaK_tagm_fcal_layer[ilayer-1]->Sumw2();
		
		h_xy_fcal_layer[ilayer-1] = new TH2F(Form("xy_fcal_layer_%d",ilayer), 
			Form("FCAL Y vs. X (Layer %d)",ilayer), 500, -100., 100., 500, -100., 100.);
	}
	
	for(int ilayer=1; ilayer<=5; ilayer++) {
		
		h_deltaK_tagh_ccal_layer[ilayer-1] = new TH2F(Form("deltaK_tagh_ccal_layer_%d",ilayer), 
			Form("#DeltaK (CCAL Layer %d)",ilayer), 274, 0.5, 274.5, 1000, -8.0, 8.0);
		h_deltaK_tagm_ccal_layer[ilayer-1] = new TH2F(Form("deltaK_tagm_ccal_layer_%d",ilayer), 
			Form("#DeltaK (CCAL Layer %d)",ilayer), 102, 0.5, 102.5, 1000, -8.0, 8.0);
		
		h_deltaK_tagh_ccal_layer[ilayer-1]->Sumw2();
		h_deltaK_tagm_ccal_layer[ilayer-1]->Sumw2();
		
		h_xy_ccal_layer[ilayer-1] = new TH2F(Form("xy_ccal_layer_%d",ilayer), 
			Form("CCAL Y vs. X (Layer %d)",ilayer), 500,  -13.,  13., 500,  -13.,  13.);
	}
	
	return;
}

void ComptonAna::resetHistograms_systematics() {
	
	h_mgg_vs_deltaK->Reset();
	h_mgg_vs_deltaK_cut->Reset();
	
	//-----------------------------------------//
	// Vary minimum FCAL shower energy cut:
	
	for(int ihist=0; ihist<m_n_hists_fcalE; ihist++) {
		h_deltaK_tagh_fcalE[ihist]->Reset();
		h_deltaK_tagm_fcalE[ihist]->Reset();
	}
	
	//-----------------------------------------//
	// Vary minimum CCAL shower energy cut:
	
	for(int ihist=0; ihist<m_n_hists_ccalE; ihist++) {
		h_deltaK_tagh_ccalE[ihist]->Reset();
		h_deltaK_tagm_ccalE[ihist]->Reset();
	}
	
	//-----------------------------------------//
	// Vary minimum FCAL shower timing cut:
	
	for(int ihist=0; ihist<m_n_hists_fcalT; ihist++) {
		h_deltaK_tagh_fcalT[ihist]->Reset();
		h_deltaK_tagm_fcalT[ihist]->Reset();
	}
	
	//-----------------------------------------//
	// Vary minimum CCAL shower timing cut:
	
	for(int ihist=0; ihist<m_n_hists_ccalT; ihist++) {
		h_deltaK_tagh_ccalT[ihist]->Reset();
		h_deltaK_tagm_ccalT[ihist]->Reset();
	}
	
	//-----------------------------------------//
	// Vary size of square fiducial beam-hole cuts:
	
	for(int icut=0; icut<m_n_fid_cuts; icut++) {
		h_deltaK_tagh_fcalfid[icut]->Reset();
		h_deltaK_tagm_fcalfid[icut]->Reset();
	}
	
	for(int icut=0; icut<m_n_fid_cuts; icut++) {
		h_deltaK_tagh_ccalfid[icut]->Reset();
		h_deltaK_tagm_ccalfid[icut]->Reset();
	}
	
	//-----------------------------------------//
	// Vary width of DeltaE cut:
	
	for(int ihist=0; ihist<m_n_hists_deltaE; ihist++) {
		h_deltaK_tagh_sigE[ihist]->Reset();
		h_deltaK_tagm_sigE[ihist]->Reset();
	}
	
	//-----------------------------------------//
	// Vary width of DeltaPhi cut:
	
	for(int ihist=0; ihist<m_n_hists_deltaPhi; ihist++) {
		h_deltaK_tagh_sigPhi[ihist]->Reset();
		h_deltaK_tagm_sigPhi[ihist]->Reset();
	}
	
	//-----------------------------------------//
	// Measure cross-section in different regions of FCAL and CCAL:
	
	for(int iphi=0; iphi<8; iphi++) {
		h_deltaK_tagh_fcal_phi[iphi]->Reset();
		h_deltaK_tagm_fcal_phi[iphi]->Reset();
		h_xy_fcal_phi[iphi]->Reset();
	}
	
	for(int iphi=0; iphi<8; iphi++) {
		h_deltaK_tagh_ccal_phi[iphi]->Reset();
		h_deltaK_tagm_ccal_phi[iphi]->Reset();
		h_xy_ccal_phi[iphi]->Reset();
	}
	
	for(int ilayer=1; ilayer<=8; ilayer++) {
		h_deltaK_tagh_fcal_layer[ilayer-1]->Reset();
		h_deltaK_tagm_fcal_layer[ilayer-1]->Reset();
		h_xy_fcal_layer[ilayer-1]->Reset();
	}
	
	for(int ilayer=1; ilayer<=5; ilayer++) {
		h_deltaK_tagh_ccal_layer[ilayer-1]->Reset();
		h_deltaK_tagm_ccal_layer[ilayer-1]->Reset();
		h_xy_ccal_layer[ilayer-1]->Reset();
	}
	
	return;
}

void ComptonAna::writeHistograms_systematics() {
	
	cout << "writing histograms to: " << m_output_fname << "..." << endl;
	
	TFile *fOut = new TFile(m_output_fname.c_str(), "RECREATE");
	fOut->cd();
	
	h_mgg_vs_deltaK->Write();
	h_mgg_vs_deltaK_cut->Write();
	
	//-----------------------------------------//
	// Vary minimum FCAL shower energy cut:
	
	TDirectory *dir_fcalE = new TDirectoryFile("fcalE", "fcalE");
	dir_fcalE->cd();
	for(int ihist=0; ihist<m_n_hists_fcalE; ihist++) {
		h_deltaK_tagh_fcalE[ihist]->Write();
		h_deltaK_tagm_fcalE[ihist]->Write();
	}
	dir_fcalE->cd("../");
	
	//-----------------------------------------//
	// Vary minimum CCAL shower energy cut:
	
	TDirectory *dir_ccalE = new TDirectoryFile("ccalE", "ccalE");
	dir_ccalE->cd();
	for(int ihist=0; ihist<m_n_hists_ccalE; ihist++) {
		h_deltaK_tagh_ccalE[ihist]->Write();
		h_deltaK_tagm_ccalE[ihist]->Write();
	}
	dir_ccalE->cd("../");
	
	//-----------------------------------------//
	// Vary minimum FCAL shower timing cut:
	
	TDirectory *dir_fcalT = new TDirectoryFile("fcalT", "fcalT");
	dir_fcalT->cd();
	for(int ihist=0; ihist<m_n_hists_fcalT; ihist++) {
		h_deltaK_tagh_fcalT[ihist]->Write();
		h_deltaK_tagm_fcalT[ihist]->Write();
	}
	dir_fcalT->cd("../");
	
	//-----------------------------------------//
	// Vary minimum CCAL shower timing cut:
	
	TDirectory *dir_ccalT = new TDirectoryFile("ccalT", "ccalT");
	dir_ccalT->cd();
	for(int ihist=0; ihist<m_n_hists_ccalT; ihist++) {
		h_deltaK_tagh_ccalT[ihist]->Write();
		h_deltaK_tagm_ccalT[ihist]->Write();
	}
	dir_ccalT->cd("../");
	
	//-----------------------------------------//
	// Vary size of square fiducial beam-hole cuts:
	
	TDirectory *dir_fcalfid = new TDirectoryFile("fcalfid", "fcalfid");
	dir_fcalfid->cd();
	for(int icut=0; icut<m_n_fid_cuts; icut++) {
		h_deltaK_tagh_fcalfid[icut]->Write();
		h_deltaK_tagm_fcalfid[icut]->Write();
	}
	dir_fcalfid->cd("../");
	
	TDirectory *dir_ccalfid = new TDirectoryFile("ccalfid", "ccalfid");
	dir_ccalfid->cd();
	for(int icut=0; icut<m_n_fid_cuts; icut++) {
		h_deltaK_tagh_ccalfid[icut]->Write();
		h_deltaK_tagm_ccalfid[icut]->Write();
	}
	dir_ccalfid->cd("../");
	
	//-----------------------------------------//
	// Vary width of DeltaE cut:
	
	TDirectory *dir_deltaE = new TDirectoryFile("DeltaE", "DeltaE");
	dir_deltaE->cd();
	for(int ihist=0; ihist<m_n_hists_deltaE; ihist++) {
		h_deltaK_tagh_sigE[ihist]->Write();
		h_deltaK_tagm_sigE[ihist]->Write();
	}
	dir_deltaE->cd("../");
	
	//-----------------------------------------//
	// Vary width of DeltaPhi cut:
	
	TDirectory *dir_deltaPhi = new TDirectoryFile("DeltaPhi", "DeltaPhi");
	dir_deltaPhi->cd();
	for(int ihist=0; ihist<m_n_hists_deltaPhi; ihist++) {
		h_deltaK_tagh_sigPhi[ihist]->Write();
		h_deltaK_tagm_sigPhi[ihist]->Write();
	}
	dir_deltaPhi->cd("../");
	
	//-----------------------------------------//
	// Measure cross-section in different regions of FCAL and CCAL:
	
	TDirectory *dir_fcal_phi = new TDirectoryFile("fcal_phi", "fcal_phi");
	dir_fcal_phi->cd();
	for(int iphi=0; iphi<8; iphi++) {
		h_deltaK_tagh_fcal_phi[iphi]->Write();
		h_deltaK_tagm_fcal_phi[iphi]->Write();
		h_xy_fcal_phi[iphi]->Write();
	}
	dir_fcal_phi->cd("../");
	
	TDirectory *dir_ccal_phi = new TDirectoryFile("ccal_phi", "ccal_phi");
	dir_ccal_phi->cd();
	for(int iphi=0; iphi<8; iphi++) {
		h_deltaK_tagh_ccal_phi[iphi]->Write();
		h_deltaK_tagm_ccal_phi[iphi]->Write();
		h_xy_ccal_phi[iphi]->Write();
	}
	dir_ccal_phi->cd("../");
	
	TDirectory *dir_fcal_layer = new TDirectoryFile("fcal_layer", "fcal_layer");
	dir_fcal_layer->cd();
	for(int ilayer=1; ilayer<=8; ilayer++) {
		h_deltaK_tagh_fcal_layer[ilayer-1]->Write();
		h_deltaK_tagm_fcal_layer[ilayer-1]->Write();
		h_xy_fcal_layer[ilayer-1]->Write();
	}
	dir_fcal_layer->cd("../");
	
	TDirectory *dir_ccal_layer = new TDirectoryFile("ccal_layer", "ccal_layer");
	dir_ccal_layer->cd();
	for(int ilayer=1; ilayer<=5; ilayer++) {
		h_deltaK_tagh_ccal_layer[ilayer-1]->Write();
		h_deltaK_tagm_ccal_layer[ilayer-1]->Write();
		h_xy_ccal_layer[ilayer-1]->Write();
	}
	dir_ccal_layer->cd("../");
	
	fOut->Write();
	
	return;
}
