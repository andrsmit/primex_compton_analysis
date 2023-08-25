// $Id$
//
//    File: JEventProcessor_compton_simulation_detailed.h
// Created: Thu Feb 10 21:33:24 EST 2022
// Creator: andrsmit (on Linux ifarm1802.jlab.org 3.10.0-1160.11.1.el7.x86_64 x86_64)
//

#ifndef _JEventProcessor_compton_simulation_detailed_
#define _JEventProcessor_compton_simulation_detailed_

#include <JANA/JEventProcessor.h>
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>

#include "TTree.h"
#include "TH1.h"
#include "TH2.h"

using namespace jana;
using namespace std;

#include <CCAL/DCCALShower.h>
#include <CCAL/DCCALGeometry.h>
#include <FCAL/DFCALShower.h>
#include <FCAL/DFCALGeometry.h>
#include <PID/DBeamPhoton.h>
#include <PID/DEventRFBunch.h>

#include <TRACKING/DMCThrown.h>
#include <PID/DMCReaction.h>
#include <TRIGGER/DL1Trigger.h>
#include <TRIGGER/DTrigger.h>
#include <DAQ/Df250WindowRawData.h>
#include <DAQ/Df250PulseData.h>

#include <HDGEOMETRY/DGeometry.h>

#include "units.h"
#include "DLorentzVector.h"
#include "DVector3.h"
#include "TRandom3.h"

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <thread>
#include <mutex>

#include <CompCand.h>

class JEventProcessor_compton_simulation_detailed:public jana::JEventProcessor{
	public:
		JEventProcessor_compton_simulation_detailed(){};
		~JEventProcessor_compton_simulation_detailed(){};
		const char* className(void){return "JEventProcessor_compton_simulation_detailed";}

	private:
		jerror_t init(void);
		jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);
		jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);
		jerror_t erun(void);
		jerror_t fini(void);
		
		int fcal_fiducial_cut(DVector3 pos, DVector3 vertex);
		int ccal_fiducial_cut(DVector3 pos, DVector3 vertex);
		
		void fill_histograms(vector<ComptonCandidate_t> Comp_Candidates, DVector3 vertex, 
			bool is_double, vector<const DMCThrown*> mc_thrown);
		void set_cuts(int32_t runnumber);
		
		double get_event_weight(double vertex_z);
		
		TRandom3 *rndm_gen;
		
		//----------   Constants   ----------//
		
		double m_beamX, m_beamY, m_beamZ;
		double m_fcalX, m_fcalY, m_fcalZ;
		double m_ccalX, m_ccalY, m_ccalZ;
		
		double m_target_length;
		double m_target_density;
		double m_atten;
		
		double m_fcalX_new, m_fcalY_new;
		double m_ccalX_new, m_ccalY_new;
		
		DVector3 fcal_correction, ccal_correction;
		
		const double c    =  29.9792458;
		const double m_e  =  0.510998928e-3;
		
		const double FCAL_RF_time_cut = 3.0;
		const double CCAL_RF_time_cut = 3.0;
		const double BEAM_RF_time_cut = 2.004;
		
		const double FCAL_min_energy_cut = 0.35;
		const double CCAL_min_energy_cut = 3.0;
		const double BEAM_min_energy_cut = 3.0;
		
		int bfield_val = 0;
		int phase_val  = 0;
		
		//----------   Histograms   ----------//
		
		TH1F *h_theta_thrown_photon[2], *h_theta_thrown_electron[2], *h_theta_thrown_secondary[2];
		TH1F *h_phi12_thrown[2],        *h_phi13_thrown[3],          *h_phi23_thrown[3];
		
		TH2F *h_deltaE_vs_deltaK[3];
		TH1F *h_deltaE[3], *h_deltaPhi[3], *h_deltaK[3];
		
		TH1F *h_deltax_photon[2],   *h_deltay_photon[2],   *h_deltaE_photon[2];
		TH1F *h_deltax_electron[2], *h_deltay_electron[2], *h_deltaE_electron[2];
		
		TH1F *h_vertex,          *h_vertex_cut; 
		TH1F *h_vertex_weight,   *h_vertex_weight_cut;
		TH1F *h_vertex_accepted, *h_vertex_accepted_cut;
		
		TH1F *h_fcal_rf_dt, *h_ccal_rf_dt, *h_beam_rf_dt;
		TH1F *h_beam_rf_dt_cut;
		TH1F *h_beam;
		TH1F *h_tagh_flux, *h_tagm_flux;
		TH1F *h_beam_energy;
		
		TH1F *h_nccal, *h_nfcal, *h_n_ccal_fcal, *h_n_ccal_fcal_cut;
		TH2F *h_ccalE,     *h_fcalE;
		TH2F *h_ccalE_cut, *h_fcalE_cut;
		
		TH2F *h_deltaE_tagh,     *h_deltaE_smeared_tagh;
		TH2F *h_deltaE_tagm,     *h_deltaE_smeared_tagm;
		TH2F *h_deltaPhi_tagh,   *h_deltaPhi_smeared_tagh;
		TH2F *h_deltaPhi_tagm,   *h_deltaPhi_smeared_tagm;
		
		TH2F *h_deltaK_tagh,                      *h_deltaK_tagm;
		TH2F *h_deltaK_shifted_tagh,              *h_deltaK_shifted_tagm;
		TH2F *h_deltaK_smeared_tagh,              *h_deltaK_smeared_tagm;
		TH2F *h_deltaK_tagh_cut,                  *h_deltaK_tagm_cut;
		TH2F *h_deltaK_shifted_tagh_cut,          *h_deltaK_shifted_tagm_cut;
		TH2F *h_deltaK_smeared_tagh_cut,          *h_deltaK_smeared_tagm_cut;
		
		TH2F *h_deltaK_tagh_cut_main,             *h_deltaK_tagm_cut_main;
		
		TH2F *h_fcal_xy, *h_ccal_xy;
		
		//----------  Cuts  ---------//
		
		TF1 *f_deltaE_mu_mc,   *f_deltaE_sig_mc;
		TF1 *f_deltaE_mu_data, *f_deltaE_sig_data;
		
		double deltaE_mu_p0_mc,    deltaE_mu_p1_mc,   
		       deltaE_mu_p2_mc,    deltaE_mu_p3_mc;
		double deltaE_mu_p0_data,  deltaE_mu_p1_data, 
		       deltaE_mu_p2_data,  deltaE_mu_p3_data;
		double deltaE_sig_p0_mc,   deltaE_sig_p1_mc,   deltaE_sig_p2_mc;
		double deltaE_sig_p0_data, deltaE_sig_p1_data, deltaE_sig_p2_data;
		
		
		TF1 *f_deltaPhi_mu_mc,   *f_deltaPhi_sig_mc;
		TF1 *f_deltaPhi_mu_data, *f_deltaPhi_sig_data;
		
		double deltaPhi_mu_p0_mc,    deltaPhi_mu_p1_mc,    
		       deltaPhi_mu_p2_mc,    deltaPhi_mu_p3_mc;
		double deltaPhi_mu_p0_data,  deltaPhi_mu_p1_data,  
		       deltaPhi_mu_p2_data,  deltaPhi_mu_p3_data;
		double deltaPhi_sig_p0_mc,   deltaPhi_sig_p1_mc,   
		       deltaPhi_sig_p2_mc,   deltaPhi_sig_p3_mc;
		double deltaPhi_sig_p0_data, deltaPhi_sig_p1_data, 
		       deltaPhi_sig_p2_data, deltaPhi_sig_p3_data;
		
		
		TF1 *f_deltaK_mu_mc,   *f_deltaK_sig_mc;
		TF1 *f_deltaK_mu_data, *f_deltaK_sig_data;
		
		double deltaK_mu_p0_mc,    deltaK_mu_p1_mc,  
		       deltaK_mu_p2_mc,    deltaK_mu_p3_mc;
		double deltaK_mu_p0_data,  deltaK_mu_p1_data, 
		       deltaK_mu_p2_data,  deltaK_mu_p3_data;
		double deltaK_sig_p0_mc,   deltaK_sig_p1_mc, 
		       deltaK_sig_p2_mc,   deltaK_sig_p3_mc;
		double deltaK_sig_p0_data, deltaK_sig_p1_data, 
		       deltaK_sig_p2_data, deltaK_sig_p3_data;
};

#endif // _JEventProcessor_compton_simulation_detailed_

