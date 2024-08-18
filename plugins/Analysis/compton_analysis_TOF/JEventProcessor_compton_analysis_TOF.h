// $Id$
//
//    File: JEventProcessor_compton_analysis_TOF.h
// Created: Wed May 22 03:17:27 PM EDT 2024
// Creator: andrsmit (on Linux ifarm180302.jlab.org 5.14.0-362.13.1.el9_3.x86_64 x86_64)
//

#ifndef _JEventProcessor_compton_analysis_TOF_
#define _JEventProcessor_compton_analysis_TOF_

#include <JANA/JEventProcessor.h>
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>

#include "TTree.h"
#include "TH1.h"
#include "TH2.h"

using namespace jana;
using namespace std;

#include <CCAL/DCCALShower.h>
#include <FCAL/DFCALShower.h>
#include <PID/DBeamPhoton.h>
#include <PID/DEventRFBunch.h>

#include <TOF/DTOFPoint.h>

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

class JEventProcessor_compton_analysis_TOF:public jana::JEventProcessor{
	
	public:
		JEventProcessor_compton_analysis_TOF();
		~JEventProcessor_compton_analysis_TOF(){};
		const char* className(void){return "JEventProcessor_compton_analysis_TOF";}
	
	private:
		jerror_t init(void);
		jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);
		jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);
		jerror_t erun(void);
		jerror_t fini(void);
		
		int fcal_fiducial_cut(DVector3 pos, DVector3 vertex, double layer_cut);
		int ccal_fiducial_cut(DVector3 pos, DVector3 vertex);
		
		void set_cuts(int32_t runnumber);
		
		double get_vertex_weight(double vertex_z);
		
		void check_TOF_match(DVector3 pos, double rfTime, DVector3 vertex, 
			vector<const DTOFPoint*> tof_points, double &dx_min, double &dy_min, 
			double &dt_min, double rf_time_cut);
		
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
		
		DVector3 m_fcal_correction, m_ccal_correction;
		
		const double c    =  29.9792458;
		const double m_e  =  0.510998928e-3;
		
		const double FCAL_RF_time_cut = 3.0;
		const double CCAL_RF_time_cut = 3.0;
		const double BEAM_RF_time_cut = 2.004;
		
		const double FCAL_min_energy_cut = 0.5;
		const double CCAL_min_energy_cut = 3.0;
		const double BEAM_min_energy_cut = 6.0;
		
		int bfield_val = 0;
		int phase_val  = 0;
		
		//----------   Histograms   ----------//
		
		TH1F *h_beam;
		TH2F *h_vertex_xy;
		TH1F *h_vertex;
		TH1F *h_vertex_weight;
		TH1F *h_vertex_accepted;
		TH1F *h_reaction_weight;
		
		TH1F *h_mass_thrown;
		TH2F *h_rec_vs_thrown_mass, *h_rec_vs_thrown_mass_cut;
		TH2F *h_deltaE_vs_thrown_mass;
		TH2F *h_deltaK_vs_thrown_mass;
		
		static const int n_cuts = 4;
		TH2F *h_deltaE_tagh[n_cuts],     *h_deltaE_tagm[n_cuts];
		TH2F *h_deltaK_tagh[n_cuts],     *h_deltaK_tagm[n_cuts];
		TH2F *h_deltaK_tagh_cut[n_cuts], *h_deltaK_tagm_cut[n_cuts];
		TH2F *h_elas_vs_deltaE[n_cuts];
		TH2F *h_deltaK_vs_deltaE[n_cuts];
		TH2F *h_deltaCCAL_vs_deltaE[n_cuts], *h_deltaCCAL_vs_deltaK[n_cuts], *h_deltaCCAL_vs_deltaK_cut[n_cuts];
		TH2F *h_deltaFCAL_vs_deltaE[n_cuts], *h_deltaFCAL_vs_deltaK[n_cuts], *h_deltaFCAL_vs_deltaK_cut[n_cuts];
		TH2F *h_mgg_vs_deltaE[n_cuts], *h_mgg_vs_deltaK[n_cuts];
		TH2F *h_ccal_xy[n_cuts], *h_fcal_xy[n_cuts];
		
		TH1F *h_n_showers,          *h_n_showers_ccal,          *h_n_showers_fcal;
		TH1F *h_n_showers_cut,      *h_n_showers_ccal_cut,      *h_n_showers_fcal_cut;
		TH1F *h_n_good_showers,     *h_n_good_showers_ccal,     *h_n_good_showers_fcal;
		TH1F *h_n_good_showers_cut, *h_n_good_showers_ccal_cut, *h_n_good_showers_fcal_cut;
		
		TH1F *h_tof_match_case, *h_tof_match_case_cut;
		
		TH1F *h_extra_fcal_shower_energy;
		TH1F *h_extra_fcal_shower_distance;
		TH1F *h_extra_fcal_shower_deltaPhi;
		TH1F *h_extra_fcal_shower_deltaPhi_ccal;
		TH1F *h_extra_fcal_shower_elasticity;
		TH1F *h_extra_fcal_shower_elasticity_new;
		TH1F *h_extra_fcal_shower_elasticity_corr;
		TH2F *h_extra_fcal_shower_xy;
		
		TH1F *h_extra_fcal_shower_energy_cut;
		TH1F *h_extra_fcal_shower_distance_cut;
		TH1F *h_extra_fcal_shower_deltaPhi_cut;
		TH1F *h_extra_fcal_shower_deltaPhi_ccal_cut;
		TH1F *h_extra_fcal_shower_elasticity_cut;
		TH1F *h_extra_fcal_shower_elasticity_new_cut;
		TH1F *h_extra_fcal_shower_elasticity_corr_cut;
		TH2F *h_extra_fcal_shower_xy_cut;
		
		TH1F *h_extra_ccal_shower_energy;
		TH1F *h_extra_ccal_shower_distance;
		TH1F *h_extra_ccal_shower_deltaPhi;
		TH1F *h_extra_ccal_shower_deltaPhi_fcal;
		TH1F *h_extra_ccal_shower_elasticity;
		TH1F *h_extra_ccal_shower_elasticity_new;
		TH1F *h_extra_ccal_shower_elasticity_corr;
		TH2F *h_extra_ccal_shower_xy;
		
		TH1F *h_extra_ccal_shower_energy_cut;
		TH1F *h_extra_ccal_shower_distance_cut;
		TH1F *h_extra_ccal_shower_deltaPhi_cut;
		TH1F *h_extra_ccal_shower_deltaPhi_fcal_cut;
		TH1F *h_extra_ccal_shower_elasticity_cut;
		TH1F *h_extra_ccal_shower_elasticity_new_cut;
		TH1F *h_extra_ccal_shower_elasticity_corr_cut;
		TH2F *h_extra_ccal_shower_xy_cut;
		
		TH1F *h_extra_ccal_shower_energy_e;
		TH1F *h_extra_ccal_shower_distance_e;
		TH1F *h_extra_ccal_shower_deltaPhi_e;
		TH1F *h_extra_ccal_shower_deltaPhi_fcal_e;
		TH1F *h_extra_ccal_shower_elasticity_e;
		TH1F *h_extra_ccal_shower_elasticity_new_e;
		TH1F *h_extra_ccal_shower_elasticity_corr_e;
		TH2F *h_extra_ccal_shower_xy_e;
		
		TH1F *h_extra_ccal_shower_energy_e_cut;
		TH1F *h_extra_ccal_shower_distance_e_cut;
		TH1F *h_extra_ccal_shower_deltaPhi_e_cut;
		TH1F *h_extra_ccal_shower_deltaPhi_fcal_e_cut;
		TH1F *h_extra_ccal_shower_elasticity_e_cut;
		TH1F *h_extra_ccal_shower_elasticity_new_e_cut;
		TH1F *h_extra_ccal_shower_elasticity_corr_e_cut;
		TH2F *h_extra_ccal_shower_xy_e_cut;
		
		//----------      Cuts      ---------//
		
		int    m_USE_REACTION_WEIGHT;
		double m_REACTION_CUT_WEIGHT;
		
		TF1 *f_deltaE_mu_mc,   *f_deltaE_sig_mc;
		TF1 *f_deltaE_mu_data, *f_deltaE_sig_data;
		
		double deltaE_mu_p0_mc,    deltaE_mu_p1_mc,    deltaE_mu_p2_mc,    deltaE_mu_p3_mc;
		double deltaE_mu_p0_data,  deltaE_mu_p1_data,  deltaE_mu_p2_data,  deltaE_mu_p3_data;
		
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
		
		double deltaK_mu_p0_mc,    deltaK_mu_p1_mc,    deltaK_mu_p2_mc,    deltaK_mu_p3_mc;
		double deltaK_mu_p0_data,  deltaK_mu_p1_data,  deltaK_mu_p2_data,  deltaK_mu_p3_data;
		
		double deltaK_sig_p0_mc,   deltaK_sig_p1_mc,   deltaK_sig_p2_mc,   deltaK_sig_p3_mc;
		double deltaK_sig_p0_data, deltaK_sig_p1_data, deltaK_sig_p2_data, deltaK_sig_p3_data;
		
		//-----------------------------------//
};

#endif // _JEventProcessor_compton_analysis_TOF_

