// $Id$
//
//    File: JEventProcessor_compton_simulation_systematics.h
// Created: Mon Nov 23 16:24:15 EST 2020
// Creator: andrsmit (on Linux ifarm1901.jlab.org 3.10.0-1062.4.1.el7.x86_64 x86_64)
//

#ifndef _JEventProcessor_compton_simulation_systematics_
#define _JEventProcessor_compton_simulation_systematics_

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

#include <CCAL/DCCALGeometry.h>
#include <FCAL/DFCALGeometry.h>

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

class JEventProcessor_compton_simulation_systematics:public jana::JEventProcessor{
	
	public:
		JEventProcessor_compton_simulation_systematics();
		~JEventProcessor_compton_simulation_systematics() {};
		const char* className(void){return "JEventProcessor_compton_simulation_systematics";}

	private:
		jerror_t init(void);
		jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);
		jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);
		jerror_t erun(void);
		jerror_t fini(void);
		
		int fcal_fiducial_cut(DVector3 pos, DVector3 vertex);
		int ccal_fiducial_cut(DVector3 pos, DVector3 vertex);
		
		void fill_histograms(vector< ComptonCandidate_t > Comp_Candidates, DVector3 vertex);
		void set_cuts(int32_t runnumber);
		
		double get_vertex_weight(double vertex_z);
		
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
		
		const double c   = 29.9792458;
		const double m_e = 0.510998928e-3;
		
		const double FCAL_RF_time_cut = 3.0;
		const double CCAL_RF_time_cut = 3.0;
		const double BEAM_RF_time_cut = 2.004;
		
		const double FCAL_min_energy_cut = 0.35;
		const double CCAL_min_energy_cut = 3.0;
		const double BEAM_min_energy_cut = 3.0;
		
		int bfield_val = 0;
		int phase_val  = 0;
		
		//----------   Histograms   ----------//
		
		TH1F *h_beam;
		TH1F *h_tagh_flux, *h_tagm_flux;
		TH1F *h_double_compton;
		
		TH2F *h_vertex_xy;
		TH1F *h_vertex;
		TH1F *h_vertex_weight;
		TH1F *h_vertex_accepted;
		TH1F *h_reaction_weight;
		
		TH2F *h_deltaK_tagh_fcalE[20],      *h_deltaK_tagm_fcalE[20];
		TH2F *h_deltaK_tagh_cut_fcalE[20],  *h_deltaK_tagm_cut_fcalE[20];
		
		TH2F *h_deltaK_tagh_ccalE[13],      *h_deltaK_tagm_ccalE[13];
		TH2F *h_deltaK_tagh_cut_ccalE[13],  *h_deltaK_tagm_cut_ccalE[13];
		
		double cut_sigmas[16] = {
		1.0, 
		2.0, 
		3.0, 
		3.5, 
		4.0, 
		4.5, 
		5.0, 
		5.5, 
		6.0, 
		6.5, 
		7.0, 
		8.0, 
		9.0, 
		10.0, 
		15.0, 
		20.0
		};
		
		TH2F *h_deltaK_tagh_sigE[16],       *h_deltaK_tagm_sigE[16];
		TH2F *h_deltaK_tagh_cut_sigE[16],   *h_deltaK_tagm_cut_sigE[16];
		
		TH2F *h_deltaK_tagh_sigPhi[16],     *h_deltaK_tagm_sigPhi[16];
		TH2F *h_deltaK_tagh_cut_sigPhi[16], *h_deltaK_tagm_cut_sigPhi[16];
		
		TH2F *h_deltaK_tagh_sigK[16],       *h_deltaK_tagm_sigK[16];
		
		TH2F *h_deltaK_tagh_fcalT[10],      *h_deltaK_tagm_fcalT[10];
		TH2F *h_deltaK_tagh_cut_fcalT[10],  *h_deltaK_tagm_cut_fcalT[10];
		
		TH2F *h_deltaK_tagh_ccalT[10],      *h_deltaK_tagm_ccalT[10];
		TH2F *h_deltaK_tagh_cut_ccalT[10],  *h_deltaK_tagm_cut_ccalT[10];
		
		static const int N_FID_CUTS = 20;
		TH2F *h_deltaK_tagh_fcalfid[N_FID_CUTS];
		TH2F *h_deltaK_tagm_fcalfid[N_FID_CUTS];
		TH2F *h_deltaK_tagh_ccalfid[N_FID_CUTS];
		TH2F *h_deltaK_tagm_ccalfid[N_FID_CUTS];
		
		TH2F *h_deltaK_tagh_fcal_phi[8],   *h_deltaK_tagm_fcal_phi[8];
		TH2F *h_deltaK_tagh_ccal_phi[8],   *h_deltaK_tagm_ccal_phi[8];
		TH2F *h_deltaK_tagh_fcal_layer[8], *h_deltaK_tagm_fcal_layer[8];
		TH2F *h_deltaK_tagh_ccal_layer[5], *h_deltaK_tagm_ccal_layer[5];
		
		TH2F *h_deltaK_tagh_cut_fcal_phi[8],   *h_deltaK_tagm_cut_fcal_phi[8];
		TH2F *h_deltaK_tagh_cut_ccal_phi[8],   *h_deltaK_tagm_cut_ccal_phi[8];
		TH2F *h_deltaK_tagh_cut_fcal_layer[8], *h_deltaK_tagm_cut_fcal_layer[8];
		TH2F *h_deltaK_tagh_cut_ccal_layer[5], *h_deltaK_tagm_cut_ccal_layer[5];
		
		TH2F *h_xy_fcal_phi[8],   *h_xy_ccal_phi[8];
		TH2F *h_xy_fcal_layer[8], *h_xy_ccal_layer[5];
		
		TH2F *h_deltaK_tagh_single,    *h_deltaK_tagm_single;
		
		//----------      Cuts      ---------//
		
		int    m_USE_REACTION_WEIGHT;
		double m_REACTION_CUT_WEIGHT;
		
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
		
		//-----------------------------------//
};

#endif // _JEventProcessor_compton_simulation_systematics_

