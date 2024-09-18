// $Id$
//
//    File: JEventProcessor_compton_analysis.h
// Created: Thu Feb 10 21:33:24 EST 2022
// Creator: andrsmit (on Linux ifarm1802.jlab.org 3.10.0-1160.11.1.el7.x86_64 x86_64)
//

#ifndef _JEventProcessor_compton_analysis_
#define _JEventProcessor_compton_analysis_

#include <JANA/JEventProcessor.h>
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>

#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TRandom3.h"
#include "TSystem.h"

using namespace jana;
using namespace std;

#include "CCAL/DCCALShower.h"
#include "FCAL/DFCALShower.h"
#include "PID/DBeamPhoton.h"
#include "PID/DEventRFBunch.h"

#include "TRACKING/DMCThrown.h"
#include "PID/DMCReaction.h"
#include "TRIGGER/DL1Trigger.h"
#include "TRIGGER/DTrigger.h"
#include "DAQ/Df250WindowRawData.h"
#include "DAQ/Df250PulseData.h"

#include "HDGEOMETRY/DGeometry.h"

#include "particleType.h"
#include "units.h"
#include "DLorentzVector.h"
#include "DVector3.h"

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <thread>
#include <mutex>

#include <CompCand.h>

class JEventProcessor_compton_analysis:public jana::JEventProcessor{
	public:
		JEventProcessor_compton_analysis();
		~JEventProcessor_compton_analysis(){};
		const char* className(void){return "JEventProcessor_compton_analysis";}

	private:
		jerror_t init(void);
		jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);
		jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);
		jerror_t erun(void){return NOERROR;};
		jerror_t fini(void){return NOERROR;};
		
		//---------------------------------------//
		// Unique function declarations:
		
		int fcal_fiducial_cut(DVector3 pos);
		int ccal_fiducial_cut(DVector3 pos);
		
		double get_acc_scaling_factor(double eb);
		
		void fill_histograms(vector<ComptonCandidate_t> Comp_Candidates);
		void set_cuts(int32_t runnumber);
		
		int cut_deltaE(  double deltaE,   double eb, double n_sigma_left, double n_sigma_right);
		int cut_deltaPhi(double deltaPhi, double eb, double n_sigma_left, double n_sigma_right);
		int cut_deltaK(  double deltaK,   double eb, double n_sigma_left, double n_sigma_right);
		
		//---------------------------------------//
		// Physical Constants
		
		const double m_c  = 29.9792458;
		const double m_e  = 0.510998928e-3;
		
		//---------------------------------------//
		// Geometry
		
		Particle_t m_Target;
		
		DVector3 m_beamSpot;
		DVector3 m_fcalFace, m_ccalFace;
		DVector3 m_fcal_correction, m_ccal_correction;
		
		double m_target_length;
		double m_target_density;
		double m_atten;
		
		int m_bfield_val = 0;
		int m_phase_val  = 0;
		
		//---------------------------------------//
		// Cuts (defaults set inside constructor)
		
		double m_cut_fcalrfdt, m_cut_ccalrfdt, m_cut_beamrfdt;
		double m_cut_fcalE,    m_cut_ccalE,    m_cut_beamE;
		double m_cut_deltaE,   m_cut_deltaPhi, m_cut_deltaK;
		double m_cut_fcal_layers;
		int    m_beam_bunches_acc;
		
		double m_deltaE_mu_pars[4],   m_deltaE_sigma_pars[3];
		double m_deltaPhi_mu_pars[4], m_deltaPhi_sigma_pars[4];
		double m_deltaK_mu_pars[4],   m_deltaK_sigma_pars[4];
		
		double m_HodoscopeHiFactor    = 1.0;
		double m_HodoscopeHiFactorErr = 1.0;
		double m_HodoscopeLoFactor    = 1.0;
		double m_HodoscopeLoFactorErr = 1.0;
		double m_MicroscopeFactor     = 1.0;
		double m_MicroscopeFactorErr  = 1.0;
		double m_TAGMEnergyBoundHi    = 1.0;
		double m_TAGMEnergyBoundLo    = 1.0;
		
		//---------------------------------------//
		// Histograms
		
		// Trigger:
		TH1F *h_trig, *h_fptrig;
		
		// Timing histograms:
		TH1F *h_fcal_rf_dt, *h_ccal_rf_dt, *h_beam_rf_dt;
		TH1F *h_beam_rf_dt_cut;
		TH2F *h_beam_rf_dt_tagh, *h_beam_rf_dt_tagm;
		
		// Multiplicity histograms:
		TH1F *h_nccal, *h_nfcal, *h_n_ccal_fcal, *h_n_ccal_fcal_cut;
		TH2F *h_ccalE,     *h_fcalE;
		TH2F *h_ccalE_cut, *h_fcalE_cut;
		
		// Compton-scattering distributions:
		TH2F *h_deltaE_tagh,      *h_deltaE_tagm;
		TH2F *h_deltaPhi_tagh,    *h_deltaPhi_tagm;
		TH2F *h_deltaK_tagh,      *h_deltaK_tagm;
		TH2F *h_deltaK_tagh_cut,  *h_deltaK_tagm_cut;
		TH2F *h_deltaK_tagh_main, *h_deltaK_tagm_main;
		TH2F *h_deltaK_tagh_acc,  *h_deltaK_tagm_acc;
		TH2F *h_deltaK_vs_deltaE;
		TH2F *h_fcal_xy, *h_ccal_xy;
};

#endif // _JEventProcessor_compton_analysis_

