// $Id$
//
//    File: JEventProcessor_CCAL_ComptonGains.h
// Created: Mon Dec 03 16:04:16 EST 2019
// Creator: andrsmit (on Linux ifarm1402)
//

#ifndef _JEventProcessor_CCAL_ComptonGains_
#define _JEventProcessor_CCAL_ComptonGains_

#include <JANA/JEventProcessor.h>
#include <JANA/JApplication.h>

#include "TTree.h"
#include "TH1.h"
#include "TH2.h"

#include <PID/DBeamPhoton.h>
#include <PID/DEventRFBunch.h>

#include <CCAL/DCCALShower.h>
#include <FCAL/DFCALShower.h>
#include <BCAL/DBCALShower.h>
#include <TOF/DTOFPoint.h>

#include <TRIGGER/DL1Trigger.h>
#include <TRIGGER/DTrigger.h>
#include <DAQ/Df250WindowRawData.h>
#include <DAQ/Df250PulseData.h>

#include <HDGEOMETRY/DGeometry.h>

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

using namespace jana;
using namespace std;


class JEventProcessor_CCAL_ComptonGains:public JEventProcessor {
 	
	public:
		
		JEventProcessor_CCAL_ComptonGains();
  		~JEventProcessor_CCAL_ComptonGains(){};
  		const char* className(void){return "JEventProcessor_CCAL_ComptonGains";}
		
 	private:
		
		jerror_t init(void);
		jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);
		jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);
		jerror_t erun(void){return NOERROR;};
		jerror_t fini(void){return NOERROR;};
		
		int fcal_fiducial_cut(DVector3 pos, DVector3 vertex, double cut_layer);
		int ccal_fiducial_cut(DVector3 pos, DVector3 vertex);
		
		void check_TOF_match(DVector3 pos, double rfTime, DVector3 vertex, 
			vector<const DTOFPoint*> tof_points, double &dx_min, double &dy_min, 
			double &dt_min, double rf_time_cut);
		
		//----- Cut Parameters -----//
		
		const double RF_TIME_CUT = 6.012; 
		const int n_bunches_acc  = 2;
		
		bool BYPASS_TRIG;
		
		double MIN_FCAL_ENERGY_CUT;
		double MIN_CCAL_ENERGY_CUT;
		double MIN_BEAM_ENERGY_CUT;
		
		double COPL_CUT;
		double DELTA_K_CUT;
		
		double MIN_CALIB_ENERGY, MAX_CALIB_ENERGY;
		
		//----- Physical Constants -----//
		
		const double c    =  29.9792458;
		const double m_e  =  0.510998928e-3;
		
		double m_beamX, m_beamY, m_beamZ;
		double m_fcalX, m_fcalY, m_fcalZ;
		double m_ccalX, m_ccalY, m_ccalZ;
		
		double m_fcalX_new, m_fcalY_new;
		double m_ccalX_new, m_ccalY_new;
		
		DVector3 m_fcal_correction, m_ccal_correction;
		
		//----- Histograms -----//
		
		TH1F *h_trig, *h_fptrig;
		
		TH1F *h_tof_rf_dt,     *h_fcal_tof_dt;
		TH1F *h_tof_rf_dt_cut, *h_fcal_tof_dt_cut;
		
		TH1F *h_fcal_rf_dt, *h_ccal_rf_dt;
		TH1F *h_beam_rf_dt, *h_ccal_beam_dt, *h_fcal_beam_dt;
		
		TH1F *h_deltaPhi, *h_deltaT;
		TH1F *h_deltaE,   *h_deltaK;
		
		static const int n_cuts = 3;
		
		TH2F *h_CompRatio[n_cuts], *h_CompRatio_g[n_cuts], *h_CompRatio_e[n_cuts];
		TH2F *h_ElasRatio[n_cuts], *h_ElasRatio_g[n_cuts], *h_ElasRatio_e[n_cuts];
		
		TH2F *h_CompRatio_cut[n_cuts], *h_CompRatio_g_cut[n_cuts], *h_CompRatio_e_cut[n_cuts];
		TH2F *h_ElasRatio_cut[n_cuts], *h_ElasRatio_g_cut[n_cuts], *h_ElasRatio_e_cut[n_cuts];
		
		TH2F *h_CompRatio_vs_E[n_cuts], *h_CompRatio_vs_E_g[n_cuts], *h_CompRatio_vs_E_e[n_cuts];
		TH2F *h_ElasRatio_vs_E[n_cuts], *h_ElasRatio_vs_E_g[n_cuts], *h_ElasRatio_vs_E_e[n_cuts];
		
		TH2F *h_xy_fcal, *h_xy_ccal;
};

#endif // _JEventProcessor_CCAL_ComptonGains_

