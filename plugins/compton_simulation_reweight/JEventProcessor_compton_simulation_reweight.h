// $Id$
//
//    File: JEventProcessor_compton_simulation_reweight.h
// Created: Mon Nov 23 16:24:15 EST 2020
// Creator: andrsmit (on Linux ifarm1901.jlab.org 3.10.0-1062.4.1.el7.x86_64 x86_64)
//

#ifndef _JEventProcessor_compton_simulation_reweight_
#define _JEventProcessor_compton_simulation_reweight_

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

class JEventProcessor_compton_simulation_reweight:public jana::JEventProcessor{
	
	public:
		JEventProcessor_compton_simulation_reweight() {};
		~JEventProcessor_compton_simulation_reweight() {};
		const char* className(void){return "JEventProcessor_compton_simulation_reweight";}

	private:
		jerror_t init(void);
		jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);
		jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);
		jerror_t erun(void);
		jerror_t fini(void);
		
		int ccalLayer(int row, int col);
		int fcalLayer(int row, int col);
		
		void fill_histograms( vector< ComptonCandidate_t > Comp_Candidates );
		void set_cuts( int32_t runnumber );
		
		
		//----------   Constants   ----------//
		
		double m_beamX, m_beamY, m_beamZ;
		double m_fcalX, m_fcalY, m_fcalZ;
		double m_ccalX, m_ccalY, m_ccalZ;
		
		double m_ccalX_new, m_ccalY_new;
		double m_beamX_new, m_beamY_new;
		
		const double c    =  29.9792458;
		const double m_e  =  0.510998928e-3;
		
		const double FCAL_RF_time_cut = 3.0;
		const double CCAL_RF_time_cut = 3.0;
		const double BEAM_RF_time_cut = 2.004;
		
		const double FCAL_min_energy_cut = 0.5;
		const double CCAL_min_energy_cut = 1.0;
		const double BEAM_min_energy_cut = 6.0;
		
		
		//----------   Histograms   ----------//
		
		TH1F *h_fcal_rf_dt,  *h_ccal_rf_dt,    *h_beam_rf_dt;
		
		TH1F *h_beam;
		TH1F *h_tagh_flux,   *h_tagm_flux;
		TH1F *h_double_compton;
		
		TH1F *h_vertex,      *h_vertex_weighted;
		TH1F *h_vertex_cut,  *h_vertex_cut_weighted;
		
		TH2F *h_deltaE_tagh,     *h_deltaE_tagh_weighted;
		TH2F *h_deltaE_tagm,     *h_deltaE_tagm_weighted;
		
		TH2F *h_deltaPhi_tagh,   *h_deltaPhi_tagh_weighted;
		TH2F *h_deltaPhi_tagm,   *h_deltaPhi_tagm_weighted;
		
		TH2F *h_deltaK_tagh,     *h_deltaK_tagh_weighted;
		TH2F *h_deltaK_tagm,     *h_deltaK_tagm_weighted;
		
		TH2F *h_deltaK_tagh_cut, *h_deltaK_tagh_cut_weighted;
		TH2F *h_deltaK_tagm_cut, *h_deltaK_tagm_cut_weighted;
		
		TH2F *h_fcal_xy, *h_fcal_xy_weighted;
		TH2F *h_ccal_xy, *h_ccal_xy_weighted;
		
		
		//----------      Cuts      ---------//
		
		// DeltaE mu function is 3rd order polynomial
		// DeltaE sig/E function is [0] + [1]/sqrt(x) + [2]/x
		
		TF1 *f_deltaE_mu,   *f_deltaE_sig;
		
		double deltaE_mu_p0,  deltaE_mu_p1,  deltaE_mu_p2,  deltaE_mu_p3;
		double deltaE_sig_p0, deltaE_sig_p1, deltaE_sig_p2;
		
		// DeltaPhi mu function is 3rd order polynomial
		// DeltaPhi sig function is 3rd order polynomial
		
		TF1 *f_deltaPhi_mu, *f_deltaPhi_sig;
		
		double deltaPhi_mu_p0,  deltaPhi_mu_p1,  deltaPhi_mu_p2,  deltaPhi_mu_p3;
		double deltaPhi_sig_p0, deltaPhi_sig_p1, deltaPhi_sig_p2, deltaPhi_sig_p3;
		
		// DeltaK mu function is 3rd order polynomial
		// DeltaK sig function is 3rd order polynomial
		
		TF1 *f_deltaK_mu,   *f_deltaK_sig;
		
		double deltaK_mu_p0,  deltaK_mu_p1,  deltaK_mu_p2,  deltaK_mu_p3;
		double deltaK_sig_p0, deltaK_sig_p1, deltaK_sig_p2, deltaK_sig_p3;
		
		
		//-----------------------------------//
		
};

#endif // _JEventProcessor_compton_simulation_reweight_

