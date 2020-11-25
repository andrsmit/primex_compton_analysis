// $Id$
//
//    File: JEventProcessor_compton_simulation.h
// Created: Mon Nov 23 16:24:15 EST 2020
// Creator: andrsmit (on Linux ifarm1901.jlab.org 3.10.0-1062.4.1.el7.x86_64 x86_64)
//

#ifndef _JEventProcessor_compton_simulation_
#define _JEventProcessor_compton_simulation_

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

class JEventProcessor_compton_simulation:public jana::JEventProcessor{
	
	public:
		JEventProcessor_compton_simulation();
		~JEventProcessor_compton_simulation() {};
		const char* className(void){return "JEventProcessor_compton_simulation";}

	private:
		jerror_t init(void);
		jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);
		jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);
		jerror_t erun(void);
		jerror_t fini(void);
		
		int ccalLayer(int row, int col);
		int fcalLayer(int row, int col);
		
		void fill_histograms( vector< ComptonCandidate_t > Comp_Candidates );
		
		
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
		
		TH2F *h_deltaPhi[2], *h_deltaPhi_e[2], *h_deltaPhi_ek[2], *h_deltaPhi_ekp[2];
		TH2F *h_deltaT[2],   *h_deltaT_e[2],   *h_deltaT_ep[2],   *h_deltaT_epk[2];
		TH2F *h_deltaR[2],   *h_deltaR_e[2],   *h_deltaR_ep[2],   *h_deltaR_epk[2];
		TH2F *h_deltaE[2],   *h_deltaE_p[2],   *h_deltaE_pk[2],   *h_deltaE_pke[2];
		TH2F *h_deltaK[2],   *h_deltaK_e[2],   *h_deltaK_ep[2],   *h_deltaK_epk[2];
		TH2F *h_deltaK2[2],  *h_deltaK2_e[2],  *h_deltaK2_ep[2],  *h_deltaK2_epk[2];
		
		TH2F *h_fcal_xy,     *h_ccal_xy;
		
		
		//----------      Cuts      ---------//
		
		// DeltaE mu function is 3rd order polynomial
		// DeltaE sig/E function is [0] + [1]/sqrt(x) + [2]/x
		
		TF1 *f_deltaE_mu,   *f_deltaE_sig;
		
		double deltaE_mu_p0  = -5.19163e-02;
		double deltaE_mu_p1  = -3.00848e-03;
		double deltaE_mu_p2  = -3.80583e-04;
		double deltaE_mu_p3  =  6.22919e-05;
		
		double deltaE_sig_p0 =  2.11133e-03;
		double deltaE_sig_p1 =  5.01150e-02;
		double deltaE_sig_p2 = -2.65049e-02;
		
		
		
		// DeltaPhi mu function is 3rd order polynomial
		// DeltaPhi sig function is 3rd order polynomial
		
		TF1 *f_deltaPhi_mu, *f_deltaPhi_sig;
		
		double deltaPhi_mu_p0  =  1.81175e+02;
		double deltaPhi_mu_p1  = -4.54542e-01;
		double deltaPhi_mu_p2  =  4.61759e-02;
		double deltaPhi_mu_p3  = -1.66371e-03;
		
		double deltaPhi_sig_p0 =  8.47957e+00;
		double deltaPhi_sig_p1 = -6.26961e-01;
		double deltaPhi_sig_p2 =  5.07839e-02;
		double deltaPhi_sig_p3 = -1.29749e-03;
		
		
		
		// DeltaK mu function is 3rd order polynomial
		// DeltaK sig function is 3rd order polynomial
		
		TF1 *f_deltaK_mu,   *f_deltaK_sig;
		
		double deltaK_mu_p0  =  2.66494e-01;
		double deltaK_mu_p1  = -8.37391e-02;
		double deltaK_mu_p2  =  7.70817e-03;
		double deltaK_mu_p3  = -2.79719e-04;
		
		double deltaK_sig_p0 =  1.81033e-01;
		double deltaK_sig_p1 = -3.39539e-02;
		double deltaK_sig_p2 =  6.40301e-03;
		double deltaK_sig_p3 = -2.56990e-04;
		
		
		//-----------------------------------//
		
};

#endif // _JEventProcessor_compton_simulation_

