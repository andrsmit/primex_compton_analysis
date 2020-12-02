// $Id$
//
//    File: JEventProcessor_compton_analysis2.h
// Created: Mon Nov 23 16:24:15 EST 2020
// Creator: andrsmit (on Linux ifarm1901.jlab.org 3.10.0-1062.4.1.el7.x86_64 x86_64)
//

#ifndef _JEventProcessor_compton_analysis2_
#define _JEventProcessor_compton_analysis2_

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

class JEventProcessor_compton_analysis2:public jana::JEventProcessor{
	
	public:
		JEventProcessor_compton_analysis2();
		~JEventProcessor_compton_analysis2() {};
		const char* className(void){return "JEventProcessor_compton_analysis2";}

	private:
		jerror_t init(void);
		jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);
		jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);
		jerror_t erun(void);
		jerror_t fini(void);
		
		int ccalLayer(int row, int col);
		int fcalLayer(int row, int col);
		
		void fill_histograms( vector< ComptonCandidate_t > Comp_Candidates );
		void read_cuts();
		
		int  RUN_GROUP;
		char cut_pathName[256];
		
		
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
		
		TH1F *hTrig,         *hfpTrig;
		TH1F *h_fcal_rf_dt,  *h_ccal_rf_dt,    *h_beam_rf_dt;
		
		TH2F *h_deltaPhi[2], *h_deltaPhi_t[2], *h_deltaPhi_te[2], *h_deltaPhi_tek[2], *h_deltaPhi_tekp[2];
		TH2F *h_deltaT[2],   *h_deltaT_e[2],   *h_deltaT_ep[2],   *h_deltaT_epk[2],   *h_deltaT_epkt[2];
		TH2F *h_deltaR[2],   *h_deltaR_t[2],   *h_deltaR_te[2],   *h_deltaR_tep[2],   *h_deltaR_tepk[2];
		TH2F *h_deltaE[2],   *h_deltaE_t[2],   *h_deltaE_tp[2],   *h_deltaE_tpk[2],   *h_deltaE_tpke[2];
		TH2F *h_deltaK[2],   *h_deltaK_t[2],   *h_deltaK_te[2],   *h_deltaK_tep[2],   *h_deltaK_tepk[2];
		TH2F *h_deltaK2[2],  *h_deltaK2_t[2],  *h_deltaK2_te[2],  *h_deltaK2_tep[2],  *h_deltaK2_tepk[2];
		
		TH2F *h_fcal_xy,     *h_ccal_xy;
		
		
		//----------      Cuts      ---------//
		
		double   deltaE_mu_tagh[274],   deltaE_sig_tagh[274];
		double   deltaE_mu_tagm[102],   deltaE_sig_tagm[102];
		
		double deltaPhi_mu_tagh[274], deltaPhi_sig_tagh[274];
		double deltaPhi_mu_tagm[102], deltaPhi_sig_tagm[274];
		
		double   deltaK_mu_tagh[274],   deltaK_sig_tagh[274];
		double   deltaK_mu_tagm[102],   deltaK_sig_tagm[102];
		
		
		// DeltaT mu function is 2nd order polynomial
		// DeltaT sig function is 2nd order polynomial
		
		TF1 *f_deltaT_mu,   *f_deltaT_sig;
		
		double deltaT_mu_p0  =  1.41009e-01;
		double deltaT_mu_p1  = -2.72818e-02;
		double deltaT_mu_p2  =  1.16503e-03;
		
		double deltaT_sig_p0 =  5.43319e-01;
		double deltaT_sig_p1 = -4.34696e-02;
		double deltaT_sig_p2 =  2.32939e-03;
		
		//-----------------------------------//
		
};

#endif // _JEventProcessor_compton_analysis2_

