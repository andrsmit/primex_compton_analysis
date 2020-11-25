// $Id$
//
//    File: JEventProcessor_compton_analysis.h
// Created: Mon Nov 23 16:24:15 EST 2020
// Creator: andrsmit (on Linux ifarm1901.jlab.org 3.10.0-1062.4.1.el7.x86_64 x86_64)
//

#ifndef _JEventProcessor_compton_analysis_
#define _JEventProcessor_compton_analysis_

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

class JEventProcessor_compton_analysis:public jana::JEventProcessor{
	
	public:
		JEventProcessor_compton_analysis();
		~JEventProcessor_compton_analysis() {};
		const char* className(void){return "JEventProcessor_compton_analysis";}

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
		
		TH2F *h_deltaPhi[2], *h_deltaPhi_e[2], *h_deltaPhi_ek[2], *h_deltaPhi_ekp[2];
		TH2F *h_deltaT[2],   *h_deltaT_e[2],   *h_deltaT_ep[2],   *h_deltaT_epk[2];
		TH2F *h_deltaR[2],   *h_deltaR_e[2],   *h_deltaR_ep[2],   *h_deltaR_epk[2];
		TH2F *h_deltaE[2],   *h_deltaE_p[2],   *h_deltaE_pk[2],   *h_deltaE_pke[2];
		TH2F *h_deltaK[2],   *h_deltaK_e[2],   *h_deltaK_ep[2],   *h_deltaK_epk[2];
		TH2F *h_deltaK2[2],  *h_deltaK2_e[2],  *h_deltaK2_ep[2],  *h_deltaK2_epk[2];
		
		TH2F *h_fcal_xy,     *h_ccal_xy;
		
		
		//----------      Cuts      ---------//
		
		double   deltaE_mu_tagh[274],   deltaE_sig_tagh[274];
		double   deltaE_mu_tagm[102],   deltaE_sig_tagm[102];
		
		double deltaPhi_mu_tagh[274], deltaPhi_sig_tagh[274];
		double deltaPhi_mu_tagm[102], deltaPhi_sig_tagm[274];
		
		double   deltaK_mu_tagh[274],   deltaK_sig_tagh[274];
		double   deltaK_mu_tagm[102],   deltaK_sig_tagm[102];
		
		//-----------------------------------//
		
};

#endif // _JEventProcessor_compton_analysis_

