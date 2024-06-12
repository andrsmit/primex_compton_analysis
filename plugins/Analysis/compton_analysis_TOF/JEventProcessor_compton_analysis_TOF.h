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

#include <CompCand.h>

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
  
  void fill_histograms(vector<ComptonCandidate_t> Comp_Candidates, DVector3 vertex, 
		double rfTime, vector<const DTOFPoint*> tof_points);
	
  void set_cuts(int32_t runnumber);
	
	double get_vertex_weight(double vertex_z);
  
  int check_TOF_match(DVector3 pos, double rfTime, DVector3 vertex, 
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
  
  DVector3 fcal_correction, ccal_correction;
  
  const double c    =  29.9792458;
  const double m_e  =  0.510998928e-3;
  
  const double FCAL_RF_time_cut = 3.0;
  const double CCAL_RF_time_cut = 3.0;
  const double BEAM_RF_time_cut = 2.004;
  
  const double FCAL_min_energy_cut = 0.35;
  const double CCAL_min_energy_cut = 3.0;
  const double BEAM_min_energy_cut = 6.0;
  
  int bfield_val = 0;
  int phase_val  = 0;
  
  //----------   Histograms   ----------//
  
	TH1F *h_beam;
	TH1F *h_tagh_flux, *h_tagm_flux;
	
	TH2F *h_vertex_xy;
	TH1F *h_vertex;
	TH1F *h_vertex_weight;
	TH1F *h_vertex_accepted;
	TH1F *h_reaction_weight;
	
  TH2F *h_tof_match_egam,  *h_tof_match_egam_cut;
	TH2F *h_tof_match_theta, *h_tof_match_theta_cut;
	TH1F *h_tof_rf_dt,       *h_tof_rf_dt_cut;
	TH1F *h_tof_fcal_dt;

  TH2F *h_deltaE_tagh,            *h_deltaE_tagm;
  TH2F *h_deltaE_tagh_nomatch,    *h_deltaE_tagm_nomatch;
  TH2F *h_deltaE_tagh_tofmatch,   *h_deltaE_tagm_tofmatch;
	TH2F *h_deltaPhi_tagh,          *h_deltaPhi_tagm;
  TH2F *h_deltaPhi_tagh_nomatch,  *h_deltaPhi_tagm_nomatch;
  TH2F *h_deltaPhi_tagh_tofmatch, *h_deltaPhi_tagm_tofmatch;
  TH2F *h_deltaK_tagh,            *h_deltaK_tagm;
  TH2F *h_deltaK_tagh_nomatch,    *h_deltaK_tagm_nomatch;
  TH2F *h_deltaK_tagh_tofmatch,   *h_deltaK_tagm_tofmatch;
  
  TH2F *h_deltaK_match_vs_energy, *h_deltaK_nomatch_vs_energy;
  TH2F *h_deltaK_match_vs_theta,  *h_deltaK_nomatch_vs_theta;
	
	TH2F *h_deltaE_vs_deltaK_unique;
	TH2F *h_deltaE_vs_deltaK_unique_nomatch;
	TH2F *h_deltaE_vs_deltaK_unique_tofmatch;
	
  TH2F *h_deltaE_vs_deltaK,          *h_deltaE_vs_deltaK_e,          *h_deltaE_vs_deltaK_g;
  TH2F *h_deltaE_vs_deltaK_nomatch,  *h_deltaE_vs_deltaK_e_nomatch,  *h_deltaE_vs_deltaK_g_nomatch;
  TH2F *h_deltaE_vs_deltaK_tofmatch, *h_deltaE_vs_deltaK_e_tofmatch, *h_deltaE_vs_deltaK_g_tofmatch;

  TH2F *h_deltaK_e_vs_deltaK, *h_deltaK_g_vs_deltaK;
	
	TH2F *h_elas_vs_deltaE, *h_elas_vs_deltaE_nomatch, *h_elas_vs_deltaE_tofmatch;
	TH2F *h_mgg_vs_deltaE, *h_mgg_vs_deltaE_nomatch, *h_mgg_vs_deltaE_tofmatch;
	
	//-----------------------------------------//
	// Vary width of DeltaE cut:
	
	static const int n_hists_deltaE = 14;
	double deltaE_cuts[n_hists_deltaE] = {
		1.0, 2.0, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 9.0, 10.0
	};
	TH2F *h_deltaK_tagh_sigE_nomatch[n_hists_deltaE];
	TH2F *h_deltaK_tagh_cut_sigE_nomatch[n_hists_deltaE];
	TH2F *h_deltaK_tagm_sigE_nomatch[n_hists_deltaE];
	TH2F *h_deltaK_tagm_cut_sigE_nomatch[n_hists_deltaE];
  
	TH2F *h_deltaK_tagh_sigE_tofmatch[n_hists_deltaE];
	TH2F *h_deltaK_tagh_cut_sigE_tofmatch[n_hists_deltaE];
	TH2F *h_deltaK_tagm_sigE_tofmatch[n_hists_deltaE];
	TH2F *h_deltaK_tagm_cut_sigE_tofmatch[n_hists_deltaE];
  
	TH2F *h_deltaK_tagh_sigE_nomatch_cut[n_hists_deltaE];
	TH2F *h_deltaK_tagh_cut_sigE_nomatch_cut[n_hists_deltaE];
	TH2F *h_deltaK_tagm_sigE_nomatch_cut[n_hists_deltaE];
	TH2F *h_deltaK_tagm_cut_sigE_nomatch_cut[n_hists_deltaE];
  
	TH2F *h_deltaK_tagh_sigE_tofmatch_cut[n_hists_deltaE];
	TH2F *h_deltaK_tagh_cut_sigE_tofmatch_cut[n_hists_deltaE];
	TH2F *h_deltaK_tagm_sigE_tofmatch_cut[n_hists_deltaE];
	TH2F *h_deltaK_tagm_cut_sigE_tofmatch_cut[n_hists_deltaE];
	
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

