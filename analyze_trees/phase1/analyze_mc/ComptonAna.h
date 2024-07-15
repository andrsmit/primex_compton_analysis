#ifndef _ComptonAna_
#define _ComptonAna_

#define MAX_BEAM 200
#define MAX_FCAL 20
#define MAX_CCAL 20
#define MAX_TOF  50
#define MAX_MC   5

using namespace std;

#include <stdio.h>
#include <fstream>
#include <math.h>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>

#include "TVector3.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TSystem.h"
#include "TDirectory.h"
#include "TString.h"
#include "TRandom3.h"

class ComptonAna {
	private:
		
		double m_runNumber;
		TRandom3 *m_random;
		
		// Geometry:
		
		int    m_phase_val;
		double m_target_length, m_target_density, m_target_atten;
		
		TVector3 m_fcal_face,       m_ccal_face;
		TVector3 m_fcal_correction, m_ccal_correction;
		TVector3 m_vertex;
		
		// Cut parameters:
		
		double m_deltaE_mu_pars[2][4],   m_deltaE_sigma_pars[2][4];
		double m_deltaPhi_mu_pars[2][4], m_deltaPhi_sigma_pars[2][4];
		double m_deltaK_mu_pars[2][4],   m_deltaK_sigma_pars[2][4];
		
		// For a second-layer fiducial cut:
		
		double m_deltaE_mu_pars_two[2][4],   m_deltaE_sigma_pars_two[2][4];
		double m_deltaPhi_mu_pars_two[2][4], m_deltaPhi_sigma_pars_two[2][4];
		double m_deltaK_mu_pars_two[2][4],   m_deltaK_sigma_pars_two[2][4];
		
		// Constants:
		
		static constexpr double m_c = TMath::C() * 1.e-7; // [cm/ns]
		static constexpr double m_e = 0.510998928e-3;     // [GeV]
		
		static constexpr double m_fcal_block_size = 4.0157;
		static constexpr double m_ccal_block_size = 2.09;
		
		// Variables:
		
		string m_output_fname;
		
		TFile *m_infile;
		TTree *m_tree;
		
		int m_event;
		
		int loadTree();
		int loadCutParameters();
		void setGeometry();
		void readEvent();
		void comptonAnalysis();
		
		int cut_deltaE(  double deltaE,   double eb, double n_sigma);
		int cut_deltaPhi(double deltaPhi, double eb, double n_sigma);
		int cut_deltaK(  double deltaK,   double eb, double n_sigma);
		
		int cut_deltaE_two(  double deltaE,   double eb, double n_sigma);
		int cut_deltaPhi_two(double deltaPhi, double eb, double n_sigma);
		int cut_deltaK_two(  double deltaK,   double eb, double n_sigma);
		
		double smear_deltaE(double deltaE, double eb);
		double smear_deltaPhi(double deltaPhi, double eb);
		double smear_deltaK(double deltaK, double eb);
		
		double smear_deltaE_two(double deltaE, double eb);
		double smear_deltaPhi_two(double deltaPhi, double eb);
		double smear_deltaK_two(double deltaK, double eb);
		
		int acceptRejectEvent();
		
		TVector3 getFCALPosition(int index);
		TVector3 getCCALPosition(int index);
		
		int fcal_fiducial_cut(TVector3 pos, double cut_layer);
		int ccal_fiducial_cut(TVector3 pos);
		
		void check_TOF_match(TVector3 pos, double &dx_min, double &dy_min, double &dz_min, double rf_time_cut);
		
		// TTree Variables:
		int    m_eventNum;
		double m_rfTime;
		int    m_nbeam;
		int    m_tag_sys[MAX_BEAM];
		int    m_tag_counter[MAX_BEAM];
		double m_beam_e[MAX_BEAM];
		double m_beam_t[MAX_BEAM];
		int    m_nfcal;
		double m_fcal_e[MAX_FCAL];
		double m_fcal_x[MAX_FCAL];
		double m_fcal_y[MAX_FCAL];
		double m_fcal_z[MAX_FCAL];
		double m_fcal_t[MAX_FCAL];
		int    m_fcal_nblocks[MAX_FCAL];
		int    m_nccal;
		double m_ccal_e[MAX_CCAL];
		double m_ccal_x[MAX_CCAL];
		double m_ccal_y[MAX_CCAL];
		double m_ccal_z[MAX_CCAL];
		double m_ccal_t[MAX_CCAL];
		int    m_ccal_idmax[MAX_CCAL];
		int    m_ccal_nblocks[MAX_CCAL];
		int    m_ntof;
		double m_tof_x[MAX_TOF];
		double m_tof_y[MAX_TOF];
		double m_tof_z[MAX_TOF];
		double m_tof_t[MAX_TOF];
		int    m_nmc;
		int    m_mc_pdgtype[MAX_MC];
		double m_mc_x[MAX_MC];
		double m_mc_y[MAX_MC];
		double m_mc_z[MAX_MC];
		double m_mc_t[MAX_MC];
		double m_mc_e[MAX_MC];
		double m_mc_p[MAX_MC];
		double m_mc_theta[MAX_MC];
		double m_mc_phi[MAX_MC];
		int    m_reaction_type;
		double m_reaction_weight;
		double m_reaction_energy;
		
		// Histograms:
		
		TH1F *h_reaction_weight, *h_reaction_weight_double;
		TH1F *h_vertex, *h_vertex_accepted;
		
		TH1F *h_fcal_rf_dt, *h_ccal_rf_dt, *h_beam_rf_dt, *h_beam_rf_dt_cut;
		
		static const int m_n_cuts = 12;
		TH2F *h_deltaE_tagh[m_n_cuts],   *h_deltaE_tagm[m_n_cuts];
		TH2F *h_deltaPhi_tagh[m_n_cuts], *h_deltaPhi_tagm[m_n_cuts];
		TH2F *h_deltaK_tagh[m_n_cuts],   *h_deltaK_tagm[m_n_cuts];
		TH2F *h_elas_vs_deltaE[m_n_cuts];
		TH2F *h_deltaK_vs_deltaE[m_n_cuts];
		TH2F *h_deltaCCAL_vs_deltaE[m_n_cuts], *h_deltaCCAL_vs_deltaK[m_n_cuts];
		TH2F *h_deltaFCAL_vs_deltaE[m_n_cuts], *h_deltaFCAL_vs_deltaK[m_n_cuts];
		TH2F *h_ccal_xy[m_n_cuts],             *h_fcal_xy[m_n_cuts];
		TH2F *h_opangle_tagh[m_n_cuts], *h_opangle_tagh_ecut[m_n_cuts], *h_opangle_tagh_ekcut[m_n_cuts];
		TH2F *h_opangle_tagm[m_n_cuts], *h_opangle_tagm_ecut[m_n_cuts], *h_opangle_tagm_ekcut[m_n_cuts];
		
		TH1F *h_ccal_nblocks,       *h_fcal_nblocks;
		TH1F *h_ccal_nblocks_cut,   *h_fcal_nblocks_cut;
		TH2F *h_ccal_xy_highdeltaE, *h_fcal_xy_highdeltaE;
		TH2F *h_ccal_xy_lowdeltaE,  *h_fcal_xy_lowdeltaE;
		
		static const int m_n_hists_deltaE = 16;
		double m_deltaE_cuts[m_n_hists_deltaE] = {
			1.0, 2.0, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 9.0, 10.0, 15.0, 20.0
		};
		TH2F *h_deltaK_tagh_sigE[m_n_hists_deltaE], *h_deltaK_tagm_sigE[m_n_hists_deltaE];
		
	public:
		
		// Cuts:
		
		double m_cut_fcalE;
		double m_cut_ccalE;
		double m_cut_fcalrfdt;
		double m_cut_ccalrfdt;
		double m_cut_beamrfdt;
		double m_cut_deltaE;
		double m_cut_deltaK;
		double m_cut_deltaPhi;
		
		double m_beam_bunches_main = 1.0;
		double m_beam_bunches_acc  = 5.0; // how many bunches to use on each side for accidental subtraction
		
		ComptonAna();
		~ComptonAna(){};
		
		void initHistograms();
		void resetHistograms();
		void writeHistograms();
		
		void runAnalysis(TString infname);
		
		void setRunNumber(int runNum);
		void setOutputFileName(string name);
};

#endif
