#ifndef _ComptonAna_
#define _ComptonAna_

#define MAX_BEAM 1000
#define MAX_FCAL 100
#define MAX_CCAL 100
#define MAX_TOF  200
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
		
		int m_runNumber;
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
		void comptonAnalysis_systematics();
		
		int cut_deltaE(  double deltaE,   double eb, double n_sigma_left, double n_sigma_right);
		int cut_deltaPhi(double deltaPhi, double eb, double n_sigma_left, double n_sigma_right);
		int cut_deltaK(  double deltaK,   double eb, double n_sigma_left, double n_sigma_right);
		
		int cut_deltaE_two(  double deltaE,   double eb, double n_sigma_left, double n_sigma_right);
		int cut_deltaPhi_two(double deltaPhi, double eb, double n_sigma_left, double n_sigma_right);
		int cut_deltaK_two(  double deltaK,   double eb, double n_sigma_left, double n_sigma_right);
		
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
		int ccal_fiducial_cut(TVector3 pos, double cut_layer=1.0);
		
		void check_TOF_match(TVector3 pos, double &dx_min, double &dy_min, double &dz_min, double rf_time_cut);
		
		// TTree Variables:
		int    m_eventNum;
		double m_rfTime;
		int    m_nbeam;
		int    m_tag_sys[MAX_BEAM];
		int    m_tag_counter[MAX_BEAM];
		double m_beam_e[MAX_BEAM];
		double m_beam_t[MAX_BEAM];
		double m_acc_scale_factor[MAX_BEAM];
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
		double m_mc_pdgtype[MAX_MC];
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
		
		TH1F *h_thrown_photon_angle;
		TH1F *h_thrown_photon_angle_accepted;
		
		TH1F *h_reaction_weight, *h_reaction_weight_double;
		TH1F *h_vertex, *h_vertex_accepted;
		TH1F *h_thrown_mass;
		TH2F *h_deltaE_vs_thrown_mass, *h_deltaE_vs_rec_mass;
		TH2F *h_deltaK_vs_thrown_mass, *h_deltaK_vs_rec_mass;
		TH2F *h_elas_low_mass, *h_elas_high_mass;
		TH2F *h_rec_vs_thrown_mass;
		
		TH1F *h_fcal_rf_dt, *h_ccal_rf_dt, *h_beam_rf_dt, *h_beam_rf_dt_cut;
		
		TH2F *h_deltaE_pair[6], *h_deltaK_pair[6], *h_deltaPhi_pair[6];
		TH2F *h_deltaE_triplet[6], *h_deltaK_triplet[6], *h_deltaPhi_triplet[6];
		
		static const int m_n_cuts = 12;
		
		TH2F   *h_deltaE_tagh[m_n_cuts],   *h_deltaE_tagm[m_n_cuts];
		TH2F *h_deltaPhi_tagh[m_n_cuts], *h_deltaPhi_tagm[m_n_cuts];
		TH2F   *h_deltaK_tagh[m_n_cuts],   *h_deltaK_tagm[m_n_cuts];
		
		TH2F         *h_elas_vs_deltaE[m_n_cuts];
		TH2F       *h_deltaK_vs_deltaE[m_n_cuts];
		TH2F    *h_deltaCCAL_vs_deltaE[m_n_cuts],        *h_deltaCCAL_vs_deltaK[m_n_cuts];
		TH2F    *h_deltaFCAL_vs_deltaE[m_n_cuts],        *h_deltaFCAL_vs_deltaK[m_n_cuts];
		TH2F *h_deltaFCAL_vs_deltaCCAL[m_n_cuts], *h_deltaFCAL_vs_deltaCCAL_cut[m_n_cuts];
		
		TH2F      *h_ccal_xy[m_n_cuts], *h_fcal_xy[m_n_cuts];
		TH2F *h_opangle_tagh[m_n_cuts], *h_opangle_tagh_ecut[m_n_cuts], *h_opangle_tagh_ekcut[m_n_cuts];
		TH2F *h_opangle_tagm[m_n_cuts], *h_opangle_tagm_ecut[m_n_cuts], *h_opangle_tagm_ekcut[m_n_cuts];
		
		TH1F *h_ccal_nblocks,       *h_fcal_nblocks;
		TH1F *h_ccal_nblocks_cut,   *h_fcal_nblocks_cut;
		TH2F *h_ccal_xy_highdeltaE, *h_fcal_xy_highdeltaE;
		TH2F *h_ccal_xy_lowdeltaE,  *h_fcal_xy_lowdeltaE;
		
		TH2F *h_sumPhi_vs_deltaPhi;
		
		//----------------------------------------------------------------------------------------------//
		// Systematics:
		
		TH2F *h_mgg_vs_deltaK, *h_mgg_vs_deltaK_cut;
		
		// Vary minimum FCAL shower energy cut:
		
		static const int m_n_hists_fcalE = 20;
		TH2F *h_deltaK_tagh_fcalE[m_n_hists_fcalE], *h_deltaK_tagm_fcalE[m_n_hists_fcalE];
		
		// Vary minimum FCAL shower energy cut:
		
		static const int m_n_hists_ccalE = 13;
		TH2F *h_deltaK_tagh_ccalE[m_n_hists_ccalE], *h_deltaK_tagm_ccalE[m_n_hists_ccalE];
		
		// Vary minimum FCAL shower timing cut:
		
		static const int m_n_hists_fcalT = 12;
		TH2F *h_deltaK_tagh_fcalT[m_n_hists_fcalT], *h_deltaK_tagm_fcalT[m_n_hists_fcalT];
		
		// Vary minimum CCAL shower timing cut:
		
		static const int m_n_hists_ccalT = 12;
		TH2F *h_deltaK_tagh_ccalT[m_n_hists_ccalT], *h_deltaK_tagm_ccalT[m_n_hists_ccalT];
		
		// Vary size of square fiducial beam-hole cuts:
		
		static const int m_n_fid_cuts = 20;
		TH2F *h_deltaK_tagh_fcalfid[m_n_fid_cuts], *h_deltaK_tagm_fcalfid[m_n_fid_cuts];
		TH2F *h_deltaK_tagh_ccalfid[m_n_fid_cuts], *h_deltaK_tagm_ccalfid[m_n_fid_cuts];
		
		// Vary width of DeltaE cut:
		
		static const int m_n_hists_deltaE = 16;
		double m_deltaE_cuts[m_n_hists_deltaE] = {
			1.0, 2.0, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 9.0, 10.0, 15.0, 20.0
		};
		TH2F *h_deltaK_tagh_sigE[m_n_hists_deltaE], *h_deltaK_tagm_sigE[m_n_hists_deltaE];
		
		// Vary width of DeltaPhi cut:
		
		static const int m_n_hists_deltaPhi = 16;
		double m_deltaPhi_cuts[m_n_hists_deltaPhi] = {
			1.0, 2.0, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 9.0, 10.0, 15.0, 20.0
		};
		TH2F *h_deltaK_tagh_sigPhi[m_n_hists_deltaPhi], *h_deltaK_tagm_sigPhi[m_n_hists_deltaPhi];
		
		// Measure cross-section in different regions of FCAL and CCAL:
		
		TH2F *h_deltaK_tagh_fcal_phi[8],   *h_deltaK_tagm_fcal_phi[8];
		TH2F *h_deltaK_tagh_ccal_phi[8],   *h_deltaK_tagm_ccal_phi[8];
		TH2F *h_deltaK_tagh_fcal_layer[8], *h_deltaK_tagm_fcal_layer[8];
		TH2F *h_deltaK_tagh_ccal_layer[5], *h_deltaK_tagm_ccal_layer[5];
		
		// Plot distribution of showers on FCAL and CCAL:
		
		TH2F *h_xy_fcal_phi[8],   *h_xy_ccal_phi[8];
		TH2F *h_xy_fcal_layer[8], *h_xy_ccal_layer[5];
		
		//----------------------------------------------------------------------------------------------//
		
	public:
		
		int m_SMEAR_DISTRIBUTIONS = 0;
		int m_SHIFT_DISTRIBUTIONS = 0;
		
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
		
		int getPrimexPhase(int run_number);
		
		// default analysis:
		void initHistograms();
		void runAnalysis(TString infname);
		void resetHistograms();
		void writeHistograms();
		
		// full systematics analysis:
		void  initHistograms_systematics();
		void     runAnalysis_systematics(TString infname);
		void resetHistograms_systematics();
		void writeHistograms_systematics();
		
		void setRunNumber(int runNum);
		void setOutputFileName(string name);
};

#endif
