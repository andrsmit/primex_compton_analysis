#ifndef _CompCand_
#define _CompCand_

typedef struct {
  
  int bunch_val;
  
  double deltaT;	
  double deltaPhi;
  double dTheta;
  double deltaR;
  double deltaE;
  double deltaK;
  double ecomp1, ecomp2;
	double invmass;
	
  double e1;
	double t1;
  double x1;
  double y1;
  double z1;
  
  double e2;
	double t2;
  double x2;
  double y2;
  double z2;
  
  int ring;
  int pair_cut;
  int phi_slice;
  int unique_val;
	
  double eb, brfdt;
  int tag_counter;
  int tag_sys;
  
  int tof_match;
  
  int is_mc;
  double event_weight;
  
} ComptonCandidate_t;

#endif
