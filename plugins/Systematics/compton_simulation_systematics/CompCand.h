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
	double deltaK2;
	
	double rfTime;
	
	double e1, t1, x1, y1, z1;
	double e2, t2, x2, y2, z2;
	
	int ring;
	int pair_cut;
	int phi_slice;
	
	int fcal_phi_sect;
	int ccal_phi_sect;
	
	int fcal_layer;
	int ccal_layer;
	
	double eb;
	int tag_counter;
	int tag_sys;
	
	double vz;
	double event_weight;
	
} ComptonCandidate_t;

#endif
