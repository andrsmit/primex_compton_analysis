#ifndef _PHASE1ENERGY_
#define _PHASE1ENERGY_

#include "compton_cs.h"

void get_phase1_energy_bin(double eb, int &old_tag_counter, int &old_tag_sys) {
	
	double loc_min_e = 12.0;
	old_tag_sys      = -1;
	old_tag_counter  =  0;
	
	double old_eb = 0.;
	
	for(int tagh_counter=1; tagh_counter<=274; tagh_counter++) {
		if(tagh_counter>126 && tagh_counter<179) continue;
		double loc_deltaE = fabs(tagh_en_phase1[tagh_counter-1] - eb);
		if(loc_deltaE < loc_min_e) {
			loc_min_e = loc_deltaE;
			old_tag_sys = 0;
			old_tag_counter = tagh_counter;
			old_eb = tagh_en_phase1[tagh_counter-1];
		}
	}
	for(int tagm_counter=1; tagm_counter<=102; tagm_counter++) {
		double loc_deltaE = fabs(tagm_en_phase1[tagm_counter-1] - eb);
		if(loc_deltaE < loc_min_e) {
			loc_min_e = loc_deltaE;
			old_tag_sys = 1;
			old_tag_counter = tagm_counter;
			old_eb = tagm_en_phase1[tagm_counter-1];
		}
	}
	
	return;
}

#endif
