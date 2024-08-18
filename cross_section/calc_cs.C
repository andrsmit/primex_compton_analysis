
void calc_cs(int tag_sys, int counter, double &loc_cs, double &loc_csE) {
	
	loc_cs = 0., loc_csE = 0.;
	
	double loc_eb;
	double loc_yield,  loc_flux,  loc_acc,  loc_fabs;
	double loc_yieldE, loc_fluxE, loc_accE, loc_fabsE;
	
	if(tag_sys==0) {
		
		loc_eb     = tagh_en[counter-1];
		loc_yield  = tagh_yield[counter-1];
		loc_yieldE = tagh_yieldE[counter-1];
		loc_flux   = tagh_flux[counter-1];
		loc_fluxE  = tagh_fluxE[counter-1];
		if(USE_F_ACC) loc_acc = f_acc->Eval(loc_eb);
		else {
			double test_acc = tagh_acc[counter-1];
			if(test_acc==0.) loc_acc = f_acc->Eval(tagh_en[counter-1]);
			else loc_acc = test_acc;
		}
	} else {
		
		loc_eb     = tagm_en[counter-1];
		loc_yield  = tagm_yield[counter-1];
		loc_yieldE = tagm_yieldE[counter-1];
		loc_flux   = tagm_flux[counter-1];
		loc_fluxE  = tagm_fluxE[counter-1];
		if(USE_F_ACC) loc_acc = f_acc->Eval(loc_eb);
		else {
			double test_acc = tagm_acc[counter-1];
			if(test_acc==0.) loc_acc = f_acc->Eval(tagm_en[counter-1]);
			else loc_acc = test_acc;
		}
	}
	
	if(loc_flux<=0. || loc_acc<=0. || n_e<=0. || mb<=0. || loc_yield<=0.) {
		loc_cs  = 0.;
		loc_csE = 0.;
		return;
	}
	
	loc_accE = 0.;
	if(loc_acc <= 0. || loc_flux <= 0.) {
		loc_cs  = 0.;
		loc_csE = 0.;
	} else {
		loc_cs  = loc_yield / (loc_flux * loc_acc * n_e * mb);
		loc_csE = sqrt(
			pow(loc_yieldE / (loc_flux * loc_acc * n_e * mb), 2.0) + 
			pow(loc_fluxE * loc_yield / (loc_flux * loc_flux * loc_acc * n_e * mb), 2.0)
		);
	}
	
	loc_cs  /= f_abs;
	loc_csE /= f_abs;
	
	return;
}
