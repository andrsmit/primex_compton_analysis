
void plot_YieldRatio()
{
	
	int tag_sys = 0;
	
	const bool SAVE_PLOTS = true;
	
	int MAX_COUNTERS;
	if( tag_sys==0 ) MAX_COUNTERS=274;
	else             MAX_COUNTERS=102;
	
	
	vector<int> sigma_vec = {10, 20, 30, 35, 40, 45, 50, 55, 60, 70, 80};
	
	int n_ens = (int)sigma_vec.size();
	
	double tagh_en[274];
	double tagh_yield[n_ens][274];
	double tagh_yieldE[n_ens][274];
	double tagh_acc[n_ens][274];
	double tagh_accE[n_ens][274];
	double tagh_cs[n_ens][274];
	double tagh_csE[n_ens][274];
	
	double tagm_en[102];
	double tagm_yield[n_ens][102];
	double tagm_yieldE[n_ens][102];
	double tagm_acc[n_ens][102];
	double tagm_accE[n_ens][102];
	double tagm_cs[n_ens][102];
	double tagm_csE[n_ens][102];
	
	for( int ie = 0; ie < n_ens; ie++ ) {
		
		int a; double b, c, d, e, f, g, h;
		
		ifstream inf_tagh( Form("tagh_sigE_%02d.txt",sigma_vec[ie]) );
		for( int ic=0; ic<274; ic++ ) {
			inf_tagh >> a >> b >> c >> d >> e >> f >> g >> h;
			if( ie==0 ) tagh_en[ic] = b;
			tagh_cs[ie][ic]     = c;
			tagh_csE[ie][ic]    = d;
			tagh_yield[ie][ic]  = e;
			tagh_yieldE[ie][ic] = f;
			tagh_acc[ie][ic]    = g;
			tagh_accE[ie][ic]   = h;
		}
		inf_tagh.close();
		
		ifstream inf_tagm( Form("tagm_sigE_%02d.txt",sigma_vec[ie]) );
		for( int ic=0; ic<102; ic++ ) {
			inf_tagm >> a >> b >> c >> d >> e >> f >> g >> h;
			if( ie==0 ) tagm_en[ic] = b;
			tagm_cs[ie][ic]     = c;
			tagm_csE[ie][ic]    = d;
			tagm_yield[ie][ic]  = e;
			tagm_yieldE[ie][ic] = f;
			tagm_acc[ie][ic]    = g;
			tagm_accE[ie][ic]   = h;
		}
		inf_tagm.close();
		
	}
	
	
	
	
	TCanvas *canvas1 = new TCanvas( "canvas1", "canvas1", 1000, 600 );
	canvas1->SetTickx(); canvas1->SetTicky();
	
	TCanvas *canvas2 = new TCanvas( "canvas2", "canvas2", 1000, 600 );
	canvas2->SetTickx(); canvas2->SetTicky();
	
	
	
	int loc_counter = 0;
	
	for( int ic = 1; ic < MAX_COUNTERS; ic++ ) 
	{
		
		
		if( tag_sys==0 ) {
			if( tagh_yield[6][ic-1] <= 0. ) continue;
			if(   tagh_cs[6][ic-1]  <= 0. ) continue;
		} else {
			if( tagm_yield[6][ic-1] <= 0. ) continue;
			if(   tagm_cs[6][ic-1]  <= 0. ) continue;
		}
		
		if( tag_sys==0 ) {
			if( tagh_yield[0][ic-1] <= 0. ) continue;
		} else {
			if( tagm_yield[0][ic-1] <= 0. ) continue;
		}
		
		
		if( tag_sys==0 ) cout << "TAGH Counter " << ic << endl;
		else             cout << "TAGM Counter " << ic << endl;
		
		
		double *sigmas       = new double[n_ens];
		double *zeros        = new double[n_ens];
		double *yield_ratio  = new double[n_ens];
		double *yield_ratioE = new double[n_ens];
		double *acc_ratio    = new double[n_ens];
		double *acc_ratioE   = new double[n_ens];
		
		double *cs  = new double[n_ens];
		double *csE = new double[n_ens];
		
		for( int ie = 0; ie < n_ens; ie++ ) {
			
			sigmas[ie]      = 0.1 * (double)sigma_vec[ie];
			zeros[ie]       = 0.;
			
			double loc_yield_ratio, loc_yield_ratioE;
			double loc_acc_ratio, loc_acc_ratioE;
			double loc_cs, loc_csE;
			
			if( !tag_sys ) {
				loc_yield_ratio  = tagh_yield[ie][ic-1] / tagh_yield[0][ic-1];
				loc_yield_ratioE = sqrt( 
					pow( tagh_yieldE[ie][ic-1]/tagh_yield[0][ic-1], 2.0) + 
					pow( tagh_yield[ie][ic-1]*tagh_yieldE[0][ic-1]
					/(tagh_yield[0][ic-1]*tagh_yield[0][ic-1]), 2.0) );
				loc_acc_ratio    = tagh_acc[ie][ic-1] / tagh_acc[0][ic-1];
				loc_acc_ratioE   = sqrt( 
					pow( tagh_accE[ie][ic-1]/tagh_acc[0][ic-1], 2.0) + 
					pow( tagh_acc[ie][ic-1]*tagh_accE[0][ic-1]
					/(tagh_acc[0][ic-1]*tagh_acc[0][ic-1]), 2.0) );
				loc_cs  = tagh_cs[ie][ic-1] / tagh_cs[6][ic-1];
				loc_csE = tagh_csE[ie][ic-1] / tagh_cs[6][ic-1];
			} else {
				loc_yield_ratio  = tagm_yield[ie][ic-1] / tagm_yield[0][ic-1];
				loc_yield_ratioE = sqrt( 
					pow( tagm_yieldE[ie][ic-1]/tagm_yield[0][ic-1], 2.0) + 
					pow( tagm_yield[ie][ic-1]*tagm_yieldE[0][ic-1]
					/(tagm_yield[0][ic-1]*tagm_yield[0][ic-1]), 2.0) );
				loc_acc_ratio    = tagm_acc[ie][ic-1] / tagm_acc[0][ic-1];
				loc_acc_ratioE   = sqrt( 
					pow( tagm_accE[ie][ic-1]/tagm_acc[0][ic-1], 2.0) + 
					pow( tagm_acc[ie][ic-1]*tagm_accE[0][ic-1]
					/(tagm_acc[0][ic-1]*tagm_acc[0][ic-1]), 2.0) );
				loc_cs  = tagm_cs[ie][ic-1] / tagm_cs[6][ic-1];
				loc_csE = tagm_csE[ie][ic-1] / tagm_cs[6][ic-1];
			}
			
			yield_ratio[ie]  = loc_yield_ratio;
			yield_ratioE[ie] = loc_yield_ratioE;
			
			acc_ratio[ie]    = loc_acc_ratio;
			acc_ratioE[ie]   = loc_acc_ratioE;
			
			cs[ie]  = loc_cs;
			csE[ie] = loc_csE;
			
			
		}
		
		
		TGraphErrors *gYieldRatio = new TGraphErrors( n_ens, sigmas, yield_ratio, 
			zeros, yield_ratioE );
		gYieldRatio->GetXaxis()->SetTitle( "#DeltaE Cut Width [#sigma]" );
		gYieldRatio->GetYaxis()->SetTitle( "Relative Yield" );
		if( !tag_sys ) {
			gYieldRatio->SetTitle( Form("TAGH Counter %d (E_{#gamma} = %.4f GeV)", 
				ic, tagh_en[ic-1]) );
		} else {
			gYieldRatio->SetTitle( Form("TAGM Counter %d (E_{#gamma} = %.4f GeV)", 
				ic, tagm_en[ic-1]) );
		}
		gYieldRatio->SetMarkerStyle(23);
		gYieldRatio->SetMarkerColor(kBlue);
		
		
		TGraphErrors *gAccRatio = new TGraphErrors( n_ens, sigmas, acc_ratio, 
			zeros, acc_ratioE );
		gAccRatio->GetXaxis()->SetTitle( "#DeltaE Cut Width [#sigma]" );
		gAccRatio->GetYaxis()->SetTitle( "Relative Yield" );
		if( !tag_sys ) {
			gAccRatio->SetTitle( Form("TAGH Counter %d (E_{#gamma} = %.4f GeV)", 
				ic, tagh_en[ic-1]) );
		} else {
			gAccRatio->SetTitle( Form("TAGM Counter %d (E_{#gamma} = %.4f GeV)", 
				ic, tagm_en[ic-1]) );
		}
		gAccRatio->SetMarkerStyle(23);
		gAccRatio->SetMarkerColor(kRed);
		
		double loc_max = 0.;
		for( int ie=0; ie < n_ens; ie++ ) {
			if( yield_ratio[ie] > 2.0 ) continue;
			if( yield_ratio[ie] > loc_max ) loc_max = yield_ratio[ie];
		}
		
		double loc_min = loc_max;
		for( int ie=0; ie < n_ens; ie++ ) {
			if( yield_ratio[ie] <= 0.2 ) continue;
			if( yield_ratio[ie] < loc_min ) loc_min = yield_ratio[ie];
		}
		
		gYieldRatio->GetYaxis()->SetRangeUser( 0.95*loc_min, 1.05*loc_max );
		
		canvas1->cd();
		gYieldRatio->Draw("AP");
		gAccRatio->Draw("P same");
		
		
		TGraphErrors *gCS = new TGraphErrors( n_ens, sigmas, cs, zeros, csE );
		gCS->GetXaxis()->SetTitle( "#DeltaE Cut Width [#sigma]" );
		gCS->GetYaxis()->SetTitle( "Cross Section Difference" );
		if( !tag_sys ) {
			gCS->SetTitle( Form("TAGH Counter %d (E_{#gamma} = %.4f GeV)", 
				ic, tagh_en[ic-1]) );
		} else {
			gCS->SetTitle( Form("TAGM Counter %d (E_{#gamma} = %.4f GeV)", 
				ic, tagm_en[ic-1]) );
		}
		gCS->SetMarkerStyle(23);
		gCS->SetMarkerColor(kBlue);
		
		gCS->GetYaxis()->SetRangeUser( 0.85, 1.15 );
		
		canvas2->cd();
		gCS->Draw("AP");
		
		
		if( SAVE_PLOTS ) {
		if( loc_counter==0 ) {
			if( tag_sys==0 ) {
				canvas1->Print( "tagh_deltaE_yieldRatio.pdf(",   "pdf" );
				canvas2->Print( "tagh_deltaE_crossSection.pdf(", "pdf" );
			} else {
				canvas1->Print( "tagm_deltaE_yieldRatio.pdf(",   "pdf" );
				canvas2->Print( "tagm_deltaE_crossSection.pdf(",  "pdf" );
			}
		} else {
			if( tag_sys==0 ) {
				canvas1->Print( "tagh_deltaE_yieldRatio.pdf",    "pdf" );
				canvas2->Print( "tagh_deltaE_crossSection.pdf",  "pdf" );
			} else {
				canvas1->Print( "tagm_deltaE_yieldRatio.pdf",    "pdf" );
				canvas2->Print( "tagm_deltaE_crossSection.pdf",  "pdf" );
			}
		}
		}
		
		loc_counter++;
		
	}
	
	if( SAVE_PLOTS ) {
	if( tag_sys==0 ) {
		canvas1->Print( "tagh_deltaE_yieldRatio.pdf)",   "pdf" );
		canvas2->Print( "tagh_deltaE_crossSection.pdf)", "pdf" );
	} else {
		canvas1->Print( "tagm_deltaE_yieldRatio.pdf)",   "pdf" );
		canvas2->Print( "tagm_deltaE_crossSection.pdf)", "pdf" );
	}
	}
	
	return;
}
