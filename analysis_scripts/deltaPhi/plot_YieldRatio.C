
void plot_YieldRatio()
{
	
	
	int tag_sys =   1;
	int counter =  44;
	
	if( counter < 1 || counter > 274 ) return;
	if( tag_sys ==1 && counter > 102 ) return;
	
	
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
		
		ifstream inf_tagh( Form("tagh_sigPhi_%02d.txt",sigma_vec[ie]) );
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
		
		ifstream inf_tagm( Form("tagm_sigPhi_%02d.txt",sigma_vec[ie]) );
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
	
	
	
	double *sigmas       = new double[n_ens];
	double *zeros        = new double[n_ens];
	double *yield_ratio  = new double[n_ens];
	double *yield_ratioE = new double[n_ens];
	double *acc_ratio    = new double[n_ens];
	double *acc_ratioE   = new double[n_ens];
	
	double *cs  = new double[n_ens];
	double *csE = new double[n_ens];
	
	for( int ie = 0; ie < n_ens; ie++ ) {
		
		sigmas[ie]  = (double)sigma_vec[ie];
		zeros[ie]   = 0.;
		
		double loc_yield_ratio, loc_yield_ratioE;
		double loc_acc_ratio, loc_acc_ratioE;
		double loc_cs, loc_csE;
		
		if( !tag_sys ) {
			loc_yield_ratio  = tagh_yield[ie][counter-1] / tagh_yield[0][counter-1];
			loc_yield_ratioE = sqrt( 
				pow( tagh_yieldE[ie][counter-1]/tagh_yield[0][counter-1], 2.0) + 
				pow( tagh_yield[ie][counter-1]*tagh_yieldE[0][counter-1]
				/(tagh_yield[0][counter-1]*tagh_yield[0][counter-1]), 2.0) );
			loc_acc_ratio    = tagh_acc[ie][counter-1] / tagh_acc[0][counter-1];
			loc_acc_ratioE   = sqrt( 
				pow( tagh_accE[ie][counter-1]/tagh_acc[0][counter-1], 2.0) + 
				pow( tagh_acc[ie][counter-1]*tagh_accE[0][counter-1]
				/(tagh_acc[0][counter-1]*tagh_acc[0][counter-1]), 2.0) );
			loc_cs  = tagh_cs[ie][counter-1];
			loc_csE = tagh_csE[ie][counter-1];
		} else {
			loc_yield_ratio  = tagm_yield[ie][counter-1] / tagm_yield[0][counter-1];
			loc_yield_ratioE = sqrt( 
				pow( tagm_yieldE[ie][counter-1]/tagm_yield[0][counter-1], 2.0) + 
				pow( tagm_yield[ie][counter-1]*tagm_yieldE[0][counter-1]
				/(tagm_yield[0][counter-1]*tagm_yield[0][counter-1]), 2.0) );
			loc_acc_ratio    = tagm_acc[ie][counter-1] / tagm_acc[0][counter-1];
			loc_acc_ratioE   = sqrt( 
				pow( tagm_accE[ie][counter-1]/tagm_acc[0][counter-1], 2.0) + 
				pow( tagm_acc[ie][counter-1]*tagm_accE[0][counter-1]
				/(tagm_acc[0][counter-1]*tagm_acc[0][counter-1]), 2.0) );
			loc_cs  = tagm_cs[ie][counter-1];
			loc_csE = tagm_csE[ie][counter-1];
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
	gYieldRatio->GetXaxis()->SetTitle( "#Delta#phi Cut Width [#sigma]" );
	gYieldRatio->GetYaxis()->SetTitle( "Relative Yield" );
	if( !tag_sys ) {
		gYieldRatio->SetTitle( Form("TAGH Counter %d (E_{#gamma} = %.4f GeV)", 
			counter, tagh_en[counter-1]) );
	} else {
		gYieldRatio->SetTitle( Form("TAGM Counter %d (E_{#gamma} = %.4f GeV)", 
			counter, tagm_en[counter-1]) );
	}
	gYieldRatio->SetMarkerStyle(23);
	gYieldRatio->SetMarkerColor(kBlue);
	
	TGraphErrors *gAccRatio = new TGraphErrors( n_ens, sigmas, acc_ratio, 
		zeros, acc_ratioE );
	gAccRatio->GetXaxis()->SetTitle( "#Delta#phi Cut Width [#sigma]" );
	gAccRatio->GetYaxis()->SetTitle( "Relative Yield" );
	if( !tag_sys ) {
		gAccRatio->SetTitle( Form("TAGH Counter %d (E_{#gamma} = %.4f GeV)", 
			counter, tagh_en[counter-1]) );
	} else {
		gAccRatio->SetTitle( Form("TAGM Counter %d (E_{#gamma} = %.4f GeV)", 
			counter, tagm_en[counter-1]) );
	}
	gAccRatio->SetMarkerStyle(23);
	gAccRatio->SetMarkerColor(kRed);
	
	TCanvas *canvas = new TCanvas( "canvas", "canvas", 1000, 600 );
	canvas->SetTickx(); canvas->SetTicky();
	gYieldRatio->Draw("AP");
	gAccRatio->Draw("P same");
	
	
	
	TGraphErrors *gCS = new TGraphErrors( n_ens, sigmas, cs, zeros, csE );
	gCS->GetXaxis()->SetTitle( "#Delta#phi Cut Width [#sigma]" );
	gCS->GetYaxis()->SetTitle( "Cross Section [mb / electron]" );
	if( !tag_sys ) {
		gCS->SetTitle( Form("TAGH Counter %d (E_{#gamma} = %.4f GeV)", 
			counter, tagh_en[counter-1]) );
	} else {
		gCS->SetTitle( Form("TAGM Counter %d (E_{#gamma} = %.4f GeV)", 
			counter, tagm_en[counter-1]) );
	}
	gCS->SetMarkerStyle(23);
	gCS->SetMarkerColor(kBlue);
	
	TCanvas *canvas2 = new TCanvas( "canvas2", "canvas2", 1000, 600 );
	canvas2->SetTickx(); canvas2->SetTicky();
	gCS->Draw("AP");
	
	return;
}
