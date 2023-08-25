
void plot_weight()
{
	TRandom3 *rand = new TRandom3(0);
	
	TF1 *f1 = new TF1("f1","exp(-0.0082*0.122*x)", 0.01, 30.);
	
	f1->Draw();
	
	/*
	TH1F *h_vertex = new TH1F("vertex", "Vertex Z position", 300, 0., 30.);
	
	for(int i=0; i<1000000; i++) {
		double x = 0.;
		int found_val = 0;
		while(x<29.5) {
			x += 0.05;
			double P = f1->Eval(x);
			double pp = rand->Uniform(0.,1.);
			if(P>pp) {
				found_val++;
				break;
			}
		}
		if(found_val) {
			h_vertex->Fill(x);
		}
	}
	
	h_vertex->Draw();
	*/
	return;
}
