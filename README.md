# primex_compton_analysis
Code used for the Compton Cross section measurement in PrimEx-Eta experiment.

Plugins:

1) compton_analysis

2) compton_analysis2
	- Same as 'compton_analysis', but with different timing cuts applied. In the first plugin,
	the FCAL and CCAL showers were required to be within +/- 3 ns of the DEventRFBunch time and the
	DBeamPhoton was required to be within +/- 2.004 ns of the DEventRFBunch time. 
	However, in the data there are tails in the fcal-rf timing distribution as well as tagm-rf 
	timing and ccal-rf timing distributions that aren't present in the simulation. This means 
	we might be excluding events from analysis which are not also excluded when calculating 
	acceptance.
	
	In this plugin we don't place a cut on the fcal-rf time. Instead we can place a cut on the 
	CCAL-RF time and select FCAL showers based on a cut between FCAL-CCAL time. 

3) compton_analysis3
	- Same as compton_analysis2, but with CCAL-RF Cut set to +/- 6.012 ns instead of +/- 3 ns.

4) compton_analysis4
	- Same as compton_analysis2, but with including 3 rf-beam bunches in the main bunch (still 
	using 2 sidebands for accidental subtraction, though). 
