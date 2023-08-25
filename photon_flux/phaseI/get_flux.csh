#!/usr/bin/tcsh
#

set prog = "/work/halld/home/andrsmit/primex_compton_analysis/photon_flux/primex_flux.py"

# Be target + empty with FIELD OFF:
set runlist = "Be_target_runs.dat"

foreach run ( `cat ${runlist}` )
	echo " "
	echo "$run"
	python ${prog} -b ${run} -e ${run}
	mv ${run}_flux.root rootFiles/
	mv ${run}_tagh_ps_acc_cor.txt txtFiles/
	mv ${run}_tagm_ps_acc_cor.txt txtFiles/
end

