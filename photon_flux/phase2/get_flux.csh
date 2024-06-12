#!/usr/bin/tcsh
#
set prog = "/work/halld/home/andrsmit/primex_compton_analysis/photon_flux/primex_flux.py"

source $BUILD_SCRIPTS/gluex_env_jlab.csh $HOME/version_h2root.xml

set run_list_dir = /work/halld/home/andrsmit/run_lists

set run_lists = {"phase2_be_target_runs","phase2_be_empty_runs","phase2_he_target_runs_fieldoff","phase2_he_empty_runs_fieldoff","phase2_he_target_runs_fieldon","phase2_he_empty_runs_fieldon"}

foreach list ($run_lists) 
	foreach run ( `cat $run_list_dir/$list.dat` )
		echo " "
		echo "$run"
		python ${prog} -b ${run} -e ${run}
		mv ${run}_flux.root rootFiles_new/
		mv ${run}_tagh_ps_acc_cor.txt txtFiles/
		mv ${run}_tagm_ps_acc_cor.txt txtFiles/
	end
end
