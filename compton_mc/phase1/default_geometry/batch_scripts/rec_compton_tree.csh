#!/usr/bin/tcsh
#

# command line inputs:

set tag_sys       = $1
set tag_counter   = $2
set plugin_name   = $3
set plugin_opts   = $4
set output_subdir = $5
set writeDir      = $6
set jfile         = $7

#------------------------------------------------------------------------------------------#
#Source setup script:

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.csh
gxenv /home/andrsmit/version.xml

setenv CCDB_CONNECTION mysql://ccdb_user@hallddb-farm.jlab.org/ccdb
setenv JANA_CALIB_URL  mysql://ccdb_user@hallddb-farm.jlab.org/ccdb

# use mc ccdb variation:
setenv JANA_CALIB_CONTEXT "variation=mc"

set bin = ${HALLD_RECON_HOME}/Linux_Alma9-x86_64-gcc11.4.1/bin/hd_root

echo "Environment variables set... Now running hd_root..."

#****************************************************************#
# Analyze simulation with random background:

if ( -f $writeDir/recRootTrees_new/${tag_sys}_${tag_counter}.root ) then
	echo "Output ROOT file already exists - nothing else to do."
else
	
	if ( -f $writeDir/recHddmFiles_new/${tag_sys}_${tag_counter}.hddm ) then
		
		echo "cp $writeDir/recHddmFiles_new/${tag_sys}_${tag_counter}.hddm output.hddm"
		cp $writeDir/recHddmFiles_new/${tag_sys}_${tag_counter}.hddm output.hddm
		
		echo "$bin -PPLUGINS=${plugin_name} ${plugin_opts} output.hddm"
		$bin -PPLUGINS=${plugin_name} ${plugin_opts} output.hddm
		
		echo "mv primex_compton.root ${writeDir}/recRootTrees_new/${tag_sys}_${tag_counter}.root"
		mv primex_compton.root ${writeDir}/recRootTrees_new/${tag_sys}_${tag_counter}.root
		
		#echo "mv hd_root.root ${writeDir}/${output_subdir}/${tag_sys}_${tag_counter}.root"
		#mv hd_root.root ${writeDir}/${output_subdir}/${tag_sys}_${tag_counter}.root
		
		echo "rm hd_root.root"
		rm hd_root.root
		
		echo "rm output.hddm"
		rm output.hddm
		
	endif
	
	#--------------------------------------------------------------#
	# Now do the same for the hddm files without random background:
	
	#if ( -f $writeDir/recHddmFiles_noBkgd/${tag_sys}_${tag_counter}.hddm ) then
	#	
	#	echo "cp $writeDir/recHddmFiles_noBkgd/${tag_sys}_${tag_counter}.hddm output.hddm"
	#	cp $writeDir/recHddmFiles_noBkgd/${tag_sys}_${tag_counter}.hddm output.hddm
	#	
	#	echo "$bin -PPLUGINS=${plugin_name} ${plugin_opts} output.hddm"
	#	$bin -PPLUGINS=${plugin_name} ${plugin_opts} output.hddm
	#	
	#	echo "mv primex_compton.root ${writeDir}/recRootTrees_noBkgd/${tag_sys}_${tag_counter}.root"
	#	mv primex_compton.root ${writeDir}/recRootTrees_noBkgd/${tag_sys}_${tag_counter}.root
	#	
	#	#echo "mv hd_root.root ${writeDir}/${output_subdir}/${tag_sys}_${tag_counter}.root"
	#	#mv hd_root.root ${writeDir}/${output_subdir}/${tag_sys}_${tag_counter}.root
	#	
	#	echo "rm hd_root.root"
	#	rm hd_root.root
	#	
	#	echo "rm output.hddm"
	#	rm output.hddm
	#	
	#endif
	
endif

#****************************************************************#
# clean local directory:

echo "rm $jfile"
rm $jfile
