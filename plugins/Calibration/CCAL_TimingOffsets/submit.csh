#!/usr/bin/tcsh
#
set dir_mss = /mss/halld/RunPeriod-2019-01/rawdata
set phase   = 1
set outdir  = /work/halld/home/andrsmit/primex_compton_analysis/plugins/Calibration/CCAL_TimingOffsets
set script  = ${outdir}/run_rec.csh
set ij      = 0

set submit_runs = go

#=========================#
#Job Resources:

set workflow   = compton_ana_tof_veto
set account    = halld
set partition  = production
set ram        = 12GB
set cores      = 4
set time       = 240min
set disk       = 23GB
set constraint = el9

#=========================#

set run_lists = /work/halld/home/andrsmit/run_lists

#foreach run ( `cat $run_lists/phase1_be_target_runs.dat $run_lists/phase1_be_empty_runs.dat` )
set run = 61321
	
	# Convert run number into 6 digit variable:
	
	set runnumber = $run
	if ($runnumber<"100000") then
		set runnumber = "0$run"
	endif
	echo "run ${runnumber}"
	
	# Find the appropriate directory on tape for this run number:
	
	if ( $runnumber>"60000" && $runnumber<"69999" ) then
		set phase = 1
		set dir_mss = /mss/halld/RunPeriod-2019-01/rawdata/Run${runnumber}
	else if ( $runnumber>"80000" && $runnumber<"89999" ) then
		set phase = 2
		set dir_mss = /mss/halld/RunPeriod-2021-08/rawdata/Run${runnumber}
	else if ( $runnumber>"110000" && $runnumber<"119999" ) then
		set phase = 3
		set dir_mss = /mss/halld/RunPeriod-2022-08/rawdata/Run${runnumber}
	else 
		continue
	endif
	
	# Make sure the output directory exists for this run number:
	
	set loc_outdir = $outdir/rootFiles/$runnumber
	if( ! -d $loc_outdir ) then
		echo "Output directory: $loc_outdir does not exist."
		continue
	endif
	
	# Loop over file extensions:
	set iext = 0
	while ($iext < 250)
		
		# Convert extension into 3 digit variable:
		
		set ext = "$iext"
		if ( $iext<"10" ) then
			set ext = "00$iext"
		else if ( $iext<"100" ) then
			set ext = "0$iext"
		endif
		
		# Set up file name of output ROOT file:
		
		set of_name = $loc_outdir/${runnumber}_${ext}.root
		
		# If this file has already been processed, skip it:
		
		if( -f $of_name ) then
			set iext = `expr $iext + 1`
			continue
		endif
		
		# Set up jsub file name:
		
		set jf_name = $outdir/jsub/${run}_${ext}.jsub
		
		# If this file is currently being processed, skip it:
		
		if( -f $jf_name ) then
		set iext = `expr $iext + 1`
			continue
		endif
		
		# Set raw data file name:
		
		set loc_fname = hd_rawdata_${runnumber}_${ext}.evio
		set mss_fname = ${dir_mss}/${loc_fname}
		
		# Check if file exists on tape:
		
		if( -f $mss_fname ) then
			
			echo "  ${runnumber}_${ext}"
			
			if($submit_runs == "go") then
				set command = "swif2 add-job -workflow ${workflow}"
				set command = "$command -name ccal_hit_tw_${runnumber}_${ext}"
				set command = "$command -account ${account} -partition ${partition}"
				set command = "$command -cores ${cores} -ram ${ram} -time ${time} -disk ${disk}"
				set command = "$command -constraint ${constraint}"
				set command = "$command -input $loc_fname mss:$mss_fname"
				set command = "$command $script $loc_outdir $runnumber $ext $mss_fname $loc_fname $jf_name"
				
				echo "$command" > $jf_name
				
				$command
				
			endif
			set ij = `expr $ij + 1`
			
		endif
		
		#*******************************************************************#
		
		set iext = `expr $iext + 1`
		
	end #loop over file extensions
#end #loop over run numbers

echo "submitted $ij jobs"
