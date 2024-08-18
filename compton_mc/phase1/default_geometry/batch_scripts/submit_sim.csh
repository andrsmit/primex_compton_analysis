#!/usr/bin/tcsh
#
set outdir = /work/halld/home/andrsmit/primex_compton_analysis/compton_mc/phase1/default_geometry
set script = ${outdir}/batch_scripts/sim_compton.csh
set ij     = 0

set submit_runs = go

#=========================#
#Job Resources:

set workflow   = compton_mc_tof_veto
set account    = halld
set partition  = production
set ram        = 8GB
set cores      = 10
set time       = 500min
set disk       = 12GB
set constraint = el9

#=========================#

set run_list = {"61866"}
set tag_list = {"tagh","tagm"}

foreach rnb ($run_list)
	
	set run_number = $rnb
	if ( $run_number < 100000 ) then
		set run_number = "0$rnb"
	endif
	
	set gendir   = /work/halld/home/andrsmit/primex_compton_analysis/compton_mc/genDir/Run${run_number}
	set writedir = ${outdir}/Run${run_number}
	
	foreach tag_sys ($tag_list)
		
		set max_counter = 230
		if("$tag_sys" == "tagm") then
			set max_counter = 102
		endif
		
		# iterate from counters 1-127 + 179-274
		set counter = 1
		while ( $counter <= $max_counter )
			
			#skip the TAGM region in the TAGH:
			if( $counter > 127 && $counter < 179 ) then
				set counter = `expr $counter + 1`
				continue
			endif
			
			#setup the counter number string:
			set tag_counter = $counter
			if ( $counter < 100 ) then
				if ( $counter < 10 ) then
					set tag_counter = 00$counter
				else
					set tag_counter = 0$counter
				endif
			endif
			
			#***********************************************************#
			# submit jobs:
			
			set loc_jobname = compton_sim_rnb-${run_number}_${tag_sys}_${tag_counter}
			
			set jf = ${outdir}/jsub/${loc_jobname}.jsub
			
			# skip jobs that have already been submitted:
			if ( -f $jf ) then
				set counter = `expr $counter + 1`
				continue
			endif
			
			# skip counters which have already been processed:
			if ( -f ${writedir}/recHddmFiles_unsmeared/${tag_sys}_${tag_counter}.hddm ) then
				set counter = `expr $counter + 1`
				continue
			endif
			
			if ( -f ${gendir}/${tag_sys}_${tag_counter}.hddm ) then
				
				echo "${tag_sys}_${tag_counter}"
				
				if($submit_runs == "go") then
					
					set command = "swif2 add-job -workflow ${workflow}"
					set command = "$command -name ${loc_jobname}"
					set command = "$command -account ${account} -partition ${partition}"
					set command = "$command -cores ${cores} -ram ${ram} -time ${time} -disk ${disk}"
					set command = "$command -constraint ${constraint}"
					set command = "$command $script $rnb $tag_sys $counter $gendir $writedir $jf"
					
					echo "$command" > $jf
					
					$command
					
					set ij = `expr $ij + 1`
					
				endif
				
			endif
			
			#***********************************************************#
			
			set counter = `expr $counter + 1`
		end
		
	end
end

echo "submitted $ij jobs"
