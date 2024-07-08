#!/usr/bin/tcsh
#
set outdir = /work/halld/home/andrsmit/primex_compton_analysis/compton_mc/phase1/genDir
set script = ${outdir}/gen_compton.csh
set ij     = 0

set submit_runs = go

#=========================#
#Job Resources:

set workflow  = compton_sim
set account   = halld
set partition = production
set ram       = 4GB
set cores     = 1
set time      = 10min
set disk      = 500MB

#=========================#

set run_list = {"61321","61866"}
set tag_list = {"tagh","tagm"}

foreach rnb ($run_list)
	
	set run_number = $rnb
	if ( $rnb < 100000 ) then:
		set run_number = "0$rnb"
	endif
	
	foreach tag_sys ($tag_list)
		
		set max_counter = 274
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
			
			set loc_jobname = compton_gen_rnb-${run_number}_${tag_sys}_${tag_counter}
			
			set writedir = ${outdir}/Run$run_number
			set jf       = ${outdir}/jsub/${loc_jobname}.jsub
			
			# skip jobs that have already been submitted:
			if ( -f $jf ) then
				set counter = `expr $counter + 1`
				continue
			endif
			
			# skip counters which have already been processed:
			if ( -f ${writedir}/genRootFiles/${tag_sys}_${tag_counter}.root ) then
				set counter = `expr $counter + 1`
				continue
			endif
			
			echo "${tag_sys}_${tag_counter}"
			
			if($submit_runs == "go") then
				
				set command = "swif2 add-job -workflow ${workflow}"
				set command = "$command -name ${loc_jobname}"
				set command = "$command -account ${account} -partition ${partition}"
				set command = "$command -cores ${cores} -ram ${ram} -time ${time} -disk ${disk}"
				set command = "$command $script $rnb $tag_sys $tag_counter $writedir $jf"
				
				echo "$command" > $jf
				
				$command
				
				set ij = `expr $ij + 1`
				
			endif
			
			#***********************************************************#
			
			set counter = `expr $counter + 1`
		end
		
	end
end

echo "submitted $ij jobs"
