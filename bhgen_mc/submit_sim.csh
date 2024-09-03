#!/usr/bin/tcsh
#
set outdir = /work/halld/home/andrsmit/primex_compton_analysis/bhgen_mc
set script = ${outdir}/simulate_pairs.csh
set ij     = 0

set submit_runs = go

#=========================#
#Job Resources:

set workflow   = bhgen_test
set account    = halld
set partition  = production
set ram        = 16GB
set cores      = 10
set time       = 14hr
set disk       = 18GB
set constraint = el9

#=========================#

set run_no = 61866

set run_number = $run_no
if( $run_number < 100000 ) then
	set run_number = 0$run_no
endif

set max_seeds = 15000

set seed = 1
while ( $seed <= $max_seeds )
	
	set seed_val = $seed
	if( $seed_val < 10 ) then
		set seed_val = 0000$seed
	else if ( $seed_val < 100 ) then
		set seed_val = 000$seed
	else if ( $seed_val < 1000 ) then
		set seed_val = 00$seed
	else if ( $seed_val < 10000 ) then
		set seed_val = 0$seed
	endif
	
	#***********************************************************#
	# submit jobs:
	
	set fname = pair_sim_${run_number}_${seed_val}
	set jf    = $outdir/jsub/$fname.jsub
	
	# skip jobs that have already been submitted:
	if ( -f $jf ) then
		set seed = `expr $seed + 1`
		continue
	endif
	
	set ofname = pair_sim_${seed_val}
	
	# skip jobs that have already been analyzed:
	if ( -f $outdir/recHddmFiles/Run${run_number}/$ofname.hddm ) then
		set seed = `expr $seed + 1`
		continue
	endif
	
	echo "$fname"
	
	if($submit_runs == "go") then
		
		set command = "swif2 add-job -workflow ${workflow}"
		set command = "$command -name ${fname}"
		set command = "$command -account ${account} -partition ${partition}"
		set command = "$command -cores ${cores} -ram ${ram} -time ${time}"
		set command = "$command -disk ${disk} -constraint ${constraint}"
		set command = "$command $script $run_no $seed $jf"
		
		echo "$command" > $jf
		
		$command
		
	endif
	
	set ij   = `expr $ij   + 1`
	set seed = `expr $seed + 1`
	
end # while ( $seed < $max_seeds )

echo "submitted $ij jobs"
