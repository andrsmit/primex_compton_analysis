#!/usr/bin/tcsh
#

set run_no = $1
set seed   = $2
set jfile  = $3

# number of events to simulate in each iteration:
set n_events = 7500000

# convert run number into a 6-digit number:
set run_number = $run_no
if( $run_number < 100000 ) then
	set run_number = 0$run_no
endif

# setup 'pathname' variable for convenience:
set pathname = /work/halld/home/andrsmit/primex_compton_analysis/bhgen_mc

setenv JANA_CALIB_URL mysql://ccdb_user@hallddb-farm.jlab.org/ccdb
setenv JANA_CALIB_CONTEXT variation=mc

#==================================================================#
# Get target length from run number:

set target_length = 29.5
if( $run_number > 61320 && $run_number < 61355 ) then
	set target_length = 1.775
else if( $run_number > 81261 && $run_number < 81396 ) then
	set target_length = 1.775
else if( $run_number > 110445 && $run_number < 110622 ) then
	set target_length = 1.775
endif
echo "target_length = $target_length"

#==================================================================#
# Get random number seeds:

set rnd_max = 215

# seed: 1 --> $rnd1 

set rnd1 = `expr $seed / $rnd_max`
set num2 = `expr $rnd1 \* $rnd_max`
set rnd2 = `expr $seed - $num2`

set rnd1 = `expr $rnd1 + 1`
set rnd2 = `expr $rnd2 + 1`

#echo "RNDM $rnd1 $rnd2"

#==================================================================#
# convert seed to 5 digit number:

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

#==================================================================#

#
# Procedure:
#  1. Generate list of photons from python script simple_beam_generator.py
#     Syntax: 'python simple_beam_generator.py <n_evts> <run_number>
#     <run_number> should either be 61321 (Be) or 61866 (He) because it will look for the corresponding flux ROOT file
#     If it can't find that ROOT file <run_number>_flux.root, then it will generate photons uniformly in energy
#
#  2. Run hdgeant4 with hddm file created in above step as input using GENBEAM 'BHgen'
#  2.1 Delete initial photon beam input file
#
#  3. Run accept-reject script to equalize weights of e+e- pairs and write them to a new hddm file:
#     Syntax: 'python accept_pairs.py <input_hddm_filename>'
#     Note: When writing to a new file, we need to delete the vertex containing the thrown photon
#
#  4. Clean directory:
#     Delete full hddm file from hdgeant4 (expected to be about 17.5GB for n_evts=30M
#
#  5. Repeate steps 1-4 5 times, and add all files together
#
#  6. Use hddm file from step 4 as the input to hdgeant4 (~1M events when n_evs=30M) and simulate through entire setup
#  6.1 smear output hddm file
#  6.2 Delete unsmeared file
#  6.3 Delete smear.root BHgen_thread* BHgen_stats.astate
#  
#  7. Run hd_root over smeared hddm file from previous step
#
#  8. Move hd_root.root and smeared hddm file to output directory
#

set n_iterations = 5

set loc_seed_base = `expr $seed \* $n_iterations`

set add_cmd = "python ${pathname}/combine_hddm_files.py"
set clean_cmd = "rm "

set ctrl = control.in
set hdg4 = ${HDGEANT4_HOME}/bin/Linux-g++/hdgeant4

set ngen_val = 0
while ( $ngen_val < $n_iterations )
	
	set loc_seed = `expr $loc_seed_base + $ngen_val`
	
	set rnd3 = `expr $loc_seed / $rnd_max`
	set num3 = `expr $rnd3 \* $rnd_max`
	set rnd4 = `expr $loc_seed - $num3`
	
	set rnd3 = `expr $rnd3 + 1`
	set rnd4 = `expr $rnd4 + 1`
	
	echo "RNDM $rnd3 $rnd4"
	
	#==================================================================#
	# 1. Generate list of photons from python script:
	
	echo "python ${pathname}/simple_beam_generator.py $n_events $run_no"
	python ${pathname}/simple_beam_generator.py $n_events $run_no
	
	#==================================================================#
	# 2. Run hdgeant4 with GENBEAM 'BHgen':
	
	set pair_gen_fname = pair_gen_$ngen_val.hddm
	
	# First, write the control.in file:
	
	echo "INFILE photonbeam.hddm" > $ctrl
	echo "TRIG ${n_events}" >> $ctrl
	echo "RUNG ${run_no}" >> $ctrl
	echo "VERTEX 'beam_spot(ccdb) * ${target_length}'" >> $ctrl
	echo "BEAM 11.1666 11.1666 3.0 76.0 0.005 4.e-9 50.e-6 0.5e-3" >> $ctrl
	echo "GENBEAM 'BHgen'" >> $ctrl
	echo "OUTFILE ${pair_gen_fname}" >> $ctrl
	echo "RNDM $rnd3 $rnd4" >> $ctrl
	echo "END" >> $ctrl
	
	# run hdgeant4:
	
	echo "$hdg4 -t10"
	$hdg4 -t10
	
	# Clean up:
	
	echo "rm BHgen_thread*.astate"
	rm BHgen_thread*.astate
	
	echo "rm BHgen_stats.astate"
	rm BHgen_stats.astate
	
	echo "rm control.in"
	rm control.in
	
	echo "rm photonbeam.hddm"
	rm photonbeam.hddm
	
	#==================================================================#
	# 3. Run script to select e+e- pairs produced inside the target
	
	echo "python $pathname/accept_pairs.py $pair_gen_fname"
	python $pathname/accept_pairs.py $pair_gen_fname
	
	#==================================================================#
	# 4. Delete hddm file from first pass of hdgeant4:
	
	echo "rm $pair_gen_fname"
	rm $pair_gen_fname
	
	echo "mv BHgen_accepted.hddm $pair_gen_fname"
	mv BHgen_accepted.hddm $pair_gen_fname
	
	set add_cmd = "${add_cmd} ${pair_gen_fname}"
	set clean_cmd = "${clean_cmd} ${pair_gen_fname}"
	set ngen_val = `expr $ngen_val + 1`
end

#==================================================================#
# 5. Combine all files from above together: pair_gen_x.hddm (x: 0-n_iterations)

echo "$add_cmd"
$add_cmd

echo "$clean_cmd"
$clean_cmd

echo "mv combined_hddm.hddm BHgen_${seed_val}.hddm"
mv combined_hddm.hddm BHgen_${seed_val}.hddm

#==================================================================#
# 6. Simulate pairs through target and experimental setup:

set pair_sim_fname = pair_sim_${seed_val}.hddm

# first, write a new control.in file:

echo "INFILE BHgen_${seed_val}.hddm" > $ctrl
echo "TRIG ${n_events}" >> $ctrl
echo "RUNG ${run_no}" >> $ctrl
echo "VERTEX 'beam_spot(ccdb) * ${target_length}'" >> $ctrl
echo "BEAM 11.1666 11.1666 3.0 76.0 0.005 4.e-9 50.e-6 0.5e-3" >> $ctrl
echo "OUTFILE ${pair_sim_fname}" >> $ctrl
echo "POSTSMEAR 1" >> $ctrl
echo "DELETEUNSMEARED 0" >> $ctrl
echo "MCSMEAROPTS '-PNTHREADS=10'" >> $ctrl
echo "RNDM $rnd1 $rnd2" >> $ctrl
echo "TOFMAX 1e-5" >> $ctrl
echo "CKOV 1" >> $ctrl
echo "LABS 1" >> $ctrl
echo "END" >> $ctrl

# run hdgeant4:

echo "$hdg4 -t10"
$hdg4 -t10

# Clean up:

echo "rm BHgen_${seed_val}.hddm"
rm BHgen_${seed_val}.hddm

echo "rm BHgen_thread*.astate"
rm BHgen_thread*.astate

echo "rm BHgen_stats.astate"
rm BHgen_stats.astate

echo "rm control.in"
rm control.in

echo "rm smear.root"
rm smear.root

# Delete unsmeared file

echo "rm $pair_sim_fname"
rm $pair_sim_fname

echo "mv pair_sim_${seed_val}_smeared.hddm $pair_sim_fname"
mv pair_sim_${seed_val}_smeared.hddm $pair_sim_fname

#==================================================================#
# 7. Run reconstruction:

set hdrt = ${HALLD_RECON_HOME}/Linux_Alma9-x86_64-gcc11.4.1/bin/hd_root

echo "$hdrt -PPLUGINS=compton_tree -PNTHREADS=10 ${pair_sim_fname}"
$hdrt -PPLUGINS=compton_tree -PNTHREADS=10 ${pair_sim_fname}

echo "mv primex_compton.root $pathname/recRootTrees/Run${run_number}/trees/pair_rec_${seed_val}.root"
mv primex_compton.root $pathname/recRootTrees/Run${run_number}/trees/pair_rec_${seed_val}.root

echo "mv hd_root.root $pathname/recRootTrees/Run${run_number}/hists/pair_rec_${seed_val}.root"
mv hd_root.root $pathname/recRootTrees/Run${run_number}/hists/pair_rec_${seed_val}.root

echo "mv ${pair_sim_fname} $pathname/recHddmFiles/Run${run_number}/${pair_sim_fname}"
mv ${pair_sim_fname} $pathname/recHddmFiles/Run${run_number}/${pair_sim_fname}

#==================================================================#
# Clean up further:

echo "rm $jfile"
#rm $jfile
