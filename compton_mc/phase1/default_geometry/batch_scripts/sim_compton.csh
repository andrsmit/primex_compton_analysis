#!/usr/bin/tcsh
#

set run      = $1
set tag_sys  = $2
set counter  = $3
set gendir   = $4
set writedir = $5
set jf       = $6

set n_events  = 2000000 # just needs to be larger than number of events from input hddm file
set n_threads = 10      # This should be adjusted according to the number of cores for the submitted job

#------------------------------------------------------------------------------------------#
# Adjust geometry depending on run number:

# Set Defaults:

set endpoint_energy     = 11.6061 # I don't think this is actually used in the simulation...
set beam_spot_x         =   0.027
set beam_spot_y         =  -0.128
set beam_spot_z         =  64.935
set target_length       =  1.7755
set target_string       = "Be"
set beam_current_vals   = {200}
set random_trigger_path = "/work/halld/home/andrsmit/primex_compton_analysis/data/random_triggers/phase1"

# Adjust based on run number:

if ( $run > 60000 && $run < 69999 ) then
	
	# PrimEx-eta Phase 1
	
	#endpoint energy:
	if ( $run > 61433 ) then
		set endpoint_energy = 11.1689
	endif
	
	#target z:
	if ( $run > 61354 ) then
		set beam_spot_z   = 65.0
		set target_length = 29.5
		set target_string = "He"
		set beam_current_vals = {200,100,50}
	endif
	
	#beam spot:
	if ( $run < 61483 ) then
		set beam_spot_x =  0.027
		set beam_spot_y = -0.128
	else if ( $run < 61774 ) then
		set beam_spot_x =  0.001
		set beam_spot_y = -0.077
	else
		set beam_spot_x =  0.038
		set beam_spot_y = -0.095
	endif
	
	#random trigger background:
	#if ( $run < 61355 ) then
	#	set random_trigger_file = "Be_200nA_FIELDOFF_random.hddm"
	#else if ( $run < 61914 ) then
	#	set random_trigger_file = "He_200nA_FIELDOFF_random.hddm"
	#else if ( $run < 61949 ) then
	#	set random_trigger_file = "He_050nA_FIELDOFF_random.hddm"
	#else 
	#	set random_trigger_file = "He_200nA_FIELDOFF_random.hddm"
	#endif
	
	#set random_trigger_file = ${random_trigger_path}/${random_trigger_file}
	
endif

set beam_current_string = ""
foreach beam_current ($beam_current_vals)
	set beam_current_string = "$beam_current_string$beam_current, "
end

echo ""
echo "Simulating Compton scattering for Run $run"
echo "  ${target_string} Target"
echo "  Beam spot: (${beam_spot_x}, ${beam_spot_y}, ${beam_spot_z})"
echo "  Target length: $target_length"
echo "  Random trigger hddm path: $random_trigger_path"
echo "  Beam currents: $beam_current_string"
echo ""

#------------------------------------------------------------------------------------------#
# Setup the counter number string:

set tag_counter = $counter
if ( $counter < 100 ) then
	if ( $counter < 10 ) then
		set tag_counter = 00$counter
	else
		set tag_counter = 0$counter
	endif
endif

#------------------------------------------------------------------------------------------#
# Set random seed according to tagger counter number:

set rndm_seed = $counter
if ( $counter > 178 ) then
	set rndm_seed = `expr $rndm_seed - 51`
endif

#------------------------------------------------------------------------------------------#
# Source setup script:

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.csh
gxenv ${writedir}/version.xml

#------------------------------------------------------------------------------------------#
# Set appropriate geometry:

#setenv JANA_GEOMETRY_URL ccdb:///GEOMETRY/main_HDDS.xml
setenv JANA_GEOMETRY_URL xmlfile:///${HDDS_HOME}/main_HDDS.xml
setenv JANA_CALIB_URL mysql://ccdb_user@hallddb-farm.jlab.org/ccdb
setenv JANA_CALIB_CONTEXT "variation=mc"

cp ${gendir}/${tag_sys}_${tag_counter}.hddm compton_gen.hddm

#------------------------------------------------------------------------------------------#
# Setup the control.in file:

set ctrl_file = control.in
echo "INFILE 'compton_gen.hddm'" > $ctrl_file
echo "TRIG ${n_events}" >> $ctrl_file
echo "RUNG ${run}" >> $ctrl_file
echo "BEAM ${endpoint_energy} 0.0 3.0 76.00 0.005 10.e-9 20.e-6 1e-3 -0.0 +0.0" >> $ctrl_file
echo "VERTEX 'beam_spot(${beam_spot_x},${beam_spot_y},${beam_spot_z},0.0625,0.0,0.0625,0.0,0.0) * ${target_length}'" >> $ctrl_file
echo "OUTFILE 'compton_rec.hddm'" >> $ctrl_file
echo "POSTSMEAR 1" >> $ctrl_file
echo "DELETEUNSMEARED 0" >> $ctrl_file
echo "MCSMEAROPTS '-PNTHREADS=${n_threads}'" >> $ctrl_file
echo "RNDM ${rndm_seed}" >> $ctrl_file
echo "TOFMAX 1e-5" >> $ctrl_file
echo "CKOV 1" >> $ctrl_file
echo "LABS 1" >> $ctrl_file
echo "END" >> $ctrl_file

#------------------------------------------------------------------------------------------#
# Simulate:

set bin1 = ${HDGEANT4_HOME}/bin/Linux-g++/hdgeant4

echo "$bin1 -t${n_threads}"
$bin1 -t${n_threads}

# before continuing, copy the unsmeared hddm file to the output directory:

echo "cp compton_rec.hddm $writedir/recHddmFiles_unsmeared/${tag_sys}_${tag_counter}.hddm"
cp compton_rec.hddm $writedir/recHddmFiles_unsmeared/${tag_sys}_${tag_counter}.hddm

#------------------------------------------------------------------------------------------#
# Produce root tree for file with noBkgd:

gxenv /home/andrsmit/version.xml

set bin2 = ${HALLD_SIM_HOME}/Linux_Alma9-x86_64-gcc11.4.1/bin/mcsmear
set bin3 = ${HALLD_RECON_HOME}/Linux_Alma9-x86_64-gcc11.4.1/bin/hd_root

set plugin_name   = "compton_tree"
set plugin_opts   = "-PNTHREADS=${n_threads}"

echo "$bin3 -PPLUGINS=${plugin_name} ${plugin_opts} compton_rec_smeared.hddm"
$bin3 -PPLUGINS=${plugin_name} ${plugin_opts} compton_rec_smeared.hddm

echo "mv primex_compton.root $writedir/recRootTrees_noBkgd/${tag_sys}_${tag_counter}.root"
mv primex_compton.root $writedir/recRootTrees_noBkgd/${tag_sys}_${tag_counter}.root

echo "rm hd_root.root"
rm hd_root.root

echo "mv compton_rec_smeared.hddm $writedir/recHddmFiles_noBkgd/${tag_sys}_${tag_counter}.hddm"
mv compton_rec_smeared.hddm $writedir/recHddmFiles_noBkgd/${tag_sys}_${tag_counter}.hddm

#------------------------------------------------------------------------------------------#
# Loop over all beam currents and re-smear simulation, adding in random triggers from selected hddm file:

foreach beam_current ($beam_current_vals)
	
	set loc_beam_current = $beam_current
	if ( $beam_current < 100 ) then
		set loc_beam_current = 0$beam_current
	endif
	
	set random_trigger_file = ${random_trigger_path}/${target_string}_${loc_beam_current}nA_FIELDOFF_random.hddm
	
	if ( -f $random_trigger_file ) then
		
		#------------------------------------------------------------------------------------------#
		# Smear with selected random trigger file:
		
		echo "cp ${random_trigger_file} random.hddm"
		cp ${random_trigger_file} random.hddm
		
		echo "$bin2 -PNTHREADS=${n_threads} compton_rec.hddm random.hddm:1"
		$bin2 -PNTHREADS=${n_threads} compton_rec.hddm random.hddm:1
		
		#------------------------------------------------------------------------------------------#
		# Next, produce ROOT tree:
		
		echo "$bin3 -PPLUGINS=${plugin_name} ${plugin_opts} compton_rec_smeared.hddm"
		$bin3 -PPLUGINS=${plugin_name} ${plugin_opts} compton_rec_smeared.hddm
		
		echo "mv primex_compton.root $writedir/recRootTrees_${loc_beam_current}nA/${tag_sys}_${tag_counter}.root"
		mv primex_compton.root $writedir/recRootTrees_${loc_beam_current}nA/${tag_sys}_${tag_counter}.root
		
		echo "rm hd_root.root"
		rm hd_root.root
		
		echo "mv compton_rec_smeared.hddm $writedir/recHddmFiles_${loc_beam_current}nA/${tag_sys}_${tag_counter}.hddm"
		mv compton_rec_smeared.hddm $writedir/recHddmFiles_${loc_beam_current}nA/${tag_sys}_${tag_counter}.hddm
		
	endif
end

#------------------------------------------------------------------------------------------#
# Clean up local directory:

echo "rm compton_rec.hddm $writedir/recHddmFiles_unsmeared/${tag_sys}_${tag_counter}.hddm"
rm compton_rec.hddm $writedir/recHddmFiles_unsmeared/${tag_sys}_${tag_counter}.hddm

echo "rm smear.root"
rm smear.root

echo "rm compton_gen.hddm"
rm compton_gen.hddm

echo "rm control.in"
rm control.in

echo "rm BHgen*.astate"
rm BHgen*.astate

echo "rm random.hddm"
rm random.hddm

echo "rm $jf"
#rm $jf
