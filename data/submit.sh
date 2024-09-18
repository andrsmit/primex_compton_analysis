#!/usr/bin/sh
#

submit_runs=go

# base directory where output will be stored:
outdir=/work/halld/home/andrsmit/primex_compton_analysis/data

####################################################
# Modify these parameters:

# set which run period to analyze:
phase=1

# get list of runs we want to analyze:
run_list=${outdir}/run_list.dat

####################################################

# script that will execute the commands for each job:
script=${outdir}/run_rec.csh

# JANA configuration file:
cfg_file=${outdir}/jana.config

#-----------------------------------#

run_period="2019-01"
if [[ $phase == "1" ]]; then
	run_period="2019-01"
elif [[ $phase == "2" ]]; then
	run_period="2021-08"
elif [[ $phiase == "3" ]]; then
	run_period="2022-08"
else
	echo "Non-primex run period specified. Example syntax: ./submit.csh [1-3]"
	exit 1
fi

#-----------------------------------#

# directoy where rawdata evio files are stored:
dir_mss=/mss/halld/RunPeriod-${run_period}/rawdata

# this variable will count how many jobs are submitted:
ij=0

#=========================#
#Job Resources:

workflow=compton_tree_production_phase$phase
account=halld
partition=production
ram=12GB
cores=4
run_time=120min
disk=23GB
constraint=el9

#=========================#

# make sure config file exists:
if [ ! -f $cfg_file ]; then
	echo "JANA config file does not exist."
	exit 1
fi

while read run; do
	
	# Convert run number into 6 digit variable:
	
	run_number=$run
	if [[ $run_number -lt "100000" ]]; then
		run_number="0$run"
	fi
	echo "run ${run_number}"
	
	
	# set up output directory where rootTrees will be written:
	dir_tree=${outdir}/rootTrees/phase${phase}/${run_number}
	mkdir -p $dir_tree
	
	# set up output directory where hd_root files will be written:
	dir_rfile=${outdir}/rootFiles/phase${phase}/${run_number}
	mkdir -p $dir_rfile
	
	# set up output directory where jsub will be written:
	dir_jsub=${outdir}/jsub/phase${phase}
	mkdir -p $dir_jsub
	
	# Loop over file extensions:
	for file in ${dir_mss}/Run${run_number}/hd_rawdata_${run_number}_*.evio; do
		
		# get extension number from filename:
		cut_c1=$((${#file}-7))
		cut_c2=$((${#file}-5))
		ext=`echo "$file" | cut -c ${cut_c1}-${cut_c2}`
		
		# set up output root tree file name:
		tree_file=${dir_tree}/${run_number}_${ext}.root
		if [ -f $tree_file ]; then
			continue
		fi
		
		# set up output hd_root file name:
		root_file=${dir_rfile}/${run_number}_${ext}.root
		if [ -f $root_file ]; then
			continue
		fi
		
		# set up job submission file name:
		jsub_file=${dir_jsub}/${run_number}_${ext}.jsub
		if [ -f $jsub_file ]; then
			continue
		fi
		
		# Get the file size of the input file so we know how much disk space to request:
		sum=0
		while IFS= read -r line; do
		if [[ $line =~ size=([0-9]+) ]]; then
			size=${BASH_REMATCH[1]}
	    	# Add the number to the sum
	    	sum=$((sum + size))
		fi
    	done < "$file"
		sum_gb=$(($sum / 1000000000))
		#echo "$sum_gb"
		disk="$(($sum_gb+2))GB"
		
		# set job name:
		job_name="compton_tree_${run_number}_${ext}"
		
		command="swif2 add-job -workflow ${workflow}"
		command="$command -name ${job_name} -account ${account} -partition ${partition}"
		command="$command -cores ${cores} -ram ${ram} -time ${run_time} -disk ${disk}"
		command="$command -constraint ${constraint}"
		command="$command -input `basename $file` mss:${file}"
		command="$command $script `basename $file` $tree_file $root_file $cfg_file $jsub_file"
		
		if [ $submit_runs == "go" ]; then
			echo "$command" > $jsub_file
			$command
		fi
		ij=$(($ij+1))
	done #loop over file extensions
	
done < $run_list

echo "submitted $ij jobs"
if [ $submit_runs == "go" ]; then
	swif2 run $workflow
fi
