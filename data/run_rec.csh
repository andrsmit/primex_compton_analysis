#!/usr/bin/tcsh
#

# command line inputs:

set input_file = $1
set  tree_file = $2
set  root_file = $3
set   cfg_file = $4
set  jsub_file = $5

# source GlueX setup scripts:

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.csh
gxenv /home/andrsmit/version.xml

setenv CCDB_CONNECTION mysql://ccdb_user@hallddb-farm.jlab.org/ccdb
setenv JANA_CALIB_URL mysql://ccdb_user@hallddb-farm.jlab.org/ccdb

setenv HALLD_MY /home/andrsmit/halld_my

echo "pwd"
pwd

echo "ls -a"
ls -a

# copy JANA config file:

echo "cp ${cfg_file} jana.config"
cp ${cfg_file} jana.config

# set up hd_root:

set bin = ${HALLD_RECON_HOME}/Linux_Alma9-x86_64-gcc11.4.1/bin/hd_root

echo "$bin --config=jana.config -PNTHREADS=4 $input_file"
$bin --config=jana.config -PNTHREADS=4 $input_file

# move output root tree to appropriate directory:
if ( -f primex_compton.root ) then
	echo "mv primex_compton.root ${tree_file}"
	mv primex_compton.root ${tree_file}
endif

# move output root file to appropriate directory:
if ( -f hd_root.root ) then
	echo "mv hd_root.root ${root_file}"
	mv hd_root.root ${root_file}
endif

# delete local copies of input files and jsub file:

echo "rm -f ${input_file}"
rm -f ${input_file}

echo "rm -f ${jsub_file}"
rm -f ${jsub_file}
