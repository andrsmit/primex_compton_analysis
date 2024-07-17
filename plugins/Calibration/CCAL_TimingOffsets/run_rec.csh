#!/usr/bin/tcsh
#

# command line inputs:

set outdir    = $1
set runnumber = $2
set ext       = $3
set mss_fname = $4
set loc_fname = $5
set jf_fname  = $6

# source GlueX setup scripts:

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.csh
gxenv /home/andrsmit/version.xml

setenv CCDB_CONNECTION mysql://ccdb_user@hallddb-farm.jlab.org/ccdb
setenv JANA_CALIB_URL mysql://ccdb_user@hallddb-farm.jlab.org/ccdb

# use default ccdb variation:
#setenv JANA_CALIB_CONTEXT variation=ccal_energy_calib

setenv HALLD_MY /home/andrsmit/halld_my

# remove symlink created by swif and replace with local file:

if ( -f /cache$mss_fname ) then
	echo "copying file from cache..."
	rm $loc_fname
	cp /cache$mss_fname $loc_fname
endif

# set up hd_root:

set bin = ${HALLD_RECON_HOME}/Linux_Alma9-x86_64-gcc11.4.1/bin/hd_root

echo "$bin -PPLUGINS=CCAL_TimingOffsets -PEVIO:PARSE_F125=0 -PNTHREADS=4 $loc_fname"
$bin -PPLUGINS=CCAL_TimingOffsets -PEVIO:PARSE_F125=0 -PNTHREADS=4 $loc_fname

# move output root file to appropriate directory

mv hd_root.root $outdir/${runnumber}_${ext}.root

rm $jf_fname
rm -f $loc_fname
