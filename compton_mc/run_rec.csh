#!/usr/bin/tcsh
#

set tag_sys = $1
set counter = $2
set jfile   = $3
set outDir  = $4


set bin   = /work/halld/home/andrsmit/halld_software/halld_recon/Linux_CentOS7.7-x86_64-gcc4.8.5/bin/hd_root

setenv BMS_OSNAME Linux_CentOS7.7-x86_64-gcc4.8.5
setenv HALLD_RECON_HOME /work/halld/home/andrsmit/halld_software/halld_recon
setenv CCDB_CONNECTION mysql://ccdb_user@hallddb.jlab.org/ccdb
setenv JANA_CALIB_URL  mysql://ccdb_user@hallddb.jlab.org/ccdb
setenv JANA_CALIB_CONTEXT "variation=mc"

setenv JANA_GEOMETRY_URL xmlfile://${HDDS_HOME}/main_HDDS.xml

cp $outDir/recHddmFiles_v1/${tag_sys}_${counter}.hddm output.hddm


$bin -PPLUGINS=compton_simulation -PNTHREADS=5 -PCCAL:DO_NONLINEAR_CORRECTION=0 output.hddm

mv hd_root.root ${outDir}/recRootFiles_v1_sys/${tag_sys}_${counter}.root



# clean local directory:

rm output.hddm
rm $jfile
