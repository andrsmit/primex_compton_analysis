#!/usr/bin/tcsh
#

set run      = $1
set tag_sys  = $2
set counter  = $3
set gendir   = $4
set writedir = $5
set jf       = $6

#Source setup script:

source /group/halld/Software/build_scripts/gluex_env_jlab.csh $writedir/version.xml

#setenv JANA_GEOMETRY_URL ccdb:///GEOMETRY/main_HDDS.xml
setenv JANA_GEOMETRY_URL xmlfile:///${HDDS_HOME}/main_HDDS.xml
setenv JANA_CALIB_URL mysql://ccdb_user@hallddb-farm.jlab.org/ccdb
setenv JANA_CALIB_CONTEXT "variation=mc"

cp ${gendir}/${tag_sys}_${counter}.hddm compton_gen.hddm
cp ${writedir}/control.in .

set bin1 = ${HDGEANT4_HOME}/bin/Linux-g++/hdgeant4

echo "$bin1 -t10"
$bin1 -t10

mv compton_rec_smeared.hddm $writedir/recHddmFiles_noBkgd/${tag_sys}_${counter}.hddm

cp /work/halld/home/andrsmit/primex_compton_analysis/compton_mc/phaseI/random_hddm/random_61866-61895.hddm random.hddm

set bin2 = $HALLD_SIM_HOME/Linux_CentOS7.7-x86_64-gcc4.8.5/bin/mcsmear

echo "$bin2 compton_rec.hddm random.hddm:1"
$bin2 compton_rec.hddm random.hddm:1

mv compton_rec.hddm $writedir/recHddmFiles_unsmeared/${tag_sys}_${counter}.hddm
mv compton_rec_smeared.hddm $writedir/recHddmFiles/${tag_sys}_${counter}.hddm

rm smear.root
rm compton_gen.hddm
rm control.in
rm BHgen*.astate
rm random.hddm
rm $jf
