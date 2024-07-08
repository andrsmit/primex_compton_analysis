#!/usr/bin/tcsh
#

set run      = $1
set tag_sys  = $2
set counter  = $3
set gendir   = $4
set writedir = $5
set jf       = $6

#------------------------------------------------------------------------------------------#
#Source setup script:

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.csh
gxenv /home/andrsmit/version.xml

#------------------------------------------------------------------------------------------#
# Set appropriate geometry:

#setenv JANA_GEOMETRY_URL ccdb:///GEOMETRY/main_HDDS.xml
setenv JANA_GEOMETRY_URL xmlfile:///${HDDS_HOME}/main_HDDS.xml
setenv JANA_CALIB_URL mysql://ccdb_user@hallddb-farm.jlab.org/ccdb
setenv JANA_CALIB_CONTEXT "variation=mc"

cp ${gendir}/${tag_sys}_${counter}.hddm compton_gen.hddm
cp ${writedir}/control.in .

#------------------------------------------------------------------------------------------#
# Simulate:

set bin1 = ${HDGEANT4_HOME}/bin/Linux-g++/hdgeant4

echo "$bin1 -t10"
$bin1 -t10

#------------------------------------------------------------------------------------------#
# Save smeared file to output directory:

echo "mv compton_rec_smeared.hddm $writedir/recHddmFiles_noBkgd/${tag_sys}_${counter}.hddm"
mv compton_rec_smeared.hddm $writedir/recHddmFiles_noBkgd/${tag_sys}_${counter}.hddm

#------------------------------------------------------------------------------------------#
# Run mcsmear again and mix in random trigger background

#cp /work/halld/home/andrsmit/primex_compton_analysis/compton_mc/phaseI/random_hddm/random_61866-61895.hddm random.hddm
echo "cp /work/osgpool/halld/random_triggers/recon-2019_01-ver01/run${run}_random.hddm random.hddm"
cp /work/osgpool/halld/random_triggers/recon-2019_01-ver01/run${run}_random.hddm random.hddm

set bin2 = $HALLD_SIM_HOME/Linux_Alma9-x86_64-gcc11.4.1/bin/mcsmear

echo "$bin2 compton_rec.hddm random.hddm:1"
$bin2 compton_rec.hddm random.hddm:1

#------------------------------------------------------------------------------------------#
# Save smeared files to output directory and clean up:

echo "mv compton_rec.hddm $writedir/recHddmFiles_unsmeared/${tag_sys}_${counter}.hddm"
mv compton_rec.hddm $writedir/recHddmFiles_unsmeared/${tag_sys}_${counter}.hddm

echo "mv compton_rec_smeared.hddm $writedir/recHddmFiles/${tag_sys}_${counter}.hddm"
mv compton_rec_smeared.hddm $writedir/recHddmFiles/${tag_sys}_${counter}.hddm

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
rm $jf
