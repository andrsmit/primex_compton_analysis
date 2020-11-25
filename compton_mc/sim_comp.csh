#!/usr/bin/tcsh
#

#nominal beam spot: 
# x, y, z, var_xx, var_xy, var_yy, dxdz, dydz:
# 0.1755, -0.05027, 64.939, 0.0625, 0.0, 0.0729, 0.0, 0.0

set tag_sys = $1
set counter = $2
set jf      = $3
set odir    = $4


set pwd  = `pwd`


#Source setup script:

source $BUILD_SCRIPTS/gluex_env_jlab.csh ${odir}/version.xml

setenv JANA_GEOMETRY_URL xmlfile:///${odir}/hdds/main_HDDS.xml
setenv JANA_CALIB_URL mysql://ccdb_user@hallddb.jlab.org/ccdb
setenv JANA_CALIB_CONTEXT "variation=mc"




set genDir = /work/halld/home/andrsmit/compton_analysis/tag_sim/genHddmFiles

cp $genDir/${tag_sys}_${counter}.hddm compton_gen.hddm
cp /work/halld/home/andrsmit/compton_analysis/tag_sim/control.in .

set bin1  = ${HDGEANT4_HOME}/bin/Linux-g++/hdgeant4




echo "${bin1} -t5"
$bin1 -t5


mv compton_rec_smeared.hddm $odir/recHddmFiles_v1/${tag_sys}_${counter}.hddm
mv compton_rec.hddm $odir/recHddmFiles_unsmeared_v1/${tag_sys}_${counter}.hddm
rm compton_rec.hddm


rm smear.root
rm compton_gen.hddm
rm control.in
rm ${jf}
