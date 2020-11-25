#!/usr/bin/tcsh
#

set tag_sys = $1
set counter = $2
set jf      = $3
set outdir  = $4

set ebeam   = "11.6061"
set ngen    = "250000"


set loc_int = ${counter}
set pwd     = `pwd`


set ofile = "test.root"

if( ${counter} < "10" ) then
  set ofile = "${tag_sys}_00${counter}"
else 
  if( ${counter} < "100" ) then
    set ofile = "${tag_sys}_0${counter}"
  else 
    set ofile = "${tag_sys}_${counter}"
  endif
endif



cp ${HALLD_SIM_HOME}/src/programs/Simulation/gen_primex_compton/sd_compton/bases-init.dbf .
cp ${HALLD_SIM_HOME}/src/programs/Simulation/gen_primex_compton/sd_compton/primex_tagh.txt .
cp ${HALLD_SIM_HOME}/src/programs/Simulation/gen_primex_compton/sd_compton/primex_tagm.txt .



echo "sd_compton bases-init.dbf test ${ebeam} ${ngen} primex_${tag_sys}.txt ${counter} 1"
sd_compton bases-init.dbf test ${ebeam} ${ngen} primex_${tag_sys}.txt ${counter} 1



h2root test.hbook



rm bases-init.dbf
rm primex_tagh.txt
rm primex_tagm.txt


mv test.root ${outdir}/genRootFiles/${ofile}.root
mv test.dat ${outdir}/genCS/${ofile}.dat


rm test.bin
rm test.hbook
