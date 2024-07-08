#!/usr/bin/tcsh
#

set rnb         = $1
set tag_sys     = $2
set tag_counter = $3
set writedir    = $4
set jf          = $5

# setup environment:

source /group/halld/Software/build_scripts/gluex_env_jlab.csh /home/andrsmit/version.xml

# adjust beam energy according to run number:

set endpoint_energy = 11.6061
set endpoint_calib  = 11.6061

if ( "$rnb" == 61866 ) then
	set endpoint_energy = 11.1666
	set endpoint_calib  = 11.1689
endif
if ( "$rnb" == 61914 ) then
	set endpoint_energy = 11.1671
	set endpoint_calib  = 11.1689
endif
if ( "$rnb" == 61947 ) then
	set endpoint_energy = 11.1664
	set endpoint_calib  = 11.1689
endif

# Number of events to generate:
set ngen  = "500000"

# Output file name:
set ofile = "${tag_sys}_${tag_counter}"

# Generate events:

cp ${HALLD_SIM_HOME}/src/programs/Simulation/gen_primex_compton/sd_compton/bases-init.dbf .
cp /work/halld/home/andrsmit/primex_compton_analysis/photon_flux/phase1/primex_tagh.txt .
cp /work/halld/home/andrsmit/primex_compton_analysis/photon_flux/phase1/primex_tagm.txt .

echo "sd_compton bases-init.dbf test ${endpoint_energy} ${endpoint_calib} ${ngen} primex_${tag_sys}.txt ${tag_counter} 1"
sd_compton bases-init.dbf test ${endpoint_energy} ${endpoint_calib} ${ngen} primex_${tag_sys}.txt ${tag_counter} 1

# apparently root-6.24.04 wasn't built with h2root...

source /group/halld/Software/build_scripts/gluex_env_jlab.csh /home/andrsmit/version_h2root.xml

h2root test.hbook

rm bases-init.dbf
rm primex_tagh.txt
rm primex_tagm.txt

mv test.root ${writedir}/genRootFiles/${ofile}.root
mv test.dat ${writedir}/genCS/${ofile}.dat

rm test.bin
rm test.hbook

rm $jf
