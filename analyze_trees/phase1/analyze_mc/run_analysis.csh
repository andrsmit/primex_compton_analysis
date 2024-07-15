#!/usr/bin/tcsh
#

./ana_compton -s0 -b1 -e50 > tagh_001_050.log &
./ana_compton -s0 -b51 -e100 > tagh_051_100.log &
./ana_compton -s0 -b101 -e190 > tagh_101_190.log &
./ana_compton -s0 -b191 -e221 > tagh_191_221.log &
./ana_compton -s1 -b1 -e51 > tagm_001_051.log &
./ana_compton -s1 -b52 -e102 > tagm_052_102.log &

#./ana_compton -s1 -b1 -e25 > tagm_001_025.log &
#./ana_compton -s1 -b26 -e50 > tagm_026_050.log &
#./ana_compton -s1 -b51 -e75 > tagm_051_075.log &
#./ana_compton -s1 -b76 -e102 > tagm_076_102.log &
