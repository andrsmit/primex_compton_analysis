#!/usr/bin/tcsh
#

set run_number   = $1
set beam_current = $2

./ana_compton -r$run_number -c$beam_current -s0 -b1   -e50  > ${run_number}_${beam_current}_tagh_001_050.log &
./ana_compton -r$run_number -c$beam_current -s0 -b51  -e100 > ${run_number}_${beam_current}_tagh_051_100.log &
./ana_compton -r$run_number -c$beam_current -s0 -b101 -e190 > ${run_number}_${beam_current}_tagh_101_190.log &
./ana_compton -r$run_number -c$beam_current -s0 -b191 -e221 > ${run_number}_${beam_current}_tagh_191_221.log &
./ana_compton -r$run_number -c$beam_current -s1 -b1   -e51  > ${run_number}_${beam_current}_tagm_001_051.log &
./ana_compton -r$run_number -c$beam_current -s1 -b52  -e102 > ${run_number}_${beam_current}_tagm_052_102.log &
