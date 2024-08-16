#!/usr/bin/tcsh
#
# 28 runs split over 7 processes (left out run 61879):

./ana_compton -b61866 -e61874 > 61866_61874.log &
./ana_compton -b61875 -e61878 > 61875_61878.log &
./ana_compton -b61880 -e61883 > 61880_61883.log &
./ana_compton -b61884 -e61888 > 61884_61888.log &
./ana_compton -b61889 -e61892 > 61889_61892.log &
./ana_compton -b61893 -e61905 > 61893_61905.log &
./ana_compton -b61906 -e61910 > 61906_61910.log &

