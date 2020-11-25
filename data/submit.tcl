#!/usr/bin/tclsh
#
set dir_mss /mss/halld/RunPeriod-2019-01/rawdata
set outdir /work/halld/home/andrsmit/primex_compton_analysis/data
set script ${outdir}/runbatch.tcl
set ext 0
set icounter 0
set ij 0


if 0 {
set runlist "61321 61322 61323 61325 61327 61329 61330 61331 61332 61333 \
	     61334 61335 61336 61337 61340 61341 61343 61344 \
	     61345 61346 61347 61348 61349 61350 61352 61353 61354 \
	     
	     61852 61854 61855 61856 61857 61858 61859 61860 61861 61862 \
	     
	     61914 61915 61916 61917 61918 61930 61931 61933 61934 61935 \
	     61936 61937 61938 61939 \
	     61947 61950 61951 61952 61953 61954 61955 61956"
}

set runlist "61321 61322 61323 61325 61327 61329 61330 61331 61332 61333 \
	     61334 61335 61336 61337 61340 61341 61343 61344 \
	     61345 61346 61347 61348 61349 61350 61352 61353 61354"


foreach run $runlist {
  set files [glob ${dir_mss}/Run0${run}/hd_rawdata_0${run}_???.evio]
  foreach file $files {
    set ext [string range $file end-13 end-5]
    if [file exists ${outdir}/rootFiles/${run}/${ext}.root] {continue}
    
    set run [string range $file end-13 end-9]
    if {($run < 7000)} {continue}
    
    set fnam [string range $file end-25 end]
    set jf $outdir/jsub/${ext}.jsub

    exec -- echo "PROJECT: gluex" > $jf
    exec -- echo "JOBNAME: compton_${ext}" >> $jf
    exec -- echo "COMMAND: $script $ext $fnam $jf $run" >> $jf
    exec -- echo "OS: centos77" >> $jf
    exec -- echo "TIME: 300" >> $jf
    exec -- echo "TRACK: analysis" >> $jf
    exec -- echo "MEMORY: 12000 MB" >> $jf
    exec -- echo "CPU:  5" >> $jf
    exec -- echo "DISK_SPACE: 23 GB" >> $jf
    exec -- echo "INPUT_FILES: " >> $jf
    exec -- echo "$file" >> $jf
    exec -- echo "SINGLE_JOB: true" >> $jf
    exec -- jsub $jf
    incr ij
    puts "submitting run $ext"
  }
}

puts "run $ij batch jobs"
