#!/usr/bin/tclsh
#
set tcl_precision 15
set outdir /work/halld/home/andrsmit/compton_analysis/tag_sim
set script ${outdir}/run_rec.csh
set ij 0


set taglist "tagm"

set counterlist "001 002 003 004 005 006 007 008 009 010 \
	         011 012 013 014 015 016 017 018 019 020 \
		 021 022 023 024 025 026 027 028 029 030 \
		 031 032 033 034 035 036 037 038 039 040 \
		 041 042 043 044 045 046 047 048 049 050 \
		 051 052 053 054 055 056 057 058 059 060 \
		 061 062 063 064 065 066 067 068 069 070 \
		 071 072 073 074 075 076 077 078 079 080 \
		 081 082 083 084 085 086 087 088 089 090 \
		 091 092 093 094 095 096 097 098 099 100 \
		 101 102"


foreach tag_sys $taglist {
  foreach counter $counterlist {
    
    if [file exists ${outdir}/recHddmFiles_v1/${tag_sys}_${counter}.hddm] {
      
      set jf $outdir/jsub/${tag_sys}_${counter}_rec.jsub
      
      if [file exists ${outdir}/recRootFiles_v1_sys/${tag_sys}_${counter}.root] {continue}
      
      exec -- echo "PROJECT: gluex" > $jf
      exec -- echo "JOBNAME: comp_rec_${tag_sys}_${counter}" >> $jf
      exec -- echo "COMMAND: $script $tag_sys $counter $jf $outdir" >> $jf
      exec -- echo "OS: centos77" >> $jf
      exec -- echo "TIME: 300" >> $jf
      exec -- echo "TRACK: simulation" >> $jf
      exec -- echo "MEMORY: 8000 MB" >> $jf
      exec -- echo "CPU:  5" >> $jf
      exec -- echo "DISK_SPACE: 10 GB" >> $jf
      exec -- echo "SINGLE_JOB: true" >> $jf
      exec -- jsub $jf
      incr ij
      puts "submitting run ${tag_sys}_${counter}"
         
    }
  }
}

puts "run $ij batch jobs"
