#!/usr/bin/tclsh
#
set tcl_precision 15
set outdir /work/halld/home/andrsmit/compton_analysis/tag_sim
set script ${outdir}/sim_comp.csh
set ij 0


set genDir /work/halld/home/andrsmit/compton_analysis/tag_sim/genHddmFiles


set taglist "tagh"

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
		 101 102 103 104 105 106 107 108 109 110 \
		 111 112 113 114 115 116 117 118 119 120 \
		 121 122 123 124 125 126 127     179 180 \
		 181 182 183 184 185 186 187 188 189 190 \
		 191 192 193 194 195 196 197 198 199 200 \
		 201 202 203 204 205 206 207 208 209 210 \
		 211 212 213 214 215 216 217 218 219 220 \
		 221"


foreach tag_sys $taglist {
  foreach counter $counterlist {
    
    if [file exists ${genDir}/${tag_sys}_${counter}.hddm] {
      
      set jf $outdir/jsub/${tag_sys}_${counter}.jsub
      
      if [file exists ${outdir}/recHddmFiles_v1/${tag_sys}_${counter}.hddm] {continue}
      
      exec -- echo "PROJECT: gluex" > $jf
      exec -- echo "JOBNAME: comp_sim_${tag_sys}_${counter}" >> $jf
      exec -- echo "COMMAND: $script ${tag_sys} ${counter} $jf $outdir" >> $jf
      exec -- echo "OS: centos77" >> $jf
      exec -- echo "TIME: 400" >> $jf
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
