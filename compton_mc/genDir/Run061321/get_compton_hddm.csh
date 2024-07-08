#!/usr/bin/tcsh
#

set runnumber = "61321"

set tag_sys = $1

set counter = 1
while ( $counter <= 274 )
  
  set ofile = "test.root"
  
  set ofile = "${tag_sys}_${counter}"
  if( $counter < 10 ) then
    set ofile = "${tag_sys}_00${counter}"
  else if ( $counter < 100 ) then
    set ofile = "${tag_sys}_0${counter}"
  endif
  
  set pwd   = `pwd`
  set compf = $pwd/compton.cfg
  
  if ( -f $pwd/genRootFiles/${ofile}.root ) then
    
    echo "run: false" > $compf
    echo "shell: tcsh" >> $compf
    echo "dir: " >> $compf
    echo "file: $pwd/genRootFiles/${ofile}.root" >> $compf
    echo "out_dir: " >> $compf
    echo "workflow: " >> $compf
    
    echo "${tag_sys} counter ${counter}"
    
    echo "gen_primex_compton -e compton.cfg -hd output.hddm -r ${runnumber} -s 0"
    gen_primex_compton -e compton.cfg -hd output.hddm -r ${runnumber} -s 0
    
    mv output.hddm ${ofile}.hddm
    
    rm compton.cfg
    rm gen_primex_compton_runnb_${runnumber}.root
    
  else 
    
	echo "${tag_sys} counter ${counter} (Skipped)"
	
  endif
  
  set counter  = `expr $counter + 1`
  
end
