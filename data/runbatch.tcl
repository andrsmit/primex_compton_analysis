#!/usr/bin/tclsh

set ext  [lindex $argv 0]
set file [lindex $argv 1]
set jf   [lindex $argv 2]
set run  [lindex $argv 3]

set odir /work/halld/home/andrsmit/primex_compton_analysis/data

set env(HALLD_RECON_HOME) /work/halld/home/andrsmit/halld_software/halld_recon

set env(CCDB_CONNECTION) mysql://ccdb_user@hallddb.jlab.org/ccdb
set env(JANA_CALIB_URL)  mysql://ccdb_user@hallddb.jlab.org/ccdb
set env(JANA_CALIB_CONTEXT) variation=primex_cal

set bin /work/halld/home/andrsmit/halld_software/halld_recon/Linux_CentOS7.7-x86_64-gcc4.8.5/bin/hd_root

exec -- rm $file
exec -- cp /cache/mss/halld/RunPeriod-2019-01/rawdata/Run0${run}/${file} .


if [catch {exec -- $bin -PPLUGINS=compton_analysis -PNTHREADS=5 -PCOMPTON_ANALYSIS:RUN_GROUP=0 $file >& /dev/null} err] {
  exec -- echo $err > ${odir}/log/${ext}.err
}

if [file exists hd_root.root] {
if {([file size hd_root.root] < 5000)} {
  puts stderr "job was not finished normally for ext $ext"
} else {
  catch {exec -- mv  hd_root.root $odir/rootFiles/${run}/${ext}.root}
  catch {exec -- rm $jf}
}
} else {
    puts stderr "job was not finished normally for ext $ext"
}

exec -- rm $file
