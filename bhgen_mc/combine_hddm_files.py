import hddm_s
import sys
from pathlib import Path

# command line arguments are the list of files to combine

if len(sys.argv) < 2:
    print("check arguments")
    exit

output_fname = "combined_hddm.hddm"

loc_fout = hddm_s.ostream(output_fname)
#loc_fout.compression = hddm_s.k_z_compression

eventNo = 1

summed_files = 0

for ifile in range(1,len(sys.argv)):
    loc_fname = sys.argv[ifile]
    loc_file = Path(loc_fname)
    if loc_file.is_file():
        print("",loc_fname)
        summed_files  = summed_files + 1
        for rec in hddm_s.istream(loc_fname):
            rec.getPhysicsEvents()[0].eventNo = eventNo
            loc_fout.write(rec)
            eventNo = eventNo+1

print(" Summed",summed_files,"input hddm files into file:",output_fname)
