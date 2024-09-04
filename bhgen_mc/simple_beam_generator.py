import hddm_s
import random
import sys
import ROOT

nevts = int(sys.argv[1])
runNo = int(sys.argv[2])

DRAW_HISTS = 0

# set random seed for photon beam generator:
random.seed()

ROOT.gRandom.SetSeed(0)

flux_file_name = "/work/halld/home/andrsmit/primex_compton_analysis/bhgen_mc/photon_flux/%06d_flux.root" % runNo
print("Reading photon flux from: ",flux_file_name)

# set target according to run number:
target_type = 14
target_mass = 0.938272046
if(runNo>60000 and runNo<69999):
    if(runNo<61355):
        target_type = 64
        target_mass = 8.39479
    else:
        target_type = 47
        target_mass = 3.727379238
elif(runNo>80000 and runNo<89999):
    if(runNo<81396):
        target_type = 64
        target_mass = 8.39479
    else:
        target_type = 47
        target_mass = 3.727379238
elif(runNo>110000 and runNo<119999):
    if(runNo<110622):
        target_type = 64
        target_mass = 8.39479
    else:
        target_type = 47
        target_mass = 3.727379238

# check if file name exists:

flux_hist = ROOT.TH1F("flux", "", 12000, 0., 12.)
if ROOT.gSystem.AccessPathName(flux_file_name):
    print("Tagged flux ROOT file for run",runNo,"does not exist.")
    for ibin in range(flux_hist.GetXaxis().GetNbins()):
        loc_eb = flux_hist.GetBinCenter(ibin+1)
        if loc_eb > 6.0 and loc_eb < 11.2:
            flux_hist.SetBinContent(ibin+1, 1.0)
else:
    # get flux histogram from file;
    flux_root_file = ROOT.TFile(flux_file_name, "READ")
    flux_hist      = flux_root_file.Get("tagged_flux")
    flux_hist.SetDirectory(0)
    flux_root_file.Close()

h_generated_flux = ROOT.TH1F("generated_flux", "", 1200, 0., 12.)

fout = hddm_s.ostream("photonbeam.hddm")
for eventNo in range(nevts):
    
    if eventNo % 1000000 == 0:
        print("  processing event",eventNo)
    
    eb = flux_hist.GetRandom()
    
    # randomly generate vertex position from beam spot:
    vertex_x = random.gauss(0.027, 0.255)
    vertex_y = random.gauss(-0.128, 0.255)
    
    # be target is 1.27 cm in radius 
    if abs(vertex_x) >= 1.27:
        vertex_x = 0.027
    if abs(vertex_y) >= 1.27:
        vertex_y = -0.128
    
    # event record:
    rec = hddm_s.HDDM()
    
    # create physics event:
    pev = rec.addPhysicsEvents(1)
    pev[0].runNo = runNo
    pev[0].eventNo = eventNo+1
    
    # create reaction:
    rea = pev[0].addReactions(1)
    rea[0].weight = 1.0
    
    # create vertex:
    ver = rea[0].addVertices(1)
    
    # vertex product (photon):
    pro = ver[0].addProducts(1)
    pro[0].id = 1
    pro[0].pdgtype = 22
    pro[0].type = 1
    
    # photon momentum:
    mom = pro[0].addMomenta(1)
    mom[0].E = eb
    mom[0].px = 0
    mom[0].py = 0
    mom[0].pz = mom[0].E
    
    # photon polarization:
    pol = pro[0].addPolarizations(1)
    pol[0].Px = 0
    pol[0].Py = 0
    pol[0].Pz = 0
    
    # vertex origin:
    ori = ver[0].addOrigins(1)
    ori[0].vz = 50
    ori[0].vx = vertex_x
    ori[0].vy = vertex_y
    
    # create target:
    tar = rea[0].addTargets(1)
    
    #Be-9:
    tar[0].type = target_type
    tmo = tar[0].addMomenta(1)
    tmo[0].E = target_mass
    
    # write record:
    fout.write(rec)
    
    h_generated_flux.Fill(eb)

if DRAW_HISTS > 0:
    c0 = ROOT.TCanvas("c0","c0",500,500)
    c0.cd()
    h_generated_flux.Draw("hist")
    input()
