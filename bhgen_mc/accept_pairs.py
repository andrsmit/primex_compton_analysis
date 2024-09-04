import hddm_s
import ROOT
import random
import math
import sys

DRAW_HISTS = 0

infile_name = sys.argv[1]

h_ngen       = ROOT.TH1F("ngen", 
    "Number e+e- Pairs Generated; E_{#gamma} (GeV); N_{gen}", 100, 2, 12)
h_pair_cs    = ROOT.TH1F("pair_cs", 
    "e+e- Pair CS; E_{#gamma} (GeV); #sigma (#mub)", 100, 2, 12)
h_triplet_cs = ROOT.TH1F("triplet_cs", 
    "e+e- Triplet CS; E_{#gamma} (GeV); #sigma (#mub)", 100, 2, 12)

h_max_weight = ROOT.TH1F("max_weight", "Max Event Weight; E_{#gamma} (GeV)", 100, 2, 12)

n_gen_total = 0.0

weight_sum         = 0.0
weight_sum_pair    = 0.0
weight_sum_triplet = 0.0

n_nans = 0

#------------------------------------------------------------------------------------------------#
# Loop over all event records to get average cross section and maximum event weight for filtering:

n_event = 0

for rec in hddm_s.istream(infile_name):
    
    if n_event % 100000 == 0:
        print("    processing event",n_event)
    
    for rea in rec.getReactions():
        
        beams = rea.getBeams()
        if len(beams) != 1:
            continue
        
        E = rea.getBeams()[0].getMomenta()[0].E
        #E = rec.getMomenta()[0].E
        n_gen_total = n_gen_total + 1.0
        h_ngen.Fill(E)
        
        iE = h_ngen.FindBin(E)
        
        if math.isnan(rea.weight):
            #print("Reaction weight is not a number. Reaction type: ",rea.type)
            n_nans = n_nans+1
            continue
        
        weight_sum += rea.weight
        if rea.weight > h_max_weight.GetBinContent(iE):
            h_max_weight.SetBinContent(iE, rea.weight)
        
        if(rea.type == 221):
            weight_sum_pair += rea.weight
            h_pair_cs.Fill(E, rea.weight)
        elif(rea.type == 231):
            weight_sum_triplet += rea.weight
            h_triplet_cs.Fill(E, rea.weight)
        #else:
        #    print("reaction weight = ",rea.weight)
    
    n_event = n_event+1

print("  total cs =",weight_sum/n_gen_total)
print("  pair cs =",weight_sum_pair/n_gen_total)
print("  triplet cs =",weight_sum_triplet/n_gen_total)
print("  number of nan weights = ",n_nans)

n_pairs_thrown = 0
n_pairs_kept   = 0

n_triplets_thrown = 0
n_triplets_kept   = 0

h_pair_cs.Divide(h_ngen)
h_triplet_cs.Divide(h_ngen)

if DRAW_HISTS > 0:
    cCS = ROOT.TCanvas("cCS", "CS", 1000, 500)
    cCS.Divide(2,1)
    cCS.cd(1)
    h_pair_cs.Draw()
    cCS.cd(2)
    h_triplet_cs.Draw()

#------------------------------------------------------------------------------------------------#
# Do a second loop over the events to determine whether each event should be saved:

fout = hddm_s.ostream("BHgen_accepted.hddm")

#fout_pair    = hddm_s.ostream("BHgen_pair_filtered.hddm")
#fout_triplet = hddm_s.ostream("BHgen_triplet_filtered.hddm")

# histogram to show the flux of events that survived:

gen_flux = ROOT.TH1F("gen_flux", "Generated Flux", 100, 2., 12.)
saved_flux = ROOT.TH1F("saved_flux", "Saved Flux", 100, 2., 12.)

n_events_written = 0

for rec in hddm_s.istream(infile_name):
    for rea in rec.getReactions():
        
        if math.isnan(rea.weight):
            continue
        
        beams = rea.getBeams()
        if len(beams) != 1:
            continue
        
        E = rea.getBeams()[0].getMomenta()[0].E
		#E = rec.getMomenta()[0].E
        iE = h_ngen.FindBin(E)
        wmax = h_max_weight.GetBinContent(iE)
        
        gen_flux.Fill(E)
        
        write_val = 0
        if rea.type == 221:
            n_pairs_thrown = n_pairs_thrown+1
            if rea.weight > random.random() * wmax:
                n_pairs_kept = n_pairs_kept+1
                rea.weight = h_pair_cs.GetBinContent(iE)
                write_val = 1
                saved_flux.Fill(E)
            
        elif rea.type == 231:
            n_triplets_thrown = n_triplets_thrown+1
            if rea.weight > random.random() * wmax:
                n_triplets_kept = n_triplets_kept+1
                rea.weight = h_triplet_cs.GetBinContent(iE)
                write_val = 1
                saved_flux.Fill(E)
        
        # Write output to new file, removing the primary vertex from the beam photon:
        if write_val == 0:
            continue
        
        #fout.write(rec)
        
        new_rec = hddm_s.HDDM()
        
        #-----------------------------------------------------------------------------#
        # Get objects from the old event record needed for writing new record:
        
        old_pev_list = rec.getPhysicsEvents()
        if len(old_pev_list) != 1:
            print("",len(old_pev_list),"physics event(s) found")
            continue
        old_pev = old_pev_list[0]
        
        old_vertex_list = rea.getVertices()
        if len(old_vertex_list) != 2:
            print("",len(old_vertex_list),"vertex(s) found at event number",old_pev.eventNo)
            continue
        old_pair_vertex = old_vertex_list[1]
        old_pair_origin = old_pair_vertex.getOrigins()[0]
        
        old_target_list = rea.getTargets()
        if len(old_target_list) != 1:
            print("",len(old_target_list),"target(s) found at event number",old_pev.eventNo)
            continue
        old_target = old_target_list[0].getMomenta()[0]
        
        old_beam_list = rea.getBeams()
        if len(old_beam_list) != 1:
            print("",len(old_beam_list),"beam photon(s) found at event number",old_pev.eventNo)
            continue
        old_beam_photon = old_beam_list[0].getMomenta()[0]
        
        old_seed_list = rea.getRandoms()
        if len(old_seed_list) != 1:
            print("",len(old_seed_list),"seed(s) found at event number",old_pev.eventNo)
            continue
        old_seeds = old_seed_list[0]
        
        old_product_list = old_pair_vertex.getProducts()
        
        # set target type according to run number (default is proton):
        runNo = old_pev.runNo
        target_type = 14
        if(runNo>60000 and runNo<69999):
            if(runNo<61355):
                target_type = 64
            else:
                target_type = 47
        elif(runNo>80000 and runNo<89999):
            if(runNo<81396):
                target_type = 64
            else:
                target_type = 47
        elif(runNo>110000 and runNo<119999):
            if(runNo<110622):
                target_type = 64
            else:
                target_type = 47
        
        #-----------------------------------------------------------------------------#
        # Write event:
        
        # create physics event:
        new_pev = new_rec.addPhysicsEvents(1)
        new_pev[0].runNo   = runNo
        new_pev[0].eventNo = n_events_written + 1
        
        # create reaction:
        new_rea = new_pev[0].addReactions(1)
        new_rea[0].type   = rea.type
        new_rea[0].weight = rea.weight
        
	    # create target:
        new_tar = new_rea[0].addTargets(1)
        new_tar[0].type = target_type
        new_tmo = new_tar[0].addMomenta(1)
        new_tmo[0].E  = old_target.E
        new_tmo[0].px = old_target.px
        new_tmo[0].py = old_target.py
        new_tmo[0].pz = old_target.pz
        
        # create beam photon:
        new_beam = new_rea[0].addBeams(1)
        new_beam[0].type = 1
        new_beam_mom = new_beam[0].addMomenta(1)
        new_beam_mom[0].E  = old_beam_photon.E
        new_beam_mom[0].px = old_beam_photon.px
        new_beam_mom[0].py = old_beam_photon.py
        new_beam_mom[0].pz = old_beam_photon.pz
        new_beam_pol = new_beam[0].addPolarizations(1)
        new_beam_pol[0].Px = 0
        new_beam_pol[0].Py = 0
        new_beam_pol[0].Pz = 0
        new_beam_prop = new_beam[0].addPropertiesList(1)
        new_beam_prop[0].charge = 0
        new_beam_prop[0].mass   = 0
        
        # create seeds:
        new_seeds = new_rea[0].addRandoms(1)
        new_seeds[0].seed1 = old_seeds.seed1
        new_seeds[0].seed2 = old_seeds.seed2
        new_seeds[0].seed3 = old_seeds.seed3
        new_seeds[0].seed4 = old_seeds.seed4
        
        # create vertex:
        new_ver = new_rea[0].addVertices(1)
        
        # create origin:
        new_ori = new_ver[0].addOrigins(1)
        new_ori[0].t  = old_pair_origin.t
        new_ori[0].vx = old_pair_origin.vx
        new_ori[0].vy = old_pair_origin.vy
        new_ori[0].vz = old_pair_origin.vz
        
        # create products:
        new_pro = new_ver[0].addProducts(len(old_product_list))
        for ipro in range(len(old_product_list)):
            
            # particle information:
            new_pro[ipro].decayVertex = 0
            new_pro[ipro].id = ipro+1
            new_pro[ipro].mech = old_product_list[ipro].mech
            new_pro[ipro].parentid = 1
            old_pdgtype = old_product_list[ipro].pdgtype
            new_pro[ipro].pdgtype = old_pdgtype
            
            new_g3type = 0
            if old_pdgtype == -11:
                new_g3type = 2 # positron
            elif old_pdgtype == 11:
                new_g3type = 3 # electron
            elif old_pdgtype == 1000040090:
                new_g3type = 64 # Beryllium-9
            elif old_pdgtype == 2212:
                new_g3type = 14 # Proton
            elif old_pdgtype == 2112:
                new_g3type = 13 # Neutron
            elif old_pdgtype == 22:
                new_g3type = 1 # Photon
            elif old_pdgtype == 1000020040:
                new_g3type = 47 # Helium-4
            else:
                print(" Unrecognized particle type: ",old_pdgtype,",",old_product_list[ipro].type)
            
            new_pro[ipro].type = new_g3type
            
            # add momentum for particle:
            old_pro_mom = old_product_list[ipro].getMomenta()[0]
            new_pro_mom = new_pro[ipro].addMomenta(1)
            new_pro_mom[0].E  = old_pro_mom.E
            new_pro_mom[0].px = old_pro_mom.px
            new_pro_mom[0].py = old_pro_mom.py
            new_pro_mom[0].pz = old_pro_mom.pz
        
        fout.write(new_rec)
        n_events_written = n_events_written + 1

print("Wrote",n_pairs_kept," out of",n_pairs_thrown," pairs (",(n_pairs_thrown/n_pairs_kept),")")
print("Wrote",n_triplets_kept," out of",n_triplets_thrown," triplets (",(n_triplets_thrown/n_triplets_kept),")")

if DRAW_HISTS > 0:
    cFlux = ROOT.TCanvas("cFlux", "Flux", 1000, 500)
    cFlux.Divide(2,1)
    cFlux.cd(1)
    gen_flux.Draw()
    cFlux.cd(2)
    saved_flux.Draw()
    input()
