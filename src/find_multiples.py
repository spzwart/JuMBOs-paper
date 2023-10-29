import numpy as np
from matplotlib import pyplot as plt

from amuse.units import units
from amuse.units import nbody_system
from amuse.ic.plummer import new_plummer_model
from amuse.community.ph4.interface import ph4
from amuse.plot import scatter
from amuse.io import store
from amuse.community.seba.interface import SeBa
from amuse.ext.LagrangianRadii import LagrangianRadii
from amuse import datamodel
from amuse.ext.evrard_test import uniform_unit_sphere
from amuse.lab import new_kroupa_mass_distribution
from amuse.lab import read_set_from_file
from amuse.lab import Particle, Particles, constants
from amuse.lab import write_set_to_file, read_set_from_file

from amuse.ext.orbital_elements import get_orbital_elements_from_arrays
from amuse.ext.orbital_elements import orbital_elements_from_binary

def is_this_binary_planet_orbiting_a_star(stars, planet):

    if len(planet.companions)>0:
        multiplanet = construct_center_of_mass_particle(planet)
    else:
        return []|units.au, []
        
    total_masses = multiplanet.mass + stars.mass
    rel_pos = multiplanet.position-stars.position
    rel_vel = multiplanet.velocity-stars.velocity
    sma, ecc, true_anomaly,\
        inc, long_asc_node, arg_per_mat =\
            get_orbital_elements_from_arrays(rel_pos,
                                             rel_vel,
                                             total_masses,
                                             G=constants.G)

    a_mp = [] | units.au
    e_mp = [] 
    n_starbinaryplanet = 0
    for i in range(len(sma)):
        if ecc[i]>0 and ecc[i]<1 and sma[i]<0.01|units.au:
            #print("Bound Jumbos:", sma[i].in_(units.au), ecc[i])
            #print("planet with M=", multiplanet.mass.in_(units.MJupiter))
            #print("companions:", sma[i].in_(units.au), "e= ", ecc[i])
            n_starbinaryplanet += 1
            a_mp.append(sma[i])
            e_mp.append(ecc[i])
    #print(f"Nstar-binary-planet = {n_starbinaryplanet}")
    return a_mp, e_mp

def free_floating_planets(bodies):
    single_freefloatrs = bodies[bodies.type=="planet"].copy()
    for bi in bodies:
        if len(bi.companions)>0:
            for ci in bi.companions:
                if ci in single_freefloaters:
                    if ci.semimajor_axis<0.01|units.pc:
                        single_freefloaters.remove_particle(ci)
                if bi.type=="planet" and bi in single_freefloaters:
                    single_freefloaters.remove_particle(bi)
    return single_freefloaters

def construct_center_of_mass_from_particle_pair(bi, ci):
    if bi.mass<ci.mass:
        ai = ci
        ci = bi
        bi = ai
    p = Particles(2)
    p[0].mass = bi.mass
    p[0].position = bi.position
    p[0].velocity = bi.velocity
    p[1].mass = ci.mass
    p[1].position = ci.position
    p[1].velocity = ci.velocity
    cm = Particle()
    cm.name = f"({bi.name}+{ci.name})"
    cm.type = "com"
    cm.id1 = bi.key
    cm.id2 = ci.key
    cm.ncomponents = bi.ncomponents + ci.ncomponents
    cm.mass = p.mass.sum()
    cm.q = p[1].mass/p[0].mass
    cm.position = p.center_of_mass()
    cm.velocity = p.center_of_mass_velocity()
    return cm

def construct_center_of_mass_particle(body):
    p = Particles()
    p.add_particle(body)
    p.add_particles(body.companions)
    cm = Particle()
    cm.type = "com"
    cm.ncomponents = p.ncomponents.sum()
    cm.mass = p.mass.sum()
    cm.position = p.center_of_mass()
    cm.velocity = p.center_of_mass_velocity()
    return cm

def find_companions(stars, planets):

    binaries = Particles()
    components = Particles()
    for bi in range(len(stars)):
        if stars[bi] in planets:
            plts = planets - stars[bi]
        else:
            plts = planets
        total_masses = plts.mass + stars[bi].mass
        rel_pos = plts.position-stars[bi].position
        rel_vel = plts.velocity-stars[bi].velocity
        sma, ecc, true_anomaly,\
            inc, long_asc_node, arg_per_mat =\
                get_orbital_elements_from_arrays(rel_pos,
                                                 rel_vel,
                                                 total_masses,
                                                 G=constants.G)
        e = np.ma.masked_where(ecc>1, ecc)
        e = np.ma.masked_invalid(e)
        if(len(np.ma.compressed(e))>0):
            k = np.ma.masked_array(plts.key, e.mask)
            k = np.ma.compressed(k)
            p = Particles()
            for pi in plts:
                if pi.key in k:
                    p.add_particle(pi)
            #p = planets[planets.key==k] # does not seem to work properly?
            a = np.ma.masked_array(sma, e.mask)
            a = np.ma.compressed(a)
            i = np.ma.masked_array(inc, e.mask)
            i = np.ma.compressed(i)
            e = np.ma.compressed(e)
            #sort on semi-major axis
            a, e, i, p = (list(t) for t in zip(*sorted(zip(a, e, i, p))))
            
            #print("companion found:", e, len(p))
            if len(p)>=1:
                if a[0]<1000|units.au:
                    b = Particle()
                    ai = a[0]
                    ii = i[0]
                    ei = e[0]
                    b = construct_center_of_mass_from_particle_pair(stars[bi], p[0])
                    b.semimajor_axis = ai
                    b.eccentricity = ei
                    b.inclination = ii
                    if p[0].key in components.key:
                        #print("already done")
                        pass
                    else:
                        #print(f"a={ai.in_(units.au)}, e={ei}: {stars[bi].name} ({stars[bi].mass.in_(units.MSun)}), {p[0].name} ({p[0].mass.in_(units.MSun)})")
                        components.add_particle(stars[bi])
                        components.add_particle(p[0])
                        binaries.add_particle(b)
    return binaries, components

def print_multiple_data(binaries, typename):
    #print("type,key,k1,k2,M,m,a,e,i,x,y,z,vx,vy,vz")
    for bi in binaries:
        print(f"{typename},{bi.key},{bi.id1},{bi.id2},{(bi.mass*bi.q).value_in(units.MJupiter)},{(bi.mass*(1-bi.q)).value_in(units.MJupiter)},{bi.semimajor_axis.value_in(units.au)},{bi.eccentricity},{bi.inclination.value_in(units.deg)}, {bi.x.value_in(units.pc)},{bi.y.value_in(units.pc)},{bi.z.value_in(units.pc)},{bi.vx.value_in(units.kms)},{bi.vy.value_in(units.kms)},{bi.vz.value_in(units.kms)}")

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-f", 
                      dest="filename", default = "jumbos_i0000.amuse",
                      help="input filename [%default]")
    result.add_option("-p", action='store_true',
                      dest="plot", default="False",
                      help="show figures [%default]")
    return result
    
if __name__=="__main__":
    o, arguments  = new_option_parser().parse_args()

    outfile = "processed_"+o.filename
    
    # read in all partilces
    bodies = read_set_from_file(o.filename, close_file=True)
    #for bi in bodies:
    #    bi.companions = Particles()
    #    bi.parent = Particles()

    stars = bodies.copy()
    stars.ncomponents = 1
    stars.id1 = 0
    stars.id2 = 0

    binaries, components = find_companions(stars, stars)
    #print(f"N binaries: {len(binaries)}")
    stars.remove_particles(components)
    stars.add_particles(binaries)
    print_multiple_data(binaries, "binary")
    
    triples, triple_components = find_companions(binaries, stars)
    #print(f"N triples: {len(triples)}")
    stars.remove_particles(triple_components)
    stars.add_particles(triples)
    print_multiple_data(triples, "triple")
    
    binary_binary, quadruple_components = find_companions(binaries, binaries)
    #print(f"N binary-binary: {len(binary_binary)}")
    for bi in triple_components:
        if bi.key in quadruple_components:
            print(f"star in triple as well as in quadruple: bi")
    print_multiple_data(binary_binary, "binbin")
            
    triple_single, components = find_companions(triples, stars)
    #print(f"N triples-single: {len(triple_single)}")
    print_multiple_data(triple_single, "tripsing")
    
    
