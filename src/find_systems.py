import numpy as np
from matplotlib import pyplot as plt

from scipy.stats import iqr
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

def rounded_means_quartiles(data):
    M25 =iqr(data, rng=(25, 50))
    M75 =iqr(data, rng=(50, 75))
    Mm  = np.median(data)
    return Mm, M25, M75

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

def find_multiple_complex_systems(stars):

    mplanet = 50 | units.MJupiter
    stars_done = Particles()
    
    N = {}
    mprim = {}
    msec = {}
    sma = {}
    ecc = {}
    planets = stars.copy()
    for bi in range(len(stars)):
        if stars[bi] in planets:
            planets.remove_particle(stars[bi])
        if len(planets)<=1:
            print("All combinations of stars done.")
            break
        if stars[bi] in stars_done:
            continue
        #stars_done.add_particle(stars[bi])
        total_masses = planets.mass + stars[bi].mass
        rel_pos = planets.position-stars[bi].position
        rel_vel = planets.velocity-stars[bi].velocity
        smai, ecci, true_anomaly,\
            inc, long_asc_node, arg_per_mat =\
                get_orbital_elements_from_arrays(rel_pos,
                                                 rel_vel,
                                                 total_masses,
                                                 G=constants.G)
        e = np.ma.masked_where(ecci>1, ecci)
        e = np.ma.masked_invalid(e)
        if(len(np.ma.compressed(e))>0):
            k = np.ma.masked_array(planets.key, e.mask)
            k = np.ma.compressed(k)
            p = Particles()
            for pi in planets:
                if pi.key in k:
                    p.add_particle(pi)
            a = np.ma.masked_array(smai, e.mask)
            a = np.ma.compressed(a)
            i = np.ma.masked_array(inc, e.mask)
            i = np.ma.compressed(i)
            e = np.ma.compressed(e)
            #sort on semi-major axis
            a, e, i, p = (list(t) for t in zip(*sorted(zip(a, e, i, p))))
            if stars[bi].mass>mplanet:
                name = "s"
            else:
                name = "p"
            for pi in range(len(p)):
                planets.remove_particle(p[pi])
                #stars_done.add_particle(p[pi])
                if a[pi]<1000|units.au:
                    #planets.remove_particle(p[pi])
                    print("N=", len(planets))
                    if p[pi].mass>mplanet:
                        name += "s"
                    else:
                        name += "p"
            name = "".join(sorted(name, reverse=True))
            
            if a[0]>25|units.au and a[0]<1000|units.au:
            #if a[0]<1000|units.au:
                if name in N:
                    N[name] += 1
                    sma[name].append(a[0])
                    ecc[name].append(e[0])
                    mprim[name].append(max(stars[bi].mass, p[0].mass))
                    msec[name].append(min(stars[bi].mass, p[0].mass))
                else:
                    N[name] = 1
                    sma[name] = [a[0].value_in(units.au)] | units.au
                    ecc[name] = [e[0]]
                    mprim[name] = [max(stars[bi].mass.value_in(units.MJupiter), p[0].mass.value_in(units.MJupiter))] | units.MJupiter
                    msec[name] = [min(stars[bi].mass.value_in(units.MJupiter), p[0].mass.value_in(units.MJupiter))] | units.MJupiter
            else:
                if "s" in N:
                    N["s"] += name.count("s")
                else:
                    N["s"] = name.count("s")
                if "p" in N:
                    N["p"] += name.count("p")
                else:
                    N["p"] = name.count("p")
                
    return N, sma, ecc, mprim, msec
                    
def print_multiple_data(binaries, typename):
    #print("type,key,k1,k2,M,m,a,e,i,dx,dy,dz,x,y,z,vx,vy,vz")
    for bi in binaries:
        print(f"{typename},{bi.key},{bi.id1},{bi.id2},{(bi.mass*bi.q).value_in(units.MJupiter)},{(bi.mass*(1-bi.q)).value_in(units.MJupiter)},{bi.semimajor_axis.value_in(units.au)},{bi.eccentricity},{bi.inclination.value_in(units.deg)},{bi.dx.value_in(units.au)},{bi.dy.value_in(units.au)},{bi.dz.value_in(units.au)},{bi.x.value_in(units.pc)},{bi.y.value_in(units.pc)},{bi.z.value_in(units.pc)},{bi.vx.value_in(units.kms)},{bi.vy.value_in(units.kms)},{bi.vz.value_in(units.kms)}")

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

    # read in all partilces
    bodies = read_set_from_file(o.filename, close_file=True)
    bodies.dx = 0 | units.au
    bodies.dy = 0 | units.au
    bodies.dz = 0 | units.au
    bodies.semimajor_axis = 0 | units.au
    bodies.eccentricity = 0 
    bodies.inclination = 0 | units.deg
    bodies.q = 1
    mplanet = 50 | units.MJupiter
    bodies = bodies[bodies.mass>0.6|units.MJupiter]

    stars = bodies.copy()
    stars.ncomponents = 1
    stars.id1 = 0
    stars.id2 = 0

    N, sma, ecc, mp, ms = find_multiple_complex_systems(stars)

    planets = stars[stars.mass<mplanet]
    jupiters = planets[planets.mass>0.6|units.MJupiter]
    stars = stars - planets
    Njup = len(jupiters)
    Nstr = len(stars)
        
    print(f"N={N}")
    
    for key, value in mp.items():
        print(f"{np.mean(value)} &", end=" ")
    for key, value in ms.items():
        print(f"{np.mean(value)} &", end=" ")
    for key, value in sma.items():
        print(f"{np.mean(value).value_in(units.au)} &", end=" ")
    for key, value in ecc.items():
        print(f"{np.mean(value)} &", end=" ")

    print(f"N={N}")
    for key in ["s", "ss", "sp", "p", "pp", 'spp']:
        if key not in N:
            N[key] = 0
        
    Np = Njup-N['sp']-2*N['pp']-2*N['spp']
    Ns = Nstr-N['sp']-2*N['ss']-N['spp']
    print(f"Ns={Ns}, Np={Np}")
    #print(f"{o.filename} & ", end=" ")
    n = 0
    for key, value in N.items():
        n += N[key] 
    for key in ["s", "ss", "sp", "p", "pp"]:
        n -= N[key] 
        print(N[key], "&", end="")
    print("\\\\")

    print(o.filename,"&", Ns, "&",  N['ss'], "&",  N['sp'], "&",  N['spp'], "&",  Np, "&",  N['pp'], "\\\\")

    key = 'pp'
    print(f"Orbital elements for {key} systems:")
    if N['pp']>0:
        print(f"{o.filename} & ", end=" ")
        M, dM25, dM75 = rounded_means_quartiles(mp[key].value_in(units.MJupiter))
        print(f"${M:2.2f}^{{+{dM25:2.2f}}}_{{-{dM75:2.2f}}}$ &", end=" ")
        M, dM25, dM75 = rounded_means_quartiles(ms[key].value_in(units.MJupiter))
        print(f"${M:2.2f}^{{+{dM25:2.2f}}}_{{-{dM75:2.2f}}}$ &", end=" ")
        a, da25, da75 = rounded_means_quartiles(sma[key].value_in(units.au))
        print(f"${a:2.2f}^{{+{da25:2.2f}}}_{{-{da75:2.2f}}}$ &", end=" ")
        e, de25, de75 = rounded_means_quartiles(ecc[key])
        print(f"${e:2.2f}^{{+{de25:2.2f}}}_{{-{de75:2.2f}}}$ \\\\")


    key = 'ss'
    print(f"Orbital elements for {key} systems:")
    if N['pp']>0:
        print(f"{o.filename} & ", end=" ")
        M, dM25, dM75 = rounded_means_quartiles(mp[key].value_in(units.MJupiter))
        print(f"${M:2.2f}^{{+{dM25:2.2f}}}_{{-{dM75:2.2f}}}$ &", end=" ")
        M, dM25, dM75 = rounded_means_quartiles(ms[key].value_in(units.MJupiter))
        print(f"${M:2.2f}^{{+{dM25:2.2f}}}_{{-{dM75:2.2f}}}$ &", end=" ")
        a, da25, da75 = rounded_means_quartiles(sma[key].value_in(units.au))
        print(f"${a:2.2f}^{{+{da25:2.2f}}}_{{-{da75:2.2f}}}$ &", end=" ")
        e, de25, de75 = rounded_means_quartiles(ecc[key])
        print(f"${e:2.2f}^{{+{de25:2.2f}}}_{{-{de75:2.2f}}}$ \\\\")

        
