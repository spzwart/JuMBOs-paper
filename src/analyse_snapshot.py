import numpy as np
from math import sqrt

from matplotlib import pyplot
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
                                                 
import random

from matplotlib import pyplot as plt
def orbital_elements_of_binary(primary, secondary):
    b = Particles(1)
    b[0].mass = primary.mass
    b[0].position = primary.position
    b[0].velocity = primary.velocity
    b.add_particle(secondary)
    #print(b.mass.in_(units.MSun), b.position.in_(units.au), b.velocity.in_(units.kms), b.mass.in_(units.MSun))

    from amuse.ext.orbital_elements import orbital_elements_from_binary
    #b = read_set_from_file(filename+"_binary.amuse", "hdf5", close_file=True)
    M, m, a, e, ta_out, inc, lan_out, aop_out = orbital_elements_from_binary(b, G=constants.G)
    #print("orbit=", M, m, a, e, ta_out, inc, lan_out, aop_out)
    b.semimajor_axis = a
    b.eccentricity = e
    b.ta_out = ta_out
    b.inclination = inc
    b.lan_out = lan_out
    b.aop_out = aop_out
    return a, e

def find_stellar_companion_to_jumbos(jumbo, stars):

    r = np.zeros(len(stars)) | units.pc
    for si in range(len(stars)):
        r[si] = (stars[si].position-jumbo.position).length()
    r, s = (list(t) for t in zip(*sorted(zip(r, stars))))
    nn = s[0].as_set().get_intersecting_subset_in(stars)
    a, e, = orbital_elements_of_binary(stars[si], nn)
    print(f"Stellar orbit: {a.in_(units.au)} e={e}")
    return a, e
    
def find_host_stellar_companion(jumbos, stars):

    triples = Particles()
    singles = Particles()
    sma = [] | units.au
    ecc = []
    for ji in jumbos:
        if hasattr(ji, "child1"):
            print(f"try jumbo N={1+len(sma)}")
            r = np.zeros(len(stars)) | units.pc
            for si in range(len(stars)):
                r[si] = (stars[si].position-ji.position).length()
            r, s = (list(t) for t in zip(*sorted(zip(r, stars))))
            nn = s[0].as_set().get_intersecting_subset_in(stars)[0]
            a, e, = orbital_elements_of_binary(ji, nn)
            if e<1:
                print(f"Jumbo Orbit around star: a={a.in_(units.au)}, {e}")
                print(nn)
                nn.planets.add_particle(ji)
                triples.add_particle(nn)
                sma.append(a)
                ecc.append(e)
            else:
                singles.add_particles(ji)
            print(f"N(SJ)= {len(triples)}, N(sp)={len(singles)}")
    #print("sma=", sma)
    #print(ecc)
    plt.title("triple parameters: (s,j)")
    plt.scatter(sma.value_in(units.au), ecc, c='b')
    plt.show()
    return triples, sma, ecc, singles #all_jumbos

def find_binary_planets(bodies):

    jumbos = bodies[bodies.type=="planet"]
    stars = bodies - jumbos
    all_jumbos = Particles()
    single_planets = jumbos.copy()
    #find nearest jumbo
    sma = [] | units.au
    ecc = []
    #jdone = Particles()
    for ji in jumbos:
        #jdone.add_particle(ji)
        b = (jumbos-ji).copy()
        b.position-=ji.position
        b.velocity-=ji.velocity
        r = np.zeros(len(b)) | units.pc
        for bi in range(len(b)):
            r[bi] = b[bi].position.length()
        r, b = (list(t) for t in zip(*sorted(zip(r, b))))
        nn = b[0].as_set().get_intersecting_subset_in(bodies)
        a, e, = orbital_elements_of_binary(ji, nn)
        if e<1:
            print(f"Jumbo Orbit: a={a.in_(units.au)}, {e}")
            tmp = Particles(1)
            tmp.add_particle(ji)
            tmp.add_particle(nn)
            jbb = Particle()
            jbb.mass = tmp.mass.sum()
            jbb.position = tmp.center_of_mass()
            jbb.velocity = tmp.center_of_mass_velocity()
            jbb.child1 = Particles(0)
            jbb.child2 = Particles(0)
            jbb.child1.add_particle(ji)
            jbb.child2.add_particle(nn)
            if ji in single_planets:
                single_planets.remove_particle(ji)
            if nn in single_planets:
                single_planets.remove_particle(nn)
            all_jumbos.add_particle(jbb)
            sma.append(a)
            ecc.append(e)
    print(f"N(JJ)= {len(sma)}")
    print(sma)
    print(ecc)
    plt.scatter(sma.value_in(units.au), ecc, c='b')
    plt.semilogx()
    plt.title("binary planets (Jumbos)")
    plt.show()
    return all_jumbos, single_planets

def print_planetary_orbits(bodies):

    jumbos = bodies[bodies.type=="planet"]
    stars = bodies[bodies.name=="host"]
    primaries = Particles()
    sma = [] | units.au
    ecc = []
    for si in stars:
        total_masses = jumbos.mass + si.mass
        rel_pos = jumbos.position-si.position
        rel_vel = jumbos.velocity-si.velocity
        a, e, true_anomaly,\
            inc, long_asc_node, arg_per_mat =\
                get_orbital_elements_from_arrays(rel_pos,
                                                 rel_vel,
                                                 total_masses,
                                                 G=constants.G)
        for ei in range(len(e)):
            if e[ei]<1:
                if si not in primaries:
                    si.planets = Particles()
                    primaries.add_particle(si)
                si.planets.add_particle(jumbos[ei])
    """
        e, a, true_anomaly,\
            inc, long_asc_node, arg_per_mat, j =\
        (list(t) for t in zip(*sorted(zip(e, a,
                                          true_anomaly,
                                          inc,
                                          long_asc_node,
                                          arg_per_mat,
                                          jumbos))))
        #print(a, e)
        if e[0]<1 and e[1]<1:
            sprimary.add_particle(si)
            jsecondary.add_particle(j[0])
            print(f"Jumbo Orbit: ae=({a[0].in_(units.au)}, {e[0]}), ({a[1].in_(units.au)}, {e[1]})")
        for ei in range(len(e)):
            if e[ei]<1:
                sma.append(a[ei])
                ecc.append(e[ei])
    
    print(f"N(JJ)= {len(sma)}")
    #print(sma)
    #print(ecc)
    plt.scatter(sma.value_in(units.au), ecc, c='b')
    plt.title("Planets in orbits around stars.")
    plt.semilogx()
    plt.show()
    """
    return primaries
    
def plot_snapshot(bodies):

    jumbos = bodies[bodies.name=="JuMBOs"]
    hosts = bodies[bodies.name=="host"]
    stars = bodies - hosts - jumbos

    bound_jumbos = find_jumbo_orbits(bodies.copy())
    triples = find_host_stellar_companion(bound_jumbos, stars+hosts)
    primaries = print_planetary_orbits(bodies.copy())
    
    m = 100*stars.mass/stars.mass.max()
    pyplot.scatter(stars.x.value_in(units.parsec), stars.y.value_in(units.parsec),
                   s=m, label="single stars", c='k', alpha= 0.2)
    if len(hosts)>0:
        m = 100*hosts.mass/hosts.mass.max()
        pyplot.scatter(hosts.x.value_in(units.parsec),
                       hosts.y.value_in(units.parsec),
                       s=m, c='b', label="host stars")
    if len(jumbos)>0:
        m = 10*jumbos.mass/jumbos.mass.max()
        pyplot.scatter(jumbos.x.value_in(units.parsec), jumbos.y.value_in(units.parsec),
                       s=m, c='r', label="jumbos")
    pyplot.xlabel('X [pc]')
    pyplot.ylabel('Y [pc]')
    pyplot.legend()
  
    pyplot.show()

def find_nearest_neighbor_planets(planets):

    nnr = [] | units.pc
    nnp = Particles()
    for ji in range(len(planets)-1):
        r = np.zeros(len(planets)) | units.pc
        for bi in range(len(planets)):
            if bi is not ji:
                r[bi] = (planets[bi].position-planets[ji].position).length()
        print("n=", len(r), len(planets))
        r, b = (list(t) for t in zip(*sorted(zip(r, planets))))
        print(r[0].in_(units.au))
        nnr.append(r[0])
        nnp.add_particle(b[0])

    return nnr, nnp

def find_nearest_star(planets, stars):
    nearest_star = Particles()
    for pi in planets:
        r = np.zeros(len(stars)) | units.pc
        for si in range(len(stars)):
            r[si] = (stars[si].position-pi.position).length()
        print("n=", si)
        r, s = (list(t) for t in zip(*sorted(zip(r, stars))))
        print(r[0].in_(units.au))
        nearest_star.add_particle(s[0])
    return nearest_star

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-n", 
                      dest="n", type="int", default = -1,
                      help="number of planets [%default]")
    result.add_option("-f", 
                      dest="filename", default = "jumbos_i0000.amuse",
                      help="input filename [%default]")
    return result
    
if __name__=="__main__":
    o, arguments  = new_option_parser().parse_args()

    bodies = read_set_from_file(o.filename)

    planets = bodies[bodies.type=="planet"]
    stars = bodies-planets
    stars.planets = Particles()
    host_stars = bodies[bodies.name=="host"]
    other_stars = bodies - host_stars - planets
    print(f"number of stars={len(other_stars)}, hosts={len(host_stars)}, jumbos={len(planets)}")

    if o.n>0:
        planets = planets[:o.n]
        bodies = stars + planets

    primaries = print_planetary_orbits(bodies.copy())
    if len(primaries)>0:
        total_masses = [] | units.MSun
        rel_pos = [] | units.pc
        rel_vel = [] | units.kms
        for pi in primaries:
            for pl in pi.planets:
                total_masses.append(pi.mass + pl.mass)
                rel_pos.append(pl.position-pi.position)
                rel_vel.append(pl.velocity-pi.velocity)
        sma, ecc, true_anomaly,\
            inc, long_asc_node, arg_per_mat =\
                get_orbital_elements_from_arrays(rel_pos,
                                                 rel_vel,
                                                 total_masses,
                                                 G=constants.G)
    
        print(f"N(JJ)= {len(sma)}")
        #print(sma)
        #print(ecc)
        plt.scatter(sma.value_in(units.au), ecc, c='b')
        plt.title(f"Planets orbiting stars: N(S)={len(primaries)}, N(pl)={len(sma)}")
        plt.semilogx()
        plt.show()
    
    bound_jumbos, single_planets = find_binary_planets(bodies.copy())
    write_set_to_file(bound_jumbos, "bound_jumbos.amuse", "amuse", close_file=True, version=2, overwrite_file=True)
    total_masses = [] | units.MSun
    rel_pos = [] | units.au
    rel_vel = [] | units.kms
    for ji in bound_jumbos:
        total_masses.append(ji.mass)
        rel_pos.append(ji.child1[0].position-ji.child2[0].position)
        rel_vel.append(ji.child1[0].velocity-ji.child2[0].velocity)
    rel_pos.reshape((-1,3))
    rel_vel.reshape((-1,3))
    sma, ecc, true_anomaly,\
        inc, long_asc_node, arg_per_mat =\
            get_orbital_elements_from_arrays(rel_pos,
                                             rel_vel,
                                             total_masses,
                                             G=constants.G)
    
    print(f"N(JJ)= {len(sma)}")
    #print(sma)
    #print(ecc)
    plt.scatter(sma.value_in(units.au), ecc, c='b')
    plt.title(f"Binary planets (Jumbo) N={len(bound_jumbos)}")
    plt.semilogx()
    plt.show()

    r = sorted(sma.value_in(units.au))
    f = np.linspace(0, 1, len(r))
    plt.plot(r, f, c='r')
    plt.show()

    
    nearest_star = find_nearest_star(bound_jumbos, stars)
    r = np.zeros(len(nearest_star)) | units.au
    for ji in range(len(nearest_star)):
        r[ji] = (nearest_star[ji].position-bound_jumbos[ji].position).length()
    r = sorted(r.value_in(units.au))
    f = np.linspace(0, 1, len(r))
    plt.plot(r, f, c='r')
    plt.title(f"Nearst star to jumbo N={len(r)}")
    plt.show()

    """        
    bodies = read_set_from_file(o.filename)
    planets = bodies[bodies.type=="planet"]
    stars = bodies-planets
    stars.planets = Particles()
    stars.type="star"
    bound_jumbos = read_set_from_file("bound_jumbos.amuse", "amuse", close_file=True)
    """
    
    if len(bound_jumbos)>0:
        triples, aj, ej, singles  = find_host_stellar_companion(bound_jumbos, stars)
        print(f"number of jumbos orbiting stars N(Sj)={len(triples)}, N(sp)={len(singles)}")
        plt.scatter(aj.value_in(units.au), ej, c='b')
    else:
        print("No binary planets orbiting stars")
    if len(single_planets)>0:
        psystems, ap, ep, singles = find_host_stellar_companion(single_planets, stars)
        print(f"number of planets orbiting stars: N(Sp)={len(psystems)}, N(ffp)={len(singles)}")
        plt.scatter(ap.value_in(units.au), ep, c='r')
    else:
        print("No single jumbos around planets")
    plt.title(f"binary planets orbiting stars N(SJ)={len(aj)}, N(SP)={len(ap)}")
    plt.semilogx()
    plt.show()
    
    
