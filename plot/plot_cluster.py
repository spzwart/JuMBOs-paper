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

from amuse.ext.orbital_elements import get_orbital_elements_from_arrays
                                                 
import random

def check_paired_planets(planets):

    pairs = Particles()
    for ji in planets:
        d = (planets.position-ji.position).lengths()
        d = np.sort(d)[1]
        if d<(1000| units.au):
            ji.d = d
            pairs.add_particle(ji)
    return pairs

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

    stars.position-=jumbo.position
    stars.velocity-=jumbo.velocity
    r = [] | units.pc
    for si in stars:
        r.append(si.position.length())
    r, s = (list(t) for t in zip(*sorted(zip(r, stars))))
    nn = s[0].as_set().get_intersecting_subset_in(stars)
    a, e, = orbital_elements_of_binary(si, nn)
    print(f"Stellar orbit: {a.in_(units.au)} e={e}")
    return a, e

def find_host_stellar_companion(jumbos, stars):

    triples = Particles()
    sma = [] | units.au
    ecc = []
    for ji in jumbos:
        if hasattr(ji, "child1"):
            r = [] | units.pc
            for si in stars:
                r.append((si.position-ji.position).length())
            r, s = (list(t) for t in zip(*sorted(zip(r, stars))))
            nn = s[0].as_set().get_intersecting_subset_in(stars)
            a, e, = orbital_elements_of_binary(ji, nn)
            print(f"Jumbo Orbit around star: a={a.in_(units.au)}, {e}")
            if e<1:
                triples.add_particles(nn)
                sma.append(a)
                ecc.append(e)
    print(f"N(JJ)= {len(sma)}")
    print(sma)
    print(ecc)
    plt.title("triple parameters: (s,j)")
    plt.scatter(sma.value_in(units.au), ecc, c='b')
    plt.show()
    return triples

def find_jumbo_orbits(bodies):

    #jumbos = bodies[bodies.type=="planet"]
    jumbos = bodies[bodies.name=="jumbos"]
    stars = bodies - jumbos
    all_jumbos = Particles()
    #find nearest jumbo
    sma = [] | units.au
    ecc = []
    sma_s = [] | units.au
    ecc_s = []
    for ji in jumbos:
        b = (jumbos-ji).copy()
        b.position-=ji.position
        b.velocity-=ji.velocity
        r = [] | units.pc
        for bi in b:
            r.append(bi.position.length())
        r, b = (list(t) for t in zip(*sorted(zip(r, b))))
        nn = b[0].as_set().get_intersecting_subset_in(bodies)
        a, e, = orbital_elements_of_binary(ji, nn)
        print(f"Jumbo Orbit: a={a.in_(units.au)}, {e}")
        if e<1:
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
            all_jumbos.add_particle(jbb)
            sma.append(a)
            ecc.append(e)
        else:
            #print("companion might be a star.")
            a, e = find_stellar_companion_to_jumbos(ji, stars)
            if e<1:
                sma_s.append(a)
                ecc_s.append(e)
    print(f"N(JJ)= {len(sma)}  N(SJ)={len(sma_s)}")
    print(sma)
    print(ecc)
    plt.scatter(sma.value_in(units.au), ecc, c='b')
    plt.scatter(sma_s.value_in(units.au), ecc_s, c='r')
    plt.show()
    return all_jumbos

def print_planetary_orbits(bodies):

    jumbos = bodies[bodies.name=="jumbos"]
    stars = bodies[bodies.name=="host"]
    sma = [] | units.au
    ecc = []
    for si in stars:
        primary_mass = np.zeros(len(jumbos)) | units.MSun
        primary_mass += si.mass
        secondary_mass = jumbos.mass

        total_masses = jumbos.mass + si.mass
        rel_pos = jumbos.position-si.position
        rel_vel = jumbos.velocity-si.velocity
        a, e, true_anomaly,\
            inc, long_asc_node, arg_per_mat =\
                get_orbital_elements_from_arrays(rel_pos,
                                                 rel_vel,
                                                 total_masses,
                                                 G=constants.G)

        e, a, true_anomaly,\
            inc, long_asc_node, arg_per_mat, j =\
        (list(t) for t in zip(*sorted(zip(e, a,
                                          true_anomaly,
                                          inc,
                                          long_asc_node,
                                          arg_per_mat,                                          jumbos))))
        print(a, e)
        print(f"Jumbo Orbit: ae=({a[0].in_(units.au)}, {e[0]}), ({a[1].in_(units.au)}, {e[1]})")
        for ei in range(len(e)):
            if e[ei]<1:
                sma.append(a[ei])
                ecc.append(e[ei])
    print(f"N(JJ)= {len(sma)}")
    print(sma)
    print(ecc)
    plt.scatter(sma.value_in(units.au), ecc, c='b')
    plt.show()

def plot_snapshot(bodies, verbose, time=0| units.Myr, figname=None):
    
    fig = pyplot.figure(figsize=(8,8))
    pyplot.rcParams.update({'font.size': 22})
    mplanet = 50 | units.MJupiter
    planets = bodies[bodies.mass<mplanet]
    stars = bodies - planets

    if verbose==3:
        bound_jumbos = find_jumbo_orbits(bodies.copy())
        triples = find_host_stellar_companion(bound_jumbos, stars+hosts)
        print_planetary_orbits(bodies.copy())
    
    m = 100*np.sqrt(stars.mass/stars.mass.max())
    c = stars.mass.value_in(units.MSun)**2
    pyplot.scatter(stars.x.value_in(units.parsec), stars.y.value_in(units.parsec),
                   s=m, label="stars", c=c, cmap="bone", alpha= 0.3)
    jupiters = planets[planets.mass>=1|units.MJupiter]
    m = 10*np.sqrt(jupiters.mass/jupiters.mass.max())
    earthers = planets - jupiters
    pyplot.scatter(jupiters.x.value_in(units.parsec), jupiters.y.value_in(units.parsec),
                   s=m, c='r', label="Jupiter-like planets")
    pyplot.scatter(earthers.x.value_in(units.parsec), earthers.y.value_in(units.parsec),
                   s=0.1, c='k', label="Earth-like planets")

    if verbose==2:    
        pairs = check_paired_planets(jupiters)
        print(f"number of stars={len(stars)}, jumbos={len(pairs)}, jupiters={len(jupiters)}, Earthers={len(earthers)}")

        if len(pairs)>0:
            pyplot.scatter(pairs.x.value_in(units.parsec), pairs.y.value_in(units.parsec),
                           s=30, c=pairs.d.value_in(units.au), cmap="turbo",
                           label="jumbos")
        pyplot.title(f"t={time.value_in(units.Myr):2.3f}Myr, N={len(pairs)}")
    else:
        pyplot.title(f"t={time.value_in(units.Myr):2.3f}Myr")
    
    pyplot.xlabel('X [pc]')
    pyplot.ylabel('Y [pc]')
    #pyplot.legend()
    if not figname == None:
        pyplot.savefig(figname)
    pyplot.show()
    pyplot.close()

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-f", 
                      dest="filename", default = "jumbos_i0000.amuse",
                      help="input filename [%default]")
    result.add_option("--verbose", 
                      dest="verbose", type="int", default = 1,
                      help="plotted verbosity [%default]")
    return result
    
if __name__=="__main__":
    o, arguments  = new_option_parser().parse_args()

    stars = read_set_from_file(o.filename)
    print(stars)
    plot_snapshot(stars, o.verbose)

