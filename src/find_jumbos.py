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

from process_jumbos import find_companions
from process_jumbos import is_this_binary_planet_orbiting_a_star


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
    for bi in bodies:
        bi.companions = Particles()
        bi.parent = Particles()
    #bodies.add_vector_attribute("abundance", ('H', 'He', 'Li', 'Fe'))
    planets = bodies[bodies.type=="planet"]
    stars = bodies-planets
    print(f"Nbodies={len(bodies)}, {len(stars)}, {len(planets)}")

    #find companion of stars to individual stars
    jumbos = find_companions(planets, planets)
    mmax = stars.mass.max()
    for ji in jumbos:
        if len(bi.companions)>0:
            print(f"Jumbo M= {ji.mass.in_(units.MSun)}, a= {ji.companions.semimajor_axis.in_(units.au)}, e= {ji.companions.eccentricity}")
        else:
            a, e, = is_this_binary_planet_orbiting_a_star(stars, ji)
            if len(e)>0:
                print(f"Stellar companions to planet: a= {a.in_(units.au)}, e={e}")
