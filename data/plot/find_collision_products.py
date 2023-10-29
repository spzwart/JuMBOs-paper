import numpy as np
from amuse.units import units
from amuse.plot import plot
from matplotlib import pyplot as plt
from amuse.lab import *

import pandas as pd

def find_collision_products(b1, b2):

    absent = Particles()
    for bi in range(len(b1)):
        if b1[bi].key not in b2.key:
            absent.add_particle(b1[bi])
    return absent

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("--f1", dest="f1", default = "jumbos_i0000.amuse",
                      help="input filename[%default]")
    result.add_option("--f2", dest="f2", default = "jumbos_i0010.amuse",
                      help="input filename[%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    #type,key,k1,k2,M,m,a,e,i,x,y,z,vx,vy,vz
    b1 = read_set_from_file(o.f1)
    b2 = read_set_from_file(o.f2)

    same = b1.get_intersecting_subset_in(b2)
    print(f"unaffectd N={len(same)}")
    lost = b1-same
    gain = b2-same
    #print(len(same), len(lost), len(gain))
    #print(lost.mass.in_(units.MJupiter))
    #print(gain.mass.in_(units.MJupiter))

    jj_collision = lost[lost.mass<20|units.MJupiter]
    ss_collision = lost - jj_collision
    print(f"Number of stars in collisions = {len(ss_collision)}")
    print(f"Number of jupiters in collisions = {len(jj_collision)}")
    
