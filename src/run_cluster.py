import sys
import numpy
from numpy import random
from amuse.lab import *
from amuse.io import store
from amuse.community.seba.interface import SeBa
from amuse.ext.LagrangianRadii import LagrangianRadii
from amuse.community.petar.interface import Petar
from amuse.community.fractalcluster.interface import new_fractal_cluster_model

from matplotlib import pyplot as plt
from make_jumbos import make_outer_planetary_systems
from make_jumbos import make_isolated_jumbos
from make_jumbos import make_planetplanet
from make_jumbos import make_jumbo_as_planetmoon_pair
from make_jumbos import make_singletons

def ZAMS_radius(mass):
    log_mass = numpy.log10(mass.value_in(units.MSun))
    mass_sq = (mass.value_in(units.MSun))**2
    alpha = 0.08353 + 0.0565*log_mass
    beta  = 0.01291 + 0.2226*log_mass
    gamma = 0.1151 + 0.06267*log_mass
    r_zams = pow(mass.value_in(units.MSun), 1.25) * (0.1148 + 0.8604*mass_sq) / (0.04651 + mass_sq)
    return r_zams | units.RSun

def merge_two_stars(bodies, particles_in_encounter):
    com_pos = particles_in_encounter.center_of_mass()
    com_vel = particles_in_encounter.center_of_mass_velocity()
    new_particle=Particles(1)
    new_particle.mass = particles_in_encounter.total_mass()
    new_particle.position = com_pos
    new_particle.velocity = com_vel
    new_particle.radius = particles_in_encounter.radius.sum()
    if "jumbos" in particles_in_encounter.name:
        print("Merge with jumbo: ", particles_in_encounter.name, "M=", particles_in_encounter.mass.in_(units.MJupiter))
    else:
        print("Merge with stars: ", particles_in_encounter.name, "M=", particles_in_encounter.mass.in_(units.MJupiter))

    if new_particle.mass>=0.08|units.MSun:
        new_particle.name = "star"
        new_particle.radius = ZAMS_radius(new_particle.mass)
    else:
        new_particle.name = "jumbos"
        rho = (1|units.MJupiter)/(1|units.RJupiter)**3 
        new_particle.radius = (new_particle.mass/rho)**(1./3.)
    bodies.add_particles(new_particle)
    bodies.remove_particles(particles_in_encounter)

def resolve_collision(collision_detection, gravity, bodies):
    if collision_detection.is_set():
        E_coll = gravity.kinetic_energy + gravity.potential_energy
        print("At time=", gravity.model_time.in_(units.Myr), \
              "number of encounters=", len(collision_detection.particles(0)))
        Nenc = 0
        for ci in range(len(collision_detection.particles(0))):
            print("Collision between: ", collision_detection.particles(0)[ci].key, "and", collision_detection.particles(1)[ci].key)
            particles_in_encounter \
                = Particles(particles=[collision_detection.particles(0)[ci],
                                       collision_detection.particles(1)[ci]])
            particles_in_encounter \
                = particles_in_encounter.get_intersecting_subset_in(bodies)

            merge_two_stars(bodies, particles_in_encounter)
            bodies.synchronize_to(gravity.particles)
            Nenc += 1
            print("Resolve encounter Number:", Nenc)
        dE_coll = E_coll - (gravity.kinetic_energy + gravity.potential_energy)
        print("dE_coll =", dE_coll, "N_enc=", Nenc)
    sys.stdout.flush()

def resolve_supernova(supernova_detection, bodies, time):
    if supernova_detection.is_set():
        print("At time=", time.in_(units.Myr), \
              len(supernova_detection.particles(0)), 'supernova(e) detected')

        Nsn = 0
        for ci in range(len(supernova_detection.particles(0))):
            print(supernova_detection.particles(0))
            particles_in_supernova \
                = Particles(particles=supernova_detection.particles(0))
            natal_kick_x = particles_in_supernova.natal_kick_x
            natal_kick_y = particles_in_supernova.natal_kick_y
            natal_kick_z = particles_in_supernova.natal_kick_z

            particles_in_supernova \
                = particles_in_supernova.get_intersecting_subset_in(bodies)
            particles_in_supernova.vx += natal_kick_x
            particles_in_supernova.vy += natal_kick_y
            particles_in_supernova.vz += natal_kick_z
            Nsn += 1

        print('Resolved', Nsn, 'supernova(e)')

def make_initial_cluster(Nstars, Njumbos, Rvir, Fd, jumbo_model,
                         a1=800|units.au, a2=1000|units.au,
                         x=-2.0, mmin=0.3|units.MJupiter,
                         jumbo_mass_function=True):

    if jumbo_model=="singletons":
        return make_singletons(Nstars, Njumbos, Rvir, Fd, x, mmin)

    N = Nstars
    if jumbo_model=="freefloaters":
        N = Nstars + int(Njumbos)

    Mmin = 0.08 | units.MSun
    Mmax = 100 | units.MSun
    masses = new_kroupa_mass_distribution(N, mass_min=Mmin, mass_max=Mmax)
    
    Mtot_init = masses.sum()
    converter=nbody_system.nbody_to_si(Mtot_init, Rvir)#1|units.Myr)
    if Fd>0:
        bodies = new_fractal_cluster_model(N, fractal_dimension=1.6, convert_nbody=converter)
    else:
        bodies = new_plummer_model(N, convert_nbody=converter)
    bodies.mass = masses
    bodies.name = "star"
    bodies.type = "star"
    bodies.radius = ZAMS_radius(bodies.mass)
    if jumbo_model=="freefloaters":
        JuMBOs = bodies.random_sample(Njumbos)
        JuMBOs.mass = new_salpeter_mass_distribution(Njumbos,
                                                     0.8|units.MJupiter,
                                                     14|units.MJupiter, alpha=x)
        q = numpy.sqrt(numpy.random.uniform(0.2**2, 1, Njumbos))
        JuMBOs.name = "jumbos"
        JuMBOs.type = "planet"
        JuMBOs.radius = 1 | units.RJupiter
        JuMBOs.m2 = JuMBOs.mass * (1-q)
        bodies.scale_to_standard(convert_nbody=converter)
        bodies.move_to_center()
        jumbos = make_isolated_jumbos(bodies, a1, a2)
        bodies.remove_particles(JuMBOs)
    else:
        bodies.scale_to_standard(convert_nbody=converter)
        bodies.move_to_center()
        bodies = bodies.sorted_by_attribute('mass')
        #print(bodies.mass.in_(units.MSun))
        for bi in range(len(bodies)):
            if bodies[bi].mass>0.6|units.MSun:
                break
        print(bodies[bi].mass.in_(units.MSun))
        print(bodies[bi-int(Njumbos/2)].mass.in_(units.MSun))
        print(bodies[bi+int(Njumbos/2)].mass.in_(units.MSun))
        host_stars = bodies[bi-int(Njumbos/2):bi+int(Njumbos/2)]
        host_stars.name = "host"
        nhost_stars = len(host_stars)
        print(f"Mass limit for jumbos: {host_stars.mass.min().in_(units.MSun)}, {host_stars.mass.max().in_(units.MSun)}")
        if jumbo_model=="circum_stellar":
            jumbos = make_planetplanet(bodies, a1, a2, jumbo_mass_function)
        elif jumbo_model=="planetmoon":
            jumbos = make_jumbo_as_planetmoon_pair(bodies, a1, a2, x)
        elif jumbo_model=="oligarchic":
            jumbos = make_outer_planetary_systems(bodies)
        else:
            print(f"No Jumbo model selected: {jumbo_model}")
    #print(jumbos)
    bodies.add_particles(jumbos)
    #from plot_cluster import print_planetary_orbits
    #print_planetary_orbits(bodies.copy())
    
    return bodies
        
def  run_cluster(bodies, Rvir, t_end, dt):

    stars = bodies[bodies.name=="star"]
    hosts = bodies[bodies.name=="host"]
    jumbos = bodies-stars-hosts
    print(f"Stars N={len(stars)}, hosts N={len(hosts)}, Jumbos N={len(jumbos)}")
    print(f"Initial rage in a={o.a1.in_(units.au)}, {o.a2.in_(units.au)}") 

    """
    plt.scatter(stars.x.value_in(units.pc),  stars.y.value_in(units.pc), s=1, c='k') 
    plt.scatter(hosts.x.value_in(units.pc), hosts.y.value_in(units.pc), s=3, c='y')
    plt.scatter(jumbos.x.value_in(units.pc), jumbos.y.value_in(units.pc), s=1, c='r')
    plt.show()
    """

    converter=nbody_system.nbody_to_si(bodies.mass.sum(), Rvir)
    gravity = ph4(converter, number_of_workers=3)
    #gravity = Petar(converter)#, mode="gpu")#, number_of_workers=6)
    #gravity = Petar(converter, mode="gpu")#, number_of_workers=6)
    #print(gravity.parameters)
    #gravity.parameters.epsilon_squared = (1|units.au)**2
    #gravity.parameters.timestep_parameter = 0.05
    gravity.particles.add_particles(bodies)
    collision_detection = gravity.stopping_conditions.collision_detection
    collision_detection.enable()

    channel_from_gd = gravity.particles.new_channel_to(bodies)
    channel_to_gd = bodies.new_channel_to(gravity.particles)

    index = 0
    filename = "jumbos_i{0:04}.amuse".format(index)
    write_set_to_file(bodies, filename, 'amuse',
                          close_file=True, overwrite_file=True)
    E_init = gravity.kinetic_energy + gravity.potential_energy
    
    Nsn = 0
    dE_coll = zero
    time = zero

    dt_diag = dt
    diag_time = dt_diag
    while time < t_end:

        print(f"Time steps: {dt.in_(units.Myr)} time= {time.in_(units.Myr)}")
        time += dt

        while gravity.model_time<time:
            channel_to_gd.copy()
            E_dyn = gravity.kinetic_energy  + gravity.potential_energy 
            gravity.evolve_model(time)
            dE_dyn = E_dyn - (gravity.kinetic_energy  + gravity.potential_energy)
            resolve_collision(collision_detection, gravity, bodies)
            channel_from_gd.copy()
            
        print_diagnostics(time, bodies.mass.sum(), E_dyn, dE_dyn, dE_coll)
        
        sys.stdout.flush()
        if diag_time <= time:
            diag_time += dt_diag
            print("Diagnostics:", time.in_(units.Myr), "N=", len(bodies))
            Rv = bodies.virial_radius()
            SMBH = bodies[bodies.name=="SMBH"]
            RL = LagrangianRadii(bodies-SMBH)
            print("Time=", time.in_(units.Myr), "Rv=", Rv.in_(units.pc), "Rl=",
                  RL[-5].value_in(units.pc),
                  RL[-4].value_in(units.pc),
                  RL[-3].value_in(units.pc))
            index += 1
            filename = "jumbos_i{0:04}.amuse".format(index)
            write_set_to_file(bodies, filename, 'hdf5',
                              close_file=True, overwrite_file=False)
            #write_set_to_file(bodies.savepoint(time), filename, 'hdf5',
            #                  close_file=True, overwrite_file=False)

        
    gravity.stop()

def print_diagnostics(time, Mtot, Etot, dE_dyn, dE_coll):
    print("T=", time.in_(units.Myr), end=' ') 
    print("M=", Mtot.in_(units.MSun), end=' ') 
    print("E= ", Etot, end=' ') 
    print("dE(dyn)=", dE_dyn/Etot, end=' ') 
    print("dE(coll)=", dE_coll/Etot) 
    
def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("--Nstars", dest="Nstars", type="int",default = 2500,
                      help="number of stars [%default]")
    result.add_option("--NJuMBOs", dest="Njumbos", type="int",default = 300,
                      help="number of JuMBs [%default]")
    result.add_option("-F", dest="Fd", type="float", default = -1,
                      help="fractal dimension [%default]")
    result.add_option("--dt", unit=units.Myr,
                      dest="dt", type="float",default = 0.1|units.Myr,
                      help="output timesteps [%default]")
    result.add_option("-x", 
                      dest="x", type="float",default = -2.0,
                      help="mass function slope [%default]")
    result.add_option("--mmin", unit=units.MJupiter,
                      dest="mmin", type="float",default = 0.3|units.MJupiter,
                      help="mass function minimum [%default]")
    result.add_option("-R", unit=units.parsec,
                      dest="Rvir", type="float",default = 0.5 | units.pc,
                      help="cluser virial radius [%default.value_in(units.parsec)]")
    result.add_option("--a1", unit=units.au,
                      dest="a1", type="float",default = 800 | units.au,
                      help="inner binary orbit for circum_stellar model [%default.value_in(units.parsec)]")
    result.add_option("--a2", unit=units.au,
                      dest="a2", type="float",default = 1000 | units.au,
                      help="outer binary orbit for circum_stellar model [%default.value_in(units.parsec)]")
    result.add_option("-t", unit=units.Myr,
                      dest="t_end", type="float", default = 1.0 | units.Myr,
                      help="end time of the simulation [%default.value_in(units.Myr]")
    result.add_option("--model", dest="jumbo_model", default = "classic",
                      help="select jumbo model (freefloaters, circum_stellar, planetmoon, oligarchic, singletons) [%default]")
    result.add_option("-q", action='store_true',
                      dest="jumbo_mass_function", default="True",
                      help="jumbo mass function [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()

    bodies = make_initial_cluster(o.Nstars, o.Njumbos, o.Rvir, o.Fd, o.jumbo_model,
                                  o.a1, o.a2, o.x, o.mmin, o.jumbo_mass_function)
    run_cluster(bodies, o.Rvir, o.t_end, o.dt)

