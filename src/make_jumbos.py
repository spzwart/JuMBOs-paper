import numpy as np
import os
from numpy import random

from amuse.lab import *
from amuse.ext.orbital_elements import generate_binaries
from amuse.ext.orbital_elements import orbital_elements_from_binary
from amuse.ic import make_planets_oligarch

from src.data_process import DataProcess

def new_rotation_matrix_from_euler_angles(phi, theta, chi):
    cosp=np.cos(phi)
    sinp=np.sin(phi)
    cost=np.cos(theta)
    sint=np.sin(theta)
    cosc=np.cos(chi)
    sinc=np.sin(chi)
    #see wikipedia: http://en.wikipedia.org/wiki/Rotation_matrix
    return np.array(
        [[cost*cosc, -cosp*sinc + sinp*sint*cosc, sinp*sinc + cosp*sint*cosc], 
         [cost*sinc, cosp*cosc + sinp*sint*sinc, -sinp*cosc + cosp*sint*sinc],
         [-sint,  sinp*cost,  cosp*cost]])

def rotate(position, velocity, phi, theta, psi): # theta and phi in radians
    Runit = position.unit
    Vunit = velocity.unit
    matrix = new_rotation_matrix_from_euler_angles(phi, theta, psi)
    return (np.dot(matrix, position.value_in(Runit)) | Runit,
           np.dot(matrix, velocity.value_in(Vunit)) | Vunit)

def make_outer_planetary_systems(bodies): 

    host_stars = bodies[bodies.name=="host"]
    nhost_stars = len(host_stars)
    print(f"N= {nhost_stars}")

    jumbos = Particles()
    for si in host_stars:
        mass_star = si.mass
        radius_star = si.radius
        inner_radius_disk = 10| units.au
        outer_radius_disk = 300|units.au/mass_star.value_in(units.MSun)
        mass_disk = 0.02*mass_star
        planetary_system = make_planets_oligarch.new_system(mass_star, 
                                                            radius_star,
                                                            inner_radius_disk,
                                                            outer_radius_disk,
                                                            mass_disk)
        all_planets = planetary_system.planets[0]
        outer_planets = all_planets[-2:]

        phi = np.radians(random.uniform(0.0, 90.0, 1)[0])     #rotate under x
        theta0 = np.radians((random.normal(-90.0,90.0,1)[0])) #rotate under y
        theta0 = 0
        theta_inclination = np.radians(random.normal(0, 1.0, 2)) 
        theta_inclination[0] = 0
        theta = theta0 + theta_inclination
        psi = np.radians(random.uniform(0.0, 180.0, 1))[0]
        
        for pi in range(len(outer_planets)):
            outer_planets[pi].name = "J1"
            outer_planets[pi].type = "planet"
            pos = outer_planets[pi].position
            vel = outer_planets[pi].velocity
            pos,vel = rotate(pos, vel, 0, 0, psi)       #theta and phi in radians
            pos,vel = rotate(pos, vel, 0, theta[pi], 0) #theta and phi in radians
            pos,vel = rotate(pos, vel, phi, 0, 0)       #theta and phi in radians
            outer_planets[pi].position += si.position 
            outer_planets[pi].velocity += si.velocity
            outer_planets[pi].radius = 1 | units.RJupiter

        jumbos.add_particles(outer_planets)
    return jumbos

def make_isolated_jumbos(bodies, config, data_direc):
    """Making JuMBOs.
       Mass ratio and semi-major axis are taken from dynamical range
       observed by Pearson and McCaughrean (2023)
    """

    data_manip = DataProcess(data_direc)

    JuMBOs = bodies[bodies.name=="JuMBOs"]
    njumbos = len(JuMBOs)
    
    q = np.random.uniform(np.sqrt((13.75)**-1), 1, njumbos)
    mvals = np.random.uniform(0.0006, 0.013, njumbos)
                                           
    JuMBOs.mass = mvals * (1 | units.MSun)
    mprim = JuMBOs.mass
    msec = JuMBOs.mass*q

    sma = np.random.uniform(10, 1000, njumbos) | units.au
    ecc = np.sqrt(np.random.uniform(0, np.sqrt(0.9), njumbos))
    inc = np.arccos(1-2*np.random.uniform(0,1, njumbos)) | units.rad
    loan = np.random.uniform(0, 2*np.pi, njumbos) | units.rad
    aop = np.random.uniform(0, 2*np.pi, njumbos) | units.rad
    true_anomaly = np.random.uniform(0, 2*np.pi, njumbos)
    
    primaries, secondaries = generate_binaries(primary_mass=mprim,
                                               secondary_mass=msec,
                                               semi_major_axis=sma,
                                               eccentricity=ecc,
                                               true_anomaly=true_anomaly, 
                                               inclination=inc,
                                               longitude_of_the_ascending_node=loan,
                                               argument_of_periapsis=aop,
                                               G=constants.G)
    print(primaries.position)
    STOP
    primaries.position += JuMBOs.position
    primaries.velocity += JuMBOs.velocity
    primaries.name = "JuMBOs"
    prvals = (primaries.mass/ (1 | units.MJupiter))**(1/3)
    primaries.radius = prvals * 1 | units.RJupiter

    secondaries.position += JuMBOs.position
    secondaries.velocity += JuMBOs.velocity
    secondaries.name = "JuMBOs"
    svals = (secondaries.mass/ (1 | units.MJupiter))**(1/3)
    secondaries.radius =  svals * 1 | units.RJupiter

    for prim_, second_ in zip(primaries, secondaries):
        bin_sys = Particles()
        bin_sys.add_particle(prim_)
        bin_sys.add_particle(second_)
        kepler_elements = orbital_elements_from_binary(bin_sys, G=constants.G)

        sem = kepler_elements[2]
        ecc = kepler_elements[3]
        inc = kepler_elements[4]
        arg_peri = kepler_elements[5]
        asc_node = kepler_elements[6]
        true_anm = kepler_elements[7]

        lines = ["Key1: {}".format(prim_.key), "Key2: {}".format(second_.key),
                 "M1: {} ".format(prim_.mass.in_(units.MSun)), 
                 "M2: {} ".format(second_.mass.in_(units.MSun)),
                 "Semi-major axis: {}".format(abs(sem).in_(units.au)),
                 "Eccentricity: {}".format(ecc),
                 "Inclination: {} deg".format(inc),
                 "Argument of Periapsis: {} degrees".format(arg_peri),
                 "Longitude of Asc. Node: {} degrees".format(asc_node),
                 "True Anomaly: {} degrees".format(true_anm),
                 "================================================================="]
        
        with open(os.path.join(data_manip.path+str('initial_binaries'),
                'initial_bins_'+str(config)+'.txt'), 'a') as f:
            for line_ in lines:
                f.write(line_)
                f.write('\n')

    jumbos = Particles()
    jumbos.add_particles(primaries)
    jumbos.add_particles(secondaries)
    return jumbos

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("--Njumbos", dest="Njumbos", type="int",default = 500,
                      help="number of JuMBOs [%default]")
    result.add_option("-f", unit=units.Myr,
                      dest="filename", default = "stars.amuse",
                      help="end time of the simulation [%default.value_in(units.Myr]")
    result.add_option("--mmin", unit=units.MSun, 
                      dest="mmin", type="float", default = 0.5|units.MSun,
                      help="minimum stellar mass for planets [%default.value_in(units.Myr]")
    result.add_option("--mmax", unit=units.MSun, 
                      dest="mmax", type="float", default = 3.0|units.MSun,
                      help="maximum stellar mass for planets [%default.value_in(units.Myr]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    stars = read_set_from_file(o.filename)
    if Njumbos>0:
        jumbos = make_isolated_jumbos(stars, o.Njumbos)
    else:
        jumbos = make_outer_planetary_systems(stars)
    stars.add_particles(jumbos)
