import numpy as np
import os

from amuse.lab import units, Particles, new_salpeter_mass_distribution, constants
from amuse.ext.orbital_elements import generate_binaries
from amuse.ext.orbital_elements import orbital_elements_from_binary

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

def make_isolated_jumbos(bodies, config, data_direc):
    """Making JuMBOs"""

    data_manip = DataProcess(data_direc)

    JuMBOs = bodies[bodies.name=="JuMBOs"]
    njumbos = len(JuMBOs)
    
    JMO_minm = 0.0006 | units.MSun
    JMO_maxm = 0.013 | units.MSun
    q_min = (13.75)**-1
    q = np.sqrt(np.random.uniform(np.sqrt(q_min), 1, njumbos))
    mvals = new_salpeter_mass_distribution(njumbos, JMO_minm, JMO_maxm, alpha = -1.2)
    JuMBOs.mass = mvals
    mprim = JuMBOs.mass
    msec = JuMBOs.mass*q

    Nunder = len(msec[msec < (JMO_minm)])
    mvals = np.random.uniform(JMO_minm, JMO_maxm, Nunder) 
    msec[msec < JMO_minm] = mvals * (1 | units.MSun)
    
    sma = np.random.uniform(25, 400, njumbos) | units.au
    ecc = np.sqrt(np.random.uniform(0, np.sqrt(0.9), njumbos))
    inc = np.arccos(1-2*np.random.uniform(0, 1, njumbos)) | units.rad
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
    primaries.position += JuMBOs.position
    primaries.velocity += JuMBOs.velocity
    primaries.name = "JuMBOs"
    prvals = (primaries.mass/(1 | units.MJupiter))**(1/3)
    primaries.radius = prvals * 1 | units.RJupiter

    secondaries.position += JuMBOs.position
    secondaries.velocity += JuMBOs.velocity
    secondaries.name = "JuMBOs"
    svals = (secondaries.mass/(1 | units.MJupiter))**(1/3)
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

        lines = ["Key1: {}".format(prim_.key), 
                 "Key2: {}".format(second_.key),
                 "M1: {}".format(prim_.mass.in_(units.MSun)), 
                 "M2: {}".format(second_.mass.in_(units.MSun)),
                 "Semi-major axis: {}".format(abs(sem).in_(units.au)),
                 "Eccentricity: {}".format(ecc),
                 "Inclination: {} deg".format(inc),
                 "Argument of Periapsis: {} deg".format(arg_peri),
                 "Longitude of Asc. Node: {} deg".format(asc_node),
                 "True Anomaly: {} deg".format(true_anm),
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