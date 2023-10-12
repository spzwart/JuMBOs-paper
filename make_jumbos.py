import numpy
from numpy import random
from amuse.lab import *
from amuse.ext.orbital_elements import generate_binaries
from amuse.ic import make_planets_oligarch

def new_rotation_matrix_from_euler_angles(phi, theta, chi):
    cosp=numpy.cos(phi)
    sinp=numpy.sin(phi)
    cost=numpy.cos(theta)
    sint=numpy.sin(theta)
    cosc=numpy.cos(chi)
    sinc=numpy.sin(chi)
    #see wikipedia: http://en.wikipedia.org/wiki/Rotation_matrix
    return numpy.array(
        [[cost*cosc, -cosp*sinc + sinp*sint*cosc, sinp*sinc + cosp*sint*cosc], 
         [cost*sinc, cosp*cosc + sinp*sint*sinc, -sinp*cosc + cosp*sint*sinc],
         [-sint,  sinp*cost,  cosp*cost]])

def rotate(position, velocity, phi, theta, psi): # theta and phi in radians
    Runit = position.unit
    Vunit = velocity.unit
    matrix = new_rotation_matrix_from_euler_angles(phi, theta, psi)
    return (numpy.dot(matrix, position.value_in(Runit)) | Runit,
           numpy.dot(matrix, velocity.value_in(Vunit)) | Vunit)


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

        phi = numpy.radians(random.uniform(0.0, 90.0, 1)[0])#rotate under x
        theta0 = numpy.radians((random.normal(-90.0,90.0,1)[0]))#rotate under y
        theta0 = 0
        theta_inclination = numpy.radians(random.normal(0, 1.0, 2)) 
        theta_inclination[0] = 0
        theta = theta0 + theta_inclination
        psi = numpy.radians(random.uniform(0.0, 180.0, 1))[0]
        
        print("J=",
              outer_planets.mass.in_(units.MJupiter),
              outer_planets[1].position.length().in_(units.au)-outer_planets[0].position.length().in_(units.au))
        for pi in range(len(outer_planets)):
            outer_planets[pi].name = "J1"
            outer_planets[pi].type = "planet"
            pos = outer_planets[pi].position
            vel = outer_planets[pi].velocity
            pos,vel = rotate(pos, vel, 0, 0, psi) # theta and phi in radians
            pos,vel = rotate(pos, vel, 0, theta[pi], 0)#theta and phi in radians
            pos,vel = rotate(pos, vel, phi, 0, 0) # theta and phi in radians
            outer_planets[pi].position += si.position
            outer_planets[pi].velocity += si.velocity
            outer_planets[pi].radius = 1 | units.RJupiter

        jumbos.add_particles(outer_planets)
    return jumbos

def make_isolated_jumbos(bodies):

    JuMBOs = bodies[bodies.name=="JuMBOs"]
    njumbos = len(JuMBOs)
    print(f"N= {njumbos}")
    q = numpy.random.uniform(0.5, 1, njumbos)
    mprim = JuMBOs.mass*q
    msec = JuMBOs.mass*(1-q)
    sma = numpy.random.uniform(10, 1000, njumbos) | units.au
    #sma = 10**numpy.random.uniform(1, 4, njumbos) | units.au
    #print(sorted(sma.in_(units.au)))

    ecc = numpy.sqrt(numpy.random.uniform(0, numpy.sqrt(0.9), njumbos))
    inc = numpy.arccos(1-2*numpy.random.uniform(0,1, njumbos))| units.rad
    loan = numpy.random.uniform(0, 2*numpy.pi, njumbos)| units.rad
    aop = numpy.random.uniform(0, 2*numpy.pi, njumbos)| units.rad
    true_anomaly = numpy.random.uniform(0, 2*numpy.pi, njumbos)

    #print(mprim.in_(units.MJupiter), msec.in_(units.MJupiter), sma.in_(units.au),
    #      ecc, inc, loan, aop, true_anomaly)
    primaries, secondaries = generate_binaries(
        primary_mass=mprim,
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
    primaries.radius = 1 | units.RJupiter
    secondaries.position += JuMBOs.position
    secondaries.velocity += JuMBOs.velocity
    secondaries.name = "JuMBOs"
    secondaries.radius = 1 | units.RJupiter

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
    write_set_to_file(stars, o.outfile, "amuse", close_file=True)
