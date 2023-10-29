import numpy
from numpy import random
from amuse.lab import *
from amuse.ext.orbital_elements import generate_binaries, new_binary_from_orbital_elements
from amuse.ic import make_planets_oligarch

def Hill_radius(Mstar, a, Mplanet):
    return a * (2*Mplanet/(3.0*Mstar))**(1./3.)

def sma_from_Hill_radius(Mstar, Mplanet, rH):
    a = rH/((Mplanet/(3.0*Mstar))**(1./3.))
    return a

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
    print("Jumbos as Oligarchic outer circum-stellar planets.")

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

def make_arXiv2310_06016(bodies, a1=800|units.au, a2=1000|units.au, jumbo_mass_function=False):
    print("Jumbos as outer circum-stellar planets following ArXiv 2310_06016.")
    
    host_stars = bodies[bodies.name=="host"]
    nhost_stars = len(host_stars)
    print(f"make ArXiv Jumbos N={nhost_stars}")

    jumbos = Particles()
    for si in host_stars:
        if jumbo_mass_function:
            msec = [0, 0] | units.MJupiter
            while msec[0]<0.6|units.MJupiter or msec[1]<0.6|units.MJupiter:
                Mjumbos = new_salpeter_mass_distribution(1,
                                                         0.8|units.MJupiter,
                                                         14|units.MJupiter, alpha=-1.2)[0]
                #q = numpy.sqrt(numpy.random.uniform(0.5**2, 1, 1))        
                #q = numpy.sqrt(numpy.random.uniform(0, 1, 1))[0]
                q = numpy.sqrt(numpy.random.uniform(0.2**2, 1, 1)[0])
                mprim = [si.mass.value_in(units.MSun), si.mass.value_in(units.MSun)] | units.MSun
                msec = [q, 1-q] * Mjumbos
        else:
            Mjumbos = 1.e-3 * si.mass
            mprim = [si.mass.value_in(units.MSun), si.mass.value_in(units.MSun)] | units.MSun
            msec = 1.e-3 * mprim
        #sma = [a1.value_in(units.au), a2.value_in(units.au)] | units.au
        #a1 = 10**numpy.random.uniform(1, 3, 1)[0] | units.au
        #a1 = numpy.random.uniform(30, 3000, 1)[0] | units.au
        a1 = numpy.random.uniform(20, 2000, 1)[0] | units.au
        #a1 = 1000|units.au
        rH = Hill_radius(mprim[0], a1, msec[0])
        #print("Hill radius=", a1.in_(units.au), rH.in_(units.au))
        #a2 = sma_from_Hill_radius(mprim[1], msec[1], rH)
        #a2 = a1 + 10*rH
        a2 = a1 + 5*rH
        #a2 = a1 + 2*rH
        sma = [a1.value_in(units.au), a2.value_in(units.au)] | units.au
        #print("initial semimajor_axis=", sma.in_(units.au))
        #print("initial planet masses=", msec.in_(units.MJupiter))
        ecc = numpy.sqrt(numpy.random.uniform(0, 0.02**2, 2))
        inc = numpy.arccos(1-2*numpy.random.uniform(0,1, 2)) | units.rad
        inc[1] = inc[0] + numpy.radians(numpy.random.uniform(-1,1, 1)) | units.rad
        loan = numpy.random.uniform(0, 2*numpy.pi, 2) | units.rad
        loan[1] = loan[0] + numpy.radians(numpy.random.uniform(-1, 1, 1)) | units.rad
        aop = numpy.random.uniform(0, 2*numpy.pi, 2)| units.rad
        aop[1] = aop[0] + numpy.radians(numpy.random.uniform(-1, 1, 1)) | units.rad
        ta = numpy.random.uniform(0, 2*numpy.pi, 2)| units.rad

        for bi in range(len(sma)):
            binary = new_binary_from_orbital_elements(
                mprim[bi], msec[bi], 
                semimajor_axis=sma[bi],
                eccentricity=ecc[bi],
                true_anomaly=ta[bi],
                inclination=inc[bi] ,
                longitude_of_the_ascending_node=loan[bi] ,
                argument_of_periapsis=aop[bi] ,
                G=constants.G)
            star = binary[0]
            planet = binary[1]
            planet.position -= star.position
            planet.velocity -= star.velocity
            planet.position += si.position
            planet.velocity += si.velocity
            planet.name = "jumbos"
            planet.type = "planet"
            planet.radius = 1 | units.RJupiter
            jumbos.add_particle(planet)

    return jumbos

def make_jumbo_as_planetmoon_pair(bodies): 
    print("Jumbos as outer circum-stellar planets-moon pair.")

    host_stars = bodies[bodies.name=="host"]
    nhost_stars = len(host_stars)
    print(f"make ArXiv Jumbos N={nhost_stars}")

    all_jumbos = Particles()

    for si in host_stars:
        #q = numpy.sqrt(numpy.random.uniform(0.5**2, 1, 1))
        #q = numpy.random.uniform(0.5, 1, 1)
        q = numpy.sqrt(numpy.random.uniform(0.2**2, 1, 1)[0])
        Mjumbos = new_salpeter_mass_distribution(1,
                                                 1|units.MJupiter,
                                                 20|units.MJupiter, alpha=-1.2)
        mprim = Mjumbos*q
        msec = Mjumbos*(1-q) 
        sma = numpy.random.uniform(10, 400, 1) | units.au
        ecc = numpy.sqrt(numpy.random.uniform(0, 0.02**2, 1))
        inc = numpy.arccos(1-2*numpy.random.uniform(0,1, 1)) | units.rad
        loan = numpy.random.uniform(0, 2*numpy.pi, 1) | units.rad
        aop = numpy.random.uniform(0, 2*numpy.pi, 1)| units.rad
        ta = numpy.random.uniform(0, 2*numpy.pi, 1)| units.rad

        jumbo = new_binary_from_orbital_elements(
            mprim, msec, 
            semimajor_axis=sma,
            eccentricity=ecc,
            true_anomaly=ta,
            inclination=inc ,
            longitude_of_the_ascending_node=loan ,
            argument_of_periapsis=aop ,
            G=constants.G)
        jumbo.name = "jumbos"
        jumbo.type = "planet"
        jumbo.radius = 1 | units.RJupiter

        mprim = si.mass
        msec = jumbo.mass.sum()

        rH = 3*sma
        sma = sma_from_Hill_radius(si.mass, jumbo.mass.sum(), rH)
        ecc = numpy.sqrt(numpy.random.uniform(0, 0.02**2, 1))
        inc = numpy.arccos(1-2*numpy.random.uniform(0, 1, 1)) | units.rad
        loan = numpy.random.uniform(0, 2*numpy.pi, 1) | units.rad
        aop = numpy.random.uniform(0, 2*numpy.pi, 1)| units.rad
        ta = numpy.random.uniform(0, 2*numpy.pi, 1)| units.rad
        
        binary = new_binary_from_orbital_elements(
            mprim, msec, 
            semimajor_axis=sma,
            eccentricity=ecc,
            true_anomaly=ta,
            inclination=inc ,
            longitude_of_the_ascending_node=loan ,
            argument_of_periapsis=aop ,
            G=constants.G)
        binary[1].position -= binary[0].position
        binary[1].velocity -= binary[0].velocity
        jumbo.position += binary[1].position
        jumbo.velocity += binary[1].velocity
        jumbo.position += si.position
        jumbo.velocity += si.velocity
        all_jumbos.add_particles(jumbo)            
    return all_jumbos

def make_isolated_jumbos(bodies):
    print("Jumbos as isolated free floating binary pair.")

    JuMBOs = bodies[bodies.name=="jumbos"]
    njumbos = len(JuMBOs)
    print(f"N= {njumbos}")
    q = numpy.sqrt(numpy.random.uniform(0.5**2, 1, njumbos))
    mprim = JuMBOs.mass*q
    msec = JuMBOs.mass*(1-q)
    sma = numpy.random.uniform(10, 1000, njumbos) | units.au
    #sma = 10**numpy.random.uniform(1, 4, njumbos) | units.au
    #print(sorted(sma.in_(units.au)))

    ecc = numpy.sqrt(numpy.random.uniform(0, 0.9**2, njumbos))
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
    primaries.name = "jumbos"
    primaries.type = "planet"
    primaries.radius = 1 | units.RJupiter
    secondaries.position += JuMBOs.position
    secondaries.velocity += JuMBOs.velocity
    secondaries.name = "jumbos"
    secondaries.type = "planet"
    secondaries.radius = 1 | units.RJupiter

    jumbos = Particles()
    jumbos.add_particles(primaries)
    jumbos.add_particles(secondaries)
    return jumbos


def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("--Njumbos", dest="Njumbos", type="int",default = 5,
                      help="number of JuMBOs [%default]")
    result.add_option("--model", dest="jumbo_model", type="int",default = "classic",
                      help="select jumbo model (freefloaters, circum_stellar, planetmoon, oligarchic) [%default]")
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

    masses = new_kroupa_mass_distribution(N, mass_min=Mmin, mass_max=Mmax)
    Rvir = 1|units.pc
    converter=nbody_system.nbody_to_si(masses.sum(), Rvir)
    stars = new_plummer_model(10, convert_nbody=converter)
    bodies.name = "star"
    
    if Njumbos>0:
        jumbos = make_isolated_jumbos(stars, o.Njumbos)
    else:
        jumbos = bodies.random_sample(o.Njumbos)
        q = numpy.sqrt(numpy.random.uniform(0.5**2, 1, o.Njumbos))
        #q = numpy.random.uniform(0.5, 1, o.Njumbos)
        jumbo.mass = new_salpeter_mass_distribution(o.Njumbos,
                                                 1|units.MJupiter,
                                                 20|units.MJupiter, alpha=-1.2)
        jumbos.name = "JuMBOs"
        jumbos.radius = 1 | units.RJupiter
        if o.arxive==0:        
            jumbos = make_outer_planetary_systems(stars)
        else:
            jumbos = make_arXiv2310_06016(stars)
    
    stars.add_particles(jumbos)
    write_set_to_file(stars, o.outfile, "amuse", close_file=True)
