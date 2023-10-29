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
    print(f"Nstar-binary-planet = {n_starbinaryplanet}")
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
    p = Particles(2)
    p[0].mass = bi.mass
    p[0].position = bi.position
    p[0].velocity = bi.velocity
    p[1].mass = ci.mass
    p[1].position = ci.position
    p[1].velocity = ci.velocity
    cm = Particle()
    cm.mass = p.mass.sum()
    cm.position = p.center_of_mass()
    cm.velocity = p.center_of_mass_velocity()
    return cm

def construct_center_of_mass_particle(body):
    p = Particles()
    p.add_particle(body)
    p.add_particles(body.companions)
    cm = Particle()
    cm.mass = p.mass.sum()
    cm.position = p.center_of_mass()
    cm.velocity = p.center_of_mass_velocity()
    return cm

def find_companions(stars, planets):

    for bi in range(len(stars)):
        if stars[bi] in planets:
            plts = planets - stars[bi]
        else:
            plts = planets
        total_masses = plts.mass + stars[bi].mass
        rel_pos = plts.position-stars[bi].position
        rel_vel = plts.velocity-stars[bi].velocity
        sma, ecc, true_anomaly,\
            inc, long_asc_node, arg_per_mat =\
                get_orbital_elements_from_arrays(rel_pos,
                                                 rel_vel,
                                                 total_masses,
                                                 G=constants.G)
        e = np.ma.masked_where(ecc>1, ecc)
        e = np.ma.masked_invalid(e)
        if(len(np.ma.compressed(e))>0):
            k = np.ma.masked_array(plts.key, e.mask)
            k = np.ma.compressed(k)
            p = Particles()
            for pi in plts:
                if pi.key in k:
                    p.add_particle(pi)
            #p = planets[planets.key==k] # does not seem to work properly?
            a = np.ma.masked_array(sma, e.mask)
            a = np.ma.compressed(a)
            e = np.ma.compressed(e)
            #sort on semi-major axis
            a, e, p = (list(t) for t in zip(*sorted(zip(a, e, p))))
            
            #print("companion found:", e, len(p))
            for pi in range(len(p)):
                if hasattr(p[pi], "semimajor_axis") and p[pi].semimajor_axis>0|units.au:
                    print(f"object already a binary member: a={p[pi].semimajor_axis.in_(units.au)}, new_a={a[pi].in_(units.au)}")
                else:
                    p[pi].semimajor_axis = a[pi]
                    p[pi].eccentricity = e[pi]
                    stars[bi].companions.add_particle(p[pi])
    """
    for si in stars:
        if len(stars.companions)>0:
            for ci in si.companions:
                print(ci.eccentricity)
    """
    return stars
            
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
    binaries = find_companions(stars, stars)
    mmax = stars.mass.max()
    for bi in binaries:
        if len(bi.companions)>0:
            #print("star with M=", bi.mass.in_(units.MSun))
            #print("stellar companions:", bi.companions.semimajor_axis.in_(units.au),
            #      "e=", bi.companions.eccentricity)
            s = 300*bi.mass/mmax
            plt.scatter(bi.companions.semimajor_axis.value_in(units.au),
                        bi.companions.eccentricity, marker="*", s=s)
            #label="star-planet pair")

            
    single_freefloaters = Particles()
    single_freefloaters = planets.copy()
    #find companion planets to individual stars
    stars = find_companions(stars, planets)
    mmax = stars.mass.max()
    n_starplanet = []
    a_starplanet = []| units.au
    a_starsplanet = []| units.au
    a_starbplanet = []| units.au
    for bi in stars:
        if len(bi.companions)>0:
            #print("star with M=", bi.mass.in_(units.MSun))
            #print("companions:", bi.companions.semimajor_axis.in_(units.au),
            #      "e= ", bi.companions.eccentricity)
            s = 300*bi.mass/mmax
            plt.scatter(bi.companions.semimajor_axis.value_in(units.au),
                        bi.companions.eccentricity, marker="^", s=s)
            #label="star-planet pair")
            nsp = 0
            for ci in bi.companions:
                if ci.semimajor_axis<0.01|units.pc:
                    a_starplanet.append(ci.semimajor_axis)
                    if len(ci.parent)==0:
                        ci.parent.add_particle(bi)
                    nsp+=1
            if nsp>0:
                if nsp==1:
                    a_starsplanet.append(a_starplanet[-1])
                elif nsp==2:
                    a_starbplanet.append(a_starplanet[-2])
                n_starplanet.append(nsp)
            
    #find companion planets to individual planets
    planets = find_companions(planets, planets)
    mmax = planets.mass.max()
    a_jumbos = [] | units.au
    e_jumbos = []
    a_multiplanet = [] | units.au
    e_multiplanet = []
    q_jumbos = []
    for bi in planets:
        if len(bi.companions)>0:
            if len(bi.parent)==0:
                #print("planet with M=", bi.mass.in_(units.MJupiter))
                #print("companions:", bi.companions.semimajor_axis.in_(units.au),
                #      "e= ", bi.companions.eccentricity)
                a_jumbos.append(bi.companions.semimajor_axis)
                e_jumbos.append(bi.companions.eccentricity)
                for ci in bi.companions:
                    com = construct_center_of_mass_from_particle_pair(bi, ci)
                    print("M,m,a,e,x,y,z,vx,vy,vz")
                    print(f"JuMBO: {bi.mass.value_in(units.MJupiter)},{ci.mass.value_in(units.MJupiter)},{ci.semimajor_axis.value_in(units.au)},{ci.eccentricity},{com.x.value_in(units.pc)},{com.y.value_in(units.pc)},{com.z.value_in(units.pc)},{com.vx.value_in(units.kms)},{com.vy.value_in(units.kms)},{com.vz.value_in(units.kms)}")
                    
                    if ci.semimajor_axis<0.01|units.pc:
                        if bi.mass>ci.mass:
                            q_jumbos.append(ci.mass/bi.mass)
                        else:
                            q_jumbos.append(bi.mass/ci.mass)
                        s = 100*bi.mass/mmax
                        plt.scatter(bi.companions.semimajor_axis.value_in(units.au),
                                    bi.companions.eccentricity, marker="o", s=s)
            else:
                print("Binary Planet in orbit around star:")
                print("planet with M=", bi.mass.in_(units.MJupiter))
                print("companions:", bi.companions.mass.in_(units.MJupiter),
                      bi.companions.semimajor_axis.in_(units.au),
                      "e= ", bi.companions.eccentricity)
                print(f"Planet already claimed by a star at d={(bi.parent.position-bi.position).length().in_(units.au)}")
                print(f"distance of other planets: d={(bi.position-bi.companions.position).lengths().in_(units.au)}, companion: a={bi.companions.semimajor_axis.in_(units.au)}, e={bi.companions.eccentricity}")       
                print(f"distance of other companions to star: d={(bi.parent.position-bi.companions.position).lengths().in_(units.au)}")
                a, e = is_this_binary_planet_orbiting_a_star(stars, bi)
                print(f"Orbiting star at a={a.in_(units.au)}, e={e}")
                for i in range(len(a)):
                    a_multiplanet.append(a[i])
                    e_multiplanet.append(e[i])
    print("Multiplanets:", len(a_multiplanet.value_in(units.au)), len(e_multiplanet))
    plt.scatter(a_multiplanet.value_in(units.au), e_multiplanet, s=100, marker="v", c="r")

    print(f"Number of JuMBOs N={len(a_jumbos)}, <a>={np.mean(a_jumbos)}+-{np.std(a_jumbos)}, <e>= {np.mean(e_jumbos)}+-{np.std(e_jumbos)}")
    single_freefloaters = free_floating_planets(bodies)
    print(f"Number of try free floating planets N= {len(single_freefloaters)}")

    n_singlep = n_starplanet.count(1)
    n_doublep = n_starplanet.count(2)
    
    #amin amax, n, a, da, e, de, nff
    if "+" in o.filename:
        dirname = o.filename
    else:
        import os
        dirname = os.getcwd()
    #dirname = os.path.basename(path)
    print(dirname)
    if "+" in dirname:
        amin = float(dirname.split("+")[1])
        amax = float(dirname.split("+")[2])
        print(f"{amin},{amax},{int(len(a_jumbos))},{np.mean(a_jumbos).value_in(units.au)},{np.std(a_jumbos).value_in(units.au)},{np.mean(e_jumbos)},{np.std(e_jumbos)},{len(single_freefloaters)},{len(n_starplanet)},{n_singlep},{n_doublep}")
    else:
        print(f"{30},{3000},{int(len(a_jumbos))},{np.mean(a_jumbos).value_in(units.au)},{np.std(a_jumbos).value_in(units.au)},{np.mean(e_jumbos)},{np.std(e_jumbos)},{len(single_freefloaters)},{len(n_starplanet)},{n_singlep},{n_doublep}")
        
    print("Star-planet systems:", len(n_starplanet))
    #print(len(n_starplanet[n_starplanet==1]), len(n_starplanet[n_starplanet==2]))
    #print(n_starplanet)

    """
    stars_and_planets = Particles()
    stars_and_planets.add_particles(stars)
    stars_and_planets.add_particles(planets)
    write_set_to_file(stars_and_planets, "reduced.amuse", "amuse", version=2, close_file=True, overwrite_file=True, append_to_file=False)
    """
    
    if o.plot==True:
        plt.xlabel("a [au]")
        plt.ylabel("e")
        plt.xlim(1, 1.e+7)
        plt.semilogx()
        plt.legend()
        plt.show()

    if o.plot==True:
        asp = sorted(a_starplanet.value_in(units.au))
        f = np.linspace(0, 1, len(asp))
        plt.plot(asp, f, c='k', label="planetary systems")
        assp = sorted(a_starsplanet.value_in(units.au))
        f = np.linspace(0, 1, len(assp))
        plt.plot(assp, f, c='r', label="star single-planet")
        asbp = sorted(a_starbplanet.value_in(units.au))
        f = np.linspace(0, 1, len(asbp))
        plt.plot(asbp, f, c='b', label="star binary-planet")
        plt.ylim(0, 1)
        plt.xlabel("a [au]")
        plt.title(f"star-planet orbital distribution N={len(asp)}: Nsp={n_singlep}, Nbp={n_doublep}")
        plt.semilogx()
        plt.show()
        
    if o.plot==True:
        a = sorted(a_jumbos.value_in(units.au))
        f = np.linspace(0, 1, len(a))
        plt.plot(a, f, c='k')
        plt.ylim(0, 1)
        plt.xlabel("a [au]")
        plt.title(f"jumbo orbital distribution N={len(a)}")
        plt.show()

    if o.plot==True:
        plt.hist(q_jumbos, 50)
        plt.xlabel("q")
        plt.show()
        
    if o.plot==True:
        q = sorted(q_jumbos)
        f = np.linspace(0, 1, len(q))
        plt.plot(q, f, c='k')
        plt.ylim(0, 1)
        plt.xlabel("q")
        plt.title(f"jumbo mass ratio distribution N={len(q)}")
        plt.show()
    
    if o.plot==True:
        e = sorted(e_jumbos)
        f = np.linspace(0, 1, len(e))
        plt.plot(e, f, c='k')
        plt.xlim(0, 1)
        plt.ylim(0, 1)
        plt.xlabel("e")
        plt.title(f"jumbo eccentricity distribution N={len(e)}")
        plt.show()

    multiplanets = Particles()
    for pi in planets:
        if len(pi.companions)>0:
            multiplanets.add_particle(construct_center_of_mass_particle(pi))

    a_mp = [] |units.au
    e_mp = []
    mmax = stars.mass.max()
    n_starbinaryplanet = 0
    for mi in range(len(multiplanets)):
        total_masses = multiplanets[mi].mass + stars.mass
        rel_pos = multiplanets[mi].position-stars.position
        rel_vel = multiplanets[mi].velocity-stars.velocity
        sma, ecc, true_anomaly,\
            inc, long_asc_node, arg_per_mat =\
                get_orbital_elements_from_arrays(rel_pos,
                                                 rel_vel,
                                                 total_masses,
                                                 G=constants.G)

        for i in range(len(sma)):
            if ecc[i]<1:# and sma[i]<0.01|units.au:
                print("Bound Jumbos:", sma[i].in_(units.au), ecc[i])
                print("planet with M=", multiplanets[mi].mass.in_(units.MJupiter))
                print("companions:", sma[i].in_(units.au), "e= ", ecc[i])
                n_starbinaryplanet += 1
                
                a_mp.append(sma[i])
                e_mp.append(ecc[i])

    s = 100
    print(f"MPS={a_mp.value_in(units.au)}, {e_mp}")
    if o.plot==True:
        plt.scatter(a_mp.value_in(units.au), e_mp, marker="o", s=s)
        plt.xlabel("a [au]")
        plt.ylabel("e")
        plt.xlim(1, 1.e+7)
        plt.semilogx()
        plt.legend()
        plt.show()
    print("star binary-planet systems:", int(n_starbinaryplanet/2), len(e_mp))

    #print(f"{30},{3000},{int(len(a_jumbos))},{np.mean(a_jumbos).value_in(units.au)},{np.std(a_jumbos).value_in(units.au)},{np.mean(e_jumbos)},{np.std(e_jumbos)},{np.mean(q_jumbos)},{np.std(q_jumbos)},{len(single_freefloaters)},{len(n_starplanet)},{n_singlep},{n_doublep}, {int(n_starbinaryplanet/2)}")

    exit(0)
    xxx
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
    
