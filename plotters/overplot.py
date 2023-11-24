import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia
import csv
import glob
from itertools import combinations
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import natsort
import numpy as np
import os

from amuse.datamodel import Particles
from amuse.ext.orbital_elements import orbital_elements_from_binary
from amuse.io.base import read_set_from_file
from amuse.units import units, constants

from plotter_setup import PlotterSetup


clean_plot = PlotterSetup()
JuMBO_max_mass = 15 | units.MJupiter
bound_threshold = 2000 | units.au

############OBSERVATIONAL DATA####################
query = """
        SELECT *
        FROM gaiadr2.gaia_source
        WHERE CONTAINS(POINT(ra, dec), CIRCLE(83.8, -5.23, 0.2)) = 1
        """

job = Gaia.launch_job(query)
result = job.get_results()
or_stars = result[(1000/result["parallax"] > 385) & (1000/result["parallax"] < 395)]

colors = ["black", "red", "darkviolet", "brown", "cyan", "blue"]
labels = ["Star", "FF", "Star-Star", "Star-Jup.", r"$N\geq 3$", "JuMBO"]

xObs = [ ]
yObs = [ ]
for ra_, dec_ in zip(or_stars["ra"], or_stars["dec"]):
    coord = SkyCoord(ra=ra_, dec=dec_, unit=(u.deg, u.deg), 
                     distance = 390*u.pc)
    xObs.append(coord.cartesian.x.value)
    yObs.append(coord.cartesian.y.value)
center = SkyCoord(ra="83.8", dec="-5 23 14.45", 
                  unit=(u.deg, u.deg), distance = 390*u.pc)

file = os.path.join("data/observations/src/obs_data.txt")
q_obs = [ ]
mprim_obs = [ ]
xJmb = [ ]
yJmb = [ ]

ra = [ ]
dec = [ ]
with open(file, 'r', newline ='') as file:
    csv_reader = csv.reader(file)
    for row_ in csv_reader:
        ra.append(float(row_[1]))
        dec.append(float(row_[2]))
        coord = SkyCoord(ra=float(row_[1]), dec=float(row_[2]), 
                         unit=(u.deg, u.deg), distance = 390*u.pc)
        xJmb.append(coord.cartesian.x.value)
        yJmb.append(coord.cartesian.y.value)

        mass1 = float(row_[3])
        mass2 = float(row_[5])
        q = min(mass1/mass2, mass2/mass1)
        q_obs.append(q)
        max_mass = max(mass1, mass2)
        mprim_obs.append(30*((1 | units.MSun)/(1 | units.MJupiter)*max_mass)**0.5)

xJmb -= center.cartesian.x.value
yJmb -= center.cartesian.y.value
xObs -= center.cartesian.x.value
yObs -= center.cartesian.y.value

################NUMERICAL DATA####################
model = "Fractal_rvir0.5_FF"
path = "data/Simulation_Data/"+str(model)
output_file = "plotters/figures/model_obs_snapshot.pdf"
configs = natsort.natsorted(glob.glob(os.path.join(str(path+"/simulation_snapshot/")+"*")))[1]
final_dt = natsort.natsorted(glob.glob(os.path.join(configs+"/*")))[-1]

parti_data = read_set_from_file(final_dt, "hdf5")
parti_data.move_to_center()
parti_data.name = "star"
parti_data[parti_data.mass <= JuMBO_max_mass].name = "FF"
components = parti_data.connected_components(threshold = bound_threshold)

q_jmb = [ ]
for c in components:
    if len(c) > 1:
        multi_syst = 0
        bin_combo = list(combinations(c, 2)) #All possible binaries in system
        keys = [ ]
        for bin_ in bin_combo:
            bin_sys = Particles()  
            bin_sys.add_particle(bin_[0])
            bin_sys.add_particle(bin_[1])
            kepler_elements = orbital_elements_from_binary(bin_sys, G=constants.G)
            semimajor = kepler_elements[2]
            eccentric = kepler_elements[3]

            if (eccentric < 1) and (semimajor < 0.5*bound_threshold):
                multi_syst += 1
                if max(bin_[0].mass, bin_[1].mass) <= JuMBO_max_mass:
                    q = min(bin_[1].mass/bin_[0].mass, bin_[0].mass/bin_[1].mass)
                    parti_data[parti_data.key == bin_[0].key].q = q
                    parti_data[parti_data.key == bin_[1].key].q = q
                    parti_data[parti_data.key == bin_[0].key].name = "JuMBO"
                    parti_data[parti_data.key == bin_[1].key].name = "JuMBO"
                elif min(bin_[0].mass, bin_[1].mass) > JuMBO_max_mass:
                    parti_data[parti_data.key == bin_[0].key].name = "Star-Star"
                    parti_data[parti_data.key == bin_[1].key].name = "Star-Star"
                else:
                    parti_data[parti_data.key == bin_[0].key].name = "Star-JMO"
                    parti_data[parti_data.key == bin_[1].key].name == "Star-JMO"
            keys = np.concatenate((keys, [bin_[0].key, bin_[1].key]), axis = None)

        if multi_syst > 1:
            for key_ in keys:
                parti_data[parti_data.key == key_].name = "Multi"

stars = parti_data[parti_data.name=="star"]
FF = parti_data[parti_data.name=="FF"]
SS = parti_data[parti_data.name=="Star-Star"]
SJ = parti_data[parti_data.name=="Star-Jupiter"]
multi = parti_data[parti_data.name=="Multi"]
jumbo = parti_data[parti_data.name=="JuMBO"]

theta = -2*np.pi/5
rotate = np.array([
    [np.cos(theta), -np.sin(theta)],
    [np.sin(theta),  np.cos(theta)]])

xy = [[ ], [ ]]
coords_obs = np.vstack((xObs, yObs))
coords_obs = np.dot(rotate, coords_obs)
rotated_coord = np.dot(rotate, coords_obs)
rot_obsx = rotated_coord[0,:]
rot_obsy = rotated_coord[1,:]
coords = np.vstack((parti_data.x.value_in(units.pc), 
                    parti_data.y.value_in(units.pc)))
rotated_coord = np.dot(rotate, coords)
rot_x = rotated_coord[0,:]
rot_y = rotated_coord[1,:]

ndot_size = 30*(jumbo.mass.value_in(units.MJupiter))**0.5
fig, ax = plt.subplots()
ax.hist2d(rot_x, rot_y, bins=100, range=([-2,2],[-2,2]), cmap='viridis')
ax.scatter(jumbo.x.value_in(units.pc), jumbo.y.value_in(units.pc),
           c = jumbo.q, edgecolors="red", s = ndot_size)
ax.scatter(xObs, yObs, color = "gold", alpha = 0.8)
colour_axes = ax.scatter(xJmb, yJmb, c = q_obs, s = mprim_obs, edgecolors="black")
cbar = plt.colorbar(colour_axes, ax=ax)
cbar.set_label(label = r'$q$', fontsize =  16)
ax.add_patch(Rectangle((-0.625, -0.4), 1.25, 0.8, alpha=1,  color = "white", fill=None))
clean_plot.tickers(ax, "hist")
ax.set_xlim(-2,2)
ax.set_ylim(-2,2)
ax.set_xlabel(r'$x$ [pc]', fontsize = 16)
ax.set_ylabel(r'$y$ [pc]', fontsize = 16)
plt.savefig("plotters/figures/overplotting_HEATMAP.pdf", 
            dpi=700, bbox_inches='tight')