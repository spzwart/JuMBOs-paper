import glob
import matplotlib.pyplot as plt
import natsort
import numpy as np
import os
import pandas as pd
from scipy import stats
from itertools import combinations

from amuse.datamodel import Particles
from amuse.ext.LagrangianRadii import LagrangianRadii
from amuse.ext.orbital_elements import orbital_elements_from_binary
from amuse.io.base import read_set_from_file
from amuse.units import units, constants

from plotter_setup import PlotterSetup

class SystemAnimations(object):
    def __init__(self):
        self.clean_plot = PlotterSetup()
        self.bound_threshold = 1000 | units.au
        self.JuMBO_max_mass = 0.013 | units.MSun

        self.image_dir = "plotters/figures/system_evolution/output_movie"
        self.models = ["Fractal_rvir0.5_HighRes"]#, "Fractal_rvir0.5_FF_10Myr"]
        self.leg_label = r"Fractal, $R_{\mathrm{vir}} = 0.5$ pc with FF"
        self.run_choice = 0

    def system_evolution(self):
        model_iter = 0
        for model_ in self.models:

            path = "data/Simulation_Data/"+str(model_)
            output_file = "plotters/figures/system_evolution/movie"+str(model_)+".mp4"
            init_snap =  natsort.natsorted(glob.glob(os.path.join(str(path+"/simulation_snapshot/")+"*")))
            system_dir = init_snap[self.run_choice]
            system_run = natsort.natsorted(glob.glob(os.path.join(system_dir+"/*")))
            print("Making movie for: ", model_, " Config: ", init_snap[self.run_choice])

            scale_min = 0 | units.pc
            for snapshot_ in system_run:
                parti_data = read_set_from_file(snapshot_, "hdf5")
                parti_data.move_to_center()
                scale = LagrangianRadii(parti_data)[-2]
                scale_min = max(scale, scale_min)
            scale_min = scale_min.value_in(units.pc)

            if model_ == "Plummer_rvir0.5" \
                or model_ == "Plummer_rvir0.5_FF" \
                    or model_ == "Plummer_rvir1.0":
                scale_min *= 2
            
            parti_data_init = read_set_from_file(system_run[0], "hdf5")
            parti_data[parti_data.mass <= self.JuMBO_max_mass].name = "FF"
            key_set = set(parti_data_init.key)
            for i, key_fin in enumerate(parti_data.key):
                if key_fin not in key_set:
                    parti_data[i].name = "Merger"
                else:
                    if parti_data[i].mass > self.JuMBO_max_mass:
                        parti_data[i].name = "Star"

            components = parti_data.connected_components(threshold = self.bound_threshold)
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
                        if (eccentric < 1) and semimajor < self.bound_threshold:
                            multi_syst += 1
                            if max(bin_[0].mass, bin_[1].mass) <= self.JuMBO_max_mass:
                                parti_data[parti_data.key == bin_[0].key].name = "JuMBO"
                                parti_data[parti_data.key == bin_[1].key].name = "JuMBO"
                            elif min(bin_[0].mass, bin_[1].mass) > self.JuMBO_max_mass:
                                parti_data[parti_data.key == bin_[0].key].name = "Star-Star"
                                parti_data[parti_data.key == bin_[1].key].name = "Star-Star"
                            else:
                                parti_data[parti_data.key == bin_[0].key].name = "Star-Jupiter"
                                parti_data[parti_data.key == bin_[1].key].name == "Star-Jupiter"
                        keys = np.concatenate((keys, [bin_[0].key, bin_[1].key]), axis = None)

                    if multi_syst > 1:
                        for key_ in keys:
                            parti_data[parti_data.key == key_].name = "Multi"
                
            print("# Star: ", len(parti_data[parti_data.name == "Star"]))
            print("# FF: ", len(parti_data[parti_data.name == "FF"]))
            print("# Merger: ", len(parti_data[parti_data.name == "Merger"]))
            print("# Ghost: ", len(parti_data[parti_data.name == "Ghost"]))
            print("# JuMBO: ", len(parti_data[parti_data.name == "JuMBO"]))
            print("# S-S: ", len(parti_data[parti_data.name == "Star-Star"]))
            print("# J-S: ", len(parti_data[parti_data.name == "Star-Jupiter"]))
            print("# N>3: ", len(parti_data[parti_data.name == "Multi"]))

            dt = 1000/len(system_run)
            if "10Myr" in model_:
                dt *= 10

            snap_ = 0
            for snapshot_ in system_run:
                print("Reading iter: ", snap_, "/", len(system_run))
                fig, ax = plt.subplots(figsize=(8, 6))

                dt_snapshot = read_set_from_file(snapshot_, "hdf5")
                for parti_ in dt_snapshot:
                    parti_exist = parti_data[parti_data.key == parti_.key]
                    if len(parti_exist) == 0:# or len(parti_exist[parti_exist.Nej == 1]) != 0:
                        parti_.name = "Ghost"
                    if len(parti_exist) > 0:
                        parti_.name = parti_data[parti_data.key == parti_.key].name
                dt_snapshot.move_to_center()
                
                ax.scatter(1000, 1000, s=30, edgecolor = "Black", color = "gold", label = "Star")
                ax.scatter(1000, 1000, s=30, edgecolor = "Black", color = "red", label = "FF")
                ax.scatter(1000, 1000, s=30, edgecolor = "Black", color = "blue", label = "JuMBO")
                ax.scatter(1000, 1000, s=30, edgecolor = "Black", color = "darkviolet", label = "Merger")
                ax.scatter(1000, 1000, s=30, edgecolor = "Black", color = "black", label = "Ghost")
                ax.scatter(1000, 1000, s=30, edgecolor = "Black", color = "chocolate", label = "Star-Jupiter")
                ax.scatter(1000, 1000, s=30, edgecolor = "Black", color = "darkorange", label = "Star-Star")
                ax.scatter(1000, 1000, s=30, edgecolor = "Black", color = "cyan", label = r"$N\geq3$")

                if len(dt_snapshot[dt_snapshot.name == "Star"]) > 0:
                    ax.scatter(dt_snapshot[dt_snapshot.name == "Star"].x.value_in(units.pc), 
                               dt_snapshot[dt_snapshot.name == "Star"].y.value_in(units.pc), 
                               s=30*(dt_snapshot[dt_snapshot.name == "Star"].mass/dt_snapshot.mass.max())**0.25, 
                               edgecolor = "black", linewidth = 0.1, color = "gold")
                
                if len(dt_snapshot[dt_snapshot.name == "FF"]) > 0:
                    ax.scatter(dt_snapshot[dt_snapshot.name == "FF"].x.value_in(units.pc), 
                               dt_snapshot[dt_snapshot.name == "FF"].y.value_in(units.pc), 
                               s=60*(dt_snapshot[dt_snapshot.name == "FF"].mass/dt_snapshot.mass.max())**0.25, 
                               edgecolor = "black", linewidth = 0.1, color = "red")

                if len(dt_snapshot[dt_snapshot.name == "Star-Star"]) > 0:
                    ax.scatter(dt_snapshot[dt_snapshot.name == "Star-Star"].x.value_in(units.pc), 
                               dt_snapshot[dt_snapshot.name == "Star-Star"].y.value_in(units.pc), 
                               s=30*(dt_snapshot[dt_snapshot.name == "Star-Star"].mass/dt_snapshot.mass.max())**0.25, 
                               edgecolor = "black", linewidth = 0.1, color = "darkorange")

                if len(dt_snapshot[dt_snapshot.name == "Star-Jupiter"]) > 0:
                    ax.scatter(dt_snapshot[dt_snapshot.name == "Star-Jupiter"].x.value_in(units.pc), 
                               dt_snapshot[dt_snapshot.name == "Star-Jupiter"].y.value_in(units.pc), 
                               s=30*(dt_snapshot[dt_snapshot.name == "Star-Jupiter"].mass/dt_snapshot.mass.max())**0.25, 
                               edgecolor = "black", linewidth = 0.1, color = "chocolate")
                    
                if len(dt_snapshot[dt_snapshot.name == "Multi"]) > 0:
                    ax.scatter(dt_snapshot[dt_snapshot.name == "Multi"].x.value_in(units.pc), 
                               dt_snapshot[dt_snapshot.name == "Multi"].y.value_in(units.pc), 
                               s=30*(dt_snapshot[dt_snapshot.name == "Multi"].mass/dt_snapshot.mass.max())**0.25, 
                               edgecolor = "black", linewidth = 0.1, color = "cyan")

                if len(dt_snapshot[dt_snapshot.name == "Ghost"]) > 0:
                    ax.scatter(dt_snapshot[dt_snapshot.name == "Ghost"].x.value_in(units.pc), 
                               dt_snapshot[dt_snapshot.name == "Ghost"].y.value_in(units.pc), 
                               s=30*(dt_snapshot[dt_snapshot.name == "Ghost"].mass/dt_snapshot.mass.max())**0.25, 
                               edgecolor = "black", linewidth = 0.1, color = "black")
                            
                if len(dt_snapshot[dt_snapshot.name == "Merger"]) > 0:
                    ax.scatter(dt_snapshot[dt_snapshot.name == "Merger"].x.value_in(units.pc), 
                               dt_snapshot[dt_snapshot.name == "Merger"].y.value_in(units.pc), 
                               s=30*(dt_snapshot[dt_snapshot.name == "Merger"].mass/dt_snapshot.mass.max())**0.25, 
                               edgecolor = "black", linewidth = 0.1, color = "darkviolet")

                if len(dt_snapshot[dt_snapshot.name == "JuMBO"]) > 0:
                    ax.scatter(dt_snapshot[dt_snapshot.name == "JuMBO"].x.value_in(units.pc), 
                               dt_snapshot[dt_snapshot.name == "JuMBO"].y.value_in(units.pc), 
                               s=60*(dt_snapshot[dt_snapshot.name == "JuMBO"].mass/dt_snapshot.mass.max())**0.25, 
                               edgecolor = "black", linewidth = 0.1, color = "blue")

                ax.text(0, self.clean_plot.top, 
                        self.leg_label, 
                        horizontalalignment='left',
                        verticalalignment='bottom',
                        fontsize = self.clean_plot.axlabel_size,
                        transform=ax.transAxes)
                ax.text(0, self.clean_plot.top, 
                        r"$t = {:.1f}$ kyr".format(float(snap_*dt)),
                        horizontalalignment='left',
                        verticalalignment='top',
                        fontsize = self.clean_plot.axlabel_size,
                        transform=ax.transAxes)
                ax.set_xlim(-scale_min, scale_min)
                ax.set_ylim(-scale_min, scale_min)
                self.clean_plot.tickers(ax, "plot")
                ax.set_xlabel(r"$x$ [pc]", fontsize = self.clean_plot.axlabel_size)
                ax.set_ylabel(r"$y$ [pc]", fontsize = self.clean_plot.axlabel_size)
                ax.legend(prop={'size': self.clean_plot.axlabel_size}, 
                          bbox_to_anchor=(1.02, 0.5), loc='center left', fancybox=True)
                fig.savefig(self.image_dir+"/system_"+str(snap_)+".png", dpi=300, bbox_inches='tight')
                snap_ += 1
                plt.close(fig)
            if model_ != "Fractal_rvir0.5_HighRes":
                os.system(f"ffmpeg -r 20 -i {self.image_dir}/system_%d.png -c:v libx264 -preset slow -crf 18 -vf 'scale=1920:1080' -y {output_file}")
            else:
                os.system(f"ffmpeg -r 170 -i {self.image_dir}/system_%d.png -c:v libx264 -preset slow -crf 18 -vf 'scale=1920:1080' -y {output_file}")
            files = glob.glob(self.image_dir+"/*")
            
            for f in files:
                os.remove(f)
            model_iter += 1

    def JuMBO_evolution(self):
        model_iter = 0
        for model_ in self.models:
            dir = "data/Simulation_Data/"+str(model_)+"/"
            bin_keys = [ ]
            bin_choice = 10

            init_bins = natsort.natsorted(glob.glob(os.path.join(str(dir)+"initial_binaries/*")))[self.run_choice]
            with open(init_bins) as f:
                line = f.readlines()
                bin_idx = bin_choice*11
                bin_keys.append(float(line[bin_idx][6:]))
                bin_keys.append(float(line[bin_idx+1][6:]))
                    
            path = "data/Simulation_Data/"+str(model_)
            output_file = "plotters/figures/system_evolution/movie_JuMBO"+str(model_)+".mp4"
            init_snap =  natsort.natsorted(glob.glob(os.path.join(str(path+"/simulation_snapshot/")+"*")))
            system_dir = init_snap[self.run_choice]
            system_run = natsort.natsorted(glob.glob(os.path.join(system_dir+"/*")))
            snap_ = 0
            dt = 1000/len(system_run)

            for snapshot_ in system_run:
                print("Reading iter: ", snap_, "/", len(system_run))
                fig, ax = plt.subplots()
                left, width = .45, .5
                bottom, height = .58, .5
                right = left + width
                top = bottom + height

                parti_data = read_set_from_file(snapshot_, "hdf5")
                parti_data.move_to_center()

                stars = parti_data[parti_data.name == "star"]
                jmb = parti_data[parti_data.name == "JuMBOs"]
                focus_1 = jmb[jmb.key == bin_keys[0]]
                focus_2 = jmb[jmb.key == bin_keys[1]]
                jmb = jmb-focus_1-focus_2
                FF = parti_data[parti_data.name == "FF"]

                starx = (stars.x-focus_1.x).value_in(units.au)
                stary = (stars.y-focus_1.y).value_in(units.au)
                jmbx = (jmb.x-focus_1.x).value_in(units.au)
                jmby = (jmb.y-focus_1.y).value_in(units.au)
                focus_1x = (focus_1.x-focus_1.x).value_in(units.au)
                focus_1y = (focus_1.y-focus_1.y).value_in(units.au)
                focus_2x = (focus_2.x-focus_1.x).value_in(units.au)
                focus_2y = (focus_2.y-focus_1.y).value_in(units.au)
                FFx = (FF.x-focus_1.x).value_in(units.au)
                FFy = (FF.y-focus_1.y).value_in(units.au)

                ax.scatter(focus_1x, focus_1y, s=100*(focus_1.mass/parti_data.mass.max())**0.25, 
                           label = "JuMBO #1", edgecolor = "black", linewidth = 0.1, color = "red")
                ax.scatter(focus_2x, focus_2y, s=100*(focus_2.mass/parti_data.mass.max())**0.25, 
                           label = "JuMBO #2", edgecolor = "black", linewidth = 0.1, color = "blue")
                ax.scatter(starx, stary, s=100*(stars.mass/parti_data.mass.max())**0.25, 
                           edgecolor = "black", linewidth = 0.1, color = "gold")
                ax.scatter(jmbx, jmby, s=100*(jmb.mass/parti_data.mass.max())**0.25, 
                           edgecolor = "black", linewidth = 0.1, color = "brown")
                if len(FFx) > 0:
                    ax.scatter(FFx[0], FFy[0], s=100*(FF.mass/parti_data.mass.max())**0.25, label = "FF", 
                               edgecolor = "black", linewidth = 0.1, color = "cyan")
                    ax.scatter(FFx, FFy, s=100*(FF.mass/parti_data.mass.max())**0.25, edgecolor = "black", 
                               linewidth = 0.1, color = "cyan")
                ax.text(self.clean_plot.right, self.clean_plot.top, 
                        self.leg_label, 
                        horizontalalignment='right',
                        verticalalignment='bottom',
                        fontsize = self.clean_plot.axlabel_size,
                        transform=ax.transAxes)
                if "10Myr" in model_:
                    ax.text(self.clean_plot.right, self.clean_plot.top, 
                            r"$t = {:.1f}$ Myr".format(float(snap_*dt)/10),
                            horizontalalignment='right',
                            verticalalignment='top',
                            fontsize = self.clean_plot.axlabel_size,
                            transform=ax.transAxes)
                else:
                    ax.text(self.clean_plot.right, self.clean_plot.top, 
                            r"$t = {:.1f}$ kyr".format(float(snap_*dt)),
                            horizontalalignment='right',
                            verticalalignment='top',
                            fontsize = self.clean_plot.axlabel_size,
                            transform=ax.transAxes)
                ax.set_xlim(-1300, 1300)
                ax.set_ylim(-1300, 1300)
                self.clean_plot.tickers(ax, "plot")
                ax.set_xlabel(r"$x$ [au]", fontsize = self.clean_plot.axlabel_size)
                ax.set_ylabel(r"$y$ [au]", fontsize = self.clean_plot.axlabel_size)
                ax.legend(prop={'size': self.clean_plot.axlabel_size}, loc=3)
                fig.savefig(self.image_dir+"/system_JuMBO_"+str(snap_)+".png", dpi=300, bbox_inches='tight')
                snap_ += 1
                plt.close(fig)
            if model_ != "Fractal_rvir0.5_HighRes":
                os.system(f"ffmpeg -r 20 -i {self.image_dir}/system_%d.png -c:v libx264 -preset slow -crf 18 -vf 'scale=1920:1080' -y {output_file}")
            else:
                os.system(f"ffmpeg -r 70 -i {self.image_dir}/system_%d.png -c:v libx264 -preset slow -crf 18 -vf 'scale=1920:1080' -y {output_file}")
            files = glob.glob(self.image_dir+"/*")
            
            for f in files:
                os.remove(f)
            model_iter += 1

    def mix_sem_ecc(self):
        models = ["Fractal_rvir0.5", "Fractal_rvir0.5_FF", "Fractal_rvir1.0",
                            "Plummer_rvir0.5", "Plummer_rvir0.5_FF", "Plummer_rvir1.0",
                            "Fractal_rvir0.5_FF_10Myr", 
                            "Fractal_rvir0.5_Obs", "Fractal_rvir0.5_Obs_Circ",
                            "Fractal_rvir0.5_FF_Obs", "Plummer_rvir0.5_FF_Obs"]

        mass_img = "plotters/figures/system_evolution/mass_fluctuate/"
        orb_img = "plotters/figures/system_evolution/sem_ecc_fluctuate/"

        for model_ in models:
            fname = "data/Simulation_Data/"+model_+"/Processed_Data/Track_JuMBO/mixed_sys_data"
            events = pd.read_hdf(fname, 'Data')
            mass_fname = "plotters/figures/system_evolution/"+model_+"mass_evol.mp4"
            orb_fname = "plotters/figures/system_evolution/"+model_+"orb_evol.mp4"

            sim_run = [ ]
            mprim_arr = [ ]
            q_arr = [ ]
            semi_arr = [ ]
            ecc_arr = [ ]
            for col_ in events:
                sim_run.append(events.iloc[0][col_])
                mprim_arr.append(events.iloc[1][col_].value_in(units.MSun))
                q_arr.append(np.log10(events.iloc[2][col_]))
                semi_arr.append(events.iloc[3][col_].value_in(units.au))
                ecc_arr.append(events.iloc[4][col_])
            mprim_arr = np.asarray(mprim_arr)
            q_arr = np.asarray(q_arr)
            semi_arr = np.asarray(semi_arr)
            ecc_arr = np.asarray(ecc_arr)

            uq_runs = np.unique(sim_run)
            dt = 0
            for run_ in uq_runs:
                mprim_vals = mprim_arr[sim_run == run_]
                q_vals = q_arr[sim_run == run_]
                semi_vals = semi_arr[sim_run == run_]
                ecc_vals = ecc_arr[sim_run == run_]

                values = np.vstack([mprim_vals, q_vals])
                xx, yy = np.mgrid[0:10:200j, -4:0:200j]
                positions = np.vstack([xx.ravel(), yy.ravel()])
                kernel = stats.gaussian_kde(values, bw_method = "silverman")
                f = np.reshape(kernel(positions).T, xx.shape)

                fig, ax = plt.subplots()
                cfset = ax.contourf(xx, yy, f, cmap="Blues", levels = 7, zorder = 1)
                cset = ax.contour(xx, yy, f, colors = "k", levels = 7, zorder = 2)
                ax.clabel(cset, inline=1, fontsize=10)
                ax.set_xlabel(r"$M_{\mathrm{prim}}$", fontsize = self.clean_plot.axlabel_size)
                ax.set_ylabel(r"$\log_{10}q$", fontsize = self.clean_plot.axlabel_size)
                ax.set_xlim(0,10)
                ax.set_ylim(-4,0)
                self.clean_plot.tickers(ax, "hist")
                ax.text(self.clean_plot.right, self.clean_plot.top, 
                        r"$t = {:.1f}$ kyr".format(float(dt*10)),
                        horizontalalignment='right',
                        verticalalignment='top',
                        fontsize = self.clean_plot.axlabel_size,
                        transform=ax.transAxes)
                fig.savefig(mass_img+"system_"+str(dt)+".png", dpi=300, bbox_inches='tight')
                plt.close()

                values = np.vstack([semi_vals, ecc_vals])
                xx, yy = np.mgrid[0:1000:200j, 0:1:200j]
                positions = np.vstack([xx.ravel(), yy.ravel()])
                kernel = stats.gaussian_kde(values, bw_method = "silverman")
                f = np.reshape(kernel(positions).T, xx.shape)

                fig, ax = plt.subplots()
                cfset = ax.contourf(xx, yy, f, cmap="Blues", levels = 7, zorder = 1)
                cset = ax.contour(xx, yy, f, colors = "k", levels = 7, zorder = 2)
                ax.clabel(cset, inline=1, fontsize=10)
                ax.set_xlabel(r"$a$ [au]", fontsize = self.clean_plot.axlabel_size)
                ax.set_ylabel(r"$e$", fontsize = self.clean_plot.axlabel_size)
                ax.set_xlim(0,1000)
                ax.set_ylim(0,1)
                self.clean_plot.tickers(ax, "hist")
                ax.text(self.clean_plot.right, self.clean_plot.top, 
                        r"$t = {:.1f}$ kyr".format(float(dt*10)),
                        horizontalalignment='right',
                        verticalalignment='top',
                        fontsize = self.clean_plot.axlabel_size,
                        transform=ax.transAxes)
                fig.savefig(orb_img+"system_"+str(dt)+".png", dpi=300, bbox_inches='tight')
                plt.close()
                dt += 1

            for fname, img in zip([mass_fname, orb_fname], [mass_img, orb_img]):
                os.system(f"ffmpeg -r 10 -i {img}system_%d.png -c:v libx264 -preset slow -crf 18 -vf 'scale=1920:1080' -y {fname}")
                files = glob.glob(img+"/*")
                for f in files:
                    os.remove(f)

def sim_checker():
    """Check for correct initial conditions of each simulation"""

    models = ["Fractal_rvir0.5_FF_Obs", "Plummer_rvir0.5_FF_Obs"]

    for model_ in models:
        path = "data/Simulation_Data/"+str(model_)
        init_snap =  natsort.natsorted(glob.glob(os.path.join(str(path+"/simulation_snapshot/")+"*")))
        prev_key = -0
        for dir_ in init_snap:
            system_run = natsort.natsorted(glob.glob(os.path.join(dir_+"/*")))
            parti_data = read_set_from_file(system_run[0], "hdf5")
            print(dir_)
            print("nStar:", len(parti_data[parti_data.name == "star"]), "nJMO:", len(parti_data[parti_data.name != "star"]), "nJuMBO: ", len(parti_data[parti_data.name == "JuMBOs"]))
            print("m ~Jup:", len(parti_data[parti_data.mass <= (15 | units.MJupiter)]), "m ~star:", len(parti_data[parti_data.mass > (50 | units.MJupiter)]))
            Q = abs(parti_data.kinetic_energy()/parti_data.potential_energy())
            if Q < 0.4 or Q > 0.6:
                print("Error: Initialised Q wrong. Q = ", Q)
                STOP
            Rvir = parti_data.virial_radius().value_in(units.pc)
            if Rvir < 0.4 or Rvir > 1.1:
                print("Error: Initialised Rvir wrong. Rvir = ", Rvir, " pc")
                STOP
            if parti_data.key[0] == prev_key:
                print("Error: Double counting same simulation.")
                print(parti_data.key[0], prev_key)
                STOP
            prev_key = parti_data.key[0]
            print(Q, Rvir)
            """plt.scatter(parti_data.x.value_in(units.pc), parti_data.y.value_in(units.pc))
            plt.scatter(parti_data[parti_data.mass <= (14 | units.MJupiter)].x.value_in(units.pc), parti_data[parti_data.mass  <= (14 | units.MJupiter)].y.value_in(units.pc))
            plt.show()"""

animate = SystemAnimations()
animate.system_evolution()
animate.JuMBO_evolution()