import glob
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
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

from read_data import ReadData
from plotter_setup import PlotterSetup

class SystemAnimations(object):
    def __init__(self):
        self.clean_plot = PlotterSetup()
        self.bound_threshold = 1000 | units.au
        self.JuMBO_max_mass = 0.013 | units.MSun

        self.image_dir = "plotters/figures/system_evolution/output_movie"
        self.models = ["Fractal_rvir0.5_FF_10Myr", "Fractal_rvir0.5_HighRes"]
        self.leg_label = r"Fractal, $R_{\mathrm{vir}} = 0.5$ pc with FF"
        self.run_choice = 2

    def config_init_final(self):
        models = ["Fractal_rvir0.5", "Plummer_rvir0.5", "Fractal_rvir0.5_FF_10Myr"]

        for model_ in models:
            for dt_ in [0, -1]:
                path = "data/Simulation_Data/"+str(model_)
                fname = "plotters/figures/system_evolution/simulation_snap"+str(model_)+"_dt="+str(dt_)+".pdf"
                init_snap =  natsort.natsorted(glob.glob(os.path.join(str(path+"/simulation_snapshot/")+"*")))
                system_dir = init_snap[self.run_choice]
                system_run = natsort.natsorted(glob.glob(os.path.join(system_dir+"/*")))

                parti_data = read_set_from_file(system_run[dt_], "hdf5")
                parti_data_init = read_set_from_file(system_run[0], "hdf5")
                key_set = set(parti_data_init.key)
                for i, key_fin in enumerate(parti_data.key):
                    if key_fin not in key_set:
                        parti_data[i].name = "Merger"
                    else:
                        if parti_data[i].mass > self.JuMBO_max_mass:
                            parti_data[i].name = "Star"
                            
                if dt_ == -1:
                    parti_data[parti_data.mass <= self.JuMBO_max_mass].name = "FF"
                parti_data[parti_data.mass > self.JuMBO_max_mass].name = "Star"

                components = parti_data.connected_components(threshold = self.bound_threshold)
                for c in components:
                    if len(c) > 1:
                        multi_syst = 0

                        bin_combo = list(combinations(c, 2))
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
                                if dt_ == -1:
                                    if max(bin_[0].mass, bin_[1].mass) <= self.JuMBO_max_mass:
                                        parti_data[parti_data.key == bin_[0].key].name = "JuMBOs"
                                        parti_data[parti_data.key == bin_[1].key].name = "JuMBOs"
                            keys = np.concatenate((keys, [bin_[0].key, bin_[1].key]), axis = None)

                        if multi_syst > 1:
                            for key_ in keys:
                                parti_data[parti_data.key == key_].name = "Multi"

                stars = parti_data[parti_data.name == "Star"]
                FF = parti_data[parti_data.name == "FF"]
                merger = parti_data[parti_data.name == "Merger"]
                jumbo = parti_data[parti_data.name == "JuMBOs"]

                fig, ax = plt.subplots(figsize=(8, 6))
                if len(stars) > 0:
                    ax.scatter(stars.x.value_in(units.pc), stars.y.value_in(units.pc), 
                            s=150*(stars.mass/parti_data.mass.max())**0.5, 
                            edgecolor = "black", linewidth = 0.1, color = "gold", alpha =0.75)
                
                if len(FF) > 0:
                    ax.scatter(FF.x.value_in(units.pc), FF.y.value_in(units.pc), 
                            edgecolor = "black", linewidth = 0.1, color = "red",
                            s = 20, alpha = 0.75)
                            
                if len(merger) > 0:
                    ax.scatter(merger.x.value_in(units.pc), merger.y.value_in(units.pc), 
                            s=150*(merger.mass/parti_data.mass.max())**0.5, 
                            edgecolor = "black", linewidth = 0.1, color = "darkviolet")

                if len(jumbo) > 0:
                    ax.scatter(jumbo.x.value_in(units.pc), jumbo.y.value_in(units.pc), 
                            edgecolor = "black", linewidth = 0.1, 
                            s = 20, color = "blue")

                xmax = max(abs(parti_data.x))
                ymax = max(abs(parti_data.y))
                lims = max(xmax, ymax)

                ax.set_xlim(-1.02*lims.value_in(units.pc), 1.02*lims.value_in(units.pc))
                ax.set_ylim(-1.02*lims.value_in(units.pc), 1.02*lims.value_in(units.pc))
                ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
                ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
                self.clean_plot.tickers(ax, "plot")
                ax.set_xlabel(r"$x$ [pc]", fontsize=self.clean_plot.axlabel_size)
                ax.set_ylabel(r"$y$ [pc]", fontsize=self.clean_plot.axlabel_size)
                fig.savefig(fname, dpi=300, bbox_inches='tight')
                plt.close(fig)

    def system_evolution(self):
        model_iter = 0
        for model_ in self.models:
            
            path = "data/Simulation_Data/"+str(model_)
            output_file = "plotters/figures/system_evolution/movie"+str(model_)+".mp4"
            init_snap =  natsort.natsorted(glob.glob(os.path.join(str(path+"/simulation_snapshot/")+"*")))
            print("Making movie for: ", model_, " Config: ", init_snap[self.run_choice])
            system_dir = init_snap[self.run_choice]
            system_run = natsort.natsorted(glob.glob(os.path.join(system_dir+"/*")))

            scale_min = 0 | units.pc
            for snapshot_ in system_run:
                parti_data = read_set_from_file(snapshot_, "hdf5")
                parti_data.move_to_center()
                scale = LagrangianRadii(parti_data)[-2]
                scale_min = max(scale, scale_min)
            scale_min = scale_min.value_in(units.pc)

            if "Plummer" in model_:
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

                    bin_combo = list(combinations(c, 2))
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
                                parti_data[parti_data.key == bin_[0].key].name = "JuMBOs"
                                parti_data[parti_data.key == bin_[1].key].name = "JuMBOs"
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
            print("# JuMBO: ", len(parti_data[parti_data.name == "JuMBOs"]))
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
                
                ax.scatter(np.nan, np.nan, s=30, edgecolor = "Black", 
                          color = "gold", label = "Star")
                ax.scatter(np.nan, np.nan, s=30, edgecolor = "Black", 
                           color = "red", label = "FF")
                ax.scatter(np.nan, np.nan, s=30, edgecolor = "Black", 
                           color = "blue", label = "JuMBOs")
                ax.scatter(np.nan, np.nan, s=30, edgecolor = "Black", 
                           color = "darkviolet", label = "Merger")
                ax.scatter(np.nan, np.nan, s=30, edgecolor = "Black", 
                           color = "black", label = "Ghost")
                ax.scatter(np.nan, np.nan, s=30, edgecolor = "Black", 
                           color = "chocolate", label = "Star-Jupiter")
                ax.scatter(np.nan, np.nan, s=30, edgecolor = "Black", 
                           color = "darkorange", label = "Star-Star")
                ax.scatter(np.nan, np.nan, s=30, edgecolor = "Black", 
                           color = "cyan", label = r"$N\geq3$")
                
                star = dt_snapshot[dt_snapshot.name == "Star"]
                FF = dt_snapshot[dt_snapshot.name == "FF"]
                star_star = dt_snapshot[dt_snapshot.name == "Star-Star"]
                star_jmo = dt_snapshot[dt_snapshot.name == "Star-Jupiter"]
                multi_sys = dt_snapshot[dt_snapshot.name == "Multi"]
                ghost = dt_snapshot[dt_snapshot.name == "Ghost"]
                merger = dt_snapshot[dt_snapshot.name == "Merger"]
                jumbo = dt_snapshot[dt_snapshot.name == "JuMBOs"]
                if len(star) > 0:
                    ax.scatter(star.x.value_in(units.pc), star.y.value_in(units.pc), 
                               s=30*(star.mass/dt_snapshot.mass.max())**0.25, 
                               edgecolor = "black", linewidth = 0.1, color = "gold")
                
                if len(FF) > 0:
                    ax.scatter(FF.x.value_in(units.pc), FF.y.value_in(units.pc), 
                               s=60*(FF.mass/dt_snapshot.mass.max())**0.25, 
                               edgecolor = "black", linewidth = 0.1, color = "red")

                if len(star_star) > 0:
                    ax.scatter(star_star.x.value_in(units.pc), star_star.y.value_in(units.pc), 
                               s=30*(star_star.mass/dt_snapshot.mass.max())**0.25, 
                               edgecolor = "black", linewidth = 0.1, color = "darkorange")

                if len(star_jmo) > 0:
                    ax.scatter(star_jmo.x.value_in(units.pc), star_jmo.y.value_in(units.pc), 
                               s=30*(star_jmo.mass/dt_snapshot.mass.max())**0.25, 
                               edgecolor = "black", linewidth = 0.1, color = "chocolate")
                    
                if len(multi_sys) > 0:
                    ax.scatter(multi_sys.x.value_in(units.pc), multi_sys.y.value_in(units.pc), 
                               s=30*(multi_sys.mass/dt_snapshot.mass.max())**0.25, 
                               edgecolor = "black", linewidth = 0.1, color = "cyan")

                if len(ghost) > 0:
                    ax.scatter(ghost.x.value_in(units.pc), ghost.y.value_in(units.pc), 
                               s=30*(ghost.mass/dt_snapshot.mass.max())**0.25, 
                               edgecolor = "black", linewidth = 0.1, color = "black")
                            
                if len(merger) > 0:
                    ax.scatter(merger.x.value_in(units.pc), merger.y.value_in(units.pc), 
                               s=30*(merger.mass/dt_snapshot.mass.max())**0.25, 
                               edgecolor = "black", linewidth = 0.1, color = "darkviolet")

                if len(jumbo) > 0:
                    ax.scatter(jumbo.x.value_in(units.pc), jumbo.y.value_in(units.pc), 
                               s=60*(jumbo.mass/dt_snapshot.mass.max())**0.25, 
                               edgecolor = "black", linewidth = 0.1, color = "blue")

                ax.text(0, self.clean_plot.top, self.leg_label, 
                        horizontalalignment='left',
                        verticalalignment='bottom',
                        fontsize=self.clean_plot.axlabel_size,
                        transform=ax.transAxes)

                if model_ == "Fractal_rvir0.5_FF_10Myr":
                    ax.text(0, self.clean_plot.top, 
                            r"$t = {:.1f}$ Myr".format(float(snap_*dt*10**-3)),
                            horizontalalignment='left',
                            verticalalignment='top',
                            fontsize=self.clean_plot.axlabel_size,
                            transform=ax.transAxes)
                else:
                    ax.text(0, self.clean_plot.top, 
                            r"$t = {:.1f}$ kyr".format(float(snap_*dt)),
                            horizontalalignment='left',
                            verticalalignment='top',
                            fontsize=self.clean_plot.axlabel_size,
                            transform=ax.transAxes)

                ax.set_xlim(-scale_min, scale_min)
                ax.set_ylim(-scale_min, scale_min)
                self.clean_plot.tickers(ax, "plot")
                ax.set_xlabel(r"$x$ [pc]", fontsize=self.clean_plot.axlabel_size)
                ax.set_ylabel(r"$y$ [pc]", fontsize=self.clean_plot.axlabel_size)
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

    def JuMBO_evolution(self, model):
        def be_calc(m1,m2,a):
            return -(constants.G*m1*m2)/(2*abs(a))

        path = "data/Simulation_Data/"+str(model)
        output_file = "plotters/figures/system_evolution/movie_JuMBOs"+str(model)+".mp4"
        init_snap =  natsort.natsorted(glob.glob(os.path.join(str(path+"/simulation_snapshot/")+"*")))
        system_dir = init_snap[self.run_choice]
        system_run = natsort.natsorted(glob.glob(os.path.join(system_dir+"/*")))
        snap_ = 0
        dt = 1000/len(system_run)
        
        parti_data = read_set_from_file(system_run[-1], "hdf5")
        parti_data = parti_data[parti_data.mass <= self.JuMBO_max_mass]
        components = parti_data.connected_components(threshold = self.bound_threshold)
        no_det = 0
        semimajor = 0 | units.au
        for c in components:
            if len(c) > 1 and semimajor <= 80 | units.au:
                bin_combo = list(combinations(c, 2))
                for bin_ in bin_combo:
                    bin_sys = Particles()  
                    bin_sys.add_particle(bin_[0])
                    bin_sys.add_particle(bin_[1])

                    kepler_elements = orbital_elements_from_binary(bin_sys, G=constants.G)
                    semimajor = kepler_elements[2]
                    eccentric = kepler_elements[3]
                    print(semimajor.value_in(units.au))
                    if (eccentric < 1) and semimajor < self.bound_threshold:
                        bin_keys = [bin_[0].key, bin_[1].key]
                        no_det += 1

        max_BE = 0
        min_BE = -10**100
        BE_arr = [ ]
        for snapshot_ in system_run:
            parti_data = read_set_from_file(snapshot_, "hdf5")
            parti_data.move_to_center()

            stars = parti_data[parti_data.name == "star"]
            jmb = parti_data[parti_data.name == "JuMBOs"]
            focus_1 = jmb[jmb.key == bin_keys[0]]
            focus_2 = jmb[jmb.key == bin_keys[1]]

            bin_sys = Particles()
            bin_sys.add_particle(focus_1)
            bin_sys.add_particle(focus_2)
            kepler_elements = orbital_elements_from_binary(bin_sys, G=constants.G)
            semimajor = kepler_elements[2]
            BE = be_calc(focus_1.mass, focus_2.mass, semimajor)[0]
            BE = BE.value_in(units.J)
            BE_arr.append(BE)
            max_BE = min(max_BE, BE)
            min_BE = max(min_BE, BE)
        normalise = plt.Normalize(min_BE/max_BE, max_BE/max_BE)

        for snapshot_ in system_run:
            print("Reading iter: ", snap_, "/", len(system_run))

            parti_data = read_set_from_file(snapshot_, "hdf5")
            parti_data.move_to_center()

            stars = parti_data[parti_data.name == "star"]
            jmb = parti_data[parti_data.name == "JuMBOs"]
            if jmb[jmb.key == bin_keys[0]].mass > jmb[jmb.key == bin_keys[1]].mass:
                focus_1 = jmb[jmb.key == bin_keys[0]]
                focus_2 = jmb[jmb.key == bin_keys[1]]
            else:
                focus_1 = jmb[jmb.key == bin_keys[1]]
                focus_2 = jmb[jmb.key == bin_keys[0]]

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

            fig, ax = plt.subplots()
            ax.scatter(focus_1x, focus_1y, s=200*(focus_1.mass/parti_data.mass.max())**0.25, 
                       edgecolor = "black", norm = normalise,
                       linewidth = 0.1, c = BE_arr[snap_]/max_BE)
            colour_axes = ax.scatter(focus_2x, focus_2y, s=200*(focus_2.mass/parti_data.mass.max())**0.25, 
                                     edgecolor = "black", norm = normalise, 
                                     linewidth = 0.1, c = BE_arr[snap_]/max_BE)
            cbar = plt.colorbar(colour_axes, ax=ax)
            cbar.set_label(label = r'$B_E/B_{E,\ \mathrm{max}}$', fontsize=self.clean_plot.axlabel_size)
            ax.scatter(starx[0], stary[0], s=100*(stars[0].mass/parti_data.mass.max())**0.25, 
                       label = "Star", edgecolor = "black", linewidth = 0.1, color = "gold")
            ax.scatter(starx[1:], stary[1:], s=100*(stars[1:].mass/parti_data.mass.max())**0.25, 
                       edgecolor = "black", linewidth = 0.1, color = "gold")
            ax.scatter(jmbx[0], jmby[0], s=200*(jmb[0].mass/parti_data.mass.max())**0.25, 
                       label = "JuMBO", edgecolor = "black", linewidth = 0.1, color = "blue")
            ax.scatter(jmbx[1:], jmby[1:], s=200*(jmb[1:].mass/parti_data.mass.max())**0.25, 
                       edgecolor = "black", linewidth = 0.1, color = "blue")
            if len(FFx) > 0:
                ax.scatter(FFx[0], FFy[0], s=200*(FF[0].mass/parti_data.mass.max())**0.25, label = "FF", 
                            edgecolor = "black", linewidth = 0.1, color = "red")
                ax.scatter(FFx[1:], FFy[1:], s=200*(FF[1:].mass/parti_data.mass.max())**0.25, 
                           edgecolor = "black", linewidth = 0.1, color = "red")
            ax.text(self.clean_plot.right, self.clean_plot.top, 
                    self.leg_label, 
                    horizontalalignment='right',
                    verticalalignment='bottom',
                    fontsize=self.clean_plot.axlabel_size,
                    transform=ax.transAxes)
            if "10Myr" in model:
                ax.text(self.clean_plot.right, self.clean_plot.top, 
                        r"$t = {:.1f}$ Myr".format(float(snap_*dt)/100),
                        horizontalalignment='right',
                        verticalalignment='top',
                        fontsize=self.clean_plot.axlabel_size,
                        transform=ax.transAxes)
            else:
                ax.text(self.clean_plot.right, self.clean_plot.top, 
                        r"$t = {:.1f}$ kyr".format(float(snap_*dt)),
                        horizontalalignment='right',
                        verticalalignment='top',
                        fontsize=self.clean_plot.axlabel_size,
                        transform=ax.transAxes)
            ax.set_xlim(-400,400)#(-1300, 1300)
            ax.set_ylim(-400,400)#(-1300, 1300)
            self.clean_plot.tickers(ax, "plot")
            ax.set_xlabel(r"$x$ [au]", fontsize=self.clean_plot.axlabel_size)
            ax.set_ylabel(r"$y$ [au]", fontsize=self.clean_plot.axlabel_size)
            ax.legend(prop={'size': self.clean_plot.axlabel_size}, loc=3)
            fig.savefig(self.image_dir+"/system_JuMBO_"+str(snap_)+".png", dpi=300, bbox_inches='tight')
            snap_ += 1
            plt.close(fig)
        if model != "Fractal_rvir0.5_HighRes":
            os.system(f"ffmpeg -r 20 -i {self.image_dir}/system_JuMBO_%d.png -c:v libx264 -preset slow -crf 18 -vf 'scale=1920:1080' -y {output_file}")
        else:
            os.system(f"ffmpeg -r 70 -i {self.image_dir}/system_JuMBO_%d.png -c:v libx264 -preset slow -crf 18 -vf 'scale=1920:1080' -y {output_file}")
        files = glob.glob(self.image_dir+"/*")
        
        print("Removing JuMBO files")
        for f in files:
            os.remove(f)

    def mix_sem_ecc(self):
        models = ["Fractal_rvir0.5", "Fractal_rvir0.5_FF",
                  "Plummer_rvir0.5", "Plummer_rvir0.5_FF"]

        mass_img = "plotters/figures/system_evolution/mass_fluctuate/"
        orb_img = "plotters/figures/system_evolution/sem_ecc_fluctuate/"

        for model_ in models:
            path = "data/Simulation_Data/"+str(model_)
            traj_files = glob.glob(os.path.join(str(path+"/simulation_snapshot/")+"*"))
            NSims = len(traj_files)

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

            uq_runs, Nsysts = np.unique(sim_run, return_counts=True)
            
            dt = 0
            for run_, Nsyst in zip(uq_runs, Nsysts):
                mprim_vals = mprim_arr[sim_run == run_]
                q_vals = q_arr[sim_run == run_]
                semi_vals = semi_arr[sim_run == run_]
                ecc_vals = ecc_arr[sim_run == run_]
                
                mprim_vals = mprim_vals[np.isfinite(q_vals)]
                q_vals = q_vals[np.isfinite(q_vals)]

                values = np.vstack([mprim_vals, q_vals])
                xx, yy = np.mgrid[0:10:200j, -4:0:200j]
                positions = np.vstack([xx.ravel(), yy.ravel()])
                kernel = stats.gaussian_kde(values, bw_method = "silverman")
                f = np.reshape(kernel(positions).T, xx.shape)

                fig, ax = plt.subplots()
                cfset = ax.contourf(xx, yy, f, cmap="Blues", levels = 7, zorder = 1)
                cset = ax.contour(xx, yy, f, colors = "k", levels = 7, zorder = 2)
                ax.clabel(cset, inline=1, fontsize=10)
                ax.set_xlabel(r"$M_{\mathrm{prim}} [M_{\mathrm{\odot}}]$", fontsize=self.clean_plot.axlabel_size)
                ax.set_ylabel(r"$\log_{10}q$", fontsize=self.clean_plot.axlabel_size)
                ax.set_xlim(0,10)
                ax.set_ylim(-4,0)
                self.clean_plot.tickers(ax, "hist")
                ax.text(self.clean_plot.right, self.clean_plot.top-1, 
                        r"$\langle N_{{\mathrm{{systs}}}}\rangle = {:.1f}$".format(float(Nsyst/NSims)),
                        horizontalalignment='right',
                        verticalalignment='top',
                        fontsize=self.clean_plot.axlabel_size,
                        transform=ax.transAxes)
                ax.text(self.clean_plot.right, self.clean_plot.top-0.95, 
                        r"$t = {:.1f}$ kyr".format(float(dt*10)),
                        horizontalalignment='right',
                        verticalalignment='top',
                        fontsize=self.clean_plot.axlabel_size,
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
                ax.set_xlabel(r"$a$ [au]", fontsize=self.clean_plot.axlabel_size)
                ax.set_ylabel(r"$e$", fontsize=self.clean_plot.axlabel_size)
                ax.set_xlim(0,1000)
                ax.set_ylim(0,1)
                self.clean_plot.tickers(ax, "hist")
                ax.text(self.clean_plot.right, self.clean_plot.top-1, 
                        r"$\langle N_{{\mathrm{{systs}}}}\rangle = {:.1f}$".format(float(Nsyst/NSims)),
                        horizontalalignment='right',
                        verticalalignment='top',
                        fontsize=self.clean_plot.axlabel_size,
                        transform=ax.transAxes)
                ax.text(self.clean_plot.right, self.clean_plot.top-0.95, 
                        r"$t = {:.1f}$ kyr".format(float(dt*10)),
                        horizontalalignment='right',
                        verticalalignment='top',
                        fontsize=self.clean_plot.axlabel_size,
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
            print("nStar:", len(parti_data[parti_data.name == "star"]), end=" ")
            print("nJMO:", len(parti_data[parti_data.name != "star"]), end=" ")
            print("nJuMBO: ", len(parti_data[parti_data.name == "JuMBOs"]))
            print("m ~Jup:", len(parti_data[parti_data.mass <= (15 | units.MJupiter)]), end=" ")
            print("m ~star:", len(parti_data[parti_data.mass > (50 | units.MJupiter)]))
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



animate = SystemAnimations()
#animate.config_init_final()
#animate.mix_sem_ecc()
animate.JuMBO_evolution(animate.models[1])
animate.system_evolution()