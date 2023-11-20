import csv
import glob
import itertools
import matplotlib.pyplot as plt
import natsort
import numpy as np
import os
import pandas as pd
import statsmodels.api as sm
import warnings
from scipy import stats

from amuse.datamodel import Particles
from amuse.ext.orbital_elements import orbital_elements_from_binary
from amuse.io.base import read_set_from_file
from amuse.units import units, constants

from plotter_func import PlotterFunctions
from plotter_setup import PlotterSetup
from read_data import ReadData

warnings.filterwarnings("ignore")


class FinalInitialProperties(object):
    def __init__(self):
        self.clean_plot = PlotterSetup()
        self.plot_func = PlotterFunctions()
        self.JuMBO_max_mass = 0.013 | units.MSun
        self.Star_min_mass = 0.08 | units.MSun

        self.models =  ["Fractal_rvir0.5", "Fractal_rvir0.5_FF", "Fractal_rvir1.0",
                        "Plummer_rvir0.5", "Plummer_rvir0.5_FF", "Plummer_rvir1.0",
                        "Fractal_rvir0.5_FF_10Myr", "Fractal_rvir0.5_FFOnly",
                        "Fractal_rvir0.5_Obs", "Fractal_rvir0.5_Obs_Circ",
                        "Fractal_rvir0.5_FF_Obs", "Plummer_rvir0.5_FF_Obs"]
                      
        self.leg_label = [r"Fractal, $R_{\mathrm{vir}} = 0.5$ pc", 
                          r"Fractal, $R_{\mathrm{vir}} = 0.5$ pc with FF", 
                          r"Fractal, $R_{\mathrm{vir}} = 1.0$ pc", 
                          r"Plummer, $R_{\mathrm{vir}} = 0.5$ pc", 
                          r"Plummer, $R_{\mathrm{vir}} = 0.5$ pc with FF", 
                          r"Plummer, $R_{\mathrm{vir}} = 1.0$ pc",
                          r"Fractal, $R_{\mathrm{vir}} = 0.5$ pc with FF",
                          r"Fractal, $R_{\mathrm{vir}} = 0.5$ pc",
                          r"Fractal, $R_{\mathrm{vir}} = 0.5$ pc",
                          r"Fractal, $R_{\mathrm{vir}} = 0.5$ pc with FF",
                          r"Plummer, $R_{\mathrm{vir}} = 0.5$ pc with FF", ]
        self.linestyles = ["dashed", "dash dot dotted", 
                           "dotted", (0, (3, 1, 1, 1, 1, 1))]
        self.colours = ["red", "blue", "dodgerblue", "pink"]
        
    def initialise(self, crop_):
        """Initialise complete data set"""

        self.file_count = [ ]

        self.init_bkeys = [[ ] for i in self.models]
        self.init_bsemi = [[ ] for i in self.models]
        self.init_bmass = [[ ] for i in self.models]
        self.init_becce = [[ ] for i in self.models]
        self.init_bincl = [[ ] for i in self.models]
        
        self.jmb_idx = [[ ] for i in self.models]
        self.fin_bsim_iter = [[ ] for i in self.models]
        self.fin_bkeys = [[ ] for i in self.models]
        self.fin_bsemi = [[ ] for i in self.models]
        self.fin_bproj = [[ ] for i in self.models]
        self.fin_bmass = [[ ] for i in self.models]
        self.fin_becce = [[ ] for i in self.models]
        self.fin_bincl = [[ ] for i in self.models]
        self.fin_btype = [[ ] for i in self.models]

        self.fin_msim_iter = [[ ] for i in self.models]
        self.fin_mkeys = [[ ] for i in self.models]
        self.fin_mtype = [[ ] for i in self.models]
        self.fin_mmass = [[ ] for i in self.models]
        self.msyst_pops = [[ ] for i in self.models]
        self.msyst_types = [[ ] for i in self.models]
        self.mkeys = [[ ] for i in self.models]

        model_iter = 0
        for model_ in self.models:
            dir = "data/Simulation_Data/"+str(model_)+"/"
            self.file_count.append(len(glob.glob(os.path.join(dir+"initial_binaries/*"))))
            if self.file_count[model_iter] >= 0:
                data = ReadData()
                data.read_init_data(dir)
                self.init_bkeys[model_iter] = data.init_bkeys
                self.init_bsemi[model_iter] = data.init_bsemi
                self.init_bmass[model_iter] = data.init_bmass
                self.init_becce[model_iter] = data.init_becce
                self.init_bincl[model_iter] = data.init_bincl

                data.read_final_data(dir, crop_)
                self.jmb_idx[model_iter] = data.jmb_idx
                self.fin_bsim_iter[model_iter] = data.fin_bsim_iter
                self.fin_bkeys[model_iter] = data.fin_bkeys
                self.fin_bsemi[model_iter] = data.fin_bsemi
                self.fin_bproj[model_iter] = data.fin_bproj
                self.fin_bmass[model_iter] = data.fin_bmass
                self.fin_becce[model_iter] = data.fin_becce
                self.fin_bincl[model_iter] = data.fin_bincl
                self.fin_btype[model_iter] = data.fin_btype

                data.read_final_multi_data(dir, crop_)
                self.fin_msim_iter[model_iter] = data.fin_msim_iter
                self.fin_mkeys[model_iter] = data.fin_mkeys
                self.fin_mtype[model_iter] = data.fin_mtype
                self.fin_mmass[model_iter] = data.fin_mmass
                moutput = self.plot_func.multi_sys_flag(self.fin_mkeys[model_iter], 
                                                        self.fin_mmass[model_iter],
                                                        self.fin_mtype[model_iter])
                self.msyst_pops[model_iter] = moutput[0]
                self.msyst_types[model_iter] = moutput[1] 
                self.mkeys[model_iter] = moutput[2]
                
            model_iter += 1

    def process_final_data(self, dt_choice):
        data = ReadData()
        if (dt_choice):
            self.models =  ["Fractal_rvir0.5", "Fractal_rvir0.5_FF", "Fractal_rvir1.0",
                        "Plummer_rvir0.5", "Plummer_rvir0.5_FF", "Plummer_rvir1.0",
                        "Fractal_rvir0.5_FF_10Myr", "Fractal_rvir0.5_FFOnly",
                        "Fractal_rvir0.5_Obs", "Fractal_rvir0.5_Obs_Circ",
                        "Fractal_rvir0.5_FF_Obs"]
        
        for model_ in ["Fractal_rvir0.5"]:#self.models:
            print("Processing data for ", model_)
            data.process_final_snapshot(model_, dt_choice)

    def event_statistics(self, dt_crop):    
        """Extract statistics on mergers + ejections"""

        for model_ in self.models:
            if (dt_crop):
                tname = model_+'_events_crop.txt'
            else:
                tname = model_+'_events.txt'
                
            path = "data/Simulation_Data/"+str(model_)
            traj_files = natsort.natsorted(glob.glob(os.path.join(str(path+"/simulation_snapshot/")+"*")))

            ejec_jmo = np.zeros(len(traj_files))
            ejec_star = np.zeros(len(traj_files))
            ejec_jj = np.zeros(len(traj_files))
            ejec_js = np.zeros(len(traj_files))
            ejec_ss = np.zeros(len(traj_files))
            merge_jj = np.zeros(len(traj_files))
            merge_js = np.zeros(len(traj_files))
            merge_ss = np.zeros(len(traj_files))

            fiter = 0
            if model_ == "Fractal_rvir0.5_FF_10Myr": #No event_files saved
                for config_ in traj_files:
                    sim_snapshot = natsort.natsorted(glob.glob(os.path.join(config_+"/*")))
                    
                    prev_set = read_set_from_file(sim_snapshot[0], "hdf5")
                    prev_len = len(prev_set)
                    for dt in sim_snapshot:
                        print("Processing event statistics for ", dt)
                        parti_set = read_set_from_file(dt, "hdf5")
                        dN = prev_len - len(parti_set) 
                        if dN != 0:
                            ghost_parti = Particles()
                            for key_ in prev_set.key:
                                if len(parti_set[parti_set.key == key_]) == 0: #Add missing particle
                                    ghost_parti.add_particle(prev_set[prev_set.key == key_]) 

                            if len(ghost_parti) == 1:
                                if ghost_parti.mass <= self.Star_min_mass:
                                    ejec_jmo[fiter] += 1
                                else:
                                    ejec_star[fiter] += 1
                            elif len(ghost_parti) == 2:
                                if dN == 2:
                                    kepler_elements = orbital_elements_from_binary(ghost_parti, G=constants.G)
                                    semimajor = kepler_elements[2]
                                    eccentric = kepler_elements[3]
                                    if semimajor < 1000 | units.au and eccentric < 1:
                                        if max(ghost_parti.mass) < self.Star_min_mass: #Account for runaway JMO mergers
                                            ejec_jj[fiter] += 1
                                        elif min(ghost_parti.mass) >= self.Star_min_mass:
                                            ejec_ss[fiter] += 1
                                        else:
                                            ejec_js[fiter] += 1

                                else:
                                    if max(ghost_parti.mass) < self.Star_min_mass:
                                        merge_jj[fiter] += 1
                                    elif min(ghost_parti.mass) >= self.Star_min_mass:
                                        merge_ss[fiter] += 1
                                    else:
                                        merge_js[fiter] += 1
                            else:
                                emerge_parti = Particles()
                                for key_ in parti_set.key:
                                    if len(prev_set[prev_set.key == key_]) == 0:  #Clump all mergers
                                        emerge_parti.add_particle(parti_set[parti_set.key == key_])
                                
                                bin_combo = list(itertools.combinations(ghost_parti, 2)) #All possible binaries in system
                                det_keys_g = [ ]
                                det_keys_e = [ ]
                                for bin_ in bin_combo:
                                    for emerge in emerge_parti:
                                        bin = Particles()  
                                        bin.add_particle(bin_[0])
                                        bin.add_particle(bin_[1])
                                        if abs(bin.mass.sum() - emerge.mass) <= (1 | units.MEarth):
                                            if max(bin.mass) < self.Star_min_mass:
                                                merge_jj[fiter] += 1
                                            elif min(bin.mass) >= self.Star_min_mass:
                                                merge_ss[fiter] += 1
                                            else:
                                                merge_js[fiter] += 1
                                            det_keys_g.append(bin_[0].key)
                                            det_keys_g.append(bin_[1].key)
                                            det_keys_e.append(emerge.key)
                                            break

                                for key_ in det_keys_e: #Remove all remnants who have their product found
                                    emerge_parti -= emerge_parti[emerge_parti.key == key_]
                                for key_ in det_keys_g: #Remove all products who have their merger found
                                    ghost_parti -= ghost_parti[ghost_parti.key == key_]

                                if len(ghost_parti) != 0:
                                    components = ghost_parti.connected_components(threshold = 
                                                                            (20000 |units.au))
                                    det_keys_g = [ ]
                                    det_keys_e = [ ]
                                    for c in components: #Looks at possibility for N > 2 merging scenarios
                                        for emerge in emerge_parti:
                                            if abs(c.mass.sum() - emerge.mass) <= (1 | units.MEarth):
                                                det_keys_g = np.concatenate((det_keys_g, c.key), axis=None)
                                                det_keys_e.append(emerge_parti.key)
                                                if len(c) == 2:
                                                    if max(c.mass) < self.Star_min_mass:
                                                        merge_jj[fiter] += 1
                                                    elif min(c.mass) >= self.Star_min_mass:
                                                        merge_ss[fiter] += 1
                                                    else:
                                                        merge_js[fiter] += 1
                                                else:       
                                                    if len(c[c.mass >= self.Star_min_mass]) == len(c):
                                                        merge_ss[fiter] += len(c) - 1
                                                    elif len(c[[c.mass < self.Star_min_mass]]) == len(c):
                                                        merge_jj[fiter] += len(c) - 1
                                                    else:
                                                        tkey = [ ]
                                                        for comp_ in c:
                                                            dist_temp = 10 | units.pc
                                                            for part_ in c[c.key != comp_.key]:
                                                                distance = (comp_.position - part_.position).length()
                                                                if distance < dist_temp:
                                                                    pkey = part_.key
                                                            tkey.append(sorted([comp_.key, pkey]))
                                                        
                                                        couple, counts = np.unique(tkey, return_counts=True)
                                                        p1 = ghost_parti[ghost_parti.key == couple[counts > 1][0]]
                                                        p2 = ghost_parti[ghost_parti.key == couple[counts > 1][1]]
                                                        if max(p1.mass, p2.mass) < self.Star_min_mass:
                                                            merge_jj[fiter] += 1
                                                        elif min(p1.mass, p2.mass) >= self.Star_min_mass:
                                                            merge_ss[fiter] += 1
                                                        else:
                                                            merge_js[fiter] += 1
                                
                                    for key_ in det_keys_e:
                                        emerge_parti -= emerge_parti[emerge_parti.key == key_]
                                    for key_ in det_keys_g:
                                        ghost_parti -= ghost_parti[ghost_parti.key == key_]
                                    
                                    if len(ghost_parti) != 0: #Look at ejection events for remainder
                                        components = ghost_parti.connected_components(threshold = (20000 |units.au))
                                        ejec_jmo[fiter] += len(ghost_parti[ghost_parti.mass < self.Star_min_mass])
                                        ejec_star[fiter] += len(ghost_parti[ghost_parti.mass > self.Star_min_mass])
                                        if len(components) != 0:
                                            for c in components:
                                                if len(c) > 1:
                                                    bin_sys = Particles()
                                                    bin_sys.add_particle(c[0])
                                                    bin_sys.add_particle(c[1])

                                                    kepler_elements = orbital_elements_from_binary(bin_sys, G=constants.G)
                                                    semimajor = kepler_elements[2]
                                                    eccentric = kepler_elements[3]

                                                    if semimajor < 1000 | units.au and eccentric < 1:
                                                        if max(bin.mass) < self.Star_min_mass:
                                                            ejec_jj[fiter] += 1
                                                            ejec_jmo[fiter] -= 2
                                                        elif min(bin.mass) >= self.Star_min_mass:
                                                            ejec_ss[fiter] += 1
                                                            ejec_star[fiter] -= 2
                                                        else:
                                                            ejec_js[fiter] += 1
                                                            ejec_jmo[fiter] -= 1
                                                            ejec_star[fiter] -= 1
                                                else:
                                                    if c.mass < self.Star_min_mass:
                                                        ejec_jmo[fiter] += 1
                                                    elif c.mass >= self.Star_min_mass:
                                                        ejec_star[fiter] += 1
                        prev_set = read_set_from_file(dt, "hdf5")
                    fiter += 1

            else:
                event_files = natsort.natsorted(glob.glob(os.path.join(str(path+"/event_data/")+"*")))
                for config_ in traj_files:
                    events = pd.read_hdf(event_files[fiter], 'dfEvent')
                    print("Processing statistics for ", config_)
                    tracked_key = np.asarray([ ])
                    for col_ in range(np.shape(events)[1]):
                        event_key = events.iloc[1][col_]
                        event_type = events.iloc[2][col_]
                        event_mass = events.iloc[3][col_]
                        proc = True
                        if "Merger" in event_type:
                            if len(tracked_key[tracked_key == event_key[0]]) != 0:
                                proc = False
                        else:
                            if len(tracked_key[tracked_key == event_key]) != 0:
                                proc = False
                        
                        if (proc):
                            if "Merger" in event_type:
                                if max(event_mass) < self.Star_min_mass:
                                    merge_jj[fiter] += 1
                                elif min(event_mass) > self.Star_min_mass:
                                    merge_ss[fiter] += 1
                                else:
                                    merge_js[fiter] += 1
                            else:
                                if event_mass < self.Star_min_mass:
                                    ejec_jmo[fiter] += 1
                                else:
                                    ejec_star[fiter] += 1
                            tracked_key = np.asarray(np.concatenate((tracked_key, 
                                                                     event_key), 
                                                                     axis=None))

                        """if len(events.iloc[0][events.iloc[0] == event_time]) > 1 and "Merger" not in event_type:
                            for key_ in events.iloc[1][events.iloc[0] == event_time]:
                                if key_ not in tracked_key:
                                    masses = events.iloc[3][(events.iloc[0] == event_time) & (events.iloc[2] == event_type)]
                                    dt = int(np.floor(event_time.value_in(units.kyr)/10))
                                    parti_set = read_set_from_file(sim_snapshot[dt], "hdf5")
                                    for mass_ in masses[masses != event_mass]:
                                        components = parti_set.connected_components(threshold = 
                                                                                    2000 |units.au)
                                        for c in components:
                                            if len(c) > 1:
                                                bin_combo = list(itertools.combinations(c, 2)) #All possible binaries in system
                                                for bin_ in bin_combo:
                                                    bin_sys = Particles()
                                                    bin_sys.add_particle(bin_[0])
                                                    bin_sys.add_particle(bin_[1])

                                                    kepler_elements = orbital_elements_from_binary(bin_sys, G=constants.G)
                                                    semimajor = kepler_elements[2]
                                                    eccentric = kepler_elements[3]

                                                    if semimajor < 1000 | units.au and eccentric < 1:
                                                        if max(bin_sys.mass) < self.Star_min_mass:
                                                            ejec_jj[fiter] += 1
                                                            ejec_jmo[fiter] -= 2
                                                        elif min(bin_sys.mass) >= self.Star_min_mass:
                                                            ejec_ss[fiter] += 1
                                                            ejec_star[fiter] -= 2
                                                        else:
                                                            ejec_js[fiter] += 1
                                                            ejec_star[fiter] -= 1
                                                            ejec_jmo[fiter] -= 1
                                                        
                                                        tracked_key = np.asarray(np.concatenate((tracked_key, 
                                                                                 [bin_[0].key, bin_[1].key]), 
                                                                                 axis=None))"""

                    fiter += 1

            complete_data = [merge_jj, merge_js, merge_ss, ejec_jmo, 
                             ejec_star, ejec_jj, ejec_js, ejec_ss]
            complete_string = ["JJ Merge", "JS Merge", "SS Merge", 
                               "JMO Ejec", "Star Ejec", "JuMBO Ejec", 
                               "JS Ejec", "SS Ejec"]

            liter = 0
            with open(os.path.join("plotters/figures/system_evolution/outputs/", 
                      tname), 'w') as f:
                for val_, prop_ in zip(complete_data, complete_string):
                    liter += 1
                    print("Property: ", prop_, val_)
                    if len(val_) > 0:
                        median = np.median(val_)
                        IQR_low, IQR_high = np.percentile(val_, [25, 75])
                        q1 = median - IQR_low
                        q3 = IQR_high - median
                    
                        lines = ["Raw {} Data: {}".format(prop_, val_),
                                "Median: {}".format(median),
                                "IQR Low: {}, IQR High: {}".format(q1, q3)]
                        for line_ in lines:
                            f.write(line_+'\n')
                        if liter == 4:
                            f.write('\n')
                        f.write('\n')
                    
    def population_statistics(self, model_iter, dt_crop):
        """Extract population statistics of systems present"""
        
        directory = "data/Simulation_Data/"+str(self.models[model_iter])+"/"
        dir_configs = natsort.natsorted(glob.glob(os.path.join(str(directory)+"simulation_snapshot/**")))
        Nsims = len(dir_configs)

        if (dt_crop):
            fname = "plotters/figures/system_evolution/"+str(self.models[model_iter])+"_sem_ecc_mixed_systs_crop.pdf"
            tname = self.models[model_iter]+'_syst_info_crop.txt'
        else:
            fname = "plotters/figures/system_evolution/"+str(self.models[model_iter])+"_sem_ecc_mixed_systs.pdf"
            tname = self.models[model_iter]+'_syst_info.txt'

        flatten_mkeys = self.plot_func.flatten_arr(self.fin_mkeys, model_iter)
        done_keys = flatten_mkeys

        #Process N > 2 Systems
        mpops = [ ]
        mtype = [ ]
        run_temp = [ ]
        mpops_run = [[ ] for i in range(Nsims+1)]
        for system in range(len(self.mkeys[model_iter])):
            temp_key = [ ]
            temp_mtype = [ ]
            for key_ in range(len(self.mkeys[model_iter][system])):
                key = self.mkeys[model_iter][system][key_]
                if key not in temp_key:
                    temp_key.append(key)
                    mass = self.fin_mmass[model_iter][system][0][key_]
                    if mass <= self.JuMBO_max_mass:
                        temp_mtype.append("JuMBOs")
                    else:
                        temp_mtype.append("star")
            if len(temp_key) > 2:
                mpops.append(len(temp_key))
                mtype.append(sorted(temp_mtype))
                run_idx = self.fin_msim_iter[model_iter][system] - 1
                run_temp.append(run_idx)
                run_idx = len(np.unique(run_temp)) - 1
                mpops_run[run_idx].append(len(temp_key))

        upop = [ ]
        cpop = [ ]
        for pops in mpops_run:
            unique_pop_run, pop_counts_run = np.unique(pops, return_counts=True)
            upop.append(unique_pop_run)
            cpop.append(pop_counts_run)

        flat_upop = [ ]
        flat_cpop = [ ]
        for system, counts in zip(upop, cpop):
            flat_upop = np.concatenate((flat_upop, system), axis=None)
            flat_cpop = np.concatenate((flat_cpop, counts), axis=None)

        pops = [ ]
        med_cpop = [ ]
        IQR_cpop = [ ]
        for val_ in np.unique(flat_upop):
            pops.append(val_)
            med_cpop.append(np.median(flat_cpop[flat_upop == val_]))
            q11, q31 = np.percentile(np.median(flat_cpop[flat_upop == val_]), [25, 75])
            q11 = np.median(flat_cpop[flat_upop == val_]) - q11
            q31 -= np.median(flat_cpop[flat_upop == val_])
            IQR_cpop.append([q11, q31])
        
        unique_pops, pops_counts = np.unique(mpops, return_counts=True)
        if len(mtype)>1:
            unique_type, type_counts = np.unique(mtype, return_counts=True)
            type_counts = [i/Nsims for i in type_counts]
        else:
            unique_type = mtype
            type_counts = [[ ] for i in mtype]
        pops_counts = [i/Nsims for i in pops_counts]

        red_array = [ ]
        for type_, count_ in zip(unique_type, type_counts):
            red_array.append([type_, count_])

        jmb_bin = np.zeros(Nsims)
        str_bin = np.zeros(Nsims)
        mix_bin = np.zeros(Nsims)
        mix_sem = [ ]
        mix_ecc = [ ]
        mix_masses = [ ]

        final_jmb_keys = [ ]
        final_jmb_runs = [ ]

        #Processing N = 2 Systems (Only Preliminary Filtered JuMBO)
        for jmb_idx_ in range(len(self.jmb_idx[model_iter])):
            proc = True
            for key_ in self.fin_bkeys[model_iter][jmb_idx_]:
                if key_ in done_keys:
                    proc = False
                else:
                    done_keys = np.concatenate((done_keys, key_), axis = None)
            if (proc):
                masses = self.fin_bmass[model_iter][jmb_idx_]
                run_idx = self.fin_bsim_iter[model_iter][jmb_idx_]-1
                if max(masses)*(1 | units.MSun) <= self.JuMBO_max_mass:
                    jmb_bin[run_idx] += 1
                    final_jmb_keys.append(sorted(self.fin_bkeys[model_iter][jmb_idx_]))
                    final_jmb_runs.append(run_idx)
                elif min(masses)*(1 | units.MSun) > self.JuMBO_max_mass:
                    str_bin[run_idx] += 1
                else:
                    mix_bin[run_idx] += 1
                    mix_sem.append(self.fin_bsemi[model_iter][jmb_idx_])
                    mix_ecc.append(self.fin_becce[model_iter][jmb_idx_])
                    mix_masses.append(max(masses))

        #Secondary JuMBO Filters
        for syst_ in range(len(self.fin_becce[model_iter])):
            proc = True
            for key_ in self.fin_bkeys[model_iter][syst_]:
                if key_ in done_keys:
                    proc = False
                else:
                    done_keys = np.concatenate((done_keys, key_), axis = None)
            if (proc):
                masses = self.fin_bmass[model_iter][syst_]
                run_idx = self.fin_bsim_iter[model_iter][syst_]-1
                if max(masses * (1 | units.MSun)) <= self.JuMBO_max_mass:
                    jmb_bin[run_idx] += 1
                    final_jmb_keys.append(sorted(self.fin_bkeys[model_iter][syst_]))
                    final_jmb_runs.append(run_idx)
                elif min(masses * (1 | units.MSun)) > self.JuMBO_max_mass:
                    str_bin[run_idx] += 1
                else:
                    mix_bin[run_idx] += 1
                    mix_sem.append(self.fin_bsemi[model_iter][syst_])
                    mix_ecc.append(self.fin_becce[model_iter][syst_])
                    mix_masses.append(max(masses))

        #Compute persistent binaries
        same_init = np.zeros(Nsims)
        for fkey_, run_idx in zip(final_jmb_keys, final_jmb_runs):
            for ikey_ in self.init_bkeys[model_iter]:
                sorted_ikey = sorted(ikey_)
                fkey1 = str(fkey_[0])[0] + "." + str(fkey_[0])[1:]
                fkey2 = str(fkey_[1])[0] + "." + str(fkey_[1])[1:]

                if fkey1[:14] == str(sorted_ikey[0])[:14] \
                    and fkey2[:14] == str(sorted_ikey[1])[:14]:
                    same_init[run_idx] += 1
                    break

        new_jmb = jmb_bin - same_init      
        with open(os.path.join("plotters/figures/system_evolution/outputs/", 
                  tname), 'w') as f:
            complete_data = [jmb_bin, new_jmb, str_bin, mix_bin, 
                            mix_sem, mix_ecc, mix_masses, unique_pops]
            complete_string = ["JuMBO Systs", "New JuMBOs", "SS Systs", 
                               "JS Systs", "JS Semi-major", "JS Eccentricity", 
                               "JS masses", "N>2 Populations"]
            for val_, prop_ in zip(complete_data, complete_string):
                if len(val_) > 0:
                    median = np.median(val_)
                    IQR_low, IQR_high = np.percentile(val_, [25,75])
                    q1 = median - IQR_low
                    q3 = IQR_high - median
                else:
                    median = None
                    q1 = None
                    q3 = None

                lines = ["Raw data {}: {}".format(prop_, val_),
                         "Median: {}".format(median),
                         "IQR Low: {}, IQR High: {}".format(q1, q3)]
                if prop_ == "JuMBO Systs":
                    init_J = len(self.init_bkeys[model_iter])/Nsims
                    surv_rateJmed = median/init_J
                    lines.append("Med. Survival Rate: {}".format(surv_rateJmed))

                for line_ in lines:
                    f.write(line_+'\n')
                f.write("\n")
            f.write("\n===============================\nSystem Type & Avg. Detections:")
            for arr_ in red_array:
                f.write(str(arr_)+'\n')

        xx, yy = np.mgrid[0:1000:200j, 0:1:200j]
        positions = np.vstack([xx.ravel(), yy.ravel()])
        values = np.vstack([mix_sem, mix_ecc])
        if len(mix_sem) > 5:
            kernel = stats.gaussian_kde(values, bw_method = "silverman")
            f = np.reshape(kernel(positions).T, xx.shape)

            fig, ax = plt.subplots()
            cfset = ax.contourf(xx, yy, f, cmap="Blues", levels = 7, zorder = 1)
            cset = ax.contour(xx, yy, f, colors = "black", levels = 7, zorder = 2)
            ax.clabel(cset, inline=1, fontsize=10)
            ax.set_xlabel(r"$a$ [au]", fontsize = self.clean_plot.axlabel_size)
            ax.set_ylabel(r"$e$", fontsize = self.clean_plot.axlabel_size)
            ax.set_xlim(0,1000)
            ax.set_ylim(0,1)
            self.clean_plot.tickers(ax, "hist")
            fig.savefig(fname, dpi=700, bbox_inches='tight')
            plt.clf()
        
    def mass_CDF(self, model_choices, dt_crop):
        """Function to plot the CDF of primary masses"""

        file = os.path.join("data/observations/src/obs_data.txt")
        q_obs = [ ]
        mprim_obs = [ ]
        with open(file, 'r', newline ='') as file:
            csv_reader = csv.reader(file)
            for row_ in csv_reader:
                mass1 = float(row_[3]) * (1 | units.MSun)/(1 | units.MJupiter)
                mass2 = float(row_[5]) * (1 | units.MSun)/(1 | units.MJupiter)
                q = min(mass1/mass2, mass2/mass1)
                q_obs.append(q)
                max_mass = max(mass1, mass2)
                mprim_obs.append(max_mass)
        obs_sort = np.sort(mprim_obs)
        obs_iter = np.asarray([i for i in enumerate(obs_sort)])
        obs_iter = obs_iter[:,0]
        obs_iter /= max(obs_iter)

        fig, ax = plt.subplots()
        ax.plot(obs_sort, obs_iter, color = "black", 
                linewidth = 2.5, label = "Observation")
        labels, plot_labels, extra_str = self.clean_plot.model_layout(model_choices)
        if (dt_crop):
            fname = "plotters/figures/binary_evolution/"+str(extra_str)+"mprim_vs_obs_crop.pdf"
        else:
            fname = "plotters/figures/binary_evolution/"+str(extra_str)+"mprim_vs_obs.pdf"

        iloop = 0
        for model_iter in model_choices:
            print("Processing ", self.models[model_iter])
            flatten_mkeys = self.plot_func.flatten_arr(self.fin_mkeys, model_iter)
            done_keys = flatten_mkeys

            if self.models[model_iter] != "Fractal_rvir0.5_FFOnly":
                init_pmass = [ ]
                for syst_ in range(len(self.init_bmass[model_iter])):
                    mass1 = self.init_bmass[model_iter][syst_][0] * (1 | units.MSun)/(1 | units.MJupiter)
                    mass2 = self.init_bmass[model_iter][syst_][1] * (1 | units.MSun)/(1 | units.MJupiter)
                    max_mass = max(mass1, mass2)
                    init_pmass.append(max_mass)
                simi_sort = np.sort(init_pmass)
                simi_iter = np.asarray([i for i in enumerate(simi_sort)])
                simi_iter = simi_iter[:,0]
                simi_iter /= max(simi_iter)
                if model_iter == model_choices[0]:
                    ax.plot(simi_sort, simi_iter, label = "Initial", 
                            color = self.colours[iloop])
                else:
                    ax.plot(simi_sort, simi_iter, color = self.colours[iloop])

            fin_pmass = [ ]
            for idx_ in self.jmb_idx[model_iter]:
                proc = True
                for key_ in self.fin_bkeys[model_iter][idx_]:
                    if key_ in done_keys:
                        proc = False
                    else:
                        done_keys = np.concatenate((done_keys, key_), axis = None)
                if (proc):
                    mass1 = self.fin_bmass[model_iter][idx_][0] * (1 | units.MSun)/(1 | units.MJupiter)
                    mass2 = self.fin_bmass[model_iter][idx_][1] * (1 | units.MSun)/(1 | units.MJupiter)
                    if (max(mass1, mass2)) <= self.JuMBO_max_mass.value_in(units.MJupiter):
                        fin_pmass.append(max(mass1, mass2))
                
            for syst_ in range(len(self.fin_bkeys[model_iter])):
                proc = True
                for key_ in self.fin_bkeys[model_iter][syst_]:
                    if key_ in done_keys:
                        proc = False
                    else:
                        done_keys = np.concatenate((done_keys, key_), axis = None)
                if (proc):
                    mass1 = self.fin_bmass[model_iter][syst_][0] * (1 | units.MSun)/(1 | units.MJupiter)
                    mass2 = self.fin_bmass[model_iter][syst_][1] * (1 | units.MSun)/(1 | units.MJupiter)
                    if (max(mass1, mass2)) <= self.JuMBO_max_mass.value_in(units.MJupiter):
                        fin_pmass.append(max(mass1, mass2))
                        
            simf_sort = np.sort(fin_pmass)
            simf_iter = np.asarray([i for i in enumerate(simf_sort)])
            simf_iter = simf_iter[:,0]
            simf_iter /= max(simf_iter)
            ax.plot(simf_sort, simf_iter, label = labels[iloop], 
                    color = self.colours[iloop], linestyle = "-.")
            iloop += 1
            
        self.clean_plot.tickers(ax, "plot")
        ax.set_xlabel(r"$M_{\mathrm{prim}} [\mathrm{M}_{\mathrm{Jup}}]$", 
                      fontsize = self.clean_plot.axlabel_size)
        ax.set_ylabel(r"CDF", fontsize = self.clean_plot.axlabel_size)
        ax.set_ylim(0, 1.02)
        ax.legend(prop={'size': self.clean_plot.axlabel_size})
        fig.savefig(fname, dpi=700, bbox_inches='tight')
        plt.clf()

    def mass_params(self, model_iter, dt_crop):
        """Plot the evolution of the mass ratio and raw masses"""

        if (dt_crop):
            fname = "plotters/figures/binary_evolution/"+str(self.models[model_iter])+"_mass_distr_crop.pdf"
            tname = self.models[model_iter]+'mass_params_crop.txt'
        else:
            fname = "plotters/figures/binary_evolution/"+str(self.models[model_iter])+"_mass_distr.pdf"
            tname = self.models[model_iter]+'mass_params.txt'

        file = os.path.join("data/observations/src/obs_data.txt")
        q_obs = [ ]
        mprim_obs = [ ]
        with open(file, 'r', newline ='') as file:
            csv_reader = csv.reader(file)
            for row_ in csv_reader:
                mass1 = float(row_[3]) * (1 | units.MSun)/(1 | units.MJupiter)
                mass2 = float(row_[5]) * (1 | units.MSun)/(1 | units.MJupiter)
                q = min(mass1/mass2, mass2/mass1)
                q_obs.append(q)
                max_mass = max(mass1, mass2)
                mprim_obs.append(max_mass)

        flatten_mkeys = self.plot_func.flatten_arr(self.fin_mkeys, model_iter)
        done_keys = flatten_mkeys

        fin_q = [ ]
        fin_pmass = [ ]
        fin_smass = [ ]
        for idx_ in self.jmb_idx[model_iter]:
            proc = True
            for key_ in self.fin_bkeys[model_iter][idx_]:
                if key_ in done_keys:
                    proc = False
                else:
                    done_keys = np.concatenate((done_keys, key_), axis=None)
            if (proc):
                mass1 = self.fin_bmass[model_iter][idx_][0] * (1 | units.MSun)/(1 | units.MJupiter)
                mass2 = self.fin_bmass[model_iter][idx_][1] * (1 | units.MSun)/(1 | units.MJupiter)
                if (max(mass1, mass2)) <= self.JuMBO_max_mass.value_in(units.MJupiter):
                    q = min(mass1/mass2, mass2/mass1)
                    fin_q.append(q)
                    fin_pmass.append(max(mass1, mass2))
                    fin_smass.append(min(mass1, mass2))
            
        for syst_ in range(len(self.fin_bkeys[model_iter])):
            proc = True
            for key_ in self.fin_bkeys[model_iter][syst_]:
                if key_ in done_keys:
                    proc = False
                else:
                    done_keys = np.concatenate((done_keys, key_), axis = None)
            if (proc):
                mass1 = self.fin_bmass[model_iter][syst_][0] * (1 | units.MSun)/(1 | units.MJupiter)
                mass2 = self.fin_bmass[model_iter][syst_][1] * (1 | units.MSun)/(1 | units.MJupiter)
                if (max(mass1, mass2)) <= self.JuMBO_max_mass.value_in(units.MJupiter):
                    q = min(mass1/mass2, mass2/mass1)
                    fin_q.append(q)
                    fin_pmass.append(max(mass1, mass2))
                    fin_smass.append(min(mass1, mass2))
            
        xlims = [0, 14]
        ylims = [0, 1.001]

        xx, yy = np.mgrid[xlims[0]:xlims[1]:200j, ylims[0]:ylims[1]:200j]
        positions = np.vstack([xx.ravel(), yy.ravel()])
        values = np.vstack([fin_pmass, fin_q])
        kernel = stats.gaussian_kde(values, bw_method = "silverman")
        f = np.reshape(kernel(positions).T, xx.shape)

        fig, ax = plt.subplots()
        cfset = ax.contourf(xx, yy, f, cmap="Blues", levels = 7, zorder = 1)
        cset = ax.contour(xx, yy, f, colors = "k", levels = 7, zorder = 2)
        ax.clabel(cset, inline=1, fontsize=10)
        ax.set_xlabel(r'$M_{\mathrm{prim}} [M_{\mathrm{Jup}}]$', 
                      fontsize = self.clean_plot.axlabel_size)
        ax.set_ylabel(r'$q$', fontsize = self.clean_plot.axlabel_size)
        ax.scatter(mprim_obs, q_obs, marker = "X", color = "red", 
                   edgecolor = "black", zorder = 3)
        self.clean_plot.tickers(ax, "hist")
        ax.set_xlim(0, 14)
        ax.set_ylim(0, 1.001)
        fig.savefig(fname, dpi=700, bbox_inches='tight')
        plt.close() 

        q11, q31 = np.percentile(fin_q, [25, 75])
        q11 = np.median(fin_q) - q11
        q31 -= np.median(fin_q)

        q12, q32 = np.percentile(fin_pmass, [25, 75])
        q12 = np.median(fin_pmass) - q12
        q32 -= np.median(fin_pmass)

        q12s, q32s = np.percentile(fin_smass, [25, 75])
        q12s = np.median(fin_smass) - q12s
        q32s -= np.median(fin_smass)
        lines = ["Model "+str(self.models[model_iter]),
                 "Median q: "+str(np.median(fin_q)),
                 "IQR: "+str(q11)+" "+str(q31),
                 "Median Mprim: "+str(np.median(fin_pmass)),
                 "IQR Mprim: "+str(q12)+" "+str(q32),
                 "Median Msec: "+str(np.median(fin_smass)),
                 "IQR Msec: "+str(q12s)+" "+str(q32s)]
        with open(os.path.join("plotters/figures/binary_evolution/outputs/", 
                  tname), 'w') as f:
            for line_ in lines:
                f.write(line_+'\n')

    def nearest_neigh_process(self):
        """Function to track the demographic and statistics of bound systems
           Also plots the distribution for nearest neighbours of each particle
        """

        for model_iter in range(len(self.models[:3])-1):
            nn_star = [ ]
            nn_jmb = [ ]
            semi_jmb = [ ]
            proj_sep_jmb = [ ] 
            nn_FF = [ ]

            print("Processing ", self.models[model_iter])
            directory = "data/Simulation_Data/"+str(self.models[model_iter])+"/"
            dir_configs = natsort.natsorted(glob.glob(os.path.join(str(directory)+"simulation_snapshot/**")))
            nmulti_jmb = 0 
            nmulti = 0
            
            flatten_mkeys = self.plot_func.flatten_arr(self.fin_mkeys, model_iter)
            flatten_keys = self.plot_func.flatten_arr(self.fin_bkeys, model_iter)

            for config_ in dir_configs:
                if len(glob.glob(config_+"/*")) > 100:
                    file_name = natsort.natsorted(glob.glob(config_+"/*"))[-1]
                    print("Reading Data: ", file_name)
                    data = read_set_from_file(file_name, "hdf5")
                    temp_keys = [ ]
                    for parti_ in data:
                        if parti_.key in flatten_mkeys:
                            for system in self.fin_mkeys[model_iter]:
                                if parti_.key in system and parti_.key not in temp_keys:
                                    syst = Particles()
                                    for key_ in system:
                                        if key_ not in temp_keys:
                                            temp_keys.append(key_)
                                            syst.add_particle(data[data.key == key_])
                                    com_pos = syst.center_of_mass()
                                    sim_syst = data - syst
                                    if max(syst.mass) <= (0.013 | units.MSun):
                                        nn_dist = min(np.sqrt((com_pos.x - sim_syst.x)**2 + 
                                                              (com_pos.y - sim_syst.y)**2 + 
                                                              (com_pos.z - sim_syst.z)**2))
                                        nn_jmb.append(np.log10(nn_dist.value_in(units.au)))
                                        nmulti_jmb += 1

                                        proj_sep = np.sqrt((syst[0].x - syst[1].x)**2 
                                                          +(syst[0].y - syst[1].y)**2)
                                        proj_sep_jmb.append(proj_sep.value_in(units.au))

                                        bin_sys = Particles()
                                        bin_sys.add_particle(syst[0])
                                        bin_sys.add_particle(syst[1])
                                        kepler_elements = orbital_elements_from_binary(bin_sys, G=constants.G)
                                        semi_jmb.append(kepler_elements[2].value_in(units.au))

                                    else:
                                        nn_dist = min(np.sqrt((com_pos.x - sim_syst.x)**2 + 
                                                              (com_pos.y - sim_syst.y)**2 + 
                                                              (com_pos.z - sim_syst.z)**2))
                                        nn_star.append(np.log10(nn_dist.value_in(units.au)))
                                        nmulti += 1

                        elif parti_.key in flatten_keys:
                            for system in self.fin_bkeys[model_iter]:
                                if parti_.key in system and parti_.key not in temp_keys:
                                    for key_ in system:
                                        temp_keys.append(key_)
                                    body1 = data[(data.key == system[0])]
                                    body2 = data[(data.key == system[1])]
                                    if body1.mass <= (0.013 | units.MSun) \
                                        and body2.mass <= (0.013 | units.MSun):
                                        syst = Particles()
                                        syst.add_particle(body1)
                                        syst.add_particle(body2)
                                        com_pos = syst.center_of_mass()
                                        sim_syst = data - body1 - body2
                                        
                                        nn_dist = min(np.sqrt((com_pos.x - sim_syst.x)**2 + 
                                                              (com_pos.y - sim_syst.y)**2 + 
                                                              (com_pos.z - sim_syst.z)**2))
                                        nn_jmb.append(np.log10(nn_dist.value_in(units.au)))

                                        proj_sep = np.sqrt((syst[0].x - syst[1].x)**2 
                                                          +(syst[0].y - syst[1].y)**2)
                                        proj_sep_jmb.append(proj_sep.value_in(units.au))

                                        bin_sys = Particles()
                                        bin_sys.add_particle(syst[0])
                                        bin_sys.add_particle(syst[1])
                                        kepler_elements = orbital_elements_from_binary(bin_sys, G=constants.G)
                                        semi_jmb.append(kepler_elements[2].value_in(units.au))
                                
                        else:
                            if parti_.mass < (0.013 | units.MSun):
                                sim_syst = data - parti_
                                nn_dist = min(np.sqrt((parti_.x - sim_syst.x)**2 + 
                                                      (parti_.y - sim_syst.y)**2 + 
                                                      (parti_.z - sim_syst.z)**2))
                                nn_FF.append(np.log10(nn_dist.value_in(units.au)))

                            else:
                                sim_syst = data - parti_
                                nn_dist = min(np.sqrt((parti_.x - sim_syst.x)**2 + 
                                                      (parti_.y - sim_syst.y)**2 + 
                                                      (parti_.z - sim_syst.z)**2))
                                nn_star.append(np.log10(nn_dist.value_in(units.au)))
            file_name = "nn_output_"+self.models[model_iter]+".h5"
            array = pd.DataFrame([nn_jmb, semi_jmb, proj_sep_jmb, nn_star, nn_FF])
            array.to_hdf(os.path.join(file_name), key ='df', mode = 'w')

    def nn_plotter(self, model):
        
        for model_iter in model:
            data = pd.read_hdf('nn_output_'+self.models[model_iter]+'.h5', 'df')
            
            nn_jmb = data.iloc[0]
            nn_jmb = nn_jmb[~np.isnan(nn_jmb)]
            nn_star = data.iloc[3]
            nn_star = nn_star[~np.isnan(nn_star)]
            nn_FF = data.iloc[4]
            nn_FF = nn_FF[~np.isnan(nn_FF)]

            kde_jmb = sm.nonparametric.KDEUnivariate((nn_jmb))
            kde_jmb.fit()
            kde_jmb.density = (kde_jmb.density/max(kde_jmb.density))

            kde_star = sm.nonparametric.KDEUnivariate((nn_star))
            kde_star.fit()
            kde_star.density = (kde_star.density/max(kde_star.density))

            kde_FF = sm.nonparametric.KDEUnivariate((nn_FF))
            kde_FF.fit()
            kde_FF.density = (kde_FF.density/max(kde_FF.density))

            fig, ax = plt.subplots()
            ax.plot(kde_star.support, kde_star.density, 
                    color = "black", label = "Star")
            ax.fill_between(kde_star.support, kde_star.density, 
                            alpha = 0.33, color = "black")
            ax.plot(kde_jmb.support, kde_jmb.density, 
                    color = "blue", label = "JuMBO")
            ax.fill_between(kde_jmb.support, kde_jmb.density, 
                            alpha = 0.33, color = "blue")
            ax.plot(kde_FF.support, kde_FF.density, 
                    color = "red", label = "FF")
            ax.fill_between(kde_FF.support, kde_FF.density, 
                            alpha = 0.33, color = "red")

            ax.legend(prop={'size': self.clean_plot.axlabel_size})
            self.clean_plot.tickers(ax, "plot")
            ax.set_xlabel(r"$\log_{10}r_{\mathrm{NN}}$ [au]", 
                          fontsize = self.clean_plot.axlabel_size)
            ax.set_ylabel(r"$\log_{10}(\rho/\rho_{\mathrm{max}})$", 
                          fontsize = self.clean_plot.axlabel_size)
            ax.set_ylim(0, 1.05)
            plt.savefig("plotters/figures/binary_evolution/"+str(self.models[model_iter])+"_nn_hist.pdf", 
                        dpi=700, bbox_inches='tight')
            plt.close()
            
    def orb_cdf_plot(self, model_choices, dt_crop):
        """Plot the evolution of JuMBO orbital parameters"""

        print("Plotting ", self.models[model_choices[0]], " and ", self.models[model_choices[-1]])
        file = os.path.join("data/observations/src/obs_data.txt")

        proj_sep = [ ]
        mprim_obs = [ ]
        with open(file, 'r', newline ='') as file:
            csv_reader = csv.reader(file)
            for row_ in csv_reader:
                proj_sep.append(float(row_[7]))
                mass1 = float(row_[3]) * (1 | units.MSun)/(1 | units.MJupiter)
                mass2 = float(row_[5]) * (1 | units.MSun)/(1 | units.MJupiter)
                max_mass = max(mass1, mass2)
                mprim_obs.append(max_mass)

        obs_sort = np.sort(proj_sep)
        obs_iter = np.asarray([i for i in enumerate(obs_sort)])
        obs_iter = obs_iter[:,0]
        obs_iter /= max(obs_iter)

        xlabel = [r'$e$', r'$a$ [au]', r'$r_{ij}$ [au]']
        file_name = ["eccentricity", "sem_axis", "proj_sep"]
        leg_txt, plot_label, extra_str = self.clean_plot.model_layout(model_choices)

        citer = 0            
        iparams = [[ ], [ ], [ ], [ ]]
        fparams = [[ ], [ ], [ ], [ ]]
        for config_ in model_choices:
            iparams[citer] = [self.init_becce[config_], 
                              self.init_bsemi[config_],
                              [None], 
                              self.init_bincl[config_]]
            fparams[citer] = [self.fin_becce[config_], 
                              self.fin_bsemi[config_],
                              self.fin_bproj[config_], 
                              self.fin_bincl[config_]]
            citer += 1

        for plot_ in range(len(file_name)):
            direc = "plotters/figures/binary_evolution/"
            if len(model_choices) > 2:
                mcomp = 2
            else:
                mcomp = 1

            if (dt_crop):
                fname = str(direc)+str(extra_str)+str(file_name[plot_])+"_crop.pdf"
                tname = str(self.models[model_choices[0]])+'_'+str(file_name[plot_])\
                        +'_'+self.models[model_choices[mcomp]]+'_crop.txt'
            else:
                fname = str(direc)+str(extra_str)+str(file_name[plot_])+".pdf"
                tname = str(self.models[model_choices[0]])+'_'+str(file_name[plot_])\
                        +'_'+self.models[model_choices[mcomp]]+'.txt'

            fig, ax = plt.subplots()  
            ax.set_xlabel(xlabel[plot_], fontsize = self.clean_plot.axlabel_size)
            ax.set_ylabel(r"CDF", fontsize = self.clean_plot.axlabel_size)
            if file_name[plot_] == "sem_axis":
                ax.plot(obs_sort, obs_iter, color = "black", 
                        linewidth = 2.5, label = "Observation")
            if file_name[plot_] != "proj_sep" and file_name[plot_] != "sem_axis":
                direc = "plotters/figures/binary_evolution/backup/"

            fvals = [[ ], [ ], [ ], [ ]]
            nsamp = [[ ], [ ], [ ], [ ]]
            ncrop = [[ ], [ ], [ ], [ ]]
            for k_ in range(len(model_choices)):
                flatten_mkeys = self.plot_func.flatten_arr(self.fin_mkeys, model_choices[k_])
                done_keys = flatten_mkeys
                isort, ipop = self.plot_func.cdf_plotter(iparams[k_][plot_])
                fin_jmb = [ ]
                fin_jmb_sem = [ ]
                fin_jmb_cropped = [ ]
                
                dir = "data/Simulation_Data/"+str(self.models[model_choices[k_]])+"/"
                Nsims = (len(glob.glob(os.path.join(dir+"simulation_snapshot/*"))))
                nsamples = np.zeros(Nsims)
                nsamples_cropped = np.zeros(Nsims)
                for idx_ in self.jmb_idx[model_choices[k_]]:
                    proc = True
                    for key_ in self.fin_bkeys[model_choices[k_]][idx_]:
                        if key_ in done_keys:
                            proc = False
                        else:
                            done_keys = np.concatenate((done_keys, key_), axis=None)
                    if (proc):
                        masses = self.fin_bmass[model_choices[k_]][idx_]
                        if max(masses * (1 | units.MSun)) <= self.JuMBO_max_mass:
                            fin_jmb.append(fparams[k_][plot_][idx_])
                            if file_name[plot_] == "proj_sep" or file_name[plot_] == "sem_axis":
                                fin_jmb_sem.append(fparams[k_][plot_-1][idx_])
                                run_idx = self.fin_bsim_iter[model_choices[k_]][idx_] - 1
                                nsamples[run_idx] += 1
                                if fparams[k_][plot_][idx_] > 25:
                                    fin_jmb_cropped.append(fparams[k_][plot_][idx_])
                                    nsamples_cropped[run_idx] += 1

                for syst_ in range(len(self.fin_bkeys[model_choices[k_]])):
                    proc = True
                    for key_ in self.fin_bkeys[model_choices[k_]][syst_]:
                        if key_ in done_keys:
                            proc = False
                        else:
                            done_keys = np.concatenate((done_keys, key_), axis = None)
                    if (proc):
                        masses = self.fin_bmass[model_choices[k_]][syst_]
                        if max(masses * (1 | units.MSun)) <= self.JuMBO_max_mass:
                            fin_jmb.append(fparams[k_][plot_][syst_])
                            if file_name[plot_] == "proj_sep" or file_name[plot_] == "sem_axis":
                                fin_jmb_sem.append(fparams[k_][plot_[1]][idx_])
                                run_idx = self.fin_bsim_iter[model_choices[k_]][idx_] - 1
                                nsamples[run_idx] += 1
                                if fparams[k_][plot_][idx_] > 25:
                                    fin_jmb_cropped.append(fparams[k_][plot_][idx_])
                                    nsamples_cropped[run_idx] += 1
                        
                fsort, fpop = self.plot_func.cdf_plotter(fin_jmb)
                if file_name[plot_] == "proj_sep":
                    fsorta, fpopa = self.plot_func.cdf_plotter(fin_jmb_sem)
                    if k_ == 0:
                        ax.plot(fsorta, fpopa, color = self.colours[k_], label = r"$a$", alpha = 0.5, zorder = 0)
                    else:
                        ax.plot(fsorta, fpopa, color = self.colours[k_], alpha = 0.5)

                if k_ == 0 and file_name[plot_] != "proj_sep":
                    ax.plot(isort, ipop, color = self.colours[k_], 
                            label = "Initial")
                else:
                    ax.plot(isort, ipop, color = self.colours[k_])

                ax.plot(fsort, fpop, label = "Final, "+leg_txt[k_], 
                        color = self.colours[k_], linestyle = "-.")
                if file_name[plot_] == "sem_axis":
                    fsortc, fpopc = self.plot_func.cdf_plotter(fin_jmb_cropped)
                    ax.plot(fsortc, fpopc, color = self.colours[k_], linestyle = ":")
                fvals[k_].append(fin_jmb)
                nsamp[k_].append(nsamples)
                ncrop[k_].append(nsamples_cropped)
            self.clean_plot.tickers(ax, "plot")
            ax.legend(prop={'size': self.clean_plot.axlabel_size})
            fig.savefig(fname, dpi=300, bbox_inches='tight')
            plt.close()

            output_path = "plotters/figures/binary_evolution/outputs/"
            ks_test = stats.ks_2samp(fvals[0][0], fvals[mcomp][0])
            q11, q31 = np.percentile(fvals[0], [25, 75])
            q11 = np.median(fvals[0]) - q11
            q31 -= np.median(fvals[0])

            q12, q32 = np.percentile(fvals[mcomp], [25, 75])
            q12 = np.median(fvals[mcomp]) - q12
            q32 -= np.median(fvals[mcomp])
            lines = ["Orbital parameter: "+str(file_name[plot_]),
                     "Model "+str(self.models[model_choices[0]])+" vs "+str(self.models[model_choices[mcomp]]),
                     "KS value: "+str(ks_test)+"\n",
                     "Median "+str(self.models[model_choices[0]])+": "+str(np.median(fvals[0][0])),
                     "IQR "+str(self.models[model_choices[0]])+": "+str(q11)+", "+str(q31)+"\n",
                     "Median "+str(self.models[model_choices[mcomp]])+": "+str(np.median(fvals[mcomp][0])),
                     "IQR "+str(self.models[model_choices[mcomp]])+": "+str(q12)+", "+str(q32),"\n"]

            if file_name[plot_] == "proj_sep":
                frac_bel0 = [i/j for i, j in zip(ncrop[0], nsamp[0])]
                frac0 = np.median(frac_bel0)
                q120, q320 = np.percentile(frac_bel0, [25, 75])
                q120 = np.median(frac0) - q120
                q320 -= np.median(frac0)

                frac_bel1 = np.asarray([i/j for i, j in zip(ncrop[mcomp], nsamp[mcomp])])
                frac1 = np.median(frac_bel1[~np.isnan(frac_bel1)])
                q121, q321 = np.percentile(frac_bel1[~np.isnan(frac_bel1)], [25, 75])
                q121 = np.median(frac1) - q121
                q321 -= np.median(frac1)

                lines = np.concatenate((lines, [str(self.models[model_choices[0]])+": Fraction > 25au Sep: "+str(frac0), "IQR "+str(q120)+", "+str(q320),
                                                str(self.models[model_choices[mcomp]])+"Fraction > 25au Sep: "+str(frac1), "IQR "+str(q121)+", "+str(q321)]),
                                      axis=None)
            with open(os.path.join(output_path, tname), 'w') as f:
                for line_ in lines:
                    f.write(line_+'\n')


    def time_evol_nJumbo(self, model_choices):
        """Plot evolution of nJuMBO in time"""

        leg_txt, plot_label, extra_str = self.clean_plot.model_layout(model_choices)
        fig, ax = plt.subplots()

        fiter = 0
        idx_frac = [ ]
        miny = 100
        for model_iter in model_choices:
            print("Processing ", self.models[model_iter])
            data = ReadData()
            data.proc_time_evol_JuMBO(self.models[model_iter], self.JuMBO_max_mass)

            file_path = "data/Simulation_Data/"+str(self.models[model_iter])+"/Processed_Data/Track_JuMBO/"
            
            med_fJuMBO = pd.read_hdf(file_path+"frac_JuMBO", 'Data')
            nearest = abs(med_fJuMBO.iloc[:,0] - 0.09)
            nearest_idx = nearest.idxmin()
            idx_frac.append(nearest_idx)

            IQR_low = pd.read_hdf(file_path+"IQR_low_fJuMBO", 'Data')
            IQR_high = pd.read_hdf(file_path+"IQR_high_fJuMBO", 'Data')

            time = np.linspace(0, 10**2, len(med_fJuMBO.iloc[:,0]))
            time = [10*i for i in time]
            vals = np.log10(IQR_low.iloc[:,0])
            vals = min(vals[np.isfinite(vals)])
            print(vals)
            miny = min(vals, miny)
            
            ax.plot(time, np.log10(med_fJuMBO.iloc[:,0]), color = self.colours[fiter], label = leg_txt[fiter])
            ax.plot(time, np.log10(IQR_high.iloc[:,0]), color = self.colours[fiter], alpha = 0.3)
            ax.plot(time, np.log10(IQR_low.iloc[:,0]), color = self.colours[fiter], alpha = 0.3)
            ax.plot([time[nearest_idx], time[nearest_idx]], [1.05*miny, np.log10(med_fJuMBO.iloc[nearest_idx])], 
                    color = self.colours[fiter])
            ax.fill_between(time, np.log10(IQR_low.iloc[:,0]), np.log10(IQR_high.iloc[:,0]), color = self.colours[fiter], alpha = 0.3)
            ax.scatter(time[nearest_idx], np.log10(med_fJuMBO.iloc[nearest_idx]), color = self.colours[fiter], zorder = 5)
            fiter += 1
        ax.set_xlabel(r"$t$ [kyr]", fontsize = self.clean_plot.axlabel_size)
        ax.set_ylabel(r"$\log_{10}(\langle N_{\mathrm{JuMBO}}\rangle/\langle N_{\mathrm{JMO}}\rangle)$", fontsize = self.clean_plot.axlabel_size)
        ax.set_xlim(0,1000)
        ax.set_ylim(1.05*miny, 0)
        self.clean_plot.tickers(ax, "plot")
        ax.legend(prop={'size': self.clean_plot.axlabel_size})
        plt.savefig("plotters/figures/system_evolution/"+str(extra_str)+"fJuMBO_evol.pdf", dpi=700, bbox_inches='tight')
        fig.clear()

        fig, ax = plt.subplots()
        fiter = 0
        for model_iter in model_choices:
            file_path = "data/Simulation_Data/"+str(self.models[model_iter])+"/Processed_Data/Track_JuMBO/"
            
            med_val= pd.read_hdf(file_path+"SemiMajor", 'Data')
            IQR_low = pd.read_hdf(file_path+"SemiMajor_IQRL", 'Data')
            IQR_high = pd.read_hdf(file_path+"SemiMajor_IQRH", 'Data')
            time = np.linspace(0, 10**3, len(med_val))
            ax.plot(time, med_val, color = self.colours[fiter], label = leg_txt[fiter])
            ax.plot(time, IQR_high.iloc[:,0], color = self.colours[fiter], alpha = 0.3)
            ax.plot(time, IQR_low.iloc[:,0], color = self.colours[fiter], alpha = 0.3)
            ax.fill_between(time, IQR_low.iloc[:,0], IQR_high.iloc[:,0], color = self.colours[fiter], alpha = 0.3)
            if model_iter != 3 and model_iter != 4 and model_iter != 5:
                frac_dt = idx_frac[fiter]
            else:
                frac_dt = -1
            ax.scatter(time[frac_dt], med_val.iloc[frac_dt], color = self.colours[fiter], zorder = 5)
            with open(os.path.join("plotters/figures/system_evolution/outputs/"+self.models[model_iter]+'_evol.txt'), 'w') as f:
                f.write(r"t = {:.1f} kyr, a = {:.1f}, IQR [{:.1f}, {:.1f}]".format(time[frac_dt+1], 
                        med_val.iloc[frac_dt,0], IQR_low.iloc[frac_dt,0], 
                        IQR_high.iloc[frac_dt,0]))
            fiter += 1
        ax.set_xlabel(r"$t$ [kyr]", fontsize = self.clean_plot.axlabel_size)
        ax.set_ylabel(r"$a$ [au]", fontsize = self.clean_plot.axlabel_size)
        ax.set_xlim(10,990)
        self.clean_plot.tickers(ax, "plot")
        ax.legend(prop={'size': self.clean_plot.axlabel_size})
        plt.savefig("plotters/figures/system_evolution/"+str(extra_str)+"sem_evol.pdf", dpi=700, bbox_inches='tight')
        fig.clear()

    def two_point_correlation(self, model_iter, dt_crop):
        """Function plotting the two-point correlation function
           Input:
           model_iter: The model wishing to analyse
           dt_crop:    Boolean to investigate final snapshot, or when fJuMBO = 9%
        """
        
        directory = "data/Simulation_Data/"+str(self.models[model_iter])+"/"
        dir_configs = glob.glob(os.path.join(str(directory)+"simulation_snapshot/**"))
        if (dt_crop):
            file = os.path.join("plotters/figures/system_evolution/outputs/"+str(self.models[model_iter])+"_evol.txt")
            with open(file) as f:
                output = csv.reader(f)
                for row_ in output:
                    time = float(row_[0][4:-3])
                    file_idx = int(np.floor(time/10))
            fname = "plotters/figures/system_evolution/two_point_corr"+str(self.models[model_iter])+"_crop.pdf"
        else:
            file_idx = -1
            fname = "plotters/figures/system_evolution/two_point_corr"+str(self.models[model_iter])+".pdf"

        jmb_jmo_arr = [ ]
        jmb_jmb_arr = [ ]
        str_str_arr = [ ]
        str_jmb_arr = [ ]
        str_jmo_arr = [ ]

        citer = 1
        file_no = 0
        for config_ in dir_configs:
            file_no += 1
            print("Reading file #", file_no)

            sim_snapshot = natsort.natsorted(glob.glob(os.path.join(config_+"/*")))[file_idx]
            mass = [ ]
            key = [ ]
            mkey = [ ]
            iter = 0 
            for run_ in self.fin_bsim_iter[model_iter]:
                if run_ == citer:
                    mass = np.concatenate((mass, self.fin_bmass[model_iter][iter]))
                    key = np.concatenate((key, self.fin_bkeys[model_iter][iter]))
                iter += 1

            iter = 0
            for run_ in self.fin_msim_iter[model_iter]:
                if run_ == citer:
                    mkey = np.concatenate((mkey, self.fin_mkeys[model_iter][iter]))
                iter += 1
                
            data = read_set_from_file(sim_snapshot, "hdf5")
            star = data[data.mass >= self.Star_min_mass]
            jmo = data[data.mass < self.Star_min_mass]
            for dmass_, mass_ in itertools.product(jmo.key, key):
                dmass_str = str(dmass_)[1:9]
                mass_str = str(mass_)[2:10]
                if dmass_str == mass_str:
                    idx = np.where(key == mass_)[0]
                    if idx%2 == 0:
                        pidx = idx+1
                    else:
                        pidx = idx-1
                    pmass = mass[pidx] * (1 | units.MSun)
                    if pmass < self.Star_min_mass:
                        jmo[jmo.key == dmass_].jmb = 1

            if len([jmo.jmb == 1]) != 0:
                jumbo = jmo[jmo.jmb == 1]
            else:
                print("Curious")
                jumbo = data - data
            FF = data - star - jumbo
            
            jmb_jmo = [ ]
            jmb_jmb = [ ]
            str_str = [ ]
            str_jmb = [ ]
            str_jmo = [ ]

            pset = [data, jumbo, FF]
            str_dist_arr = [str_str, str_jmb, str_jmo]
            jmb_dist_arr = [jmb_jmb, jmb_jmo]
            for star_ in star:
                sdata = star - star_
                pset[0] = sdata
                for type_, darr_ in zip(pset, str_dist_arr):
                    if len(type_) > 0:
                        dx = type_.x - star_.x
                        dy = type_.y - star_.y
                        dz = type_.z - star_.z
                        dist = np.sqrt(dx**2+dy**2+dz**2)
                        darr_.append(min(dist).value_in(units.au))
            
            pset = [data, jumbo, FF]
            for jmb_ in jumbo:
                pset[1] = jumbo - jmb_
                for type_, darr_ in zip(pset, jmb_dist_arr):
                    if len(type_) > 0:
                        dx = type_.x - jmb_.x
                        dy = type_.y - jmb_.y
                        dz = type_.z - jmb_.z
                        dist = np.sqrt(dx**2+dy**2+dz**2)
                        dist_jmb = np.sort(jmb_dist_arr)[1]
                        if len(type_) == len(jumbo):
                            darr_.append(min(dist_jmb).value_in(units.au))
                        else:
                            darr_.append(min(dist).value_in(units.au))
                            
            jmb_jmo_arr = np.concatenate((jmb_jmo_arr, jmb_jmo), axis=None)
            jmb_jmb_arr = np.concatenate((jmb_jmb_arr, jmb_jmb), axis=None)
            str_str_arr = np.concatenate((str_str_arr, str_str), axis=None)
            str_jmb_arr = np.concatenate((str_jmb_arr, str_jmb), axis=None)
            str_jmo_arr = np.concatenate((str_jmo_arr, str_jmo), axis=None)

            citer += 1
        data_arr = [jmb_jmo_arr, jmb_jmb_arr, str_str_arr, str_jmb_arr, str_jmo_arr]

        label_arr = ["JuMBO-JMO", "JuMBO-JuMBO", "Star-Star", "Star-JuMBO", "Star-JMO"]
        colours_arr = ["Black", "Red", "Blue", "pink", "dodgerblue"]

        fig, ax = plt.subplots()
        ax.set_xlabel(r"$\log_{10} r_{ij}$ [au]", fontsize = self.clean_plot.axlabel_size)
        ax.set_ylabel(r"CDF", fontsize = self.clean_plot.axlabel_size)
        for data_, label_, colour_ in zip(data_arr, label_arr, colours_arr):
            dsort = np.sort(data_)
            diter = np.asarray([i for i in enumerate(dsort)])
            diter = diter[:,0]
            diter /= max(diter)
            ax.plot(np.log10(dsort), diter, label = label_, color = colour_)
        ax.set_ylim(0,1)
        ax.legend(prop={'size': self.clean_plot.axlabel_size})
        self.clean_plot.tickers(ax, "plot")
        fig.savefig(fname, dpi=700, bbox_inches='tight')
        plt.show()
        plt.clf()