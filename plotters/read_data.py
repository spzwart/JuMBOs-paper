import csv
import glob
import natsort
import numpy as np
import os
import pandas as pd
import pickle as pkl
from itertools import combinations

from amuse.datamodel import Particles
from amuse.ext.orbital_elements import orbital_elements_from_binary
from amuse.io.base import read_set_from_file
from amuse.units import units, constants

class ReadData(object):
    def process_final_snapshot(self, model, dt_crop):
        """Process final snapshot and store all detected binaries
        
           Inputs:
           model:   Model to process data of
           dt_crop: Boolean (1 dt(Njmb/NJMO = 0.09) || 0 Final dt)
        """

        dir_path = "data/Simulation_Data/"+str(model)+"/"
        dir_configs = glob.glob(os.path.join(dir_path+"simulation_snapshot/**"))
        
        if (dt_crop):
            file = os.path.join("plotters/figures/system_evolution/outputs/"+str(model)+"_evol.txt")
            with open(file) as f:
                output = csv.reader(f)
                for row_ in output:
                    time = float(row_[0][4:-3])
                    file_idx = int(np.floor(time/10))
            dir = dir_path+"Processed_Data/Final_Binaries_crop/"
            dir_multi = dir_path+"Processed_Data/Final_MultiSysts_crop/"

        else:
            file_idx = -1
            dir = dir_path+"Processed_Data/Final_Binaries/"
            dir_multi = dir_path+"Processed_Data/Final_MultiSysts/"

        bound_threshold = 2000 | units.au #1000 | units.au
        no_files = 0
        no_sim_rn = 0
        for config_ in dir_configs:
            file_name = natsort.natsorted(glob.glob(config_+"/*"))[file_idx]
            msys_file = 0
            no_sim_rn += 1
            if "DONE" not in file_name:
                print("Processing Data: ", file_name)
                syst = 0
                data = read_set_from_file(file_name, "hdf5")
                components = data.connected_components(threshold = 
                                                       bound_threshold)
                for c in components:
                    if len(c) > 1:
                        syst += 1
                        multi_syst = 0

                        bin_combo = list(combinations(c, 2)) #Permuting through connected components
                        keys = [ ]
                        msyst_type = [ ]
                        msyst_keys = [ ]
                        msyst_mass = [ ]
                        for bin_ in bin_combo:
                            bin_sys = Particles()  
                            bin_sys.add_particle(bin_[0])
                            bin_sys.add_particle(bin_[1])
                            bin_sys.move_to_center()

                            KE1 = 0.5*bin_sys[0].mass*bin_sys[0].velocity.length()**2
                            KE2 = 0.5*bin_sys[1].mass*bin_sys[1].velocity.length()**2
                            PE = bin_sys.potential_energy()
                            if (KE1+PE < (0 | units.J)) and (KE2+PE < (0 | units.J)):
                                kepler_elements = orbital_elements_from_binary(bin_sys, G=constants.G)
                                semimajor = kepler_elements[2]
                                eccentric = kepler_elements[3]
                                if (eccentric<1) and semimajor<0.5*bound_threshold:
                                    no_files += 1
                                    multi_syst += 1

                                    msyst_keys.append([int(bin_sys[0].key), int(bin_sys[1].key)])
                                    msyst_mass = np.concatenate((msyst_mass, 
                                                        [bin_sys[0].mass, bin_sys[1].mass]), 
                                                        axis=None)
                                    msyst_type = np.concatenate((msyst_type, 
                                                        [bin_sys[0].name , bin_sys[1].name]), 
                                                        axis=None)

                                    keys = [bin_[0].key, bin_[1].key]
                                    bin_type = [bin_[0].name, bin_[1].name]
                                    bin_mass = [bin_[0].mass.value_in(units.MSun), 
                                                bin_[1].mass.value_in(units.MSun)]
                                    
                                    proj_xy = np.sqrt((bin_[0].x - bin_[1].x)**2 
                                                    + (bin_[0].y-bin_[1].y)**2)[0]
                                    proj_xz = np.sqrt((bin_[0].x - bin_[1].x)**2 
                                                    + (bin_[0].z-bin_[1].z)**2)[0]
                                    proj_yz = np.sqrt((bin_[0].y - bin_[1].y)**2 
                                                    + (bin_[0].z-bin_[1].z)**2)[0]
                                    mproj = (proj_xy+proj_xz+proj_yz)/3

                                    inclinate = kepler_elements[4]
                                    arg_peri = kepler_elements[5]
                                    asc_node = kepler_elements[6]
                                    true_anom = kepler_elements[7]

                                    fname = "Final_Binary_no"+str(no_sim_rn)+"_syst"+str(no_files)+".pkl"
                                    df_arr = pd.DataFrame()
                                    df_vals = pd.Series({"Simulation": no_sim_rn, 
                                                         "Keys ": keys, 
                                                         "Type" : bin_type, 
                                                         "Mass" : bin_mass, 
                                                         "Semi-major" : semimajor, 
                                                         "Inclinate" : inclinate, 
                                                         "Arg. of Pericenter" : arg_peri, 
                                                         "Ascending Node" : asc_node, 
                                                         "True anomaly" : true_anom, 
                                                         "Eccentricity" : eccentric,
                                                         "Mean Proj. Sep" : mproj})
                                    df_arr = df_arr._append(df_vals, ignore_index=True)
                                    df_arr.to_pickle(os.path.join(dir, fname))  
                        
                        if multi_syst > 1:
                            all_keys = [ ]
                            for sublist_ in msyst_keys:
                                for item_ in sublist_:
                                    all_keys.append(item_)
                                    
                            uq_keys, counts = np.unique(all_keys, return_counts=True)
                            uq_keys = uq_keys[counts > 1]
                            nhost = [[] for i in range(len(uq_keys))]
                            syst_key = [[] for i in range(len(uq_keys))]
                            host_keys = [[] for i in range(len(uq_keys))]

                            hidx = 0
                            tracked = [ ]
                            for host_ in uq_keys:
                                if host_ not in tracked:
                                    nhost[hidx] = 1
                                    tracked = np.concatenate((tracked, host_), axis=None)
                                    syst_key[hidx] = np.concatenate((syst_key[hidx], host_), axis=None)
                                    host_keys[hidx] = np.concatenate((host_keys[hidx], host_), axis=None)
                                    for sublist_ in msyst_keys:
                                        sublist_ = np.asarray(sublist_)
                                        if host_ in sublist_:
                                            partner = int(sublist_[sublist_ != host_])
                                            syst_key[hidx] = np.concatenate((syst_key[hidx], partner), axis=None)
                                            if partner in uq_keys: #Check if partner also host multiple systems
                                                nhost[hidx] += 1
                                                host_keys[hidx] = np.concatenate((host_keys[hidx], partner), axis=None)
                                                for sublist_ in msyst_keys:
                                                    sublist_ = np.asarray(sublist_)
                                                    if partner in sublist_:
                                                        partner2 = int(sublist_[sublist_ != partner])
                                                        tracked = np.concatenate((tracked, partner), axis=None)
                                                        if partner2 not in syst_key[hidx]:
                                                            syst_key[hidx] = np.concatenate((syst_key[hidx], partner2), axis=None)
                                hidx += 1
                            
                            idx = 0
                            for syst_ in syst_key:
                                msys_file += 1
                                if len(syst_) > 0:
                                    pset = Particles()
                                    for key_ in syst_:
                                        if len(pset[pset.key == key_]) == 0:
                                            pset.add_particle(data[data.key == key_])
                                    print("Nhosts: ", nhost[idx], "Nsyst: ", len(pset), "Host: ", host_keys[idx])
                                    print(pset)

                                    df_arr = pd.DataFrame()
                                    df_vals = pd.Series({"Simulation": no_sim_rn, 
                                                         "Keys": pset.key, 
                                                         "Type": pset.name,
                                                         "Mass": [pset.mass],
                                                         "Host": host_keys[idx],
                                                         "nHosts": nhost[idx]})
                                    df_arr = df_arr._append(df_vals, ignore_index=True)

                                    fname = "Final_Binary_no"+str(no_sim_rn)+"_syst"+str(msys_file)+".pkl"
                                    df_arr.to_pickle(os.path.join(dir_multi, fname))

                                idx += 1

                if "DONE" not in file_name:
                    os.rename(file_name, file_name+"DONE")

    def proc_time_evol_JuMBO(self, model, JMO_mass):
        """Process final snapshot and store all detected binaries
        
           Input:
           model:    Model processing data of
           JMO_mass: Maximum JMO mass
        """

        dir_configs = glob.glob(os.path.join("data/Simulation_Data/"+str(model)+"/simulation_snapshot/**"))
        file_path = "data/Simulation_Data/"+str(model)+"/Processed_Data/Track_JuMBO/"
        file_checker = glob.glob(os.path.join(file_path+"*"))
        if len(file_checker) > 0:
            print("Data analysis already done for ", model)
            None
        else:
            snap_dt = 102
            snap_frac = 1
            bound_threshold = 2000 | units.au
            snapshot = natsort.natsorted(glob.glob(dir_configs[0]+"/*"))[0]
            if model == "Fractal_rvir0.5_FF_10Myr":
                nsnap = 1002
            else:
                nsnap = 101

            semi_arr = [[ ] for i in range(nsnap)]
            ecc_arr = [[ ] for i in range(nsnap)]
            no_JuMBO_all = [ ]
            no_JMO_all = [ ]
            idx_snap = [ ]

            no_configs = 0
            for config_ in dir_configs:
                print("Reading #", no_configs, ": ", config_, len(snapshot))
                no_configs += 1
                files = natsort.natsorted(glob.glob(config_+"/*"))

                no_JuMBO = np.zeros(len(files)+1)
                no_JMO = np.zeros(len(files)+1)

                dt = 0
                data = read_set_from_file(files[0], "hdf5")
                data = data[data.mass <= JMO_mass]
                no_JuMBO[dt] = len(data[data.name == "JuMBOs"])
                no_JMO[dt] = len(data)
                for snap_ in files:
                    print("File: ", snap_)
                    dt += 1
                    data = read_set_from_file(snap_, "hdf5")
                    data = data[data.mass <= JMO_mass]
                    no_JMO[dt] = len(data)
                    ecc_dt = [ ]
                    sem_dt = [ ]

                    components = data.connected_components(threshold = bound_threshold)
                    tracked_keys = [ ]
                    for c in components:
                        if len(c) > 1:
                            c.move_to_center()
                            multi_syst = 0
                            bin_combo = list(combinations(c, 2))
                            for bin_ in bin_combo:
                                bin_sys = Particles()  
                                bin_sys.add_particle(bin_[0])
                                bin_sys.add_particle(bin_[1])

                                kepler_elements = orbital_elements_from_binary(bin_sys, G=constants.G)
                                semimajor = kepler_elements[2]
                                eccentric = kepler_elements[3]
                                if (eccentric<1) and semimajor<(1000 | units.AU):
                                    multi_syst += 1
                                    no_JuMBO[dt] += 2

                                    ecc_dt.append(eccentric)
                                    sem_dt.append(semimajor.value_in(units.au))
                                    
                                    tracked_keys.append(bin_[0].key)
                                    tracked_keys.append(bin_[1].key)

                    uq_keys, counts = np.unique(tracked_keys, return_counts=True)
                    fiter = 0
                    del_idx = [ ]
                    for key_ in tracked_keys:
                        if key_ in uq_keys[counts > 1]:
                            del_idx.append(int(np.floor(fiter/2)))
                        fiter += 1
                    ecc_dt = np.delete(ecc_dt, del_idx)
                    sem_dt = np.delete(sem_dt, del_idx)
                    
                    ecc_arr[dt-1] = np.concatenate((ecc_arr[dt-1], ecc_dt), axis=None)
                    semi_arr[dt-1] = np.concatenate((semi_arr[dt-1], sem_dt), axis=None)
                    for count_ in counts[counts>1]:
                        no_JuMBO[dt] -= (count_-1)
                    no_JuMBO[dt] = min(no_JuMBO[dt], no_JMO[dt])

                    dt_frac = no_JuMBO[dt]/no_JMO[dt]
                    val1 = abs(0.09-snap_frac)
                    val2 = abs(0.09-dt_frac)
                    if val2 < val1:
                        snap_dt = dt
                        snap_frac = dt_frac
                no_JMO_all.append(no_JMO)
                no_JuMBO_all.append(no_JuMBO)
                idx_snap.append(snap_dt)

            print("...writing files...")
            with open(os.path.join(file_path+"snap_idx.txt"), 'w') as f:
                f.write(str(idx_snap))

            med_JuMBO = [np.median(x) for x in zip(*no_JuMBO_all)]
            IQR_JB_low = [np.percentile(x, 25) for x in zip(*no_JuMBO_all)]
            IQR_JB_hig = [np.percentile(x, 75) for x in zip(*no_JuMBO_all)]
            med_JMO = [np.median(x) for x in zip(*no_JMO_all)]

            frac_JuMBO = [i/j for i, j in zip(med_JuMBO, med_JMO)]
            IQR_high = [i/j for i, j in zip(IQR_JB_hig, med_JMO)]
            IQR_low = [i/j for i, j in zip(IQR_JB_low, med_JMO)]

            file_name = ["frac_JuMBO", "IQR_low_fJuMBO", "IQR_high_fJuMBO"]
            file_data = [frac_JuMBO, IQR_high, IQR_low]
            for fname_, fdata_ in zip(file_name, file_data):
                data_arr = pd.DataFrame(fdata_)
                data_arr.to_hdf(os.path.join(file_path+fname_), key ='Data', mode = 'w')

            for data_, fname_ in zip([semi_arr, ecc_arr], ["SemiMajor", "Eccentric"]):
                median = [ ]
                IQR_low = [ ]
                IQR_high = [ ]
                for subdata_ in data_[1:]:
                    median.append(np.median(subdata_))
                    if not np.isnan(np.median(subdata_)):
                        IQR_low.append(np.percentile(subdata_, 25))
                        IQR_high.append(np.percentile(subdata_, 75))
                    else:
                        IQR_low.append(np.nan)
                        IQR_high.append(np.nan)
                
                data_arr_med = pd.DataFrame(median)
                data_arr_IQRL = pd.DataFrame(IQR_low)
                data_arr_IQRH = pd.DataFrame(IQR_high)

                data_arr_med.to_hdf(os.path.join(file_path+fname_), 
                                    key ='Data', mode = 'w')
                data_arr_IQRL.to_hdf(os.path.join(file_path+fname_+"_IQRL"), 
                                     key ='Data', mode = 'w')
                data_arr_IQRH.to_hdf(os.path.join(file_path+fname_+"_IQRH"), 
                                     key ='Data', mode = 'w')

    def proc_time_evol_mixed(self, model, JMO_mass):
        """Process final snapshot and store semi + ecc of detected mixed systems"""

        dir_configs = glob.glob(os.path.join("data/Simulation_Data/"+str(model)+"/simulation_snapshot/**"))
        file_path = "data/Simulation_Data/"+str(model)+"/Processed_Data/Track_JuMBO/"
        
        bound_threshold = 2000 | units.au
        semi_arr = [ ]
        ecc_arr = [ ]
        mprim_arr = [ ]
        q_arr = [ ]
        file_iter = []
        
        no_configs = 0
        for config_ in dir_configs:
            print("Reading: ", config_)
            no_configs += 1
            files = natsort.natsorted(glob.glob(config_+"/*"))

            dt = 0
            for snap_ in files[1:]:
                dt += 1
                data = read_set_from_file(snap_, "hdf5")

                components = data.connected_components(threshold = bound_threshold)
                for c in components:
                    if len(c) > 1:
                        c.move_to_center()
                        bin_combo = list(combinations(c, 2))
                        for bin_ in bin_combo:
                            if (min(bin_[0].mass, bin_[1].mass) <= JMO_mass) \
                                and (max(bin_[0].mass, bin_[1].mass) > 10*JMO_mass): 
                                bin_sys = Particles()  
                                bin_sys.add_particle(bin_[0])
                                bin_sys.add_particle(bin_[1])

                                kepler_elements = orbital_elements_from_binary(bin_sys, G=constants.G)
                                semimajor = kepler_elements[2]
                                eccentric = kepler_elements[3]
                                if (eccentric<1) and semimajor<(1000 | units.AU):
                                    semi_arr.append(semimajor)
                                    ecc_arr.append(eccentric)
                                    mprim_arr.append(max(bin_[0].mass, 
                                                         bin_[1].mass))
                                    q_arr.append(min(bin_[0].mass/bin_[1].mass, 
                                                     bin_[1].mass/bin_[0].mass))
                                    file_iter.append(dt)
        
        data_arr = pd.DataFrame([file_iter, mprim_arr, q_arr, semi_arr, ecc_arr])
        data_arr.to_hdf(os.path.join(file_path+"mixed_sys_data"), key ='Data', mode = 'w')

    def read_init_data(self, directory):
        """Store data from simulations"""

        self.init_bkeys = [ ]
        self.init_bmass = [ ]      #Raw data given in MSun
        self.init_bsemi = [ ]      #Raw data given in au
        self.init_becce = [ ]
        self.init_bincl = [ ]      #Raw data given in deg
        self.init_barg_peri = [ ]  #Raw data given in deg
        self.init_basc_node = [ ]  #Raw data given in deg
        self.init_banom = [ ]      #Raw data given in deg
        self.no_JuMBOs = 0

        dir_init_bin = glob.glob(os.path.join(str(directory)+"initial_binaries/*"))
        for file_ in dir_init_bin:
            with open(file_) as f:
                line = f.readlines()
                temp_mass = [ ]
                row = 0
                for data_ in line:
                    row += 1
                    if row%11 == 0:
                        self.no_JuMBOs += 1
                    elif row%11 == 1:
                        temp_keys = [ ]
                        temp_keys.append(float(data_[5:-2]))
                    elif row%11 == 2:
                        temp_keys.append(float(data_[5:-2]))
                        self.init_bkeys.append(temp_keys)
                    elif row%11 == 3:
                        temp_mass = [ ]
                        temp_mass.append(float(data_[4:-6]))
                    elif row%11 == 4:
                        temp_mass.append(float(data_[4:-6]))
                        self.init_bmass.append(temp_mass)
                    elif row%11 == 5:
                        self.init_bsemi.append(float(data_[17:-3]))
                    elif row%11 == 6:
                        self.init_becce.append(float(data_[14:]))
                    elif row%11 == 7:
                        self.init_bincl.append(float(data_[13:-4]))
                    elif row%11 == 8:
                        self.init_barg_peri.append(float(data_[23:-8]))
                    elif row%11 == 9:
                        self.init_basc_node.append(float(data_[23:-8]))
                    elif row%11 == 10:
                        self.init_banom.append(float(data_[13:-8]))

    def read_final_data(self, directory, dt_crop):
        """Read properties of binaries existing"""

        print("...Reading Final Binaries...")
        dir_configs = glob.glob(os.path.join(str(directory)+"Processed_Data/Final_Binaries/*"))
        if (dt_crop):
            dir_configs = glob.glob(os.path.join(str(directory)+"Processed_Data/Final_Binaries_crop/*"))

        self.jmb_idx = [ ]
        self.fin_bsim_iter = [ ]
        self.fin_bkeys = [ ]
        self.fin_btype = [ ]
        self.fin_bmass = [ ] 
        self.fin_bsemi = [ ] 
        self.fin_bproj = [ ]
        self.fin_becce = [ ]
        self.fin_bincl = [ ]
        self.fin_barg_peri = [ ] 
        self.fin_basc_node = [ ]
        self.fin_banom = [ ]

        file_iter = 0 
        for file_ in dir_configs:
            with open(file_, 'rb') as input_file:
                file_iter += 1
                if file_iter%500 == 0:
                    print("Reading File #", file_iter, ": ", file_)
                data_file = pkl.load(input_file)
                
                self.fin_bsim_iter.append(data_file.iloc[0][0])
                self.fin_bkeys.append(data_file.iloc[0][1])
                self.fin_btype.append(data_file.iloc[0][2])
                self.fin_bmass.append([data_file.iloc[0][3][0], 
                                       data_file.iloc[0][3][1]])
                self.fin_bsemi.append(data_file.iloc[0][4].value_in(units.au))
                self.fin_bincl.append(data_file.iloc[0][5])
                self.fin_barg_peri.append(data_file.iloc[0][6])
                self.fin_basc_node.append(data_file.iloc[0][7])
                self.fin_banom.append(data_file.iloc[0][8])
                self.fin_becce.append(data_file.iloc[0][9])
                self.fin_bproj.append(data_file.iloc[0][10].value_in(units.au))

                if max(self.fin_bmass[-1]) <= 0.013:
                    self.jmb_idx.append(len(self.fin_bkeys)-1)
    
    def read_final_multi_data(self, directory, dt_crop):
        """Read properties of N>2 systems"""

        dir_configs = glob.glob(os.path.join(str(directory)+"Processed_Data/Final_MultiSysts/*"))
        if (dt_crop):
            dir_configs = glob.glob(os.path.join(str(directory)+"Processed_Data/Final_MultiSysts_crop/*"))

        self.fin_msim_iter = [ ]
        self.fin_mkeys = [ ]
        self.fin_mtype = [ ]
        self.fin_midx = [ ]
        self.fin_mmass = [ ]

        for file_ in dir_configs:
            with open(file_, 'rb') as input_file:
                data_file = pkl.load(input_file)
                self.fin_msim_iter.append(data_file.iloc[0][0])
                self.fin_mkeys.append(data_file.iloc[0][1])
                self.fin_mtype.append(data_file.iloc[0][2])
                self.fin_mmass.append(data_file.iloc[0][3])
                