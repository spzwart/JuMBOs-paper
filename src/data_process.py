import fnmatch
import glob
import os
import pandas as pd
import warnings

from amuse.ext.LagrangianRadii import LagrangianRadii as LagRad
from amuse.ext.orbital_elements import orbital_elements_from_binary
from amuse.units import units, constants

warnings.simplefilter(action = 'ignore', category = pd.errors.PerformanceWarning)
def file_reset():
    """Function to remove all files from a directory"""

    dir = r"data/Ph4"
    filelist = glob.glob(os.path.join(dir, "*"))
    filelist_1 = filelist[filelist != "data/Ph4/simulation_snapshot"]
    filelist_2 = filelist[filelist == "data/Ph4/simulation_snapshot"]
    
    print("!!!Warning: Deleting all files!!!")
    ans = input("Are you sure? (Y|N)")
    filelist = glob.glob(os.path.join(dir+str(filelist_1), "*"))
    if ans.lower() == "y":
        for f in filelist_1:
            os.remove(f)
    filelist = glob.glob(os.path.join(dir+str(filelist_2[0]), "*"))
    for f in filelist_2:
        os.remove(f)
    STOP
#file_reset()

class DataProcess(object):
    def __init__(self, data_direc):
        self.sim_prefix = "Simulation_Data"
        self.path = data_direc
        self.file_count = len(fnmatch.filter(os.listdir(self.path+str("initial_binaries/")), '*'))
        self.read_files = False

    def init_arr(self, grav_code, pset, rcluster):
        #Initialising data storage
        self.time_arr = [grav_code.model_time.value_in(units.kyr)]
        self.lag25_arr = [LagRad(pset)[5].value_in(units.pc)]
        self.lag50_arr = [LagRad(pset)[6].value_in(units.pc)]
        self.lag75_arr = [LagRad(pset)[7].value_in(units.pc)]
        self.cejec_lim = [rcluster.value_in(units.pc)]

        self.parti_key = [pset.key]
        self.parti_mass = [pset.mass.value_in(units.MSun)]
        
        self.dE = [ ]
        self.dA = [ ]
    
    def append_data(self, grav_code, pset, rcluster, dE, dA):
        """Function to append data files"""

        self.time_arr.append(grav_code.model_time.value_in(units.kyr))
        self.lag25_arr.append(LagRad(pset)[5].value_in(units.pc))
        self.lag50_arr.append(LagRad(pset)[6].value_in(units.pc))
        self.lag75_arr.append(LagRad(pset)[7].value_in(units.pc))
        self.cejec_lim.append(rcluster.value_in(units.pc))
        
        self.parti_key.append(pset.key)
        self.parti_mass.append(pset.mass.value_in(units.MSun))
        self.dE.append(dE)
        self.dA.append(dA)

    def energy_data(self, time, dE, dA, config):
        """Initialise data set to track energy and angular momentum error
        
           Inputs:
           time:       Array with simulation times      (yr)
           d(E/A):     Energy/ang. momentum error
           config:     Simulation configuration
        """
        
        file_name = "energy_data/energy_sim_config_"+str(config)+"_2.0.h5"
        event_array = pd.DataFrame([time, dE, dA])
        event_array.to_hdf(os.path.join(self.path+file_name), key ='dfEnergy', mode = 'w')

        if (self.read_files):
            print("READING ENERGY")
            print(pd.read_hdf('data/'+str(self.sim_prefix)+'/energy_data/energy_sim0.h5', 'dfEnergy'))

    def event_data(self, event_time, event_key, event_type, event_mass, config):
        """Initialise data set tracking dissolution, mergers, ejections

           Inputs:
           event_time: Time which events occur   (yr)
           event_key:  Key of relevant particles
           event_type: Flag denoting event type (Cluster)
           config:     Simulation configuration
        """
           
        file_name = "event_data/events_sim_"+str(config)+"_2.0.h5"
        event_array = pd.DataFrame([event_time, event_key, event_type, event_mass])
        event_array.to_hdf(os.path.join(self.path+file_name), key ='dfEvent', mode = 'w')

        if (self.read_files):
            print("READING EVENT")
            print(pd.read_hdf('data/'+str(self.sim_prefix)+'/event_data/events_sim0.h5', 'dfEvent'))

    def ini_parti_data(self, pset, config):
        """Function storing initial particle data"""
        file_name = "initial_particle_data/init_parti_sim"+str(config)+".h5"
        ic_parti_array = pd.DataFrame([pset.key, pset.mass, pset.sub_worker_radius])
        ic_parti_array.to_hdf(os.path.join(self.path+file_name), key ='dfPartiIC', mode = 'w')
        if (self.read_files):
            print("READING EVENT")
            print(pd.read_hdf('data/'+str(self.sim_prefix)+'/initial_particle_data/init_parti_sim0.h5', 'dfPartiIC'))

    def lagrangian_data(self, time, lag25_data, lag50_data, lag75_data, cluster_lim, config):
        """Store Lagrangian radius evolution
        
           Inputs:
           time:        Array with simulation times      (yr)
           lagXYZ_data: Array with XYZ lagrangian radius (pc)
           cluster_lim: Array with radius defining cluster threshold
           config:      Simulation configuration
        """

        file_name = "lagrangian_data/lagrangian_sim"+str(config)+"_2.0.h5"
        LG_array = pd.DataFrame([time, lag25_data, lag50_data, lag75_data, cluster_lim])
        LG_array.to_hdf(os.path.join(self.path+file_name), key ='dfLag', mode='w')
        
        if (self.read_files):
            print("READING LAGRANGIAN")
            print(pd.read_hdf('data/'+str(self.sim_prefix)+'/lagrangian_data/lagrangian_sim0.h5', 'dfLag'))
