import numpy as np
from amuse.units import units


class PlotterFunctions(object):
    def cdf_plotter(self, data):
        """Sort arrays for CDF plots

           Inputs:
           data:   The data array
           unit:   Boolean denoting if input has units
        """

        data_sort = np.sort(data)
        data_iter = np.arange(0, len(data_sort), 1)/(len(data_sort)-1)
        return data_sort, data_iter

    def flatten_arr(self, arr, iter):
        """Flatten nested arrays"""

        temp_arr = [ ]
        for sublist in arr[iter]:
            for item_ in sublist:
                if item_ not in temp_arr:
                    temp_arr.append(item_)
        return temp_arr

    def multi_sys_flag(self, key_array, mass_array, type_array):
        """Function grouping together all multi-systems from their
           permutation combinations (i.e [1,2], [1,3] --> [1,2,3])
           
           Inputs:
           key_array:   Array hosting particle keys
           mass_array:  Array hosting system masses
           type_array:  Array hosting particle types
        """
        
        syst_pops = [ ]
        syst_types = [ ]
        multi_keys = [ ]
        for syst_ in range(len(key_array)):
            temp_keys = [ ]
            temp_type = [ ]
            for indiv_ in range(len(key_array[syst_])):
                if key_array[syst_][indiv_] not in temp_keys:
                    temp_keys.append(key_array[syst_][indiv_])
                    if mass_array[syst_][0][indiv_] <= (0.013 | units.MSun) \
                        and type_array[syst_][indiv_] != "star":
                        temp_type.append("JuMBOs")
                    else:
                        temp_type.append("star")
                        
            Npop = len(np.unique(key_array[syst_]))
            syst_pops.append(Npop)
            multi_keys.append(temp_keys)
            syst_types.append(temp_type)
            
        return [syst_pops, syst_types, multi_keys]