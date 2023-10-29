import numpy as np
from amuse.units import units
from amuse.plot import plot
from matplotlib import pyplot as plt

import pandas as pd

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-f", dest="filename", default = "Plummer_N2500n600csmf_R1pcQ05_R1.csv",
                      help="input filename[%default]")
    result.add_option("-q", action='store_true',
                      dest="mass_ratio", default="False",
                      help="mass_ratio [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()

    #data = pd.read_csv("Plummer_N2500n300cs_R05pcQ05_R1.csv")
    data = pd.read_csv(o.filename)
    print(data)
    amin = data["amin"] 
    print(amin)
    amax = data["amax"] 
    n_jumbos = data["n"]
    a_jumbos = data["a"] 
    da_jumbos = data["da"] 
    e_jumbos = data["e"]
    de_jumbos = data["de"]
    if o.mass_ratio==True:
        q_jumbos = data["q"]
        dq_jumbos = data["dq"]
    n_ffps = data["nff"]

    for i in range(len(amin)):
        plt.plot([amin[i], amax[i]],
                 [n_jumbos[i], n_jumbos[i]], linestyle='-', marker='o')
        #plot([amin[i], amax[i]], [n_jumbos[i], n_jumbos[i]])
    plt.show()    

    for i in range(len(amin)):
        ai = (amin[i]+amax[i])/2
        dai = (amax[i]-amin[i])/2
        aj = a_jumbos[i]
        daj = da_jumbos[i]
        plt.errorbar(ai, aj, xerr=dai, yerr=daj)
    plt.semilogy()    
    plt.show()

    if o.mass_ratio==True:
        for i in range(len(amin)):
            ai = (amin[i]+amax[i]/2)
            dai = (amax[i]-amin[i]/2)
            qj = (q_jumbos[i]+q_jumbos[i]/2)
            dqj = (dq_jumbos[i]/2)
            plt.errorbar(ai, qj, xerr=dai, yerr=dqj)
        plt.semilogy()    
        plt.show()

    f_jumbos = np.array(n_jumbos)/(2*np.array(n_jumbos)+np.array(n_ffps))
    for i in range(len(amin)):
        plt.plot([amin[i], amax[i]],
                 [f_jumbos[i], f_jumbos[i]],
                 linestyle='-', marker='o')
    plt.show()


