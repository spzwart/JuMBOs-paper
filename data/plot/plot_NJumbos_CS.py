import numpy as np
from amuse.units import units
from amuse.plot import plot
from matplotlib import pyplot as plt

import pandas as pd

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-f", dest="filename", default = "PlN2500n600cs_R0.25pcQ05_JuMBOs.csv",
                      help="input filename[%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()

    d = pd.read_csv(o.filename)
    M1 = d["M"] 
    M2 = d["m"] 
    for i in range(len(M1)):
        print(M1[i])
        if M1[i]<M2[i]:
            m = M1[i]
            M1[i] = M2[i]
            M2[i] = m
    #data = d[d['M'] > 0.6]
    #data = data[data['m'] > 0.6]
    #data = data[data['a'] < 10000]
    data=d
    print(data)
    M1 = data["M"] 
    M2 = data["m"] 
    a = data["a"] 
    e = data["e"] 
    x = data["x"] 
    y = data["y"] 
    z = data["z"] 
    vx = data["vx"] 
    vy = data["vy"] 
    vz = data["vz"]
    q = M2/M1

    print("F,N,M,m,a,e")
    print(f"{o.filename},{len(M1)},{np.mean(M1)},{np.mean(M2)},{np.mean(q)},{np.mean(a)},{np.mean(e)}")
    
    plt.scatter(a, e)
    plt.semilogx()
    plt.xlabel("a [au]")
    plt.ylabel("e")
    plt.show()

    a = sorted(np.log10(a))
    f = np.linspace(0, 1, len(a))
    plt.plot(a, f, label="a/au")
    e = sorted(e)
    f = np.linspace(0, 1, len(e))
    plt.plot(e, f, label="e")
    plt.legend()
    plt.show()
    plt.scatter(M1, M2)
    plt.xlabel("M")
    plt.ylabel("m")
    plt.show()
    plt.scatter(M1, q)
    plt.xlabel("M")
    plt.ylabel("q")
    plt.show()
    q = sorted(q)
    f = np.linspace(0, 1, len(q))
    plt.plot(q, f, label="q")
    plt.legend()
    plt.show()
    plt.scatter(x, y)
    #for i in range(len(data)):
    #    plt.arrow(x[i], y[i], dx=vx[i], dy=vy[i])
    plt.xlabel("x [pc]")
    plt.ylabel("y [pc]")
    plt.show()
