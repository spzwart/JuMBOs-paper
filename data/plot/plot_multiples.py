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
    #type,key,k1,k2,M,m,a,e,i,x,y,z,vx,vy,vz
    data = pd.read_csv(o.filename)
    """
    M1 = d["M"] 
    M2 = d["m"]
    for i in range(len(M1)):
        if M1[i]<M2[i]:
            m = M1[i]
            M1[i] = M2[i]
            M2[i] = m
    """

    #data = d

    binaries = data[data['type'] == 'binary']
    triples = data[data['type'] == 'triple']
    tripsing = data[data['type'] == 'tripsing']
    binbin = data[data['type'] == 'binbin']

    binkk = binaries['key']
    bink1 = binaries['k1']
    bink2 = binaries['k2']
    tripk1 = triples['k1']
    tripk2 = triples['k2']
    bbk1 = binbin['k1']
    bbk2 = binbin['k2']
    tsk1 = tripsing['k1']
    tsk2 = tripsing['k2']

    #jumbos should not be member of a triple
    non_jumbos = binaries[binaries.key.isin(triples.k1) | binaries.key.isin(triples.k2)]
    jumbos = binaries.drop(non_jumbos.index)
    print(jumbos)

    #jumbos should not be member of a quadruple
    non_jumbos = jumbos[jumbos.key.isin(tripsing.k1) |jumbos.key.isin(tripsing.k2)]
    jumbos = jumbos.drop(non_jumbos.index)

    #jumbos should not be member of a quadruple
    non_jumbos = jumbos[jumbos.key.isin(binbin.k1) |jumbos.key.isin(binbin.k2)]
    jumbos = jumbos.drop(non_jumbos.index)
    print(jumbos)
    
    print(f"Njumbos={len(jumbos)}")
    #jumbos = binbin
    #print(jumbos)

    jumbos = jumbos[jumbos['M'] > 0.6]
    jumbos = jumbos[jumbos['M'] < 20]
    jumbos = jumbos[jumbos['m'] > 0.6]
    jumbos = jumbos[jumbos['m'] < 20]
    jumbos = jumbos[jumbos['a'] < 10000]

    #jumbos = binbin
    
    type = jumbos["type"] 
    key = jumbos["key"] 
    k1 = jumbos["k1"] 
    k2 = jumbos["k2"] 
    M1 = jumbos["M"] 
    M2 = jumbos["m"] 
    a = jumbos["a"] 
    e = jumbos["e"] 
    i = jumbos["i"] 
    x = jumbos["x"] 
    y = jumbos["y"] 
    z = jumbos["z"] 
    vx = jumbos["vx"] 
    vy = jumbos["vy"] 
    vz = jumbos["vz"]
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
    #for i in range(len(jumbos)):
    #    plt.arrow(x[i], y[i], dx=vx[i], dy=vy[i])
    plt.xlabel("x [pc]")
    plt.ylabel("y [pc]")
    plt.show()
