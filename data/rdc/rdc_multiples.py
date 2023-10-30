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
    try:
        amin = int(o.filename.split("+")[1])
        amax = int(o.filename.split("+")[2])
    except:
        amin = 0
        amax = 1

    """
    for idx, *row in data.itertuples():
        print(idx)
        print(row)
        print(row[6])
        if row[6]<row[8]:
            m = row[6]
            row[6] = row[8]
            row[8] = m   
    """
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

    #mplanet = 0.08
    mplanet = 50
    
    sx = binaries[binaries.M>=mplanet]
    xs = binaries[binaries.m>=mplanet]
    ss = sx[sx.m>=mplanet]
    sp = sx[sx.m<mplanet]
    px = binaries[binaries.M<mplanet]
    xp = binaries[binaries.m<mplanet]
    pp = px[px.m<mplanet]
    ps = px[px.m>mplanet]

    print(f"stars: {len(sx)}, {len(xs)}, ss={len(ss)}, sp={len(sp)}")
    print(f"planets: {len(px)}, {len(xp)},pp= {len(pp)}, ps={len(ps)}")
    n_s = len(sp) + len(ps) + 2*len(ss)
    n_p = len(sp) + len(ps) + 2*len(pp)
    print(f"N(sp + ps)= {len(ps)}+{len(sp)}={len(ps)+len(sp)}, N(s)={n_s}, N(p)={n_p}")
    n_sp = len(ps)+len(sp)

    nstars = 2500
    nplanets = 1200
    nrun = 10

    print("Numbers per cluster")
    print(f"number of single stars: {nstars-n_s/nrun}")
    print(f"number of single planets: {nplanets-n_p/nrun}")
    print(f"number of (s,s): {len(ss)/nrun}")
    print(f"number of (s,p): {n_sp/nrun}")
    print(f"number of (p,p): {len(pp)/nrun}")
    print(f"checksum stars: {nstars-n_s/nrun + 2*len(ss)/nrun + n_sp/nrun}")
    print(f"checksum planets: {nplanets-n_p/nrun + 2*len(pp)/nrun + n_sp/nrun}")

    nss = len(ss)/nrun
    nsp = n_sp/nrun
    npp = len(pp)/nrun

    tbap = triples[triples.k1.isin(binaries.key)]
    n_bs = len(tbap[tbap.m>=mplanet])/nrun
    n_jjs = len(tbap[tbap.M<mplanet])/nrun
    n_bp = len(tbap[tbap.m<mplanet])/nrun
 
    tbas = triples[triples.k2.isin(binaries.key)]
    try:
        n_sb = len(tbas[tbap.M>=mplanet])/nrun
        n_pb = len(tbas[tbap.M<mplanet])/nrun
    except:
        n_sb = len(tbas)/nrun
        n_pb = 0.0

    nbs = n_bs+n_sb
    nbp = n_bp+n_pb
    print(f"stars in triples n(bs)= {nbs}")
    print(f"planets in triples n(bp)= {nbp}")
        
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
    
    print(f"Njumbos (before mass and semimajor axis selection) ={len(jumbos)}")
    #jumbos = binbin
    #print(jumbos)

    jumbos = jumbos[jumbos['M'] > 0.6]
    jumbos = jumbos[jumbos['M'] < 20]
    jumbos = jumbos[jumbos['m'] > 0.6]
    jumbos = jumbos[jumbos['m'] < 20]
    jumbos = jumbos[jumbos['a'] < 10000]

    print(f"Njumbos ={len(jumbos)}")
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

    print("F,amin,amax,Nbs,Nbp,Nss,Nsp,Npp,M,dM,m,dm,q,dq,a,da,e,de")
    print(f"{o.filename},{amin},{amax},{nbs},{nbp},{nss},{nsp},{npp},{np.mean(M1)},{np.std(M1)},{np.mean(M2)},{np.std(M2)},{np.mean(q)},{np.std(q)},{np.mean(a)},{np.std(a)},{np.mean(e)},{np.std(e)}")
