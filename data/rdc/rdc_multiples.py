import numpy as np
from scipy.stats import iqr
from amuse.units import units
from amuse.plot import plot
from matplotlib import pyplot as plt

import pandas as pd


# find binaries with at least one planet
def split_binaries(binaries, mplanet=50):
    sx = binaries[binaries.M>=mplanet]
    ss = sx[sx.m>=mplanet]
    px = binaries[binaries.M<mplanet]
    pp = px[px.m<mplanet]
    sp = sx[sx.m<mplanet]
    ps = px[px.m>=mplanet]
    n = len(binaries)
    nb = len(ss) + len(sp) + len(ps) + len(pp)
    print("split:", len(ss), len(sp), len(ps), len(pp), "::", n, nb)
    return ss, sp, ps, pp

def find_planet_pairs_in_quadruple_systems(binaries, binbin):
    print("find quadruples.")
    ss,sp,ps,pp = split_binaries(binaries)
    print(len(pp))

    #find pp in quadruples
    ppxx = binbin[binbin.k1.isin(pp.key)]
    xxpp = binbin[binbin.k2.isin(pp.key)]
    
    #find the planet pairs which are member of a quadruple
    ppiqk1 = binaries[binaries.key.isin(ppxx.k1)]
    ppiqk2 = binaries[binaries.key.isin(xxpp.k2)]
    overlap = ppiqk1[ppiqk1.key.isin(ppiqk2.key)]
    print(overlap)
    print("NMultiples=", len(overlap))
    #remove overlap
    ppiqk2 = ppiqk2.drop(overlap.index)
    biq = pd.concat([ppiqk1, ppiqk2])
    print("nq=", len(biq))
    
    pppp = ppxx[ppxx.k2.isin(pp.key)]
    binaries = binaries.drop(ppiqk1.index)
    try:
        binaries = binaries.drop(ppiqk2.index)
    except:
        print("tried")
    print("Nbin=", len(binaries)+len(pppp)+len(ppiqk1)+len(ppiqk2))

    return binaries, biq, xxpp, ppxx, pppp

def find_planet_pairs_in_tripsing_systems(binaries, tripsing):
    print("find hierarchical quadruples.")
    ss,sp,ps,pp = split_binaries(binaries)
    print(len(pp))

    #find pp in quadruples
    ppxx = tripsing[tripsing.k1.isin(pp.key)]
    xxpp = tripsing[tripsing.k2.isin(pp.key)]
    
    #find the planet pairs which are member of a quadruple
    ppiqk1 = binaries[binaries.key.isin(ppxx.k1)]
    ppiqk2 = binaries[binaries.key.isin(xxpp.k2)]
    overlap = ppiqk1[ppiqk1.key.isin(ppiqk2.key)]
    print(overlap)
    print("NMultiples=", len(overlap))
    #remove overlap
    ppiqk2 = ppiqk2.drop(overlap.index)
    biq = pd.concat([ppiqk1, ppiqk2])
    
    pppp = ppxx[ppxx.k2.isin(pp.key)]
    #binaries = binaries.drop(ppiqk1.index)

    print("NN+", len(pppp))
    print("pppp=", len(pppp), pppp.index, pppp.key, pppp.k1, pppp.k2)
    doubles = binaries[binaries.key.isin(pppp.k2)]
    binaries = binaries.drop(doubles.index)
    doubles = binaries[binaries.key.isin(ppiqk1.key)]
    binaries = binaries.drop(doubles.index)
    doubles = binaries[binaries.key.isin(ppiqk2.key)]
    binaries = binaries.drop(doubles.index)
    
    print("Nbin=", len(binaries)+len(pppp)+len(ppiqk1)+len(ppiqk2))

    return binaries, biq, xxpp, ppxx, pppp

def find_planet_pairs_in_triple_systems(binaries, triples):
    print("find triples.")

    ss,sp,ps,pp = split_binaries(binaries)
    print("pp=",len(pp))
    
    #find pp in quadruples
    ppx = triples[triples.k1.isin(pp.key)]
    xpp = triples[triples.k2.isin(pp.key)]

    #find the planet pairs which are member of a triple
    ppitk1 = binaries[binaries.key.isin(ppx.k1)]
    ppitk2 = binaries[binaries.key.isin(xpp.k2)]
    overlap = ppitk1[ppitk1.key.isin(ppitk2.key)]
    print(overlap)
    print("NMultiples=", len(overlap))
    #remove overlap
    ppitk2 = ppitk2.drop(overlap.index)
    
    bit = pd.concat([ppitk1, ppitk2])

    print("Nbinaries in triples=", len(bit))
    
    ppp = ppx[ppx.k2.isin(pp.key)]
    print("NN+", len(ppp))
    print("ppp=", len(ppp), ppp.index, ppp.key, ppp.k1, ppp.k2)
    doubles = binaries[binaries.key.isin(ppp.k2)]
    binaries = binaries.drop(doubles.index)
    doubles = binaries[binaries.key.isin(ppitk1.key)]
    binaries = binaries.drop(doubles.index)
    doubles = binaries[binaries.key.isin(ppitk2.key)]
    binaries = binaries.drop(doubles.index)
    #binaries = binaries.drop(ppp.index)
    #binaries = binaries.drop(ppitk1.index)
    #binaries = binaries.drop(ppitk2.index)
    print("Nbin=", len(binaries)+len(ppp)+len(ppitk1)+len(ppitk2))
    #print(binaries)

    return binaries, bit, xpp, ppx, ppp

def rounded_means_quartiles(data):
    M25 =iqr(data, rng=(25, 50))
    M75 =iqr(data, rng=(50, 75))
    Mm  = np.median(data)
    return Mm, M25, M75

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-f", dest="filename", default = "PlN2500n600cs_R0.25pcQ05_JuMBOs.csv",
                      help="input filename[%default]")
    result.add_option("-n", dest="nrun", type="int", default = "1",
                      help="number of runs [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    #type,key,k1,k2,M,m,a,e,i,dx,dy,dz,x,y,z,vx,vy,vz

    f = open(o.filename)
    r = f.readlines()
    nsnapshot = 0
    for ri in r:
        if "AMUSE" in ri:
            nsnapshot += 1
    f.close()
    nsnapshot = int(nsnapshot/2)
    print("N-snapshot=", nsnapshot)
    if o.nrun != nsnapshot:
        print(f"input nrun={o.nrun} !- {nsnapshot}=nsnapshots")
        print("stop")
        exit(-1)
    
    data = pd.read_csv(o.filename)
    data = data[data["M"]>0.6]
    
    """
    data = data[data["x"]>-0.5]
    data = data[data["x"]<0.5]
    data = data[data["y"]>-0.5]
    data = data[data["y"]<0.5]
    """

    try:
        amin = int(o.filename.split("+")[1])
        amax = int(o.filename.split("+")[2])
    except:
        amin = 0
        amax = 1

    all_binaries = data[data['type'] == 'binary']
    triples = data[data['type'] == 'triple']
    tripsing = data[data['type'] == 'tripsing']
    binbin = data[data['type'] == 'binbin']
    singleton = data[data['type'] == 'singleton']

    print("N=", len(singleton), len(all_binaries), len(triples), len(tripsing), len(binbin))
    
    mplanet = 50 # objects less massive than 50 jupiters is considered a planet

    binaries, biq, xxpp, ppxx, pppp = find_planet_pairs_in_quadruple_systems(all_binaries, binbin)
    print("Binary-binaries:", len(binaries), len(xxpp), len(ppxx), len(pppp))

    binaries, biq, xxpph, ppxxh, pppph = find_planet_pairs_in_tripsing_systems(binaries, tripsing)
    print("Tripsing:", len(binaries), len(xxpph), len(ppxxh), len(pppph))
    
    binaries, bit, xpp, ppx, ppp = find_planet_pairs_in_triple_systems(binaries, triples)
    print("triples:", len(binaries), len(xpp), len(ppx), len(ppp))

    #now analyze the left-over binaries
    ss,sp,ps,pp = split_binaries(binaries)
    print("binaries:", len(ss), len(sp), len(ps), len(pp))

    #true singles
    tss = singleton[singleton['M'] >= mplanet]
    tsp = singleton[singleton['M'] < mplanet]
    tsp = tsp[tsp['M'] >= 1]
    n_s = len(tss)/o.nrun
    n_p = len(tsp)/o.nrun
    
    nspp = 0
    if len(xpp)>0:
        spp = xpp[xpp.M>=mplanet]
        nspp = len(spp)
        print(spp)
    if len(ppx)>0:
        pps = ppx[ppx.m>=mplanet]
        nspp += len(pps)
        print(pps)
    npppph = len(pppph)
    npppp = len(pppp)
    nppp = len(ppp)
    nss = len(ss)/o.nrun
    nsp = (len(sp) + len(ps))/o.nrun
    print("Npp=", len(pp), o.nrun)
    npp = len(pp)/o.nrun

    M1 = pp["M"]
    M2 = pp["m"]
    m1 = []
    m2 = []
    for i in pp.index:
        if M1[i]<M2[i]:
            m1.append(M2[i])
            m2.append(M1[i])
        else:
            m1.append(M1[i])
            m2.append(M2[i])
    
    Mm, M25, M75 = rounded_means_quartiles(m1)
    mm, m25, m75 = rounded_means_quartiles(m2)
    am, a25, a75 = rounded_means_quartiles(pp["a"])
    em, e25, e75 = rounded_means_quartiles(pp["e"])
    #dx, dx25, dx75 = rounded_means_quartiles(pp["dx"])
    #dy, dy25, dy75 = rounded_means_quartiles(pp["dy"])
    #dz, dz25, dz75 = rounded_means_quartiles(pp["dz"])
    dx = pp['dx']
    dy = pp['dy']
    dz = pp['dz']
    dxy = np.sqrt(dx**2+dy**2)
    dxz = np.sqrt(dx**2+dz**2)
    dyz = np.sqrt(dy**2+dz**2)
    print(len(dxy), len(dxz), len(dyz))
    dd = dxy
    dd = np.append(dd, dxz)
    dd = np.append(dd, dyz)
    dd, d25, d75 = rounded_means_quartiles(dd)

    print(f"Nss={len(tss)}, Nst={len(tsp)}")
    print("Table 1")
    print(f"{o.filename} & {n_s} & {nss} & {nsp} & {n_p} & {npp} & {nspp} & {nppp} & {npppph} & {npppp} \\\\")
    print(f"${Mm:2.3f}^{{+{M75:2.3f}}}_{{-{M25:2.3f}}}$ & ${mm:2.3f}^{{+{m75:2.3f}}}_{{-{m25:2.3f}}}$ & ${am:2.3f}^{{+{a75:2.3f}}}_{{-{a25:2.3f}}}$ & ${em:2.3f}^{{+{e75:2.3f}}}_{{-{e25:2.3f}}}$& ${dd:2.3f}^{{+{d75:2.3f}}}_{{-{d25:2.3f}}}$  \\\\")

    print(f"Model {o.filename},{n_s},{nss},{nsp},{n_p},{npp},{nspp},{nppp},{npppph},{npppp},{Mm:2.3f},{M75:2.3f},{mm:2.3f},{m75:2.3f},{am:2.3f},{a75:2.3f},${em:2.3f},{e75:2.3f},{dd:2.3f},{d75:2.3f}")
    

    
