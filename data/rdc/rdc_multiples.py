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
    biq = ppiqk1.append(ppiqk2)
    
    pppp = ppxx[ppxx.k2.isin(pp.key)]
    binaries = binaries.drop(ppiqk1.index)
    binaries = binaries.drop(ppiqk2.index)
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
    bit = ppitk1.append(ppitk2)
    #print(bit)
    
    ppp = ppx[ppx.k2.isin(pp.key)]
    binaries = binaries.drop(ppp.index)
    binaries = binaries.drop(ppitk1.index)
    binaries = binaries.drop(ppitk2.index)
    print("Nbin=", len(binaries)+len(ppp)+len(ppitk1)+len(ppitk2))
    #print(binaries)

    return binaries, bit, xpp, ppx, ppp

def rounded_means_quartiles(data):
    M25 =iqr(data, rng=(25, 50))
    M75 =iqr(data, rng=(50, 75))
    Mm  = np.median(data)
    M25 = round(M25,2)
    M75 = round(M75,2)
    Mm = round(Mm,2)
    return Mm, M25, M75

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-f", dest="filename", default = "PlN2500n600cs_R0.25pcQ05_JuMBOs.csv",
                      help="input filename[%default]")
    result.add_option("-n", dest="nrun", type="int", default = "4",
                      help="number of runs [%default]")
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

    all_binaries = data[data['type'] == 'binary']
    triples = data[data['type'] == 'triple']
    tripsing = data[data['type'] == 'tripsing']
    binbin = data[data['type'] == 'binbin']

    mplanet = 50 # objects less massive than 50 jupiters is considered a planet

    binaries, biq, xxpp, ppxx, pppp = find_planet_pairs_in_quadruple_systems(all_binaries, binbin)
    print(len(binaries), len(xxpp), len(ppxx), len(pppp))
    
    binaries, bit, xpp, ppx, ppp = find_planet_pairs_in_triple_systems(binaries, triples)
    print(len(binaries), len(xpp), len(ppx), len(ppp))

    #now analyze the left-over binaries
    ss,sp,ps,pp = split_binaries(binaries)
    print(len(ss), len(sp), len(ps), len(pp))

    nspp = 0
    if len(xpp)>0:
        spp = xpp[xpp.M>=mplanet]
        nspp = len(spp)
        print(spp)
    if len(ppx)>0:
        pps = ppx[ppx.m>=mplanet]
        nspp += len(pps)
        print(pps)
    npppp = len(ppp)
    nppp = len(ppp)
    nss = len(ss)/o.nrun
    nsp = (len(sp) + len(ps))/o.nrun
    npp = len(pp)/o.nrun
    n_p = 1200 - (3*nppp +2*npp +nsp)

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
    
    print("Table 1")
    print(f"{o.filename} & {nss} & {nsp} & {n_p} & {npp} & {nspp} & {nppp} & {npppp} \\\\")
    print(f"${Mm}^{{+{M75}}}_{{-{M25}}}$ & ${mm}^{{+{m75}}}_{{-{m25}}}$ & ${am}^{{+{a75}}}_{{-{a25}}}$ & ${em}^{{+{e75}}}_{{-{e25}}}$ \\\\")

