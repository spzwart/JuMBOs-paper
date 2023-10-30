import numpy as np
from amuse.units import units
from amuse.plot import plot
from matplotlib import pyplot as plt

import pandas as pd

if __name__ in ('__main__', '__plot__'):
    filename = "../reduced/Mode_PP_sma_range.csv"

    #print("F,amin,amax,Nss,Nsp,Nj,M,dM,m,dm,q,dq,a,da,e,de")
    all_data = pd.read_csv(filename)
    d = all_data
    #data = d[d['F'].str.contains("R0.25")]
    
    for idx, *row in d.itertuples():
        print(idx)
        print(row)
        print(row[6])
        if row[6]<row[8]:
            m = row[6]
            row[6] = row[8]
            row[8] = m
    data = d
    print(data)
    fname = data["F"] 
    n_jumbos = data["Nj"] 
    n_ffps = data["Nsp"]
    f = data["F"]
    amin = data["amin"]
    amax = data["amax"]
    M1 = data["M"] 
    M2 = data["m"]
    dM1 = data["dM"] 
    dM2 = data["dm"]
    a_jumbos = data["a"] 
    da_jumbos = data["da"] 
    e_jumbos = data["e"] 
    de_jumbos = data["de"] 
    q_jumbos = data["q"] 
    dq_jumbos = data["dq"] 

    print(amin)
    print(amin.index)
    #for idx, *row in d.itertuples():    
    for i in d.index:
        plt.plot([amin[i], amax[i]],
                 [n_jumbos[i], n_jumbos[i]])
    #plot([amin[i], amax[i]], [n_jumbos[i], n_jumbos[i]])
    plt.xlabel("a [au]")    
    plt.ylabel("$n_{jumbos}$")    
    plt.show()    

    for i in d.index:
        ai = (amin[i]+amax[i]/2)
        dai = (amax[i]-amin[i]/2)
        print(i, f[i], a_jumbos[i])
        aj = (a_jumbos[i]+a_jumbos[i])/2
        daj = (da_jumbos[i]/2)
        plt.errorbar(ai, aj, xerr=dai, yerr=daj)
    plt.show()

    figure = plt.figure(figsize=(6, 4))
    ax = figure.gca()
    #fig, ax = plt.subplots()
    
    f_jumbos = np.array(n_jumbos)/(2*np.array(n_jumbos)+np.array(n_ffps))
    loga = np.log10(1+a_jumbos)
    loga[0] = 2
    loga[1] = 2
    """
    c = []
    for i in range(len(n_jumbos)):
        if a_jumbos[i]>=25 and a_jumbos[i]<380:
            c.append("r")
        else:
            c.append("k")
    """
    from matplotlib import cm
    print(a_jumbos)

    c = (all_data['a']-140)/(400-140)
    #c = a_jumbos

    s = 100*(np.array(n_jumbos)/35)**2
    colors = plt.cm.winter_r(c)

    #import matplotlib as m
    #cm = m.colors.LinearSegmentedColormap('my_colormap', cdict, 1024)

    for i in d.index:
        print(i, len(a_jumbos))
        if "R0.25" in fname[i]:
            m = 'o'
        elif "R0.5" in fname[i]:
            m = '^'
        else:
            m = 's'
        plt.plot([amin[i], amax[i]],
                 [f_jumbos[i], f_jumbos[i]], marker = m,
                 linestyle="-", c=cm.winter_r(c[i]))
    print("print final jumbo semi-major axis: N=", len(n_jumbos))
    ar = []
    da = []
    n = 0
    for i in d.index:
        ar.append((amax[i]-amin[i])*n_jumbos[i])
        da.append(a_jumbos[i]*n_jumbos[i])
        n+=n_jumbos[i]
        print(f"A=[{amax[i]-amin[i]}] = {a_jumbos[i]}+-{da_jumbos[i]} (N={n_jumbos[i]})")
    print(f"semimajor-axis distribution: {np.sum(ar)/n}, {np.sum(da)/n}, N={n}") 
    print(f"semimajor-axis dispersion: {np.std(ar)}, {np.std(da)}, N={n}") 

        
        
    #cmap = ax.scatter(amin, f_jumbos, s=s, c=colors, vmin=160, vmax=400)
    #cmap = ax.scatter(amax, f_jumbos, s=s, c=colors)

    #[f_jumbos[i], f_jumbos[i]], linestyle="-", c=c[i])
    colormap = plt.cm.get_cmap('winter_r') # 'plasma' or 'viridis'
    sm = plt.cm.ScalarMappable(cmap=colormap)
    sm.set_clim(vmin=160, vmax=400)
    cb = plt.colorbar(sm)
    cb.set_label(r'$\langle a/au \rangle$', labelpad=-15, y=0.48)
    plt.xlabel("a [au]")    
    plt.ylabel("$f_{jumbos}$")
    #plt.colorbar(cb,cm.hot)
    #pyplot.colorbar(ticks=[1,5,10,20,50], format=formatter)
    #cbar = plt.colorbar()#cbar)#, fraction=0.024, pad=0.02)
    #, ticks=[1,5,10,20,50]) #ax=ax,
    #orientation='vertical',
    #norm=norm,
    #ticks=bounds)

    #import matplotlib as mpl    
    #cbar = mpl.colorbar.ColorbarBase(cmap=cmap, ax=ax)
    #                       norm=mpl.colors.Normalize(vmin=0, vmax=35))

    plt.xlim(0, 3500)
    plt.ylim(-0.002, 0.03)
    plt.savefig("fig_fjumbos_from_psystems.pdf")
    plt.show()
