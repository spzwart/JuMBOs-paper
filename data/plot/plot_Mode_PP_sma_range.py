import numpy as np
from amuse.units import units
from amuse.plot import plot
from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib import colormaps
import matplotlib as mpl
import pandas as pd

if __name__ in ('__main__', '__plot__'):
    filename = "../reduced/Model_SPP_sma_range.csv"

    #print("F,amin,amax,Nss,Nsp,Nj,M,dM,m,dm,q,dq,a,da,e,de")
    all_data = pd.read_csv(filename)
    data = all_data
    #data = d[d['F'].str.contains("R0.25")]
    
    print(data)
    fname = data["F"] 
    n_jumbos = data["Npp"]
    n_ffp = data["Np"]
    n_ffsp = data["Nsp"]
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
    #q_jumbos = data["q"] 
    #dq_jumbos = data["dq"] 

    print(a_jumbos)
    print(amin.index)
    for i in data.index:
        plt.plot([amin[i], amax[i]],
                 [n_jumbos[i], n_jumbos[i]])
    plt.xlabel("a [au]")    
    plt.ylabel("$n_{jumbos}$")    
    plt.show()    

    for i in data.index:
        ai = (amin[i]+amax[i]/2)
        dai = (amax[i]-amin[i]/2)
        print(i, fname[i], a_jumbos[i])
        aj = (a_jumbos[i]+a_jumbos[i])/2
        daj = (da_jumbos[i]/2)
        plt.errorbar(ai, aj, xerr=dai, yerr=daj)
    plt.show()

    figure = plt.figure(figsize=(6, 4))
    ax = figure.gca()
    #fig, ax = plt.subplots()

    ax.tick_params(which='both', width=1)
    ax.tick_params(which='major', length=7)
    ax.tick_params(which='minor', length=4)
    
    f_jumbos = np.array(n_jumbos)/(np.array(n_jumbos)+np.array(n_ffp))#+np.array(n_ffsp))
    loga = np.log10(1+a_jumbos)
    loga[0] = 2
    loga[1] = 2
    
    print(a_jumbos)

    c = (all_data['a']-40)/(400-40)
    #c = a_jumbos

    #s = 100*(np.array(n_jumbos)/35)**2
    s = np.array(n_jumbos)
    colors = plt.cm.winter_r(c)

    #import matplotlib as m
    #cm = m.colors.LinearSegmentedColormap('my_colormap', cdict, 1024)

    for i in data.index:
        print(i, len(a_jumbos))
        if "R0.25" in fname[i]:
            m = 'o'
        elif "R0.5" in fname[i]:
            m = '^'
        else:
            m = 's'
        plt.plot([amin[i], amax[i]],
                 [f_jumbos[i], f_jumbos[i]], marker = m, lw=1+s[i]/3,
                 linestyle="-", c=cm.winter_r(c[i]), markersize=s[i])
    print("print final jumbo semi-major axis: N=", len(n_jumbos))
    ar = []
    da = []
    n = 0
    for i in data.index:
        ar.append((amax[i]-amin[i])*n_jumbos[i])
        da.append(a_jumbos[i]*n_jumbos[i])
        n+=n_jumbos[i]
        print(f"A=[{amax[i]-amin[i]}] = {a_jumbos[i]}+-{da_jumbos[i]} (N={n_jumbos[i]})")
    print(f"semimajor-axis distribution: {np.sum(ar)/n}, {np.sum(da)/n}, N={n}") 
    print(f"semimajor-axis dispersion: {np.std(ar)}, {np.std(da)}, N={n}") 
        
    #cmap = ax.scatter(amin, f_jumbos, s=s, c=colors, vmin=160, vmax=400)
    #cmap = ax.scatter(amax, f_jumbos, s=s, c=colors)
    
    #[f_jumbos[i], f_jumbos[i]], linestyle="-", c=c[i])
    #colormap = plt.cm.get_cmap('winter_r') # 'plasma' or 'viridis'
    colormap = mpl.colormaps['winter_r']
    sm = plt.cm.ScalarMappable(cmap=colormap)
    sm.set_clim(vmin=40, vmax=400)
    cb = plt.colorbar(sm, ax=ax)
    cb.set_label(r'$\langle a/au \rangle$', labelpad=0, y=0.48)

    """
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = plt.colorbar(cax=cax)
    cbar.set_label(r'$\langle a/au \rangle$')
    """

    #from mpl_toolkits.axes_grid1 import make_axes_locatable
    #divider = make_axes_locatable(ax)
    #cax = divider.append_axes("right", size="5%", pad=0.05)
    #plt.colorbar(sm, cax=cax)
    
    ax.set_xlabel("a [au]")    
    ax.set_ylabel("$f_{jumbos}$")
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

    #plt.xlim(0, 1000)
    #plt.ylim(-0.001, 0.01)
    plt.savefig("fig_fjumbos_from_psystems.pdf")
    plt.show()
