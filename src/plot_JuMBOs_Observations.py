from amuse.lab import units
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

def plot_mass_function(masses, ximf, color, label):
    import math
    mass_min = masses.min()
    mass_max = masses.max()
    lm = math.log10(mass_min.value_in(units.MSun))
    lM = math.log10(mass_max.value_in(units.MSun))
    bins = 10**np.linspace(lm, lM, 10)
    print(bins)
    bin_number, bin_edges = np.histogram(
        masses.value_in(units.MSun), bins=bins
    )
    y = bin_number / (bin_edges[1:] - bin_edges[:-1])
    x = (bin_edges[1:] + bin_edges[:-1]) / 2.0
    for i in range(len(y)):
        y[i] = max(y[i], 1.e-10)

    plt.scatter(x, y, s=50, c=color, lw=0, label=label)

    if ximf<0:
        c = (
            (mass_max.value_in(units.MSun)**(ximf+1))
            - (mass_min.value_in(units.MSun)**(ximf+1))
        ) / (ximf+1)
        plt.plot(x, len(masses)/c * (x**ximf), c=color)


# data from https://arxiv.org/pdf/2310.01231.pdf
def read_JuMBOs_Observations(filename="JuMBOs_Observations.csv"):
    data = pd.read_csv("JuMBOs_Observations.csv")
    print(data)
    d = data["ProjSep"]
    m1 = data["MPri"]
    m2 = data["MSec"]
    d, m1, m2 = zip(*sorted(zip(d, m1, m2)))

    #print(pd.DataFrame(d, columns = ['ProjSep']))
    #d = np.sort(d)
    #scale = 25|units.au/d[0]
    d = d | units.au
    m1 = m1 | units.MSun
    m2 = m2 | units.MSun
    return d, m1, m2

if __name__=="__main__":
    d, m1, m2 = read_JuMBOs_Observations(filename="JuMBOs_Observations.csv")
    f = np.arange(0, 1, 1./len(d))
    plt.plot(d.value_in(units.au), f)
    plt.xlabel("d [au]")
    plt.show()

    print(f"Mean separation d= {np.mean(d.value_in(units.au))}, {np.std(d.value_in(units.au))}")
    print(f"Mean mass d= {np.mean(m1.value_in(units.MJupiter))}, {np.std(m1.value_in(units.MJupiter))}")
    print(f"Mean secondary mass d= {np.mean(m2.value_in(units.MJupiter))}, {np.std(m2.value_in(units.MJupiter))}")


    plt.hist(m2/m1)
    plt.xlabel("q")
    plt.show()

    plt.hist((m2+m1).value_in(units.MJupiter))
    plt.xlabel("M+m")
    plt.show()

    plot_mass_function(m1+m2, -1.2, 'r', "JMF")
    plt.loglog()
    plt.legend()
    plt.show()
    
    m1 = sorted(m1.value_in(units.MSun))
    f = np.arange(0, 1, 1./len(m1))
    plt.plot(m1, f, c='b', label="m")
    m2 = sorted(m2.value_in(units.MSun))
    f = np.arange(0, 1, 1./len(m2))
    plt.plot(m2, f, c='r', label="m")
    plt.legend()
    plt.show()
    
