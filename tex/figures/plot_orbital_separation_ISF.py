import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import iqr

sma_ISF_Pl_R050 = [36.56657529,  58.02021065,  61.60197848,  72.04547175,  91.58728582
,  92.31469155,  95.6128096, 99.07622804, 103.11351672, 113.80453782
, 124.6990304,  125.82308741, 128.55473301, 140.18360905, 148.09731408
, 190.92117824, 198.30238643, 206.72865147, 229.15349937, 233.39782822
, 258.39657966, 266.08351559, 272.20186237, 280.36160748, 312.43800673
, 321.28821853, 324.68875108, 360.35482219, 367.66210772, 374.54544409
, 378.95082429, 412.90975329, 419.41339564, 434.62205627, 464.98785845
, 476.45823847, 516.36611098, 568.50645258, 592.51284362, 609.2956579
, 639.68115048, 645.86700783, 692.75336978, 700.68425096, 706.74326507
, 758.10210933, 952.92346468, 989.53223176]

sma_ISF_Fr_t020Myr = [25.13205738, 28.64257288,  31.83956996,  34.38640628,  34.58028167
,  40.27892813,  41.10085002,  42.33359957,  45.48206101,  49.97169816
,  52.89809377,  54.87410959,  57.38573183,  65.25434526,  65.90330072
,  66.10346271,  71.87536154,  76.67255995,  76.67307729,  77.24709588
,  79.9936153, 88.42910723,  88.57655303,  97.36543678, 101.16196131
, 101.36953312, 101.78958185, 102.92135051, 105.08462145, 106.30826295
, 108.95878271, 118.07598438, 118.50326231, 120.08353802, 127.31427931
, 135.91358431, 138.33679293, 139.91693276, 140.1243508,  146.9814593
, 148.58509957, 152.07763788, 152.22258227, 166.63567843, 167.47704658
, 170.85948254, 175.99719898, 183.04979846, 183.35224921, 191.17020877
, 191.93044962, 193.7966243,  196.96741143, 197.89171734, 197.90795757
, 201.01087469, 208.02569426, 209.92506137, 210.86995,   218.44777112
, 224.04457253, 226.76921756, 229.06796919, 230.17775743, 230.21054105
, 231.89382436, 233.25176484, 233.88386456, 237.38952414, 239.56455441
, 246.55393031, 252.32431127, 257.71608185, 282.16198379, 282.22895114
, 286.0297906,  286.95240854, 290.34320419, 294.68131282, 294.78045174
, 295.74130703, 309.02678538, 309.87978887, 318.97157329, 321.13570285
, 337.27544923, 337.46237425, 342.25638984, 344.72706091, 348.71676545
, 355.48014736, 356.02824146, 356.16681375, 356.42001793, 361.11964563
, 366.79070079, 382.60830286, 386.39319862, 392.79234014, 582.8429529 ]

from scipy.stats import ks_2samp


# data from https://arxiv.org/pdf/2310.01231.pdf
def read_JuMBOs_Observations(filename="JuMBOs_Observations.csv"):
    data = pd.read_csv("../observations/JuMBOs_Observations.csv")
    print(data)
    d = data["ProjSep"]
    m1 = data["MPri"]
    m2 = data["MSec"]
    d, m1, m2 = zip(*sorted(zip(d, m1, m2)))

    #print(pd.DataFrame(d, columns = ['ProjSep']))
    #d = np.sort(d)
    #scale = 25|units.au/d[0]
    d = d #| units.au
    m1 = m1 #| units.MSun
    m2 = m2 #| units.MSun
    return d, m1, m2

if __name__=="__main__":

    plt.figure(figsize=(6, 5))
    plt.rcParams.update({'font.size': 14})

    f = np.linspace(0, 1, len(sma_ISF_Fr_t020Myr))
    plt.plot(sma_ISF_Fr_t020Myr, f, c='r', lw=4, label="ISF_Fr_R050 at t=0.2Myr")

    f = np.linspace(0, 1, len(sma_ISF_Pl_R050))
    plt.plot(sma_ISF_Pl_R050, f, c='b', lw=2, label="ISF_Pl_R050")

    d, m1, m2 = read_JuMBOs_Observations(filename="JuMBOs_Observations.csv")
    f = np.arange(0, 1, 1./len(d))
    plt.plot(d, f, c='k', lw=4, label="Oberved JuMBOs")


    print("Fr=", ks_2samp(sma_ISF_Fr_t020Myr, d))
    print("Pl=", ks_2samp(sma_ISF_Pl_R050, d))
    #print(ks_2samp(sma_ISF_Pl_R050,sma_ISF_Fr_t020Myr))

    models = ["../Erwan/Fractal_rvir0.5", "../Erwan/Fractal_rvir0.5_FF",
              "../Erwan/Plummer_rvir0.5", "../Erwan/Plummer_rvir0.5_FF"]

    #fname = "../Erwan/semi_major_Fr_050.h5"
    fname = "../Erwan/semi_major_Fractal_rvir0.5.h5"
    #fname = "../Erwan/semi_major_Plummer_rvir0.5.h5"
    #fname = "../Erwan/semi_major_Fractal_rvir0.5_FF.h5"
    data = pd.read_hdf(fname)
    semi_axis = sorted(data[0])
    a = []
    for ai in range(len(semi_axis)):
        if semi_axis[ai]>25:
            break
    semi_axis = semi_axis[ai:]
    #print(sorted(semi_axis)) #Values are in au    

    
    f = np.arange(0, 1, 1./len(semi_axis))
    plt.plot(semi_axis, f, c='orange', lw=4, label="ISF_Pl_R050")
    
    lg = plt.legend(loc="upper left")
    lg.get_frame().set_alpha(0)
    plt.xlabel("d [au]")
    plt.ylabel("$f_{<d}$")
    plt.semilogx()
    plt.xlim(10, 1000)
    plt.savefig("fig_orbital_separation_ISF.pdf")
    plt.show()

