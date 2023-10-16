from amuse.lab import units
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

def read_JuMBOs_Observations(filename="JuMBOs_Observations.csv"):
    data = pd.read_csv("JuMBOs_Observations.csv")
    print(data)
    d = data["ProjSep"]
    #print(pd.DataFrame(d, columns = ['ProjSep']))
    d = np.sort(d)
    scale = 25|units.au/d[0]
    d *= scale
    return d

if __name__=="__main__":
    d = read_JuMBOs_Observations(filename="JuMBOs_Observations.csv")
    f = np.arange(0, 1, 1./len(d))
    plt.plot(d.value_in(units.au), f)
    plt.show()
