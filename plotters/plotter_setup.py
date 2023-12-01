from amuse.lab import *
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

class PlotterSetup(object):
    def __init__(self):
        plt.rcParams["font.family"] = "Times New Roman"
        plt.rcParams["mathtext.fontset"] = "cm"
        self.axlabel_size = 14
        self.tick_size = 14
        
        left = .45
        width = .5
        bottom = .58
        height = .5
        self.right = left + width
        self.top = bottom + height

    def tickers(self, ax, ptype):
        """Function to setup axis"""

        ax.yaxis.set_ticks_position('both')
        ax.xaxis.set_ticks_position('both')
        ax.xaxis.set_minor_locator(mtick.AutoMinorLocator())
        ax.yaxis.set_minor_locator(mtick.AutoMinorLocator())
        if ptype == "hist":
            ax.tick_params(axis="y", labelsize = self.tick_size)
            ax.tick_params(axis="x", labelsize = self.tick_size)
            return ax
        else:
            ax.tick_params(axis="y", which = 'both', direction="in", labelsize = self.tick_size)
            ax.tick_params(axis="x", which = 'both', direction="in", labelsize = self.tick_size)
            return ax

    def model_layout(self, model_choices):
        if model_choices == [0,1,2]:
            labels = ["$\mathcal{ISF}$_Fr050", "$\mathcal{ISF}$_Fr050ff", "$\mathcal{ISF}$_Fr100"]
            extra_str = "Fractal_General_"
            
        elif model_choices == [0,3]:
            labels = ["$\mathcal{ISF}$_Fr050", "$\mathcal{ISF}$_Pl050"]
            extra_str = "Distr_0.5_noFF_"

        elif model_choices == [1,10]:
            labels = ["$\mathcal{ISF}$_Fr050ff", "$\mathcal{ISF}$_Fr050ffO"]
            extra_str = "Frac_FF_GenObs_"

        elif model_choices == [3,4,5]:
            labels = ["$\mathcal{ISF}$_Pl050", "$\mathcal{ISF}$_Pl050ff", "$\mathcal{ISF}$_Pl100"]
            extra_str = "Plummer_General_"

        elif model_choices == [0,6]:
            labels = ["$\mathcal{ISF}$_Fr050", "$\mathcal{ISF}$_Fr050ffL"]
            extra_str = "SimTime_"

        elif model_choices == [1,7,10]:
            labels = ["$\mathcal{ISF}$_Fr050ff", "$\mathcal{ISF}$_Fr050OC",
                      "$\mathcal{ISF}$_Fr050ffO"]
            extra_str = "Fractal_FF_"

        elif model_choices == [1,7]:
            labels = ["$\mathcal{ISF}$_Fr050ff", "$\mathcal{ISF}$_Fr050OC"]
            extra_str = "Fractal_FFvsJFF_"

        elif model_choices == [1,4,10,11]:
            labels = ["$\mathcal{ISF}$_Fr050ff", "$\mathcal{ISF}$_Pl050ff",
                      "$\mathcal{ISF}$_Fr050ffO", "$\mathcal{ISF}$_Pl050ffO"]
            extra_str = "Distr_FF_GenObs_"
        elif model_choices == [0,8,9]:
            labels = ["$\mathcal{ISF}$_Fr050", "$\mathcal{ISF}$_Fr050O", "$\mathcal{ISF}$_Fr050OC"]
            extra_str = "Fractal_GenObs_"
        else:
            labels = ["$\mathcal{ISF}$_Fr050ffL"]
            extra_str = "Frac_FF_10Myr"
        return labels, extra_str