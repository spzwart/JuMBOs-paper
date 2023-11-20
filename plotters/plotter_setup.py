from amuse.lab import *
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

class PlotterSetup(object):
    def __init__(self):
        plt.rcParams["font.family"] = "Times New Roman"
        plt.rcParams["mathtext.fontset"] = "cm"
        self.axlabel_size = 16
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
            labels = [r"$R_{\mathrm{vir}} = 0.5$ pc", 
                      r"$R_{\mathrm{vir}} = 0.5$ pc, w/ FF", 
                      r"$R_{\mathrm{vir}} = 1.0$ pc"]
            plot_label = r"Fractal"
            extra_str = "Fractal_General_"
            
        elif model_choices == [0,3]:
            labels = [r"Fractal", r"Plummer"]
            plot_label = r"$R_{\mathrm{vir}} = 0.5$ pc"
            extra_str = "Distr_0.5_noFF_"

        elif model_choices == [1,10]:
            labels = [r"General", r"Constrained"]
            plot_label = r"$R_{\mathrm{vir}} = 0.5$ pc"
            extra_str = "FF_GenObs_"

        elif model_choices == [3,4,5]:
            labels = [r"$R_{\mathrm{vir}} = 0.5$ pc", 
                      r"$R_{\mathrm{vir}} = 0.5$ pc, w/ FF", 
                      r"$R_{\mathrm{vir}} = 1.0$ pc"]
            plot_label = r"Plummer"
            extra_str = "Plummer_General_"

        elif model_choices == [0,6]:
            labels = [r"$t_{\mathrm{end}} = 1$ Myr", 
                      r"$t_{\mathrm{end}} = 10$ Myr"]
            plot_label = r"Fractal, $R_{\mathrm{vir}} = 0.5$ pc, w/FF"
            extra_str = "SimTime_"

        elif model_choices == [0,8,9]:
            labels = [r"General", r"Constrained",
                      r"Constrained + Circularised"]
            plot_label = r"Fractal, $R_{\mathrm{vir}} = 0.5$ pc"
            extra_str = "Fractal_noFF_"

        elif model_choices == [1,7,10]:
            labels = [r"JuMBOs+FF", r"FF Only",
                      r"JuMBOs+FF Constrained"]
            plot_label = r"Fractal, $R_{\mathrm{vir}} = 0.5$ pc, w/FF"
            extra_str = "Fractal_FF_"

        elif model_choices == [1,7]:
            labels = [r"JuMBOs+FF", r"FF Only"]
            plot_label = r"Fractal, $R_{\mathrm{vir}} = 0.5$ pc, w/FF"
            extra_str = "Fractal_FFvsJFF_"

        elif model_choices == [1,4,10,11]:
            labels = [r"Fractal General", r"Plummer General",
                      r"Fractal Constrained", r"Plummer Constrained"]
            plot_label = r"Fractal, $R_{\mathrm{vir}} = 0.5$ pc, w/FF"
            extra_str = "Distr_FF_GenObs_"
        else:
            labels = ["Dummy"]
            plot_label = ["Dummy"]
            extra_str = ["Dummy"]
        return labels, plot_label, extra_str