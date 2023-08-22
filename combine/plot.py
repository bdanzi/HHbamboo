import numpy as np
import json
import numpy as np

from matplotlib import pyplot as plt
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Tahoma']
from matplotlib import pyplot as plt

def getArrays(json_dict):
    exp_0  = []
    exp_p1 = []
    exp_m1 = []
    exp_p2 = []
    exp_m2 = []
    #obs    = []

    keys = np.array([float(key) for key in json_dict.keys()])
    keys.sort()
    #keys = keys[:-1]
    for kl in keys:
        limit = json_dict[str(kl)]
        exp_0.append(limit['exp0'])
        exp_p1.append(limit['exp+1'])
        exp_m1.append(limit['exp-1'])
        exp_p2.append(limit['exp+2'])
        exp_m2.append(limit['exp-2'])
        #obs.append(limit['obs'])

    return  exp_0, exp_p1, exp_m1, exp_p2, exp_m2, keys


def klscan_plotting(year, filename):

    lumi = {'2016': '35.8', '2017': '41.5', '2018': '59.7', 'Run2':'59.7'}

    with open(filename) as json_file:
        json_dict = json.load(json_file)
        exp_0, exp_p1, exp_m1, exp_p2, exp_m2, kl = getArrays(json_dict)

        xsec_limit    = exp_0#*xsec_BR
        xsec_limit_p2 = exp_p2#*xsec_BR
        xsec_limit_p1 = exp_p1#*xsec_BR
        xsec_limit_m1 = exp_m1#*xsec_BR
        xsec_limit_m2 = exp_m2#*xsec_BR

        fig = plt.figure(figsize=(10,10))
        ax  = fig.add_subplot(1,1,1)

        print(kl)
        print(xsec_limit)

        ax.fill_between(kl,xsec_limit_p2, color="gold", label=r"2 sigma")
        ax.fill_between(kl,xsec_limit_p1, color="forestgreen", label=r"1 sigma")
        ax.fill_between(kl,xsec_limit_m1, color="gold")
        ax.fill_between(kl,xsec_limit_m2, color="white")

        ax.plot(kl, xsec_limit, color='k', ls='--', label="Expected")



        ax.set_ylim(0,2)
        #ax.set_xlim(-5,10)
        ax.set_ylabel(r"95% CL on $\sigma / \sigma_{SM}$", fontsize="23",fontweight='normal')
        ax.set_xlabel(r"$m_{A}$ (GeV)", fontsize="23", fontweight='normal')

        ax.set_title("%s $fb^{-1}$ (13 TeV)"%(lumi[year]), loc="right", fontdict={'fontsize': "18", 'fontweight':'normal'})


        ax.tick_params(direction='out', length=6, width=1, colors='k',
                       which='major')
        ax.tick_params(direction='out', length=3, width=1, colors='k',
                       which='minor')

        from matplotlib.lines import Line2D
        from matplotlib.patches import Patch

        custom_lines = [Line2D([0], [0], color="red", lw=1, ls="-", label='Theory Prediction'),
                        Line2D([0], [0], color="black", lw=1, ls="--", label='Expected 95% CL Limit'),
                        Patch(facecolor="forestgreen", edgecolor='k',label=r'Expected $\pm$1 s.d.'),
                        Patch(facecolor="gold", edgecolor='k', label=r'Expected $\pm$2 s.d.')]

        ax.legend(handles=custom_lines,fontsize=16)
        plt.grid(ls='--',color='grey')


        plt.savefig("mH_mA_limts.pdf".format(year))


limits_files = {'Run2': 'hza_limits_default.json'}

for year, filename in limits_files.items():
    print(year,filename)
    klscan_plotting(year, filename)
