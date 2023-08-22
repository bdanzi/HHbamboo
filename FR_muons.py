import pyarrow.parquet as pq
import pandas as pd


PLOT_FOLDER = "/eos/user/a/atalierc/skim_bamboo_FR_mu_ele_tau/results/"

Mumu = pq.read_table(PLOT_FOLDER+"Mu_normal_FR.parquet")
df_Mumu = Mumu.to_pandas()

print(df_Mumu.process)
print("finito pandas")


df_data = df_Mumu.loc[df_Mumu.process.str.contains("data")]

df_WZ = df_Mumu.loc[df_Mumu.process.str.contains("VV")]


df_data_all = df_data.loc[df_data.mt < 40]
df_wz_all = df_WZ.loc[df_WZ.mt < 40]


df_data_all_tight = df_data_all.loc[(df_data_all.lepton_loose_id) & (df_data_all.lepton_iso < 0.25)]
#add bkg


input_vars_=[
    
            ["lepton_pt",  "lep $p_T$ [GeV]", (0, 200), 50],    
]


import matplotlib
import matplotlib.pyplot as plt

fig, ax = plt.subplots()

list_edge = [0,25,50,75,100,150,200]

list_data_all_loose = []
list_data_all_tight = []

for var in input_vars_:

    n_, bin_, _ = plt.hist(df_data_all[var[0]], range=var[2], alpha=1,log=False,bins=list_edge, histtype=("step"),  label = 'mu loose')
    n_tight, bin_tight, _ = plt.hist(df_data_all_tight[var[0]], range=var[2], alpha=1,log=False,bins=list_edge, histtype=("step"),  label = 'mu tight')
    
    print(n_)
    print(n_tight)
    
    list_data_all_loose = n_
    list_data_all_tight = n_tight

#do the same for mc
plt.clf()

num = list_data_all_tight# - mc
den = list_data_all_loose# - mc
import numpy as np


list_edge = [0,25,50,75,100,150,200]
x = [ (list_edge[i+1]-list_edge[i])/2+list_edge[i] for i in range(0, len(list_edge)-1)]
xerr = [ (list_edge[i]-list_edge[i-1])/2 for i in range(1, len(list_edge))]
print(x)



import math
xerr_list=[0,0,0,0,0,0]

fr_loose, bin_loose, _ = plt.errorbar(x, num/den, xerr=xerr, yerr=num/den*(np.sqrt(num)/num + np.sqrt(den)/den), fmt='bo', barsabove=True, ms=0.005, elinewidth=1, label='FR mu loose WP')


plt.xlabel(f'%s'%var[1])
plt.ylabel(f'Events')
plt.legend(loc='upper right')

plt.text(0.0, 1.05,'CMS', transform=ax.transAxes, fontweight = 'bold', fontsize=12)
plt.text(0.1, 1.05,'Run2', transform=ax.transAxes,fontstyle = 'italic', fontsize=12)
plt.text(1, 1.05,'59 $fb^{-1}$ (13 TeV)', transform=ax.transAxes,ha = 'right', fontsize=10)
plt.show()

fr_loose.get_ydata()





