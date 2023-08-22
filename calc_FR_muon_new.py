import pyarrow.parquet as pq
import pandas as pd


PLOT_FOLDER = "/eos/home-a/atalierc/skim_CR_FR_good/results/"

Mumu = pq.read_table(PLOT_FOLDER+"1mu1tau_CR_FR_mu.parquet")
df_Mumu = Mumu.to_pandas()

print(df_Mumu.process)
print("finito pandas")


for i in df_Mumu.process.unique():
    print(i)


#data and bkg
#SM ZH, WWZ, WZZ, ZZZ, ttZ and ZZ
df_data = df_Mumu.loc[(df_Mumu.process == 'data') & (df_Mumu.mt < 40)]
df_bkg = df_Mumu.loc[(df_Mumu.process.str.contains('HZA_signal')) | (df_Mumu.process == 'VVV') | (df_Mumu.process == 'TTZ') | (df_Mumu.process == 'GluGluToContinToZZ') | (df_Mumu.process == 'ZZTo4L')]
df_bkg = df_bkg.loc[(df_bkg.mt < 40)]


#apply FR
def apply_FR(lept_pt,cat):
    weight = -99.
    if('muon' in cat):
        y = [0.17021336, 0.03746735, 0.02450056, 0.01748893, 0.01278937, 0.00402972]
        x = [25,50,75,100,150,200]    
        for i in range(len(x)):
            if(lept_pt <= x[i]):
                weight = y[i]#/(1-y[i])
                break
            if(lept_pt >= 200):
                weight = y[5]#/(1-y[5])
                break
    if('tau' in cat):
        y = [0.05170972, 0.0509316 , 0.05002239, 0.04904342, 0.05169376, 0.05243137]
        x = [25,50,75,100,150,200]    
        for i in range(len(x)):
            if(lept_pt <= x[i]):
                weight = y[i]#/(1-y[i])
                break
            if(lept_pt >= 200):
                weight = y[5]#/(1-y[5])
                break
    if('ele' in cat):
        y = [0.07231289, 0.04937679, 0.046803, 0.03732532, 0.03487372, 0.03024623]
        x = [25,50,75,100,150,200]    
        for i in range(len(x)):
            if(lept_pt <= x[i]):
                weight = y[i]/(1-y[i])
                break
            if(lept_pt >= 200):
                weight = y[5]/(1-y[5])
                break
    if(weight == -99.):
        print('NOT WORKING')
    else:
        return weight


#data with muon good and tau fake -> first contribution in the 3P1F region

df_good_mu_tau_fake = df_data.loc[(df_data.muon_loose_id) & (df_data.muonn_iso < 0.25)]
df_good_mu_tau_fake = df_good_mu_tau_fake.loc[(df_good_mu_tau_fake.lepton_jet_loose) < 8]
df_good_mu_tau_fake['weight_FR'] = df_good_mu_tau_fake.apply(lambda x: apply_FR(x['tau_pt'], 'tau'), axis=1)

#mc with muon good and tau fake -> second contribution in the 3P1F region

df_good_mu_tau_fake_bkg = df_bkg.loc[(df_bkg.muon_loose_id) & (df_bkg.muonn_iso < 0.25)]
df_good_mu_tau_fake_bkg = df_good_mu_tau_fake_bkg.loc[df_good_mu_tau_fake_bkg.lepton_jet_loose < 8]
df_good_mu_tau_fake_bkg['weight_FR'] = df_good_mu_tau_fake_bkg.apply(lambda x: apply_FR(x['tau_pt'],'tau')*x['weight'], axis=1)

CR_true_mu_false_tau = df_good_mu_tau_fake['weight_FR'].sum()-df_good_mu_tau_fake_bkg['weight_FR'].sum()


#add data with muon fake and tau good -> 3 contribution in the 3P1F region
#add mc with muon fake and tau good -> 3 contribution in the 3P1F region
#add 2P2F region

#ZX yield
#CR_true_mu_false_tau+CR_true_tau_false_mu-CR_2false


#plot in CR10

df_CR_10 = df_Mumu.loc[(df_Mumu.muon_loose_id) & (df_Mumu.muonn_iso < 0.25) & (df_Mumu.mt < 40)]
df_CR_10 = df_CR_10.loc[(df_CR_10.lepton_jet_loose) < 8]


df_CR_10['weight_FR'] = df_CR_10.apply(lambda x: (apply_FR(x['muon_pt'], 'muon'))*(1-apply_FR(x['tau_pt'], 'tau'))*x['weight'], axis=1)

import numpy as np

vgamma = np.array(df_CR_10.loc[df_CR_10.process.str.contains("VGamma")]['muon_pt'])
vgamma_w = np.array(df_CR_10.loc[df_CR_10.process.str.contains("VGamma")]['weight_FR'])

wjets = np.array(df_CR_10.loc[df_CR_10.process.str.contains("WJets")]['muon_pt'])
wjets_w = np.array(df_CR_10.loc[df_CR_10.process.str.contains("WJets")]['weight_FR'])

ttw = np.array(df_CR_10.loc[df_CR_10.process.str.contains("TTW")]['muon_pt'])
ttw_w = np.array(df_CR_10.loc[df_CR_10.process.str.contains("TTW")]['weight_FR'])

ttz = np.array(df_CR_10.loc[df_CR_10.process.str.contains("TTZ")]['muon_pt'])
ttz_w = np.array(df_CR_10.loc[df_CR_10.process.str.contains("TTZ")]['weight_FR'])

ttj = np.array(df_CR_10.loc[df_CR_10.process.str.contains("TTJets")]['muon_pt'])
ttj_w = np.array(df_CR_10.loc[df_CR_10.process.str.contains("TTJets")]['weight_FR'])

dy = np.array(df_CR_10.loc[df_CR_10.process.str.contains("DY")]['muon_pt'])
dy_w = np.array(df_CR_10.loc[df_CR_10.process.str.contains("DY")]['weight_FR'])

h = np.array(df_CR_10.loc[df_CR_10.process.str.contains("GluGluToContinToZZ") | df_CR_10.process.str.contains("ZZTo4L")]['muon_pt']) 
h_w = np.array(df_CR_10.loc[df_CR_10.process.str.contains("GluGluToContinToZZ") | df_CR_10.process.str.contains("ZZTo4L")]['weight_FR'])

st = np.array(df_CR_10.loc[df_CR_10.process.str.contains("SingleT")]['muon_pt'])
st_w = np.array(df_CR_10.loc[df_CR_10.process.str.contains("SingleT")]['weight_FR'])

vv = np.array(df_CR_10.loc[df_CR_10.process.str.contains("VV")]['muon_pt']) 
vv_w = np.array(df_CR_10.loc[df_CR_10.process.str.contains("VV")]['weight_FR']) 

vvv = np.array(df_CR_10.loc[df_CR_10.process.str.contains("VVV")]['muon_pt']) 
vvv_w = np.array(df_CR_10.loc[df_CR_10.process.str.contains("VVV")]['weight_FR']) 

data = np.array(df_CR_10.loc[df_CR_10.process.str.contains("data")]['muon_pt'])
data_w = np.array(df_CR_10.loc[df_CR_10.process.str.contains("data")]['weight_FR'])


import numpy as np
import matplotlib
import matplotlib.pyplot as plt


fig, ax = plt.subplots()

range_ = (0, 1000)
bins_ = 1


plt.hist([ttw,st,wjets,ttz,h,vgamma,vv,vvv,ttj,dy], weights = [ttw_w,st_w,wjets_w,ttz_w,h_w,vgamma_w,vv_w,vvv_w,ttj_w,dy_w],label=('TTW','STop','Wjets','TTZ','Higgs','Vgamma','VV','VVV','ttj','DY'), histtype=("stepfilled"),stacked = True, alpha = 0.7, bins=bins_, range=range_);
plt.hist(df_CR_10.loc[df_CR_10.process == 'data']['muon_pt'], weights = df_CR_10.loc[df_CR_10.process == 'data']['weight_FR'], label='data', histtype=("step"), bins=bins_, range=range_, linewidth=1.5);


plt.legend(loc = 'best')
#plt.xlim(100.,180.)
plt.ylim(10, 1e+3)
plt.ylabel('Events / 25 GeV')
plt.xlabel(r'Z mass [GeV]')
#plt.yscale('log')
plt.text(0.0, 1.01,'CMS', transform=ax.transAxes, fontweight = 'bold', fontsize=12)
plt.text(0.1, 1.01,'HL-LHC', transform=ax.transAxes,fontstyle = 'italic', fontsize=12)
plt.text(1, 1.01,'3000 $fb^{-1}$ (14 TeV)', transform=ax.transAxes,ha = 'right', fontsize=10)
plt.show()
