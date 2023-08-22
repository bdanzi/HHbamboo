from sklearn.model_selection import train_test_split
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
import pyarrow.parquet as pq
import matplotlib.pyplot as plt
import math
import numpy as np
import ROOT
from root_numpy import fill_hist

# Import the library
import argparse
# Create the parser
parser = argparse.ArgumentParser()
# Add an argument
parser.add_argument("--cat", type = str, help="Categories: 1mu1tau, 1ele1tau, 2tau", required=True)

# Parse the argument
args = parser.parse_args()
cat =  args.cat


PLOT_FOLDER = "./HZA_plots_dnn/"+cat+"/"

skim = pq.read_table("./parquet/"+cat+"_looseCuts.parquet")
df1 = skim.to_pandas()


df1["label"] = 0
df1.loc[df.process.str.contains("HZA_signal_600_100"), ['label']] = 1
df1.loc[df.process.str.contains("HZA_signal_600_150"), ['label']] = 2
df1.loc[df.process.str.contains("HZA_signal_600_200"), ['label']] = 3
df1.loc[df.process.str.contains("HZA_signal_600_250"), ['label']] = 4
df1.loc[df.process.str.contains("HZA_signal_600_300"), ['label']] = 5
df1.loc[df.process.str.contains("HZA_signal_600_350"), ['label']] = 6
df1.loc[df.process.str.contains("HZA_signal_600_400"), ['label']] = 7


print(df1.columns)

#add your input vars
input_vars = ["met_pt"]

print("INPUT VAR for DNN")
print(input_vars)

#SIGNAL REGION
#add your signal region
sel_tightZ = (df1.lep1_Zcand_pfRelIso04_all < 0.15 ) & (df1.lep1_Zcand_tightId) & (df1.lep2_Zcand_pfRelIso04_all < 0.15 ) & (df1.lep2_Zcand_mediumId) 
sel_tau2   = (df1.tau2_Hcand_iddeepTauVSmu >= 1 )  & (df1.tau2_Hcand_iddeepTauVSe >= 1 ) & (df1.tau2_Hcand_iddeepTauVSjet >= 1 )
if cat == "1mu1tau": 
   sel_tau1 = (df1.tau1_Hcand_mediumId ) & (df1.tau1_Hcand_pfRelIso04_all < 0.25)
elif cat == "1ele1tau":
   sel_tau1 = (df1.tau1_Hcand_mvaFall17V2noIso_WP90) #& (df1.tau1_Hcand_pfRelIso03_all < 0.25)
elif cat == "2tau":
   sel_tau1 = (df1.tau1_Hcand_iddeepTauVSmu >= 1 )  & (df1.tau1_Hcand_iddeepTauVSe >= 1 ) &  (df1.tau1_Hcand_iddeepTauVSjet >= 1 )

OS = df1.tau1_Hcand_ch != df1.tau2_Hcand_ch
onlyMC = (df1.process != "data")

selection = sel_tightZ & sel_tau1 & sel_tau2 & onlyMC & OS

df= df1.loc[selection, :]

onlyData = (df1.process == "data")
if cat == "1mu1tau":
   sel_tau1_NOT = (df1.tau1_Hcand_mediumId ==False) & (df1.tau1_Hcand_pfRelIso04_all > 0.25)
elif cat == "1ele1tau":
   sel_tau1_NOT = (df1.tau1_Hcand_mvaFall17V2noIso_WP90 == False) #& (df1.tau1_Hcand_pfRelIso03_all > 0.25)
elif cat == "2tau":
   sel_tau1_NOT = (df1.tau1_Hcand_iddeepTauVSmu >= 1 )  & (df1.tau1_Hcand_iddeepTauVSe >= 1 ) &  (df1.tau1_Hcand_iddeepTauVSjet < 1 )

selectionNOT = sel_tightZ & sel_tau1_NOT & sel_tau2 & onlyData & OS
df_data = df1.loc[selectionNOT,:]

#import models
model_dnn = keras.models.load_model(PLOT_FOLDER+'dnn_multiclass')
model_dnn.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])

df['class_bkg'] = 1-model_dnn.predict(df[input_vars])[:, 0]
df['class_signal_600_100'] = 1-model_dnn.predict(df[input_vars])[:, 1]
df['class_signal_600_150'] = 1-model_dnn.predict(df[input_vars])[:, 2]
df['class_signal_600_200'] = 1-model_dnn.predict(df[input_vars])[:, 3]
df['class_signal_600_250'] = 1-model_dnn.predict(df[input_vars])[:, 4]
df['class_signal_600_300'] = 1-model_dnn.predict(df[input_vars])[:, 5]
df['class_signal_600_350'] = 1-model_dnn.predict(df[input_vars])[:, 6]
df['class_signal_600_400'] = 1-model_dnn.predict(df[input_vars])[:, 7]


#####DATA
df_data['class_bkg'] = 1-model_dnn.predict(df_data[input_vars])[:, 0]

for s in df.loc[df.process.str.contains("HZA_signal_600_100") == True,:].process.unique():
    plt.hist(df.loc[df.process == s,:]["class_signal_600_100"], histtype=("step"), range=(0.,1.),  bins=25, label = s )

plt.hist(df.loc[df.process.str.contains("HZA_signal") == False,:]["class_bkg"], histtype=("stepfilled"), range=(0., 1.),  bins=25, label = "Background", alpha = 0.5 )
plt.yscale('log')
plt.legend(loc='best')
plt.savefig(PLOT_FOLDER+'h_dnn_score.pdf')
plt.clf()

#Optimize the choice of the cut on the DNN score

x_max = 0.6

print("Optimal cut on dnn score for division in regions: ",x_max)

##SAVE HISTOGRAMS 
EDGES_H  = np.arange(0,500, 25, dtype=float)
EDGES_ZH = np.arange(0,700, 25, dtype=float)
 

list =[ 
    ['HZA_signal_600_100','class_signal_600_100'],
    ['HZA_signal_600_150','class_signal_600_150'],
    ['HZA_signal_600_200','class_signal_600_200'],
    ['HZA_signal_600_250','class_signal_600_250'],
    ['HZA_signal_600_300','class_signal_600_300'],
    ['HZA_signal_600_350','class_signal_600_350'],
    ['HZA_signal_600_400','class_signal_600_400'],
]

########FR DICT
#add your FR values
fr_dict = {"1mu1tau" : 52,
           "1ele1tau": 12,
           "2tau"    : 99}


##LOW PURITY
out_file = ROOT.TFile(PLOT_FOLDER+'/FR/h_mH_mZH_lowp'+str(cat)+'.root',"RECREATE")

for element in list:
    histo_H = ROOT.TH1D("Hcand_m_"+str(element[0]) ,"Hcand_m_"+str(element[0]),len(np.array(EDGES_H))-1,np.array(EDGES_H))
    histo_H.Sumw2()
    fill_hist(histo_H, (df.loc[(df.process.str.contains(element[0])) & ( df[element[1]]<= x_max) ,:]['Hcand_m']).to_numpy(),
                  weights=(df.loc[(df.process.str.contains(element[0])) & ( df[element[1]]<= x_max) ,:]['weight']).to_numpy())
    for ibin in range(1,histo_H.GetNbinsX()+1):
        if histo_H.GetBinContent(ibin) < 0:
            print("sono entrato")
            histo_H.SetBinContent(ibin,1e-3)

    histo_ZH = ROOT.TH1D("ZHcand_m_"+str(element[0]) ,"ZHcand_m_"+str(element[0]),len(np.array(EDGES_ZH))-1,np.array(EDGES_ZH))
    histo_ZH.Sumw2()
    fill_hist(histo_ZH, (df.loc[(df.process.str.contains(element[0])) & ( df[element[1]]<= x_max) ,:]['ZHcand_m']).to_numpy(),
              weights=(df.loc[(df.process.str.contains(element[0])) & ( df[element[1]]<= x_max) ,:]['weight']).to_numpy())
    for ibin in range(1,histo_ZH.GetNbinsX()+1):
        if histo_ZH.GetBinContent(ibin) < 0:
            print("sono entrato")
            histo_ZH.SetBinContent(ibin,1e-3)
    histo_H.Write()
    histo_ZH.Write()


for item in df.process.unique():
    if('signal' not in item):
        histo_H = ROOT.TH1D("Hcand_m_"+str(item) ,"Hcand_m_"+str(item),len(np.array(EDGES_H))-1,np.array(EDGES_H))
        histo_H.Sumw2()
        fill_hist(histo_H, (df.loc[(df.process == item) & ( df.class_bkg<= x_max) ,:]['Hcand_m']).to_numpy(),
                  weights=(df.loc[(df.process == item) & ( df.class_bkg<= x_max) ,:]['weight']).to_numpy())
        for ibin in range(1,histo_H.GetNbinsX()+1):
            if histo_H.GetBinContent(ibin) < 0:
                print("sono entrato")
                histo_H.SetBinContent(ibin,1e-3)

        histo_ZH = ROOT.TH1D("ZHcand_m_"+str(item) ,"ZHcand_m_"+str(item),len(np.array(EDGES_ZH))-1,np.array(EDGES_ZH))
        histo_ZH.Sumw2()
        fill_hist(histo_ZH, (df.loc[(df.process == item) & ( df.class_bkg<= x_max) ,:]['ZHcand_m']).to_numpy(),
                  weights=(df.loc[(df.process == item) & ( df.class_bkg<= x_max) ,:]['weight']).to_numpy())
        for ibin in range(1,histo_ZH.GetNbinsX()+1):
            if histo_ZH.GetBinContent(ibin) < 0:
                print("sono entrato")
                histo_ZH.SetBinContent(ibin,1e-3)
        histo_H.Write()
        histo_ZH.Write()

        #apply FR
        if('ZZTo4L' in item):
            histo_H = ROOT.TH1D("Hcand_m_FR" ,"Hcand_m_FR",len(np.array(EDGES_H))-1,np.array(EDGES_H))
            histo_H.Sumw2()
            fill_hist(histo_H, (df_data.loc[( df_data.class_bkg<= x_max) ,:]['Hcand_m']).to_numpy(),
                      weights=(df_data.loc[ ( df_data.class_bkg<= x_max) ,:]['weight']).to_numpy())
            for ibin in range(1,histo_H.GetNbinsX()+1):
                if histo_H.GetBinContent(ibin) < 0:
                    print("sono entrato")
                    histo_H.SetBinContent(ibin,1e-3)
                    
            histo_ZH = ROOT.TH1D("ZHcand_m_FR" ,"ZHcand_m_FR",len(np.array(EDGES_ZH))-1,np.array(EDGES_ZH))
            histo_ZH.Sumw2()
            fill_hist(histo_ZH, (df_data.loc[(df_data.class_bkg<= x_max) ,:]['ZHcand_m']).to_numpy(),
                      weights=(df_data.loc[(df_data.class_bkg<= x_max) ,:]['weight']).to_numpy())
            for ibin in range(1,histo_ZH.GetNbinsX()+1):
                if histo_ZH.GetBinContent(ibin) < 0:
                    print("sono entrato")
                    histo_ZH.SetBinContent(ibin,1e-3)

            histo_H.Scale(fr_dict[cat]/histo_H.Integral())
            histo_ZH.Scale(fr_dict[cat]/histo_ZH.Integral())

            histo_H.Write()
            histo_ZH.Write()



data_obs = ROOT.TH1D( "data_obs" ,"data_obs",1,0.,1.)
data_obs.Fill(1)
data_obs.Write()
out_file.Close()


##HIGH PURITY
out_file2 = ROOT.TFile(PLOT_FOLDER+'/FR/h_mH_mZH_highp'+str(cat)+'.root',"RECREATE")

for element in list:
    histo_H_2 = ROOT.TH1D("Hcand_m_"+str(element[0]) ,"Hcand_m_"+str(element[0]),len(np.array(EDGES_H))-1,np.array(EDGES_H))
    histo_H_2.Sumw2()
    fill_hist(histo_H_2, (df.loc[(df.process.str.contains(element[0])) & ( df[element[1]]> x_max) ,:]['Hcand_m']).to_numpy(),
                  weights=(df.loc[(df.process.str.contains(element[0])) & ( df[element[1]]> x_max) ,:]['weight']).to_numpy())
    for ibin in range(1,histo_H_2.GetNbinsX()+1):
        if histo_H_2.GetBinContent(ibin) < 0:
            print("sono entrato")
            histo_H_2.SetBinContent(ibin,1e-3)

    histo_ZH_2 = ROOT.TH1D("ZHcand_m_"+str(element[0]) ,"ZHcand_m_"+str(element[0]),len(np.array(EDGES_ZH))-1,np.array(EDGES_ZH))
    histo_ZH_2.Sumw2()
    fill_hist(histo_ZH_2, (df.loc[(df.process.str.contains(element[0])) & ( df[element[1]]> x_max) ,:]['ZHcand_m']).to_numpy(),
              weights=(df.loc[(df.process.str.contains(element[0])) & ( df[element[1]]> x_max) ,:]['weight']).to_numpy())
    for ibin in range(1,histo_ZH_2.GetNbinsX()+1):
        if histo_ZH_2.GetBinContent(ibin) < 0:
            print("sono entrato")
            histo_ZH_2.SetBinContent(ibin,1e-3)
    histo_H_2.Write()
    histo_ZH_2.Write()

for item in df.process.unique():
    if('signal' not in str(item)):
        histo_H_2 = ROOT.TH1D("Hcand_m_"+str(item) ,"Hcand_m_"+str(item),len(np.array(EDGES_H))-1,np.array(EDGES_H))
        histo_H_2.Sumw2()
        fill_hist(histo_H_2, (df.loc[(df.process == item) &( df.class_bkg> x_max) ,:]['Hcand_m']).to_numpy(),
                  weights=(df.loc[(df.process == item) & (df.class_bkg> x_max),:]['weight']).to_numpy())
        for ibin in range(1,histo_H_2.GetNbinsX()+1):
            if histo_H_2.GetBinContent(ibin) < 0:
                print("sono entrato")
                histo_H_2.SetBinContent(ibin,1e-3)

        histo_ZH_2 = ROOT.TH1D("ZHcand_m_"+str(item) ,"ZHcand_m_"+str(item),len(np.array(EDGES_ZH))-1,np.array(EDGES_ZH))
        histo_ZH_2.Sumw2()
        fill_hist(histo_ZH_2, (df.loc[(df.process == item) & (df.class_bkg> x_max),:]['ZHcand_m']).to_numpy(),
                  weights=(df.loc[(df.process == item) & (df.class_bkg> x_max),:]['weight']).to_numpy())
        for ibin in range(1,histo_ZH_2.GetNbinsX()+1):
            if histo_ZH_2.GetBinContent(ibin) < 0:
                print("sono entrato")
                histo_ZH_2.SetBinContent(ibin,1e-3)

        histo_H_2.Write()
        histo_ZH_2.Write()

        #apply FR
        if('ZZTo4L' in item):
            histo_H = ROOT.TH1D("Hcand_m_FR" ,"Hcand_m_FR",len(np.array(EDGES_H))-1,np.array(EDGES_H))
            histo_H.Sumw2()
            fill_hist(histo_H, (df_data.loc[(df_data.class_bkg> 0) ,:]['Hcand_m']).to_numpy(),
                      weights=(df_data.loc[( df_data.class_bkg> 0) ,:]['weight']).to_numpy())
            for ibin in range(1,histo_H.GetNbinsX()+1):
                if histo_H.GetBinContent(ibin) < 0:
                    print("sono entrato")
                    histo_H.SetBinContent(ibin,1e-3)

            histo_ZH = ROOT.TH1D("ZHcand_m_FR" ,"ZHcand_m_FR",len(np.array(EDGES_ZH))-1,np.array(EDGES_ZH))
            histo_ZH.Sumw2()
            fill_hist(histo_ZH, (df_data.loc[( df_data.class_bkg>0) ,:]['ZHcand_m']).to_numpy(),
                      weights=(df_data.loc[( df_data.class_bkg> 0) ,:]['weight']).to_numpy())
            for ibin in range(1,histo_ZH.GetNbinsX()+1):
                if histo_ZH.GetBinContent(ibin) < 0:
                    print("sono entrato")
                    histo_ZH.SetBinContent(ibin,1e-3)

            histo_H.Scale(fr_dict[cat]/histo_H.Integral())
            histo_ZH.Scale(fr_dict[cat]/histo_ZH.Integral())

            histo_H.Write()
            histo_ZH.Write()


data_obs_2 = ROOT.TH1D( "data_obs" ,"data_obs",1,0.,1.)
data_obs_2.Fill(1)
data_obs_2.Write()
out_file2.Close()

