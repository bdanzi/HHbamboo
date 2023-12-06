#!/usr/bin/env python
# 2016MC_PUscenario taken from https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_20_patchX/SimGeneral/MixingModule/python/mix_2016_25ns_Moriond17MC_PoissonOOTPU_cfi.py#L25
# 2017MC_PUscenario taken from https://github.com/cms-sw/cmssw/blob/CMSSW_9_4_X/SimGeneral/MixingModule/python/mix_2017_25ns_WinterMC_PUScenarioV1_PoissonOOTPU_cfi.py#L13
# 2018MC_PUscenatio taken from https://github.com/cms-sw/cmssw/blob/CMSSW_10_4_X/SimGeneral/MixingModule/python/mix_2018_25ns_JuneProjectionFull18_PoissonOOTPU_cfi.py#L11


# To obtain DataPileup root histogram run pileupCalc.py:

# For 2016 data (Moriond 2017 setup:)
# pileupCalc.py -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/pileup_latest.txt --calcMode true --minBiasXsec 69200 --maxPileupBin 75 --numPileupBins 75 DataPileupHistogram2016_69200_75bins.root
# pileupCalc.py -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/pileup_latest.txt --calcMode true --minBiasXsec 66017 --maxPileupBin 75 --numPileupBins 75 DataPileupHistogram2016_66017_75bins.root
# pileupCalc.py -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/pileup_latest.txt --calcMode true --minBiasXsec 72383 --maxPileupBin 75 --numPileupBins 75 DataPileupHistogram2016_72383_75bins.root

#Up and down variation histograms are created varying the minBiasXsec by 4.6% (72383, 66017)

# For 2017 data:
#pileupCalc.py -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PileUp/pileup_latest.txt --calcMode true --minBiasXsec 69200 --maxPileupBin 100 --numPileupBins 100 DataPileupHistogram2017_69200_100bins.root
#pileupCalc.py -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PileUp/pileup_latest.txt --calcMode true --minBiasXsec 66017 --maxPileupBin 100 --numPileupBins 100 DataPileupHistogram2017_66017_100bins.root
#pileupCalc.py -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PileUp/pileup_latest.txt --calcMode true --minBiasXsec 72383 --maxPileupBin 100 --numPileupBins 100 DataPileupHistogram2017_72383_100bins.root

# For 2018 data: 
#pileupCalc.py -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/ReReco/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PileUp/pileup_latest.txt --calcMode true --minBiasXsec 69200 --maxPileupBin 100 --numPileupBins 100 DataPileupHistogram2018_69200_100bins.root
#pileupCalc.py -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/ReReco/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PileUp/pileup_latest.txt --calcMode true --minBiasXsec 66017 --maxPileupBin 100 --numPileupBins 100 DataPileupHistogram2018_66017_100bins.root
#pileupCalc.py -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/ReReco/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PileUp/pileup_latest.txt --calcMode true --minBiasXsec 72383 --maxPileupBin 100 --numPileupBins 100 DataPileupHistogram2018_72383_100bins.root

############## UltraLegacy ###############
#2018UL
#pileupCalc.py -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/UltraLegacy/pileup_latest_2018UL.txt --calcMode true --minBiasXsec 69200 --maxPileupBin 80 --numPileupBins 80 DataPileupHistogram2018UL_69200_80bins.root
#pileupCalc.py -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/UltraLegacy/pileup_latest_2018UL.txt --calcMode true --minBiasXsec 66017 --maxPileupBin 80 --numPileupBins 80 DataPileupHistogram2018UL_66017_80bins.root
#pileupCalc.py -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/UltraLegacy/pileup_latest_2018UL.txt --calcMode true --minBiasXsec 72383 --maxPileupBin 80 --numPileupBins 80 DataPileupHistogram2018UL_72383_80bins.root

#pileupCalc.py -i Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt --inputLumiJSON pileup_latest_2018UL.txt --calcMode true --minBiasXsec 69200 --maxPileupBin 80 --numPileupBins 80 DataPileupHistogram2018UL_69200_80bins.root
#pileupCalc.py -i Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt --inputLumiJSON pileup_latest_2018UL.txt --calcMode true --minBiasXsec 66017 --maxPileupBin 80 --numPileupBins 80 DataPileupHistogram2018UL_66017_80bins.root
#pileupCalc.py -i Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt --inputLumiJSON pileup_latest_2018UL.txt --calcMode true --minBiasXsec 72383 --maxPileupBin 80 --numPileupBins 80 DataPileupHistogram2018UL_72383_80bins.root

#2017UL
#pileupCalc.py -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/Legacy_2017/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PileUp/UltraLegacy/pileup_latest_2017UL.txt --calcMode true --minBiasXsec 69200 --maxPileupBin 80 --numPileupBins 80 DataPileupHistogram2017UL_69200_80bins.root
#pileupCalc.py -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/Legacy_2017/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PileUp/UltraLegacy/pileup_latest_2017UL.txt --calcMode true --minBiasXsec 66017 --maxPileupBin 80 --numPileupBins 80 DataPileupHistogram2017UL_66017_80bins.root
#pileupCalc.py -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/Legacy_2017/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PileUp/UltraLegacy/pileup_latest_2017UL.txt --calcMode true --minBiasXsec 72383 --maxPileupBin 80 --numPileupBins 80 DataPileupHistogram2017UL_72383_80bins.root

#pileupCalc.py -i Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt --inputLumiJSON pileup_latest_2017UL.txt --calcMode true --minBiasXsec 69200 --maxPileupBin 80 --numPileupBins 80 DataPileupHistogram2017UL_69200_80bins.root
#pileupCalc.py -i Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt --inputLumiJSON pileup_latest_2017UL.txt --calcMode true --minBiasXsec 66017 --maxPileupBin 80 --numPileupBins 80 DataPileupHistogram2017UL_66017_80bins.root
#pileupCalc.py -i Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt --inputLumiJSON pileup_latest_2017UL.txt --calcMode true --minBiasXsec 72383 --maxPileupBin 80 --numPileupBins 80 DataPileupHistogram2017UL_72383_80bins.root


#2016UL
#pileupCalc.py -i Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/UltraLegacy/pileup_latest_2016UL.txt --calcMode true --minBiasXsec 69200 --maxPileupBin 80 --numPileupBins 80 DataPileupHistogram2016UL_69200_80bins.root
#pileupCalc.py -i Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/UltraLegacy/pileup_latest_2016UL.txt --calcMode true --minBiasXsec 66017 --maxPileupBin 80 --numPileupBins 80 DataPileupHistogram2016UL_66017_80bins.root
#pileupCalc.py -i Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/UltraLegacy/pileup_latest_2016UL.txt --calcMode true --minBiasXsec 72383 --maxPileupBin 80 --numPileupBins 80 DataPileupHistogram2016UL_72383_80bins.root

#pileupCalc.py -i Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt --inputLumiJSON pileup_latest_2016UL.txt --calcMode true --minBiasXsec 69200 --maxPileupBin 80 --numPileupBins 80 DataPileupHistogram2016UL_69200_80bins.root
#pileupCalc.py -i Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt --inputLumiJSON pileup_latest_2016UL.txt --calcMode true --minBiasXsec 66017 --maxPileupBin 80 --numPileupBins 80 DataPileupHistogram2016UL_66017_80bins.root
#pileupCalc.py -i Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt --inputLumiJSON pileup_latest_2016UL.txt --calcMode true --minBiasXsec 72383 --maxPileupBin 80 --numPileupBins 80 DataPileupHistogram2016UL_72383_80bins.root

# Afetrwards run this script to produce root file which contains PU weights


import ROOT as rt

### MC pu scenario to be used
#puMCscenario = puMC['2016MC_PUscenario']
#puMCscenario = puMC['2017MC_PUscenario']
len_mc = 80

#--- 2016 data
data_file_name = 'DataPileupHistogram2016UL_69200_80bins.root'
data_file_name_varUp = 'DataPileupHistogram2016UL_72383_80bins.root'
data_file_name_varDn = 'DataPileupHistogram2016UL_66017_80bins.root'
mc_file_name = 'pileup_UL_2016.root'
#--- 2017 data
#data_file_name       = 'DataPileupHistogram2017UL_69200_80bins.root'
#data_file_name_varUp = 'DataPileupHistogram2017UL_72383_80bins.root'
#data_file_name_varDn = 'DataPileupHistogram2017UL_66017_80bins.root'
#mc_file_name = 'pileup_UL_2017.root'
#--- 2018 data
#data_file_name       = 'DataPileupHistogram2018UL_69200_80bins.root'
#data_file_name_varUp = 'DataPileupHistogram2018UL_72383_80bins.root'
#data_file_name_varDn = 'DataPileupHistogram2018UL_66017_80bins.root'
#mc_file_name = 'pileup_UL_2018.root'

rt.TH1.SetDefaultSumw2(True)

h_d = rt.TH1F('Data', '', len_mc , 0, len_mc) 
h_d_varUp = rt.TH1F('Data_varUp', '', len_mc , 0, len_mc) 
h_d_varDn = rt.TH1F('Data_varDn', '', len_mc , 0, len_mc) 


fpu = rt.TFile.Open(data_file_name,'read')
h_din = fpu.Get('pileup')
for i in range(1, len_mc + 1) :
    h_d.SetBinContent(i, h_din.GetBinContent(i))
h_d.Scale(1./h_d.Integral())
fpu.Close()

fpu = rt.TFile.Open(data_file_name_varUp,'read')
h_din = fpu.Get('pileup')
for i in range(1, len_mc + 1) :
    h_d_varUp.SetBinContent(i, h_din.GetBinContent(i))
h_d_varUp.Scale(1./h_d_varUp.Integral())
fpu.Close()

fpu = rt.TFile.Open(data_file_name_varDn,'read')
h_din = fpu.Get('pileup')
for i in range(1, len_mc + 1) :
    h_d_varDn.SetBinContent(i, h_din.GetBinContent(i))
h_d_varDn.Scale(1./h_d_varDn.Integral())
fpu.Close()


#h_mc = rt.TH1F('MC_out_of_the_box', ';true number of interactions;normalized to unity', len_mc , 0, len_mc)
#for ipu in range(len(puMCscenario)) :
#    puMCscenario[ipu]
#    h_mc.SetBinContent(ipu + 1, puMCscenario[ipu])

#h_mc.Scale(1./h_mc.Integral())

fpu = rt.TFile.Open(mc_file_name,'read')
h_mc = fpu.Get('MC_out_of_the_box')
fpu.Close()

h_w = h_d.Clone('weights')
h_w.Divide(h_mc)


h_w_varUp = h_d_varUp.Clone('weights_varUp')
h_w_varUp.Divide(h_mc)

h_w_varDn = h_d_varDn.Clone('weights_varDn')
h_w_varDn.Divide(h_mc)

h_mc_rw = h_mc.Clone('MC_reweighted')
h_mc_rw_varUp = h_mc.Clone('MC_up')
h_mc_rw_varUp.SetTitle("MC reweighted +1#sigma")
h_mc_rw_varDn = h_mc.Clone('MC_reweighted')
h_mc_rw_varDn.SetTitle("MC reweighted -1#sigma")


for i in range(1, len_mc + 1) :
    h_mc_rw.SetBinContent(i, h_mc.GetBinContent(i)*h_w.GetBinContent(i))
    h_mc_rw_varUp.SetBinContent(i, h_mc.GetBinContent(i)*h_w_varUp.GetBinContent(i))
    h_mc_rw_varDn.SetBinContent(i, h_mc.GetBinContent(i)*h_w_varDn.GetBinContent(i))

can = rt.TCanvas('can', 'can', 400, 400)
h_mc.SetLineColor(rt.kBlue)
h_mc.Draw()
#h_mc.GetYaxis().SetRangeUser
h_mc.SetMaximum(0.12)
h_mc.Draw()
h_d.SetLineColor(rt.kRed-2)
h_d.SetFillColor(rt.kRed-2)
h_d.SetFillStyle(3004)
h_d.Draw('HISTSAME')
h_mc_rw.SetLineColor(rt.kGreen)
h_mc_rw.Draw('SAME')

h_mc_rw_varUp.SetLineColor(rt.kGreen + 2)
h_mc_rw_varUp.SetLineStyle(2)
h_mc_rw_varUp.Draw('SAME')

h_mc_rw_varDn.SetLineColor(rt.kGreen + 2)
h_mc_rw_varDn.SetLineStyle(3)
h_mc_rw_varDn.Draw('SAME')

leg = can.BuildLegend()
leg.Draw('SAME')

#f_out = rt.TFile.Open('pu_weights_2016.root', 'recreate')
#f_out = rt.TFile.Open('pu_weights_2017.root', 'recreate')
f_out = rt.TFile.Open('pu_weights_2018UL.root', 'recreate')

h_w.Write()
h_w_varDn.Write()
h_w_varUp.Write()

h_mc.Write()
h_d.Write()
h_mc_rw.Write()
can.Write()
f_out.Close()
