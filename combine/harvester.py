import CombineHarvester.CombineTools.ch as ch

import ROOT as R
import glob
import os

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--reg", type=str, required=True, choices=['highp', 'lowp'], help="region in mbb" )
parser.add_argument("--tag", type=str, required=True, help="Datacard name tag" )
parser.add_argument("--input", type=str, required=True, help="Root files input path" )
parser.add_argument("--signal", type=str, required=True, help="signal")


args = parser.parse_args()

cb = ch.CombineHarvester()

reg  = args.reg
tag   = args.tag
input = args.input
sig = args.signal
#/lustre/home/taliercio/SL7/CombineNew/CMSSW_10_2_13/src/HiggsAnalysis/CombinedLimit/DM_limits
auxiliaries  = os.environ['CMSSW_BASE'] + '/src/HiggsAnalysis/CombinedLimit/HZA_limits'
aux_shapes   = auxiliaries +'/{inputFolder}'.format(inputFolder=input)

print auxiliaries
print aux_shapes
 #['2016', '2017', '2018']

chns = ['1mu1tau', '1ele1tau', '2tau']


bkg_procs = {
    '1mu1tau'   : ['ZHcand_m_DY', 'ZHcand_m_VVV', 'ZHcand_m_VV', 'ZHcand_m_VGamma', 'ZHcand_m_SingleT', 'ZHcand_m_TTJets', 'ZHcand_m_TTW', 'ZHcand_m_TTZ', 'ZHcand_m_ZZTo4L', 'ZHcand_m_GluGluToContinToZZ', 'ZHcand_m_WJets'],
    '1ele1tau' : ['ZHcand_m_DY', 'ZHcand_m_VVV', 'ZHcand_m_VV', 'ZHcand_m_VGamma', 'ZHcand_m_SingleT', 'ZHcand_m_TTJets', 'ZHcand_m_TTW', 'ZHcand_m_TTZ', 'ZHcand_m_ZZTo4L', 'ZHcand_m_GluGluToContinToZZ', 'ZHcand_m_WJets'],
    '2tau'   : ['ZHcand_m_DY', 'ZHcand_m_VVV', 'ZHcand_m_VV', 'ZHcand_m_VGamma', 'ZHcand_m_SingleT', 'ZHcand_m_TTJets', 'ZHcand_m_TTW', 'ZHcand_m_TTZ', 'ZHcand_m_ZZTo4L', 'ZHcand_m_GluGluToContinToZZ', 'ZHcand_m_WJets']
 }

sig_procs = [str(sig)]#['ZHcand_m_HZA_signal_600_200']

cats = {
  '1mu1tau' : [
    (0, '1mu1tau'),
  ],
  '1ele1tau' : [
    (0, '1ele1tau'),
 ],
  '2tau' : [
    (0, '2tau')
  ]
}

print '>> Creating processes and observations...'

for chn in chns:
    cb.AddObservations(  ['*'],  ['HZA'], ['13TeV'], [chn],                 cats[chn]      )
    cb.AddProcesses(     ['*'],  ['HZA'], ['13TeV'], [chn], bkg_procs[chn], cats[chn], False  )
    cb.AddProcesses(     [''], ['HZA'], ['13TeV'], [chn], sig_procs,      cats[chn], True   )#mod 125



print '>> Adding systematic uncertainties...'
signal = cb.cp().signals().process_set()

MC_Backgrouds = ['ZHcand_m_DY', 'ZHcand_m_VVV', 'ZHcand_m_VV', 'ZHcand_m_VGamma', 'ZHcand_m_SingleT', 'ZHcand_m_TTJets', 'ZHcand_m_TTW', 'ZHcand_m_TTZ', 'ZHcand_m_ZZTo4L', 'ZHcand_m_GluGluToContinToZZ', 'ZHcand_m_WJets']


cb.cp().process(signal+MC_Backgrouds).AddSyst(cb, "lumi_14TeV_HL_LHC", "lnN", ch.SystMap()([1.01,0.99]))
cb.cp().process(signal).AddSyst(cb, "lep_ID_mu", "lnN", ch.SystMap("channel")(['mu'] , [1.01,0.99]))

cb.cp().process(signal).AddSyst(cb, "lep_ID_e", "lnN", ch.SystMap("channel")(['e'] , [1.01,0.99]))

cb.cp().process(signal).AddSyst(cb, "tau_ID_mu", "lnN", ch.SystMap("channel")(['mu'] , [1.05,0.95]))

cb.cp().process(signal).AddSyst(cb, "tau_ID_e", "lnN", ch.SystMap("channel")(['e'] , [1.05,0.95]))
                                                                                        
cb.cp().process(signal).AddSyst(cb, "tau_ID_tau", "lnN", ch.SystMap("channel")(['tau'] , [1.05,0.95]))

cb.cp().process(signal+MC_Backgrouds).AddSyst(cb, "JES", "lnN", ch.SystMap()([1.01,0.99]))
cb.cp().process(signal+MC_Backgrouds).AddSyst(cb, "btag", "lnN", ch.SystMap()([1.02,0.98]))
#cb.cp().process(['mtt_TT_inclusive']).AddSyst(cb, "QCD_tt_tautau", "lnN", ch.SystMap()([1.012,0.982]))
#cb.cp().process(['mtt_TT_inclusive']).AddSyst(cb, "pdf_tt", "lnN", ch.SystMap()([1.021,0.979]))
cb.cp().process(signal).AddSyst(cb, "QCD_sig_tautau", "lnN", ch.SystMap()([1.021,0.951]))
cb.cp().process(signal).AddSyst(cb, "pdf_sig", "lnN", ch.SystMap()([1.03,0.97]))
cb.cp().process(signal).AddSyst(cb, "mtop_un_sig", "lnN", ch.SystMap()([1.027,0.973]))

#h_mH_mZH_highp_2tau.root
print '>> Extracting histograms from input root files...'
for chn in chns:
    file = aux_shapes + "/h_mH_mZH_"+reg+""+chn+".root"
    cb.cp().channel([chn]).era(['13TeV']).backgrounds().ExtractShapes(
        file, '$PROCESS', '$PROCESS_$SYSTEMATIC')
    cb.cp().channel([chn]).era(['13TeV']).signals().ExtractShapes(
        file, '$PROCESS', '$PROCESS$MASS_$SYSTEMATIC')

print '>> Setting standardised bin names...'
ch.SetStandardBinNames(cb)

writer = ch.CardWriter('LIMITS/$TAG/$ANALYSIS_$CHANNEL_$BINID_$ERA.txt',
                       'LIMITS/$TAG/$ANALYSIS_$CHANNEL.input.root')

print(writer)
#HZA_1ele1tau_0_13TeV.txt
writer.SetVerbosity(2)
writer.WriteCards('{}/{}'.format(tag,reg), cb) ## the first argument is the $TAG
# for chn in chns: writer.WriteCards(chn,cb.cp().channel([chn]))
#HHbbtautau_tau_0_13TeV
#LIMITS/600_200/highp/HZA_2tau_0_13TeV.txt
#LIMITS/600_200/lowp/HZA_1ele1tau_0_13TeV.txt
for chn in chns:
     with open("./LIMITS/"+tag+"/"+reg+"/HZA_"+chn+"_0_13TeV.txt",'a') as f:
         #         print("reg ", reg)
         f.write("* autoMCStats 0  1  1")

print '>> Done!'
