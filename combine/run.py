import os
list = [
    ['600_100_scr', 'ZHcand_m_HZA_signal_600_100', '100'],
    ['600_150_scr', 'ZHcand_m_HZA_signal_600_150', '150'],
    ['600_200_scr', 'ZHcand_m_HZA_signal_600_200', '200'],
    ['600_250_scr', 'ZHcand_m_HZA_signal_600_250', '250'],
    ['600_300_scr', 'ZHcand_m_HZA_signal_600_300', '300'],
    ['600_350_scr', 'ZHcand_m_HZA_signal_600_350', '350'],
    ['600_400_scr', 'ZHcand_m_HZA_signal_600_400', '400']
]


for element in list:
    os.system('python harvester.py --input input_folder_multiclass/ --reg lowp --tag '+element[0]+' --signal '+element[1]+'')
    os.system('python harvester.py --input input_folder_multiclass/ --reg highp --tag '+element[0]+' --signal '+element[1]+'')

for element in list:
    os.system('combineCards.py LIMITS/'+element[0]+'/lowp/*txt LIMITS/'+element[0]+'/highp/*txt > full_'+element[0]+'.txt')

for element in list:
    os.system('combine -M AsymptoticLimits --run blind --mass '+element[2]+' full_'+element[0]+'.txt')

#collect into a json file
os.system('combineTool.py -M CollectLimits higgsCombineTest.AsymptoticLimits.mH*.root --use-dirs -o hza_limits.json')
