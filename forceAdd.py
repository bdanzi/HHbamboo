###------------------------------------------------------------------------------------###
###  This code needs envConfig cern.ini                                                ###
###  It forces merging of bath outputs without using --distributed=finalize            ###
###  Author: Brunella D'Anzi                                                           ###
###  Date: 20/08/2023                                                                  ###
###------------------------------------------------------------------------------------###

import os
import argparse
import yaml

parser = argparse.ArgumentParser()
parser.add_argument("--outDir", type = str, help="Output directory of bambooRun", required=True)
parser.add_argument("--yml", type = str, help="yml sample file", required=True)
args   = parser.parse_args()
outDir = args.outDir
yml    = args.yml

DIR    = f"./{outDir}/batch/output/"
ADDIR  = f"./{outDir}/results/"
with open(yml, "r") as stream:
    try:
        dictYML =  yaml.safe_load(stream)
        samples =  dictYML["samples"].keys()
        print(samples)
    except yaml.YAMLError as exc:
        print(exc)

for s in samples:
    print(s)
    print(f'find {DIR} -name "{s}.root" | xargs hadd -f {ADDIR}{s}.root') 
    os.system(f'find {DIR} -name "{s}.root" | xargs hadd -f {ADDIR}{s}.root')
