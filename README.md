![](https://img.shields.io/github/v/tag/bdanzi/HHbbmumu)
[![Run test](https://github.com/bdanzi/HHbbmumu/actions/workflows/python_test.yml/badge.svg)](https://github.com/bdanzi/HHbbmumu/actions/workflows/python_test.yml)
![](https://img.shields.io/badge/CMS-Run2-blue)

# HH(bbmumu) Run-2/3 analysis
---- Work In Progress ----

This repository uses the **bamboo analysis framework**, you can install it via the instructions here: https://bamboo-hep.readthedocs.io/en/latest/install.html#fresh-install

Then clone this repository in the parent directory containing the bamboo installation:

```bash
git clone https://github.com/bdanzi/HHbbmumu.git && cd HHbbmumu
```

Execute these each time you start from a clean shell on lxplus or any other machine with an cvmfs:
```bash
source /cvmfs/sft.cern.ch/lcg/views/LCG_102/x86_64-centos7-gcc11-opt/setup.sh
source (path to your bamboo installation)/bamboovenv/bin/activate
export PYTHONPATH="${PYTHONPATH}:${PWD}/python/"
```

and the followings before submitting jobs to the batch system (HTCondor, Slurm, Dask and Spark are supported):

```bash
voms-proxy-init --voms cms -rfc --valid 192:00 
export X509_USER_PROXY=$(voms-proxy-info -path)
```
if you encounter problems with accessing files when using batch, the following lines may solve your problem

```bash
voms-proxy-init --voms cms -rfc --valid 192:00  --out ~/private/gridproxy/x509
export X509_USER_PROXY=$HOME/private/gridproxy/x509
```

Then plot various control regions via the following command line using batch (you can pass `--maxFiles 1` to use only 1 file from each sample for a quick test):

```bash
bambooRun -m  BaseNanoHHtobbmumu.py:BaseNanoHHtobbmumu samples_2018UL_all.yml --envConfig=../cern.ini -o test --distributed=finalize
```
Instead of passing everytime `--envConfig config/cern.ini`, you can copy the content of that file to `~/.config/bamboorc`.

using the `parquet` output file that contains skims, you can perform machine learning applications
