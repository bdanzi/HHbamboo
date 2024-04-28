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
mkdir bamboodev
cd bamboodev
# make a virtualenv
source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc11-opt/setup.sh
python -m venv bamboovenv
source bamboovenv/bin/activate
# clone and install bamboo
git clone -o upstream https://gitlab.cern.ch/cp3-cms/bamboo.git
pip install ./bamboo
# clone and install plotIt
git clone -o upstream https://github.com/cp3-llbb/plotIt.git
mkdir build-plotit
cd build-plotit
cmake -DCMAKE_INSTALL_PREFIX=$VIRTUAL_ENV ../plotIt
make -j2 install
pip install git+https://gitlab.cern.ch/cp3-cms/pyplotit.git
```

and the followings before submitting jobs to the batch system (HTCondor, Slurm, Dask and Spark are supported):

```bash
voms-proxy-init --voms cms -rfc --valid 192:00  --out ~/private/gridproxy/x509
export X509_USER_PROXY=$HOME/private/gridproxy/x509
```

Use condor_submit -spool on /eos/.

Each time you open a new session:

```bash
source bamboovenv/bin/activate
source /cvmfs/cms.cern.ch/cmsset_default.sh
export X509_USER_PROXY=$HOME/private/gridproxy/x509
```

Then plot various control regions via the following command line using batch (you can pass `--maxFiles 1` to use only 1 file from each sample for a quick test):

```bash
bambooRun -m  BaseNanoHHtobbmumu.py:BaseNanoHHtobbmumu samples_2018UL_all.yml --envConfig=../cern.ini -o test --distributed=driver
```
To merge all the ouput files and produce plots:

```bash
bambooRun -m  BaseNanoHHtobbmumu_*.py:BaseNanoHHtobbmumu samples_2018UL_all.yml --envConfig=../cern.ini -o test --distributed=finalize
```
Instead of passing everytime `--envConfig config/cern.ini`, you can copy the content of that file to `~/.config/bamboorc`.

using the `parquet` output file that contains skims, you can perform machine learning applications
